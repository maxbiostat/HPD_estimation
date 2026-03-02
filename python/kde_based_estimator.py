from dataclasses import dataclass
from typing import Callable, Dict

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
from icecream import ic
from jax import random
from jax.scipy.optimize import minimize
from jaxtyping import Float, Array
from quadax import quadgk
from scipy.integrate import quad
from scipy.optimize import root, root_scalar, fsolve
from scipy.stats import gaussian_kde
from scipy.stats import Normal


# plt.rcParams.update({
#     "text.usetex": True,
# })
jax.config.update("jax_enable_x64", True)


@dataclass
class RootFinderResult:
    x: Float[Array, "..."]
    success: bool
    message: str | None

@dataclass
class HPDEstimate:
    a_hat: float
    a_L: float
    a_U: float
    b_hat: float
    b_L: float
    b_U: float


def empirical_hpd(data: Float[Array, "n_samples"], alpha: Float[Array, ""]) -> Float[Array, ""]:
    n = data.shape[0]
    ecdf = lambda x: jnp.sum(data <= x) / n
    lower_ecdf = jax.vmap(ecdf)(data)
    upper_ecdf = jax.vmap(lambda p: jnp.min(jnp.array([p + alpha, 1])))(lower_ecdf)
    lower_endpoints = data
    upper_endpoints = jnp.quantile(data, upper_ecdf, method="lower")
    interval_lengths = (upper_endpoints - lower_endpoints) + jnp.inf * (upper_ecdf - lower_ecdf < alpha)
    idx = jnp.argmin(interval_lengths)
    return jnp.array([lower_endpoints[idx], upper_endpoints[idx]])


def compute_bandwidth(data: Float[Array, ""], method="rule of thumb") -> Float[Array, ""]:
    n = data.shape[0]
    if method == "rule of thumb": # Assumes kernel is gaussian
        return jnp.power(4/3, 1/5) / jnp.power(n, 1/4) * jnp.std(data)
    

def KDE(
    data: Float[Array, "n_samples"],
    bandwidth: Float[Array, ""],
    kernel: Callable[[Float[Array, ""]], Float[Array, ""]],
) -> Callable[[Float[Array, ""]], Float[Array, ""]]:
    return lambda x: 1/bandwidth * jnp.mean(kernel((x - data)/bandwidth))


def TransformedKDE(
    kde: Callable[[Float[Array, ""]], Float[Array, ""]],
    phi: Callable[[Float[Array, ""]], Float[Array, ""]],
) -> Callable[[Float[Array, "x"]], Float[Array, "x"]]:
    """
    kde is an estimate of the density of phi(X) and TransformedKDE is an estimate
    of the density of X.
    """
    return lambda x: kde(phi(x)) * jax.grad(phi)(x)


def empirical_Psi(
    interval: Float[Array, "2"],
    phi: Callable[[Float[Array, ""]], Float[Array, ""]],
    transformed_kde: Callable[[Float[Array, ""]], Float[Array, ""]],
    alpha: Float[Array, ""],
) -> Float[Array, "2"]:
    a, b = interval[0], interval[1]
    delta_transformed_kde = transformed_kde(b) - transformed_kde(a)
    return jnp.array([
        quadgk(kde, [phi(a), phi(b)])[0] - alpha, delta_transformed_kde
    ])


def jac_empirical_Psi(
    interval: Float[Array, "2"],
    phi: Callable[[Float[Array, ""]], Float[Array, ""]],
    transformed_kde: Callable[[Float[Array, ""]], Float[Array, ""]],
) -> Float[Array, "2"]:
    a, b = interval[0], interval[1]
    transformed_kde_derivative = jax.grad(transformed_kde)
    result = jnp.array([
        [
            -transformed_kde(a),
            transformed_kde(b)
        ],
        [
            -transformed_kde_derivative(a),
            transformed_kde_derivative(b)
        ]
    ])
    return result


def HPD_estimator(
    data: Float[Array, "n_samples"],
    alpha: Float[Array, ""],
    beta: Float[Array, ""], # confidence level of hpd endpoints
    kernel: Callable[[Float[Array, "..."]], Float[Array, "..."]],
    phi: Callable[[Float[Array, "..."]], Float[Array, "..."]],
    kernel_squared_norm: Float[Array, ""],
    x0: Float[Array, "2"] = None,
    maxiter: int = 1000,
    key: jax.random.PRNGKey = 42,
    single_rate: bool = False,
) -> HPDEstimate:
    n_samples = data.shape[0]
    transformed_data = phi(data)
    bandwidth = compute_bandwidth(transformed_data)
    kde_ess = n_samples * bandwidth
    kde = KDE(transformed_data, bandwidth, kernel)
    transformed_kde = TransformedKDE(kde, phi)
    transformed_kde_deriv = jax.grad(transformed_kde)
    phi_grad = jax.grad(phi)
    normal_dist = Normal(mu=0.0, sigma=1.0)
    q = normal_dist.icdf((1 + jnp.sqrt(beta))/2)

    fun = lambda interval: empirical_Psi(
        interval, phi, transformed_kde,alpha
    )
    jac = lambda interval: jac_empirical_Psi(
        interval, phi, transformed_kde
    )
    if x0 is None:
        x0 = empirical_hpd(data, alpha)
    x1 = x0
    hpd_sol = None
    it = 0
    while (hpd_sol is None) and (it <= maxiter):
        it += 1
        hpd_sol = root(fun, x0=x1, jac=jac)
        if not hpd_sol.success :
            hpd_sol = None
            key, _ = random.split(key)
            x1 = x0 + jnp.array([jnp.quantile(data, (1 - alpha)/2), jnp.quantile(data, (1 + alpha)/2)]) * random.normal(key, (2))
    if hpd_sol is None:
        raise Exception(f"Reached max number of tries: {maxiter}. Root finder did not converge")
    a_hat, b_hat = hpd_sol.x

    if single_rate:
        q = normal_dist.icdf((1 + beta)/2)
        std_err = (
            transformed_kde(a_hat) * jnp.sqrt(
                (kernel_squared_norm * (transformed_kde(a_hat) * phi_grad(a_hat) + transformed_kde(b_hat) * phi_grad(b_hat)))
                / (transformed_kde(b_hat) * transformed_kde_deriv(a_hat) - transformed_kde(a_hat) * transformed_kde_deriv(b_hat))**2
            )
        ) / jnp.sqrt(kde_ess)
        a_L = a_hat - std_err * q
        a_U = a_hat + std_err * q
        b_L = b_hat - std_err * q
        b_U = b_hat + std_err * q
    else:
        det = (
            transformed_kde(b_hat) * transformed_kde_deriv(a_hat)
            - transformed_kde(a_hat) * transformed_kde_deriv(b_hat)
        )
        inv_psi_dot = 1 / det * jnp.array([
            [transformed_kde_deriv(b_hat), - transformed_kde(b_hat)],
            [transformed_kde_deriv(a_hat), - transformed_kde(a_hat)]
        ])
        inv_R_n = jnp.array([
            [1/jnp.sqrt(n_samples), 0],
            [0, 1/jnp.sqrt(kde_ess)]
        ])
        root_cov = jnp.sqrt(jnp.array([
            [alpha * (1 - alpha), 0],
            [0, kernel_squared_norm * (transformed_kde(a_hat) * phi_grad(a_hat) + transformed_kde(b_hat) * phi_grad(b_hat))]
        ]))
        gamma = inv_psi_dot @ inv_R_n @ root_cov

        vertices = gamma @ jnp.stack(
            [
                jnp.array([q, q]),
                jnp.array([-q, q]),
                jnp.array([q, -q]),
                jnp.array([-q, -q])
            ],
            axis=1
        )
        a_L = a_hat + jnp.min(vertices[0, :])
        a_U = a_hat + jnp.max(vertices[0, :])
        b_L = b_hat + jnp.min(vertices[1, :])
        b_U = b_hat + jnp.max(vertices[1, :])
    # return jnp.array([a_hat, a_L, a_U, b_hat, b_L, b_U])
    return HPDEstimate(a_hat, a_L, a_U, b_hat, b_L, b_U)


if __name__ == "__main__":
    kernel = lambda x: 1/jnp.sqrt(2*jnp.pi)*jnp.exp(-x**2/2)
    phi = lambda x: jnp.log(x)
    kernel_squared_norm = 1 / (2*jnp.sqrt(jnp.pi))
    ALPHA = 0.95
    BETA = 0.95
    SIGMA = 1
    MU = 0.0
    N_SAMPLES = int(1E4)
    key = random.key(42)

    def log_normal_pdf(
        x: Float[Array, ""],
    ) -> Float[Array, ""]:
        return jnp.exp(-(jnp.log(x) - MU)**2 / (2 * SIGMA**2)) / (x * SIGMA * jnp.sqrt(2*jnp.pi))

    true_hpd_sol = root(
        fun=lambda x: jnp.array([quad(log_normal_pdf, x[0], x[1])[0] - ALPHA, log_normal_pdf(x[0]) - log_normal_pdf(x[1])]),
        x0=[0.1, 2.0],
        jac=lambda x: jnp.array([
                [-log_normal_pdf(x[0]), log_normal_pdf(x[1])],
                [jax.grad(log_normal_pdf)(x[0]), -jax.grad(log_normal_pdf)(x[1])]
        ]),
        method="hybr",
    )
    print(true_hpd_sol.success)
    print(true_hpd_sol.message)
    print(true_hpd_sol.x)
    assert true_hpd_sol.success
    true_hpd = true_hpd_sol.x
    a, b = true_hpd[0], true_hpd[1]
    ic(quadgk(log_normal_pdf, [a, b]))
    ic(jnp.abs(log_normal_pdf(b) - log_normal_pdf(a)))

    key, _ = random.split(key)
    data = jnp.exp(MU + SIGMA * random.normal(key, shape=(N_SAMPLES,)))
    transformed_data = phi(data)
    bandwidth = compute_bandwidth(transformed_data)
    kde = KDE(transformed_data, bandwidth, kernel)
    transformed_kde = TransformedKDE(kde, phi)

    hpd_estimator = HPD_estimator(
        data, ALPHA, BETA, kernel, phi, kernel_squared_norm
    )
    a_hat = hpd_estimator.a_hat
    a_L = hpd_estimator.a_L
    a_U = hpd_estimator.a_U
    b_hat = hpd_estimator.b_hat
    b_L = hpd_estimator.b_L
    b_U = hpd_estimator.b_U
    ic(empirical_Psi([a_hat, b_hat], phi, transformed_kde, ALPHA))
    ic(a_hat)
    ic(a_L)
    ic(a_U)
    ic(b_hat)
    ic(b_L)
    ic(b_U)
    ic(jnp.abs(transformed_kde(b_hat) - transformed_kde(a_hat)))

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)
    ax_a, ax_b = axs

    # left endpoint plot
    x_a = jnp.linspace(0.5*a, 1.5*a, num=100)
    ax_a.plot(x_a, log_normal_pdf(x_a), linewidth=1, color="black", label=r"$f(x)$")
    ax_a.plot(x_a, jax.vmap(transformed_kde)(x_a), linewidth=1, color="blue", label=r"$f_n(x)$")
    ax_a.scatter(x=a, y=0, color="black", s=10, label=r"$a$")
    ax_a.scatter(x=a_hat, y=0, color="blue", s=10, label=r"$\hat{a}$")
    ax_a.fill_between([a_L, a_U], -0.1*transformed_kde(a_hat), 0.1*transformed_kde(a_hat), color='skyblue', alpha=0.4, label="95% IC")
    ax_a.set_xlim(xmin=0.5*a, xmax=1.5*a)
    ax_a.legend()

    # right endpoint plot
    x_b = jnp.linspace(0.5*b, 1.5*b, num=100)
    ax_b.plot(x_b, log_normal_pdf(x_b), linewidth=1, color="black", label=r"$f(x)$")
    ax_b.plot(x_b, jax.vmap(transformed_kde)(x_b), linewidth=1, color="blue", label=r"$f_n(x)$")
    ax_b.scatter(x=b, y=0, color="black", s=10, label=r"$b$")
    ax_b.scatter(x=b_hat, y=0, color="blue", s=10, label=r"$\hat{b}$")
    ax_b.fill_between([b_L, b_U], -0.1*transformed_kde(b_hat), 0.1*transformed_kde(b_hat), color='skyblue', alpha=0.4, label="95% IC")
    ax_b.set_xlim(xmin=0.5*b, xmax=1.5*b)
    ax_b.legend()

    fig.tight_layout()

    plt.show()

    fig_, ax = plt.subplots()
    x = jnp.linspace(0.0, 5.0, num=300)
    ax.plot(x, log_normal_pdf(x), linewidth=1, color="black", label=r"$f(x)$")
    ax.plot(x, jax.vmap(transformed_kde)(x), linewidth=1, color="blue", label=r"$f_n(x)$")
    plt.show()