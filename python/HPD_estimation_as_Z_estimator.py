from collections.abc import Callable
from dataclasses import dataclass
from typing import Dict

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
from icecream import ic
from jax import random
from jaxtyping import Float, Array
from scipy.integrate import quad
from scipy.optimize import root
from scipy.stats import gaussian_kde
from scipy.stats import Normal


# plt.rcParams.update({
#     "text.usetex": True,
# })

def compute_bandwidth(data: Float[Array, ""], method="rule of thumb") -> Float[Array, ""]:
    n = data.shape[0]
    if method == "rule of thumb": # Assumes kernel is gaussian
        return jnp.power(4/3, 1/5) / jnp.power(n, 1/5) * jnp.std(data)
    

def KDE(
    data: Float[Array, "n_samples"],
    bandwidth: Float[Array, ""],
    kernel: Callable[[Float[Array, "x"]], Float[Array, "x"]],
) -> Callable[[Float[Array, ""]], Float[Array, ""]]:
    return lambda x: 1/bandwidth * jnp.mean(kernel((x - data)/bandwidth))


def Psi(
    interval: Float[Array, "2"],
    kde: Callable[[Float[Array, ""]], Float[Array, ""]],
    kernel_second_moment: Float[Array, ""], # mu_2(K) in the document
    bandwidth: Float[Array, ""],
    alpha: Float[Array, ""],
) -> Float[Array, "2"]:
    a, b = interval[0], interval[1]
    kde_second_derivative = jax.grad(jax.grad(kde))
    delta_kde = kde(b) - kde(a)
    delta_kde_curvature = kde_second_derivative(b) - kde_second_derivative(a)
    return jnp.array([
        quad(kde, a, b)[0] - alpha,
        delta_kde - 0.5 * kernel_second_moment * bandwidth**2 * delta_kde_curvature
    ])


def jac_Psi(
    interval: Float[Array, "2"],
    kde: Callable[[Float[Array, ""]], Float[Array, ""]],
    kernel_second_moment: Float[Array, ""], # mu_2(K) in the document
    bandwidth: Float[Array, ""],
    alpha: Float[Array, ""],
) -> Float[Array, "2"]:
    a, b = interval[0], interval[1]
    kde_derivative = jax.grad(kde)
    kde_third_derivative = jax.grad(jax.grad(kde_derivative))
    return jnp.array([
        [
            -kde(a),
            kde(b)
        ],
        [
            -kde_derivative(a) + 0.5 * kernel_second_moment * bandwidth**2 * kde_third_derivative(a),
            kde_derivative(b) - 0.5 * kernel_second_moment * bandwidth**2 * kde_third_derivative(b),
        ]
    ])


def HPD_estimator_with_confidence(
    data: Float[Array, "n_samples"],
    alpha: Float[Array, ""],
    beta: Float[Array, ""], # confidence level of hpd endpoints
    kernel: Callable[[Float[Array, "..."]], Float[Array, "..."]],
    kernel_second_moment: Float[Array, ""],
    kernel_squared_norm: Float[Array, ""],
) -> Float[Array, "6"]:
    n_samples = data.shape[0]
    bandwidth = compute_bandwidth(data)
    kde_ess = n_samples * bandwidth
    kde = KDE(data, bandwidth, kernel)
    normal_dist = Normal(mu=0.0, sigma=1.0)
    q = normal_dist.icdf(BETA)

    fun = lambda interval: Psi(
        interval, kde, kernel_second_moment, bandwidth, alpha
    )
    jac = lambda interval: jac_Psi(
        interval, kde, kernel_second_moment, bandwidth, alpha
    )
    hpd = root(fun, [0.0, 4.0], jac=jac, method="hybr").x
    a_hat, b_hat = hpd
    kde_deriv = jax.grad(kde)
    gamma_hat = jnp.sqrt(
        (
            kernel_squared_norm * (kde(a_hat) + kde(b_hat))
        )
        / (
            kde(b_hat) * kde_deriv(a_hat) - kde(a_hat) * kde_deriv(b_hat)
        )
    )
    a_L = a_hat - gamma_hat * kde(b_hat) / jnp.sqrt(kde_ess) * q
    a_U = a_hat + gamma_hat * kde(b_hat) / jnp.sqrt(kde_ess) * q
    b_L = b_hat - gamma_hat * kde(a_hat) / jnp.sqrt(kde_ess) * q
    b_U = b_hat + gamma_hat * kde(a_hat) / jnp.sqrt(kde_ess) * q
    return jnp.array([a_hat, a_L, a_U, b_hat, b_L, b_U])


if __name__ == "__main__":
    kernel = lambda x: 1/jnp.sqrt(2*jnp.pi)*jnp.exp(-x**2/2)
    kernel_second_moment = 1
    kernel_squared_norm = 1 / (2*jnp.sqrt(jnp.pi))
    ALPHA = 0.9
    BETA = 0.95
    SIGMA = 1.0
    MU = 0
    N_SAMPLES = int(1E6)
    key = random.key(0)

    def log_normal_pdf(
        x: Float[Array, ""],
    ) -> Float[Array, ""]:
        return jnp.exp(-(jnp.log(x) - MU)**2 / (2 * SIGMA**2)) / (x * SIGMA * jnp.sqrt(2*jnp.pi))

    true_hpd_sol = root(
        fun=lambda x: jnp.array([quad(log_normal_pdf, x[0], x[1])[0] - ALPHA, log_normal_pdf(x[0]) - log_normal_pdf(x[1])]),
        x0=[0.1, 0.4],
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

    key, _ = random.split(key)
    data = jnp.exp(random.normal(key, shape=(N_SAMPLES,)))
    bandwidth = compute_bandwidth(data)
    kde = KDE(data, bandwidth, kernel)

    hpd_estimator = HPD_estimator_with_confidence(
        data, ALPHA, BETA, kernel, kernel_second_moment, kernel_squared_norm
    )
    a_hat, a_L, a_U, b_hat, b_L, b_U = hpd_estimator
    ic(a_hat)
    ic(a_L)
    ic(a_U)
    ic(b_hat)
    ic(b_L)
    ic(b_U)

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)
    ax_a, ax_b = axs

    # left endpoint plot
    x_a = jnp.linspace(0.5*a, 1.5*a, num=100)
    ax_a.plot(x_a, log_normal_pdf(x_a), linewidth=1, color="black", label=r"$f(x)$")
    ax_a.plot(x_a, jax.vmap(kde)(x_a), linewidth=1, color="blue", label=r"$f_n(x)$")
    ax_a.scatter(x=a, y=0, color="blue", s=15, label=r"$a$")
    ax_a.set_xlim(xmin=0.5*a, xmax=1.5*a)
    ax_a.legend()

    # right endpoint plot
    x_b = jnp.linspace(0.5*b, 1.5*b, num=100)
    ax_b.plot(x_b, log_normal_pdf(x_b), linewidth=1, color="black", label=r"$f(x)$")
    ax_b.plot(x_b, jax.vmap(kde)(x_b), linewidth=1, color="blue", label=r"$f_n(x)$")
    ax_b.scatter(x=b, y=0, color="blue", s=15, label=r"$b$")
    ax_b.set_xlim(xmin=0.5*b, xmax=1.5*b)
    ax_b.legend()

    fig.tight_layout()

    plt.show()

    fig_, ax = plt.subplots()
    x = jnp.linspace(0.0, 5.0, num=300)
    ax.plot(x, log_normal_pdf(x), linewidth=1, color="black", label=r"$f(x)$")
    ax.plot(x, jax.vmap(kde)(x), linewidth=1, color="blue", label=r"$f_n(x)$")
    plt.show()