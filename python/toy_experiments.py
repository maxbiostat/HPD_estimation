import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from jax import random
from jax.scipy.optimize import minimize
from jaxtyping import Float, Array
from quadax import quadgk
from scipy.integrate import quad
from scipy.optimize import root
from tqdm import tqdm
    
from python.kde_based_estimator import HPD_estimator

KERNEL = lambda x: 1/jnp.sqrt(2*jnp.pi)*jnp.exp(-x**2/2)
PHI = lambda x: jnp.log(x)
KERNEL_SQUARED_NORM = 1 / (2*jnp.sqrt(jnp.pi))
KERNEL_SECOND_MOMENT = 1.0
ALPHA = 0.95
BETA = 0.95
SIGMA = 1.0
MU = 0.0
GAMMA_ALPHA = 2
GAMMA_SCALE = 1
N_SAMPLES = int(1E3)
N_REPLICATIONS = int(5E2)

def log_normal_pdf(
    x: Float[Array, ""],
) -> Float[Array, ""]:
    return jnp.exp(-(jnp.log(x) - MU)**2 / (2 * SIGMA**2)) / (x * SIGMA * jnp.sqrt(2*jnp.pi))

def gamma_pdf(
    x: Float[Array, ""],
) -> Float[Array, ""]:
    return 1 / (jax.scipy.special.gamma(GAMMA_ALPHA) * jnp.pow(GAMMA_SCALE, GAMMA_ALPHA)) * jnp.pow(x, GAMMA_ALPHA - 1) * jnp.exp(- x / GAMMA_SCALE)

pdf = log_normal_pdf

true_hpd_sol = root(
    fun=lambda x: jnp.array([quad(pdf, x[0], x[1])[0] - ALPHA, pdf(x[0]) - pdf(x[1])]),
    x0=[0.02, 5],
    jac=lambda x: jnp.array([
            [-pdf(x[0]), pdf(x[1])],
            [jax.grad(pdf)(x[0]), -jax.grad(pdf)(x[1])]
    ]),
    method="hybr",
)
print(true_hpd_sol.message)
assert true_hpd_sol.success
true_hpd = true_hpd_sol.x
a, b = true_hpd[0], true_hpd[1]
print(a)
print(b)


if __name__ == "__main__":
    key = random.key(100)
    results = np.empty(shape=(N_REPLICATIONS, 6), dtype=np.float32)
    a_L_less_than_zero = 0
    for i in tqdm(range(N_REPLICATIONS)):
        key, _ = random.split(key)
        data = jnp.exp(MU + SIGMA * random.normal(key=key, shape=(N_SAMPLES,)))
        # data = GAMMA_SCALE * random.gamma(key=key, a=GAMMA_ALPHA, shape=(N_SAMPLES,))
        key, _ = random.split(key)
        # hpd_estimator = HPD_estimator_with_confidence(
        #     data, ALPHA, BETA, KERNEL, PHI, KERNEL_SQUARED_NORM, KERNEL_SECOND_MOMENT, key=key
        # )
        hpd_estimator = HPD_estimator(
            data, ALPHA, BETA, KERNEL, PHI, KERNEL_SQUARED_NORM, key=key, single_rate=True
        )
        # a_hat, a_L, a_U, b_hat, b_L, b_U = hpd_estimator
        a_hat = hpd_estimator.a_hat
        a_L = hpd_estimator.a_L
        a_U = hpd_estimator.a_U
        b_hat = hpd_estimator.b_hat
        b_L = hpd_estimator.b_L
        b_U = hpd_estimator.b_U
        a_L_less_than_zero += (a_L <= 0)
        results[i, :] = np.array([a_hat, a_L, a_U, b_hat, b_L, b_U])
    a_L = results[:, 1]
    a_U = results[:, 2]
    b_L = results[:, 4]
    b_U = results[:, 5]
    a_coverage, b_coverage = np.mean(np.logical_and(a_L < a, a < a_U), axis=0), jnp.mean(np.logical_and(b_L < b, b < b_U), axis=0)
    average_a_interval_length = np.mean(a_U - a_L, axis=0)
    average_b_interval_length = np.mean(b_U - b_L, axis=0)
    print(f"coverage for a: {100*a_coverage:.2f}%")
    print(f"coverage for b: {100*b_coverage:.2f}%")
    print(f"average length of interval for a: {average_a_interval_length:.2f}")
    print(f"average length of interval for b: {average_b_interval_length:.2f}")
    a_hat = results[:, 0]
    b_hat = results[:, 3]
    delta_a = a_hat - a
    delta_b = b_hat - b
    reference = 2 * delta_a * delta_b / (delta_a**2 + delta_b**2)
    ic(jnp.max(reference))
    plt.hist(reference, bins="auto")
    plt.show()
    print(f"% of a_L less than 0: {100*a_L_less_than_zero/N_REPLICATIONS:.2f}%")