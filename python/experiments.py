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
    
from HPD_estimation_as_Z_estimator_modified import HPD_estimator_with_confidence

KERNEL = lambda x: 1/jnp.sqrt(2*jnp.pi)*jnp.exp(-x**2/2)
PHI = lambda x: jnp.log(x)
KERNEL_SQUARED_NORM = 1 / (2*jnp.sqrt(jnp.pi))
ALPHA = 0.90
BETA = 0.95
SIGMA = 1
MU = 0.0
N_SAMPLES = int(1E4)
N_REPLICATIONS = int(1E2)

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
print(true_hpd_sol.message)
assert true_hpd_sol.success
true_hpd = true_hpd_sol.x
a, b = true_hpd[0], true_hpd[1]
print(a)
print(b)


if __name__ == "__main__":
    key = random.key(0)
    coverage = np.empty(shape=(N_REPLICATIONS, 2), dtype=jnp.int32)
    for i in tqdm(range(N_REPLICATIONS)):
        key, _ = random.split(key)
        data = jnp.exp(MU + SIGMA * random.normal(key=key, shape=(N_SAMPLES,)))
        hpd_estimator = HPD_estimator_with_confidence(
            data, ALPHA, BETA, KERNEL, PHI, KERNEL_SQUARED_NORM,
        )
        a_hat, a_L, a_U, b_hat, b_L, b_U = hpd_estimator
        coverage[i, 0] = a_L < a < a_U
        coverage[i, 1] = b_L < b < b_U
    a_coverage, b_coverage = jnp.mean(coverage, axis=0)
    print(f"coverage for a: {100*a_coverage:.2f}%")
    print(f"coverage for b: {100*b_coverage:.2f}%")