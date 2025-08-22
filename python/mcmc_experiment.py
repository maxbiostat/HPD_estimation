import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
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
N_REPLICATIONS = int(5E2)

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
    data = pl.read_csv("python/for_Caio.csv")
    data = jnp.array(data["x"])
    print(data.min(), data.max())
    hpd_estimator = HPD_estimator_with_confidence(
        data, ALPHA, BETA, KERNEL, PHI, KERNEL_SQUARED_NORM,
    )
    a_hat, a_L, a_U, b_hat, b_L, b_U = hpd_estimator

    ic(a_hat, b_hat)
    print()

    ic(a_L)
    ic(a)
    ic(a_U)
    ic(a_L <= a <= a_U)

    ic(b_L)
    ic(b)
    ic(b_U)
    ic(b_L <= b <= b_U)
