from typing import Callable

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from icecream import ic
from jax import random
from jaxtyping import Float, Array
from tqdm import tqdm
    
from HPD_estimation_as_multiscale_Z_estimator import HPD_estimator_with_confidence

HPD_LEVEL = 0.95
CI_LEVEL = 0.95
KERNEL = lambda x: 1/jnp.sqrt(2*jnp.pi)*jnp.exp(-x**2/2)
PHI = lambda x: jnp.log(x)
KERNEL_SQUARED_NORM = 1 / (2*jnp.sqrt(jnp.pi))



def estimate_hpd_accross_replicates(data: pl.DataFrame):
    rows = []
    for replicate in data.iter_rows():
        rows.append(
            {
                "quantity":
                "point_estimate": a_hat,
                "ci_lower": a_L,
                "ci_upper": a_U,
                "method": "multiscale_z_estimator"
                "target":
                "sample_min":
                "sample_max":
                "hpd_target_coverage": HPD_LEVEL,
                "ci_target_coverage": CI_LEVEL,
            }
        )
        rows.append(
            {
                "quantity":
                "point_estimate": b_hat,
                "ci_lower": b_L,
                "ci_upper": b_U,
                "method": "multiscale_z_estimator"
                "target":
                "sample_min":
                "sample_max":
                "hpd_target_coverage": HPD_LEVEL,
                "ci_target_coverage": CI_LEVEL,
            }
        )
    