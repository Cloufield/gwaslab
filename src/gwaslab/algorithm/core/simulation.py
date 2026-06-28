"""Simulation and z-score distribution helpers.
"""

from __future__ import annotations

from math import sqrt

import numpy as np
import scipy.stats as stats
from scipy.special import erfc


def norm_sf(x: np.ndarray) -> np.ndarray:
    """Compute the standard normal survival function P(Z > x).

Parameters
----------
x : numpy.ndarray
    Z-scores.
Returns
-------
numpy.ndarray
    Survival probabilities.

Notes
-----
    Uses ``erfc`` for vectorized evaluation. Orchestration lives in
    ``util.util_in_simulate``.
"""
    x = np.asarray(x, dtype=np.float64)
    x_clipped = np.clip(x, -37.0, 37.0)
    result = 0.5 * erfc(x_clipped / sqrt(2.0))
    return np.where(np.isfinite(result), result, 0.0)


def p_from_z(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to two-sided P-values.

Parameters
----------
z : numpy.ndarray
    Z-scores.
Returns
-------
numpy.ndarray
    Two-sided P-values.
"""
    return 2.0 * norm_sf(np.abs(np.asarray(z, dtype=np.float64)))


def z_to_mlog10p(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to -log10(P) using log-space for precision.

Parameters
----------
z : numpy.ndarray
    Z-scores.
Returns
-------
numpy.ndarray
    Minus log10 two-sided P-values.

Notes
-----
    Equivalent to ``core.conversions.z_to_mlog10p``; kept here for the
    simulation module import path.
"""
    z_arr = np.asarray(z, dtype=np.float64)
    log_pvalue = np.log(2.0) + stats.norm.logsf(np.abs(z_arr))
    return -log_pvalue / np.log(10.0)
