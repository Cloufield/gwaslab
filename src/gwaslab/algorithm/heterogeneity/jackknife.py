"""Jackknife resampling statistics.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def jackknife_correlation_se(x: np.ndarray, y: np.ndarray) -> float:
    """Jackknife standard error of Pearson correlation.

Parameters
----------
x : numpy.ndarray
    First variable.
y : numpy.ndarray
    Second variable, same length as ``x``.
Returns
-------
float
    Jackknife SE of the correlation coefficient.

Notes
-----
    Implementation note: leave-one-out Pearson r with standard jackknife SE.
    Orchestration (pandas columns, logging) lives in ``viz.viz_plot_*``.

References
----------
    Miller, R. G. (1974). The jackknife—a review. Biometrics, 30(1), 1-15.
"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n < 3:
        return float("nan")

    rs = np.empty(n, dtype=float)
    for i in range(n):
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        rs[i] = stats.pearsonr(x[keep], y[keep]).statistic

    r_se = np.sqrt((n - 1) / n * np.sum((rs - np.mean(rs)) ** 2))
    return float(r_se)
