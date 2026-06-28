"""Genomic inflation factor (lambda GC) calculations.
"""

from __future__ import annotations

import numpy as np
import scipy.stats as stats

DEFAULT_QUANTILE = 0.5


def _expected_chi2(quantile: float) -> float:
    return float(stats.chi2.ppf(quantile, df=1))


def lambda_gc_from_p(
    p: np.ndarray,
    quantile: float = DEFAULT_QUANTILE,
) -> float:
    """Compute genomic inflation factor lambda from P-values.

    Transforms P-values to chi-squared statistics and compares the observed
    quantile to the expected quantile under the null.

Parameters
----------
p : numpy.ndarray
    Two-sided P-values in (0, 1]. NaNs are ignored.
quantile : float, default 0.5
    Quantile for the chi-squared transform (median).
Returns
-------
float
    Genomic inflation factor (lambda GC), or ``nan`` if no valid values.

Notes
-----
    Orchestration (CHR filtering, column auto-detect) lives in
    ``util.util_in_calculate_gc._lambda_GC``.

References
----------
    Devlin, B., & Roeder, K. (1999). Genomic control for association studies.
    Biometrics, 55(4), 964-975.
"""
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return float("nan")
    observed = float(stats.chi2.isf(np.nanmedian(p), df=1))
    expected = _expected_chi2(quantile)
    return observed / expected


def lambda_gc_from_mlog10p(
    mlog10p: np.ndarray,
    quantile: float = DEFAULT_QUANTILE,
) -> float:
    """Compute lambda GC from -log10(P) values.

Parameters
----------
mlog10p : numpy.ndarray
    Minus log10 P-values. NaNs are ignored.
quantile : float, default 0.5
    Quantile for comparison (median).
Returns
-------
float
    Genomic inflation factor, or ``nan`` if no valid values.

References
----------
    Devlin, B., & Roeder, K. (1999). Genomic control for association studies.
    Biometrics, 55(4), 964-975.
"""
    mlog10p = np.asarray(mlog10p, dtype=float)
    mlog10p = mlog10p[np.isfinite(mlog10p)]
    if mlog10p.size == 0:
        return float("nan")
    p = np.power(10.0, -mlog10p)
    return lambda_gc_from_p(p, quantile=quantile)


def lambda_gc_from_z(
    z: np.ndarray,
    quantile: float = DEFAULT_QUANTILE,
) -> float:
    """Compute lambda GC from Z-scores.

Parameters
----------
z : numpy.ndarray
    Z-scores. NaNs are ignored.
quantile : float, default 0.5
    Quantile for comparison (median).
Returns
-------
float
    Genomic inflation factor, or ``nan`` if no valid values.

References
----------
    Devlin, B., & Roeder, K. (1999). Genomic control for association studies.
    Biometrics, 55(4), 964-975.
"""
    z = np.asarray(z, dtype=float)
    z = z[np.isfinite(z)]
    if z.size == 0:
        return float("nan")
    observed = float(np.median(z ** 2))
    expected = _expected_chi2(quantile)
    return observed / expected


def lambda_gc_from_chisq(
    chisq: np.ndarray,
    quantile: float = DEFAULT_QUANTILE,
) -> float:
    """Compute lambda GC from chi-squared statistics.

Parameters
----------
chisq : numpy.ndarray
    Chi-squared statistics (df=1). NaNs are ignored.
quantile : float, default 0.5
    Quantile for comparison (median).
Returns
-------
float
    Genomic inflation factor, or ``nan`` if no valid values.

References
----------
    Devlin, B., & Roeder, K. (1999). Genomic control for association studies.
    Biometrics, 55(4), 964-975.
"""
    chisq = np.asarray(chisq, dtype=float)
    chisq = chisq[np.isfinite(chisq)]
    if chisq.size == 0:
        return float("nan")
    observed = float(np.median(chisq))
    expected = _expected_chi2(quantile)
    return observed / expected
