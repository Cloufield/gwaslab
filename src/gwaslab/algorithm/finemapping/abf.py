"""Approximate Bayesian factor fine-mapping."""

from __future__ import annotations

import numpy as np


def log_abf_from_summary(
    beta: np.ndarray,
    se: np.ndarray,
    w: float = 0.2,
) -> np.ndarray:
    """
    Compute log approximate Bayes factors from summary statistics.

    Parameters
    ----------
    beta : numpy.ndarray
        Effect size estimates.
    se : numpy.ndarray
        Standard errors.
    w : float, optional
        Prior standard deviation on effect size. Default is 0.2 (binary traits).

    Returns
    -------
    numpy.ndarray
        Log ABF values per variant.

    Notes
    -----
    Typical priors: ``w=0.2`` binary, ``w=0.15`` quantitative. Column attachment
    and region extraction live in ``util.util_abf_finemapping``.

    References
    ----------
    Wakefield, J. (2007). A Bayesian measure of the probability of false
    discovery in genetic epidemiology studies. American Journal of Human
    Genetics, 81(2), 208-227.
    """
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    omega = w ** 2
    v = se ** 2
    r = omega / (omega + v)
    z = beta / se
    return 0.5 * (np.log(1.0 - r) + (r * z ** 2))


def pip_from_log_abf(log_abf: np.ndarray) -> np.ndarray:
    """
    Compute posterior inclusion probabilities from log ABF values.

    Parameters
    ----------
    log_abf : numpy.ndarray
        Log approximate Bayes factors for variants in a region.

    Returns
    -------
    numpy.ndarray
        PIP values summing to 1 over finite inputs.

    References
    ----------
    Wakefield, J. (2007). A Bayesian measure of the probability of false
    discovery in genetic epidemiology studies. American Journal of Human
    Genetics, 81(2), 208-227.
    """
    log_abf = np.asarray(log_abf, dtype=float)
    max_log = np.nanmax(log_abf)
    log_total = np.log(np.nansum(np.exp(log_abf - max_log))) + max_log
    log_pip = log_abf - log_total
    return np.exp(log_pip)
