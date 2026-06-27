"""Fixed-effects meta-analysis helpers."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.stats import norm


def fixed_effect_meta(
    beta: np.ndarray,
    se: np.ndarray,
) -> Tuple[float, float, float, float]:
    """
    Inverse-variance fixed-effects meta-analysis for one variant.

    Parameters
    ----------
    beta : numpy.ndarray
        Study-specific effect estimates.
    se : numpy.ndarray
        Study-specific standard errors (must be positive).

    Returns
    -------
    beta_fe : float
        Combined effect estimate.
    se_fe : float
        Standard error of the combined estimate.
    z_fe : float
        Z-score of the combined estimate.
    p_fe : float
        Two-sided P-value.
    """
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    weights = 1.0 / (se ** 2)
    w_sum = weights.sum()
    beta_fe = float((beta * weights).sum() / w_sum)
    se_fe = float(np.sqrt(1.0 / w_sum))
    z_fe = beta_fe / se_fe if se_fe > 0 else float("nan")
    p_fe = float(norm.sf(abs(z_fe)) * 2.0)
    return beta_fe, se_fe, z_fe, p_fe


def cochran_q_from_studies(
    beta1: np.ndarray,
    se1: np.ndarray,
    beta2: np.ndarray,
    se2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute Cochran Q, heterogeneity P, and I-squared for two-study comparisons.

    Parameters
    ----------
    beta1, se1 : numpy.ndarray
        Effect and SE for study 1.
    beta2, se2 : numpy.ndarray
        Effect and SE for study 2.

    Returns
    -------
    q : numpy.ndarray
        Cochran Q statistic (df=1 per row).
    het_p : numpy.ndarray
        Heterogeneity P-values.
    i2 : numpy.ndarray
        I-squared statistics in [0, 1].

    References
    ----------
    Cochran, W. G. (1954). The combination of estimates from different
    experiments. Biometrics, 10(1), 101-129.

    Huedo-Medina, T. B., Sánchez-Meca, J., Marín-Martínez, F., & Botella, J.
    (2006). Assessing heterogeneity in meta-analysis: Q statistic or I² index?
    Psychological Methods, 11(2), 193-206.
    """
    beta1 = np.asarray(beta1, dtype=float)
    se1 = np.asarray(se1, dtype=float)
    beta2 = np.asarray(beta2, dtype=float)
    se2 = np.asarray(se2, dtype=float)

    w1 = 1.0 / (se1 ** 2)
    w2 = 1.0 / (se2 ** 2)
    beta_fe = (w1 * beta1 + w2 * beta2) / (w1 + w2)
    q = w1 * (beta1 - beta_fe) ** 2 + w2 * (beta2 - beta_fe) ** 2

    from scipy import stats as sp_stats

    het_p = sp_stats.chi2.sf(q, 1)
    i2 = np.where(q <= 0, 0.0, (q - 1.0) / q)
    i2 = np.clip(i2, 0.0, None)
    return q, het_p, i2
