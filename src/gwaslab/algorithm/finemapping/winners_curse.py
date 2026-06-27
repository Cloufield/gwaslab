"""Winner's curse bias correction."""

from __future__ import annotations

from typing import Union

import numpy as np
import scipy.optimize as optimize
import scipy.stats as stats

DEFAULT_SIG_LEVEL = 5e-8


def _selection_threshold(sig_level: float) -> float:
    return float(np.sqrt(stats.chi2.ppf(1.0 - sig_level, df=1)))


def _bias_equation(beta_true: float, beta_obs: float, se: float, c: float) -> float:
    z = beta_true / se
    numerator = stats.norm.pdf(z - c) - stats.norm.pdf(-z - c)
    denominator = stats.norm.cdf(z - c) + stats.norm.cdf(-z - c)
    return beta_true + se * numerator / denominator - beta_obs


def winners_curse_correct(
    beta: Union[float, np.ndarray],
    se: Union[float, np.ndarray],
    sig_level: float = DEFAULT_SIG_LEVEL,
) -> Union[float, np.ndarray]:
    """
    Apply winner's curse bias correction to observed effect sizes.

    Solves the Zhong-Prentice bias equation using Brent's method on each
    variant.

    Parameters
    ----------
    beta : float or numpy.ndarray
        Observed effect estimates.
    se : float or numpy.ndarray
        Standard errors of ``beta``.
    sig_level : float, optional
        GWAS significance threshold. Default is 5e-8.

    Returns
    -------
    float or numpy.ndarray
        Bias-corrected effect estimates.

    Notes
    -----
    Orchestration and logging live in ``util.util_in_correct_winnerscurse``.

    References
    ----------
    Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and
    confidence intervals for odds ratios in genome-wide association studies.
    Biostatistics, 9(4), 621-634.
    """
    c = _selection_threshold(sig_level)

    def _solve_scalar(b_obs: float, s: float) -> float:
        return float(
            optimize.brentq(
                lambda x: _bias_equation(x, b_obs, s, c),
                a=-100.0,
                b=100.0,
                maxiter=1000,
            )
        )

    if np.ndim(beta) == 0 and np.ndim(se) == 0:
        return _solve_scalar(float(beta), float(se))

    beta_arr = np.asarray(beta, dtype=float)
    se_arr = np.asarray(se, dtype=float)
    out = np.empty_like(beta_arr, dtype=float)
    for idx in np.ndindex(beta_arr.shape):
        out[idx] = _solve_scalar(beta_arr[idx], se_arr[idx])
    return out
