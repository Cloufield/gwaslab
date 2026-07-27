"""GWAS summary-statistic conversions (BETA, SE, Z, P, OR, MAF).
"""

from __future__ import annotations

import re

import numpy as np
import scipy.stats as stats
from scipy.special import erfcinv
from scipy.stats import norm

Z_CRITICAL_95 = float(norm.ppf(0.975))


def betase_to_z(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to two-sided Z-scores.

Parameters
----------
beta : numpy.ndarray
    Effect size estimates.
se : numpy.ndarray
    Standard errors.
Returns
-------
numpy.ndarray
    Two-sided Z-scores.
"""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        return beta / se


def betase_to_p(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to two-sided P-values.

Parameters
----------
beta : numpy.ndarray
    Effect size estimates.
se : numpy.ndarray
    Standard errors.
Returns
-------
numpy.ndarray
    Two-sided P-values.
"""
    z = betase_to_z(beta, se)
    return 2.0 * norm.sf(np.abs(z))


def z_to_p(z: np.ndarray) -> np.ndarray:
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
    z = np.asarray(z, dtype=np.float64)
    return 2.0 * norm.sf(np.abs(z))


def z_to_mlog10p(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to -log10(P).

Parameters
----------
z : numpy.ndarray
    Z-scores.
Returns
-------
numpy.ndarray
    Minus log10 two-sided P-values.
"""
    z = np.asarray(z, dtype=np.float64)
    log_pvalue = np.log(2.0) + norm.logsf(np.abs(z))
    return -log_pvalue / np.log(10.0)


def betase_to_mlog10p(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to -log10(P).

Parameters
----------
beta : numpy.ndarray
    Effect size estimates.
se : numpy.ndarray
    Standard errors.
Returns
-------
numpy.ndarray
    Minus log10 two-sided P-values.
"""
    return z_to_mlog10p(betase_to_z(beta, se))


def p_to_chisq(p: np.ndarray) -> np.ndarray:
    """Convert P-values to chi-squared statistics (df=1).

Parameters
----------
p : numpy.ndarray
    Two-sided P-values.
Returns
-------
numpy.ndarray
    Chi-squared statistics with one degree of freedom.
"""
    return stats.chi2.isf(np.asarray(p, dtype=np.float64), 1)


def z_to_chisq(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to chi-squared statistics (df=1).

Parameters
----------
z : numpy.ndarray
    Z-scores.
Returns
-------
numpy.ndarray
    Chi-squared statistics with one degree of freedom.
"""
    z = np.asarray(z, dtype=np.float64)
    return z ** 2


def chisq_to_p(chisq: np.ndarray) -> np.ndarray:
    """Convert chi-squared statistics to P-values (df=1).

Parameters
----------
chisq : numpy.ndarray
    Chi-squared statistics.
Returns
-------
numpy.ndarray
    Two-sided P-values.
"""
    return stats.chi2.sf(np.asarray(chisq, dtype=np.float64), 1)


def log10_p_to_mantissa_exponent(log10_p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split log10(P) into mantissa and exponent (P = mantissa * 10**exponent).

Parameters
----------
log10_p : numpy.ndarray
    Base-10 logarithm of P-values.
Returns
-------
tuple[numpy.ndarray, numpy.ndarray]
    Mantissa in [1, 10) and integer-scaled exponent (floor of log10_p).
"""
    log10_p = np.asarray(log10_p, dtype=np.float64)
    exponent = np.floor(log10_p)
    mantissa = np.power(10.0, log10_p - exponent)
    return mantissa, exponent


def _parse_e_precision(fmt: str) -> int:
    """Extract fractional precision from a float format string like ``'{:.4e}'``."""
    match = re.search(r"\.(\d+)e", fmt, re.IGNORECASE)
    return int(match.group(1)) if match else 4


def mantissa_exponent_to_p_format_str(
    mantissa: float,
    exponent: float,
    fmt: str = "{:.4e}",
) -> str:
    """Format P as a scientific-notation string from mantissa and exponent.

    ``P = mantissa * 10**exponent`` with mantissa in [1, 10).

Parameters
----------
mantissa : float
    Significand in [1, 10).
exponent : float
    Base-10 exponent (integer-valued in practice).
fmt : str, optional
    Float format string used to infer fractional precision (e.g. ``'{:.4e}'``).
Returns
-------
str
    P formatted as ``mantissa e±exp`` (e.g. ``1.2345e-300``).
    """
    prec = _parse_e_precision(fmt)
    return f"{float(mantissa):.{prec}f}e{int(exponent):+03d}"


def mlog10p_to_p_format_str(mlog10p: float, fmt: str = "{:.4e}") -> str:
    """Format P as scientific notation from -log10(P) without float64 underflow.

Parameters
----------
mlog10p : float
    Minus log10 P-value.
fmt : str, optional
    Float format string used to infer fractional precision (e.g. ``'{:.4e}'``).
Returns
-------
str
    P formatted as ``mantissa e±exp`` without float64 underflow.
    """
    mantissa, exponent = mlog10p_to_mantissa_exponent(
        np.atleast_1d(np.asarray(mlog10p, dtype=np.float64))
    )
    return mantissa_exponent_to_p_format_str(
        float(np.ravel(mantissa)[0]), float(np.ravel(exponent)[0]), fmt
    )


def mlog10p_to_mantissa_exponent(mlog10p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split -log10(P) into mantissa and exponent without float64 underflow.

Parameters
----------
mlog10p : numpy.ndarray
    Minus log10 P-values.
Returns
-------
tuple[numpy.ndarray, numpy.ndarray]
    Mantissa and exponent such that P = mantissa * 10**exponent.
"""
    return log10_p_to_mantissa_exponent(-np.asarray(mlog10p, dtype=np.float64))


def mlog10p_to_p(mlog10p: np.ndarray) -> np.ndarray:
    """Convert -log10(P) to P-values.

Parameters
----------
mlog10p : numpy.ndarray
    Minus log10 P-values.
Returns
-------
numpy.ndarray
    P-values.
"""
    return np.power(10.0, -np.asarray(mlog10p, dtype=np.float64))


def p_to_mlog10p(p: np.ndarray) -> np.ndarray:
    """Convert P-values to -log10(P).

Parameters
----------
p : numpy.ndarray
    P-values.
Returns
-------
numpy.ndarray
    Minus log10 P-values.
"""
    return -np.log10(np.asarray(p, dtype=np.float64))


def or_to_beta(odds_ratio: np.ndarray) -> np.ndarray:
    """Convert odds ratios to log-odds (BETA).

Parameters
----------
odds_ratio : numpy.ndarray
    Odds ratios.
Returns
-------
numpy.ndarray
    Log-odds effect sizes.
"""
    return np.log(np.asarray(odds_ratio, dtype=np.float64))


def beta_to_or(beta: np.ndarray) -> np.ndarray:
    """Convert log-odds (BETA) to odds ratios.

Parameters
----------
beta : numpy.ndarray
    Log-odds effect sizes.
Returns
-------
numpy.ndarray
    Odds ratios.
"""
    return np.exp(np.asarray(beta, dtype=np.float64))


def betase_to_or_95l(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to lower bound of 95% OR CI.

Parameters
----------
beta : numpy.ndarray
    Log-odds effect sizes.
se : numpy.ndarray
    Standard errors.
Returns
-------
numpy.ndarray
    Lower bounds of 95% odds-ratio confidence intervals.
"""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    return np.exp(beta - Z_CRITICAL_95 * se)


def betase_to_or_95u(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to upper bound of 95% OR CI.

Parameters
----------
beta : numpy.ndarray
    Log-odds effect sizes.
se : numpy.ndarray
    Standard errors.
Returns
-------
numpy.ndarray
    Upper bounds of 95% odds-ratio confidence intervals.
"""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    return np.exp(beta + Z_CRITICAL_95 * se)


def betap_to_se(beta: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Convert BETA and P to SE.

Parameters
----------
beta : numpy.ndarray
    Effect size estimates.
p : numpy.ndarray
    Two-sided P-values.
Returns
-------
numpy.ndarray
    Standard errors.
"""
    beta = np.asarray(beta, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    abs_z = np.sqrt(2.0) * erfcinv(p)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(beta / abs_z)


def or_or95u_to_se(odds_ratio: np.ndarray, or_95u: np.ndarray) -> np.ndarray:
    """Convert OR and OR upper CI bound to SE.

Parameters
----------
odds_ratio : numpy.ndarray
    Odds ratios.
or_95u : numpy.ndarray
    Upper bounds of 95% odds-ratio confidence intervals.
Returns
-------
numpy.ndarray
    Standard errors on the log-odds scale.
"""
    odds_ratio = np.asarray(odds_ratio, dtype=np.float64)
    or_95u = np.asarray(or_95u, dtype=np.float64)
    return (np.log(or_95u) - np.log(odds_ratio)) / Z_CRITICAL_95


def or_or95l_to_se(odds_ratio: np.ndarray, or_95l: np.ndarray) -> np.ndarray:
    """Convert OR and OR lower CI bound to SE.

Parameters
----------
odds_ratio : numpy.ndarray
    Odds ratios.
or_95l : numpy.ndarray
    Lower bounds of 95% odds-ratio confidence intervals.
Returns
-------
numpy.ndarray
    Standard errors on the log-odds scale.
"""
    odds_ratio = np.asarray(odds_ratio, dtype=np.float64)
    or_95l = np.asarray(or_95l, dtype=np.float64)
    return (np.log(odds_ratio) - np.log(or_95l)) / Z_CRITICAL_95


def eaf_to_maf(eaf: np.ndarray) -> np.ndarray:
    """Convert effect allele frequency to minor allele frequency.

Parameters
----------
eaf : numpy.ndarray
    Effect allele frequencies.
Returns
-------
numpy.ndarray
    Minor allele frequencies.
"""
    eaf = np.asarray(eaf, dtype=np.float64)
    return np.minimum(eaf, 1.0 - eaf)
