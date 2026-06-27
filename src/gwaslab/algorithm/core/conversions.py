"""GWAS summary-statistic conversions (BETA, SE, Z, P, OR, MAF)."""

from __future__ import annotations

import numpy as np
import scipy.stats as stats
from scipy.special import erfcinv
from scipy.stats import norm

Z_CRITICAL_95 = float(norm.ppf(0.975))


def betase_to_z(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to two-sided Z-scores."""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        return beta / se


def betase_to_p(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to two-sided P-values."""
    z = betase_to_z(beta, se)
    return 2.0 * norm.sf(np.abs(z))


def z_to_p(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to two-sided P-values."""
    z = np.asarray(z, dtype=np.float64)
    return 2.0 * norm.sf(np.abs(z))


def z_to_mlog10p(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to -log10(P)."""
    z = np.asarray(z, dtype=np.float64)
    log_pvalue = np.log(2.0) + norm.logsf(np.abs(z))
    return -log_pvalue / np.log(10.0)


def betase_to_mlog10p(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to -log10(P)."""
    return z_to_mlog10p(betase_to_z(beta, se))


def p_to_chisq(p: np.ndarray) -> np.ndarray:
    """Convert P-values to chi-squared statistics (df=1)."""
    return stats.chi2.isf(np.asarray(p, dtype=np.float64), 1)


def z_to_chisq(z: np.ndarray) -> np.ndarray:
    """Convert Z-scores to chi-squared statistics (df=1)."""
    z = np.asarray(z, dtype=np.float64)
    return z ** 2


def chisq_to_p(chisq: np.ndarray) -> np.ndarray:
    """Convert chi-squared statistics to P-values (df=1)."""
    return stats.chi2.sf(np.asarray(chisq, dtype=np.float64), 1)


def mlog10p_to_p(mlog10p: np.ndarray) -> np.ndarray:
    """Convert -log10(P) to P-values."""
    return np.power(10.0, -np.asarray(mlog10p, dtype=np.float64))


def p_to_mlog10p(p: np.ndarray) -> np.ndarray:
    """Convert P-values to -log10(P)."""
    return -np.log10(np.asarray(p, dtype=np.float64))


def or_to_beta(odds_ratio: np.ndarray) -> np.ndarray:
    """Convert odds ratios to log-odds (BETA)."""
    return np.log(np.asarray(odds_ratio, dtype=np.float64))


def beta_to_or(beta: np.ndarray) -> np.ndarray:
    """Convert log-odds (BETA) to odds ratios."""
    return np.exp(np.asarray(beta, dtype=np.float64))


def betase_to_or_95l(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to lower bound of 95% OR CI."""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    return np.exp(beta - Z_CRITICAL_95 * se)


def betase_to_or_95u(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Convert BETA and SE to upper bound of 95% OR CI."""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    return np.exp(beta + Z_CRITICAL_95 * se)


def betap_to_se(beta: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Convert BETA and P to SE."""
    beta = np.asarray(beta, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    abs_z = np.sqrt(2.0) * erfcinv(p)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(beta / abs_z)


def or_or95u_to_se(odds_ratio: np.ndarray, or_95u: np.ndarray) -> np.ndarray:
    """Convert OR and OR upper CI bound to SE."""
    odds_ratio = np.asarray(odds_ratio, dtype=np.float64)
    or_95u = np.asarray(or_95u, dtype=np.float64)
    return (np.log(or_95u) - np.log(odds_ratio)) / Z_CRITICAL_95


def or_or95l_to_se(odds_ratio: np.ndarray, or_95l: np.ndarray) -> np.ndarray:
    """Convert OR and OR lower CI bound to SE."""
    odds_ratio = np.asarray(odds_ratio, dtype=np.float64)
    or_95l = np.asarray(or_95l, dtype=np.float64)
    return (np.log(odds_ratio) - np.log(or_95l)) / Z_CRITICAL_95


def eaf_to_maf(eaf: np.ndarray) -> np.ndarray:
    """Convert effect allele frequency to minor allele frequency."""
    eaf = np.asarray(eaf, dtype=np.float64)
    return np.minimum(eaf, 1.0 - eaf)
