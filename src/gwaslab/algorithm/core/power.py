"""GWAS statistical power calculations."""

from __future__ import annotations

import numpy as np
import scipy.stats as stats


def or_to_rr(genotype_or: float, prevalence: float) -> float:
    """
    Convert odds ratio to genotype relative risk using population prevalence.

    Parameters
    ----------
    genotype_or : float
        Genotype odds ratio.
    prevalence : float
        Disease prevalence in the population.

    Returns
    -------
    float
        Genotype relative risk.

    References
    ----------
    Zhang, J., & Kai, F. Y. (1998). What's the relative risk? A method of
    correcting the odds ratio in cohort studies of common outcomes. JAMA,
    280(19), 1690-1691.
    """
    return genotype_or / ((1.0 - prevalence) + (genotype_or * prevalence))


def power_quantitative(
    beta: float,
    eaf: float,
    n: int,
    sig_level: float = 5e-8,
    phenotypic_variance: float = 1.0,
) -> float:
    """
    Compute GWAS power for a quantitative trait.

    Parameters
    ----------
    beta : float
        Effect size on the log/linear scale.
    eaf : float
        Effect allele frequency.
    n : int
        Sample size.
    sig_level : float, optional
        Significance threshold. Default is 5e-8.
    phenotypic_variance : float, optional
        Phenotypic variance. Default is 1.0.

    Returns
    -------
    float
        Statistical power in [0, 1].

    References
    ----------
    Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). Joint
    analysis is more efficient than replication-based analysis for two-stage
    genome-wide association studies. Nature Genetics, 38(2), 209-213.
    """
    c = stats.chi2.isf(sig_level, df=1)
    ncp = n * 2.0 * eaf * (1.0 - eaf) * (beta ** 2) / phenotypic_variance
    return float(1.0 - stats.ncx2.cdf(c, df=1, nc=ncp))


def power_binary(
    genotype_rr: float,
    daf: float,
    ncase: int,
    ncontrol: int,
    prevalence: float,
    sig_level: float = 5e-8,
) -> float:
    """
    Compute GWAS power for a binary trait under an additive model.

    Parameters
    ----------
    genotype_rr : float
        Genotype relative risk for the effect allele.
    daf : float
        Derived allele frequency.
    ncase : int
        Number of cases.
    ncontrol : int
        Number of controls.
    prevalence : float
        Disease prevalence.
    sig_level : float, optional
        Significance threshold. Default is 5e-8.

    Returns
    -------
    float
        Statistical power in [0, 1].

    References
    ----------
    Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). Joint
    analysis is more efficient than replication-based analysis for two-stage
    genome-wide association studies. Nature Genetics, 38(2), 209-213.
    """
    aaf = daf ** 2
    abf = 2.0 * daf * (1.0 - daf)
    bbf = (1.0 - daf) ** 2
    x = [2.0 * genotype_rr - 1.0, genotype_rr, 1.0]

    denom = x[0] * aaf + x[1] * abf + x[2] * bbf
    aap = x[0] * prevalence / denom
    abp = x[1] * prevalence / denom

    pcase = (aap * aaf + abp * abf * 0.5) / prevalence
    pcontrol = ((1.0 - aap) * aaf + (1.0 - abp) * abf * 0.5) / (1.0 - prevalence)

    vcase = pcase * (1.0 - pcase)
    vcontrol = pcontrol * (1.0 - pcontrol)
    num = pcase - pcontrol
    for_sqrt = (vcase / ncase + vcontrol / ncontrol) * 0.5
    if np.iterable(for_sqrt):
        for_sqrt = np.asarray(for_sqrt, dtype=float)
        for_sqrt[for_sqrt < 0] = np.nan
    den = np.sqrt(for_sqrt)
    u = num / den

    c = stats.norm.isf(sig_level / 2.0)
    return float(1.0 - stats.norm.cdf(c - u) + stats.norm.cdf(-c - u))
