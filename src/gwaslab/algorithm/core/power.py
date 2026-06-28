"""GWAS statistical power calculations.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple, Union

import numpy as np
import scipy.stats as stats


def or_to_rr(genotype_or: float, prevalence: float) -> float:
    """Convert odds ratio to genotype relative risk using population prevalence.

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
    genotype_or = float(genotype_or)
    prevalence = float(prevalence)
    return genotype_or / ((1.0 - prevalence) + (genotype_or * prevalence))


def resolve_genotype_rr(
    genotype_rr: Optional[float] = None,
    genotype_or: Optional[float] = None,
    beta: Optional[float] = None,
    prevalence: float = 0.15,
    use_or_to_rr: bool = False,
) -> float:
    """Resolve a binary-trait genotype relative risk from OR, RR, or log-OR beta.

Parameters
----------
genotype_rr : float, optional
    Genotype relative risk when already known.
genotype_or : float, optional
    Genotype odds ratio.
beta : float, optional
    Log odds ratio; exponentiated when OR/RR are missing.
prevalence : float, default 0.15
    Disease prevalence for Zhang–Kai OR→RR conversion.
use_or_to_rr : bool, default False
    When True, convert OR to RR via :func:`or_to_rr`.

Returns
-------
float
    Genotype relative risk used by :func:`power_binary`.
"""
    if genotype_rr is not None:
        return float(genotype_rr)
    if genotype_or is None:
        if beta is None:
            return 0.1
        genotype_or = float(np.exp(beta))
    if use_or_to_rr:
        return or_to_rr(float(genotype_or), float(prevalence))
    return float(genotype_or)


def genotype_disease_probabilities(
    genotype_rr: float,
    daf: float,
    prevalence: float,
) -> Tuple[float, float, float, float, float, float]:
    """Genotype counts and disease probabilities under an additive risk model.

Parameters
----------
genotype_rr : float
    Genotype relative risk for the effect allele.
daf : float
    Derived allele frequency.
prevalence : float
    Population disease prevalence.

Returns
-------
tuple of float
    ``(aaf, abf, bbf, p_aa, p_ab, p_bb)`` genotype frequency and disease
    probability terms for AA, AB, and BB genotypes.
"""
    genotype_rr = float(genotype_rr)
    daf = float(daf)
    prevalence = float(prevalence)
    aaf = daf ** 2
    abf = 2.0 * daf * (1.0 - daf)
    bbf = (1.0 - daf) ** 2
    x = (2.0 * genotype_rr - 1.0, genotype_rr, 1.0)
    denom = x[0] * aaf + x[1] * abf + x[2] * bbf
    aap = x[0] * prevalence / denom
    abp = x[1] * prevalence / denom
    bbp = x[2] * prevalence / denom
    return aaf, abf, bbf, aap, abp, bbp


def power_quantitative(
    beta: float,
    eaf: float,
    n: int,
    sig_level: float = 5e-8,
    phenotypic_variance: float = 1.0,
) -> float:
    """Compute GWAS power for a quantitative trait.

Parameters
----------
beta : float
    Effect size on the log/linear scale.
eaf : float
    Effect allele frequency.
n : int
    Sample size.
sig_level : float, default 5e-8
    Significance threshold.
phenotypic_variance : float, default 1.0
    Phenotypic variance.

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
    beta = float(beta)
    eaf = float(eaf)
    n = int(n)
    c = stats.chi2.isf(sig_level, df=1)
    ncp = n * 2.0 * eaf * (1.0 - eaf) * (beta ** 2) / float(phenotypic_variance)
    return float(1.0 - stats.ncx2.cdf(c, df=1, nc=ncp))


def power_binary(
    genotype_rr: float,
    daf: float,
    ncase: int,
    ncontrol: int,
    prevalence: float,
    sig_level: float = 5e-8,
) -> float:
    """Compute GWAS power for a binary trait under an additive model.

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
sig_level : float, default 5e-8
    Significance threshold.

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
    genotype_rr = float(genotype_rr)
    daf = float(daf)
    ncase = int(ncase)
    ncontrol = int(ncontrol)
    prevalence = float(prevalence)
    aaf = daf ** 2
    abf = 2.0 * daf * (1.0 - daf)
    bbf = (1.0 - daf) ** 2
    x = (2.0 * genotype_rr - 1.0, genotype_rr, 1.0)

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


def power_grid_quantitative(
    eaf_range: Tuple[float, float],
    beta_range: Tuple[float, float],
    n: int,
    sig_level: float = 5e-8,
    phenotypic_variance: float = 1.0,
    n_matrix: int = 500,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a quantitative-trait power grid over MAF and beta.

Parameters
----------
eaf_range : tuple of float
    ``(min, max)`` effect allele frequency window.
beta_range : tuple of float
    ``(min, max)`` beta window.
n : int
    Sample size.
sig_level : float, default 5e-8
    Significance threshold.
phenotypic_variance : float, default 1.0
    Phenotypic variance.
n_matrix : int, default 500
    Grid resolution per axis.

Returns
-------
grid : ndarray
    Power values with shape ``(n_matrix, n_matrix)``.
eaf_values : ndarray
    EAF coordinates (high to low, matching legacy trumpet plots).
beta_values : ndarray
    Beta coordinates (low to high).
"""
    eaf_values = np.linspace(eaf_range[1], eaf_range[0], n_matrix)
    beta_values = np.linspace(beta_range[0], beta_range[1], n_matrix)
    eaf_col = eaf_values[:, np.newaxis]
    beta_row = beta_values[np.newaxis, :]
    c = stats.chi2.isf(sig_level, df=1)
    ncp = (
        int(n)
        * 2.0
        * eaf_col
        * (1.0 - eaf_col)
        * (beta_row ** 2)
        / float(phenotypic_variance)
    )
    grid = 1.0 - stats.ncx2.cdf(c, df=1, nc=ncp)
    return np.asarray(grid, dtype=float), eaf_values, beta_values


def power_grid_binary(
    eaf_range: Tuple[float, float],
    beta_range: Tuple[float, float],
    ncase: int,
    ncontrol: int,
    prevalence: float,
    sig_level: float = 5e-8,
    n_matrix: int = 500,
    use_or_to_rr: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a binary-trait power grid over DAF and log-OR beta.

Parameters
----------
eaf_range : tuple of float
    ``(min, max)`` derived allele frequency window.
beta_range : tuple of float
    ``(min, max)`` beta window (log OR scale).
ncase, ncontrol : int
    Case and control counts.
prevalence : float
    Disease prevalence.
sig_level : float, default 5e-8
    Significance threshold.
n_matrix : int, default 500
    Grid resolution per axis.
use_or_to_rr : bool, default False
    Convert exp(beta) from OR to RR before power evaluation.

Returns
-------
grid : ndarray
    Power values with shape ``(n_matrix, n_matrix)``.
eaf_values : ndarray
    DAF coordinates (high to low).
beta_values : ndarray
    Beta coordinates (low to high).
"""
    eaf_values = np.linspace(eaf_range[1], eaf_range[0], n_matrix)
    beta_values = np.linspace(beta_range[0], beta_range[1], n_matrix)
    grid = np.zeros((n_matrix, n_matrix), dtype=float)
    prevalence = float(prevalence)
    c = stats.norm.isf(sig_level / 2.0)
    for i, daf in enumerate(eaf_values):
        daf = float(daf)
        aaf = daf ** 2
        abf = 2.0 * daf * (1.0 - daf)
        bbf = (1.0 - daf) ** 2
        genotype_or = np.exp(beta_values)
        if use_or_to_rr:
            genotype_rr = genotype_or / ((1.0 - prevalence) + (genotype_or * prevalence))
        else:
            genotype_rr = genotype_or
        x0 = 2.0 * genotype_rr - 1.0
        x1 = genotype_rr
        x2 = np.ones_like(genotype_rr)
        denom = x0 * aaf + x1 * abf + x2 * bbf
        aap = x0 * prevalence / denom
        abp = x1 * prevalence / denom
        pcase = (aap * aaf + abp * abf * 0.5) / prevalence
        pcontrol = ((1.0 - aap) * aaf + (1.0 - abp) * abf * 0.5) / (1.0 - prevalence)
        vcase = pcase * (1.0 - pcase)
        vcontrol = pcontrol * (1.0 - pcontrol)
        num = pcase - pcontrol
        for_sqrt = (vcase / ncase + vcontrol / ncontrol) * 0.5
        for_sqrt = np.asarray(for_sqrt, dtype=float)
        for_sqrt[for_sqrt < 0] = np.nan
        den = np.sqrt(for_sqrt)
        u = num / den
        grid[i, :] = 1.0 - stats.norm.cdf(c - u) + stats.norm.cdf(-c - u)
    return grid, eaf_values, beta_values


def power_contour_at_target(
    grid: np.ndarray,
    eaf_values: np.ndarray,
    beta_values: np.ndarray,
    target_power: float,
) -> list[Tuple[float, float]]:
    """Extract EAF–beta pairs on a power threshold contour.

Parameters
----------
grid : ndarray
    Power grid from :func:`power_grid_quantitative` or :func:`power_grid_binary`.
eaf_values : ndarray
    Row coordinates aligned with ``grid``.
beta_values : ndarray
    Column coordinates aligned with ``grid``.
target_power : float
    Target statistical power in ``(0, 1)``.

Returns
-------
list of tuple of float
    ``(eaf, beta)`` pairs with power at or above ``target_power``.
"""
    target_power = float(target_power)
    if target_power <= 0:
        return []
    n_matrix = grid.shape[0]
    contour: list[Tuple[float, float]] = []
    i, j = 1, 1
    while i < n_matrix - 1 and j < n_matrix - 1:
        if grid[i, j] < target_power:
            j += 1
        else:
            i += 1
            contour.append((float(eaf_values[i]), float(beta_values[j])))
    return contour


def infer_maf_range(maf_min: float) -> Tuple[float, float]:
    """Infer a default MAF axis window for trumpet plots from observed data.

Parameters
----------
maf_min : float
    Minimum observed minor allele frequency in the plot input.

Returns
-------
tuple of float
    ``(lower, upper)`` MAF limits with upper fixed at ``0.5``.
"""
    maf_min = float(maf_min)
    maf_min_power = np.floor(-np.log10(maf_min)) + 1
    lower = min(float(10.0 ** (-maf_min_power)), 1e-4)
    return lower, 0.5


def infer_beta_range(beta_max: float) -> Tuple[float, float]:
    """Infer a default beta axis window for trumpet plots.

Parameters
----------
beta_max : float
    Maximum absolute observed effect size.

Returns
-------
tuple of float
    ``(lower, upper)`` beta limits.
"""
    beta_max = float(beta_max)
    if beta_max > 3:
        return 0.0001, beta_max
    return 0.0001, 3.0
