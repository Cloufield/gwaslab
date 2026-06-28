"""Heritability scale conversion and per-SNP R-squared.
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

import numpy as np
from scipy.stats import norm


def h2_obs_to_liab(
    h2_obs: float,
    p_sample: float,
    k_pop: float,
    se_obs: Optional[float] = None,
) -> Union[float, Tuple[float, float]]:
    """Convert observed-scale heritability to liability-scale heritability.

Parameters
----------
h2_obs : float
    Heritability on the observed scale in an ascertained sample.
p_sample : float
    Phenotype prevalence in the sample, in (0, 1).
k_pop : float
    Phenotype prevalence in the population, in (0, 1).
se_obs : float, optional
    Standard error of ``h2_obs``. When provided, returns liability-scale SE.
Returns
-------
float or tuple of float
    Liability-scale heritability, or ``(h2_liab, se_liab)`` when ``se_obs`` is set.

Notes
-----
    Adapted from LDSC. Orchestration lives in ``util.util_in_convert_h2``.

References
----------
    Lee, S. H., Wray, N. R., Goddard, M. E., & Visscher, P. M. (2011).
    Estimating missing heritability for disease from genome-wide association
    studies. American Journal of Human Genetics, 89(6), 294-302.
"""
    if np.isnan(p_sample) and np.isnan(k_pop):
        return h2_obs
    if k_pop <= 0 or k_pop >= 1:
        raise ValueError("k_pop must be in the range (0, 1)")
    if p_sample <= 0 or p_sample >= 1:
        raise ValueError("p_sample must be in the range (0, 1)")

    t = norm.isf(k_pop)
    z = norm.pdf(t)
    conversion_factor = (k_pop ** 2 * (1 - k_pop) ** 2) / (p_sample * (1 - p_sample) * z ** 2)

    if se_obs is not None:
        return h2_obs * conversion_factor, conversion_factor * se_obs
    return h2_obs * conversion_factor


def h2_se_to_p(h2: float, se: float) -> float:
    """Convert heritability estimate and SE to a one-sided P-value.

Parameters
----------
h2 : float
    Heritability estimate.
se : float
    Standard error of ``h2``.
Returns
-------
float
    One-sided P-value under a standard normal null.
"""
    z = h2 / se
    return float(norm.sf(abs(z)))


def per_snp_r2_quantitative(
    beta: np.ndarray,
    eaf: np.ndarray,
    phenotypic_variance: float,
) -> np.ndarray:
    """Per-SNP R-squared for quantitative traits.

Parameters
----------
beta : numpy.ndarray
    Effect sizes.
eaf : numpy.ndarray
    Effect allele frequencies.
phenotypic_variance : float
    Total phenotypic variance Var(Y).
Returns
-------
numpy.ndarray
    Per-variant R-squared values.

References
----------
    Shim, H., et al. (2015). A multivariate genome-wide association analysis of
    10 LDL subfractions. PLoS ONE, 10(4), e0120758.
"""
    beta = np.asarray(beta, dtype=float)
    eaf = np.asarray(eaf, dtype=float)
    var_betax = 2.0 * (beta ** 2) * eaf * (1.0 - eaf)
    return var_betax / phenotypic_variance


def per_snp_r2_from_se(
    beta: np.ndarray,
    se: np.ndarray,
    n: np.ndarray,
    eaf: np.ndarray,
) -> np.ndarray:
    """Estimate per-SNP R-squared using SE, N, and EAF.

Parameters
----------
beta : numpy.ndarray
    Effect sizes.
se : numpy.ndarray
    Standard errors.
n : numpy.ndarray
    Sample sizes.
eaf : numpy.ndarray
    Effect allele frequencies.
Returns
-------
numpy.ndarray
    Per-variant R-squared values.
"""
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    n = np.asarray(n, dtype=float)
    eaf = np.asarray(eaf, dtype=float)
    var_betax = 2.0 * (beta ** 2) * eaf * (1.0 - eaf)
    sigma2 = (se ** 2) * 2.0 * n * eaf * (1.0 - eaf)
    return var_betax / (sigma2 + var_betax)


def population_allele_frequency(
    af: float,
    prop: float,
    odds_ratio: float,
    prevalence: float,
) -> float:
    """Back-calculate population allele frequency from case-control summary data.

Parameters
----------
af : float
    Reported allele frequency in the GWAS sample.
prop : float
    Case proportion in the sample.
odds_ratio : float
    Allele odds ratio.
prevalence : float
    Population disease prevalence.
Returns
-------
float
    Estimated population allele frequency, or ``nan`` if no valid solution.

Notes
-----
    Used for binary-trait per-SNP R-squared. Orchestration in
    ``util.util_in_convert_h2._get_per_snp_r2``.
"""
    a_coef = odds_ratio - 1.0
    b_coef = (af + prop) * (1.0 - odds_ratio) - 1.0
    c_coef = odds_ratio * af * prop
    d = (b_coef ** 2) - (4.0 * a_coef * c_coef)
    if d < 0:
        return float("nan")
    sol1 = (-b_coef - np.sqrt(d)) / (2.0 * a_coef)
    sol2 = (-b_coef + np.sqrt(d)) / (2.0 * a_coef)

    for sol in (sol1, sol2):
        a_cell = sol
        c_cell = prop - a_cell
        b_cell = af - a_cell
        d_cell = 1.0 + a_cell - af - prop
        if min(a_cell, b_cell, c_cell, d_cell) >= 0:
            af_case = a_cell / (a_cell + c_cell) if (a_cell + c_cell) > 0 else 0.0
            af_control = b_cell / (b_cell + d_cell) if (b_cell + d_cell) > 0 else 0.0
            return af_case * prevalence + (1.0 - prevalence) * af_control
    return float("nan")
