from typing import Optional, Union, Tuple
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log

def get_power(
              mode: str = "b",
              genotype_rr: Optional[float] = None,
              genotype_or: Optional[float] = None,
              beta: float = 0.3,
              eaf: float = 0.1,
              n: int = 10000,
              ncase: int = 2000,
              ncontrol: int = 15000,
              prevalence: float = 0.15,
              or_to_rr: bool = False,
              daf: float = 0.2,
              sig_level: float = 5e-8,
              vary: float = 1,
              log: Log = Log(),
              verbose: bool = True
             ) -> Union[float, np.ndarray]:
    """Calculate statistical power for genetic association studies.
    
Parameters
----------
mode : {'b', 'q'}, optional
    Calculation mode:
    - 'b' for binary traits (case-control)
    - 'q' for quantitative traits
genotype_rr : float, optional
    Genotype relative risk (GRR) for risk allele for binary traits
genotype_or : float, optional
    Genotype odds ratio (OR) for risk allele for binary traits
beta : float, optional
    Effect size (log scale) for quantitative traits
eaf : float, optional
    Effect allele frequency for quantitative traits
n : int, optional
    Total sample size for quantitative traits
ncase : int, optional
    Number of cases for binary mode
ncontrol : int, optional
    Number of controls for binary mode
prevalence : float, optional
    Disease prevalence in population
or_to_rr : bool, optional
    Convert OR to GRR using prevalence (Zhang & Kai method)
daf : float, optional
    Derived allele frequency
sig_level : float, optional
    Significance threshold (default: 5e-8)
vary : float, optional
    Phenotype variance for quantitative traits
Returns
-------
float or array-like
    Calculated statistical power (0-1)
    
References
----------
    Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). 
    Joint analysis is more efficient than replication-based analysis 
    for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
    
    Zhang, J., & Kai, F. Y. (1998). What's the relative risk? 
    A method of correcting the odds ratio in cohort studies of common outcomes. 
    Jama, 280(19), 1690-1691.
"""
    log.write(" Start to calculate statistical power...", verbose=verbose)
    if mode=="b":
        from gwaslab.algorithm.core.power import (
            genotype_disease_probabilities,
            power_binary,
            resolve_genotype_rr,
        )

        log.write(" -Input settings (b mode):", verbose=verbose)
        log.write("  -Number of cases:{}".format(ncase), verbose=verbose)
        log.write("  -Number of controls:{}".format(ncontrol), verbose=verbose)
        if genotype_rr is not None:
            log.write("  -Risk allele RR:{}".format(genotype_rr), verbose=verbose)
        elif genotype_or is not None:
            log.write("  -Risk allele OR:{}".format(genotype_or), verbose=verbose)
        elif beta is not None:
            log.write("  -Risk allele beta:{}".format(beta), verbose=verbose)
        else:
            log.write("  -Risk allele RR:{}".format(0.1), verbose=verbose)

        log.write("  -Disease prevalence:{}".format(prevalence), verbose=verbose)
        log.write("  -Risk allele frequency: {}".format(daf), verbose=verbose)
        log.write("  -Significance level: {}".format(sig_level), verbose=verbose)

        if genotype_rr is None:
            if or_to_rr is False:
                log.write(" -Alogorithm: Skol, Andrew D., et al. Nature genetics 38.2 (2006): 209-213....", verbose=verbose)
                log.write(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....", verbose=verbose)
            else:
                log.write(" -OR is converted to GRR using base prevalence: {}".format(prevalence), verbose=verbose)
                log.write(" -Alogorithm: Zhang, J., & Kai, F. Y. (1998). What's the relative risk?: A method of correcting the odds ratio in cohort studies of common outcomes. Jama, 280(19), 1690-1691.....", verbose=verbose)

        resolved_rr = resolve_genotype_rr(
            genotype_rr=genotype_rr,
            genotype_or=genotype_or,
            beta=beta,
            prevalence=prevalence,
            use_or_to_rr=or_to_rr,
        )
        _, _, _, aap, abp, bbp = genotype_disease_probabilities(
            resolved_rr, daf, prevalence
        )
        log.write("Probability of disease :", verbose=verbose)
        log.write(" - Individuals with AA genotype: {}".format(aap), verbose=verbose)
        log.write(" - Individuals with AB genotype: {}".format(abp), verbose=verbose)
        log.write(" - Individuals with BB genotype: {}".format(bbp), verbose=verbose)

        power = power_binary(
            genotype_rr=resolved_rr,
            daf=daf,
            ncase=ncase,
            ncontrol=ncontrol,
            prevalence=prevalence,
            sig_level=sig_level,
        )
        log.write("Expected power: {}".format(power), verbose=verbose)

    elif mode=="q":
        from gwaslab.algorithm.core.power import power_quantitative

        if beta is None:
            beta = 0.1

        log.write(" -Input settings (q mode):", verbose=verbose)
        log.write("  -Significance level: {}".format(sig_level), verbose=verbose)
        log.write("  -EAF: {}".format(eaf), verbose=verbose)
        log.write("  -BETA: {}".format(beta), verbose=verbose)
        log.write("  -N: {}".format(n), verbose=verbose)
        log.write("  -SNPR2: {}".format(2*eaf*(1-eaf)*(beta**2)), verbose=verbose)
        power = power_quantitative(beta=beta, eaf=eaf, n=n, sig_level=sig_level, phenotypic_variance=vary)
    log.write("Finished calculating statistical power.", verbose=verbose)
    return power

def get_beta(
              mode: str = "q",
              eaf_range: Tuple[float, float] = (0.0001, 0.5),
              beta_range: Tuple[float, float] = (0.0001, 10),
              t: float = 0,
              n: Optional[int] = None,
              sig_level: float = 5e-8,
              vary: float = 1,
              log: Log = Log(),
              verbose: bool = True,
              n_matrix: int = 500
             ) -> pd.DataFrame:
    """Calculate effect size (beta) values that achieve a target statistical power.
    
Parameters
----------
mode : {'q', 'b'}, optional
    Calculation mode:
    - 'q' for quantitative traits (default)
    - 'b' for binary traits (case-control)
eaf_range : tuple of float, optional
    Range of effect allele frequencies to evaluate (min, max)
beta_range : tuple of float, optional
    Range of beta values to evaluate (min, max)
t : float, optional
    Target power value (0-1) to find corresponding eaf-beta combinations
n : int, optional
    Sample size for quantitative traits
sig_level : float, optional
    Significance threshold (default: 5e-8)
vary : float, optional
    Phenotype variance for quantitative traits
log : gwaslab.g_Log.Log object, optional
    Logging object for output messages
verbose : bool, optional
    If True, writes progress messages to log
n_matrix : int, optional
    Size of the eaf-beta matrix for calculation resolution
Returns
-------
pandas.DataFrame
    DataFrame containing eaf-beta combinations that achieve the target power
    
References
----------
    Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). 
    Joint analysis is more efficient than replication-based analysis 
    for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
"""
    if mode=="q":
        if t >0:
            from gwaslab.algorithm.core.power import (
                power_contour_at_target,
                power_grid_quantitative,
            )

            log.write(" -Updating eaf-beta matrix...", verbose=verbose)
            grid, eafs, betas = power_grid_quantitative(
                eaf_range=eaf_range,
                beta_range=beta_range,
                n=n,
                sig_level=sig_level,
                phenotypic_variance=vary,
                n_matrix=n_matrix,
            )

            log.write(" -Extracting eaf-beta combinations with power = {}...".format(t), verbose=verbose)
            eaf_beta = power_contour_at_target(grid, eafs, betas, t)
        else:
            eaf_beta = []
        return pd.DataFrame(eaf_beta)

def get_beta_binary(
              prevalence: Optional[float] = None,
              ncase: Optional[int] = None,
              ncontrol: Optional[int] = None,
              eaf_range: Tuple[float, float] = (0.0001, 0.5),
              beta_range: Tuple[float, float] = (0.0001, 10),
              t: float = 0,
              sig_level: float = 5e-8,
              vary: float = 1,
              log: Log = Log(),
              verbose: bool = True,
              n_matrix: int = 500,
              or_to_rr: bool = False
             ) -> pd.DataFrame:
    """Calculate effect size (beta) values that achieve a target statistical power for binary traits.
    
Parameters
----------
prevalence : float, optional
    Disease prevalence in population
ncase : int, optional
    Number of cases for binary traits
ncontrol : int, optional
    Number of controls for binary traits
eaf_range : tuple of float, optional
    Range of derived allele frequencies to evaluate (min, max)
beta_range : tuple of float, optional
    Range of beta values to evaluate (min, max)
t : float, optional
    Target power value (0-1) to find corresponding eaf-beta combinations
sig_level : float, optional
    Significance threshold (default: 5e-8)
vary : float, optional
    Phenotype variance for quantitative traits
log : gwaslab.g_Log.Log object, optional
    Logging object for output messages
verbose : bool, optional
    If True, writes progress messages to log
n_matrix : int, optional
    Size of the eaf-beta matrix for calculation resolution
or_to_rr : bool, optional
    Convert OR to GRR using prevalence (Zhang & Kai method)
Returns
-------
pandas.DataFrame
    DataFrame containing eaf-beta combinations that achieve the target power
    
References
----------
    Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). 
    Joint analysis is more efficient than replication-based analysis 
    for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
    
    Zhang, J., & Kai, F. Y. (1998). What's the relative risk? 
    A method of correcting the odds ratio in cohort studies of common outcomes. 
    Jama, 280(19), 1690-1691.
"""
    if t >0:
        from gwaslab.algorithm.core.power import (
            power_contour_at_target,
            power_grid_binary,
        )

        log.write(" -Updating eaf-beta matrix...", verbose=verbose)
        if or_to_rr is False:
            log.write(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....", verbose=verbose)
        else:
            log.write(" -OR is converted to GRR using base prevalence: {}".format(prevalence), verbose=verbose)

        grid, eafs, betas = power_grid_binary(
            eaf_range=eaf_range,
            beta_range=beta_range,
            ncase=ncase,
            ncontrol=ncontrol,
            prevalence=prevalence,
            sig_level=sig_level,
            n_matrix=n_matrix,
            use_or_to_rr=or_to_rr,
        )

        log.write(" -Extracting eaf-beta combinations with power = {}...".format(t), verbose=verbose)
        eaf_beta = power_contour_at_target(grid, eafs, betas, t)
    else:
        eaf_beta = []
    return pd.DataFrame(eaf_beta)
