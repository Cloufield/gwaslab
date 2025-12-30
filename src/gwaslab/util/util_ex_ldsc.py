from typing import TYPE_CHECKING, Tuple, Optional, Dict, Any, Union

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
import copy
import numpy as np
import pandas as pd
import scipy.stats as ss
from gwaslab.info.g_Log import Log

from gwaslab.extension.ldsc.ldsc_sumstats import estimate_h2
from gwaslab.extension.ldsc.ldsc_sumstats import estimate_rg
from gwaslab.extension.ldsc.ldsc_sumstats import cell_type_specific


from gwaslab.io.io_read_ldsc import parse_ldsc_summary
from gwaslab.io.io_read_ldsc import parse_partitioned_ldsc_summary

from gwaslab.util.util_in_filter_value import _filter_values
from gwaslab.util.util_in_filter_value import _filter_palindromic
from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.util.util_in_filter_value import _exclude_sexchr


class ARGS():
    def __init__(self, kwargs=None):
        
        self.out = "ldsc"

        if "bfile" in kwargs.keys():
            self.bfile = kwargs["bfile"]
        else:
            self.bfile = None 
        
        if "l2" in kwargs.keys():
            self.l2 = kwargs["l2"]
        else:
            self.l2 = None 

        if "extract" in kwargs.keys():
            self.extract = kwargs["extract"]
        else:
            self.extract = None 

        if "keep" in kwargs.keys():
            self.keep = kwargs["keep"]
        else:
            self.keep = None 

        if "ld_wind_snps" in kwargs.keys():
            self.ld_wind_snps = kwargs["ld_wind_snps"]
        else:
            self.ld_wind_snps = None 

        if "ld_wind_kb" in kwargs.keys():
            self.ld_wind_kb = kwargs["ld_wind_kb"]
        else:
            self.ld_wind_kb = None 

        if "ld_wind_cm" in kwargs.keys():
            self.ld_wind_cm = kwargs["ld_wind_cm"]
        else:
            self.ld_wind_cm = None 

        if "print_snps" in kwargs.keys():
            self.print_snps = kwargs["print_snps"]
        else:
            self.print_snps = None 

        if "annot" in kwargs.keys():
            self.annot = kwargs["annot"]
        else:
            self.annot = None 

        if "thin_annot" in kwargs.keys():
            self.thin_annot = kwargs["thin_annot"]
        else:
            self.thin_annot = None 

        if "cts_bin" in kwargs.keys():
            self.cts_bin = kwargs["cts_bin"]
        else:
            self.cts_bin = None 

        if "cts_breaks" in kwargs.keys():
            self.cts_breaks = kwargs["cts_breaks"]
        else:
            self.cts_breaks = None 

        if "cts_names" in kwargs.keys():
            self.cts_names = kwargs["cts_names"]
        else:
            self.cts_names = None 

        if "per_allele" in kwargs.keys():
            self.per_allele = kwargs["per_allele"]
        else:
            self.per_allele = None 

        if "pq_exp" in kwargs.keys():
            self.pq_exp = kwargs["pq_exp"]
        else:
            self.pq_exp = None 

        if "no_print_annot" in kwargs.keys():
            self.no_print_annot = kwargs["no_print_annot"]
        else:
            self.no_print_annot = None 

        if "h2" in kwargs.keys():
            self.h2 = kwargs["h2"]
        else:
            self.h2 = None

        if "h2_cts" in kwargs.keys():
            self.h2_cts = kwargs["h2_cts"]
        else:
            self.h2_cts = None

        if "rg" in kwargs.keys():
            self.rg = kwargs["rg"]
        else:
            self.rg = None

        if "ref_ld" in kwargs.keys():
            self.ref_ld = kwargs["ref_ld"]
        else:
            self.ref_ld = None

        if "ref_ld_chr" in kwargs.keys():
            self.ref_ld_chr = kwargs["ref_ld_chr"]
        else:
            self.ref_ld_chr = None

        if "w_ld" in kwargs.keys():
            self.w_ld = kwargs["w_ld"]
        else:
            self.w_ld = None
        
        if "w_ld_chr" in kwargs.keys():
            self.w_ld_chr = kwargs["w_ld_chr"]
        else:
            self.w_ld_chr = None     

        if "overlap_annot" in kwargs.keys():
            self.overlap_annot = kwargs["overlap_annot"]
        else:
            self.overlap_annot = None 

        if "print_coefficients" in kwargs.keys():
            self.print_coefficients = kwargs["print_coefficients"]
        else:
            self.print_coefficients = "ldsc" 

        if "frqfile" in kwargs.keys():
            self.frqfile = kwargs["frqfile"]
        else:
            self.frqfile = None 

        if "frqfile_chr" in kwargs.keys():
            self.frqfile_chr = kwargs["frqfile_chr"]
        else:
            self.frqfile_chr = None 

        if "no_intercept" in kwargs.keys():
            self.no_intercept = kwargs["no_intercept"]
        else:
            self.no_intercept = None 

        if "intercept_h2" in kwargs.keys():
            self.intercept_h2 = kwargs["intercept_h2"]
        else:
            self.intercept_h2 = None 

        if "intercept_gencov" in kwargs.keys():
            self.intercept_gencov = kwargs["intercept_gencov"]
        else:
            self.intercept_gencov = None 

        if "M" in kwargs.keys():
            self.M = kwargs["M"]
        else:
            self.M = None 

        if "two_step" in kwargs.keys():
            self.two_step = kwargs["two_step"]
        else:
            self.two_step = None 

        if "chisq_max" in kwargs.keys():
            self.chisq_max = kwargs["chisq_max"]
        else:
            self.chisq_max= None 

        if "ref_ld_chr_cts" in kwargs.keys():
            self.ref_ld_chr_cts = kwargs["ref_ld_chr_cts"]
        else:
            self.ref_ld_chr_cts = None

        if "print_all_cts" in kwargs.keys():
            self.print_all_cts = kwargs["print_all_cts"]
        else:
            self.print_all_cts = False 

        if "print_cov" in kwargs.keys():
            self.print_cov = kwargs["print_cov"]
        else:
            self.print_cov = None 

        self.print_delete_vals = False
        if "print_delete_vals" in kwargs.keys():
            self.print_delete_vals = kwargs["print_delete_vals"]
        else:
            self.print_delete_vals = False 

        if "chunk_size" in kwargs.keys():
            self.chunk_size = kwargs["chunk_size"]
        else:
            self.chunk_size = 50 

        if "pickle" in kwargs.keys():
            self.pickle = kwargs["pickle"]
        else:
            self.pickle = False 

        if "yes_really" in kwargs.keys():
            self.yes_really = kwargs["yes_really"]
        else:
            self.yes_really = False 

        if "invert_anyway" in kwargs.keys():
            self.invert_anyway = kwargs["invert_anyway"]
        else:
            self.invert_anyway = False 

        if "n_blocks" in kwargs.keys():
            self.n_blocks = kwargs["n_blocks"]
        else:
            self.n_blocks = 200 

        if "not_M_5_50" in kwargs.keys():
            self.not_M_5_50 = kwargs["not_M_5_50"]
        else:
            self.not_M_5_50 = False 

        if "no_check_alleles" in kwargs.keys():
            self.no_check_alleles = kwargs["no_check_alleles"]
        else:
            self.no_check_alleles = False 
        
        if "return_silly_things" in kwargs.keys():
            self.return_silly_things = kwargs["return_silly_things"]
        else:
            self.return_silly_things = False 

        if "samp_prev" in kwargs.keys():
            self.samp_prev = kwargs["samp_prev"]
        else:
            self.samp_prev = None 
        
        if "pop_prev" in kwargs.keys():
            self.pop_prev = kwargs["pop_prev"]
        else:
            self.pop_prev = None 


####################################################################################################################
from gwaslab.qc.qc_decorator import with_logging
@with_logging(
        start_to_msg="run LD score regression",
        finished_msg="running LD score regression",
        start_cols=["CHR","POS","EA","NEA"],
        start_function="estimate_h2_by_ldsc"
)
def _estimate_h2_by_ldsc(
    insumstats: Union['Sumstats', pd.DataFrame], 
    log: Log, 
    meta: Optional[Dict[str, Any]] = None,
    verbose: bool = True, 
    munge: bool = False, 
    munge_kwargs: Optional[Dict[str, Any]] = None, 
    **raw_kwargs: Any
) -> Tuple[float, Optional[pd.DataFrame]]:
    """
    Estimate SNP heritability using LD score regression.

    Parameters
    ----------
    verbose : bool, default=True
        If True, print detailed progress and status messages during execution.
    munge : bool, default=False
        If True, apply standard munging procedures (e.g., filtering, harmonization, and QC)
        to the input summary statistics prior to analysis.
    **raw_kwargs
        Additional keyword arguments forwarded directly to the underlying LDSC call.
        Common options include:
            ref_ld_chr : str or path-like, required
                Path to reference LD score files (directory or specific file prefix).
            w_ld_chr : str or path-like, optional
                Path to LD weight scores. Often the same as `ref_ld`.
            samp_prev : float, optional
                Sample prevalence (case proportion) for case–control summary statistics.
            pop_prev : float, optional
                Population prevalence for case–control traits.

    Returns
    -------
    tuple
        A tuple containing:
        - parsed_summary (pd.DataFrame): Heritability estimate and related statistics
        - results_table (pd.DataFrame or None): Coefficient results if available; otherwise None

    Notes
    -----
    This function wraps the LDSC implementation from Bulik-Sullivan et al. (2015).
    Requires input columns: CHR, POS, EA, NEA.
    For case-control studies, provide samp_prev and pop_prev via meta or raw_kwargs.
    """
    sumstats = insumstats
    kwargs = copy.deepcopy(raw_kwargs)
    
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].fillna(sumstats["N"].median()).apply("int64")

    if munge:
        if munge_kwargs is None:
            munge_kwargs={}
        log.write("Start to munge sumstats.")
        sumstats = _munge_sumstats(sumstats, log=log, verbose=verbose,**munge_kwargs)
        log.write("Finished munging sumstats.")

    log.write(" -Run single variate LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.", verbose=verbose)
    
    if meta["gwaslab"]["sample_prevalence"]!="Unknown" and meta["gwaslab"]["population_prevalence"]!="Unknown" :
        if "samp_prev" not in kwargs.keys():
            kwargs["samp_prev"] = "{}".format(meta["gwaslab"]["sample_prevalence"])
        if "pop_prev" not in kwargs.keys():
            kwargs["pop_prev"] =  "{}".format(meta["gwaslab"]["population_prevalence"])

    log.write(" -Arguments:", verbose=verbose)
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_kwargs = ARGS(kwargs = kwargs)

    # Create Z if not present (munging may have created it from P)
    if "Z" not in sumstats.columns:
        if "BETA" in sumstats.columns and "SE" in sumstats.columns:
            sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
        else:
            raise ValueError("Cannot create Z: need either Z column, or BETA and SE columns")

    # Rename columns to LDSC format if not already renamed (munging may have done this)
    rename_dict = {}
    if "A1" not in sumstats.columns and "EA" in sumstats.columns:
        rename_dict["EA"] = "A1"
    if "A2" not in sumstats.columns and "NEA" in sumstats.columns:
        rename_dict["NEA"] = "A2"
    if "SNP" not in sumstats.columns and "rsID" in sumstats.columns:
        rename_dict["rsID"] = "SNP"
    if rename_dict:
        sumstats = sumstats.rename(columns=rename_dict)

    log.write(" -LDSC log:", verbose=verbose)
    summary = estimate_h2(sumstats, args = default_kwargs, log = log)
    
    results_table = None
    if type(summary) is tuple:
        results_table = summary[1]
        summary = summary[0]
        log.write(" -Coefficient results have been stored in .ldsc_h2_results", verbose=verbose)
        

    log.write(" -Results have been stored in .ldsc_h2", verbose=verbose)
    return parse_ldsc_summary(summary), results_table


####################################################################################################################
@with_logging(
        start_to_msg="run LD score regression",
        finished_msg="running LD score regression",
        start_cols=["CHR","POS","EA","NEA"],
        start_function="estimate_partitioned_h2_by_ldsc"
)
def _estimate_partitioned_h2_by_ldsc(insumstats,  log,  meta=None,verbose=True, **raw_kwargs):
    sumstats = insumstats.copy()
    kwargs = copy.deepcopy(raw_kwargs)
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].fillna(sumstats["N"].median()).apply("int64")
    log.write(" -Run partitioned LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)
    
    if meta["gwaslab"]["sample_prevalence"]!="Unknown" and meta["gwaslab"]["population_prevalence"]!="Unknown" :
        if "samp_prev" not in kwargs.keys():
            kwargs["samp_prev"] = "{}".format(meta["gwaslab"]["sample_prevalence"])
        if "pop_prev" not in kwargs.keys():
            kwargs["pop_prev"] =  "{}".format(meta["gwaslab"]["population_prevalence"])

    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_kwargs = ARGS(kwargs = kwargs)

    # Create Z if not present (munging may have created it from P)
    if "Z" not in sumstats.columns:
        if "BETA" in sumstats.columns and "SE" in sumstats.columns:
            sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
        else:
            raise ValueError("Cannot create Z: need either Z column, or BETA and SE columns")

    # Rename columns to LDSC format if not already renamed (munging may have done this)
    rename_dict = {}
    if "A1" not in sumstats.columns and "EA" in sumstats.columns:
        rename_dict["EA"] = "A1"
    if "A2" not in sumstats.columns and "NEA" in sumstats.columns:
        rename_dict["NEA"] = "A2"
    if "SNP" not in sumstats.columns and "rsID" in sumstats.columns:
        rename_dict["rsID"] = "SNP"
    if rename_dict:
        sumstats = sumstats.rename(columns=rename_dict)
    
    log.write(" -LDSC log:", verbose=verbose)
    summary,results = estimate_h2(sumstats, default_kwargs, log)
    
    log.write(" -Results have been stored in .ldsc_h2", verbose=verbose)
    return parse_partitioned_ldsc_summary(summary), results


####################################################################################################################


@with_logging(
        start_to_msg="run LD score regression for genetic correlation",
        finished_msg="running LD score regression for genetic correlation",
        start_cols=["CHR","POS","EA","NEA"],
        start_function="estimate_rg_by_ldsc"
)
def _estimate_rg_by_ldsc(
    insumstats: Union['Sumstats', pd.DataFrame], 
    other_traits: Union['Sumstats', pd.DataFrame], 
    log: Log, 
    meta: Optional[Dict[str, Any]] = None, 
    verbose: bool = True, 
    **raw_kwargs: Any
) -> pd.DataFrame:
    sumstats = insumstats.copy()
    kwargs = copy.deepcopy(raw_kwargs)
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].fillna(sumstats["N"].median()).apply("int64")

    log.write(" -Run cross-trait LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)

    
    samp_prev_string=""
    pop_prev_string=""

    if meta["gwaslab"]["sample_prevalence"]!="Unknown" and meta["gwaslab"]["population_prevalence"]!="Unknown" :

        if "samp_prev" not in kwargs.keys():
            samp_prev_string =  "{}".format(meta["gwaslab"]["sample_prevalence"])
        if "pop_prev" not in kwargs.keys():
            pop_prev_string =  "{}".format(meta["gwaslab"]["population_prevalence"])
    
    if "rg" in kwargs.keys():
        alias = kwargs["rg"].split(",")[1:]
    else:
        alias=[]
        for index, each_other_sumstats in enumerate(other_traits):
            alias.append(each_other_sumstats.meta["gwaslab"]["study_name"])
        kwargs["rg"]=",".join([meta["gwaslab"]["study_name"]]+alias)
    
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)

    default_kwargs = ARGS(kwargs = kwargs)

    # Create Z if not present (munging may have created it from P)
    if "Z" not in sumstats.columns:
        if "BETA" in sumstats.columns and "SE" in sumstats.columns:
            sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
        else:
            raise ValueError("Cannot create Z: need either Z column, or BETA and SE columns")

    # Rename columns to LDSC format if not already renamed (munging may have done this)
    rename_dict = {}
    if "A1" not in sumstats.columns and "EA" in sumstats.columns:
        rename_dict["EA"] = "A1"
    if "A2" not in sumstats.columns and "NEA" in sumstats.columns:
        rename_dict["NEA"] = "A2"
    if "SNP" not in sumstats.columns and "rsID" in sumstats.columns:
        rename_dict["rsID"] = "SNP"
    if rename_dict:
        sumstats = sumstats.rename(columns=rename_dict)

    other_traits_to_use = []

    for index, each_other_sumstats in enumerate(other_traits):
        log.write(" -Processing sumstats with alias {} ({})".format(alias[index], each_other_sumstats.meta["gwaslab"]["study_name"]))
        if "rsID" not in each_other_sumstats.data.columns:
            to_append = each_other_sumstats.filter_hapmap3(verbose=False).data
        else:        
            to_append = each_other_sumstats.data.copy()
        
        # Rename columns if needed
        rename_dict_other = {}
        if "A1" not in to_append.columns and "EA" in to_append.columns:
            rename_dict_other["EA"] = "A1"
        if "A2" not in to_append.columns and "NEA" in to_append.columns:
            rename_dict_other["NEA"] = "A2"
        if "SNP" not in to_append.columns and "rsID" in to_append.columns:
            rename_dict_other["rsID"] = "SNP"
        if rename_dict_other:
            to_append = to_append.rename(columns=rename_dict_other)
        
        # Create Z if not present
        if "Z" not in to_append.columns:
            if "BETA" in to_append.columns and "SE" in to_append.columns:
                to_append["Z"] = to_append["BETA"]/to_append["SE"]
            else:
                raise ValueError("Cannot create Z for other trait: need either Z column, or BETA and SE columns")

        other_traits_to_use.append(to_append[["SNP","A1","A2","Z","N"]])    
        
        if each_other_sumstats.meta["gwaslab"]["sample_prevalence"]!="Unknown" and each_other_sumstats.meta["gwaslab"]["population_prevalence"]!="Unknown" :
                samp_prev_string +=  ",{}".format(meta["gwaslab"]["sample_prevalence"])
                pop_prev_string += ",{}".format(meta["gwaslab"]["population_prevalence"])

    if len(pop_prev_string.split(",")) == len(other_traits)+1 and len(samp_prev_string.split(",")) == len(other_traits)+1:
        if "samp_prev" not in kwargs.keys():
            log.write(" -{}:{}".format("samp_prev", samp_prev_string), verbose=verbose)
            default_kwargs.samp_prev = samp_prev_string
        if "pop_prev" not in kwargs.keys():
            log.write(" -{}:{}".format("pop_prev", pop_prev_string), verbose=verbose)
            default_kwargs.pop_prev =  pop_prev_string    

    log.write(" -LDSC log:", verbose=verbose)
    summary = estimate_rg(sumstats[["SNP","A1","A2","Z","N"]], other_traits_to_use, default_kwargs, log)[1]
    
    log.write(" -Results have been stored in .ldsc_rg", verbose=verbose)
    return summary


####################################################################################################################
@with_logging(
        start_to_msg="run LD score regression",
        finished_msg="running LD score regression",
        start_cols=["CHR","POS","EA","NEA"],
        start_function="estimate_h2_cts_by_ldsc"
)
def _estimate_h2_cts_by_ldsc(
    insumstats: Union['Sumstats', pd.DataFrame], 
    log: Log, 
    verbose: bool = True, 
    **raw_kwargs: Any
) -> Tuple[Any, Any]:
    sumstats = insumstats.copy()
    kwargs = copy.deepcopy(raw_kwargs)
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].fillna(sumstats["N"].median()).apply("int64")
    log.write(" -Run cell type specific LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Finucane, H. K., Reshef, Y. A., Anttila, V., Slowikowski, K., Gusev, A., Byrnes, A., ... & Price, A. L. (2018). Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature genetics, 50(4), 621-629.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)
    
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_kwargs = ARGS(kwargs = kwargs)

    # Create Z if not present (munging may have created it from P)
    if "Z" not in sumstats.columns:
        if "BETA" in sumstats.columns and "SE" in sumstats.columns:
            sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
        else:
            raise ValueError("Cannot create Z: need either Z column, or BETA and SE columns")

    # Rename columns to LDSC format if not already renamed (munging may have done this)
    rename_dict = {}
    if "A1" not in sumstats.columns and "EA" in sumstats.columns:
        rename_dict["EA"] = "A1"
    if "A2" not in sumstats.columns and "NEA" in sumstats.columns:
        rename_dict["NEA"] = "A2"
    if "SNP" not in sumstats.columns and "rsID" in sumstats.columns:
        rename_dict["rsID"] = "SNP"
    if rename_dict:
        sumstats = sumstats.rename(columns=rename_dict)
    
    log.write(" -LDSC log:", verbose=verbose)
    summary = cell_type_specific(sumstats, default_kwargs, log)
    
    log.write(" -Results have been stored in .ldsc_partitioned_h2", verbose=verbose)
    return summary



def _munge_sumstats(
    sumstats: pd.DataFrame, 
    log: Log, 
    info: float = 0.9, 
    maf: float = 0.01, 
    n: Optional[str] = None, 
    nopalindromic: bool = True,
    exclude_hla: bool = True, 
    exclude_sexchr: bool = True,  
    verbose: bool = True, 
    **kwargs: Any
) -> pd.DataFrame:
    """
    Munge (filter and harmonize) summary statistics following LDSC workflow.
    
    Based on https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame
    log : Log
        Logging object
    info : float, default=0.9
        Minimum INFO score threshold
    maf : float, default=0.01
        Minimum minor allele frequency threshold
    n : float or None, default=None
        Minimum sample size. If None, uses 90th percentile / 1.5
    nopalindromic : bool, default=True
        If True, remove palindromic SNPs
    exclude_hla : bool, default=True
        If True, exclude HLA region
    exclude_sexchr : bool, default=True
        If True, exclude sex chromosomes
    verbose : bool, default=True
        If True, print progress messages
    **kwargs
        Additional arguments (ignored for compatibility)
    
    Returns
    -------
    pd.DataFrame
        Munge'd summary statistics
    """
    sumstats = sumstats.copy()
    original_len = len(sumstats)
    
    # Step 1: Map column names to LDSC standard (EA->A1, NEA->A2, EAF->FRQ)
    rename_dict = {}
    if "EA" in sumstats.columns and "A1" not in sumstats.columns:
        rename_dict["EA"] = "A1"
    if "NEA" in sumstats.columns and "A2" not in sumstats.columns:
        rename_dict["NEA"] = "A2"
    if "EAF" in sumstats.columns and "FRQ" not in sumstats.columns:
        rename_dict["EAF"] = "FRQ"
    if "rsID" in sumstats.columns and "SNP" not in sumstats.columns:
        rename_dict["rsID"] = "SNP"
    if rename_dict:
        sumstats = sumstats.rename(columns=rename_dict)
        log.write(" -Renamed columns: {}".format(rename_dict), verbose=verbose)
    
    # Step 2: Filter P-values (must be > 0 and <= 1)
    if "P" in sumstats.columns:
        valid_p = (sumstats["P"] > 0) & (sumstats["P"] <= 1) & sumstats["P"].notna()
        bad_p = (~valid_p).sum()
        if bad_p > 0:
            log.write(" -WARNING: {} SNPs had P outside of (0,1]. The P column may be mislabeled.".format(bad_p), verbose=verbose)
        sumstats = sumstats.loc[valid_p].copy()
        log.write(" -After P-value filtering: {} SNPs remain".format(len(sumstats)), verbose=verbose)
    
    # Step 3: Filter INFO (>= info_min, warn if outside [0,1.5])
    if "INFO" in sumstats.columns:
        if isinstance(sumstats["INFO"], pd.Series):
            bad_info = ((sumstats["INFO"] > 2.0) | (sumstats["INFO"] < 0)) & sumstats["INFO"].notna()
            valid_info = sumstats["INFO"] >= info
        else:
            # Multiple INFO columns (DataFrame)
            bad_info = (((sumstats["INFO"] > 2.0) & sumstats["INFO"].notna()).any(axis=1) | 
                       ((sumstats["INFO"] < 0) & sumstats["INFO"].notna()).any(axis=1))
            valid_info = (sumstats["INFO"].sum(axis=1) >= info * len(sumstats["INFO"].columns))
        
        bad_info_count = bad_info.sum()
        if bad_info_count > 0:
            log.write(" -WARNING: {} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled.".format(bad_info_count), verbose=verbose)
        
        sumstats = sumstats.loc[valid_info].copy()
        log.write(" -After INFO filtering (>= {}): {} SNPs remain".format(info, len(sumstats)), verbose=verbose)
    
    # Step 4: Filter FRQ/MAF (>= maf_min, warn if outside [0,1], convert to MAF)
    if "FRQ" in sumstats.columns:
        bad_frq = (sumstats["FRQ"] < 0) | (sumstats["FRQ"] > 1)
        bad_frq_count = bad_frq.sum()
        if bad_frq_count > 0:
            log.write(" -WARNING: {} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled.".format(bad_frq_count), verbose=verbose)
        
        # Convert to MAF (minor allele frequency)
        sumstats["FRQ"] = np.minimum(sumstats["FRQ"], 1 - sumstats["FRQ"])
        valid_frq = (sumstats["FRQ"] > maf) & ~bad_frq
        sumstats = sumstats.loc[valid_frq].copy()
        log.write(" -After MAF filtering (>= {}): {} SNPs remain".format(maf, len(sumstats)), verbose=verbose)
    
    # Step 5: Filter alleles (strand-unambiguous SNPs only)
    # Valid SNPs: A/T, C/G, A/C, A/G, T/C, T/G
    # Use A1/A2 if available, otherwise EA/NEA
    allele1_col = "A1" if "A1" in sumstats.columns else ("EA" if "EA" in sumstats.columns else None)
    allele2_col = "A2" if "A2" in sumstats.columns else ("NEA" if "NEA" in sumstats.columns else None)
    
    if allele1_col and allele2_col:
        valid_snps = {'AT', 'TA', 'CG', 'GC', 'AC', 'CA', 'AG', 'GA', 'TC', 'CT', 'TG', 'GT'}
        # Create allele pair string
        allele_pair = (sumstats[allele1_col].astype(str) + sumstats[allele2_col].astype(str)).str.upper()
        valid_alleles = allele_pair.isin(valid_snps)
        sumstats = sumstats.loc[valid_alleles].copy()
        log.write(" -After allele filtering (strand-unambiguous): {} SNPs remain".format(len(sumstats)), verbose=verbose)
    
    # Step 6: Filter palindromic SNPs if requested
    # _filter_palindromic expects EA/NEA, so use those column names
    if nopalindromic and allele1_col and allele2_col:
        # Temporarily rename back to EA/NEA for _filter_palindromic if needed
        temp_rename = {}
        if allele1_col == "A1":
            temp_rename["A1"] = "EA"
        if allele2_col == "A2":
            temp_rename["A2"] = "NEA"
        if temp_rename:
            sumstats = sumstats.rename(columns=temp_rename)
        
        sumstats = _filter_palindromic(sumstats, mode="out", verbose=verbose, log=log)
        
        # Rename back to A1/A2 if we renamed them
        if temp_rename:
            reverse_rename = {v: k for k, v in temp_rename.items()}
            sumstats = sumstats.rename(columns=reverse_rename)
        
        log.write(" -After removing palindromic SNPs: {} SNPs remain".format(len(sumstats)), verbose=verbose)
    
    # Step 7: Process N (filter on sample size)
    if "N" in sumstats.columns:
        if n is None:
            # Use 90th percentile / 1.5 as threshold (LDSC default)
            min_n = sumstats["N"].quantile(0.9) / 1.5
        else:
            min_n = n
        valid_n = sumstats["N"] >= min_n
        sumstats = sumstats.loc[valid_n].copy()
        log.write(" -After N filtering (>= {}): {} SNPs remain".format(min_n, len(sumstats)), verbose=verbose)
    
    # Step 8: Exclude HLA region if requested
    if exclude_hla and "CHR" in sumstats.columns and "POS" in sumstats.columns:
        sumstats = _exclude_hla(sumstats, verbose=verbose, log=log)
        log.write(" -After excluding HLA: {} SNPs remain".format(len(sumstats)), verbose=verbose)
    
    # Step 9: Exclude sex chromosomes if requested
    if exclude_sexchr and "CHR" in sumstats.columns:
        sumstats = _exclude_sexchr(sumstats, verbose=verbose, log=log)
        log.write(" -After excluding sex chromosomes: {} SNPs remain".format(len(sumstats)), verbose=verbose)
    
    # Step 10: Convert P to Z if P exists and Z doesn't
    # Note: If BETA and SE exist, Z should be calculated from them instead
    if "Z" not in sumstats.columns:
        if "BETA" in sumstats.columns and "SE" in sumstats.columns:
            # Prefer BETA/SE for Z calculation (more accurate)
            sumstats["Z"] = sumstats["BETA"] / sumstats["SE"]
            log.write(" -Calculated Z from BETA/SE", verbose=verbose)
        elif "P" in sumstats.columns:
            # Convert P to Z using two-sided test
            # Z = sign(BETA) * sqrt(chi2.isf(P, 1))
            # If we have BETA, use its sign; otherwise use absolute value
            if "BETA" in sumstats.columns:
                z_sign = np.sign(sumstats["BETA"])
                z_abs = ss.norm.isf(sumstats["P"] / 2.0)
                sumstats["Z"] = z_sign * np.abs(z_abs)
            elif "OR" in sumstats.columns:
                # For OR, sign is based on whether OR > 1
                z_sign = np.sign(np.log(sumstats["OR"]))
                z_abs = ss.norm.isf(sumstats["P"] / 2.0)
                sumstats["Z"] = z_sign * np.abs(z_abs)
            else:
                # No sign information, use absolute Z
                sumstats["Z"] = ss.norm.isf(sumstats["P"] / 2.0)
            log.write(" -Converted P to Z scores", verbose=verbose)
    
    # Step 11: Remove duplicates on SNP ID
    if "SNP" in sumstats.columns:
        old_len = len(sumstats)
        sumstats = sumstats.drop_duplicates(subset="SNP", keep="first").reset_index(drop=True)
        new_len = len(sumstats)
        if old_len != new_len:
            log.write(" -Removed {} SNPs with duplicated SNP IDs ({} SNPs remain)".format(old_len - new_len, new_len), verbose=verbose)
    
    # Final summary
    log.write(" -Munging complete: {} SNPs removed, {} SNPs remain (from original {})".format(
        original_len - len(sumstats), len(sumstats), original_len), verbose=verbose)
    
    return sumstats
