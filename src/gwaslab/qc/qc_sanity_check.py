
import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from liftover import get_lifter
from liftover import ChainFile
from functools import partial
from functools import wraps

from gwaslab.g_vchange_status import vchange_status
from gwaslab.g_vchange_status import status_match
from gwaslab.g_vchange_status import change_status
from gwaslab.g_Log import Log
from gwaslab.g_version import _get_version
from gwaslab.g_vchange_status import STATUS_CATEGORIES

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chain
from gwaslab.bd.bd_common_data import NA_STRINGS

from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_build import _set_build
from gwaslab.qc.qc_decorator import with_logging

from gwaslab.util.util_in_fill_data import _convert_betase_to_mlog10p
from gwaslab.util.util_in_fill_data import _convert_betase_to_p
from gwaslab.util.util_in_fill_data import _convert_mlog10p_to_p


def check_na_columns(sumstats, coltocheck=None, log=Log(), verbose=True):
    """
    Check for columns with all NA values and drop them.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame.
    log : Log object
        Logging object for output messages.
    verbose : bool
        If True, write detailed messages.
    
    Returns
    -------
    pd.DataFrame
        Modified DataFrame with all-NA columns removed.
    """
    if coltocheck is None:
        coltocheck = []
    log.write(f" -Checking if any columns are empty...", verbose=verbose)
    for col in coltocheck:
        if col in sumstats.columns:
            if sumstats[col].isna().all():
                log.write(f" -Dropping column '{col}' as all values are NA.", verbose=verbose)
                sumstats = sumstats.drop(columns=[col])
    return sumstats

def add_tolerence(stats, float_tolerance, mode):
    if "l" in mode:
        stats = (stats[0] - float_tolerance if stats[0]!=float("Inf") else float("Inf"), stats[1])
    if "r" in mode:
        stats = (stats[0] , stats[1] + float_tolerance if stats[0]!=float("Inf") else float("Inf"))
    return stats

def check_range(sumstats, var_range, header, coltocheck, cols_to_check, log, verbose, dtype="Int64"):
    pre_number=len(sumstats)
    if header in coltocheck and header in sumstats.columns:
        cols_to_check.append(header)
        if header=="STATUS": 
            log.write(" -Checking STATUS and converting STATUS to Int64....", verbose=verbose) 
            # Convert STATUS to integer (remove Categorical if present)
            if sumstats[header].dtype.name == 'category':
                sumstats[header] = sumstats[header].astype(str).astype(int)
            elif sumstats[header].dtype not in ['int64', 'Int64', 'int32', 'Int32']:
                sumstats[header] = sumstats[header].astype(int)
            sumstats[header] = sumstats[header].astype('Int64')
            return sumstats
        
        if dtype in ["Int64","Int32","int","int32","in64"]:
            log.write(" -Checking if {} <= {} <= {} ...".format( var_range[0] ,header, var_range[1]), verbose=verbose) 
            # Remove commas for string representations of numbers before conversion
            if sumstats[header].dtype == 'object':
                sumstats[header] = sumstats[header].astype(str).str.replace(',', '')
            sumstats[header] = np.floor(pd.to_numeric(sumstats[header], errors='coerce')).astype(dtype)
            is_valid = (sumstats[header]>=var_range[0]) & (sumstats[header]<=var_range[1])
        elif dtype in ["Float64","Float32","float","float64","float32"]:
            log.write(" -Checking if {} < {} < {} ...".format( var_range[0] ,header, var_range[1]),verbose=verbose) 
            sumstats[header] = pd.to_numeric(sumstats[header], errors='coerce').astype(dtype)
            is_valid = (sumstats[header]>var_range[0]) & (sumstats[header]<var_range[1])
        is_valid = is_valid.fillna(False)

        if header=="P":
            is_low_p =  sumstats["P"] == 0 
            if sum(is_low_p) >0:
                log.warning("Extremely low P detected (P=0 or P < minimum positive value of float64) : {}".format(sum(is_low_p)))
                log.warning("Please consider using MLOG10P instead.")
        
        if header=="INFO":
            is_high_info =  sumstats["INFO"]>1 
            if sum(is_high_info) >0:
                log.warning("High INFO detected (INFO>1) : {}".format(sum(is_high_info)))
                log.warning("max(INFO): {}".format(sumstats["INFO"].max()))
                log.warning("Please check if this is as expected.")

        if sum(~is_valid)>0:
            try:
                if "SNPID" in sumstats.columns:
                    id_to_use = "SNPID"
                elif "rsID" in sumstats.columns:
                    id_to_use = "rsID"
                invalid_ids = sumstats.loc[~is_valid, id_to_use].head().astype("string")
                invalid_values = sumstats.loc[~is_valid, header].head().astype("string").fillna("NA")
                log.write("  -Examples of invalid variants({}): {} ...".format(id_to_use, ",".join(invalid_ids.to_list()) ), verbose=verbose) 
                log.write("  -Examples of invalid values ({}): {} ...".format(header, ",".join(invalid_values.to_list()) ), verbose=verbose) 
            except:
                pass
        
        # If all values are invalid, drop the column instead of filtering rows
        if not is_valid.any():
            sumstats = sumstats.drop(columns=[header])
            log.write(" -Dropped column {} as all values are invalid.".format(header), verbose=verbose)
        else:
            sumstats = sumstats.loc[is_valid,:]
            after_number=len(sumstats)
            log.write(" -Removed {} variants with bad/na {}.".format(pre_number - after_number, header), verbose=verbose) 
    return sumstats

@with_logging(
        start_to_msg="perform sanity check for statistics",
        finished_msg="sanity check for statistics",
        start_function=".check_sanity()"
)
def _sanity_check_stats(sumstats,
                     coltocheck=None,
                     n=(0,2**31-1),
                     ncase=(0,2**31-1),
                     ncontrol=(0,2**31-1),
                     eaf=(0,1),
                     mac=(0,2**31-1),
                     maf=(0,0.5),
                     chisq=(0,float("Inf")),
                     z=(-9999,9999),
                     t=(-99999,99999),
                     f=(0,float("Inf")),
                     p=(0,1),
                     mlog10p=(0,99999),
                     beta=(-100,100),
                     se=(0,float("Inf")),
                     OR=(0,100),
                     OR_95L=(0,float("Inf")),
                     OR_95U=(0,float("Inf")),
                     HR=(0,100),
                     HR_95L=(0,float("Inf")),
                     HR_95U=(0,float("Inf")),
                     info=(0,2),
                     float_tolerance = 1e-7,
                     verbose=True,
                     log=Log()):
    '''
    Check whether numerical summary statistics fall within valid ranges.

    This function validates commonly used GWAS fields (sample sizes, allele
    frequencies, effect sizes, test statistics, etc.) against expected numeric
    ranges. Columns not present in the input are ignored.

    Parameters
    ----------
    n : tuple of (float, float), optional
        Valid range for sample size (N).
    ncase : tuple of (float, float), optional
        Valid range for number of cases (N_CASE).
    ncontrol : tuple of (float, float), optional
        Valid range for number of controls (N_CONTROL).
    eaf : tuple of (float, float), optional
        Valid range for effect allele frequency (EAF).
    mac : tuple of (float, float), optional
        Valid range for minor allele count (MAC).
    maf : tuple of (float, float), optional
        Valid range for minor allele frequency (MAF).
    chisq : tuple of (float, float), optional
        Valid range for chi-square statistics (CHISQ).
    z : tuple of (float, float), optional
        Valid range for z-scores (Z).
    t : tuple of (float, float), optional
        Valid range for t-statistics (T).
    f : tuple of (float, float), optional
        Valid range for F-statistics (F).
    p : tuple of (float, float), optional
        Valid range for p-values (P).
    mlog10p : tuple of (float, float), optional
        Valid range for negative log10 p-values (MLOG10P).
    beta : tuple of (float, float), optional
        Valid range for effect size estimates (BETA).
    se : tuple of (float, float), optional
        Valid range for standard errors (SE).
    OR : tuple of (float, float), optional
        Valid range for odds ratios (OR).
    OR_95L : tuple of (float, float), optional
        Valid range for lower bound of 95% CI for OR.
    OR_95U : tuple of (float, float), optional
        Valid range for upper bound of 95% CI for OR.
    HR : tuple of (float, float), optional
        Valid range for hazard ratios (HR).
    HR_95L : tuple of (float, float), optional
        Valid range for lower bound of 95% CI for HR.
    HR_95U : tuple of (float, float), optional
        Valid range for upper bound of 95% CI for HR.
    info : tuple of (float, float), optional
        Valid range for imputation info score (INFO).
    float_tolerance : float, default 0.0
        Numerical tolerance applied when comparing floating-point values.
    verbose : bool, default False
        If True, print progress and warnings.

    Note:
    Sanity check ranges (default; v3.4.33):
        N:      Int32    , N>0 , 
        EAF:    float32  , 0 <= EAF <=1, 
        P:      float64  , 0 <= P <= 1, 
        BETA:   float64  , abs(BETA) <100
        SE:     float64  , SE > 0
        OR:     float64  , np.exp(-100) <OR < np.exp(100)
        OR_95L: float64  , OR_95L >0
        OR_95U: float64  , OR_95L >0
        HR:     float64  , np.exp(-100) <log(HR) <np.exp(100)
        HR_95L: float64  , HR_95L >0
        HR_95U: float64  , HR_95L >0
        INFO:   float32  , 1>=INFO>0
        Z       float64  , -9999 < Z < 9999
        T       float64  , -99999 < T < 99999
        F       float64  , F > 0 

    Returns:
        pd.DataFrame: Modified sumstats with invalid variants removed.
    '''
    log.write(" -Comparison tolerance for floats: {}".format(float_tolerance), verbose=verbose) 
    eaf = add_tolerence(eaf, float_tolerance, "lr")
    maf = add_tolerence(maf, float_tolerance, "lr")
    beta = add_tolerence(beta, float_tolerance, "lr")
    se = add_tolerence(se, float_tolerance, "lr")
    mlog10p = add_tolerence(mlog10p, float_tolerance, "lr")
    OR = add_tolerence(OR, float_tolerance, "lr")
    OR_95L = add_tolerence(OR_95L, float_tolerance, "lr")
    OR_95U = add_tolerence(OR_95U, float_tolerance, "lr")
    HR = add_tolerence(HR, float_tolerance, "lr")
    HR_95L = add_tolerence(HR_95L, float_tolerance, "lr")
    HR_95U = add_tolerence(HR_95U, float_tolerance, "lr")
    info = add_tolerence(info, float_tolerance, "lr")
    z = add_tolerence(z, float_tolerance, "lr")
    p = add_tolerence(p, float_tolerance, "lr")
    f = add_tolerence(f, float_tolerance, "lr")
    chisq = add_tolerence(chisq, float_tolerance, "lr")
    ############################################################################################
    ## add direction
    if coltocheck is None:
        coltocheck = ["P","MLOG10P","INFO","Z","BETA","SE","EAF","CHISQ","F","N","N_CASE","N_CONTROL","OR","OR_95L","OR_95U","HR","HR_95L","HR_95U","STATUS"]
    
    cols_to_check=[]
    oringinal_number=len(sumstats)
    sumstats = sumstats.copy()
    sumstats = check_na_columns(sumstats, coltocheck = coltocheck, log=log, verbose=verbose)
    
    ###Int64 ################################################################################################################################################
    sumstats = check_range(sumstats, var_range=n, header="N", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")
    sumstats = check_range(sumstats, var_range=ncase, header="N_CASE", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")
    sumstats = check_range(sumstats, var_range=ncontrol, header="N_CONTROL", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")

    ###float32 ################################################################################################################################################
    sumstats = check_range(sumstats, var_range=eaf, header="EAF", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float32")
    sumstats = check_range(sumstats, var_range=maf, header="MAF", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float32")
    sumstats = check_range(sumstats, var_range=info, header="INFO", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float32")

    ###float64 ################################################################################################################################################
    sumstats = check_range(sumstats, var_range=chisq, header="CHISQ", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=z, header="Z", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=t, header="T", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=f, header="F", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=p, header="P", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=mlog10p, header="MLOG10P", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=beta, header="BETA", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=se, header="SE", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=OR, header="OR", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=OR_95L, header="OR_95L", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=OR_95U, header="OR_95U", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=HR, header="HR", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=HR_95L, header="HR_95L", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=HR_95U, header="HR_95U", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    ###STATUS ###############################################################################################################################################
    sumstats = check_range(sumstats, var_range=None, header="STATUS", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")

    after_number=len(sumstats)
    log.write(" -Removed "+str(oringinal_number - after_number)+" variants with bad statistics in total.",verbose=verbose) 
    log.write(" -Data types for each column:",verbose=verbose)
    
    return sumstats

### check consistency #############################################################################################################################################

@with_logging(
        start_to_msg="check data consistency across columns",
        finished_msg="checking data consistency across columns",
        start_function=".check_data_consistency()"
)
def _check_data_consistency(sumstats, beta="BETA", se="SE", p="P",mlog10p="MLOG10P",rtol=1e-3, atol=1e-3, equal_nan=True, verbose=True,log=Log()):
    '''
    Check consistency between related statistical values. Minor inconsistencies are likely due to rounding.

    Parameters
    ----------
    rtol : float, default 1e-05
        Relative tolerance for comparisons (fractional).
    atol : float, default 1e-08
        Absolute tolerance for comparisons (absolute threshold).
    equal_nan : bool, default False
        If True, treat NaN values as equal during consistency checks.
    verbose : bool, default False
        If True, print progress or warning messages.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table including annotations for detected inconsistencies.
    '''

    log.write(" -Tolerance: {} (Relative) and {} (Absolute)".format(rtol, atol),verbose=verbose)
    check_status = 0
    
    if "SNPID" in sumstats.columns:
        id_to_use = "SNPID"
    elif "rsID" in sumstats.columns:
        id_to_use = "rsID"
    else:
        log.write(" -SNPID/rsID not available...SKipping",verbose=verbose)
        log.write("Finished checking data consistency across columns.",verbose=verbose) 
        return 0
    
    
    if "BETA" in sumstats.columns and "SE" in sumstats.columns:
        if "MLOG10P" in sumstats.columns:
            log.write(" -Checking if BETA/SE-derived-MLOG10P is consistent with MLOG10P...",verbose=verbose)
            betase_derived_mlog10p =  _convert_betase_to_mlog10p(sumstats["BETA"], sumstats["SE"])
            is_close = np.isclose(betase_derived_mlog10p, sumstats["MLOG10P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
            diff = betase_derived_mlog10p - sumstats["MLOG10P"]
            if sum(~is_close)>0:
                log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
                log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
            else:
                log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
            check_status=1
        
        if "P" in sumstats.columns:
            log.write(" -Checking if BETA/SE-derived-P is consistent with P...",verbose=verbose)
            betase_derived_p =  _convert_betase_to_p(sumstats["BETA"], sumstats["SE"])
            is_close = np.isclose(betase_derived_p, sumstats["P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
            diff = betase_derived_p - sumstats["P"]
            if sum(~is_close)>0:
                log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
                log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
            else:
                log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
            check_status=1
    
    if "MLOG10P" in sumstats.columns and "P" in sumstats.columns:
        log.write(" -Checking if MLOG10P-derived-P is consistent with P...",verbose=verbose)
        mlog10p_derived_p = _convert_mlog10p_to_p(sumstats["MLOG10P"])
        is_close = np.isclose(mlog10p_derived_p, sumstats["P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
        diff = mlog10p_derived_p - sumstats["P"]
        if sum(~is_close)>0:
            log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
            log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
        else:
            log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
        check_status=1

    if "N" in sumstats.columns and "N_CONTROL" in sumstats.columns and "N_CASE" in sumstats.columns:
        log.write(" -Checking if N is consistent with N_CASE + N_CONTROL ...", verbose=verbose) 
        is_close = sumstats["N"] == sumstats["N_CASE"] + sumstats["N_CONTROL"] 
        #is_close = np.isclose(sumstats["N"], sumstats["N_CASE"] + sumstats["N_CONTROL"] , rtol=rtol, atol=atol, equal_nan=equal_nan)
        diff = abs(sumstats["N"] - (sumstats["N_CASE"] + sumstats["N_CONTROL"] ))
        if sum(~is_close)>0:
            log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
            log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
        else:
            log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
        check_status=1
        
    if check_status==1:
        log.write(" -Note: if the max difference is greater than expected, please check your original sumstats.",verbose=verbose)
    else:
        log.write(" -No availalbe columns for data consistency checking...Skipping...",verbose=verbose)
