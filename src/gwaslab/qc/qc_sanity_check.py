from typing import TYPE_CHECKING, Union, Optional, List, Tuple
import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from functools import partial
from functools import wraps

from gwaslab.info.g_vchange_status import vchange_status
from gwaslab.info.g_vchange_status import status_match
from gwaslab.info.g_vchange_status import change_status
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version
from gwaslab.info.g_vchange_status import STATUS_CATEGORIES

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
from gwaslab.qc.qc_reserved_headers import get_default_sanity_ranges

from gwaslab.util.util_in_fill_data import _convert_betase_to_mlog10p
from gwaslab.util.util_in_fill_data import _convert_betase_to_p
from gwaslab.util.util_in_fill_data import _convert_mlog10p_to_p


# ----- Internal Helper Functions -----
from gwaslab.io.io_input_type import _get_id_column

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def check_na_columns(sumstats: pd.DataFrame, coltocheck: Optional[List[str]] = None, log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
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

def add_tolerence(stats: Tuple[float, float], float_tolerance: float, mode: str) -> Tuple[float, float]:
    if "l" in mode:
        stats = (stats[0] - float_tolerance if stats[0] != float("Inf") else float("Inf"), stats[1])
    if "r" in mode:
        # Bug fix: Check stats[1] for right side, not stats[0]
        stats = (stats[0], stats[1] + float_tolerance if stats[1] != float("Inf") else float("Inf"))
    return stats

def check_range(sumstats: pd.DataFrame, var_range: Optional[Tuple[Optional[float], Optional[float]]], header: str, coltocheck: List[str], cols_to_check: List[str], log: Log, verbose: bool, dtype: str = "Int64") -> pd.DataFrame:
    pre_number=len(sumstats)
    # Performance optimization: Cache column existence check
    has_header = header in sumstats.columns
    if header in coltocheck and has_header:
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
        
        # Skip range checking if var_range is None or contains None values
        if var_range is None or var_range[0] is None or var_range[1] is None:
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
            is_low_p = sumstats["P"] == 0 
            # Performance optimization: Use .any() instead of sum() > 0
            if is_low_p.any():
                low_p_num = is_low_p.sum()
                log.warning("Extremely low P detected (P=0 or P < minimum positive value of float64) : {}".format(low_p_num))
                log.warning("Please consider using MLOG10P instead.")
            
            # Check for bounded P values (many identical minimum values)
            non_na_p = sumstats["P"].dropna()
            if len(non_na_p) > 0:
                min_p = non_na_p.min()
                is_min_p = (sumstats["P"] == min_p) & sumstats["P"].notna()
                if is_min_p.any():
                    min_p_count = is_min_p.sum()
                    min_p_proportion = min_p_count / len(non_na_p)
                    # Warn if more than 1% of values are at the minimum (excluding the P=0 case already handled above)
                    if min_p_proportion > 0.01 and min_p > 0:
                        log.warning("Bounded P values detected: {} ({:.2%}) variants have P = {} (minimum value)".format(
                            min_p_count, min_p_proportion, min_p))
                        log.warning("This suggests P values may be truncated/capped at a lower limit.")
        
        if header=="MLOG10P":
            # Check for bounded MLOG10P values (many identical maximum values)
            non_na_mlog10p = sumstats["MLOG10P"].dropna()
            if len(non_na_mlog10p) > 0:
                max_mlog10p = non_na_mlog10p.max()
                is_max_mlog10p = (sumstats["MLOG10P"] == max_mlog10p) & sumstats["MLOG10P"].notna()
                if is_max_mlog10p.any():
                    max_mlog10p_count = is_max_mlog10p.sum()
                    max_mlog10p_proportion = max_mlog10p_count / len(non_na_mlog10p)
                    # Warn if more than 1% of values are at the maximum
                    if max_mlog10p_proportion > 0.01:
                        log.warning("Bounded MLOG10P values detected: {} ({:.2%}) variants have MLOG10P = {} (maximum value)".format(
                            max_mlog10p_count, max_mlog10p_proportion, max_mlog10p))
                        log.warning("This suggests MLOG10P values may be truncated/capped at an upper limit.")
        
        if header=="INFO":
            is_high_info = sumstats["INFO"] > 1 
            # Performance optimization: Use .any() instead of sum() > 0
            if is_high_info.any():
                high_info_num = is_high_info.sum()
                log.warning("High INFO detected (INFO>1) : {}".format(high_info_num))
                log.warning("max(INFO): {}".format(sumstats["INFO"].max()))
                log.warning("Please check if this is as expected.")

        # Performance optimization: Use .any() instead of sum() > 0 for boolean check
        invalid = ~is_valid
        if invalid.any():
            try:
                id_to_use = _get_id_column(sumstats)
                # Performance optimization: Reuse invalid mask instead of recomputing
                invalid_ids = sumstats.loc[invalid, id_to_use].head().astype("string")
                invalid_values = sumstats.loc[invalid, header].head().astype("string").fillna("NA")
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
def _sanity_check_stats(sumstats_obj: Union['Sumstats', pd.DataFrame],
                     coltocheck: Optional[List[str]] = None,
                     n: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     ncase: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     ncontrol: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     eaf: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     mac: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     maf: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     chisq: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     z: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     t: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     f: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     p: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     mlog10p: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     beta: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     se: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     OR: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     OR_95L: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     OR_95U: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     HR: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     HR_95L: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     HR_95U: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     info: Optional[Tuple[Optional[float], Optional[float]]] = None,
                     float_tolerance: float = 1e-7,
                     verbose: bool = True,
                     log: Log = Log()) -> pd.DataFrame:
    '''
    Check whether numerical summary statistics fall within valid ranges.

    This function validates commonly used GWAS fields (sample sizes, allele
    frequencies, effect sizes, test statistics, etc.) against expected numeric
    ranges. Columns not present in the input are ignored.

    Parameters
    ----------
    n : tuple of (float, float), optional
        Valid range for sample size (N). Default from qc_researved_header_json: (0, 2147483647).
    ncase : tuple of (float, float), optional
        Valid range for number of cases (N_CASE). Default from qc_researved_header_json: (0, 2147483647).
    ncontrol : tuple of (float, float), optional
        Valid range for number of controls (N_CONTROL). Default from qc_researved_header_json: (0, 2147483647).
    eaf : tuple of (float, float), optional
        Valid range for effect allele frequency (EAF). Default from qc_researved_header_json: (0, 1).
    mac : tuple of (float, float), optional
        Valid range for minor allele count (MAC). Default: (0, 2**31-1).
    maf : tuple of (float, float), optional
        Valid range for minor allele frequency (MAF). Default from qc_researved_header_json: (0, 0.5).
    chisq : tuple of (float, float), optional
        Valid range for chi-square statistics (CHISQ). Default from qc_researved_header_json: (0, inf).
    z : tuple of (float, float), optional
        Valid range for z-scores (Z). Default from qc_researved_header_json: (-9999, 9999).
    t : tuple of (float, float), optional
        Valid range for t-statistics (T). Default from qc_researved_header_json: (-99999, 99999).
    f : tuple of (float, float), optional
        Valid range for F-statistics (F). Default from qc_researved_header_json: (0, inf).
    p : tuple of (float, float), optional
        Valid range for p-values (P). Default from qc_researved_header_json: (0, 1).
    mlog10p : tuple of (float, float), optional
        Valid range for negative log10 p-values (MLOG10P). Default from qc_researved_header_json: (0, 99999).
    beta : tuple of (float, float), optional
        Valid range for effect size estimates (BETA). Default from qc_researved_header_json: (-100, 100).
    se : tuple of (float, float), optional
        Valid range for standard errors (SE). Default from qc_researved_header_json: (0, inf).
    OR : tuple of (float, float), optional
        Valid range for odds ratios (OR). Default from qc_researved_header_json: (0, 100).
    OR_95L : tuple of (float, float), optional
        Valid range for lower bound of 95% CI for OR. Default from qc_researved_header_json: (0, inf).
    OR_95U : tuple of (float, float), optional
        Valid range for upper bound of 95% CI for OR. Default from qc_researved_header_json: (0, inf).
    HR : tuple of (float, float), optional
        Valid range for hazard ratios (HR). Default from qc_researved_header_json: (0, 100).
    HR_95L : tuple of (float, float), optional
        Valid range for lower bound of 95% CI for HR. Default from qc_researved_header_json: (0, inf).
    HR_95U : tuple of (float, float), optional
        Valid range for upper bound of 95% CI for HR. Default from qc_researved_header_json: (0, inf).
    info : tuple of (float, float), optional
        Valid range for imputation info score (INFO). Default from qc_researved_header_json: (0, 2).
    float_tolerance : float, default 1e-7
        Numerical tolerance applied when comparing floating-point values.
    verbose : bool, default True
        If True, print progress and warnings.

    Note:
    Default sanity check ranges are loaded from qc_researved_header_json (single source of truth).
    All range parameters default to None and are automatically loaded from the JSON file if not provided.
    Users can override any range by explicitly passing a tuple value. 

    Returns:
        pd.DataFrame: Modified sumstats with invalid variants removed.
    '''
    # Load default ranges from JSON file (single source of truth)
    default_ranges = get_default_sanity_ranges()
    
    # Use provided values or defaults from JSON
    n = n if n is not None else default_ranges.get('N', (0, 2**31-1))
    ncase = ncase if ncase is not None else default_ranges.get('N_CASE', (0, 2**31-1))
    ncontrol = ncontrol if ncontrol is not None else default_ranges.get('N_CONTROL', (0, 2**31-1))
    eaf = eaf if eaf is not None else default_ranges.get('EAF', (0, 1))
    mac = mac if mac is not None else (0, 2**31-1)  # MAC not in JSON, use default
    maf = maf if maf is not None else default_ranges.get('MAF', (0, 0.5))
    chisq = chisq if chisq is not None else default_ranges.get('CHISQ', (0, float("Inf")))
    z = z if z is not None else default_ranges.get('Z', (-9999, 9999))
    t = t if t is not None else default_ranges.get('T', (-99999, 99999))
    f = f if f is not None else default_ranges.get('F', (0, float("Inf")))
    p = p if p is not None else default_ranges.get('P', (0, 1))
    mlog10p = mlog10p if mlog10p is not None else default_ranges.get('MLOG10P', (0, 99999))
    beta = beta if beta is not None else default_ranges.get('BETA', (-100, 100))
    se = se if se is not None else default_ranges.get('SE', (0, float("Inf")))
    OR = OR if OR is not None else default_ranges.get('OR', (0, 100))
    OR_95L = OR_95L if OR_95L is not None else default_ranges.get('OR_95L', (0, float("Inf")))
    OR_95U = OR_95U if OR_95U is not None else default_ranges.get('OR_95U', (0, float("Inf")))
    HR = HR if HR is not None else default_ranges.get('HR', (0, 100))
    HR_95L = HR_95L if HR_95L is not None else default_ranges.get('HR_95L', (0, float("Inf")))
    HR_95U = HR_95U if HR_95U is not None else default_ranges.get('HR_95U', (0, float("Inf")))
    info = info if info is not None else default_ranges.get('INFO', (0, 2))
    
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
    
    # Handle both DataFrame and Sumstats object
    import pandas as pd
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    cols_to_check=[]
    original_number=len(sumstats)
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
    log.write(" -Removed "+str(original_number - after_number)+" variants with bad statistics in total.",verbose=verbose) 
    log.write(" -Data types for each column:",verbose=verbose)
    
    # Performance optimization: is_dataframe already determined above, no need to recompute
    
    if not is_dataframe:
        # Assign filtered dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
    
    # Update QC status
    try:
        from gwaslab.info.g_meta import _update_qc_step
        sanity_kwargs = {
            'coltocheck': coltocheck, 'n': n, 'ncase': ncase, 'ncontrol': ncontrol, 'eaf': eaf,
            'mac': mac, 'maf': maf, 'chisq': chisq, 'z': z, 't': t, 'f': f, 'p': p, 'mlog10p': mlog10p,
            'beta': beta, 'se': se, 'OR': OR, 'OR_95L': OR_95L, 'OR_95U': OR_95U, 'HR': HR,
            'HR_95L': HR_95L, 'HR_95U': HR_95U, 'info': info, 'float_tolerance': float_tolerance
        }
        _update_qc_step(sumstats_obj, "sanity", sanity_kwargs, True)
    except:
        pass
    
    if not is_dataframe:
        return sumstats_obj.data
    else:
        return sumstats

### check consistency #############################################################################################################################################

@with_logging(
        start_to_msg="check data consistency across columns",
        finished_msg="checking data consistency across columns",
        start_function=".check_data_consistency()"
)
def _check_data_consistency(sumstats_obj: Union['Sumstats', pd.DataFrame], beta: str = "BETA", se: str = "SE", p: str = "P", mlog10p: str = "MLOG10P", rtol: float = 1e-3, atol: float = 1e-3, equal_nan: bool = True, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    '''
    Check consistency between related statistical values. Minor inconsistencies are likely due to rounding.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to check.
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
    sumstats = sumstats_obj.data

    log.write(" -Tolerance: {} (Relative) and {} (Absolute)".format(rtol, atol),verbose=verbose)
    check_status = 0
    
    try:
        id_to_use = _get_id_column(sumstats)
    except (KeyError, AttributeError):
        log.write(" -SNPID/rsID not available...SKipping",verbose=verbose)
        log.write("Finished checking data consistency across columns.",verbose=verbose) 
        return 0
    
    # Performance optimization: Cache column existence checks
    has_beta = "BETA" in sumstats.columns
    has_se = "SE" in sumstats.columns
    has_mlog10p = "MLOG10P" in sumstats.columns
    has_p = "P" in sumstats.columns
    has_n = "N" in sumstats.columns
    has_n_case = "N_CASE" in sumstats.columns
    has_n_control = "N_CONTROL" in sumstats.columns
    
    if has_beta and has_se:
        if has_mlog10p:
            log.write(" -Checking if BETA/SE-derived-MLOG10P is consistent with MLOG10P...",verbose=verbose)
            betase_derived_mlog10p =  _convert_betase_to_mlog10p(sumstats["BETA"], sumstats["SE"])
            is_close = np.isclose(betase_derived_mlog10p, sumstats["MLOG10P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
            # Performance optimization: Use .any() instead of sum() > 0 for boolean check
            inconsistent = ~is_close
            if inconsistent.any():
                inconsistent_num = inconsistent.sum()
                diff = betase_derived_mlog10p - sumstats["MLOG10P"]
                log.write("  -Potentially inconsistent (likely due to rounding): {} variant(s)".format(inconsistent_num),verbose=verbose)
                log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
            else:
                log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
            check_status=1
        
        if has_p:
            log.write(" -Checking if BETA/SE-derived-P is consistent with P...",verbose=verbose)
            betase_derived_p =  _convert_betase_to_p(sumstats["BETA"], sumstats["SE"])
            
            # For P-values, use fold change as the primary consistency metric
            # since P-values span many orders of magnitude (e.g., 1e-300 to 1)
            p_values = sumstats["P"].copy()
            # Ensure derived_p_values is a pandas Series with same index
            if isinstance(betase_derived_p, pd.Series):
                derived_p_values = betase_derived_p.copy()
            else:
                derived_p_values = pd.Series(betase_derived_p, index=sumstats.index)
            
            # Handle NaN values - treat as consistent if both are NaN
            both_na = p_values.isna() & derived_p_values.isna()
            one_na = p_values.isna() | derived_p_values.isna()
            valid_mask = ~one_na & (p_values > 0) & (derived_p_values > 0)
            
            # Initialize consistency check: NaN pairs are considered consistent
            is_close = pd.Series(both_na, index=sumstats.index)
            
            # Initialize fold_change Series for all indices
            fold_change = pd.Series(index=sumstats.index, dtype=float)
            
            if valid_mask.any():
                # Calculate fold change: max(derived_p/p, p/derived_p) for symmetric measure
                # 1.0 = identical, >1.0 = different
                valid_p = p_values[valid_mask]
                valid_derived_p = derived_p_values[valid_mask]
                fold_change[valid_mask] = np.maximum(valid_derived_p / valid_p, valid_p / valid_derived_p)
                
                # Convert rtol to fold change threshold (e.g., rtol=0.001 → 1.001x acceptable)
                fold_threshold = 1.0 + rtol
                
                # Check consistency based on fold change
                is_close[valid_mask] = fold_change[valid_mask] <= fold_threshold
                
                # For variants where fold change check fails, also verify with np.isclose
                # as a secondary check (handles edge cases with very small values)
                fold_inconsistent = valid_mask & ~is_close
                if fold_inconsistent.any():
                    np_close_secondary = np.isclose(
                        derived_p_values[fold_inconsistent], 
                        p_values[fold_inconsistent], 
                        rtol=rtol, 
                        atol=atol, 
                        equal_nan=equal_nan
                    )
                    is_close[fold_inconsistent] = np_close_secondary
            
            inconsistent = ~is_close & ~one_na  # Exclude single NaN cases from inconsistent count
            if inconsistent.any():
                inconsistent_num = inconsistent.sum()
                # Calculate fold change for reporting (only for valid pairs)
                inconsistent_valid = inconsistent & valid_mask
                if inconsistent_valid.any():
                    inconsistent_fold_change = fold_change[inconsistent_valid]
                    max_fold_idx = inconsistent_fold_change.idxmax()
                    max_fold_change = inconsistent_fold_change.max()
                    log.write("  -Potentially inconsistent (likely due to rounding): {} variant(s)".format(inconsistent_num),verbose=verbose)
                    log.write("  -Variant {} with max fold change: {} with {:.2f}x".format(id_to_use, sumstats.loc[max_fold_idx,id_to_use], max_fold_change),verbose=verbose)
                else:
                    # Fallback: report absolute difference for edge cases
                    diff = (derived_p_values - p_values).abs()
                    max_diff_idx = diff[inconsistent].idxmax()
                    max_diff = diff[inconsistent].max()
                    log.write("  -Potentially inconsistent (likely due to rounding): {} variant(s)".format(inconsistent_num),verbose=verbose)
                    log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[max_diff_idx,id_to_use], max_diff),verbose=verbose)
            else:
                log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
            check_status=1
    
    if has_mlog10p and has_p:
        log.write(" -Checking if MLOG10P-derived-P is consistent with P...",verbose=verbose)
        mlog10p_derived_p = _convert_mlog10p_to_p(sumstats["MLOG10P"])
        
        # For P-values, use fold change as the primary consistency metric
        # since P-values span many orders of magnitude (e.g., 1e-300 to 1)
        p_values = sumstats["P"].copy()
        # Ensure derived_p_values is a pandas Series with same index
        if isinstance(mlog10p_derived_p, pd.Series):
            derived_p_values = mlog10p_derived_p.copy()
        else:
            derived_p_values = pd.Series(mlog10p_derived_p, index=sumstats.index)
        
        # Handle NaN values - treat as consistent if both are NaN
        both_na = p_values.isna() & derived_p_values.isna()
        one_na = p_values.isna() | derived_p_values.isna()
        valid_mask = ~one_na & (p_values > 0) & (derived_p_values > 0)
        
        # Initialize consistency check: NaN pairs are considered consistent
        is_close = pd.Series(both_na, index=sumstats.index)
        
        # Initialize fold_change Series for all indices
        fold_change = pd.Series(index=sumstats.index, dtype=float)
        
        if valid_mask.any():
            # Calculate fold change: max(derived_p/p, p/derived_p) for symmetric measure
            # 1.0 = identical, >1.0 = different
            valid_p = p_values[valid_mask]
            valid_derived_p = derived_p_values[valid_mask]
            fold_change[valid_mask] = np.maximum(valid_derived_p / valid_p, valid_p / valid_derived_p)
            
            # Convert rtol to fold change threshold (e.g., rtol=0.001 → 1.001x acceptable)
            fold_threshold = 1.0 + rtol
            
            # Check consistency based on fold change
            is_close[valid_mask] = fold_change[valid_mask] <= fold_threshold
            
            # For variants where fold change check fails, also verify with np.isclose
            # as a secondary check (handles edge cases with very small values)
            fold_inconsistent = valid_mask & ~is_close
            if fold_inconsistent.any():
                np_close_secondary = np.isclose(
                    derived_p_values[fold_inconsistent], 
                    p_values[fold_inconsistent], 
                    rtol=rtol, 
                    atol=atol, 
                    equal_nan=equal_nan
                )
                is_close[fold_inconsistent] = np_close_secondary
        
        inconsistent = ~is_close & ~one_na  # Exclude single NaN cases from inconsistent count
        if inconsistent.any():
            inconsistent_num = inconsistent.sum()
            # Calculate fold change for reporting (only for valid pairs)
            inconsistent_valid = inconsistent & valid_mask
            if inconsistent_valid.any():
                inconsistent_fold_change = fold_change[inconsistent_valid]
                max_fold_idx = inconsistent_fold_change.idxmax()
                max_fold_change = inconsistent_fold_change.max()
                log.write("  -Potentially inconsistent (likely due to rounding): {} variant(s)".format(inconsistent_num),verbose=verbose)
                log.write("  -Variant {} with max fold change: {} with {:.2f}x".format(id_to_use, sumstats.loc[max_fold_idx,id_to_use], max_fold_change),verbose=verbose)
            else:
                # Fallback: report absolute difference for edge cases
                diff = (derived_p_values - p_values).abs()
                max_diff_idx = diff[inconsistent].idxmax()
                max_diff = diff[inconsistent].max()
                log.write("  -Potentially inconsistent (likely due to rounding): {} variant(s)".format(inconsistent_num),verbose=verbose)
                log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[max_diff_idx,id_to_use], max_diff),verbose=verbose)
        else:
            log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
        check_status=1

    if has_n and has_n_control and has_n_case:
        log.write(" -Checking if N is consistent with N_CASE + N_CONTROL ...", verbose=verbose) 
        n_expected = sumstats["N_CASE"] + sumstats["N_CONTROL"]
        is_close = sumstats["N"] == n_expected
        # Performance optimization: Use .any() instead of sum() > 0, and only compute diff when needed
        inconsistent = ~is_close
        if inconsistent.any():
            inconsistent_num = inconsistent.sum()
            # Performance optimization: Only compute diff when inconsistencies are found
            diff = abs(sumstats["N"] - n_expected)
            log.write("  -Potentially inconsistent: {} variant(s)".format(inconsistent_num),verbose=verbose)
            log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
        else:
            log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
        check_status=1
        
    if check_status==1:
        log.write(" -Note: Minor differences are typically due to rounding.",verbose=verbose)
        log.write("  -If the max difference is greater than expected, please check your original sumstats.",verbose=verbose)
    else:
        log.write(" -No availalbe columns for data consistency checking...Skipping...",verbose=verbose)
