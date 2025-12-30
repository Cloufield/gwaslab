from typing import TYPE_CHECKING, Union, Optional, List, Tuple, Any
import pandas as pd
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
from scipy import stats
from gwaslab.info.g_Log import Log
import gc
#from gwaslab.qc_fix_sumstats import sortcolumn
from gwaslab.info.g_version import _get_version
from gwaslab.qc.qc_check_datatype import check_datatype
from scipy.special import erfcinv

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _fill_data( 
    insumstats: Union['Sumstats', pd.DataFrame],
    to_fill: Optional[Union[str, List[str]]] = None,
    df: Optional[pd.DataFrame] = None,
    overwrite: bool = False,
    verbose: bool = True,
    only_sig: bool = False,
    sig_level: float = 5e-8,
    extreme: bool = False,
    log: Log = Log()
    ) -> pd.DataFrame:
    """
    Fill missing statistical values in genetic summary statistics from available columns.
    
    This function systematically derives missing statistical values using relationships
    between different statistical measures (e.g., converting beta/SE to Z-scores, or ORs
    to betas). It handles multiple conversion pathways and maintains data consistency.
    
    Parameters
    ----------
    insumstats : pandas.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object.
    to_fill : str or list of str
        Column name(s) to fill. Common values include "OR","OR_95L","OR_95U","BETA","SE","P","Z","CHISQ","MLOG10P","MAF", etc.
    overwrite : bool, optional, default=False
        If True, overwrite existing values in target columns.
    verbose : bool, optional, default=True
        Whether to display progress messages.
    extreme : bool, optional, default=False
        If True, use extreme value calculations for -log10(P). Helpful when P<1e-300 (float64 datatype limits).
    
    Returns
    -------
    pandas.DataFrame
        Modified summary statistics DataFrame with filled values

    Less used parameters
    ----------------
    df : str, optional
        Column name containing degrees of freedom for chi-square tests. Only used when CHISQ 
    only_sig : bool, optional, default=False
        If True, only update values for significant variants.
    sig_level : float, optional, default=5e-8
        Significance threshold for P-value filtering.
    """
    # Handle both DataFrame and Sumstats object
    import pandas as pd
    if isinstance(insumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats_obj = insumstats
        insumstats = sumstats_obj.data
        is_dataframe = False
    
    # Normalize input
    if type(to_fill) is str:
        to_fill = [to_fill]
    if to_fill is None:
        to_fill = []
    
    sumstats = insumstats.copy()
    log.write("Start filling data using existing columns...{}".format(_get_version()), verbose=verbose)
    check_datatype(sumstats, verbose=verbose, log=log)
    
    # Filter out columns that already exist (unless overwrite=True)
    valid_targets = {"OR", "OR_95L", "OR_95U", "BETA", "SE", "P", "Z", "CHISQ", "MLOG10P", "MAF", "SIG"}
    to_fill = [col for col in to_fill if col in valid_targets]
    
    if not to_fill:
        log.write(" -No valid columns to fill.", verbose=verbose)
        log.write("Finished filling data using existing columns.", verbose=verbose)
        return sumstats
    
    if not overwrite:
        existing_cols = [col for col in to_fill if col in sumstats.columns]
        to_fill = [col for col in to_fill if col not in sumstats.columns]
        if existing_cols:
            log.write(" -Skipping existing columns: {}".format(existing_cols), verbose=verbose)
    
    if not to_fill:
        log.write(" -All target columns already exist.", verbose=verbose)
        log.write("Finished filling data using existing columns.", verbose=verbose)
        return sumstats
    
    # Perform iterative filling
    log.write(" -Target columns: {}".format(to_fill), verbose=verbose)
    fill_iteratively(sumstats, to_fill, log, only_sig, df, extreme, verbose, sig_level)
    
    # Final summary
    still_missing = [col for col in to_fill if col not in sumstats.columns]
    if still_missing:
        log.write(" -Warning: Could not fill: {}".format(still_missing), verbose=verbose)
    else:
        log.write(" -Successfully filled all requested columns.", verbose=verbose)
    
    gc.collect()
    log.write("Finished filling data using existing columns.", verbose=verbose)
    
    # Update only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        return sumstats_obj.data
    else:
        return sumstats
    
##########################################################################################################################    
    
def fill_p(sumstats: pd.DataFrame, log: Log, df: Optional[str] = None, only_sig: bool = False, sig_level: float = 5e-8, overwrite: bool = False, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill P column from MLOG10P, Z, or CHISQ."""
    if "MLOG10P" in sumstats.columns:
        log.write("    Filling P from MLOG10P...", verbose=verbose)
        sumstats["P"] = np.power(10, -sumstats["MLOG10P"])
        return 1, filled_count + 1
    
    if "Z" in sumstats.columns:
        log.write("    Filling P from Z...", verbose=verbose)
        stats.chisqprob = lambda chisq, degree_of_freedom: stats.chi2.sf(chisq, degree_of_freedom)
        sumstats["P"] = ss.chisqprob(sumstats["Z"]**2, 1)
        return 1, filled_count + 1
    
    if "CHISQ" in sumstats.columns:
        log.write("    Filling P from CHISQ...", verbose=verbose)
        stats.chisqprob = lambda chisq, degree_of_freedom: stats.chi2.sf(chisq, degree_of_freedom)
        if df is None:
            if only_sig and overwrite:
                sumstats.loc[sumstats["P"] < sig_level, "P"] = stats.chisqprob(
                    sumstats.loc[sumstats["P"] < sig_level, "CHISQ"], 1)
            else:
                sumstats["P"] = stats.chisqprob(sumstats["CHISQ"], 1)
        else:
            if only_sig and overwrite:
                sumstats.loc[sumstats["P"] < sig_level, "P"] = stats.chisqprob(
                    sumstats.loc[sumstats["P"] < sig_level, "CHISQ"],
                    sumstats.loc[sumstats["P"] < sig_level, df].astype("int"))
            else:
                sumstats["P"] = stats.chisqprob(sumstats["CHISQ"], sumstats[df].astype("int"))
        return 1, filled_count + 1
    
    return 0, filled_count

def fill_z(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill Z column from BETA/SE."""
    if "BETA" in sumstats.columns and "SE" in sumstats.columns:
        log.write("    Filling Z from BETA/SE...", verbose=verbose)
        sumstats["Z"] = sumstats["BETA"] / sumstats["SE"]
        return 1, filled_count + 1
    return 0, filled_count

def fill_chisq(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill CHISQ column from Z or P."""
    if "Z" in sumstats.columns:
        log.write("    Filling CHISQ from Z...", verbose=verbose)
        sumstats["CHISQ"] = sumstats["Z"]**2
        return 1, filled_count + 1
    if "P" in sumstats.columns:
        log.write("    Filling CHISQ from P...", verbose=verbose)
        sumstats["CHISQ"] = ss.chi2.isf(sumstats["P"], 1)
        return 1, filled_count + 1
    return 0, filled_count

def fill_or(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill OR and optionally OR_95L/OR_95U from BETA/SE."""
    if "BETA" in sumstats.columns:
        log.write("    Filling OR from BETA...", verbose=verbose)
        sumstats["OR"] = np.exp(sumstats["BETA"])
        filled_count += 1
        
        if "SE" in sumstats.columns:
            log.write("    Filling OR_95L/OR_95U from BETA/SE...", verbose=verbose)
            z_critical = ss.norm.ppf(0.975)
            sumstats["OR_95L"] = np.exp(sumstats["BETA"] - z_critical * sumstats["SE"])
            sumstats["OR_95U"] = np.exp(sumstats["BETA"] + z_critical * sumstats["SE"])
            filled_count += 1
        return 1, filled_count
    return 0, filled_count
def fill_beta(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill BETA column from OR."""
    if "OR" in sumstats.columns:
        log.write("    Filling BETA from OR...", verbose=verbose)
        sumstats["BETA"] = np.log(sumstats["OR"])
        return 1, filled_count + 1
    return 0, filled_count

def fill_se(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill SE column from BETA/P, OR/OR_95U, or OR/OR_95L."""
    if "BETA" in sumstats.columns and "P" in sumstats.columns:
        log.write("    Filling SE from BETA/P...", verbose=verbose)
        abs_z = np.sqrt(2) * erfcinv(sumstats["P"])
        sumstats["SE"] = np.abs(sumstats["BETA"] / abs_z)
        return 1, filled_count + 1
    
    if "OR" in sumstats.columns and "OR_95U" in sumstats.columns:
        log.write("    Filling SE from OR/OR_95U...", verbose=verbose)
        z_critical = ss.norm.ppf(0.975)
        sumstats["SE"] = (np.log(sumstats["OR_95U"]) - np.log(sumstats["OR"])) / z_critical
        return 1, filled_count + 1
    
    if "OR" in sumstats.columns and "OR_95L" in sumstats.columns:
        log.write("    Filling SE from OR/OR_95L...", verbose=verbose)
        z_critical = ss.norm.ppf(0.975)
        sumstats["SE"] = (np.log(sumstats["OR"]) - np.log(sumstats["OR_95L"])) / z_critical
        return 1, filled_count + 1
    
    return 0, filled_count

def fill_mlog10p(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill MLOG10P column from P."""
    if "P" in sumstats.columns:
        log.write("    Filling MLOG10P from P...", verbose=verbose)
        sumstats["MLOG10P"] = -np.log10(sumstats["P"])
        return 1, filled_count + 1
    return 0, filled_count

def fill_extreme_mlog10p(sumstats: pd.DataFrame, df: Optional[str], log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill MLOG10P using extreme value methods (for very small P-values).
    
    Tries extreme value methods first, falls back to P -> MLOG10P conversion as last resort.
    """
    # Try extreme value methods first
    if "Z" in sumstats.columns:
        log.write("    Filling MLOG10P from Z (extreme)...", verbose=verbose)
        sumstats = fill_extreme_mlog10(sumstats, "Z")
        return 1, filled_count + 1
    
    if "BETA" in sumstats.columns and "SE" in sumstats.columns:
        log.write("    Filling MLOG10P from BETA/SE (extreme, via Z)...", verbose=verbose)
        sumstats["Z"] = sumstats["BETA"] / sumstats["SE"]
        sumstats = fill_extreme_mlog10(sumstats, "Z")
        return 1, filled_count + 1
    
    if "CHISQ" in sumstats.columns and "DOF" in sumstats.columns:
        log.write("    Filling MLOG10P from CHISQ/DOF (extreme)...", verbose=verbose)
        sumstats = fill_extreme_mlog10_chisq(sumstats, "CHISQ", df)
        return 1, filled_count + 1
    
    # Last resort: convert from P
    if "P" in sumstats.columns:
        log.write("    Filling MLOG10P from P (fallback)...", verbose=verbose)
        sumstats["MLOG10P"] = -np.log10(sumstats["P"])
        return 1, filled_count + 1
    
    return 0, filled_count

def fill_maf(sumstats: pd.DataFrame, log: Log, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill MAF column from EAF."""
    if "EAF" in sumstats.columns:
        log.write("    Filling MAF from EAF...", verbose=verbose)
        sumstats["MAF"] = sumstats["EAF"].apply(lambda x: min(x, 1-x) if pd.notnull(x) else np.nan)
        return 1, filled_count + 1
    return 0, filled_count

def fill_sig(sumstats: pd.DataFrame, log: Log, sig_level: float = 5e-8, verbose: bool = True, filled_count: int = 0) -> Tuple[int, int]:
    """Fill SIGNIFICANT column from P or MLOG10P."""
    if "P" in sumstats.columns:
        log.write("    Filling SIGNIFICANT from P (threshold={})...".format(sig_level), verbose=verbose)
        is_sig = sumstats["P"] < sig_level
    elif "MLOG10P" in sumstats.columns:
        log.write("    Filling SIGNIFICANT from MLOG10P (threshold={})...".format(sig_level), verbose=verbose)
        is_sig = sumstats["MLOG10P"] > np.log10(1/sig_level)
    else:
        return 0, filled_count
    
    sumstats["SIGNIFICANT"] = False
    sumstats.loc[is_sig, "SIGNIFICANT"] = True
    return 1, filled_count + 1

####################################################################################################################
def fill_extreme_mlog10(sumstats: pd.DataFrame, z: str) -> pd.DataFrame:
    log_pvalue = np.log(2) + ss.norm.logsf(np.abs(sumstats[z])) #two-sided
    log10_pvalue = log_pvalue/np.log(10)
    mantissa = 10**(log10_pvalue %1 )
    exponent = log10_pvalue // 1
    sumstats["MLOG10P"] = -log10_pvalue
    sumstats["P_MANTISSA"]= mantissa
    sumstats["P_EXPONENT"]= exponent
    return sumstats

def fill_extreme_mlog10_chisq(sumstats: pd.DataFrame, chisq: str, df: str) -> pd.DataFrame:
    #https://stackoverflow.com/a/46416222/199475
    log_pvalue = ss.chi2.logsf(sumstats[chisq], sumstats[df])

    log10_pvalue = log_pvalue/np.log(10)
    
    mantissa = 10**(log10_pvalue %1)
    exponent = log10_pvalue // 1
    sumstats["MLOG10P"] = -log10_pvalue
    sumstats["P_MANTISSA"]= mantissa
    sumstats["P_EXPONENT"]= exponent
    return sumstats

####################################################################################################################
def _try_fill_mlog10p(sumstats: pd.DataFrame, to_fill: List[str], df: Optional[str], log: Log, verbose: bool) -> Tuple[bool, List[str]]:
    """
    Attempt to fill MLOG10P column.
    
    Strategy:
    1. Try extreme value methods first (Z, BETA/SE, CHISQ/DOF)
    2. If that fails and P doesn't exist, try to fill P first, then retry
    3. Drop intermediate P column if it wasn't requested
    4. Always drop P_MANTISSA and P_EXPONENT (created by extreme methods)
    
    Returns:
        tuple: (success: bool, columns_to_remove: list)
    """
    columns_to_remove = []
    
    # Try extreme method first (includes P fallback)
    status, _ = fill_extreme_mlog10p(sumstats, df, log, verbose=verbose, filled_count=0)
    
    if status == 1:
        # Extreme method succeeded - always drop P_MANTISSA and P_EXPONENT
        # (they're created by extreme methods but not needed after MLOG10P is filled)
        columns_to_remove.extend(["P_MANTISSA", "P_EXPONENT"])
        return True, columns_to_remove
    
    # If failed and P doesn't exist, try to create P first
    if "P" not in sumstats.columns:
        fill_p(sumstats, log, verbose=verbose, filled_count=0)
        # Retry extreme method (will now use P as fallback)
        status, _ = fill_extreme_mlog10p(sumstats, df, log, verbose=verbose, filled_count=0)
        
        # If P was created as intermediate step and not requested, mark for removal
        if "P" not in to_fill and "P" in sumstats.columns:
            columns_to_remove.append("P")
        
        # Drop P_MANTISSA and P_EXPONENT if MLOG10P was filled successfully
        if status == 1:
            columns_to_remove.extend(["P_MANTISSA", "P_EXPONENT"])
        
        return status == 1, columns_to_remove
    
    return False, columns_to_remove


def fill_iteratively(sumstats: pd.DataFrame, to_fill: List[str], log: Log, only_sig: bool, df: Optional[str], extreme: bool, verbose: bool, sig_level: float) -> pd.DataFrame:
    """
    Iteratively fill columns using available conversion pathways.
    
    This function attempts to fill target columns in multiple rounds, as newly filled
    columns may enable further conversions in subsequent rounds.
    
    Args:
        sumstats: DataFrame to fill
        to_fill: Set or list of column names to fill
        log: Log object for messages
        only_sig: Only fill significant variants
        df: Degrees of freedom column name
        extreme: (deprecated, kept for compatibility)
        verbose: Verbose logging
        sig_level: Significance threshold
    """
    # Define fill functions for standard columns
    fill_functions = {
        "OR": lambda: fill_or(sumstats, log, verbose=verbose, filled_count=0),
        "BETA": lambda: fill_beta(sumstats, log, verbose=verbose, filled_count=0),
        "SE": lambda: fill_se(sumstats, log, verbose=verbose, filled_count=0),
        "P": lambda: fill_p(sumstats, log, only_sig=only_sig, df=df, sig_level=sig_level, verbose=verbose, filled_count=0),
        "Z": lambda: fill_z(sumstats, log, verbose=verbose, filled_count=0),
        "CHISQ": lambda: fill_chisq(sumstats, log, verbose=verbose, filled_count=0),
        "MAF": lambda: fill_maf(sumstats, log, verbose=verbose, filled_count=0),
        "SIG": lambda: fill_sig(sumstats, log, sig_level=sig_level, verbose=verbose, filled_count=0),
    }
    
    # Convert to set for efficient operations
    remaining = set(to_fill)
    max_rounds = len(remaining) + 1
    
    # Iterate through rounds until all columns are filled or no progress is made
    for round_num in range(1, max_rounds + 1):
        round_filled = []
        columns_to_remove = []
        
        # Try to fill each remaining column
        for col in list(remaining):
            success = False
            
            if col == "MLOG10P":
                # Special handling for MLOG10P (always uses extreme methods)
                success, cols_to_remove = _try_fill_mlog10p(sumstats, remaining, df, log, verbose)
                columns_to_remove.extend(cols_to_remove)
                    
            elif col in fill_functions:
                # Standard column filling
                status, _ = fill_functions[col]()
                success = (status == 1)
            
            # Update tracking if successful
            if success:
                remaining.remove(col)
                round_filled.append(col)
        
        # Remove intermediate columns that were created but not requested
        # Remove duplicates and only drop columns that actually exist
        columns_to_drop = [col for col in set(columns_to_remove) if col in sumstats.columns]
        if columns_to_drop:
            sumstats.drop(labels=columns_to_drop, axis=1, inplace=True)
        
        # Log round results
        _log_round_results(round_num, round_filled, remaining, log, verbose)
        
        # Check termination conditions
        if not remaining:
            # All columns filled successfully
            break
        if not round_filled:
            # No progress made this round - cannot fill remaining columns
            break


def _log_round_results(
    round_num: int, 
    round_filled: List[str], 
    remaining: List[str], 
    log: Log, 
    verbose: bool
) -> None:
    """Log the results of a filling round."""
    if round_filled:
        log.write("  [Round {}] Filled: {}".format(round_num, round_filled), verbose=verbose)
        if remaining:
            log.write("  [Round {}] Remaining: {}".format(round_num, sorted(remaining)), verbose=verbose)
        else:
            log.write("  [Round {}] All columns filled!".format(round_num), verbose=verbose)
    elif remaining:
        log.write("  [Round {}] No progress. Unable to fill: {}".format(round_num, sorted(remaining)), verbose=verbose)

         
###Base functions########################################################################################

def _convert_betase_to_z(beta: pd.Series, se: pd.Series) -> pd.Series:
    return beta/se 

def _convert_betase_to_p(beta: pd.Series, se: pd.Series) -> pd.Series:
    z = _convert_betase_to_z(beta, se)
    p = _convert_z_to_p(z)
    return p

def _convert_betase_to_mlog10p(beta: pd.Series, se: pd.Series) -> pd.Series:
    z = _convert_betase_to_z(beta, se)
    mlog10p = _convert_z_to_mlog10p(z)
    return mlog10p

def _convert_p_to_chisq(p: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return ss.chi2.isf(p, 1)

def _convert_z_to_chisq(z: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return (z)**2

def _convert_z_to_p(z: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return ss.chi2.sf(z**2,1) 

def _convert_z_to_mlog10p(z: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    log_pvalue = np.log(2) + ss.norm.logsf(np.abs(z)) #two-sided
    mlog10p = log_pvalue/np.log(10)
    return -mlog10p 

def _conver_chisq_to_p(chisq: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return ss.chi2.sf(chisq,1)

def _convert_mlog10p_to_p(mlog10p: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return np.power(10, -mlog10p) 

def _convert_or_to_beta(OR: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return np.log(OR)  

def _convert_beta_to_or(beta: Union[pd.Series, np.ndarray, float]) -> Union[pd.Series, np.ndarray, float]:
    return np.exp(beta)   

def rank_based_int(series: pd.Series, c: float = 3/8) -> pd.Series:
    """
    Perform rank-based inverse normal transformation.
    
    Parameters
    ----------
    series : pd.Series
        Input series to transform
    c : float, default 3/8
        Blom's constant for rank transformation
        
    Returns
    -------
    pd.Series
        Transformed series with normal distribution
    """
    #https://onlinelibrary.wiley.com/doi/10.1111/biom.13214
    n=sum(~series.isna())
    normalized_value = norm.ppf((series.rank()-c)/(n+1-2*c))
    return normalized_value


################################################################################################################################################################################

def _get_multi_min(
    sumstats_multi: pd.DataFrame, 
    col: str, 
    nstudy: int, 
    log: Log = Log(), 
    verbose: bool = True
) -> pd.DataFrame:
    cols =[]
    for i in range(nstudy):
        single_header = "{}_{}".format(col, i + 1)
        if single_header in sumstats_multi.columns:
            cols.append(single_header)

    combined_header = "{}_{}".format(col, "MIN")
    log.write("  -Filling {} using {}".format(combined_header,",".join(cols)), verbose=verbose)
    sumstats_multi[combined_header] = sumstats_multi[cols].min(axis=1)
    
    combined_header_index = "{}_{}_COL".format(col, "MIN")
    sumstats_multi[combined_header_index] = sumstats_multi[cols].idxmin(axis=1)
    return sumstats_multi
