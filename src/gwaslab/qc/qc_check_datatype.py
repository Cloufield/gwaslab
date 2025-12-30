from typing import TYPE_CHECKING, Optional, List, Any
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_reserved_headers import dtype_dict

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

# pandas.api.types.is_int64_dtype
# pandas.api.types.is_categorical_dtype

def _get_preferred_dtype(header: str) -> Optional[str]:
    """
    Get the preferred dtype for a column.
    
    Parameters
    ----------
    header : str
        Column name.
    
    Returns
    -------
    str
        Preferred dtype string.
    """
    if header in dtype_dict:
        return dtype_dict[header][0]
    return None

def check_datatype(sumstats: pd.DataFrame, verbose: bool = True, log: Log = Log()) -> None:
    """
    Inspect and report data types for all columns in a summary statistics
    DataFrame, comparing each against expected dtypes.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table to check.
    verbose : bool, default True
        Whether to print log messages to stdout.
    log : gwaslab.g_Log.Log
        Logger used for structured messages and warnings.

    Behavior
    --------
    - Logs aligned lists of column names, their dtypes, and verification flags
      ("T" for match, "F" for mismatch, "NA" for unknown header).
    - Emits targeted fix suggestions for CHR, POS, allele, and ID columns when
      their dtypes are incompatible.
    - For statistics columns present (e.g., BETA, SE, Z, P, OR, HR, CHISQ,
      INFO), warns to use `Sumstats.check_sanity()` only when those columns have
      incompatible dtypes.
    """
    
    try:
        headers = []
        dtypes = []
        verified = []
        raw_verified =[]

        for header,dtype in sumstats.dtypes.items():
            width = max(len(header),len(str(dtype)))
            
            header_fix_length = header + " "*(width- len(header) )
            
            dtype_fix_length  = str(dtype) + " "*(width- len(str(dtype)))
            
            verified_str = verify_datatype(header, dtype)
            verified_fix_length  = verified_str + " " *(width- len(verified_str))
            
            headers.append(format(header_fix_length))
            dtypes.append((str(dtype_fix_length)))
            verified.append(verified_fix_length)
            if verified_str == "F":
                raw_verified.append(header)

        log.write(" -Column  :", " ".join(headers), verbose=verbose)
        log.write(" -DType   :", " ".join(dtypes), verbose=verbose)
        log.write(" -Verified:", " ".join(verified), verbose=verbose)

        if len(raw_verified)>0:
            log.warning("Columns with possibly incompatible dtypes: {}".format(",".join(raw_verified)), verbose=verbose)
            try:
                if "CHR" in raw_verified:
                    log.warning("Consider using Sumstats.fix_chr() to fix CHR dtype", verbose=verbose)
                if "POS" in raw_verified:
                    log.warning("Consider using Sumstats.fix_pos() to fix POS dtype", verbose=verbose)
                if ("EA" in raw_verified) or ("NEA" in raw_verified) or ("REF" in raw_verified) or ("ALT" in raw_verified):
                    log.warning("Consider using Sumstats.fix_allele() to fix allele dtypes", verbose=verbose)
                if ("SNPID" in raw_verified) or ("rsID" in raw_verified):
                    log.warning("Consider using Sumstats.fix_id() to fix ID dtypes", verbose=verbose)
                stats_candidates = [
                    "BETA","SE","Z","P","MLOG10P","OR","OR_95L","OR_95U",
                    "HR","HR_95L","HR_95U","CHISQ","F","T","SNPR2","INFO"
                ]
                present_stats = [c for c in stats_candidates if c in sumstats.columns]
                bad_stats = [c for c in present_stats if verify_datatype(c, sumstats.dtypes[c]) == "F"]
                if len(bad_stats) > 0:
                    for c in bad_stats:
                        log.warning("Consider using Sumstats.check_sanity() to check {} statistics".format(c), verbose=verbose)
            except:
                pass
    except:
        pass

def verify_datatype(header: str, dtype: Any) -> str:
    """
    Verify a single column's dtype against expected dtypes.

    Parameters
    ----------
    header : str
        Column name to check.
    dtype : object
        The pandas dtype object or its string representation.

    Returns
    -------
    str
        "T" if dtype matches one of the expected entries, "F" if it does not,
        and "NA" if the column is not recognized.
    """

    if header in dtype_dict.keys():
        if str(dtype) in dtype_dict[header]:
            return "T"
        else:
            return "F"
    else:
        return "NA"

def quick_convert_datatype(sumstats: pd.DataFrame, log: Log, verbose: bool) -> pd.DataFrame:
    """
    Attempt to convert columns with incompatible dtypes to the preferred dtype
    for that column. Logs successes and failures.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table to convert.
    log : gwaslab.g_Log.Log
        Logger used for messages.
    verbose : bool
        Whether to print log messages to stdout.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns converted where possible.
    """
    for col in sumstats.columns:
        if col in dtype_dict.keys():
            current_dtype = str(sumstats[col].dtype)
            # Check if conversion is needed
            if current_dtype not in dtype_dict[col]:
                # Current dtype is not acceptable, convert to preferred dtype
                datatype = _get_preferred_dtype(col)
                if datatype is None:
                    datatype = dtype_dict[col][0]  # Fallback to first option
                
                log.write(" -Trying to convert datatype for {}: {} -> {}...".format(col, current_dtype, datatype), end="" ,verbose=verbose)
                try:
                    sumstats[col] = sumstats[col].astype(datatype)
                    log.write("Success",show_time=False, verbose=verbose)
                except:
                    log.write("Failed...",show_time=False,verbose=verbose)
                    pass
    return sumstats

def check_dataframe_shape(
    sumstats: pd.DataFrame,
    log: Log,
    verbose: bool,
    sumstats_obj: Optional['Sumstats'] = None
) -> None:
    """
    Log DataFrame shape and estimated memory usage in megabytes.
    Only logs if shape has changed from the last check.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table.
    log : gwaslab.g_Log.Log
        Logger used for messages.
    verbose : bool
        Whether to print log messages to stdout.
    sumstats_obj : Sumstats, optional
        Sumstats object instance. If provided, checks _last_shape attribute
        to skip logging if shape is unchanged. If None, will try to get from log._sumstats_obj.
    """
    # If sumstats_obj not provided, try to get it from log object
    if sumstats_obj is None:
        sumstats_obj = getattr(log, '_sumstats_obj', None)
    # Use the standardized logging method which handles change detection
    log.log_dataframe_shape(sumstats, verbose=verbose, sumstats_obj=sumstats_obj)
    
def check_dataframe_memory_usage(sumstats: pd.DataFrame, log: Log, verbose: bool) -> None:
    """
    Log DataFrame memory usage in megabytes.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table.
    log : gwaslab.g_Log.Log
        Logger used for messages.
    verbose : bool
        Whether to print log messages to stdout.
    """
    memory_in_mb = sumstats.memory_usage().sum()/1024/1024
    try:
        log.write(" -Current Dataframe memory usage: {:.2f} MB".format(memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Memory usage...")

def check_datatype_for_cols(
    sumstats_obj: 'Sumstats',
    cols: Optional[List[str]] = None,
    verbose: bool = True,
    log: Log = Log(),
    fix: bool = False,
    **fix_kwargs: Any
) -> None:
    """
    Verify dtypes for a specified subset of columns and emit fix suggestions.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to check.
    cols : list[str] or None, default None
        Column names to verify. If None, no columns are checked.
    verbose : bool, default True
        Whether to print log messages to stdout.
    log : gwaslab.g_Log.Log
        Logger used for messages and warnings.
    fix : bool, default False
        If True, attempt to fix incompatible dtypes automatically.
    **fix_kwargs : dict
        Additional keyword arguments to pass to fix functions.

    Behavior
    --------
    - Collects failing columns and logs targeted fix suggestions for CHR, POS,
      allele, and ID columns.
    - For present statistics columns among the specified `cols`, warns to use
      `Sumstats.check_sanity()` only when those columns have incompatible
      dtypes.
    - Raises `ValueError` listing failing columns to encourage corrective
      actions via Sumstats methods.
    """
    sumstats = sumstats_obj.data
    if cols is None:
        cols = []
    try:
        failed = []

        for header,dtype in sumstats.dtypes.items():
            if header in cols:
                verified_str = verify_datatype(header, dtype)
                if verified_str=="F":
                    failed.append(header)
    except:
        pass

    if len(failed) >0:
        log.warning("Data types were not fixed for : {} ".format(",".join(failed)))
        try:
            if fix is True:
                try:
                    from gwaslab.qc.qc_fix_sumstats import _fix_chr, _fix_pos, _fix_allele, _fix_ID
                    from gwaslab.info.g_meta import _update_qc_step
                except Exception:
                    _fix_chr = None; _fix_pos = None; _fix_allele = None; _fix_ID = None
                    _update_qc_step = None

                if "CHR" in failed and _fix_chr is not None:
                    chr_kwargs = {k:v for k,v in fix_kwargs.items() if k in ["remove","add_prefix"]}
                    sumstats_obj.data = _fix_chr(sumstats_obj, log=log, verbose=verbose, **chr_kwargs)
                    if _update_qc_step is not None:
                        _update_qc_step(sumstats_obj, "chr", chr_kwargs, True)
                if "POS" in failed and _fix_pos is not None:
                    pos_kwargs = {k:v for k,v in fix_kwargs.items() if k in ["remove","lower_limit","upper_limit","limit"]}
                    sumstats_obj.data = _fix_pos(sumstats_obj, log=log, verbose=verbose, **pos_kwargs)
                    if _update_qc_step is not None:
                        _update_qc_step(sumstats_obj, "pos", pos_kwargs, True)
                if (set(["EA","NEA","REF","ALT"]) & set(failed)) and _fix_allele is not None:
                    allele_kwargs = {k:v for k,v in fix_kwargs.items() if k in ["remove"]}
                    sumstats_obj.data = _fix_allele(sumstats_obj, log=log, verbose=verbose, **allele_kwargs)
                    if _update_qc_step is not None:
                        _update_qc_step(sumstats_obj, "allele", allele_kwargs, True)
                if ("SNPID" in failed or "rsID" in failed) and _fix_ID is not None:
                    id_kwargs = {k:v for k,v in fix_kwargs.items() if k in ["fixprefix","fixchrpos","fixid","fixeanea","fixeanea_flip","fixsep","reversea","overwrite","forcefixid"]}
                    sumstats_obj.data = _fix_ID(sumstats_obj, log=log, verbose=verbose, **id_kwargs)
                    if _update_qc_step is not None:
                        _update_qc_step(sumstats_obj, "id", id_kwargs, True)
                
                # Update sumstats reference after fixes
                sumstats = sumstats_obj.data

                # Re-verify targeted columns after attempted fixes
                re_failed = []
                for header,dtype in sumstats.dtypes.items():
                    if header in cols:
                        verified_str = verify_datatype(header, dtype)
                        if verified_str=="F":
                            re_failed.append(header)

                if len(re_failed) == 0:
                    log.write(" -Dtype fixes applied successfully for requested columns.", verbose=verbose)
                    return sumstats_obj.data
                else:
                    log.warning("Some columns still have incompatible dtypes after fixes: {}".format(",".join(re_failed)), verbose=verbose)

            # Suggest manual next steps if fix was False or fixes failed
            if "CHR" in failed:
                log.warning("Consider using Sumstats.fix_chr() to fix CHR dtype", verbose=verbose)
            if "POS" in failed:
                log.warning("Consider using Sumstats.fix_pos() to fix POS dtype", verbose=verbose)
            if ("EA" in failed) or ("NEA" in failed) or ("REF" in failed) or ("ALT" in failed):
                log.warning("Consider using Sumstats.fix_allele() to fix allele dtypes", verbose=verbose)
            if ("SNPID" in failed) or ("rsID" in failed):
                log.warning("Consider using Sumstats.fix_id() to fix ID dtypes", verbose=verbose)

            stats_candidates = [
                "BETA","SE","Z","P","MLOG10P","OR","OR_95L","OR_95U",
                "HR","HR_95L","HR_95U","CHISQ","F","T","SNPR2","INFO"
            ]
            present_stats = [c for c in stats_candidates if c in sumstats.columns]
            bad_stats = [c for c in present_stats if verify_datatype(c, sumstats.dtypes[c]) == "F"]
            if len(bad_stats) > 0:
                for c in bad_stats:
                    log.warning("Consider using Sumstats.check_sanity() to check {} statistics".format(c), verbose=verbose)
        except:
            pass

        if fix is True:
            raise ValueError("Failed to fix dtypes for requested columns: {}".format(",".join(failed)))
        else:
            raise ValueError("Please fix dtypes using Sumstats methods (fix_chr, fix_pos, fix_allele, fix_id) for: {}. If statistics columns are present, consider Sumstats.check_sanity().".format(",".join(failed)))
