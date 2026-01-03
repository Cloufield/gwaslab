from typing import Any, Optional, List
import gc
import pandas as pd
import polars as pl
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_reserved_headers import get_dtype_dict_polars
# pandas.api.types.is_int64_dtype
# pandas.api.types.is_categorical_dtype

# Get dtype_dict from reserved headers JSON (single source of truth)
dtype_dict = get_dtype_dict_polars()

def check_datatype_polars(sumstats: pl.DataFrame, verbose: bool = True, log: Log = Log()) -> None:
    
    #try:
    headers = []
    dtypes = []
    verified = []
    raw_verified =[]
    for header,dtype in  sumstats.schema.items():
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
                log.warning("Consider using Sumstatsp.fix_chr() to fix CHR dtype", verbose=verbose)
            if "POS" in raw_verified:
                log.warning("Consider using Sumstatsp.fix_pos() to fix POS dtype", verbose=verbose)
            if ("EA" in raw_verified) or ("NEA" in raw_verified) or ("REF" in raw_verified) or ("ALT" in raw_verified):
                log.warning("Consider using Sumstatsp.fix_allele() to fix allele dtypes", verbose=verbose)
            if ("SNPID" in raw_verified) or ("rsID" in raw_verified):
                log.warning("Consider using Sumstatsp.fix_id() to fix ID dtypes", verbose=verbose)
            stats_candidates = [
                "BETA","SE","Z","P","MLOG10P","OR","OR_95L","OR_95U",
                "HR","HR_95L","HR_95U","CHISQ","F","T","SNPR2","INFO"
            ]
            present_stats = [c for c in stats_candidates if c in sumstats.columns]
            bad_stats = [c for c in present_stats if verify_datatype(c, sumstats.schema[c]) == "F"]
            if len(bad_stats) > 0:
                for c in bad_stats:
                    log.warning("Consider using Sumstatsp.check_sanity() to check {} statistics".format(c), verbose=verbose)
        except:
            pass
    #except:
    #    pass

def verify_datatype(header: str, dtype: Any) -> str:
    """
    Verify a single column's dtype against expected dtypes for polars.
    
    Parameters
    ----------
    header : str
        Column name to check.
    dtype : Any
        The polars dtype object.
    
    Returns
    -------
    str
        "T" if dtype matches one of the expected entries, "F" if it does not,
        and "NA" if the column is not recognized.
    """
    if header in dtype_dict.keys():
        # Try direct comparison first (for polars type objects)
        if dtype in dtype_dict[header]:
            return "T"
        # Also try string comparison as fallback
        dtype_str = str(dtype)
        if any(str(expected_dtype) == dtype_str for expected_dtype in dtype_dict[header]):
            return "T"
        return "F"
    else:
        return "NA"

def quick_convert_datatype(sumstats: pl.DataFrame, log: Log, verbose: bool) -> pl.DataFrame:
    for col in sumstats.columns:
        if col in dtype_dict.keys():
            if sumstats[col].dtype not in dtype_dict[col]:
                datatype=dtype_dict[col][0]
                log.write(" -Trying to convert datatype for {}: {} -> {}...".format(col, str(sumstats[col].dtype), datatype), end="" ,verbose=verbose)
                try:
                    sumstats = sumstats.cast({col: datatype})
                    log.write("Success",show_time=False, verbose=verbose)
                except:
                    log.write("Failed...",show_time=False,verbose=verbose)
                    pass
    return sumstats

def check_dataframe_shape_polars(sumstats: pl.DataFrame, log: Log, verbose: bool) -> None:
    memory_in_mb = sumstats.estimated_size(unit="mb") 
    try:
        log.write(" -Current Dataframe shape : {} x {} ; Memory usage: {:.2f} MB".format(len(sumstats),len(sumstats.columns),memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Dataframe shape...")
    
def check_dataframe_memory_usage(sumstats: pl.DataFrame, log: Log, verbose: bool) -> None:
    memory_in_mb = sumstats.estimated_size(unit="mb") 
    try:
        log.write(" -Current Dataframe memory usage: {:.2f} MB".format(memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Memory usage...")

def check_datatype_for_cols_polars(
    sumstats_obj: Any,
    cols: Optional[List[str]] = None,
    verbose: bool = True,
    log: Log = Log(),
    fix: bool = False,
    **fix_kwargs: Any
) -> Optional[pl.DataFrame]:
    """
    Verify dtypes for a specified subset of columns and optionally fix them.

    Parameters
    ----------
    sumstats_obj : Sumstatsp or pl.DataFrame
        Sumstats object or polars DataFrame containing the data to check.
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

    Returns
    -------
    pl.DataFrame or None
        Returns the DataFrame if fix=True and fixes were successful, None otherwise.
    """
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pl.DataFrame):
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    if cols is None:
        cols = []
    
    try:
        failed = []

        for header, dtype in sumstats.schema.items():
            if header in cols:
                verified_str = verify_datatype(header, dtype)
                if verified_str == "F":
                    failed.append(header)
    except:
        pass

    if len(failed) > 0:
        log.warning("Data types were not fixed for : {} ".format(",".join(failed)), verbose=verbose)
        try:
            if fix is True:
                try:
                    from gwaslab.qc.qc_fix_sumstats_polars import _fix_chrp, _fix_posp
                    from gwaslab.info.g_meta import _update_qc_step
                except Exception:
                    _fix_chrp = None
                    _fix_posp = None
                    _update_qc_step = None

                if "CHR" in failed and _fix_chrp is not None:
                    chr_kwargs = {k: v for k, v in fix_kwargs.items() if k in ["remove", "add_prefix"]}
                    if is_dataframe:
                        sumstats = _fix_chrp(sumstats, log=log, verbose=verbose, **chr_kwargs)
                    else:
                        sumstats_obj.data = _fix_chrp(sumstats_obj, log=log, verbose=verbose, **chr_kwargs)
                        sumstats = sumstats_obj.data
                        if _update_qc_step is not None:
                            _update_qc_step(sumstats_obj, "chr", chr_kwargs, True)
                
                if "POS" in failed and _fix_posp is not None:
                    pos_kwargs = {k: v for k, v in fix_kwargs.items() if k in ["remove", "lower_limit", "upper_limit", "limit"]}
                    if is_dataframe:
                        sumstats = _fix_posp(sumstats, log=log, verbose=verbose, **pos_kwargs)
                    else:
                        sumstats_obj.data = _fix_posp(sumstats_obj, log=log, verbose=verbose, **pos_kwargs)
                        sumstats = sumstats_obj.data
                        if _update_qc_step is not None:
                            _update_qc_step(sumstats_obj, "pos", pos_kwargs, True)

                # Re-verify targeted columns after attempted fixes
                re_failed = []
                for header, dtype in sumstats.schema.items():
                    if header in cols:
                        verified_str = verify_datatype(header, dtype)
                        if verified_str == "F":
                            re_failed.append(header)

                if len(re_failed) == 0:
                    log.write(" -Dtype fixes applied successfully for requested columns.", verbose=verbose)
                    return sumstats
                else:
                    log.warning("Some columns still have incompatible dtypes after fixes: {}".format(",".join(re_failed)), verbose=verbose)

            # Suggest manual next steps if fix was False or fixes failed
            if "CHR" in failed:
                log.warning("Consider using Sumstatsp.fix_chr() to fix CHR dtype", verbose=verbose)
            if "POS" in failed:
                log.warning("Consider using Sumstatsp.fix_pos() to fix POS dtype", verbose=verbose)
            if ("EA" in failed) or ("NEA" in failed) or ("REF" in failed) or ("ALT" in failed):
                log.warning("Consider using Sumstatsp.fix_allele() to fix allele dtypes", verbose=verbose)
            if ("SNPID" in failed) or ("rsID" in failed):
                log.warning("Consider using Sumstatsp.fix_id() to fix ID dtypes", verbose=verbose)

            stats_candidates = [
                "BETA","SE","Z","P","MLOG10P","OR","OR_95L","OR_95U",
                "HR","HR_95L","HR_95U","CHISQ","F","T","SNPR2","INFO"
            ]
            present_stats = [c for c in stats_candidates if c in sumstats.columns]
            bad_stats = [c for c in present_stats if verify_datatype(c, sumstats.schema[c]) == "F"]
            if len(bad_stats) > 0:
                for c in bad_stats:
                    log.warning("Consider using Sumstatsp.check_sanity() to check {} statistics".format(c), verbose=verbose)
        except:
            pass

        if fix is True:
            raise ValueError("Failed to fix dtypes for requested columns: {}".format(",".join(failed)))
        else:
            raise ValueError("Please fix dtypes using Sumstatsp methods (fix_chr, fix_pos, fix_allele, fix_id) for: {}. If statistics columns are present, consider Sumstatsp.check_sanity().".format(",".join(failed)))
    
    return None

