from typing import TYPE_CHECKING, Union, Optional, Dict, Any, List
import re
import pandas as pd
import numpy as np
import polars as pl
from os import path
from pathlib import Path
from functools import wraps
from gwaslab.info.g_Log import Log
from gwaslab.info.g_vchange_status import vchange_status
from gwaslab.qc.qc_fix_sumstats import _sort_coordinate
from gwaslab.qc.qc_fix_sumstats import _process_build
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_common_data import get_high_ld
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.hm.hm_harmonize_sumstats import is_palindromic
from gwaslab.info.g_version import _get_version
import gc

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

_HAPMAP_DF_CACHE = {}
_HAPMAP_FULL_CACHE = {}

def _get_hapmap_df_polars(build: str) -> pl.DataFrame:
    """Get Hapmap3 coordinates as a Polars DataFrame for fast lookup."""
    if build in _HAPMAP_DF_CACHE:
        return _HAPMAP_DF_CACHE[build]
    
    # Read directly with Polars - no pandas conversion
    base = Path(__file__).parents[1] / "data" / "hapmap3_SNPs"
    if build == "19":
        p = base / "hapmap3_db150_hg19.snplist.gz"
    else:
        p = base / "hapmap3_db151_hg38.snplist.gz"
    
    # Read with Polars - keep as Polars DataFrame
    # File is tab-separated, so we use "\t" as separator
    try:
        hapmap_df = pl.read_csv(
            p,
            separator="\t",
            columns=["#CHROM", "POS"],
            schema={"#CHROM": pl.Int64, "POS": pl.Int64},
            null_values=["", "NA", "N/A"]
        ).rename({"#CHROM": "CHR"})
    except Exception:
        # Fallback to pandas if Polars fails, then convert to Polars
        df = pd.read_csv(p, sep="\t", usecols=["#CHROM", "POS"], dtype={"#CHROM": "Int64", "POS": "Int64"})
        df = df.rename(columns={"#CHROM": "CHR"})
        hapmap_df = pl.from_pandas(df[["CHR", "POS"]])
    
    _HAPMAP_DF_CACHE[build] = hapmap_df
    return hapmap_df

def _get_hapmap_full_polars(build: str, include_alleles: bool = True) -> pl.DataFrame:
    """Get full Hapmap3 data as a Polars DataFrame with rsid, CHR, POS, and optionally A1, A2.
    
    Parameters
    ----------
    build : str
        Genome build version ("19" or "38")
    include_alleles : bool, default=True
        If True, includes A1 and A2 columns. If False, only includes rsid, CHR, POS.
    
    Returns
    -------
    pl.DataFrame
        Polars DataFrame with hapmap3 data
    """
    cache_key = f"{build}_{include_alleles}"
    if cache_key in _HAPMAP_FULL_CACHE:
        return _HAPMAP_FULL_CACHE[cache_key]
    
    # Read directly with Polars - no pandas conversion
    base = Path(__file__).parents[1] / "data" / "hapmap3_SNPs"
    if build == "19":
        p = base / "hapmap3_db150_hg19.snplist.gz"
    else:
        p = base / "hapmap3_db151_hg38.snplist.gz"
    
    # Determine columns to read
    columns = ["#CHROM", "POS", "rsid"]
    if include_alleles:
        columns.extend(["A1", "A2"])
    
    # Read with pandas first (handles whitespace-separated files with \s+)
    # Then convert to Polars for fast operations
    # This is still much faster than reading from disk every time due to caching
    df = pd.read_csv(
        p, 
        sep=r"\s+", 
        usecols=columns, 
        dtype={"#CHROM": "Int64", "POS": "Int64", "rsid": "string"}
    )
    df = df.rename(columns={"#CHROM": "CHR"})
    
    # Convert to Polars
    hapmap_df = pl.from_pandas(df)
    
    # Ensure A1 and A2 are strings if included
    if include_alleles:
        hapmap_df = hapmap_df.with_columns([
            pl.col("A1").cast(pl.Utf8),
            pl.col("A2").cast(pl.Utf8)
        ])
    
    _HAPMAP_FULL_CACHE[cache_key] = hapmap_df
    return hapmap_df

def with_logging_filter(start_to_msg: str, finished_msg: str) -> Any:
    def decorator(func):
        @wraps(func)  # This preserves the original function's metadata including __doc__
        def wrapper(*args, **kwargs):
            import inspect
            sig = inspect.signature(func)
            bound_kwargs = sig.bind(*args, **kwargs)
            bound_kwargs.apply_defaults()
            
            log = bound_kwargs.arguments.get('log', Log())
            verbose = bound_kwargs.arguments.get('verbose', True)
            
            # Log start message
            log.log_operation_start(start_to_msg, version=_get_version(), verbose=verbose)
            
            # Try to find Sumstats object instance for shape/memory tracking
            # First check if log has a reference to Sumstats object
            sumstats_obj = getattr(log, '_sumstats_obj', None)
            
            # If not found, try to find it in arguments
            if sumstats_obj is None:
                try:
                    from gwaslab.g_Sumstats import Sumstats
                    # Check if 'self' is a Sumstats instance (for methods)
                    if args and isinstance(args[0], Sumstats):
                        sumstats_obj = args[0]
                    # Check if 'sumstats_obj' is in arguments
                    elif 'sumstats_obj' in bound_kwargs.arguments:
                        potential_obj = bound_kwargs.arguments['sumstats_obj']
                        if isinstance(potential_obj, Sumstats):
                            sumstats_obj = potential_obj
                    # Check if any kwarg is a Sumstats instance
                    else:
                        for value in bound_kwargs.arguments.values():
                            if isinstance(value, Sumstats):
                                sumstats_obj = value
                                break
                except:
                    pass
            
            # Extract sumstats DataFrame
            # Try multiple possible parameter names
            sumstats = None
            for param_name in ['sumstats', 'sumstats_obj', 'sumstats_or_dataframe']:
                if param_name in bound_kwargs.arguments:
                    potential_sumstats = bound_kwargs.arguments[param_name]
                    # Check if it's a DataFrame directly
                    if isinstance(potential_sumstats, pd.DataFrame):
                        sumstats = potential_sumstats
                        break
                    # Check if it's a Sumstats object (extract .data)
                    elif sumstats_obj is None and potential_sumstats is not None:
                        try:
                            from gwaslab.g_Sumstats import Sumstats
                            if isinstance(potential_sumstats, Sumstats):
                                sumstats_obj = potential_sumstats
                                sumstats = potential_sumstats.data
                                break
                        except:
                            pass
            
            # If not found in named arguments, check first positional argument
            if sumstats is None and args:
                first_arg = args[0]
                if isinstance(first_arg, pd.DataFrame):
                    sumstats = first_arg
                elif sumstats_obj is None:
                    try:
                        from gwaslab.g_Sumstats import Sumstats
                        if isinstance(first_arg, Sumstats):
                            sumstats_obj = first_arg
                            sumstats = first_arg.data
                    except:
                        pass
            
            # If we have a sumstats_obj but no sumstats DataFrame yet, extract it
            if sumstats is None and sumstats_obj is not None:
                try:
                    sumstats = sumstats_obj.data
                except:
                    pass
            
            # If still not found, default to empty DataFrame
            if sumstats is None:
                sumstats = pd.DataFrame()
            
            check_dataframe_shape(sumstats=sumstats, log=log, verbose=verbose, sumstats_obj=sumstats_obj)
            prenum = len(sumstats)

            # Execute the original function
            result = func(*args, **kwargs)

            afternum = len(result)
            log.log_variants_filtered(prenum - afternum, verbose=verbose)
            check_dataframe_shape(sumstats=result, log=log, verbose=verbose, sumstats_obj=sumstats_obj)
            # Log finish message
            log.write(f"Finished {finished_msg}.", verbose=verbose)
            
            return result
        return wrapper
    return decorator


"""
Filtering and value manipulation utilities for GWAS summary statistics
Provides functions to filter variants based on various criteria including 
value thresholds, genomic regions, and special variant types.
"""

@with_logging_filter("filter variants by condition...","filtering variants")

def _filter_values(sumstats_obj: Union['Sumstats', pd.DataFrame], expr: str, remove: bool = False, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Filter variants based on a query expression.
    
    Parameters:
    -----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to filter.
    expr : str
        Query expression using pandas.DataFrame.query syntax
    remove : bool, default=False
        If True, removes variants meeting the condition
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.
    
    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame directly
        sumstats = sumstats_obj
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
    log.write(" -Expression:",expr, verbose=verbose)
    
    # Extract @ variables from caller's frame and prepare local_dict for query()
    import inspect
    local_dict = {}
    if '@' in expr:
        # Find variables referenced with @ syntax (e.g., @odd_chromosomes)
        pattern = r'@(\w+)'
        matches = re.findall(pattern, expr)
        if matches:
            # Search through call stack to find variables
            frame = inspect.currentframe()
            try:
                # Start from the frame that called this function (skip current frame)
                current_frame = frame.f_back
                # Search up to 5 frames deep to find the variables
                for _ in range(5):
                    if current_frame is None:
                        break
                    frame_locals = current_frame.f_locals
                    frame_globals = current_frame.f_globals
                    # Extract variables referenced with @
                    for var_name in matches:
                        if var_name not in local_dict:  # Don't override if already found
                            if var_name in frame_locals:
                                local_dict[var_name] = frame_locals[var_name]
                            elif var_name in frame_globals:
                                local_dict[var_name] = frame_globals[var_name]
                    # Move to next frame up
                    current_frame = current_frame.f_back
            finally:
                del frame
    
    # Try standard query first
    try:
        if local_dict:
            sumstats = sumstats.query(expr, engine='python', local_dict=local_dict).copy()
        else:
            sumstats = sumstats.query(expr, engine='python').copy()
    except Exception as e:
        # Check if error is related to string slicing (which query() doesn't support)
        error_msg = str(e).lower()
        if "slice" in error_msg or "not a supported function" in error_msg or ".str" in expr:
            # String slicing or .str accessor detected - use eval() with column namespace instead
            log.write(" -String slicing or .str accessor detected, using alternative evaluation method...", verbose=verbose)
            sumstats_copy = sumstats.copy()
            
            # Detect columns used with .str accessor in the expression
            str_columns = set()
            # Pattern to match column names before .str (e.g., "STATUS.str" -> "STATUS")
            pattern = r'(\w+)\.str\b'
            matches = re.findall(pattern, expr)
            str_columns.update(matches)
            
            # Convert columns that use .str accessor to string type
            for col in str_columns:
                if col in sumstats_copy.columns:
                    # Convert to string, handling NaN values
                    sumstats_copy[col] = sumstats_copy[col].astype(str).replace('nan', '')
            
            # Create namespace with all columns accessible
            env = {col: sumstats_copy[col] for col in sumstats_copy.columns}
            env.update({'pd': pd, 'np': np})
            # Add local variables from caller (for @ variable references)
            env.update(local_dict)
            
            # Replace @variable with variable in expression for eval()
            eval_expr = expr
            for var_name in local_dict.keys():
                eval_expr = re.sub(r'@' + re.escape(var_name) + r'\b', var_name, eval_expr)
            
            # Evaluate expression using eval() which supports string slicing
            try:
                mask = eval(eval_expr, {"__builtins__": {}}, env)
                sumstats = sumstats_copy.loc[mask].copy()
            except Exception as eval_error:
                error_str = str(eval_error)
                # Check if error is about .str accessor on non-string values
                if "can only use .str accessor" in error_str.lower():
                    raise ValueError(f"Column used with .str accessor must be convertible to string. "
                                   f"Expression: {expr}. "
                                   f"Error: {error_str}. "
                                   f"Note: Ensure the column can be converted to string type.")
                else:
                    raise ValueError(f"Unable to evaluate expression with string slicing: {expr}. "
                                   f"Error: {str(eval_error)}. "
                                   f"Note: String slicing operations (e.g., .str[:2]) are not supported "
                                   f"in pandas query(). Consider using .str.slice() or creating a temporary column.")
        else:
            # Re-raise original error if it's not a slicing issue
            raise

    gc.collect()
    # Update sumstats_obj.data with filtered result only if called with Sumstats object
    if not isinstance(sumstats_obj, pd.DataFrame):
        sumstats_obj.data = sumstats
    return sumstats

@with_logging_filter("filter out variants based on threshold values...", "filtering variants")
def _filter_out(
    sumstats_obj: Union['Sumstats', pd.DataFrame], 
    interval: Dict[str, Any] = {}, 
    lt: Dict[str, Any] = {}, 
    gt: Dict[str, Any] = {}, 
    eq: Dict[str, Any] = {}, 
    remove: bool = False, 
    verbose: bool = True, 
    log: Log = Log()
) -> pd.DataFrame:
    """
    Filter out variants based on threshold values.
    
    Parameters:
    -----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to filter.
    lt : dict, default={}
        Dictionary of {column: threshold} for lower bounds (variant values < threshold will be removed)
    gt : dict, default={}
        Dictionary of {column: threshold} for upper bounds (variant values > threshold will be removed)
    eq : dict, default={}
        Dictionary of {column: value} for equality checks (variant values == value will be removed)
    remove : bool, default=False
        If True, removes variants meeting the condition
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.
    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame directly
        sumstats = sumstats_obj.copy()
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data.copy()
    for key,threshold in gt.items():
        num = len(sumstats.loc[sumstats[key]>threshold,:])
        log.write(" -Removing "+ str(num) +" variants with "+key+" > "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]<threshold,:]
    for key,threshold in lt.items():
        num = len(sumstats.loc[sumstats[key]<threshold,:])
        log.write(" -Removing "+ str(num) +" variants with "+key+" < "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]>threshold,:]
    for key,threshold in eq.items():
        num = len(sumstats.loc[sumstats[key]==threshold,:])
        log.write(" -Removing "+ str(num) +" variants with "+key+" = "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]!=threshold,:]
    gc.collect()
    # Update sumstats_obj.data with filtered result only if called with Sumstats object
    if not isinstance(sumstats_obj, pd.DataFrame):
        sumstats_obj.data = sumstats
    return sumstats

@with_logging_filter("filter in variants based on threshold values", "filtering variants")
def _filter_in(
    sumstats_obj: Union['Sumstats', pd.DataFrame], 
    lt: Dict[str, Any] = {}, 
    gt: Dict[str, Any] = {}, 
    eq: Dict[str, Any] = {}, 
    remove: bool = False, 
    verbose: bool = True, 
    log: Log = Log()
) -> pd.DataFrame:
    """
    Filter in variants based on threshold values.
    
    Parameters:
    -----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to filter.
    lt : dict, default={}
        Dictionary of {column: threshold} for lower bounds (variant values < threshold will be kept)
    gt : dict, default={}
        Dictionary of {column: threshold} for upper bounds (variant values > threshold will be kept)
    eq : dict, default={}
        Dictionary of {column: value} for equality checks (variant values == value will be kept)
    remove : bool, default=False
        If True, removes variants meeting the condition
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.
    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame directly
        sumstats = sumstats_obj.copy()
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data.copy()
    for key,threshold in gt.items():
        num = len(sumstats.loc[sumstats[key]>threshold,:])
        log.write(" -Keeping "+ str(num) +" variants with "+key+" > "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]>threshold,:]
    for key,threshold in lt.items():
        num = len(sumstats.loc[sumstats[key]<threshold,:])
        log.write(" -Keeping "+ str(num) +" variants with "+key+" < "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]<threshold,:]
    for key,threshold in eq.items():
        num = len(sumstats.loc[sumstats[key]==threshold,:])
        log.write(" -Keeping "+ str(num) +" variants with "+key+" = "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]==threshold,:]
    gc.collect()
    # Update sumstats_obj.data with filtered result only if called with Sumstats object
    if not isinstance(sumstats_obj, pd.DataFrame):
        sumstats_obj.data = sumstats
    return sumstats

@with_logging_filter("filter in variants if in intervals defined in bed files", "filtering variants")
def _filter_region_in(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], path: Optional[str] = None, chrom: str = "CHR", pos: str = "POS", high_ld: bool = False, build: str = "19", verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Keep variants located within specified genomic regions from a BED file.
    
    Parameters:
    -----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    path : str or None, default=None
        Path to BED file containing regions of interest
        Column name for position information
    high_ld : bool, default=False
        If True, uses high LD regions from the specified build
    build : str, default="19"
        Genome build version (19 or 38) for high LD regions
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table containing only variants in the specified regions
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    sumstats = _sort_coordinate(sumstats,verbose=verbose)
    if high_ld is True:
        path = get_high_ld(build=build)
        log.write(" -Loading bed format file for hg"+build, verbose=verbose)

    else:
        log.write(" -Loading bed format file: " , path, verbose=verbose)
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    dic=get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr")
    bed[0]=bed[0].map(dic).astype("string")
    bed[0]=bed[0].astype("Int64")
    sumstats = sumstats.sort_values(["CHR","POS"])
    
    if len(bed)<100:
        log.write(" -Bed file < 100 lines: using pd IntervalIndex... ", verbose=verbose)
        for i in sumstats[chrom].unique():
            if sum(bed[0]==i)>0:
                interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==i,"tuple"])
                sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
            else:
                continue
    else:
        log.write(" -Bed file > 100 lines: using two pointers, please make files are all sorted... ", verbose=verbose)
        bed_num  =0
        bed_chr   =bed.iloc[bed_num,0]
        bed_left  =bed.iloc[bed_num,1]
        bed_right =bed.iloc[bed_num,2]
        
        sum_num=0
        sum_chr_col = sumstats.columns.get_loc(chrom)
        sum_pos_col = sumstats.columns.get_loc(pos)
        sum_ind_col = sumstats.columns.get_loc("bed_indicator")
        while bed_num<len(bed) and sum_num<len(sumstats):
            #sumstats variant chr < bed chr
            if sumstats.iat[sum_num,sum_chr_col]<bed_chr:
                # next variant
                sum_num+=1
                continue
            #sumstats variant chr > bed chr
            elif sumstats.iat[sum_num,sum_chr_col]>bed_chr:
                # next bed record
                bed_num+=1
                bed_chr  =bed.iat[bed_num,0]
                bed_left  = bed.iat[bed_num,1]
                bed_right  = bed.iat[bed_num,2]
                continue
            #sumstats variant chr == bed chr
            else:
                #sumstats variant pos < bed left
                if sumstats.iat[sum_num,sum_pos_col]<bed_left:
                    # next variant
                    sum_num+=1
                    continue
                #sumstats variant pos > bed right
                elif sumstats.iat[sum_num,sum_pos_col]>bed_right:
                    # next bed record
                    bed_num+=1
                    bed_chr  =bed.iat[bed_num,0]
                    bed_left  = bed.iat[bed_num,1]
                    bed_right  = bed.iat[bed_num,2]
                # bed left  < sumstats variant  pos < bed right
                else:
                    # set to true
                    sumstats.iat[sum_num,sum_ind_col]=True
                    # next variant
                    sum_num+=1
                           
    ## in
    
    sumstats = sumstats.loc[sumstats["bed_indicator"],:]
    log.write(" -Number of variants in the specified regions to keep:",sum(sumstats["bed_indicator"]), verbose=verbose)
    sumstats = sumstats.drop(columns="bed_indicator")
    gc.collect()
    return sumstats

@with_logging_filter("filter out variants if in intervals defined in bed files", "filtering variants")
def _filter_region_out(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], path: Optional[str] = None, chrom: str = "CHR", pos: str = "POS", high_ld: bool = False, build: str = "19", verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Remove variants located within specified genomic regions from a BED file.
    
    Parameters:
    -----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    path : str or None, default=None
        Path to BED file containing regions to exclude
    high_ld : bool, default=False
        If True, excludes variants in high LD regions from the specified build
    build : str, default="19"
        Genome build version (19 or 38) for high LD regions
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.
        
    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table with variants in specified regions removed
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    sumstats = _sort_coordinate(sumstats,verbose=verbose)
    if high_ld is True:
        path = get_high_ld(build=build)
        log.write(" -Loading bed format file for hg"+build, verbose=verbose)

    else:
        log.write(" -Loading bed format file: " , path, verbose=verbose)
            
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    
    dic=get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr")
    bed[0]=bed[0].map(dic).astype("string")
    bed[0]=bed[0].astype("Int64")
    
    if len(bed)<100:
        log.write(" -Bed file < 100 lines: using pd IntervalIndex... ", verbose=verbose)
        for i in sumstats[chrom].unique():
            if sum(bed[0]==i)>0:
                interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==i,"tuple"])
                sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
            else:
                continue
    else:
        log.write(" -Bed file > 100 lines: using two pointers, please make files are all sorted... ", verbose=verbose)
        bed_num  =0
        bed_chr  =bed.iloc[bed_num,0]
        bed_left  =bed.iloc[bed_num,1]
        bed_right =bed.iloc[bed_num,2]
        
        sum_num=0
        sum_chr_col = sumstats.columns.get_loc(chrom)
        sum_pos_col = sumstats.columns.get_loc(pos)
        sum_ind_col = sumstats.columns.get_loc("bed_indicator")
        while bed_num<len(bed) and sum_num<len(sumstats):
            if sumstats.iat[sum_num,sum_chr_col]<bed_chr:
                sum_num+=1
                continue
            elif sumstats.iat[sum_num,sum_chr_col]>bed_chr:
                bed_num+=1
                bed_chr  =bed.iat[bed_num,0]
                bed_left  = bed.iat[bed_num,1]
                bed_right  = bed.iat[bed_num,2]
                continue
            else:
                if sumstats.iat[sum_num,sum_pos_col]<bed_left:
                    sum_num+=1
                    continue
                elif sumstats.iat[sum_num,sum_pos_col]>bed_right:
                    bed_num+=1
                    bed_chr  =bed.iat[bed_num,0]
                    bed_left  = bed.iat[bed_num,1]
                    bed_right  = bed.iat[bed_num,2]
                else:
                    sumstats.iat[sum_num,sum_ind_col]=True
                    sum_num+=1
                           
    ## out
    sumstats = sumstats.loc[~sumstats["bed_indicator"],:]
    sumstats = sumstats.drop(columns="bed_indicator")
    gc.collect()
    return sumstats

@with_logging("infer genome build version using hapmap3 SNPs", 
              "inferring genome build version using hapmap3 SNPs",
              start_function=".infer_build()",
              start_cols=["CHR","POS"],
              check_dtype=True)
def _infer_build(sumstats: Union['Sumstats', pd.DataFrame], status: str = "STATUS", chrom: str = "CHR", pos: str = "POS", 
               ea: str = "EA", nea: str = "NEA", build: str = "19",
               change_status: bool = True, 
               verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Infer genome build version using Hapmap3 SNPs.
    


    
    Returns
    -------
    pandas.DataFrame
        Updated summary statistics DataFrame with inferred build version.
        If called with Sumstats object, also updates self.build and self.meta["gwaslab"]["genome_build"].
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats, pd.DataFrame):
        sumstats_data = sumstats
        is_dataframe = True
    else:
        sumstats_obj = sumstats
        sumstats_data = sumstats_obj.data
        is_dataframe = False

    inferred_build = "Unknown"
    log.write(" -Loading Hapmap3 variants data...", verbose=verbose)
    hapmap3_df_19 = _get_hapmap_df_polars("19")
    hapmap3_df_38 = _get_hapmap_df_polars("38")
    log.write(" -CHR and POS will be used for matching...", verbose=verbose)
    
    # Convert sumstats to Polars for fast matching (no pandas operations)
    # with_logging decorator with check_dtype=True ensures CHR and POS are Int64
    sumstats_pl = pl.from_pandas(sumstats_data[[chrom, pos]])
    
    # Filter out null values and ensure correct types
    sumstats_pl = sumstats_pl.filter(
        pl.col(chrom).is_not_null() & pl.col(pos).is_not_null()
    ).with_columns([
        pl.col(chrom).cast(pl.Int64).alias("CHR"),
        pl.col(pos).cast(pl.Int64).alias("POS")
    ])
    
    # Use Polars semi_join for fast counting - more efficient than full join
    # semi_join only checks existence without creating full join result
    # This is highly optimized in Polars and avoids pandas conversions
    match_count_for_19 = sumstats_pl.join(
        hapmap3_df_19, on=["CHR", "POS"], how="semi"
    ).height
    match_count_for_38 = sumstats_pl.join(
        hapmap3_df_38, on=["CHR", "POS"], how="semi"
    ).height
    log.write(" -Matching variants for hg19: num_hg19 = ", match_count_for_19, verbose=verbose)
    log.write(" -Matching variants for hg38: num_hg38 = ", match_count_for_38, verbose=verbose)
    
    if max(match_count_for_19, match_count_for_38)<10000:
        log.warning("Please be cautious due to the limited number of variants.", verbose=verbose) 
    
    if match_count_for_19 > match_count_for_38:
        log.write(" -Since num_hg19 >> num_hg38, set the genome build to hg19 for the STATUS code....", verbose=verbose) 
        if change_status==True:
            sumstats_data[status] = vchange_status(sumstats_data[status],1,"9","1")
        inferred_build="19"
    elif match_count_for_19 < match_count_for_38:
        log.write(" -Since num_hg19 << num_hg38, set the genome build to hg38 for the STATUS code....", verbose=verbose) 
        if change_status==True:
            sumstats_data[status] = vchange_status(sumstats_data[status],1,"9","3")
            sumstats_data[status] = vchange_status(sumstats_data[status],2,"9","8")
        inferred_build="38"
    else:
        log.write(" -Since num_hg19 = num_hg38, unable to infer...", verbose=verbose) 
    
    # Update Sumstats object if called with one
    if not is_dataframe:
        sumstats_obj.data = sumstats_data
        # Use setter to ensure consistency between self.build and meta
        sumstats_obj.build = inferred_build
    
    return sumstats_data

@with_logging_filter("randomly select variants from the sumstats", "sampling")
def _sampling(sumstats: pd.DataFrame, n: int = 1, p: Optional[float] = None, verbose: bool = True, log: Log = Log(), **kwargs: Any) -> pd.DataFrame:
    """
    Randomly sample variants from summary statistics.
    
    Parameters:
    -----------
    n : int, default=1
        Number of variants to sample
    p : float, default=None
        Fraction of variants to sample (alternative to n)
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.
    **kwargs : dict
        Additional arguments for pandas.DataFrame.sample
    
    Returns:
    --------
    pandas.DataFrame
        Subsampled summary statistics table
    
    """
    if p is None:
        log.write(" -Number of variants selected from the sumstats:",n, verbose=verbose)
        if n > len(sumstats):
            raise ValueError("Please input a number < {}".format(len(sumstats)))
    else:
        if p>-0.00000001 and p<1.00000001:
            log.write(" -Percentage of variants selected from the sumstats: ",p, verbose=verbose)
            n = int(len(sumstats)*p)
            log.write(" -Number of variants selected from the sumstats:",n, verbose=verbose)
        else:
            raise ValueError("Please input a number in (0,1)")
    
    if "random_state" in kwargs.keys():
        log.write(" -Random state (seed): {}".format(kwargs["random_state"]), verbose=verbose)
    else:
        kwargs["random_state"] = np.random.randint(0,4294967295)
        log.write(" -Random state (seed): {}".format(kwargs["random_state"]), verbose=verbose)
    sampled = sumstats.sample(n=n,**kwargs)
    gc.collect()
    return sampled

@with_logging_filter("extract variants in the flanking regions", "extracting variants in the flanking regions")
def _get_flanking(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], snpid: str, windowsizekb: int = 500, verbose: bool = True, log: Log = Log(), **kwargs: Any) -> pd.DataFrame:
    """
    Extract variants in flanking regions around a specified variant.
    
    Parameters:
    -----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    snpid : str
        ID of the central variant (must exist in SNPID column)
    windowsizekb : int, default=500
        Size of flanking region in kilobases
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Variants in flanking regions
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    log.write(" - Central variant: {}".format(snpid))
    
    row = sumstats.loc[sumstats["SNPID"]==snpid,:]
    
    for index, row in row.iterrows():
        chrom = row["CHR"]
        left = row["POS"] - 1000 * windowsizekb
        right = row["POS"] + 1000 * windowsizekb
    
    log.write(" - Flanking regions: {}:{}-{}".format(chrom, left, right ))

    flanking = sumstats.query("CHR==@chrom & POS > @left & POS < @right ", engine='python')
    
    log.write(" - Extracted {} variants in the regions.".format(len(flanking)),verbose=verbose)
    log.write("Finished extracting variants in the flanking regions.",verbose=verbose)

    return flanking

@with_logging_filter("extract variants in the flanking regions using rsID or SNPID", "extracting variants in the flanking regions")
def _get_flanking_by_id(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], snpid: Union[str, List[str]], windowsizekb: int = 500, verbose: bool = True, log: Log = Log(), **kwargs: Any) -> pd.DataFrame:
    """
    Extract variants in flanking regions using rsID or SNPID.
    
    Parameters:
    -----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    snpid : str or list
        Variant ID(s) to use as center (searches in rsID/SNPID columns)
    windowsizekb : int, default=500
        Size of flanking region in kilobases
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Variants in flanking regions
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    log.write(" - Central variants: {}".format(snpid), verbose=verbose)
    log.write(" - Flanking windowsize in kb: {}".format(windowsizekb), verbose=verbose)

    if type(snpid) == str:
        snpid = [snpid]
    
    if "rsID" in sumstats.columns and "SNPID" not in sumstats.columns:
        is_specified = sumstats["rsID"].isin(snpid)
    elif "rsID" not in sumstats.columns and "SNPID" in sumstats.columns:
        is_specified = sumstats["SNPID"].isin(snpid)
    else:
        is_specified = sumstats["rsID"].isin(snpid) | sumstats["SNPID"].isin(snpid)

    row = sumstats.loc[is_specified,:]
    
    is_flanking = None
    for index, row in row.iterrows():
        chrom = row["CHR"]
        left =  row["POS"] - 1000 * windowsizekb
        right = row["POS"] + 1000 * windowsizekb
        is_flancking_in_this_region = (sumstats["CHR"] == chrom) & (sumstats["POS"] >= left) & (sumstats["POS"] <= right)
        
        log.write(" - Variants in flanking region {}:{}-{} : {}".format(chrom, left, right, sum(is_flancking_in_this_region) ))
        
        if is_flanking is None:
            is_flanking = is_flancking_in_this_region
        else:
            is_flanking = is_flanking | is_flancking_in_this_region
    
    flanking = sumstats.loc[is_flanking,:]
    
    log.write(" - Extracted {} variants in the regions.".format(len(flanking)),verbose=verbose)
    log.write("Finished extracting variants in the flanking regions.",verbose=verbose)

    return flanking

@with_logging_filter("extract variants in the flanking regions using CHR and POS", "extracting variants in the flanking regions")
def _get_flanking_by_chrpos(sumstats: pd.DataFrame, chrpos: Union[List[Union[int, float]], List[List[Union[int, float]]]], windowsizekb: int = 500, verbose: bool = True, log: Log = Log(), **kwargs: Any) -> pd.DataFrame:
    """
    Extract variants in flanking regions using chromosome and position.
    
    Parameters:
    -----------
    chrpos : list or list of lists
        Chromosome and position coordinate(s) to use as center
        Format: [chromosome, position]
    windowsizekb : int, default=500
        Size of flanking region in kilobases
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Variants in flanking regions
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats, pd.DataFrame):
        # Already a DataFrame, use as is
        pass
    else:
        # Assume it's a Sumstats object, extract .data
        sumstats = sumstats.data
    
    log.write(" - Central positions: {}".format(chrpos), verbose=verbose)
    log.write(" - Flanking windowsize in kb: {}".format(windowsizekb), verbose=verbose)

    if type(chrpos) == tuple or type(chrpos) == list:
        chrpos_to_check = [chrpos]
    else:
        chrpos_to_check = chrpos

    is_flanking = None
    
    for index, row in enumerate(chrpos_to_check):
        chrom = row[0]
        left =  row[1] - 1000 * windowsizekb
        right = row[1] + 1000 * windowsizekb
        is_flancking_in_this_region = (sumstats["CHR"] == chrom) & (sumstats["POS"] >= left) & (sumstats["POS"] <= right)
        
        log.write(" - Variants in flanking region {}:{}-{} : {}".format(chrom, left, right, sum(is_flancking_in_this_region) ))
        
        if is_flanking is None:
            is_flanking = is_flancking_in_this_region
        else:
            is_flanking = is_flanking | is_flancking_in_this_region
    
    flanking = sumstats.loc[is_flanking,:]
    
    log.write(" - Extracted {} variants in the regions.".format(len(flanking)),verbose=verbose)
    return flanking

@with_logging_filter("filter palindromic variants", "filtering variants")
def _filter_palindromic(sumstats: pd.DataFrame, mode: str = "in", ea: str = "EA", nea: str = "NEA", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Filter palindromic variants based on allele symmetry.
    
    Parameters:
    -----------
    mode : str, default="in"
        "in" to keep palindromic variants, "out" to remove them
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table
    """
    is_palindromic_snp = is_palindromic(sumstats[[nea,ea]],a1=nea,a2=ea)   
    
    log.write(" -Identified palindromic variants: {}".format(sum(is_palindromic_snp)),verbose=verbose)
    
    if mode=="in":
        palindromic = sumstats.loc[is_palindromic_snp,:]
    else:
        palindromic = sumstats.loc[~is_palindromic_snp,:]

    return palindromic

@with_logging_filter("filter indels", "filtering variants")
def _filter_indel(sumstats: pd.DataFrame, mode: str = "in", ea: str = "EA", nea: str = "NEA", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Filter indels based on allele length differences.
    
    Parameters:
    -----------
    mode : str, default="in"
        "in" to keep indels, "out" to remove indels
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table
    """
    is_indel = (sumstats[ea].str.len()!=sumstats[nea].str.len()) 
    
    log.write(" -Identified indels: {}".format(sum(is_indel)),verbose=verbose)
    if mode=="in":
        indel = sumstats.loc[is_indel,:]
    else:
        indel = sumstats.loc[~is_indel,:]
    return indel

@with_logging_filter("filter SNPs", "filtering variants")
def _filter_snp(sumstats: pd.DataFrame, mode: str = "in", ea: str = "EA", nea: str = "NEA", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Filter SNPs based on allele length.
    
    Parameters:
    -----------
    mode : str, default="in"
        "in" to keep SNPs, "out" to remove SNPs
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table
    """
    is_snp = (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()==1)
    
    log.write(" -Identified SNPs: {}".format(sum(is_snp)),verbose=verbose)
    if mode=="in":
        snp = sumstats.loc[is_snp,:]
    else:
        snp = sumstats.loc[~is_snp,:]
    return snp

@with_logging_filter("exclude variants in HLA regions", "filtering variants")
def _exclude_hla(sumstats: pd.DataFrame, chrom: str = "CHR", pos: str = "POS", lower: Optional[int] = None, upper: Optional[int] = None, build: Optional[str] = None, mode: str = "xmhc", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Exclude variants in HLA regions based on genomic coordinates.
    
    Parameters:
    -----------
    chrom : str, default="CHR"
        Column name for chromosome information
    pos : str, default="POS"
        Column name for position information
    lower : int or None, default=None
        Lower bound of genomic region
    upper : int or None, default=None
        Upper bound of genomic region
    build : str or None, default=None
        Genome build version (19 or 38)
    mode : str, default="xmhc"
        "xmhc" for extended MHC region, "hla" or "mhc" for classical HLA region
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table with HLA variants removed
    """
    if build is not None:
        build = _process_build(build = build,log = log,verbose = verbose)
        # xMHC : HIST1H2AA ~ 7.6mb ~ RPL12P1
        # reference: Horton, R., Wilming, L., Rand, V., Lovering, R. C., Bruford, E. A., Khodiyar, V. K., ... & Beck, S. (2004). Gene map of the extended human MHC. Nature Reviews Genetics, 5(12), 889-899.
        # hg38:  25,726,063 ~ 33,400,644
        # hg19 : 25,726,291 ~ 33,368,421

        # HLA : GABBR1 ~ 3.78mb ~ KIFC1
        # reference: Shiina, T., Hosomichi, K., Inoko, H., & Kulski, J. K. (2009). The HLA genomic loci map: expression, interaction, diversity and disease. Journal of human genetics, 54(1), 15-39. 
        # hg38:  29,602,238 ~ 33,409,896
        # hg19:  29,570,015 ~ 33,377,673

        if build == "19":
            if mode =="xmhc":
                lower=25000000
                upper=34000000
            if mode =="hla" or mode =="mhc":
                lower=29500000
                upper=33500000
        if build == "38":
            if mode =="xmhc":
                lower=25000000
                upper=34000000
            if mode =="hla" or mode =="mhc":
                lower=29500000
                upper=33500000
    else:
        # -> 25,000,000 ~ 34,000,000
        if mode =="xmhc":
            lower=25000000
            upper=34000000
        if mode =="hla" or mode =="mhc":
            lower=29500000
            upper=33500000
        
    raw_len = len(sumstats)
    
    if str(sumstats[chrom].dtype) == "string":
        is_in_hla = (sumstats[chrom].astype("string")=="6")&(sumstats[pos]>lower)&(sumstats[pos]<upper)
    else:
        is_in_hla = (sumstats[chrom]==6)&(sumstats[pos]>lower)&(sumstats[pos]<upper)
    
    sumstats = sumstats.loc[~is_in_hla, : ]
    
    after_len = len(sumstats)
    
    log.write(" -Excluded {} variants in HLA region (chr6: {}-{} )...".format(raw_len - after_len,lower,upper),verbose=verbose)
    
    return sumstats

@with_logging_filter("exclude variants on sex chromosomes", "filtering variants")
def _exclude_sexchr(sumstats: pd.DataFrame, chrom: str = "CHR", pos: str = "POS", sexchrs: List[int] = [23,24,25], log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Exclude variants on sex chromosomes.
    
    Parameters:
    -----------
    chrom : str, default="CHR"
        Column name for chromosome information
    sexchrs : list, default=[23,24,25]
        List of sex chromosome numbers to exclude
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Filtered summary statistics table with sex chromosome variants removed
    """
    raw_len = len(sumstats)
    
    if str(sumstats[chrom].dtype) == "string":
        sexchrs_string = [str(i) for i in sexchrs]
        is_in_sexchr = sumstats[chrom].astype("string").isin(sexchrs_string)
    else:
        is_in_sexchr = sumstats[chrom].isin(sexchrs)
    
    sumstats = sumstats.loc[~is_in_sexchr, : ]
    
    after_len = len(sumstats)
    
    log.write(" -Excluded {} variants on sex chromosomes ({})...".format(raw_len - after_len,sexchrs),verbose=verbose)
    
    return sumstats

@with_logging_filter("extract specific variants by ID", "extracting variants")
def _extract(sumstats: pd.DataFrame, extract: Optional[List[str]] = None, id_use: str = "SNPID", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Extract specific variants by ID.
    
    Parameters:
    -----------
    extract : list or None, default=None
        List of variant IDs to extract
    id_use : str, default="SNPID"
        Column name containing variant IDs
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Subset of summary statistics with specified variants
    """
    if extract is not None:
        log.write(" -Extracting {} variants from sumstats...".format(len(extract)),verbose=verbose)
        sumstats = sumstats.loc[sumstats[id_use].isin(extract),:]
        log.write(" -Extracted {} variants from sumstats...".format(len(sumstats)),verbose=verbose)
    return sumstats

@with_logging_filter("exclude specific variants by ID", "filtering variants")
def _exclude(sumstats: pd.DataFrame, exclude: Optional[List[str]] = None, id_use: str = "SNPID", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Exclude specific variants by ID.
    
    Parameters:
    -----------
    exclude : list or None, default=None
        List of variant IDs to exclude
    id_use : str, default="SNPID"
        Column name containing variant IDs
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Subset of summary statistics without excluded variants
    """
    if exclude is not None:
        log.write(" -Excluding {} variants from sumstats...".format(len(exclude)),verbose=verbose)
        sumstats = sumstats.loc[~sumstats[id_use].isin(exclude),:]
        log.write(" -Excluded {} variants from sumstats...".format(len(sumstats)),verbose=verbose)
    return sumstats

@with_logging_filter("filter variants within a specific genomic region", "filtering variants")
def _filter_region(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], region: Optional[List[Union[int, str]]] = None, chrom: str = "CHR", pos: str = "POS", log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Filter variants within a specific genomic region.
    
    Parameters:
    -----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    region : List
        Genomic region to filter [chr, start, end]
    verbose : bool, default=True
        If True, writes progress to log
    inplace : bool, default=False  
        If False, return a new `Sumstats` object containing the filtered results.  
        If True, apply the filter to the current object in place and return None.

    Returns:
    --------
    pandas.DataFrame
        Subset of summary statistics in the specified region
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    if region is not None:
        
        if type(region[0]) is str:
            region[0] = int(region[0])
            
        region_chr = region[0]
        region_start = region[1]
        region_end = region[2]
        
        log.write(" -Extract SNPs in region : chr{}:{}-{}...".format(region_chr, region[1], region[2]),verbose=verbose)
        
        in_region_snp = (sumstats[chrom]==region_chr) & (sumstats[pos]<region_end) & (sumstats[pos]>region_start)
        
        log.write(" -Extract SNPs in specified regions: "+str(sum(in_region_snp)),verbose=verbose)
        sumstats = sumstats.loc[in_region_snp,:]
        return sumstats.copy()    
    
def _search_variants(sumstats: pd.DataFrame, snplist: Optional[List[Union[str, List[Union[int, float]]]]] = None, 
                     snpid: str = "SNPID", rsid: str = "rsID",
                     chrom: str = "CHR", pos: str = "POS", ea: str = "EA", nea: str = "NEA",
                     log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Search for variants in summary statistics using multiple identifier formats.
    
    Parameters:
    -----------
    sumstats : pandas.DataFrame
        GWAS summary statistics table
    snplist : list or None, default=None
        List of variant identifiers to search for. Accepts multiple formats:
        - CHR:POS (e.g., '1:123456')
        - CHR-POS (e.g., '1_123456')
        - rsID (e.g., 'rs12345')
        - SNPID (e.g., '1:123456:A:G')
        - List of [CHR, POS]
        - Full variant IDs with alleles (e.g., 'chr1:123456:A:G')
    verbose : bool, default=True
        If True, writes progress to log
    
    Returns:
    --------
    pandas.DataFrame
        Subset of summary statistics containing matching variants
    
    Examples:
    ---------
    >>> variants = _search_variants(sumstats, snplist=["rs1234", "1:100500", [2, 202000], "19:45100000:C:T"])
    >>> print(variants.shape)
    """
    # create a boolean col with FALSE 
    if snplist is None:
        return sumstats.iloc[0:0,:].copy()
    
    if snpid in sumstats.columns:
        is_extract = sumstats[snpid]!=sumstats[snpid]
    else:
        is_extract = sumstats[rsid]!=sumstats[rsid]
    
    # search each variant
    for variant in snplist:        
        
        if pd.api.types.is_list_like(variant):
            # (1:1234)
            single_chrom=variant[0]
            single_pos=variant[1]
            is_extract = is_extract | ((sumstats[pos] == single_pos ) &(sumstats[chrom] == single_chrom))
        
        elif pd.api.types.is_string_dtype(type(variant)):
            # rs123
            if "rsID" in sumstats.columns:
                is_extract = is_extract | (sumstats["rsID"] == variant)

            # 1:123:A:D
            if "SNPID" in sumstats.columns:
                is_extract = is_extract | (sumstats["SNPID"] == variant)

            # 1:123:A:D -> (1:1234)
            a= re.match(r'^(chr|Chr|CHR)?(\d+)[:_-](\d+)([:_-]([ATCG]+)[:_-]([ATCG]+))?$', variant, flags=0)
            
            if a is not None:
                if a[4] is None:
                    single_chrom=int(a[2])
                    single_pos=int(a[3])
                    is_extract = is_extract | ((sumstats[pos] == single_pos ) &(sumstats[chrom] == single_chrom))
                else:
                    single_chrom = int(a[2])
                    single_pos = int(a[3])
                    single_ea = a[5]
                    single_nea = a[6]
                    a_match = ((sumstats[nea] == single_nea) & (sumstats[ea] == single_ea)) | ((sumstats[nea] == single_ea) & (sumstats[ea] == single_nea))
                    is_extract = is_extract | ((sumstats[pos] == single_pos ) &(sumstats[chrom] == single_chrom)  & a_match)
                        
    to_search =  sumstats.loc[is_extract,:].copy()
    log.write(" -Found {} variants...".format(len(to_search)),verbose=verbose)

    log.write("Finished searching variants.", verbose=verbose)
    return to_search

def _get_region_start_and_end(
    chrom: Union[str, int],
    pos: Union[int, float, str],
    windowsizekb: int = 500,
    verbose: bool = True,
    log: Optional[Log] = None,
    **kwargs: Any) -> List[Union[int, str]]:
    """
    Determine the [chr, start, end] for a region. 

    Parameters
    ----------
    chrom : str or int
        Chromosome identifier.
    pos : int or float or str
        Base-pair position.
    windowsizekb : int, default=500
        Window size in kilobases.
    verbose : bool, default=True
        Print log message.

    Returns
    -------
    list
        [chrom, start, end], where start and end are base-pair coordinates (int).
    """
    # Ensure chrom is str
    chrom = str(chrom)

    # Convert pos to int safely
    try:
        pos = int(float(pos))
    except Exception:
        raise ValueError(f"Invalid position value: {pos}")

    # Convert window to bp
    try:
        window = int(float(windowsizekb) * 1000)
    except Exception:
        raise ValueError(f"Invalid windowsizekb value: {windowsizekb}")

    start = max(1, pos - window)
    end = pos + window

    if log is not None:
        msg = f"[REGION] CHR {chrom}: {start:,} - {end:,} (±{windowsizekb}kb around {pos:,})"
        log.write(msg)
        if verbose:
            print(msg)
    try:
        chrom = int(chrom)
    except:
        pass

    return [chrom, start, end]

@with_logging_filter("filter duplicate/multiallelic variants", "filtering variants")
def _filter_dup(sumstats_obj: Union['Sumstats', pd.DataFrame], mode: str = "dm", chrom: str = "CHR", pos: str = "POS", snpid: str = "SNPID", ea: str = "EA", nea: str = "NEA", rsid: str = "rsID", 
keep: Union[str, bool] = 'first', keep_col: str = "P", remove_na: bool = False, keep_ascend: bool = True, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Filter to keep only duplicate or multiallelic variants based on user-selected criteria.
    
    This function works like remove_dup but keeps only the duplicate rows instead of removing them.
    Useful for inspecting duplicate variants before removal.

    Parameters
    ----------
    sumstats_obj : Sumstats or pandas.DataFrame
        Sumstats object or DataFrame containing the data to filter.
    mode : str, default="dm"
        String encoding the deduplication rules; may include one or more of:
        - 'ds' : Identify duplicates using SNPID.
        - 'dr' : Identify duplicates using rsID.
        - 'dc' : Identify duplicates using chromosome, position, effect allele, and non-effect allele.
        - 'm' : Identify multi-allelic variants (same chromosome + position).
    keep : {'first', 'last', False}, default 'first'
        Which duplicates to show:
        - 'first': Show all duplicates except the first occurrence (the ones that would be removed)
        - 'last': Show all duplicates except the last occurrence (the ones that would be removed)
        - False: Show all duplicate rows (all occurrences)
    keep_col : str, default="P"
        Column to sort by before filtering (for consistency with remove_dup).
    remove_na : bool, default=False
        If True, remove rows containing missing values in deduplication-relevant columns.
    keep_ascend : bool, default=True
        If True, sort in ascending order.
    verbose : bool, default=True
        If True, writes progress to log

    Returns
    -------
    pandas.DataFrame
        Summary statistics containing only duplicate and multi-allelic variants
        according to the specified mode.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame directly
        sumstats = sumstats_obj.copy()
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data.copy()
    
    log.write(" -Filtering mode: {}".format(mode), verbose=verbose)
    log.write(" -Which duplicates to show: {}".format(keep), verbose=verbose)
    
    # Sort the variants using the specified column before filtering (for consistency)
    if keep_col is not None:
        if keep_col in sumstats.columns:
            log.write(" -Sorting the sumstats using {}...".format(keep_col), verbose=verbose)
            sumstats = sumstats.sort_values(by=keep_col, ascending=keep_ascend)
        else:
            log.write(" -Column {} was not detected... skipping...".format(keep_col), verbose=verbose)
    
    total_number = len(sumstats)
    is_duplicate = pd.Series([False] * len(sumstats), index=sumstats.index)
    
    # Filter by duplicated SNPID
    if (snpid in sumstats.columns) and ("d" in mode or "s" in mode):
        log.write(" -Filtering duplicated variants based on SNPID...", verbose=verbose)
        if snpid in sumstats.columns:
            # Use keep parameter: False shows all duplicates, 'first'/'last' exclude first/last
            is_dup_snpid = sumstats.duplicated(subset=[snpid], keep=keep)
            # Exclude NA values
            is_dup_snpid = is_dup_snpid & sumstats[snpid].notna()
            is_duplicate = is_duplicate | is_dup_snpid
            log.write(" -Found {} duplicate variants based on SNPID...".format(is_dup_snpid.sum()), verbose=verbose)
    
    # Filter by duplicated rsID
    if (rsid in sumstats.columns) and ("d" in mode or "r" in mode):
        log.write(" -Filtering duplicated variants based on rsID...", verbose=verbose)
        if rsid in sumstats.columns:
            # Use keep parameter
            is_dup_rsid = sumstats.duplicated(subset=[rsid], keep=keep)
            # Exclude NA values
            is_dup_rsid = is_dup_rsid & sumstats[rsid].notna()
            is_duplicate = is_duplicate | is_dup_rsid
            log.write(" -Found {} duplicate variants based on rsID...".format(is_dup_rsid.sum()), verbose=verbose)
    
    # Filter by duplicated variants by CHR:POS:NEA:EA
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and (nea in sumstats.columns) and (ea in sumstats.columns) and ("d" in mode or "c" in mode):
        log.write(" -Filtering duplicated variants based on CHR,POS,EA and NEA...", verbose=verbose)
        # Check for any NA values in the subset columns
        subset_cols = [chrom, pos, ea, nea]
        has_na = sumstats[subset_cols].isna().any(axis=1)
        # Use keep parameter, excluding rows with NA
        is_dup_chrpos = sumstats.duplicated(subset=[chrom, pos, ea, nea], keep=keep)
        is_dup_chrpos = is_dup_chrpos & ~has_na
        is_duplicate = is_duplicate | is_dup_chrpos
        log.write(" -Found {} duplicate variants based on CHR,POS,EA and NEA...".format(is_dup_chrpos.sum()), verbose=verbose)
    
    # Filter by multiallelic variants by CHR:POS
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and "m" in mode:
        log.write(" -Filtering multiallelic variants based on CHR:POS...", verbose=verbose)
        # Check for any NA values in the subset columns
        subset_cols = [chrom, pos]
        has_na = sumstats[subset_cols].isna().any(axis=1)
        # Use keep parameter, excluding rows with NA
        is_dup_multiallelic = sumstats.duplicated(subset=[chrom, pos], keep=keep)
        is_dup_multiallelic = is_dup_multiallelic & ~has_na
        is_duplicate = is_duplicate | is_dup_multiallelic
        log.write(" -Found {} multiallelic variants based on CHR:POS...".format(is_dup_multiallelic.sum()), verbose=verbose)
    
    # Filter to keep only duplicates
    sumstats = sumstats.loc[is_duplicate, :]
    after_number = len(sumstats)
    
    log.write(" -Found {} duplicate/multiallelic variants in total...".format(after_number), verbose=verbose)
    
    # Remove NAs if requested
    if "n" in mode or remove_na == True:
        log.write(" -Removing NAs...", verbose=verbose)
        pre_number = len(sumstats)
        specified_columns = []
        if "d" in mode:
            if rsid in sumstats.columns: specified_columns.append(rsid)
            if snpid in sumstats.columns: specified_columns.append(snpid)
            if chrom in sumstats.columns: specified_columns.append(chrom)
            if pos in sumstats.columns: specified_columns.append(pos)
            if ea in sumstats.columns: specified_columns.append(ea)
            if nea in sumstats.columns: specified_columns.append(nea)
        if "r" in mode:
            if rsid in sumstats.columns: specified_columns.append(rsid)
        if "s" in mode:
            if snpid in sumstats.columns: specified_columns.append(snpid)
        if "m" in mode:
            if chrom in sumstats.columns: specified_columns.append(chrom)
            if pos in sumstats.columns: specified_columns.append(pos)
        if "c" in mode:
            if chrom in sumstats.columns: specified_columns.append(chrom)
            if pos in sumstats.columns: specified_columns.append(pos)
            if ea in sumstats.columns: specified_columns.append(ea)
            if nea in sumstats.columns: specified_columns.append(nea)
        specified_columns = list(set(specified_columns))
        if specified_columns:
            sumstats = sumstats.loc[~sumstats[specified_columns].isna().any(axis=1), :]
            after_number = len(sumstats)
            log.write(" -Removed {} variants with NA values in {}...".format(pre_number - after_number, specified_columns), verbose=verbose)
    
    # Sort the coordinates
    if keep_col is not None:
        log.write(" -Sorting the coordinates based on CHR and POS...", verbose=verbose)
        sumstats = _sort_coordinate(sumstats, verbose=False)
    
    gc.collect()
    # Update sumstats_obj.data with filtered result only if called with Sumstats object
    if not isinstance(sumstats_obj, pd.DataFrame):
        sumstats_obj.data = sumstats
    
    return sumstats
