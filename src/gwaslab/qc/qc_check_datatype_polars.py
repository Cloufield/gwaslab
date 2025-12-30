from typing import Any
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
    #except:
    #    pass

def verify_datatype(header: str, dtype: Any) -> str:

    if header in dtype_dict.keys():
        if dtype in dtype_dict[header]:
            return "T"
        else:
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

