import gc
import pandas as pd
import polars as pl
import numpy as np
from gwaslab.g_Log import Log
# pandas.api.types.is_int64_dtype
# pandas.api.types.is_categorical_dtype

dtype_dict ={
    "SNPID":[pl.String()],
    "rsID": [pl.String()],
    "CHR":  [pl.Int64()],
    "POS":  [pl.Int64()],
    "EA":   [pl.String()],  
    "NEA":[pl.String()],  
    "REF":[pl.String()],  
    "ALT":[pl.String()],  
    "BETA":[pl.Float64()],
    "BETA_95L":[pl.Float64()],
    "BETA_95U":[pl.Float64()],
    "SE":[pl.Float64()],
    "N":[pl.Int64()],
    "N_CASE":[pl.Int64()],
    "N_CONTROL":[pl.Int64()],
    "OR":[pl.Float64()],
    "OR_95L":[pl.Float64()],
    "OR_95U":[pl.Float64()],
    "HR":[pl.Float64()],
    "HR_95L":[pl.Float64()],
    "HR_95U":[pl.Float64()],
    "P":[pl.Float64()],
    "MLOG10P":[pl.Float64()],
    "Z":[pl.Float64()],
    "F":[pl.Float64()],
    "T":[pl.Float64()],
    "TEST":[pl.String()],
    "CHISQ":[pl.Float64()],
    "I2":[pl.Float64()],
    "P_HET":[pl.Float64()],
    "SNPR2":[pl.Float64()],
    "EAF":[pl.Float64()],
    "NEAF":[pl.Float64()],
    "MAF":[pl.Float64()],
    "INFO":[pl.Float64()],
    "DOF":[pl.Int64()],  
    "STATUS":[pl.String()],  
    "DIRECTION":[pl.String()],
    'PIP'               :[pl.Float64()],
    'CREDIBLE_SET_INDEX':[pl.Int64()],
    'N_SNP'             :[pl.Int64()],
    'LOCUS'             :[pl.String()],
    'STUDY'             :[pl.String()],
    'BETA_RANDOM' :[pl.Float64()],
    'SE_RANDOM' :[pl.Float64()],
    'Z_RANDOM' :[pl.Float64()],
    'P_RANDOM' :[pl.Float64()]
    }

def check_datatype(sumstats, verbose=True, log=Log()):
    
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

def verify_datatype(header, dtype):

    if header in dtype_dict.keys():
        if dtype in dtype_dict[header]:
            return "T"
        else:
            return "F"
    else:
        return "NA"

def quick_convert_datatype(sumstats, log, verbose):
    for col in sumstats.columns:
        if col in dtype_dict.keys():
            if sumstats[col].dtype not in dtype_dict[col]:
                datatype=dtype_dict[col][0]
                log.write(" -Trying to convert datatype for {}: {} -> {}...".format(col, str(sumstats[col].dtype), datatype), end="" ,verbose=verbose)
                try:
                    sumstats = sumstats.cast({col: datatype})
                    log.write("{}".format(datatype),show_time=False, verbose=verbose)
                except:
                    log.write("Failed...",show_time=False,verbose=verbose)
                    pass
    return sumstats

def check_dataframe_shape(sumstats, log, verbose):
    memory_in_mb = sumstats.estimated_size(unit="mb") 
    try:
        log.write(" -Current Dataframe shape : {} x {} ; Memory usage: {:.2f} MB".format(len(sumstats),len(sumstats.columns),memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Dataframe shape...")
    
def check_dataframe_memory_usage(sumstats, log, verbose):
    memory_in_mb = sumstats.estimated_size(unit="mb") 
    try:
        log.write(" -Current Dataframe memory usage: {:.2f} MB".format(memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Memory usage...")

