import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
# pandas.api.types.is_int64_dtype
# pandas.api.types.is_categorical_dtype

dtype_dict ={
    "SNPID":["string","object"],
    "rsID":["string","object"],
    "CHR":["Int64","int64","int32","Int32","int"],
    "POS":["int64","Int64"],
    "EA":["category"],  
    "NEA":["category"],  
    "REF":["category"],  
    "ALT":["category"],  
    "BETA":["float64"],
    "BETA_95L":["float64"],
    "BETA_95U":["float64"],
    "SE":["float64"],
    "N":["Int64","int64","int32","Int32","int"],
    "N_CASE":["Int64","int64","int32","Int32","int"],
    "N_CONTROL":["Int64","int64","int32","Int32","int"],
    "OR":["float64"],
    "OR_95L":["float64"],
    "OR_95U":["float64"],
    "HR":["float64"],
    "HR_95L":["float64"],
    "HR_95U":["float64"],
    "P":["float64"],
    "MLOG10P":["float64"],
    "Z":["float64"],
    "F":["float64"],
    "T":["float64"],
    "TEST":["string","object","category"],
    "CHISQ":["float64"],
    "I2":["float64"],
    "PHET":["float64"],
    "SNPR2":["float64"],
    "EAF":["float64","float","float32"],
    "NEAF":["float64","float","float32"],
    "MAF":["float64","float","float32"],
    "INFO":["float64","float","float32"],
    "DOF":["Int64","int64","int32","Int32","int"],  
    "STATUS":["category"],  
    "DIRECTION":["string","object"],
    'PIP'               :["float64","float","float32"],
    'CREDIBLE_SET_INDEX':["Int64","int64","int32","Int32","int"],
    'N_SNP'             :["Int64","int64","int32","Int32","int"],
    'LOCUS'             :["string","object","category"],
    'STUDY'             :["string","object","category"]
    }

def check_datatype(sumstats, verbose=True, log=Log()):
    
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
    except:
        pass

def verify_datatype(header, dtype):

    if header in dtype_dict.keys():
        if str(dtype) in dtype_dict[header]:
            return "T"
        else:
            return "F"
    else:
        return "NA"

def quick_convert_datatype(sumstats, log, verbose):
    for col in sumstats.columns:
        if col in dtype_dict.keys():
            if str(sumstats[col].dtypes) not in dtype_dict[col]:
                datatype=dtype_dict[col][0]
                log.write(" -Trying to convert datatype for {}: {} -> {}...".format(col, str(sumstats[col].dtypes), datatype), end="" ,verbose=verbose)
                try:
                    sumstats[col] = sumstats[col].astype(datatype)
                    log.write("{}".format(datatype),show_time=False, verbose=verbose)
                except:
                    log.write("Failed...",show_time=False,verbose=verbose)
                    pass
    return sumstats



def check_dataframe_shape(sumstats, log, verbose):
    memory_in_mb = sumstats.memory_usage().sum()/1024/1024
    try:
        log.write(" -Current Dataframe shape : {} x {} ; Memory usage: {:.2f} MB".format(len(sumstats),len(sumstats.columns),memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Dataframe shape...")
    
def check_dataframe_memory_usage(sumstats, log, verbose):
    memory_in_mb = sumstats.memory_usage().sum()/1024/1024
    try:
        log.write(" -Current Dataframe memory usage: {:.2f} MB".format(memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Memory usage...")

