import gc
import pandas as pd
import numpy as np
from gwaslab.Log import Log

def check_datatype(sumstats, verbose=True, log=Log()):
    try:
        headers = []
        dtypes = []

        for header,dtype in sumstats.dtypes.items():
            width = max(len(header),len(str(dtype)))
            header_fix_length = header + " "*(width- len(header) )
            dtype_fix_length  = str(dtype) + " "*(width- len(str(dtype)))
            headers.append(format(header_fix_length))
            dtypes.append((str(dtype_fix_length)))
            
        if verbose: log.write(" -Column:", " ".join(headers))
        if verbose: log.write(" -DType :", " ".join(dtypes))
    except:
        pass