import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from liftover import get_lifter
from liftover import ChainFile
from functools import partial
from gwaslab.g_vchange_status import vchange_status
from gwaslab.g_Log import Log

def _process_build(build, log, verbose):
    if str(build).lower() in ["hg19","19","37","b37","grch37"]:
        log.write(" -Genomic coordinates are based on GRCh37/hg19...", verbose=verbose)
        final_build = "19"
    elif str(build).lower() in ["hg18","18","36","b36","grch36"]:
        log.write(" -Genomic coordinates are based on GRCh36/hg18...", verbose=verbose)
        final_build = "18"
    elif str(build).lower() in ["hg38","38","b38","grch38"]:
        log.write(" -Genomic coordinates are based on GRCh38/hg38...", verbose=verbose)
        final_build = "38"
    elif str(build).lower() in ["t2t","hs1","chm13","13"]:
        log.write(" -Genomic coordinates are based on T2T-CHM13...", verbose=verbose)
        final_build = "13"
    else:
        log.warning("Version of genomic coordinates is unknown...", verbose=verbose)
        final_build = "99"
    return final_build

def _set_build(sumstats, build="99", status="STATUS",verbose=True,log=Log()):
    build = _process_build(build,log=log,verbose=verbose)
    sumstats[status] = vchange_status(sumstats[status], 1, "139",build[0]*3)
    sumstats[status] = vchange_status(sumstats[status], 2, "89",build[1]*3)
    return sumstats, build

def _check_build(target_build, build="99", status="STATUS",verbose=True,log=Log()):
    target_build = _process_build(target_build,log=log,verbose=verbose)
    build = _process_build(build,log=log,verbose=verbose)
    if build == "99":
        raise ValueError("Sumstats build is unknown. Please run infer_build() or set_build()")
    
    if target_build == "99": 
        raise ValueError("Target build is unknown.")
    
    if build!=target_build:
        raise ValueError("Please make sure sumstats build is {}".format(target_build))
    else:
        log.write(" -Sumstats build matches target build")
    
    return True

    
