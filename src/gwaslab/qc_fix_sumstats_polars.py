import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from liftover import get_lifter
from liftover import ChainFile
from functools import partial
from gwaslab.g_vchange_status_polars import vchange_statusp
from gwaslab.g_vchange_status import status_match
from gwaslab.g_vchange_status import change_status
from gwaslab.g_Log import Log
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_chr_list
from gwaslab.qc_check_datatype import check_datatype
from gwaslab.qc_check_datatype import check_dataframe_shape
from gwaslab.qc_build import _process_build
from gwaslab.qc_build import _set_build
from gwaslab.g_version import _get_version
from gwaslab.util_in_fill_data import _convert_betase_to_mlog10p
from gwaslab.util_in_fill_data import _convert_betase_to_p
from gwaslab.util_in_fill_data import _convert_mlog10p_to_p
from gwaslab.bd_common_data import get_chain
import polars as pl
###############################################################################################################
# 20220426
def get_reverse_complementary_allele(a):
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a[::-1].translate(dic)

def flip_direction(string):
    flipped_string=""
    for char in string:
        if char=="?":
            flipped_string+="?"
        elif char=="+":
            flipped_string+="-"
        elif char=="-":
            flipped_string+="+"
        else: #sometime it is 0
            flipped_string+=char
    return flipped_string

def flip_by_swap(sumstats, matched_index, log, verbose):
    if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
        log.write(" -Swapping column: NEA <=> EA...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(EA = pl.col("EA"),
                                                               NEA= pl.col("NEA"))
    return sumstats

def flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1):
    if "OR" in sumstats.columns:
        log.write(" -Flipping column: OR = 1 / OR...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(OR =  factor / pl.col("OR"))
    if "OR_95L" in sumstats.columns:
        log.write(" -Flipping column: OR_95U = 1 / OR_95L...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(OR_95L =  factor / pl.col("OR_95L"))
    if "OR_95U" in sumstats.columns:
        log.write(" -Flipping column: OR_95L = 1 / OR_95U...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(OR_95U =  factor / pl.col("OR_95U"))
    if "HR" in sumstats.columns:
        log.write(" -Flipping column: HR = 1 / HR...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(HR =  factor / pl.col("HR"))
    if "HR_95L" in sumstats.columns:
        log.write(" -Flipping column: HR_95U = 1 / HR_95L...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(HR_95L =  factor / pl.col("HR_95L"))
    if "HR_95U" in sumstats.columns:
        log.write(" -Flipping column: HR_95L = 1 / HR_95U...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(HR_95U =  factor / pl.col("HR_95U"))
    return sumstats

def flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1):
    if "BETA" in sumstats.columns:
        log.write(" -Flipping column: EAF = 1 - EAF...", verbose=verbose) 
        sumstats = sumstats.filter(matched_index).with_columns(EAF = 1 - pl.col("EAF"))
    return sumstats

def flip_by_sign(sumstats, matched_index, log, verbose, cols=None):
    if "BETA" in sumstats.columns:
        sumstats = sumstats.filter(matched_index).with_columns(BETA = - pl.col("BETA"))
    if "BETA_95L" in sumstats.columns:
        sumstats = sumstats.filter(matched_index).with_columns(BETA_95L = - pl.col("BETA_95L"))
    if "BETA_95U" in sumstats.columns:
        sumstats = sumstats.filter(matched_index).with_columns(BETA_95U = - pl.col("BETA_95U"))
    if "T" in sumstats.columns:
        sumstats = sumstats.filter(matched_index).with_columns(T = - pl.col("T"))
    if "Z" in sumstats.columns:
        sumstats = sumstats.filter(matched_index).with_columns(Z = - pl.col("Z"))
    if "DIRECTION" in sumstats.columns:
        sumstats = sumstats.filter(matched_index).with_columns(BETA = - pl.col("DIRECTION").map_batches(lambda x: pl.Series(flip_direction(x))))
    return sumstats

def flipallelestatsp(sumstats,status="STATUS",verbose=True,log=Log()):
    ##start function with col checking#########################################################

    if_stats_flipped = False
    ###################get reverse complementary####################
    pattern = r"\w\w\w\w\w[45]\w"  
    #matched_index = status_match(sumstats[status],6,[4,5]) #
    #matched_index = sumstats[status].str[5].str.match(r"4|5")

    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w\w[45]\w")

    if len(sumstats.filter(matched_index))>0:
        log.write("Start to convert alleles to reverse complement for SNPs with status xxxxx[45]x...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            log.write(" -Converting to reverse complement : EA and NEA...", verbose=verbose) 
            
            sumstats = sumstats.filter(matched_index).with_columns(
                NEA = pl.col("NEA").map_batches(lambda x: pl.Series(get_reverse_complementary_allele(x))),
                EA = pl.col("EA").map_batches(lambda x: pl.Series(get_reverse_complementary_allele(x)))
                                                                   )
        
            sumstats  = vchange_statusp(sumstats, matched_index, status,6, "4","2")
            log.write(" -Changed the status for flipped variants : xxxxx4x -> xxxxx2x", verbose=verbose)
        if_stats_flipped = True
    
    ###################flip ref####################
    pattern = r"\w\w\w\w\w[35]\w"  
    #matched_index = status_match(sumstats[status],6,[3,5]) #sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w[35]\w\w")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        
        flip_by_swap(sumstats, matched_index, log, verbose)
        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants : xxxxx[35]x -> xxxxx[12]x", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status,6, "35","12")
        if_stats_flipped = True
        
    ###################flip ref for undistingushable indels####################
    pattern = r"\w\w\w\w[123][67]6"  
    #matched_index = status_match(sumstats[status],6,[1,2,3])|status_match(sumstats[status],6,[6,7])|status_match(sumstats[status],7,6) #sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w[123][67]6")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for standardized indels with status xxxx[123][67][6]: ALT->EA , REF->NEA...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        
        flip_by_swap(sumstats, matched_index, log, verbose)
        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status, 7, "6","4")
        if_stats_flipped = True
         # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    pattern = r"\w\w\w\w\w[012]5"  
    #matched_index = status_match(sumstats[status],6,[0,1,2]) | status_match(sumstats[status],7,[5])#sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w[05|15|25]")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 

        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status,7, "5","2")
        if_stats_flipped = True

    if if_stats_flipped != True:
        log.write(" -No statistics have been changed.")
    return sumstats