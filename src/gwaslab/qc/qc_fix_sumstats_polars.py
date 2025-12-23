import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from functools import partial
from gwaslab.info.g_vchange_status_polars import vchange_statusp
from gwaslab.info.g_vchange_status import status_match
from gwaslab.info.g_vchange_status import change_status
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chain

from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_build import _set_build

from gwaslab.util.util_in_fill_data import _convert_betase_to_mlog10p
from gwaslab.util.util_in_fill_data import _convert_betase_to_p
from gwaslab.util.util_in_fill_data import _convert_mlog10p_to_p

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

        sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  pl.col("EA")  )  
                .otherwise( pl.col("NEA") )
                .alias("NEA"),

                pl.when( matched_index )  
                .then(   pl.col("NEA")  )  
                .otherwise( pl.col("EA") )
                .alias("EA"),
                )

    return sumstats

def flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1):
    for header in ["OR","OR_95L","OR_95U","HR","HR_95L","HR_95U"]:
        if header in sumstats.columns:
                log.write(" -Flipping column: {header} = 1 / {header}...".format(header = header), verbose=verbose) 
                sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  1/ pl.col(header) )  
                .otherwise( pl.col(header) )
                .alias(header)
                )
    return sumstats

def flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1):
    header="EAF"
    if header in sumstats.columns:
        log.write(" -Flipping column: EAF = 1 - EAF...", verbose=verbose) 
        sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  1 - pl.col(header) )  
                .otherwise( pl.col(header) )
                .alias(header)
                )
    return sumstats

def flip_by_sign(sumstats, matched_index, log, verbose, cols=None):
    for header in ["BETA","BETA_95L","BETA_95U","T","Z"]:
        if header in sumstats.columns:
                log.write(" -Flipping column: {header} = - {header}...".format(header = header), verbose=verbose) 
                sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  - pl.col(header) )  
                .otherwise( pl.col(header) )
                .alias(header)
                )
    
    if "DIRECTION" in sumstats.columns:
        sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  pl.col("DIRECTION").map_batches(lambda x: pl.Series(flip_direction(x))) )  
                .otherwise( pl.col("DIRECTION") )
                .alias("DIRECTION")
                )
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
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w\w[35]\w")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        
        sumstats = flip_by_swap(sumstats, matched_index, log, verbose)
        sumstats = flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        sumstats = flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        sumstats = flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
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
        
        sumstats = flip_by_swap(sumstats, matched_index, log, verbose)
        sumstats = flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        sumstats = flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        sumstats = flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status, 7, "6","4")
        if_stats_flipped = True
         # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    pattern = r"\w\w\w\w\w[012]5"  
    #matched_index = status_match(sumstats[status],6,[0,1,2]) | status_match(sumstats[status],7,[5])#sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w\w[012]5")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 

        sumstats = flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        sumstats = flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        sumstats = flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status,7, "5","2")
        if_stats_flipped = True

    if if_stats_flipped != True:
        log.write(" -No statistics have been changed.")
    return sumstats