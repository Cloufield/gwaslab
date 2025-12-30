from typing import TYPE_CHECKING, Optional, List, Tuple, Union, Dict, Any
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from pandas.api.types import CategoricalDtype
from gwaslab.info.g_vchange_status import copy_status
from gwaslab.info.g_vchange_status_polars import vchange_statusp
from gwaslab.info.g_vchange_status_polars import copy_statusp
from gwaslab.qc.qc_fix_sumstats import _flip_allele_stats
from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_reserved_headers import DEFAULT_COLUMN_ORDER
from gwaslab.util.util_in_fill_data import _fill_data
from itertools import combinations
import polars as pl

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _merge_mold_with_sumstats_by_chrposp(
    mold: pl.DataFrame,
    sumstats: pl.DataFrame,
    ref_path: Optional[str] = None,
    add_raw_index: bool = False,
    stats_cols1: Optional[List[str]] = None,
    stats_cols2: Optional[List[str]] = None,
    windowsizeb: int = 10,
    log: Log = Log(),
    suffixes: Tuple[str, str] = ("_MOLD",""),
    merge_mode: str = "full",
    merge_by_id: bool = False,
    verbose: bool = True,
    return_not_matched_mold: bool = False
) -> Union[pl.DataFrame, Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]]:
    
    log.write("Start to merge sumstats...", verbose=verbose)
    if merge_mode=="full":
        
        sumstats = sumstats.rename({
                                            "SNPID":"_SNPID_RIGHT",
                                            "rsID":"_rsID_RIGHT"
                                            }, strict=False) #,

    
    if merge_by_id==False:    
        # drop old ids
        cols_to_drop = []
        for i in sumstats.columns:
            if i in ["SNPID","rsID"]:
                cols_to_drop.append(i)    
        if len(cols_to_drop)>0:
            log.write(" -Dropping old IDs:{}".format(cols_to_drop), verbose=verbose)
            sumstats = sumstats.drop(columns=cols_to_drop)

        ##################################################################################################################

        # mold sumffix + mold 
            # add ASET
        mold = mold.with_columns(
                pl.when( pl.col("EA_1") > pl.col("NEA_1") )  
                .then(   pl.col("EA_1") + ":" + pl.col("NEA_1")  )  
                .otherwise( pl.col("NEA_1") + ":" + pl.col("EA_1") )
                .alias("ASET")
            )
        
        sumstats = sumstats.with_columns(
                pl.when( pl.col("EA") > pl.col("NEA") )  
                .then(   pl.col("EA") + ":" + pl.col("NEA")  )  
                .otherwise( pl.col("NEA") + ":" + pl.col("EA") )
                .alias("ASET"))
        
        sumstats_len = len(sumstats)
        mold_len = len(mold)
        sumstats = sumstats.unique(subset=["CHR","POS","ASET"]) 
        mold = mold.unique(subset=["CHR","POS","ASET"]) 

        log.write(f' -Left:  dropping duplicated variants based on CHR,POS,ASET: {sumstats_len - len(sumstats)}')
        log.write(f' -Right: dropping duplicated variants based on CHR,POS,ASET: {mold_len - len(mold)}')

        mold = mold.with_columns(
            pl.when( pl.col("NEA_1").str.len_chars() != pl.col("EA_1").str.len_chars() )  
            .then(   
                pl.when( pl.col("EAF_1")<0.5 ).then(
                    pl.col("ASET") + ":" + pl.col("EA_1") 
                ).otherwise( pl.col("ASET") + ":" + pl.col("NEA_1") )
                .alias("ASET")
                )  
            .otherwise( pl.col("ASET") )
            .alias("ASET")
        )
        
        sumstats = sumstats.with_columns(
            pl.when( pl.col("NEA").str.len_chars() != pl.col("EA").str.len_chars() )  
            .then(   
                pl.when( pl.col("EAF")<0.5 ).then(
                    pl.col("ASET") + ":" + pl.col("EA") 
                ).otherwise( pl.col("ASET") + ":" + pl.col("NEA") )
                .alias("ASET")
                )  
            .otherwise( pl.col("ASET"))
            .alias("ASET")
            )

        mold_sumstats = mold.join(sumstats, on=["CHR","POS","ASET"], how=merge_mode, suffix="_", coalesce=True)
    
    elif merge_by_id==True:

        sumstats = sumstats.rename({
                                            "_SNPID_RIGHT":"SNPID",
                                            }, strict=False)
        

        sumstats_len = len(sumstats)
        mold_len = len(mold)
        sumstats = sumstats.unique(subset=["SNPID","CHR","POS"])
        mold = mold.unique(subset=["SNPID","CHR","POS"]) 
        log.write(f' -Left:  dropping duplicated variants based on CHR,POS,SNPID: {sumstats_len - len(sumstats)}')
        log.write(f' -Right: dropping duplicated variants based on CHR,POS,SNPID: {mold_len - len(mold)}')
        mold_sumstats = mold.join(sumstats, on=["SNPID","CHR","POS"], how=merge_mode, suffix="_", coalesce=True)


    if merge_mode=="full":
        is_temp_na = mold_sumstats["EA_1"].is_null()
        log.write(" -Detected {} variants not in the template...".format(sum(is_temp_na)), verbose=verbose)

        for i in ["EA_1","NEA_1","EA","NEA"]:
            mold_sumstats = mold_sumstats.with_columns(pl.col(i).cast(pl.String).alias(i))

        if merge_by_id==False:
            mold_sumstats = mold_sumstats.with_columns(
            pl.when( is_temp_na )  
                .then(   pl.col("_SNPID_RIGHT")  )  
                .otherwise( pl.col("SNPID") )
                .alias("SNPID")
            )
            mold_sumstats = mold_sumstats.drop(["_SNPID_RIGHT"])

        # for variants not in template, copy snp info
        mold_sumstats = mold_sumstats.with_columns(
            pl.when( is_temp_na )  
                .then( pl.col("EA")  )  
                .otherwise( pl.col("EA_1") )
                .alias("EA_1")
        ).with_columns(
            pl.when( is_temp_na )  
                .then( pl.col("NEA")  )  
                .otherwise( pl.col("NEA_1") )
                .alias("NEA_1")
        ).with_columns(
            pl.when( is_temp_na )  
                .then( pl.col("EAF")  )  
                .otherwise( pl.col("EAF_1"))
                .alias("EAF_1")
        ).with_columns(
            pl.when( is_temp_na )  
                .then( pl.col("STATUS")  )  
                .otherwise( pl.col("STATUS_1") )
                .alias("STATUS_1")
        )

        # 
        if "_rsID_RIGHT" in mold_sumstats.columns:
            mold_sumstats = mold_sumstats.with_columns(
                pl.when( is_temp_na )  
                .then(   pl.col("_rsID_RIGHT")  )  
                .otherwise( pl.col("rsID") )
                .alias("rsID")
                )
        
        
        # for variants not in right sumstats, copy snp info
        is_temp_na_2 = mold_sumstats["EA"].is_null()
        
        mold_sumstats = mold_sumstats.with_columns(
                pl.when( is_temp_na_2 )  
                .then(   pl.col("EA_1") )  
                .otherwise( pl.col("EA") )
                .alias("EA")
                ).with_columns(
                pl.when( is_temp_na_2 )  
                .then(   pl.col("NEA_1")  )  
                .otherwise( pl.col("NEA") )
                .alias("NEA")
                )
        
        
    if merge_by_id==False:
        mold_sumstats = mold_sumstats.unique(subset=["CHR","POS","ASET"])
        log.write(" -After merging by CHR, POS and ASET:{}".format(len(mold_sumstats)), verbose=verbose)
    else:
        mold_sumstats = mold_sumstats.unique(subset=["SNPID","CHR","POS"])
        log.write(" -After merging by SNPID, CHR and POS:{}".format(len(mold_sumstats)), verbose=verbose)

    mold_sumstats = _keep_variants_with_same_allele_setp(mold_sumstats,suffixes=suffixes)

    log.write(" -Matched variants:{}".format(len(mold_sumstats)), verbose=verbose)
    
    return mold_sumstats

def _keep_variants_with_same_allele_setp(
    sumstats: pl.DataFrame,
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_MOLD","")
) -> pl.DataFrame:

    ea1="EA"+suffixes[0]
    nea1="NEA"+suffixes[0]
    ea2="EA"+suffixes[1]
    nea2="NEA"+suffixes[1]
    
    is_perfect_match = (sumstats[ea2] == sumstats[ea1]) & (sumstats[nea2] == sumstats[nea1])
    is_flipped_match = (sumstats[ea2] == sumstats[nea1]) & (sumstats[nea2] == sumstats[ea1])
    
    log.write("  -Perfect match: {}".format(sum(is_perfect_match)), verbose=verbose)
    log.write("  -Flipped match: {}".format(sum(is_flipped_match)), verbose=verbose)
    return sumstats

def _align_with_moldp(
    sumstats: pl.DataFrame,
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_MOLD","")
) -> pl.DataFrame:
    
    ea1="EA"+suffixes[0]
    nea1="NEA"+suffixes[0]
    ea2="EA"+suffixes[1]
    nea2="NEA"+suffixes[1]
    status1="STATUS"+suffixes[0]
    status2="STATUS"+suffixes[1]

    is_perfect_match = (sumstats[ea2] == sumstats[ea1]) & (sumstats[nea2] == sumstats[nea1])
    is_flipped_match = (sumstats[ea2] == sumstats[nea1]) & (sumstats[nea2] == sumstats[ea1])
    
    log.write(" -Aligning alleles with reference: ", verbose=verbose)
    log.write("  -Perfect match: {}".format(sum(is_perfect_match)), verbose=verbose)
    log.write("  -Flipped match: {}".format(sum(is_flipped_match)), verbose=verbose)
    
    log.write("  -For perfect match: copy STATUS from reference...", verbose=verbose)

    sumstats  = copy_statusp(sumstats, is_perfect_match, status1, status2, 6)

    log.write("  -For Flipped match: convert STATUS xxxxx[456789]x to xxxxx3x...", verbose=verbose)
    
    sumstats  = vchange_statusp(sumstats, is_flipped_match, status2, 6,"456789","333333")
    
    return sumstats

def _fill_missing_columnsp(
    sumstats: pl.DataFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True
) -> pl.DataFrame:
    sumstats = _fill_data(sumstats, to_fill=columns)
    return sumstats

def _renaming_colsp(
    sumstats: pl.DataFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_1","_2")
) -> pl.DataFrame:
    to_rename =["STATUS"]
    for col in columns:
        if col in sumstats.columns:
            to_rename.append(col)
    sumstats = sumstats.rename({i:i + suffixes[1] for i in to_rename})
    log.write(" -Renaming sumstats2 columns by adding suffix {}".format(suffixes[1]),verbose=verbose)
    return sumstats

def _renaming_cols_rp(
    sumstats: pl.DataFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True,
    suffix: str = ""
) -> pl.DataFrame:
    # columns: name without suffix
    to_rename =[]
    for col in columns:
        if col + suffix in sumstats.columns:
            to_rename.append(col)
    sumstats = sumstats.rename({i + suffix:i for i in to_rename})
    log.write(" -Renaming sumstats columns by removing suffix {}".format(suffix),verbose=verbose)
    return sumstats

def _sort_pair_colsp(
    molded_sumstats: pl.DataFrame,
    verbose: bool = True,
    log: Log = Log(),
    order: Optional[List[str]] = None,
    stats_order: Optional[List[str]] = None,
    suffixes: Tuple[str, str] = ("_1","_2")
) -> pl.DataFrame:
    if stats_order is None:
        # Base order: first 6 columns from DEFAULT_COLUMN_ORDER
        base_cols = ["SNPID","rsID", "CHR", "POS", "EA", "NEA"]
        # Stats order: remaining columns from DEFAULT_COLUMN_ORDER (excluding base columns)
        stats_order = [col for col in DEFAULT_COLUMN_ORDER if col not in base_cols]
        # Add non-reserved columns that may exist in sumstats
        additional_cols = ["DIRECTION", "I2"]  # I2 is an alias for I2_HET
        for col in additional_cols:
            if col not in stats_order:
                stats_order.append(col)
        order = base_cols.copy()
        
    for suffix in suffixes:
        for i in stats_order:
            order.append(i+suffix)
    
    log.write("Start to reorder the columns...",verbose=verbose)
    
    output_columns = []
    
    for i in order:
        if i in molded_sumstats.columns: 
            output_columns.append(i)
    for i in molded_sumstats.columns:
        if i not in order: 
            output_columns.append(i)
    
    log.write(" -Reordering columns to    :", ",".join(output_columns), verbose=verbose)
    molded_sumstats = molded_sumstats[ output_columns]
    log.write("Finished sorting columns successfully!", verbose=verbose)
    
    return molded_sumstats