from typing import TYPE_CHECKING, Optional, List, Tuple, Union, Dict, Any
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from pandas.api.types import CategoricalDtype
from gwaslab.info.g_vchange_status import copy_status
from gwaslab.info.g_vchange_status import vchange_status
from gwaslab.info.g_vchange_status import ensure_status_int
from gwaslab.qc.qc_fix_sumstats import _flip_allele_stats
from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_reserved_headers import DEFAULT_COLUMN_ORDER
from gwaslab.util.util_in_fill_data import _fill_data
from itertools import combinations

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _merge_mold_with_sumstats_by_chrpos(
    mold: pd.DataFrame,
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    ref_path: Optional[str] = None,
    add_raw_index: bool = False,
    stats_cols1: Optional[List[str]] = None,
    stats_cols2: Optional[List[str]] = None,
    windowsizeb: int = 10,
    log: Log = Log(),
    suffixes: Tuple[str, str] = ("_MOLD",""),
    merge_mode: str = "inner",
    verbose: bool = True,
    return_not_matched_mold: bool = False
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.DataFrame]]:
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    log.write("Start to merge sumstats...", verbose=verbose)
    if merge_mode=="outer":
        sumstats = sumstats.rename(columns={
                                            "SNPID":"_SNPID_RIGHT",
                                            "rsID":"_rsID_RIGHT"
                                            })
    
    # drop old ids
    cols_to_drop = []
    for i in sumstats.columns:
        if i in ["SNPID","rsID"]:
            cols_to_drop.append(i)    
    if len(cols_to_drop)>0:
        log.write(" -Dropping old IDs:{}".format(cols_to_drop), verbose=verbose)
        sumstats = sumstats.drop(columns=cols_to_drop)


    if add_raw_index==True:
        index1= "_INDEX" + suffixes[0]
        index2= "_INDEX" + suffixes[1]
        mold[index1] = mold.index
        sumstats[index2] =  sumstats.index
        

    if ref_path is not None :
        # index for checking removed variants
        index1= "_INDEX" + suffixes[0]
        index2= "_INDEX" + suffixes[1]
        mold[index1] = range(len(mold))
        sumstats[index2] = range(len(sumstats))
    
    #if return_not_matched_mold:
    #   mold["_IDENTIFIER_FOR_VARIANT"] = range(len(mold))
    #   sumstats["_IDENTIFIER_FOR_VARIANT2"] = range(len(sumstats))

    # mold sumffix + mold 
    mold_sumstats = pd.merge(mold, sumstats, on=["CHR","POS"], how=merge_mode,suffixes=suffixes)

    if merge_mode=="outer":
        is_temp_na = mold_sumstats["EA_1"].isna()
        log.write(" -Detected {} variants not in the template...".format(sum(is_temp_na)), verbose=verbose)
        
        mold_sumstats["EA_1"] = mold_sumstats["EA_1"].astype("string")
        mold_sumstats["NEA_1"] = mold_sumstats["NEA_1"].astype("string")
        mold_sumstats["EA"] = mold_sumstats["EA"].astype("string")
        mold_sumstats["NEA"] = mold_sumstats["NEA"].astype("string")

        # for variants not in template, copy snp info
        mold_sumstats.loc[is_temp_na, ["SNPID","EA_1","NEA_1","STATUS_1"]] = mold_sumstats.loc[is_temp_na, ["_SNPID_RIGHT","EA","NEA","STATUS"]].values
        
        # 
        if "_rsID_RIGHT" in mold_sumstats.columns:
            mold_sumstats.loc[is_temp_na, "rsID"] = mold_sumstats.loc[is_temp_na, "_rsID_RIGHT"].values
        
        
        # for variants not in right sumstats, copy snp info
        is_temp_na_2 = mold_sumstats["EA"].isna()
        mold_sumstats.loc[is_temp_na_2, ["EA","NEA"]] = mold_sumstats.loc[is_temp_na_2, ["EA_1","NEA_1"]].values
        mold_sumstats = mold_sumstats.drop(columns=["_SNPID_RIGHT"])

    log.write(" -After merging by CHR and POS:{}".format(len(mold_sumstats)), verbose=verbose)
    
    mold_sumstats = _keep_variants_with_same_allele_set(mold_sumstats,suffixes=suffixes)

    log.write(" -Matched variants:{}".format(len(mold_sumstats)), verbose=verbose)
    
    #if ref_path is not None:
    #    # match removed sumstats
    #    mold_removed = mold.loc[~mold[index1].isin(mold_sumstats[index1]),:]
    #    iron_removed = sumstats.loc[~sumstats[index2].isin(mold_sumstats[index2]),:]
    #    _match_two_sumstats(mold_removed,iron_removed,ref_path,windowsizeb=windowsizeb)
    #    mold_sumstats.drop(columns=["_INDEX",""])
    
    if return_not_matched_mold == True:

        sumstats1 = mold.loc[~mold["_RAW_INDEX_1"].isin(mold_sumstats["_RAW_INDEX_1"]),:]
        sumstats1 = sumstats1.drop(columns=["_RAW_INDEX_1"])
        sumstats1 = _renaming_cols_r(sumstats1, stats_cols1 +["EA","NEA"],suffix="_1", verbose=False)
        
        sumstats2 = sumstats.loc[~sumstats["_RAW_INDEX_2"].isin(mold_sumstats["_RAW_INDEX_2"]),:]
        sumstats2 = sumstats2.drop(columns=["_RAW_INDEX_2"])

        mold_sumstats= mold_sumstats.drop(columns=["_RAW_INDEX_1","_RAW_INDEX_2"])
        
        return mold_sumstats, sumstats1, sumstats2
    
    return mold_sumstats

def _keep_variants_with_same_allele_set(
    sumstats: pd.DataFrame,
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_MOLD","")
) -> pd.DataFrame:

    ea1="EA"+suffixes[0]
    nea1="NEA"+suffixes[0]
    ea2="EA"+suffixes[1]
    nea2="NEA"+suffixes[1]

    all_alleles = set(list(sumstats[ea1].unique())+list(sumstats[nea1].unique())+list(sumstats[ea2].unique())+list(sumstats[nea2].unique()))
    allele_type = CategoricalDtype(categories=all_alleles, ordered=False)
    sumstats[[nea1,ea1,nea2,ea2]] = sumstats[[nea1,ea1,nea2,ea2]].astype(allele_type)
    
    is_perfect_match = (sumstats[ea2] == sumstats[ea1]) & (sumstats[nea2] == sumstats[nea1])
    is_flipped_match = (sumstats[ea2] == sumstats[nea1]) & (sumstats[nea2] == sumstats[ea1])
    is_allele_set_match = is_flipped_match | is_perfect_match
    
    log.write(" -Matching alleles and keeping only variants with same allele set: ", verbose=verbose)
    log.write("  -Perfect match: {}".format(sum(is_perfect_match)), verbose=verbose)
    log.write("  -Flipped match: {}".format(sum(is_flipped_match)), verbose=verbose)
    log.write("  -Unmatched : {}".format(sum(~is_allele_set_match)), verbose=verbose)
    
    return sumstats.loc[is_allele_set_match,:]

def _align_with_mold(
    sumstats: pd.DataFrame,
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_MOLD","")
) -> pd.DataFrame:
    
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
    
    # Ensure status columns are integer type before assignment
    sumstats = ensure_status_int(sumstats, status2)
    
    log.write("  -For perfect match: copy STATUS from reference...", verbose=verbose)
    sumstats.loc[is_perfect_match,status2] = copy_status(sumstats.loc[is_perfect_match,status1], sumstats.loc[is_perfect_match,status2],6)
    
    log.write("  -For Flipped match: convert STATUS xxxxx[456789]x to xxxxx3x...", verbose=verbose)
    sumstats.loc[is_flipped_match,status2] = vchange_status(sumstats.loc[is_flipped_match,status2],6,"456789","333333")
    
    return sumstats

def _fill_missing_columns(
    sumstats: pd.DataFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True
) -> pd.DataFrame:
    sumstats = _fill_data(sumstats, to_fill=columns)
    return sumstats

def _renaming_cols(
    sumstats: pd.DataFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_1","_2")
) -> pd.DataFrame:
    to_rename =["STATUS"]
    for col in columns:
        if col in sumstats.columns:
            to_rename.append(col)
    sumstats = sumstats.rename(columns={i:i + suffixes[1] for i in to_rename})
    log.write(" -Renaming sumstats2 columns by adding suffix {}".format(suffixes[1]),verbose=verbose)
    return sumstats

def _renaming_cols_r(
    sumstats: pd.DataFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True,
    suffix: str = ""
) -> pd.DataFrame:
    # columns: name without suffix
    to_rename =[]
    for col in columns:
        if col + suffix in sumstats.columns:
            to_rename.append(col)
    sumstats = sumstats.rename(columns={i + suffix:i for i in to_rename})
    log.write(" -Renaming sumstats columns by removing suffix {}".format(suffix),verbose=verbose)
    return sumstats

def _sort_pair_cols(
    molded_sumstats: pd.DataFrame,
    verbose: bool = True,
    log: Log = Log(),
    order: Optional[List[str]] = None,
    stats_order: Optional[List[str]] = None,
    suffixes: Tuple[str, str] = ("_1","_2")
) -> pd.DataFrame:
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

def _check_daf(sumstats, log=Log(),verbose=True,suffixes=("_MOLD","")):
    sumstats["DAF"] = sumstats["EAF"+suffixes[0]] - sumstats["EAF"]
    return sumstats

def _assign_warning_code(
    sumstats: pd.DataFrame,
    threshold: float = 0.2,
    log: Log = Log(),
    verbose: bool = True
) -> pd.DataFrame:
    is_outlier = (sumstats["DAF"] > threshold) | (-sumstats["DAF"] > threshold)
    sumstats["WARNING"] = 0
    sumstats.loc[is_outlier, "WARNING"] = 1

    log.write(" -Detected variants with large DAF with threshold of {} : {}".format(threshold,sum(is_outlier)), verbose=verbose)

    return sumstats


#def _match_two_sumstats(mold,sumstats,ref_path,windowsizeb=25,verbose=True,log=Log()):
#    
#    from gwaslab.io.io_fasta import load_fasta_auto
#    records = load_fasta_auto(ref_path, as_seqrecord=True)
#
#    chromlist = list(set(mold["CHR"].values) & set(sumstats["CHR"].values))
#    
#    for record in records:
#        if len(chromlist) ==0:
#            break
#
#        if record is not None:
#            ##############################################################################
#            record_chr = int(str(record.id).strip("chrCHR").upper())
#            
#            if record_chr in chromlist:
#                log.write(record_chr," ", end="",show_time=False,verbose=verbose) 
#                chromlist.remove(record_chr)
#            else:
#                continue
#            ###############################################################################
#            mold_chr = mold.loc[mold["CHR"]==record_chr,:]
#            sumstats_chr = sumstats.loc[sumstats["CHR"]==record_chr,:]
#
#            for index, row in sumstats_chr.iterrows():
#                if len(row["EA"])>1 or len(row["NEA"])>1:
#                    is_in_variants_lista = (mold_chr["POS"] > row["POS"] - windowsizeb) & (mold_chr["POS"]< row["POS"] + windowsizeb) 
#
#                    is_in_variants_listb = (sumstats_chr["POS"] > row["POS"] - windowsizeb) & (sumstats_chr["POS"]< row["POS"] + windowsizeb) 
#                    
#                    if sum(is_in_variants_lista)>0 and sum(is_in_variants_listb)>0 and (sum(is_in_variants_lista) + sum(is_in_variants_listb) >2):
#                        variants_lista = mold.loc[is_in_variants_lista,:]
#                        variants_listb = sumstats.loc[is_in_variants_listb,:]
#                        
#                        refseq = record[row["POS"]-1 - windowsizeb: row["POS"] + windowsizeb].seq.upper()
#                        _match_single_variant(refseq, variants_lista, variants_listb, left_offset=row["POS"] - windowsizeb, windowsizeb=windowsizeb)
#
#def _match_single_variant(refseq,  variants_lista, variants_listb, left_offset,windowsizeb):
#    
#    
#    seta=set()
#    setb=set()
#
#    seta_pumutations=[]
#    for i in range(1, len(variants_lista)+1):
#        seta_pumutations+=combinations(variants_lista.index, i)
#
#    for i in seta_pumutations:
#        if _is_ref_overlap(variants_lista.loc[i,:],suffix="_MOLD"):
#            continue
#        else:
#            seta = _form_haplotype(refseq, variants_lista.loc[i,:], seta, left_offset,suffix="_MOLD")
#    
#    setb_pumutations=[]
#    for i in range(1,len(variants_listb)+1):
#        setb_pumutations+=combinations(variants_listb.index, i)
#    for i in setb_pumutations:
#        if _is_ref_overlap(variants_listb.loc[i,:],suffix=""):
#            continue
#        else:
#            setb = _form_haplotype(refseq, variants_listb.loc[i,:], setb, left_offset,suffix="")
#        
#    if len(seta & setb)>0:
#        print("-Topmed--------------------------------")  
#        print(variants_lista[["CHR","POS","NEA_MOLD","EA_MOLD","EAF_MOLD"]])
#        print("-Finngen--------------------------------")  
#        print(variants_listb[["CHR","POS","NEA","EA","EAF"]])
#        print(refseq,left_offset)
#        print("-set a--------------------------------")  
#        print(seta)
#        print("-set b---------------------------------")    
#        print(setb)   
#        print("------------------------------------")    
#        print("maybe equivalent ########################################################################")
#        a = seta & setb
#        for i in a:
#            print(i)
#
#def _is_ref_overlap(variants_list,suffix="_MOLD"):
#    previous_end = 0
#    for index, row in variants_list.iterrows():
#        if row["POS"] <= previous_end:
#            return True    
#        if row["POS"] + len(row["NEA"+suffix]) -1 > previous_end:
#            previous_end = row["POS"] + len(row["NEA"+suffix]) -1
#    return False
#
#def _form_haplotype(refseq, variants_list, haplotype_set, left_offset,suffix="_MOLD"):       
#        new_haplotype = ""
#        lastpos = 0
#        for index, row in variants_list.iterrows():
#            new_haplotype += refseq[lastpos:row["POS"] - left_offset]
#            new_haplotype += row["EA"+suffix]
#            lastpos = row["POS"] + len(row["NEA"+suffix])- left_offset
#        new_haplotype  += refseq[lastpos:]
#        haplotype_set.add(new_haplotype)
#        return haplotype_set