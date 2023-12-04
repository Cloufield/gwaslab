import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from pandas.api.types import CategoricalDtype
from gwaslab.g_vchange_status import copy_status
from gwaslab.g_vchange_status import vchange_status
from gwaslab.qc_fix_sumstats import flipallelestats
from gwaslab.util_in_fill_data import filldata
from Bio import SeqIO
from itertools import combinations

def _merge_mold_with_sumstats(mold, sumstats, ref_path=None, windowsizeb=10, log=Log(),verbose=True):
    cols_to_drop = []
    for i in sumstats.columns:
        if i in ["SNPID","rsID"]:
            cols_to_drop.append(i)
    
    if len(cols_to_drop)>0:
        log.write("Dropping old IDs:{}".format(cols_to_drop))
        sumstats = sumstats.drop(columns=cols_to_drop)
    
    mold["_INDEX"] = range(len(mold))
    sumstats["_INDEX"] = range(len(sumstats))
    
    mold_sumstats = pd.merge(mold, sumstats, on=["CHR","POS"], how="inner",suffixes=("_MOLD",""))
    log.write("After merging by CHR and POS:{}".format(len(mold_sumstats)))
    mold_sumstats = _keep_variants_with_same_allele_set(mold_sumstats)
    log.write("Matched variants:{}".format(len(mold_sumstats)))
    
    if ref_path is not None:
        mold_removed = mold.loc[~mold["_INDEX"].isin(mold_sumstats["_INDEX_MOLD"]),:]
        iron_removed = sumstats.loc[~sumstats["_INDEX"].isin(mold_sumstats["_INDEX"]),:]
        _match_two_sumstats(mold_removed,iron_removed,ref_path,windowsizeb=windowsizeb)

    return mold_sumstats

def _keep_variants_with_same_allele_set(sumstats, log=Log(),verbose=True):

    all_alleles = set(list(sumstats["EA"].unique())+list(sumstats["NEA"].unique())+list(sumstats["EA_MOLD"].unique())+list(sumstats["NEA_MOLD"].unique()))
    allele_type = CategoricalDtype(categories=all_alleles, ordered=False)
    sumstats.loc[:, ["EA","EA_MOLD","NEA","NEA_MOLD"]] = sumstats.loc[:, ["EA","EA_MOLD","NEA","NEA_MOLD"]].astype(allele_type)
    
    is_perfect_match = (sumstats["EA"] == sumstats["EA_MOLD"]) & (sumstats["NEA"] == sumstats["NEA_MOLD"])
    is_flipped_match = (sumstats["EA"] == sumstats["NEA_MOLD"]) & (sumstats["NEA"] == sumstats["EA_MOLD"])
    is_allele_set_match = is_flipped_match | is_perfect_match
    
    sumstats.loc[~is_allele_set_match,:]

    return sumstats.loc[is_allele_set_match,:]

def _align_with_mold(sumstats, log=Log(),verbose=True):
    is_perfect_match = (sumstats["EA"] == sumstats["EA_MOLD"]) & (sumstats["NEA"] == sumstats["NEA_MOLD"])
    is_flipped_match = (sumstats["EA"] == sumstats["NEA_MOLD"]) & (sumstats["NEA"] == sumstats["EA_MOLD"])
    sumstats.loc[is_perfect_match,"STATUS"] = copy_status(sumstats.loc[is_perfect_match,"STATUS_MOLD"], sumstats.loc[is_perfect_match,"STATUS"],6)
    sumstats.loc[is_flipped_match,"STATUS"] = vchange_status(sumstats.loc[is_flipped_match,"STATUS"],6,"456789","333333")
    return sumstats

def _fill_missing_columns(sumstats, columns, log=Log(),verbose=True):
    sumstats = filldata(sumstats, to_fill=columns)
    return sumstats

def _check_daf(sumstats, log=Log(),verbose=True):
    sumstats["DAF"] = sumstats["EAF_MOLD"] - sumstats["EAF"]
    return sumstats

def _assign_warning_code(sumstats, threshold=0.2, log=Log(),verbose=True):
    is_outlier = (sumstats["DAF"] > threshold) | (-sumstats["DAF"] > threshold)
    sumstats["WARNING"] = 0
    sumstats.loc[is_outlier, "WARNING"] = 1

    log.write(" -Detected variants with large DAF with threshold of {} : {}".format(threshold,sum(is_outlier)), verbose=verbose)

    return sumstats


def _match_two_sumstats(mold,sumstats,ref_path,windowsizeb=25,verbose=True,log=Log()):
    
    records = SeqIO.parse(ref_path, "fasta")

    chromlist = list(set(mold["CHR"].values) & set(sumstats["CHR"].values))
    
    for record in records:
        if len(chromlist) ==0:
            break

        if record is not None:
            ##############################################################################
            record_chr = int(str(record.id).strip("chrCHR").upper())
            
            if record_chr in chromlist:
                if verbose:  log.write(record_chr," ", end="",show_time=False) 
                chromlist.remove(record_chr)
            else:
                continue
            ###############################################################################
            mold_chr = mold.loc[mold["CHR"]==record_chr,:]
            sumstats_chr = sumstats.loc[sumstats["CHR"]==record_chr,:]

            for index, row in sumstats_chr.iterrows():
                if len(row["EA"])>1 or len(row["NEA"])>1:
                    is_in_variants_lista = (mold_chr["POS"] > row["POS"] - windowsizeb) & (mold_chr["POS"]< row["POS"] + windowsizeb) 

                    is_in_variants_listb = (sumstats_chr["POS"] > row["POS"] - windowsizeb) & (sumstats_chr["POS"]< row["POS"] + windowsizeb) 
                    
                    if sum(is_in_variants_lista)>0 and sum(is_in_variants_listb)>0 and (sum(is_in_variants_lista) + sum(is_in_variants_listb) >2):
                        variants_lista = mold.loc[is_in_variants_lista,:]
                        variants_listb = sumstats.loc[is_in_variants_listb,:]
                        
                        refseq = record[row["POS"]-1 - windowsizeb: row["POS"] + windowsizeb].seq.upper()
                        _match_single_variant(refseq, variants_lista, variants_listb, left_offset=row["POS"] - windowsizeb, windowsizeb=windowsizeb)

def _match_single_variant(refseq,  variants_lista, variants_listb, left_offset,windowsizeb):
    
    
    seta=set()
    setb=set()

    seta_pumutations=[]
    for i in range(1, len(variants_lista)+1):
        seta_pumutations+=combinations(variants_lista.index, i)

    for i in seta_pumutations:
        if _is_ref_overlap(variants_lista.loc[i,:],suffix="_MOLD"):
            continue
        else:
            seta = _form_haplotype(refseq, variants_lista.loc[i,:], seta, left_offset,suffix="_MOLD")
    
    setb_pumutations=[]
    for i in range(1,len(variants_listb)+1):
        setb_pumutations+=combinations(variants_listb.index, i)
    for i in setb_pumutations:
        if _is_ref_overlap(variants_listb.loc[i,:],suffix=""):
            continue
        else:
            setb = _form_haplotype(refseq, variants_listb.loc[i,:], setb, left_offset,suffix="")
        
    if len(seta & setb)>0:
        print("-Topmed--------------------------------")  
        print(variants_lista[["CHR","POS","NEA_MOLD","EA_MOLD","EAF_MOLD"]])
        print("-Finngen--------------------------------")  
        print(variants_listb[["CHR","POS","NEA","EA","EAF"]])
        print(refseq,left_offset)
        print("-set a--------------------------------")  
        print(seta)
        print("-set b---------------------------------")    
        print(setb)   
        print("------------------------------------")    
        print("maybe equivalent ########################################################################")
        a = seta & setb
        for i in a:
            print(i)

def _is_ref_overlap(variants_list,suffix="_MOLD"):
    previous_end = 0
    for index, row in variants_list.iterrows():
        if row["POS"] <= previous_end:
            return True    
        if row["POS"] + len(row["NEA"+suffix]) -1 > previous_end:
            previous_end = row["POS"] + len(row["NEA"+suffix]) -1
    return False

def _form_haplotype(refseq, variants_list, haplotype_set, left_offset,suffix="_MOLD"):       
        new_haplotype = ""
        lastpos = 0
        for index, row in variants_list.iterrows():
            new_haplotype += refseq[lastpos:row["POS"] - left_offset]
            new_haplotype += row["EA"+suffix]
            lastpos = row["POS"] + len(row["NEA"+suffix])- left_offset
        new_haplotype  += refseq[lastpos:]
        haplotype_set.add(new_haplotype)
        return haplotype_set