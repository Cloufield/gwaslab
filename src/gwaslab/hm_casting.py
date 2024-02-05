import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from pandas.api.types import CategoricalDtype
from gwaslab.g_vchange_status import copy_status
from gwaslab.g_vchange_status import vchange_status
from gwaslab.qc_fix_sumstats import flipallelestats
from gwaslab.qc_check_datatype import check_datatype
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.util_in_fill_data import filldata
from Bio import SeqIO
from itertools import combinations

def _merge_mold_with_sumstats_by_chrpos(mold, sumstats, ref_path=None, windowsizeb=10, log=Log(),suffixes=("_MOLD",""),verbose=True,return_not_matched_mold =False):
    

    cols_to_drop = []
    for i in sumstats.columns:
        if i in ["SNPID","rsID"]:
            cols_to_drop.append(i)

    log.write("Start to merge sumstats...", verbose=verbose)
    
    if len(cols_to_drop)>0:
        log.write(" -Dropping old IDs:{}".format(cols_to_drop), verbose=verbose)
        sumstats = sumstats.drop(columns=cols_to_drop)
    
    if ref_path is not None :
        # index for checking removed variants
        index1= "_INDEX" + suffixes[0]
        index2= "_INDEX" + suffixes[1]
        mold[index1] = range(len(mold))
        sumstats[index2] = range(len(sumstats))
    
    if return_not_matched_mold:
        mold["_IDENTIFIER_FOR_VARIANT"] = range(len(mold))

    # mold sumffix + mold 
    mold_sumstats = pd.merge(mold, sumstats, on=["CHR","POS"], how="inner",suffixes=suffixes)
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
        sumstats1 = mold.loc[~mold["_IDENTIFIER_FOR_VARIANT"].isin(mold_sumstats["_IDENTIFIER_FOR_VARIANT"]),:]
        sumstats1= sumstats1.drop(columns=["_IDENTIFIER_FOR_VARIANT"])
        mold_sumstats= mold_sumstats.drop(columns=["_IDENTIFIER_FOR_VARIANT"])
        return mold_sumstats, sumstats1
    
    return mold_sumstats

def _keep_variants_with_same_allele_set(sumstats, log=Log(),verbose=True,suffixes=("_MOLD","")):

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

def _align_with_mold(sumstats, log=Log(),verbose=True, suffixes=("_MOLD","")):
    
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
    sumstats.loc[is_perfect_match,status2] = copy_status(sumstats.loc[is_perfect_match,status1], sumstats.loc[is_perfect_match,status2],6)
    
    log.write("  -For Flipped match: convert STATUS xxxxx[456789]x to xxxxx3x...", verbose=verbose)
    sumstats.loc[is_flipped_match,status2] = vchange_status(sumstats.loc[is_flipped_match,status2],6,"456789","333333")
    
    return sumstats

def _fill_missing_columns(sumstats, columns, log=Log(),verbose=True):
    sumstats = filldata(sumstats, to_fill=columns)
    return sumstats

def _renaming_cols(sumstats, columns, log=Log(),verbose=True, suffixes=("_1","_2")):
    to_rename =["STATUS"]
    for col in columns:
        if col in sumstats.columns:
            to_rename.append(col)
    sumstats = sumstats.rename(columns={i:i + suffixes[1] for i in to_rename})
    log.write(" -Renaming sumstats2 columns by adding suffix {}".format(suffixes[1]),verbose=verbose)
    return sumstats

def _sort_pair_cols(molded_sumstats, verbose=True, log=Log(), order=None, stats_order=None,suffixes=("_1","_2")):
    if stats_order is None:
        order = ["SNPID","rsID", "CHR", "POS", "EA", "NEA"]
        stats_order = ["EAF", "MAF", "BETA", "SE","BETA_95L","BETA_95U", "Z",
        "CHISQ", "P", "MLOG10P", "OR", "OR_95L", "OR_95U","HR", "HR_95L", "HR_95U","INFO", "N","N_CASE","N_CONTROL","DIRECTION","I2","P_HET","DOF","SNPR2","STATUS"]
        
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
                log.write(record_chr," ", end="",show_time=False,verbose=verbose) 
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