import re
#import modin.pandas as pd
import pandas as pd
import numpy as np
from os import path
from gwaslab.bd_common_data import get_high_ld
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.g_Log import Log
from gwaslab.g_vchange_status import vchange_status
from gwaslab.qc_fix_sumstats import sortcoordinate
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.hm_harmonize_sumstats import is_palindromic

import gc
def filtervalues(sumstats,expr,remove=False,verbose=True,log=Log()):
    log.write("Start filtering values by condition:",expr, verbose=verbose)
    prenum = len(sumstats)
    sumstats = sumstats.query(expr,engine='python').copy()
    afternum = len(sumstats)
    log.write(" -Removing "+ str(prenum-afternum) +" variants not meeting the conditions:",expr, verbose=verbose)
    log.write("Finished filtering values.", verbose=verbose)
    gc.collect()
    return sumstats

def filterout(sumstats,interval={},lt={},gt={},eq={},remove=False,verbose=True,log=Log()):
    log.write("Start filtering values:", verbose=verbose)
    for key,threshold in gt.items():
        num = len(sumstats.loc[sumstats[key]>threshold,:])
        log.write(" -Removing "+ str(num) +" variants with "+key+" > "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]<threshold,:]
    for key,threshold in lt.items():
        num = len(sumstats.loc[sumstats[key]<threshold,:])
        log.write(" -Removing "+ str(num) +" variants with "+key+" < "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]>threshold,:]
    for key,threshold in eq.items():
        num = len(sumstats.loc[sumstats[key]==threshold,:])
        log.write(" -Removing "+ str(num) +" variants with "+key+" = "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]!=threshold,:]
    log.write("Finished filtering values.", verbose=verbose)
    gc.collect()
    return sumstats.copy()

def filterin(sumstats,lt={},gt={},eq={},remove=False,verbose=True,log=Log()):
    log.write("Start filtering values:", verbose=verbose)
    for key,threshold in gt.items():
        num = len(sumstats.loc[sumstats[key]>threshold,:])
        log.write(" -Keeping "+ str(num) +" variants with "+key+" > "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]>threshold,:]
    for key,threshold in lt.items():
        num = len(sumstats.loc[sumstats[key]<threshold,:])
        log.write(" -Keeping "+ str(num) +" variants with "+key+" < "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]<threshold,:]
    for key,threshold in eq.items():
        num = len(sumstats.loc[sumstats[key]==threshold,:])
        log.write(" -Keeping "+ str(num) +" variants with "+key+" = "+ str(threshold)+" ...", verbose=verbose)
        sumstats = sumstats.loc[sumstats[key]==threshold,:]
    log.write("Finished filtering values.", verbose=verbose)
    gc.collect()
    return sumstats.copy()

def filterregionin(sumstats,path=None, chrom="CHR",pos="POS", high_ld=False, build="19", verbose=True,log=Log()):
    sumstats = sortcoordinate(sumstats,verbose=verbose)
    log.write("Start to filter in variants if in intervals defined in bed files:", verbose=verbose)
    log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns), verbose=verbose)
    
    if high_ld is True:
        path = get_high_ld(build=build)
        log.write(" -Loading bed format file for hg"+build, verbose=verbose)

    else:
        log.write(" -Loading bed format file: " , path, verbose=verbose)
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    dic=get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr")
    bed[0]=bed[0].map(dic).astype("string")
    bed[0]=bed[0].astype("Int64")
    sumstats = sumstats.sort_values(["CHR","POS"])
    
    if len(bed)<100:
        log.write(" -Bed file < 100 lines: using pd IntervalIndex... ", verbose=verbose)
        for i in sumstats[chrom].unique():
            if sum(bed[0]==i)>0:
                interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==i,"tuple"])
                sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
            else:
                continue
    else:
        log.write(" -Bed file > 100 lines: using two pointers, please make files are all sorted... ", verbose=verbose)
        bed_num  =0
        bed_chr   =bed.iloc[bed_num,0]
        bed_left  =bed.iloc[bed_num,1]
        bed_right =bed.iloc[bed_num,2]
        
        sum_num=0
        sum_chr_col = sumstats.columns.get_loc(chrom)
        sum_pos_col = sumstats.columns.get_loc(pos)
        sum_ind_col = sumstats.columns.get_loc("bed_indicator")
        while bed_num<len(bed) and sum_num<len(sumstats):
            #sumstats variant chr < bed chr
            if sumstats.iat[sum_num,sum_chr_col]<bed_chr:
                # next variant
                sum_num+=1
                continue
            #sumstats variant chr > bed chr
            elif sumstats.iat[sum_num,sum_chr_col]>bed_chr:
                # next bed record
                bed_num+=1
                bed_chr  =bed.iat[bed_num,0]
                bed_left  = bed.iat[bed_num,1]
                bed_right  = bed.iat[bed_num,2]
                continue
            #sumstats variant chr == bed chr
            else:
                #sumstats variant pos < bed left
                if sumstats.iat[sum_num,sum_pos_col]<bed_left:
                    # next variant
                    sum_num+=1
                    continue
                #sumstats variant pos > bed right
                elif sumstats.iat[sum_num,sum_pos_col]>bed_right:
                    # next bed record
                    bed_num+=1
                    bed_chr  =bed.iat[bed_num,0]
                    bed_left  = bed.iat[bed_num,1]
                    bed_right  = bed.iat[bed_num,2]
                # bed left  < sumstats variant  pos < bed right
                else:
                    # set to true
                    sumstats.iat[sum_num,sum_ind_col]=True
                    # next variant
                    sum_num+=1
                           
    ## in
    
    sumstats = sumstats.loc[sumstats["bed_indicator"],:]
    log.write(" -Number of variants in the specified regions to keep:",sum(sumstats["bed_indicator"]), verbose=verbose)
    log.write(" -Number of variants removed:",sum(~sumstats["bed_indicator"]), verbose=verbose)
    sumstats = sumstats.drop(columns="bed_indicator")
    log.write("Finished filtering in variants.", verbose=verbose)
    gc.collect()
    return sumstats

def filterregionout(sumstats, path=None, chrom="CHR",pos="POS", high_ld=False, build="19", verbose=True,log=Log()):
    sumstats = sortcoordinate(sumstats,verbose=verbose)
    log.write("Start to filter out variants if in intervals defined in bed files:", verbose=verbose)
    log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns), verbose=verbose)
    if high_ld is True:
        path = get_high_ld(build=build)
        log.write(" -Loading bed format file for hg"+build, verbose=verbose)

    else:
        log.write(" -Loading bed format file: " , path, verbose=verbose)
            
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    
    dic=get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr")
    bed[0]=bed[0].map(dic).astype("string")
    bed[0]=bed[0].astype("Int64")
    
    if len(bed)<100:
        log.write(" -Bed file < 100 lines: using pd IntervalIndex... ", verbose=verbose)
        for i in sumstats[chrom].unique():
            if sum(bed[0]==i)>0:
                interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==i,"tuple"])
                sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
            else:
                continue
    else:
        log.write(" -Bed file > 100 lines: using two pointers, please make files are all sorted... ", verbose=verbose)
        bed_num  =0
        bed_chr  =bed.iloc[bed_num,0]
        bed_left  =bed.iloc[bed_num,1]
        bed_right =bed.iloc[bed_num,2]
        
        sum_num=0
        sum_chr_col = sumstats.columns.get_loc(chrom)
        sum_pos_col = sumstats.columns.get_loc(pos)
        sum_ind_col = sumstats.columns.get_loc("bed_indicator")
        while bed_num<len(bed) and sum_num<len(sumstats):
            if sumstats.iat[sum_num,sum_chr_col]<bed_chr:
                sum_num+=1
                continue
            elif sumstats.iat[sum_num,sum_chr_col]>bed_chr:
                bed_num+=1
                bed_chr  =bed.iat[bed_num,0]
                bed_left  = bed.iat[bed_num,1]
                bed_right  = bed.iat[bed_num,2]
                continue
            else:
                if sumstats.iat[sum_num,sum_pos_col]<bed_left:
                    sum_num+=1
                    continue
                elif sumstats.iat[sum_num,sum_pos_col]>bed_right:
                    bed_num+=1
                    bed_chr  =bed.iat[bed_num,0]
                    bed_left  = bed.iat[bed_num,1]
                    bed_right  = bed.iat[bed_num,2]
                else:
                    sumstats.iat[sum_num,sum_ind_col]=True
                    sum_num+=1
                           
    ## out
    
    sumstats = sumstats.loc[~sumstats["bed_indicator"],:]
    log.write(" -Number of variants in the specified regions to exclude:",sum(sumstats["bed_indicator"]), verbose=verbose)
    log.write(" -Number of variants left:",len(sumstats), verbose=verbose)
    sumstats = sumstats.drop(columns="bed_indicator")
    log.write("Finished filtering out variants.", verbose=verbose)
    gc.collect()
    return sumstats

def inferbuild(sumstats,status="STATUS",chrom="CHR", pos="POS", ea="EA", nea="NEA",build="19", verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "infer genome build version using hapmap3 SNPs"
    _end_line = "inferring genome build version using hapmap3 SNPs"
    _start_cols = [chrom,pos]
    _start_function = ".infer_build()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                              log=log,
                              verbose=verbose,
                              start_line=_start_line,
                              end_line=_end_line,
                              start_cols=_start_cols,
                              start_function=_start_function,
                              **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    inferred_build="Unknown"
    log.write("Start to infer genome build version using hapmap3 SNPs...", verbose=verbose)    
    data_path_19 =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db150_hg19.snplist.gz'    
    data_path_38 =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db151_hg38.snplist.gz'    
    log.write(" -Loading Hapmap3 variants data...", verbose=verbose)        
    hapmap3_ref_19 = pd.read_csv(data_path_19,sep="\s+",usecols=["#CHROM","POS"],dtype={"#CHROM":"string","POS":"string"})
    hapmap3_ref_38 = pd.read_csv(data_path_38,sep="\s+",usecols=["#CHROM","POS"],dtype={"#CHROM":"string","POS":"string"})
    
    log.write(" -CHR:POS will be used for matching...", verbose=verbose)
    raw_chrpos = sumstats[chrom].astype("string")+":"+sumstats[pos].astype("string")
    
    hapmap3_ref_19["chr:pos"] = hapmap3_ref_19["#CHROM"]+":"+hapmap3_ref_19["POS"]
    hapmap3_ref_38["chr:pos"] = hapmap3_ref_38["#CHROM"]+":"+hapmap3_ref_38["POS"]
    
    match_count_for_19 = sum(raw_chrpos.isin(hapmap3_ref_19["chr:pos"].values))
    match_count_for_38 = sum(raw_chrpos.isin(hapmap3_ref_38["chr:pos"].values))
    
    log.write(" -Matching variants for hg19: num_hg19 = ",match_count_for_19, verbose=verbose)        
    log.write(" -Matching variants for hg38: num_hg38 = ",match_count_for_38, verbose=verbose) 
    
    if max(match_count_for_19, match_count_for_38)<10000:
        log.warning("Please be cautious due to the limited number of variants.", verbose=verbose) 
    
    if match_count_for_19 > match_count_for_38:
        log.write(" -Since num_hg19 >> num_hg38, assigning genome build hg19...", verbose=verbose) 
        sumstats[status] = vchange_status(sumstats[status],1,"9","1")
        sumstats[status] = vchange_status(sumstats[status],2,"9","9")
        inferred_build="19"
    elif match_count_for_19 < match_count_for_38:
        log.write(" -Since num_hg19 << num_hg38, assigning genome build hg38...", verbose=verbose) 
        sumstats[status] = vchange_status(sumstats[status],1,"9","3")
        sumstats[status] = vchange_status(sumstats[status],2,"9","8")
        inferred_build="38"
    else:
        log.write(" -Since num_hg19 = num_hg38, unable to infer...", verbose=verbose) 
        
    finished(log,verbose,_end_line)
    return sumstats, inferred_build

def sampling(sumstats,n=1, p=None, verbose=True,log=Log(),**args):

    log.write("Start to randomly select variants from the sumstats...", verbose=verbose) 
    if p is None:
        log.write(" -Number of variants selected from the sumstats:",n, verbose=verbose)
        if n > len(sumstats):
            raise ValueError("Please input a number < {}".format(len(sumstats)))
    else:
        if p>-0.00000001 and p<1.00000001:
            log.write(" -Percentage of variants selected from the sumstats: ",p, verbose=verbose)
            n = int(len(sumstats)*p)
            log.write(" -Number of variants selected from the sumstats:",n, verbose=verbose)
        else:
            raise ValueError("Please input a number in (0,1)")
    
    if "random_state" in args.keys():
        log.write(" -Random state (seed): {}".format(args["random_state"]), verbose=verbose)
    else:
        args["random_state"] = np.random.randint(0,4294967295)
        log.write(" -Random state (seed): {}".format(args["random_state"]), verbose=verbose)
    sampled = sumstats.sample(n=n,**args)
    log.write("Finished sampling...", verbose=verbose)
    gc.collect()
    return sampled

def _get_flanking(sumstats, snpid, windowsizekb=500, verbose=True,log=Log(),**args):
    
    log.write("Start to extract variants in the flanking regions:",verbose=verbose)
    log.write(" - Central variant: {}".format(snpid))
    
    row = sumstats.loc[sumstats["SNPID"]==snpid,:]
    
    for index, row in row.iterrows():
        chrom = row["CHR"]
        left = row["POS"] - 1000 * windowsizekb
        right = row["POS"] + 1000 * windowsizekb
    
    log.write(" - Flanking regions: {}:{}-{}".format(chrom, left, right ))

    flanking = sumstats.query("CHR==@chrom & POS > @left & POS < @right ")
    
    log.write(" - Extracted {} variants in the regions.".format(len(flanking)),verbose=verbose)
    log.write("Finished extracting variants in the flanking regions.",verbose=verbose)

    return flanking

def _get_flanking_by_id(sumstats, snpid, windowsizekb=500, verbose=True,log=Log(),**args):
    
    log.write("Start to extract variants in the flanking regions using rsID or SNPID...",verbose=verbose)
    log.write(" - Central variants: {}".format(snpid), verbose=verbose)
    log.write(" - Flanking windowsize in kb: {}".format(windowsizekb), verbose=verbose)

    if type(snpid) == str:
        snpid = [snpid]
    
    if "rsID" in sumstats.columns and "SNPID" not in sumstats.columns:
        is_specified = sumstats["rsID"].isin(snpid)
    elif "rsID" not in sumstats.columns and "SNPID" in sumstats.columns:
        is_specified = sumstats["SNPID"].isin(snpid)
    else:
        is_specified = sumstats["rsID"].isin(snpid) | sumstats["SNPID"].isin(snpid)

    row = sumstats.loc[is_specified,:]
    
    is_flanking = None
    for index, row in row.iterrows():
        chrom = row["CHR"]
        left =  row["POS"] - 1000 * windowsizekb
        right = row["POS"] + 1000 * windowsizekb
        is_flancking_in_this_region = (sumstats["CHR"] == chrom) & (sumstats["POS"] >= left) & (sumstats["POS"] <= right)
        
        log.write(" - Variants in flanking region {}:{}-{} : {}".format(chrom, left, right, sum(is_flancking_in_this_region) ))
        
        if is_flanking is None:
            is_flanking = is_flancking_in_this_region
        else:
            is_flanking = is_flanking | is_flancking_in_this_region
    
    flanking = sumstats.loc[is_flanking,:]
    
    log.write(" - Extracted {} variants in the regions.".format(len(flanking)),verbose=verbose)
    log.write("Finished extracting variants in the flanking regions.",verbose=verbose)

    return flanking

def _get_flanking_by_chrpos(sumstats, chrpos, windowsizekb=500, verbose=True,log=Log(),**args):
    
    log.write("Start to extract variants in the flanking regions using CHR and POS...",verbose=verbose)
    log.write(" - Central positions: {}".format(chrpos), verbose=verbose)
    log.write(" - Flanking windowsize in kb: {}".format(windowsizekb), verbose=verbose)

    if type(chrpos) == tuple:
        chrpos_to_check = [chrpos]
    else:
        chrpos_to_check = chrpos

    is_flanking = None
    
    for index, row in enumerate(chrpos_to_check):
        chrom = row[0]
        left =  row[1] - 1000 * windowsizekb
        right = row[1] + 1000 * windowsizekb
        is_flancking_in_this_region = (sumstats["CHR"] == chrom) & (sumstats["POS"] >= left) & (sumstats["POS"] <= right)
        
        log.write(" - Variants in flanking region {}:{}-{} : {}".format(chrom, left, right, sum(is_flancking_in_this_region) ))
        
        if is_flanking is None:
            is_flanking = is_flancking_in_this_region
        else:
            is_flanking = is_flanking | is_flancking_in_this_region
    
    flanking = sumstats.loc[is_flanking,:]
    
    log.write(" - Extracted {} variants in the regions.".format(len(flanking)),verbose=verbose)
    log.write("Finished extracting variants in the flanking regions.",verbose=verbose)

    return flanking

def _filter_palindromic(sumstats, mode="in", ea="EA",nea="NEA", log=Log(),verbose=True):
    log.write("Start to filter palindromic variants...",verbose=verbose)
    is_palindromic_snp = is_palindromic(sumstats[[nea,ea]],a1=nea,a2=ea)   
    
    log.write(" -Identified palindromic variants: {}".format(sum(is_palindromic_snp)),verbose=verbose)
    
    if mode=="in":
        palindromic = sumstats.loc[is_palindromic_snp,:]
    else:
        palindromic = sumstats.loc[~is_palindromic_snp,:]

    log.write("Finished filtering palindromic variants.",verbose=verbose)
    return palindromic

def _filter_indel(sumstats, mode="in", ea="EA",nea="NEA", log=Log(),verbose=True):
    log.write("Start to filter indels...",verbose=verbose)
    is_indel = (sumstats[ea].str.len()!=sumstats[nea].str.len()) 
    
    log.write(" -Identified indels: {}".format(sum(is_indel)),verbose=verbose)
    if mode=="in":
        indel = sumstats.loc[is_indel,:]
    else:
        indel = sumstats.loc[~is_indel,:]
    log.write("Finished filtering indels.",verbose=verbose)
    return indel

def _filter_snp(sumstats, mode="in", ea="EA",nea="NEA", log=Log(),verbose=True):
    log.write("Start to filter SNPs...",verbose=verbose)
    is_snp = (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()==1)
    
    log.write(" -Identified SNPs: {}".format(sum(is_snp)),verbose=verbose)
    if mode=="in":
        snp = sumstats.loc[is_snp,:]
    else:
        snp = sumstats.loc[~is_snp,:]
    log.write("Finished filtering SNPs.",verbose=verbose)
    return snp

def _exclude_hla(sumstats, chrom="CHR", pos="POS", lower=25000000 ,upper=34000000 ,log=Log(), verbose=True):
    
    raw_len = len(sumstats)
    
    if str(sumstats[chrom].dtype) == "string":
        is_in_hla = (sumstats[chrom].astype("string")=="6")&(sumstats[pos]>lower)&(sumstats[pos]<upper)
    else:
        is_in_hla = (sumstats[chrom]==6)&(sumstats[pos]>lower)&(sumstats[pos]<upper)
    
    sumstats = sumstats.loc[~is_in_hla, : ]
    
    after_len = len(sumstats)
    
    log.write(" -Excluded {} variants in HLA region (chr6: {}-{} )...".format(raw_len - after_len,lower,upper),verbose=verbose)
    
    return sumstats

def _extract(sumstats, extract=None, id_use="SNPID", log=Log(), verbose=True ):
    if extract is not None:
        log.write(" -Extracting {} variants from sumstats...".format(len(extract)),verbose=verbose)
        sumstats = sumstats.loc[sumstats[id_use].isin(extract),:]
        log.write(" -Extracted {} variants from sumstats...".format(len(sumstats)),verbose=verbose)
    return sumstats

def _exclude(sumstats, exclude=None, id_use="SNPID", log=Log(), verbose=True ):
    if exclude is not None:
        log.write(" -Excluding {} variants from sumstats...".format(len(exclude)),verbose=verbose)
        sumstats = sumstats.loc[~sumstats[id_use].isin(exclude),:]
        log.write(" -Excluded {} variants from sumstats...".format(len(sumstats)),verbose=verbose)
    return sumstats