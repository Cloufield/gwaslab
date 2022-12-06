import re
#import modin.pandas as pd
import pandas as pd
import numpy as np
from os import path
from gwaslab.CommonData import get_high_ld
from gwaslab.CommonData import get_chr_to_number
from gwaslab.Log import Log
from gwaslab.vchangestatus import vchange_status
from gwaslab.fixdata import sortcoordinate

import gc
def filtervalues(sumstats,expr,remove=False,verbose=True,log=Log()):
    if verbose: log.write("Start filtering values by condition:",expr)
    prenum = len(sumstats)
    sumstats = sumstats.query(expr,engine='python').copy()
    afternum = len(sumstats)
    if verbose: log.write(" -Removing "+ str(prenum-afternum) +" variants not meeting the conditions:",expr)
    if verbose: log.write("Finished filtering values.")
    gc.collect()
    return sumstats

def filterout(sumstats,interval={},lt={},gt={},eq={},remove=False,verbose=True,log=Log()):
    if verbose: log.write("Start filtering values:")
    for key,threshold in gt.items():
        num = len(sumstats.loc[sumstats[key]>threshold,:])
        if verbose:log.write(" -Removing "+ str(num) +" variants with "+key+" > "+ str(threshold)+" ...")
        sumstats = sumstats.loc[sumstats[key]<threshold,:]
    for key,threshold in lt.items():
        num = len(sumstats.loc[sumstats[key]<threshold,:])
        if verbose:log.write(" -Removing "+ str(num) +" variants with "+key+" < "+ str(threshold)+" ...")
        sumstats = sumstats.loc[sumstats[key]>threshold,:]
    for key,threshold in eq.items():
        num = len(sumstats.loc[sumstats[key]==threshold,:])
        if verbose:log.write(" -Removing "+ str(num) +" variants with "+key+" = "+ str(threshold)+" ...")
        sumstats = sumstats.loc[sumstats[key]!=threshold,:]
    if verbose: log.write("Finished filtering values.")
    gc.collect()
    return sumstats.copy()

def filterin(sumstats,lt={},gt={},eq={},remove=False,verbose=True,log=Log()):
    if verbose: log.write("Start filtering values:")
    for key,threshold in gt.items():
        num = len(sumstats.loc[sumstats[key]>threshold,:])
        if verbose:log.write(" -Keeping "+ str(num) +" variants with "+key+" > "+ str(threshold)+" ...")
        sumstats = sumstats.loc[sumstats[key]>threshold,:]
    for key,threshold in lt.items():
        num = len(sumstats.loc[sumstats[key]<threshold,:])
        if verbose:log.write(" -Keeping "+ str(num) +" variants with "+key+" < "+ str(threshold)+" ...")
        sumstats = sumstats.loc[sumstats[key]<threshold,:]
    for key,threshold in eq.items():
        num = len(sumstats.loc[sumstats[key]==threshold,:])
        if verbose:log.write(" -Keeping "+ str(num) +" variants with "+key+" = "+ str(threshold)+" ...")
        sumstats = sumstats.loc[sumstats[key]==threshold,:]
    if verbose: log.write("Finished filtering values.")
    gc.collect()
    return sumstats.copy()

def filterregionin(sumstats,path=None, chrom="CHR",pos="POS", high_ld=False, build="19", verbose=True,log=Log()):
    sumstats = sortcoordinate(sumstats,verbose=verbose)
    if verbose: log.write("Start to filter in variants if in intervals defined in bed files:")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))
    
    if high_ld is True:
        path = get_high_ld(build=build)
        if verbose: log.write(" -Loading bed format file for hg"+build)

    else:
        if verbose: log.write(" -Loading bed format file: " , path)
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    dic=get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr")
    bed[0]=bed[0].map(dic).astype("string")
    bed[0]=bed[0].astype("Int64")
    sumstats = sumstats.sort_values(["CHR","POS"])
    
    if len(bed)<100:
        if verbose: log.write(" -Bed file < 100 lines: using pd IntervalIndex... ")
        for i in sumstats[chrom].unique():
            if sum(bed[0]==i)>0:
                interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==i,"tuple"])
                sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
            else:
                continue
    else:
        if verbose: log.write(" -Bed file > 100 lines: using two pointers, please make files are all sorted... ")
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
    if verbose: log.write(" -Number of variants in the specified regions to keep:",sum(sumstats["bed_indicator"]))
    if verbose: log.write(" -Number of variants removed:",sum(~sumstats["bed_indicator"]))
    sumstats = sumstats.drop(columns="bed_indicator")
    if verbose: log.write("Finished filtering in variants.")
    gc.collect()
    return sumstats

def filterregionout(sumstats, path=None, chrom="CHR",pos="POS", high_ld=False, build="19", verbose=True,log=Log()):
    sumstats = sortcoordinate(sumstats,verbose=verbose)
    if verbose: log.write("Start to filter out variants if in intervals defined in bed files:")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))
    if high_ld is True:
        path = get_high_ld(build=build)
        if verbose: log.write(" -Loading bed format file for hg"+build)

    else:
        if verbose: log.write(" -Loading bed format file: " , path)
            
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    
    dic=get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr")
    bed[0]=bed[0].map(dic).astype("string")
    bed[0]=bed[0].astype("Int64")
    
    if len(bed)<100:
        if verbose: log.write(" -Bed file < 100 lines: using pd IntervalIndex... ")
        for i in sumstats[chrom].unique():
            if sum(bed[0]==i)>0:
                interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==i,"tuple"])
                sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
            else:
                continue
    else:
        if verbose: log.write(" -Bed file > 100 lines: using two pointers, please make files are all sorted... ")
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
    if verbose: log.write(" -Number of variants in the specified regions to exclude:",sum(sumstats["bed_indicator"]))
    if verbose: log.write(" -Number of variants left:",len(sumstats))
    sumstats = sumstats.drop(columns="bed_indicator")
    if verbose: log.write("Finished filtering out variants.")
    gc.collect()
    return sumstats

def inferbuild(sumstats,status="STATUS",chrom="CHR", pos="POS", ea="EA", nea="NEA",build="19", verbose=True,log=Log()):
    inferred_build="Unknown"
    if verbose:log.write("Start to infer genome build version using hapmap3 SNPs...")    
    data_path_19 =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db150_hg19.snplist.gz'    
    data_path_38 =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db151_hg38.snplist.gz'    
    if verbose:log.write(" -Loading Hapmap3 variants data...")        
    hapmap3_ref_19 = pd.read_csv(data_path_19,sep="\s+",usecols=["#CHROM","POS"],dtype={"#CHROM":"string","POS":"string"})
    hapmap3_ref_38 = pd.read_csv(data_path_38,sep="\s+",usecols=["#CHROM","POS"],dtype={"#CHROM":"string","POS":"string"})
    
    if chrom in sumstats.columns and pos in sumstats.columns:
        if verbose: log.write(" -CHR:POS will be used for matching...")
        raw_chrpos = sumstats[chrom].astype("string")+":"+sumstats[pos].astype("string")
        
        hapmap3_ref_19["chr:pos"] = hapmap3_ref_19["#CHROM"]+":"+hapmap3_ref_19["POS"]
        hapmap3_ref_38["chr:pos"] = hapmap3_ref_38["#CHROM"]+":"+hapmap3_ref_38["POS"]
        
        match_count_for_19 = sum(raw_chrpos.isin(hapmap3_ref_19["chr:pos"].values))
        match_count_for_38 = sum(raw_chrpos.isin(hapmap3_ref_38["chr:pos"].values))
        
        if verbose:log.write(" -Matching variants for hg19: num_hg19 = ",match_count_for_19)        
        if verbose:log.write(" -Matching variants for hg38: num_hg38 = ",match_count_for_38) 
        
        if max(match_count_for_19, match_count_for_38)<10000:
            if verbose:log.write(" -Warning: please be cautious due to the limited number of variants.") 
        
        if match_count_for_19 > match_count_for_38:
            if verbose:log.write(" -Since num_hg19 >> num_hg38, assigning genome build hg19...") 
            sumstats.loc[:,status] = vchange_status(sumstats.loc[:,status],1,"9","1")
            sumstats.loc[:,status] = vchange_status(sumstats.loc[:,status],2,"9","9")
            inferred_build="19"
        elif match_count_for_19 < match_count_for_38:
            if verbose:log.write(" -Since num_hg19 << num_hg38, assigning genome build hg38...") 
            sumstats.loc[:,status] = vchange_status(sumstats.loc[:,status],1,"9","3")
            sumstats.loc[:,status] = vchange_status(sumstats.loc[:,status],2,"9","8")
            inferred_build="38"
        else:
            if verbose:log.write(" -Since num_hg19 = num_hg38, unable to infer...") 
        gc.collect()
        return sumstats, inferred_build
    else:
        gc.collect()
        raise ValueError("Not enough information to match SNPs. Please check if CHR and POS columns are in your sumstats...")

def sampling(sumstats,n,verbose=True,log=Log()):
    if verbose:log.write("Start to randomly select variants from the sumstats...") 
    if verbose:log.write(" -Number of variants selected from the sumstats:",n)
    sampled = sumstats.sample(n=n)
    if verbose:log.write("Finished sampling...")
    gc.collect()
    return sampled

