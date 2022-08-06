import re
import gwaslab as gl
import pandas as pd
import numpy as np

def filterout(sumstats,lt={},gt={},eq={},remove=False,verbose=True,log=gl.Log()):
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
    return sumstats.copy()

def filterin(sumstats,lt={},gt={},eq={},remove=False,verbose=True,log=gl.Log()):
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
    return sumstats.copy()

def filterregionin(sumstats,path=None, chrom="CHR",pos="POS", high_ld=False, build="19", verbose=True,log=gl.Log()):
    if verbose: log.write("Start filtering in variants if in intervals defined in bed files:")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))
    
    if high_ld is True:
        path = gl.get_high_ld(build=build)
        if verbose: log.write(" -Loading bed format file for hg"+build)

    else:
        if verbose: log.write(" -Loading bed format file: " , path)
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    dic=gl.get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr").map(dic).astype("Int64")
    
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
                           
    ## in
    
    sumstats = sumstats.loc[sumstats["bed_indicator"],:]
    if verbose: log.write(" -After filtering, number of variants in the specified regions:",len(sumstats))
    sumstats = sumstats.drop(columns="bed_indicator")
    return sumstats
    
def filterregionout(sumstats, path=None, chrom="CHR",pos="POS", high_ld=False, build="19", verbose=True,log=gl.Log()):
    if verbose: log.write("Start filtering out variants if in intervals defined in bed files:")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))
    if high_ld is True:
        path = gl.get_high_ld(build=build)
        if verbose: log.write(" -Loading bed format file for hg"+build)

    else:
        if verbose: log.write(" -Loading bed format file: " , path)
            
    bed = pd.read_csv(path,sep="\s+",header=None,dtype={0:"string",1:"Int64",2:"Int64"})
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    dic=gl.get_chr_to_number(out_chr=True)
    bed[0]=bed[0].str.strip("chr").map(dic).astype("Int64")
    
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
    if verbose: log.write(" -After filtering, number of variants in the specified regions:",len(sumstats))
    sumstats = sumstats.drop(columns="bed_indicator")
    return sumstats