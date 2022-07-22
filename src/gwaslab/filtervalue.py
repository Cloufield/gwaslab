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
    return sumstats

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
    return sumstats

def filterinbed(sumstats,path, chrom="CHR",pos="POS", verbose=True,log=gl.Log()):
    if verbose: log.write("Start filtering in variants if in intervals defined in bed files:")
    bed = pd.read_csv(path,"\t",header=None)
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    for i in sumstats[chrom].unique():
        if sum(bed[0]==("chr"+str(i)))>0:
            interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==("chr"+str(i)),"tuple"])
            sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
        else:
            continue
    ## in
    sumstats = sumstats.loc[sumstats["bed_indicator"] is True,:]
    sumstats = sumstats.drop(columns="bed_indicator")
    return sumstats
    
def filteroutbed(sumstats, bed, chrom="CHR",pos="POS", verbose=True,log=gl.Log()):
    if verbose: log.write("Start filtering out variants if in intervals defined in bed files:")
    bed = pd.read_csv(path,"\t",header=None)
    bed["tuple"] = bed.apply(lambda x: (x[1],x[2]),axis=1)
    sumstats["bed_indicator"]=False
    for i in sumstats[chrom].unique():
        if sum(bed[0]==("chr"+str(i)))>0:
            interval = pd.IntervalIndex.from_tuples(bed.loc[bed[0]==("chr"+str(i)),"tuple"])
            sumstats.loc[sumstats[chrom]==i,"bed_indicator"] = sumstats.loc[sumstats[chrom]==i,pos].apply(lambda x: any(interval.contains(x)))
        else:
            continue
    ## out
    sumstats = sumstats.loc[sumstats["bed_indicator"] is False,:] 
    sumstats = sumstats.drop(columns="bed_indicator")
    return sumstats