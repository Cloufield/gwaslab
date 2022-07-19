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