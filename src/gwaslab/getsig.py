import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
def getsig(insumstats,id,chrom,pos,p,windowsizekb=500,verbose=True,sig_level=5e-8):
    
    if verbose: print("  - Processing "+str(len(insumstats))+" variants...")
    sumstats=insumstats.loc[~insumstats[id].isna(),:]
    sumstats_sig = sumstats.loc[sumstats[p]<sig_level,:]
    if verbose:print("  - Found "+str(len(sumstats_sig))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")
    sumstats_sig = sumstats_sig.sort_values([chrom,pos])
    
    if len(sumstats_sig)==0:
        if verbose:print("  - No lead snps at given significance threshold!")
        return None
    
    sig_index_list=[]
    current_sig_index = False
    current_sig_p = 1
    current_sig_pos = 0
    current_sig_chr = 0
    for line_number,(index, row) in enumerate(sumstats_sig.iterrows()):
        #when finished one chr 
        if row[chrom]>current_sig_chr:
            if current_sig_index:sig_index_list.append(current_sig_index)
            current_sig_chr=row[chrom]
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
            continue
        
        #next loci
        if row[pos]>current_sig_pos + windowsizekb*1000:
            sig_index_list.append(current_sig_index)
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
            continue
        # update current pos and p
        if row[p]<current_sig_p:
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
        else:
            current_sig_pos=row[pos]
        #when last line
        if  line_number == len(sumstats_sig)-1:
            sig_index_list.append(current_sig_index)
            continue
    if verbose:print("  - Identified "+str(len(sig_index_list))+" lead variants!")
    return sumstats_sig.loc[sumstats_sig[id].isin(sig_index_list),:]