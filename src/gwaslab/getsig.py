import pandas as pd
import numpy as np
import scipy as sp
import gwaslab as gl

def getsig(insumstats,
           id,
           chrom,
           pos,
           p,
           windowsizekb=500,
           sig_level=5e-8,
           log=gl.Log(),
           verbose=True):
    
    if verbose: log.write("Start extracting lead variants...")
    if verbose: log.write(" -Processing "+str(len(insumstats))+" variants...")
    if verbose: log.write(" -Significance threshold :", sig_level)
    if verbose: log.write(" -Sliding window size:", str(windowsizekb) ," kb")
    #load data
    
    sumstats=insumstats.loc[~insumstats[id].isna(),:].copy()
    
    #convert chrom to int
    
    sumstats[chrom]=sumstats[chrom].astype("string")
    sumstats.loc[sumstats[chrom]=="X",chrom] = "23"
    sumstats.loc[sumstats[chrom]=="Y",chrom] = "24"
    sumstats.loc[sumstats[chrom]=="MT",chrom] = "25"
    
    sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    sumstats[p] = pd.to_numeric(sumstats[p], errors='coerce')
    
    #extract all significant variants
    sumstats_sig = sumstats.loc[sumstats[p]<sig_level,:]
    if verbose:log.write(" -Found "+str(len(sumstats_sig))+" significant variants in total...")
    
    #sort the coordinates
    sumstats_sig = sumstats_sig.sort_values([chrom,pos])
    
    if sumstats_sig is None:
        if verbose:log.write(" -No lead snps at given significance threshold!")
        return None
    
    sig_index_list=[]
    current_sig_index = False
    current_sig_p = 1
    current_sig_pos = 0
    current_sig_chr = 0
    
    #iterate through all significant snps
    for line_number,(index, row) in enumerate(sumstats_sig.iterrows()):
        #when finished one chr 
        if row[chrom]!=current_sig_chr:
            #add the current lead variants id to lead variant list
            if current_sig_index is not False:sig_index_list.append(current_sig_index)
            
            #update lead vairant info to the new variant
            current_sig_chr=row[chrom]
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
            
            # only one significant variant on a new chromsome and this is the last sig variant
            if  line_number == len(sumstats_sig)-1:
                sig_index_list.append(current_sig_index)
            continue
        
        # next loci : gap > windowsizekb*1000
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
            
        #when last line in sig_index_list
        if  line_number == len(sumstats_sig)-1:
            sig_index_list.append(current_sig_index)
            continue
    
    if verbose:log.write(" -Identified "+str(len(sig_index_list))+" lead variants!")
    
    chromosome_conversion_dict = {str(i):str(i) for i in range(1,23)}
    chromosome_conversion_dict["23"] = "X"
    chromosome_conversion_dict["24"] = "Y"
    chromosome_conversion_dict["25"] = "MT"
    sumstats_sig.loc[:,chrom] = sumstats_sig[chrom].astype("string")
    sumstats_sig.loc[:,chrom] = sumstats_sig.loc[:,chrom].apply(lambda x:chromosome_conversion_dict[x])
    output = sumstats_sig.loc[sumstats_sig[id].isin(sig_index_list),:].copy()
    return output