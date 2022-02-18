import pandas as pd
from os import path
def ldsc(insumstats,chrom=None,pos=None,rsid=None,verbose=True):
    if verbose:print("  - Processing "+str(len(insumstats))+" raw variants...")
    if rsid:
        sumstats = insumstats[rsid]
    elif chrom and pos:
        sumstats = insumstats.loc[:,[chrom,pos]].astype("string")
        sumstats["chr:pos"]=sumstats[chrom]+":"+sumstats[pos] 
    data_path = path.join(path.dirname(__file__), 'data/hapmap3_rs_chr_pos_a1_a2/hapmap3_db150_hg19.snplist.gz')
    hapmap3_ref = pd.read_csv(data_path,sep="\s+",dtype="string")
    
    if rsid:
        output=insumstats.loc[sumstats.isin(hapmap3_ref["rsid"].values),:]
        if verbose:print("  - Raw input contains "+str(len(output))+" hapmaps variants based on rsid...")
        return output
    elif chrom and pos:
        hapmap3_ref["chr:pos"]=hapmap3_ref["#CHROM"]+":"+hapmap3_ref["POS"]
        output=insumstats.loc[sumstats["chr:pos"].isin(hapmap3_ref["chr:pos"].values),:]
        if verbose:print("  - Raw input contains "+str(len(output))+" hapmaps variants based on chr:pos...")
        return output