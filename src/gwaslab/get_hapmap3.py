import pandas as pd
from os import path

#A unique identifier (e.g., the rs number)
#Allele 1 (effect allele)
#Allele 2 (non-effect allele)
#Sample size (which often varies from SNP to SNP)
#A P-value
#A signed summary statistic (beta, OR, log odds, Z-score, etc)

def ldsc(insumstats,chrom=None,pos=None,rsid=None,signedstats=True,verbose=True):

    if verbose:print("  - Processing "+str(len(insumstats))+" raw variants...")
    
    if signedstats:
        if verbose:print("  - Signed stats to use: BETA..")
        if signedstats is True:
            signedstats = ["BETA"]
    
    if rsid:
        sumstats = insumstats[rsid]
    
    elif chrom and pos:
        insumstats["chr:pos"]=insumstats[chrom].astype("string")+":"+insumstats[pos].astype("string")
    
    data_path = '/home/heyunye/work/gwaslab/gwaslab/src/gwaslab/data/hapmap3_rs_chr_pos_a1_a2/hapmap3_db150_hg19.snplist.gz'
    
    hapmap3_ref = pd.read_csv(data_path,sep="\s+",dtype="string")
    
    if rsid:
        output = insumstats.loc[sumstats.isin(hapmap3_ref["rsid"].values),:]
        if verbose: print("  - Raw input contains "+str(len(output))+" hapmaps variants based on rsid...")
        return output.loc[:, ["rsID","EA","NEA","N","P"]+signedstats ]
    
    elif chrom and pos:
        hapmap3_ref["chr:pos"]=hapmap3_ref["#CHROM"]+":"+hapmap3_ref["POS"]
        
        output = pd.merge(insumstats,hapmap3_ref,left_on="chr:pos",right_on="chr:pos",how="inner",suffixes=('', '_hapmap3'))
        
        if verbose: print("  - Raw input contains "+str(len(output))+" hapmaps variants based on chr:pos...")
        
        output=output.rename(columns={"rsid":"rsID"})
        output["POS"] = output["POS"].astype('int')
        return output.loc[:, ["rsID","CHR","POS","EA","NEA","N","P"]+signedstats ]