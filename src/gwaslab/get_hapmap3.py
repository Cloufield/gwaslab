import pandas as pd
from os import path

#A unique identifier (e.g., the rs number)
#Allele 1 (effect allele)
#Allele 2 (non-effect allele)
#Sample size (which often varies from SNP to SNP)
#A P-value
#A signed summary statistic (beta, OR, log odds, Z-score, etc)

def forldsc(sumstats,chrom=None,pos=None,rsid=None,signedstats=["BETA"],verbose=True):
    if verbose: print("Start fornmatting to ldsc input format...")
    ldsc_cols = ["EA","NEA","N","P","EAF"] + signedstats +["SE"]
    usecols =[]
    for i in sumstats.columns:
        if i in ldsc_cols:
            usecols.append(i)
    if rsid: 
        usecols.insert(0,"rsID")
    else:
        usecols.insert(0,"POS")
        usecols.insert(0,"CHR")
        usecols.insert(0,"MARKERNAME")
            
    insumstats = sumstats.loc[:,usecols].copy()
    
    if verbose:print("  - Processing "+str(len(insumstats))+" raw variants...")
    if verbose:print("  - Signed stats to use: ", signedstats )
    
    data_path =  path.dirname(__file__) + '/data/hapmap3_rs_chr_pos_a1_a2/hapmap3_db150_hg19.snplist.gz'
    
    hapmap3_ref = pd.read_csv(data_path,sep="\s+",dtype={"#CHROM":"string","POS":"int"})
    #rsid    A1      A2      #CHROM  POS
    #rs3094315       G       A       1       752566
    
    if rsid:
        
        output = pd.merge(insumstats,hapmap3_ref,left_on="rsID",right_on="rsid",how="inner",suffixes=('', '_hapmap3'))
        #output = insumstats.loc[sumstats.isin(hapmap3_ref["rsid"].values),:]
        if verbose: print("  - Raw input contains "+str(len(output))+" hapmaps variants based on rsid...")
        output = output.rename(columns={"#CHROM":"CHR"})
        return output.loc[:, ["rsID","CHR","POS","EA","NEA","N","P"]+signedstats ]
    
    elif chrom and pos:
        insumstats["chr:pos"]=insumstats[chrom].astype("string")+":"+insumstats[pos].astype("string")
        
        hapmap3_ref["chr:pos"]=hapmap3_ref["#CHROM"]+":"+hapmap3_ref["POS"].astype("string")
        
        output = pd.merge(insumstats,hapmap3_ref,left_on="chr:pos",right_on="chr:pos",how="inner",suffixes=('', '_hapmap3'))
        
        if verbose: print("  - Raw input contains "+str(len(output))+" hapmaps variants based on chr:pos...")
        
        output = output.rename(columns={"rsid":"rsID"})
        output["POS"] = output["POS"].astype('int')
        return output.loc[:, ["rsID","CHR","POS","EA","NEA","N","P"]+signedstats ]