import pandas as pd
from os import path
from gwaslab.Log import Log
#A unique identifier (e.g., the rs number)
#Allele 1 (effect allele)
#Allele 2 (non-effect allele)
#Sample size (which often varies from SNP to SNP)
#A P-value
#A signed summary statistic (beta, OR, log odds, Z-score, etc)

def gethapmap3(sumstats,rsid="rsID",chrom="CHR", pos="POS", ea="EA", nea="NEA",build="19", verbose=True,log=Log()):
    if verbose:log.write(" -Processing "+str(len(sumstats))+" raw variants...")

    if build=="19":
        data_path =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db150_hg19.snplist.gz'
    elif build=="38":
        data_path =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db151_hg38.snplist.gz'
    
    if verbose:log.write(" -Loading Hapmap3 variants data...")
        
    hapmap3_ref = pd.read_csv(data_path,sep="\s+",usecols=["#CHROM","POS","rsid"],dtype={"#CHROM":"string","POS":"string"})
    #rsid    A1      A2      #CHROM  POS
    #rs3094315       G       A       1       752566
    if rsid in sumstats.columns:
        output = sumstats.loc[sumstats[rsid].isin(hapmap3_ref["rsid"].values),:].copy()
        return output
    elif chrom in sumstats.columns and pos in sumstats.columns:
        if verbose: log.write(" -Since rsID not in sumstats, chr:pos( build "+build+") will be used for matching...")
        sumstats   ["chr:pos"] = sumstats[chrom].astype("string")+":"+sumstats[pos].astype("string")
        hapmap3_ref["chr:pos"] = hapmap3_ref["#CHROM"]+":"+hapmap3_ref["POS"]
        hapmap3_ref = hapmap3_ref.rename(columns={"rsid":"rsID"})
        output = pd.merge(sumstats,hapmap3_ref.loc[:,["chr:pos","rsID"]],left_on="chr:pos",right_on="chr:pos",how="inner",suffixes=('', '_hapmap3')).copy()
        output = output.drop(columns="chr:pos")
        if verbose: log.write(" -Raw input contains "+str(len(output))+" hapmaps variants based on chr:pos...")
        return output
    else:
        raise ValueError("Not enough information to match SNPs. Please check your sumstats...")