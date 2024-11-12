import pandas as pd
from os import path
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import skipped
from gwaslab.qc_fix_sumstats import finished

#A unique identifier (e.g., the rs number)
#Allele 1 (effect allele)
#Allele 2 (non-effect allele)
#Sample size (which often varies from SNP to SNP)
#A P-value
#A signed summary statistic (beta, OR, log odds, Z-score, etc)

def gethapmap3(sumstats,rsid="rsID",chrom="CHR", pos="POS", ea="EA", nea="NEA",build="19", verbose=True, match_allele= True, how="inner", log=Log()):
    ##start function with col checking##########################################################
    _start_line = "extract HapMap3 SNPs"
    _end_line = "extracting HapMap3 SNPs"
    _start_cols =[]
    _start_function = ".gethapmap3"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None

    ############################################################################################
    if build=="19":
        data_path =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db150_hg19.snplist.gz'
    elif build=="38":
        data_path =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db151_hg38.snplist.gz'
    
    log.write(" -Loading Hapmap3 variants from built-in datasets...", verbose=verbose)
    
    if match_allele:
        additional_cols= ["A1","A2"]
    else:
        additional_cols=[]
    hapmap3_ref = pd.read_csv(data_path,sep="\s+",usecols=["#CHROM","POS","rsid"]+additional_cols, dtype={"#CHROM":"string","POS":"string"})

    #rsid    A1      A2      #CHROM  POS
    #rs3094315       G       A       1       752566
    
    if rsid in sumstats.columns:
        log.write(" -rsID will be used for matching...", verbose=verbose)
        output = sumstats.loc[sumstats[rsid].isin(hapmap3_ref["rsid"].values),:].copy()
        log.write(" -Raw input contains "+str(len(output))+" Hapmap3 variants based on rsID...", verbose=verbose)
        return output
    
    elif chrom in sumstats.columns and pos in sumstats.columns:
        log.write(" -Since rsID not in sumstats, CHR:POS( build "+build+") will be used for matching...", verbose=verbose)
        sumstats   ["chr:pos"] = sumstats[chrom].astype("string")+":"+sumstats[pos].astype("string")
        hapmap3_ref["chr:pos"] = hapmap3_ref["#CHROM"]+":"+hapmap3_ref["POS"]
        hapmap3_ref = hapmap3_ref.rename(columns={"rsid":"rsID"})
        output = pd.merge(sumstats,hapmap3_ref.loc[:,["chr:pos","rsID"]+additional_cols],left_on="chr:pos",right_on="chr:pos",how=how,suffixes=('', '_hapmap3')).copy()
        if match_allele:
            log.write(" -Checking if alleles are same...")
            is_matched = ((output[ea].astype("string") == output["A1"]) & (output[nea].astype("string") == output["A2"])) \
                            | ((output[ea].astype("string") == output["A2"]) & (output[nea].astype("string") == output["A1"]))
            if how=="right":
                is_matched = ((output[ea].astype("string") == output["A1"]) & (output[nea].astype("string") == output["A2"])) \
                            | ((output[ea].astype("string") == output["A2"]) & (output[nea].astype("string") == output["A1"])) | output[ea].isna()

            log.write(" -Variants with macthed alleles: {}".format(sum(is_matched)))
            output = output.loc[is_matched,:]
        output = output.drop(columns=["chr:pos"]+additional_cols)
        log.write(" -Raw input contains "+str(len(output))+" Hapmap3 variants based on CHR:POS...", verbose=verbose)
        finished(log=log,verbose=verbose,end_line=_end_line)
        return output
    else:
        raise ValueError("Not enough information to match SNPs. Please check your sumstats...")