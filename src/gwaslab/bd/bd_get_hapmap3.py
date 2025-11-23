import pandas as pd
from os import path
from gwaslab.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
from pathlib import Path

#A unique identifier (e.g., the rs number)
#Allele 1 (effect allele)
#Allele 2 (non-effect allele)
#Sample size (which often varies from SNP to SNP)
#A P-value
#A signed summary statistic (beta, OR, log odds, Z-score, etc)
@with_logging(
        start_to_msg="extract HapMap3 SNPs",
        finished_msg="extracting HapMap3 SNPs",
        start_function=".gethapmap3"
)
def gethapmap3(sumstats,rsid="rsID",chrom="CHR", pos="POS", ea="EA", nea="NEA",build="19", verbose=True, match_allele= True, how="inner", log=Log()):
    """
    Extract HapMap3 SNPs from summary statistics based on rsID or genomic coordinates.

    Parameters
    ----------

    build : str, optional
        Genome build version ("19" or "38"). Default is "19".
    verbose : bool, optional
        Print progress messages. Default is True.
    match_allele : bool, optional
        Check allele matching. Default is True.
    how : str, optional
        Type of merge to perform. Default is "inner".


    Returns
    -------
    pd.DataFrame
        Filtered summary statistics with HapMap3 SNPs.

    Less used parameters
    --------
    sumstats : pd.DataFrame
        Input summary statistics dataframe.
    rsid : str, optional
        Column name for rsID. Default is "rsID".
    chrom : str, optional
        Column name for chromosome. Default is "CHR".
    pos : str, optional
        Column name for position. Default is "POS".
    ea : str, optional
        Column name for effect allele. Default is "EA".
    nea : str, optional
        Column name for non-effect allele. Default is "NEA".
    log : Log, optional
        Logging object. Default is Log().
    """
    if build=="19":
        #data_path =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db150_hg19.snplist.gz'
        data_path = path.join( Path(__file__).parents[1], "data","hapmap3_SNPs","hapmap3_db150_hg19.snplist.gz")
    elif build=="38":
        #data_path =  path.dirname(__file__) + '/data/hapmap3_SNPs/hapmap3_db151_hg38.snplist.gz'
        data_path = path.join( Path(__file__).parents[1], "data","hapmap3_SNPs","hapmap3_db151_hg38.snplist.gz")

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
        hapmap3_ref = hapmap3_ref.rename(columns={"rsid":rsid})
        
        output = sumstats.loc[sumstats[rsid].isin(hapmap3_ref[rsid].values),:].copy()
        
        output = pd.merge(output, hapmap3_ref, on = rsid, how=how, suffixes=('', '_hapmap3'))

        raw_rsid_count= len(output)
        log.write(f" -Raw input contains {raw_rsid_count} Hapmap3 variants based on rsID...", verbose=verbose)

        if match_allele:
            log.write(" -Checking if alleles are same...")
            is_matched = ((output[ea].astype("string") == output["A1"]) & (output[nea].astype("string") == output["A2"])) \
                            | ((output[ea].astype("string") == output["A2"]) & (output[nea].astype("string") == output["A1"]))
            if how=="right":
                is_matched = ((output[ea].astype("string") == output["A1"]) & (output[nea].astype("string") == output["A2"])) \
                            | ((output[ea].astype("string") == output["A2"]) & (output[nea].astype("string") == output["A1"])) | output[ea].isna()
            output = output.loc[is_matched,:]
            output = output.drop(columns=["#CHROM","A1","A2"] )
            log.write(f" -Filtered {raw_rsid_count - len(output)} Hapmap3 variants due to unmatech alleles...", verbose=verbose)
        
        for i in ["#CHROM","A1","A2","POS_hapmap3"]:
            todrop=[]
            if i in output.columns:
                todrop.append(i)
        output = output.drop(columns=todrop)
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
        return output
    else:
        raise ValueError("Not enough information to match SNPs. Please check your sumstats...")
