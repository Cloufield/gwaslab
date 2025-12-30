from typing import TYPE_CHECKING, Union
import pandas as pd
import polars as pl
from os import path
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
from pathlib import Path
from gwaslab.util.util_in_filter_value import _get_hapmap_full_polars

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

#A unique identifier (e.g., the rs number)
#Allele 1 (effect allele)
#Allele 2 (non-effect allele)
#Sample size (which often varies from SNP to SNP)
#A P-value
#A signed summary statistic (beta, OR, log odds, Z-score, etc)
@with_logging(
        start_to_msg="extract HapMap3 SNPs",
        finished_msg="extracting HapMap3 SNPs",
        start_function=".gethapmap3",
        required_species="homo sapiens"
)
def _get_hapmap3(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], rsid: str = "rsID", chrom: str = "CHR", pos: str = "POS", ea: str = "EA", nea: str = "NEA", build: str = "19", verbose: bool = True, match_allele: bool = True, how: str = "inner", log: Log = Log()) -> pd.DataFrame:
    """
    Extract HapMap3 SNPs from summary statistics based on rsID or genomic coordinates.

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
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

    build : str, optional
        Genome build version ("19" or "38"). Default is "19".
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
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    log.write(" -Loading Hapmap3 variants from built-in datasets...", verbose=verbose)
    
    # Use cached Polars function for fast loading
    hapmap3_ref_pl = _get_hapmap_full_polars(build, include_alleles=match_allele)
    
    # Convert to pandas for compatibility with existing merge logic
    # Only convert what we need - this is still faster than reading from disk
    hapmap3_ref = hapmap3_ref_pl.to_pandas()
    
    # Ensure CHR and POS are strings for compatibility with existing string concatenation logic
    hapmap3_ref["#CHROM"] = hapmap3_ref["CHR"].astype("string")
    hapmap3_ref["POS"] = hapmap3_ref["POS"].astype("string")
    
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
            log.write(f" -Filtered {raw_rsid_count - len(output)} Hapmap3 variants due to unmatched alleles...", verbose=verbose)
        
        for i in ["#CHROM","A1","A2","POS_hapmap3"]:
            todrop=[]
            if i in output.columns:
                todrop.append(i)
        output = output.drop(columns=todrop)
        return output
    
    elif chrom in sumstats.columns and pos in sumstats.columns:
        log.write(" -Since rsID not in sumstats, CHR:POS( build "+build+") will be used for matching...", verbose=verbose)
        
        # Use Polars for fast CHR:POS matching - much faster than string concatenation
        # Convert full sumstats to Polars to preserve all columns
        sumstats_pl = pl.from_pandas(sumstats)
        
        # Filter out null values and ensure correct types for join columns
        sumstats_pl = sumstats_pl.filter(
            pl.col(chrom).is_not_null() & pl.col(pos).is_not_null()
        ).with_columns([
            pl.col(chrom).cast(pl.Int64).alias("_CHR_join"),
            pl.col(pos).cast(pl.Int64).alias("_POS_join")
        ])
        
        # Prepare hapmap3 reference for join
        hapmap3_join = hapmap3_ref_pl.select([
            pl.col("CHR").alias("_CHR_join"),
            pl.col("POS").alias("_POS_join"),
            pl.col("rsid").alias("rsID")
        ] + (["A1", "A2"] if match_allele else []))
        
        # Map how parameter to Polars join type
        join_type_map = {"inner": "inner", "left": "left", "right": "right"}
        join_type = join_type_map.get(how, "inner")
        
        # Join on CHR and POS using Polars (fast)
        matched_pl = sumstats_pl.join(
            hapmap3_join,
            on=["_CHR_join", "_POS_join"],
            how=join_type,
            suffix="_hapmap3"
        )
        
        # Drop temporary join columns
        matched_pl = matched_pl.drop(["_CHR_join", "_POS_join"])
        
        # Convert back to pandas for allele matching and final processing
        output = matched_pl.to_pandas()
        
        if match_allele:
            log.write(" -Checking if alleles are same...")
            is_matched = ((output[ea].astype("string") == output["A1"]) & (output[nea].astype("string") == output["A2"])) \
                            | ((output[ea].astype("string") == output["A2"]) & (output[nea].astype("string") == output["A1"]))
            if how=="right":
                is_matched = ((output[ea].astype("string") == output["A1"]) & (output[nea].astype("string") == output["A2"])) \
                            | ((output[ea].astype("string") == output["A2"]) & (output[nea].astype("string") == output["A1"])) | output[ea].isna()

            log.write(" -Variants with matched alleles: {}".format(sum(is_matched)))
            output = output.loc[is_matched,:]
            output = output.drop(columns=["A1", "A2"], errors="ignore")
        
        log.write(" -Raw input contains "+str(len(output))+" Hapmap3 variants based on CHR:POS...", verbose=verbose)
        return output
    else:
        raise ValueError("Not enough information to match SNPs. Please check your sumstats...")
