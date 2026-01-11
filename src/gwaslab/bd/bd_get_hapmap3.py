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
    
    # Helper function for allele matching
    def _match_alleles(output_df, ea_col, nea_col, how_join):
        """Check if alleles match between sumstats and hapmap3."""
        ea_str = output_df[ea_col].astype("string")
        nea_str = output_df[nea_col].astype("string")
        a1_str = output_df["A1"].astype("string")
        a2_str = output_df["A2"].astype("string")
        
        is_matched = ((ea_str == a1_str) & (nea_str == a2_str)) | ((ea_str == a2_str) & (nea_str == a1_str))
        if how_join == "right":
            is_matched = is_matched | output_df[ea_col].isna()
        return is_matched
    
    # Helper function to drop temporary columns
    def _drop_temp_columns(output_df, cols_to_drop):
        """Drop temporary columns if they exist."""
        cols_to_remove = [col for col in cols_to_drop if col in output_df.columns]
        if cols_to_remove:
            output_df = output_df.drop(columns=cols_to_remove, errors="ignore")
        return output_df
    
    if rsid in sumstats.columns:
        log.write(" -rsID will be used for matching...", verbose=verbose)
        
        # Convert to pandas only for rsID path (needed for string operations)
        hapmap3_ref = hapmap3_ref_pl.to_pandas()
        hapmap3_ref["#CHROM"] = hapmap3_ref["CHR"].astype("string")
        hapmap3_ref["POS"] = hapmap3_ref["POS"].astype("string")
        hapmap3_ref = hapmap3_ref.rename(columns={"rsid": rsid})
        
        output = sumstats.loc[sumstats[rsid].isin(hapmap3_ref[rsid].values), :].copy()
        output = pd.merge(output, hapmap3_ref, on=rsid, how=how, suffixes=('', '_hapmap3'))

        raw_rsid_count = len(output)
        log.write(f" -Raw input contains {raw_rsid_count} Hapmap3 variants based on rsID...", verbose=verbose)

        if match_allele:
            log.write(" -Checking if alleles are same...")
            is_matched = _match_alleles(output, ea, nea, how)
            output = output.loc[is_matched, :]
            log.write(f" -Filtered {raw_rsid_count - len(output)} Hapmap3 variants due to unmatched alleles...", verbose=verbose)
        
        # Drop temporary columns
        output = _drop_temp_columns(output, ["#CHROM", "A1", "A2", "POS_hapmap3"])
        return output
    
    elif chrom in sumstats.columns and pos in sumstats.columns:
        log.write(" -Since rsID not in sumstats, CHR:POS( build "+build+") will be used for matching...", verbose=verbose)
        
        # Optimize: First filter by POS (more selective) for inner join, then join on CHR+POS
        # Get all POS values from hapmap3 for fast filtering
        hapmap3_pos_set = set(hapmap3_ref_pl["POS"].to_list())
        
        # Filter out null values in pandas (no Polars conversion overhead)
        sumstats_filtered = sumstats.dropna(subset=[chrom, pos])
        
        # Fast POS filtering: filter sumstats to POS values in hapmap3 (reduces dataset size significantly before join)
        # This is safe for all join types since non-matching POS values won't match anyway
        sumstats_filtered = sumstats_filtered[sumstats_filtered[pos].isin(hapmap3_pos_set)]
        
        # Prepare hapmap3 reference for join (convert to pandas only when needed)
        hapmap3_ref_pd = hapmap3_ref_pl.select([
            pl.col("CHR"),
            pl.col("POS"),
            pl.col("rsid").alias("rsID")
        ] + (["A1", "A2"] if match_allele else [])).to_pandas()
        
        # Join on CHR and POS using pandas merge (no Polars conversion overhead)
        output = pd.merge(
            sumstats_filtered,
            hapmap3_ref_pd,
            left_on=[chrom, pos],
            right_on=["CHR", "POS"],
            how=how,
            suffixes=("", "_hapmap3")
        )
        
        # Drop duplicate CHR and POS columns from hapmap3 (keep original sumstats columns)
        # The merge adds "CHR" and "POS" from hapmap3, but we keep the original columns from sumstats
        cols_to_drop = []
        # Only drop hapmap3's CHR/POS if they're different from the original column names
        if "CHR" in output.columns and chrom != "CHR":
            cols_to_drop.append("CHR")
        if "POS" in output.columns and pos != "POS":
            cols_to_drop.append("POS")
        if cols_to_drop:
            output = output.drop(columns=cols_to_drop, errors="ignore")
        
        if match_allele:
            log.write(" -Checking if alleles are same...")
            is_matched = _match_alleles(output, ea, nea, how)
            log.write(" -Variants with matched alleles: {}".format(sum(is_matched)))
            output = output.loc[is_matched, :]
            output = _drop_temp_columns(output, ["A1", "A2"])
        
        log.write(" -Raw input contains "+str(len(output))+" Hapmap3 variants based on CHR:POS...", verbose=verbose)
        return output
    else:
        raise ValueError("Not enough information to match SNPs. Please check your sumstats...")
