from typing import TYPE_CHECKING, Union, Optional
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
import gc

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

@with_logging(
        start_to_msg="calculate signal DENSITY",
        finished_msg="calculating signal DENSITY successfully!"
)
def _get_signal_density2(
    insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    snpid: str = "SNPID",
    chrom: str = "CHR",
    pos: str = "POS",
    bwindowsizekb: int = 100,
    sig_sumstats: Optional[pd.DataFrame] = None,
    log: Optional[Log] = None,
    verbose: bool = True
    ) -> pd.DataFrame:
    """
    Calculate signal density in genomic data using a sliding window approach.

    This function computes signal density by analyzing the distribution of variants
    across the genome within specified window sizes. It provides statistical summaries
    of density values including mean, median, standard deviation, and maximum values.

    Parameters
    ----------
    insumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame containing variants.
    snpid : str, optional
        Column name containing variant identifiers. Default is "SNPID".
    chrom : str, optional
        Column name containing chromosome numbers. Default is "CHR".
    pos : str, optional
        Column name containing genomic positions. Default is "POS".
    bwindowsizekb : int, optional
        Window size in kilobases for density calculation. Default is 100.
    sig_sumstats : pandas.DataFrame, optional
        Summary statistics DataFrame containing significant variants. If provided,
        density is calculated based on significant variants (conditional analysis).
        If None, density is calculated based on all variants. Default is None.
    log : Log, optional
        Log object for writing messages. Default is None.
    verbose : bool, optional
        Whether to display progress messages. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with added "DENSITY" column containing calculated density values.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data
    
    if log is None:
        class Dummy:
            def write(self, *args, **kwargs): pass
        log = Dummy()
    
    # Auto-detect ID column: prefer SNPID, fallback to rsID
    if snpid == "SNPID" and snpid not in insumstats.columns:
        if "rsID" in insumstats.columns:
            snpid = "rsID"
    
    wsize = bwindowsizekb * 1000

    sumstats = insumstats.copy()
    sumstats = sumstats.sort_values([chrom, pos], ignore_index=True)

    sumstats[pos] = pd.to_numeric(sumstats[pos], errors="coerce")

    valid_mask = sumstats[chrom].notna() & sumstats[pos].notna()

    densities = pd.Series(pd.NA, index=sumstats.index, dtype="Int32")

    # If sig_sumstats is provided, calculate density based on significant variants
    if sig_sumstats is not None:
        sig_sumstats = sig_sumstats.copy()
        sig_sumstats[pos] = pd.to_numeric(sig_sumstats[pos], errors="coerce")
        sig_valid_mask = sig_sumstats[chrom].notna() & sig_sumstats[pos].notna()
        
        # Initialize density to 0 for all variants
        densities = pd.Series(0, index=sumstats.index, dtype="Int32")
        counter = 0
        
        # For each significant variant, count how many variants in insumstats are within window
        for chrom_i, sig_chr in sig_sumstats.loc[sig_valid_mask, :].groupby(chrom, sort=False):
            sig_chr_valid = sig_chr.loc[sig_chr[pos].notna(), :]
            if sig_chr_valid.empty:
                continue
                
            # Get corresponding chromosome data from sumstats
            chr_mask = (sumstats[chrom] == chrom_i) & valid_mask
            if not chr_mask.any():
                continue
                
            chr_data = sumstats.loc[chr_mask, :]
            chr_positions = chr_data[pos].to_numpy()
            sig_positions = sig_chr_valid[pos].to_numpy()
            
            # For each significant variant, find all variants within window
            for sig_pos in sig_positions:
                counter += 1
                # Find variants within window: [sig_pos - wsize, sig_pos + wsize]
                left_idx = np.searchsorted(chr_positions, sig_pos - wsize, side="left")
                right_idx = np.searchsorted(chr_positions, sig_pos + wsize, side="right")
                # Increment density for variants in this window
                chr_indices = chr_data.index[left_idx:right_idx]
                densities.loc[chr_indices] += 1
                
                if counter % 1000 == 0:
                    log.write(f" -Processed {counter//1000}k signals", verbose=verbose)
                    gc.collect()
    else:
        # Original behavior: calculate density based on all variants
        for chrom_i, df_chr in sumstats.loc[valid_mask, :].groupby(chrom, sort=False):
            df_chr_valid = df_chr.loc[df_chr[pos].notna(), :]
            positions = df_chr_valid[pos].to_numpy()
            # Use searchsorted to find right edges efficiently
            # For each variant, find the rightmost index within +window
            right_idx = np.searchsorted(positions, positions + wsize, side="right")
            left_idx = np.searchsorted(positions, positions - wsize, side="left")
            # Count how many fall within window (excluding itself)
            density_chr = (right_idx - left_idx - 1).astype(np.int32)
            densities[df_chr_valid.index] = density_chr

    sumstats["DENSITY"] = densities

    # Basic stats
    density_valid = sumstats["DENSITY"].dropna()
    if density_valid.empty:
        bmean = pd.NA
        bmedian = pd.NA
        bsd = pd.NA
        bmax = pd.NA
        bmaxid = "NA"
    else:
        bmean = density_valid.mean()
        bmedian = density_valid.median()
        # Handle std() for single value case
        if len(density_valid) > 1:
            bsd = density_valid.std()
        else:
            bsd = 0.0
        bmax = density_valid.max()
        
        # Ensure snpid column exists for reporting
        if snpid not in sumstats.columns:
            if "rsID" in sumstats.columns:
                snpid = "rsID"
            else:
                bmaxid = "NA"
        else:
            bmaxid = sumstats.loc[density_valid.idxmax(), snpid]

    log.write(f" -Mean : {bmean:.3f} signals per {bwindowsizekb} kb", verbose=verbose)
    log.write(f" -SD : {bsd:.3f}", verbose=verbose)
    log.write(f" -Median : {bmedian:.3f} signals per {bwindowsizekb} kb", verbose=verbose)
    log.write(f" -Max : {bmax} signals per {bwindowsizekb} kb at variant {bmaxid}", verbose=verbose)

    return sumstats