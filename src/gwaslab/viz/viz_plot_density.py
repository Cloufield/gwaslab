from gwaslab.viz.viz_aux_quickfix import _get_largenumber
import pandas as pd
from gwaslab.util.util_in_get_density import _get_signal_density2

def _process_density(sumstats, mode, bwindowsizekb, chrom, pos, verbose, log):
    """Create Brisbane plot (variant density plot)
    Y axis: number of neiboring variants within a window around each variant)
    X axis: same as Manhattan plot
"""

    if "b" in mode and "DENSITY" not in sumstats.columns:
        
        sumstats =  _get_signal_density2(insumstats_or_dataframe=sumstats,
                                                snpid="SNPID",
                                                chrom=chrom,
                                                pos=pos,
                                                bwindowsizekb=bwindowsizekb,
                                                log=log,
                                                verbose=verbose)

        bmean=sumstats.drop_duplicates(subset="SNPID")["DENSITY"].mean()
        bmedian=sumstats.drop_duplicates(subset="SNPID")["DENSITY"].median()
    elif "b" in mode and "DENSITY" in sumstats.columns:
        bmean=sumstats["DENSITY"].mean()
        bmedian=sumstats["DENSITY"].median()
        log.write(" -DENSITY column exists. Skipping calculation...",verbose=verbose) 
    return   sumstats, bmean, bmedian