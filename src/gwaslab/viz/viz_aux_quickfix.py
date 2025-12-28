import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from matplotlib import ticker
import matplotlib.pyplot as plt
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from math import ceil

def _quick_fix(sumstats, chr_dict=get_chr_to_number(), scaled=False, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",log=Log(), verbose=True):
    '''
    quick sanity check for input sumstats
    '''

    # quick fix chr
    sumstats[chrom] = _quick_fix_chr(sumstats[chrom], chr_dict=chr_dict, verbose=verbose, log=log)

    # quick fix pos
    sumstats[pos] = _quick_fix_pos(sumstats[pos], verbose=verbose, log=log)

    # quick fix p
    sumstats = _quick_fix_p_value(sumstats, p=p, mlog10p=mlog10p, scaled=scaled, verbose=verbose, log=log)

    # quick fix mlog10p
    sumstats = _quick_fix_mlog10p(sumstats, p=p, mlog10p=mlog10p, scaled=scaled, verbose=verbose, log=log)

    return sumstats


def _quick_fix_p_value(sumstats, p="P", mlog10p="MLOG10P", scaled=False,verbose=True, log=Log()):
    '''
    drop variants with bad P values
    '''
    if scaled==True:
        # if scaled, add scaled P and P col 
        log.write(" -P values are already scaled...", verbose=verbose)
        log.write(" -Sumstats -log10(P) values are being converted to P...", verbose=verbose)
        sumstats["scaled_P"] = sumstats[mlog10p]
        sumstats[p]= np.power(10,-sumstats[mlog10p].astype("float64"))
        return sumstats
    # bad p : na and outside (0,1]
    bad_p_value = (sumstats[p].isna()) | (sumstats[p] > 1) | (sumstats[p] <= 0)
    bad_p_count = sum(bad_p_value)

    log.write(" -Sanity check after conversion: " + str(bad_p_count) +
                  " variants with P value outside of (0,1] will be removed...", verbose=verbose if bad_p_count > 0 else False)
    sumstats = sumstats.loc[~bad_p_value, :]
    return sumstats


def _quick_fix_mlog10p(insumstats,p="P", mlog10p="MLOG10P", scaled=False, log=Log(), verbose=True):
    '''
    drop variants with bad -log10(P) values
    '''
    if scaled != True:
        log.write(" -Sumstats P values are being converted to -log10(P)...", verbose=verbose)
        # Use pd.to_numeric for better performance and error handling
        p_values = pd.to_numeric(insumstats[p], errors='coerce')
        insumstats["scaled_P"] = -np.log10(p_values)
    elif "scaled_P" not in insumstats.columns:
        # If scaled but scaled_P doesn't exist, create it from mlog10p
        insumstats["scaled_P"] = pd.to_numeric(insumstats[mlog10p], errors='coerce')
    
    # More efficient: use np.isfinite to check for both inf and nan in one operation
    scaled_P = insumstats["scaled_P"]
    if_inf_na = ~np.isfinite(scaled_P)
    inf_na_count = if_inf_na.sum()

    log.write(" -Sanity check: "+str(inf_na_count) +
                  " na/inf/-inf variants will be removed...", verbose=verbose if inf_na_count > 0 else False)
    sumstats = insumstats.loc[~if_inf_na, :]
    return sumstats


def _quick_fix_eaf(seires,log=Log(), verbose=True):
    '''
    conversion of eaf to maf
    '''
    # Early exit: Check if already numeric
    if pd.api.types.is_numeric_dtype(seires):
        log.write(" -EAF data type is already numeric. Skipping conversion...", verbose=verbose)
        # Still need to flip values > 0.5
        flipped = seires > 0.5
        seires = seires.copy()
        seires[flipped] = 1 - seires[flipped]
        return seires
    
    seires = pd.to_numeric(seires, errors='coerce')
    flipped = seires > 0.5
    seires[flipped] = 1 - seires[flipped]
    return seires


def _quick_fix_chr(seires, chr_dict,log=Log(), verbose=True):
    '''
    conversion and check for chr
    '''
    # Early exit: Check if already numeric
    if pd.api.types.is_numeric_dtype(seires):
        log.write(" -CHR data type is already numeric. Skipping conversion...", verbose=verbose)
        return seires.astype('Int64')
    
    if pd.api.types.is_string_dtype(seires):
        # if chr is string dtype: convert using chr_dict
        seires = seires.map(chr_dict, na_action="ignore")
    seires = np.floor(pd.to_numeric(seires, errors='coerce')).astype('Int64')
    return seires


def _quick_fix_pos(seires,log=Log(), verbose=True):
    '''
    force conversion for pos
    '''
    # Early exit: Check if already numeric
    if pd.api.types.is_numeric_dtype(seires):
        log.write(" -POS data type is already numeric. Skipping conversion...", verbose=verbose)
        return np.floor(seires).astype('Int64')
    
    seires = np.floor(pd.to_numeric(seires, errors='coerce')).astype('Int64')
    return seires




from gwaslab.qc.qc_normalize_args import _normalize_region

def _dropna_in_cols(sumstats, cols, log=Log(), verbose=True):
    to_drop = sumstats[cols].isna().any(axis=1)
    log.write(" -Dropping {} variants due to missing values in {}.".format(sum(to_drop),cols))
    return sumstats.loc[~to_drop,:]


def _get_largenumber(*args,log=Log(), verbose=True):
    '''
    get a helper large number, >> max(pos)
    '''
    large_number = 1000000000
    max_number = np.nanmax(args)
    for i in range(7):
        if max_number*10 > large_number:
            large_number = int(large_number * 10)
        else:
            break
        if i == 7:
            log.warning("Max POS is too large!")
    return large_number


def _quick_add_tchrpos(sumstats, chr="chr", pos="POS", large_number=10000000000, dropchrpos=False,log=Log(), verbose=True):
    sumstats["TCHR+POS"] = sumstats["CHR"]*large_number + sumstats["POS"]
    sumstats["TCHR+POS"] = sumstats["TCHR+POS"].astype('Int64')
    if dropchrpos == True:
        sumstats = sumstats.drop(labels=["CHR", "POS"], axis=1)
    sumstats = sumstats.dropna()
    return sumstats


def _quick_merge_sumstats(sumstats1, sumstats2, log=Log(), verbose=True):
    merged_sumstats = pd.merge(sumstats1, sumstats2, on="TCHR+POS", how="outer", suffixes=('_1', '_2'))
    merged_sumstats["CHR"] = merged_sumstats["CHR_1"]
    merged_sumstats["POS"] = merged_sumstats["POS_1"]
    merged_sumstats.loc[merged_sumstats["CHR_1"].isna(), "CHR"] = merged_sumstats.loc[merged_sumstats["CHR_1"].isna(), "CHR_2"]
    merged_sumstats.loc[merged_sumstats["POS_1"].isna(), "POS"] = merged_sumstats.loc[merged_sumstats["POS_1"].isna(), "POS_2"]
    merged_sumstats = merged_sumstats.drop(labels=["CHR_1", "CHR_2", "POS_1", "POS_2"],axis=1)
    return merged_sumstats

def _quick_assign_i(sumstats, chrom="CHR",pos="POS",log=Log(), verbose=True):
    # sort by CHR an POS
    sumstats = sumstats.sort_values([chrom,pos])
    # set new id
    sumstats["_ID"]=range(len(sumstats))
    sumstats = sumstats.set_index("_ID")
    #create a df , groupby by chromosomes , and get the maximum position
    posdic = sumstats.groupby(chrom)[pos].max()
    # convert to dictionary
    posdiccul = posdic.to_dict()
    # fill empty chr with 0
    max_chr = sumstats[chrom].max()
    posdiccul = {i: posdiccul.get(i, 0) for i in range(max_chr + 1)}
    # get interval
    interval_between_chr = sumstats[pos].max()*0.05
    # cumulative sum dictionary
    for i in range(1, max_chr + 1):
        posdiccul[i] = posdiccul[i-1] + posdiccul[i] + interval_between_chr
    # convert base pair postion to x axis position using the cumulative sum dictionary
    # Use vectorized map with pre-computed mapping Series
    chrom_int = sumstats[chrom].astype(int)
    add_mapping = pd.Series(posdiccul)
    sumstats["_ADD"] = (chrom_int - 1).map(add_mapping)
    sumstats["i"]   = sumstats[pos] + sumstats["_ADD"]
    # drop add
    sumstats = sumstats.drop(labels=["_ADD"],axis=1)
    #for plot, get the chr text tick position  
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    # fix dtype for i
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
    return sumstats, chrom_df

def _quick_assign_i_with_rank(sumstats, chrpad, use_rank=False, chrom="CHR",pos="POS",drop_chr_start=False,_posdiccul=None,log=Log(), verbose=True):
    # align all variants on a single axis (i)
    sumstats = sumstats.sort_values([chrom,pos])
    if use_rank is True: 
        sumstats["_POS_RANK"] = sumstats.groupby(chrom)[pos].rank("dense", ascending=True)
        pos="_POS_RANK"
    sumstats["_ID"]=range(len(sumstats))
    sumstats=sumstats.set_index("_ID")

    #create a df , groupby by chromosomes , and get the maximum position
    if use_rank is True: 
        posdic = sumstats.groupby(chrom)["_POS_RANK"].max()
    else:
        posdic = sumstats.groupby(chrom)[pos].max()
    
    if _posdiccul is None:
        # convert to dictionary
        posdiccul = posdic.to_dict()
        
        # fill empty chr with 0
        max_chr = sumstats[chrom].max()
        posdiccul = {i: posdiccul.get(i, 0) for i in range(max_chr + 1)}

        # cumulative sum dictionary
        for i in range(1, max_chr + 1):
            posdiccul[i] = posdiccul[i-1] + posdiccul[i] + sumstats[pos].max()*chrpad
    else:
        posdiccul = _posdiccul
    # convert base pair postion to x axis position using the cumulative sum dictionary
    # Use vectorized map with pre-computed mapping Series
    chrom_int = sumstats[chrom].astype(int)
    add_mapping = pd.Series(posdiccul)
    sumstats["_ADD"] = (chrom_int - 1).map(add_mapping)
    
    if drop_chr_start==True:
            posdic_min = sumstats.groupby(chrom)[pos].min()
            posdiccul_min = posdic_min.to_dict()
            max_chr = sumstats[chrom].max()
            posdiccul_min = {i: posdiccul_min.get(i, 0) for i in range(max_chr + 1)}
            for i in range(1, max_chr + 1):
                posdiccul_min[i] = posdiccul_min[i-1] + posdiccul_min[i]
            min_mapping = pd.Series(posdiccul_min)
            sumstats["_ADD"] = sumstats["_ADD"] - (chrom_int - 1).map(min_mapping)
        
    if use_rank is True: 
        sumstats["i"]=sumstats["_POS_RANK"]+sumstats["_ADD"]
    else:
        sumstats["i"]=sumstats[pos]+sumstats["_ADD"]
    

    #for plot, get the chr text tick position      
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    #sumstats["i"] = sumstats["i"]+((sumstats[chrom].map(dict(chrom_df)).astype("int")))*0.02
    #sumstats["i"] = sumstats["i"].astype("Int64")
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
    return sumstats, chrom_df

def _quick_assign_marker_relative_size(series, sig_level = 5e-8, suggestive_sig_level=5e-6, lower_level=5e-4,log=Log(), verbose=True):
    size_series = series.copy()
    size_series[:] = 1

    is_lower_level          = series > -np.log10(lower_level)
    is_suggestive_sig_level = series > -np.log10(suggestive_sig_level)
    is_sig_level            = series > -np.log10(sig_level)
    
    size_series[is_lower_level] = 2
    size_series[is_suggestive_sig_level] = 3
    size_series[is_sig_level] = 4
    return size_series

def _quick_assign_highlight_hue(sumstats,highlight,highlight_windowkb, snpid="SNPID",chrom="CHR",pos="POS",log=Log(), verbose=True):
    # Initialize HUE column if it doesn't exist
    if "HUE" not in sumstats.columns:
        sumstats["HUE"] = None
    
    # Get highlighted variants
    to_highlight = sumstats.loc[sumstats[snpid].isin(highlight), [chrom, pos]]
    
    if len(to_highlight) == 0:
        return sumstats
    
    # Vectorized approach: create mask for all highlighted regions at once
    highlight_mask = pd.Series(False, index=sumstats.index)
    
    for target_chr, target_pos in zip(to_highlight[chrom], to_highlight[pos]):
        target_chr = int(target_chr)
        target_pos = int(target_pos)
        window = highlight_windowkb * 1000
        # Vectorized boolean mask
        mask = (sumstats[chrom] == target_chr) & \
               (sumstats[pos] > target_pos - window) & \
               (sumstats[pos] < target_pos + window)
        highlight_mask |= mask
    
    sumstats.loc[highlight_mask, "HUE"] = "0"
    return sumstats

def _quick_assign_highlight_hue_pair(sumstats, highlight1, highlight2, highlight_windowkb, chrom="CHR",pos="POS",log=Log(), verbose=True):
    #assign colors: 0 is hightlight color
    to_highlight1 = pd.DataFrame()
    to_highlight2 = pd.DataFrame()
    
    if len(highlight1) > 0:
        to_highlight1 = sumstats.loc[sumstats["TCHR+POS"].isin(highlight1), [chrom, pos]]
    if len(highlight2) > 0:
        to_highlight2 = sumstats.loc[sumstats["TCHR+POS"].isin(highlight2), [chrom, pos]]
    
    # Initialize HUE columns if they don't exist
    if "HUE1" not in sumstats.columns:
        sumstats["HUE1"] = None
    if "HUE2" not in sumstats.columns:
        sumstats["HUE2"] = None
    
    window = highlight_windowkb * 1000
    
    if len(to_highlight1) > 0:
        highlight1_mask = pd.Series(False, index=sumstats.index)
        for target_chr, target_pos in zip(to_highlight1[chrom], to_highlight1[pos]):
            target_chr = int(target_chr)
            target_pos = int(target_pos)
            mask = (sumstats["CHR"] == target_chr) & \
                   (sumstats["POS"] > target_pos - window) & \
                   (sumstats["POS"] < target_pos + window)
            highlight1_mask |= mask
        sumstats.loc[highlight1_mask, "HUE1"] = "0"
    
    if len(to_highlight2) > 0:
        highlight2_mask = pd.Series(False, index=sumstats.index)
        for target_chr, target_pos in zip(to_highlight2[chrom], to_highlight2[pos]):
            target_chr = int(target_chr)
            target_pos = int(target_pos)
            mask = (sumstats["CHR"] == target_chr) & \
                   (sumstats["POS"] > target_pos - window) & \
                   (sumstats["POS"] < target_pos + window)
            highlight2_mask |= mask
        sumstats.loc[highlight2_mask, "HUE2"] = "0"
    
    return sumstats, to_highlight1, to_highlight2

def _quick_extract_snp_in_region(sumstats, region, chrom="CHR",pos="POS",log=Log(), verbose=True):
    region_chr = region[0]
    region_start = region[1]
    region_end = region[2]
    log.write(" -Extract SNPs in region : chr"+str(region_chr)+":"+str(region[1])+"-"+str(region[2])+ "...", verbose=verbose)
    is_in_region_snp = (sumstats[chrom]==region_chr) &(sumstats[pos]<region_end) &(sumstats[pos]>region_start)
    log.write(" -Extract SNPs in specified regions: "+str(sum(is_in_region_snp)), verbose=verbose)
    sumstats = sumstats.loc[is_in_region_snp,:]
    return sumstats

def _cut(series, mode,cutfactor,cut,skip, ylabels, cut_log, verbose, lines_to_plot, log):
    """
    Shrink/compress extremely large y-axis values to prevent them from dominating the plot.
    
    This function applies a "cut" transformation to values above a threshold, compressing them
    so that very high values (e.g., extremely significant p-values) don't take up most of the
    plot space. Two modes are available:
    - Linear mode: Values above cut are shrunk by dividing by cutfactor
    - Log mode: Values above cut are compressed using logarithmic scaling
    
    Parameters
    ----------
    series : pd.Series
        The y-axis values to transform (typically -log10(P) values)
    mode : str
        Plot mode (e.g., "mqq", "b" for density plot)
    cutfactor : float
        Shrinkage factor for linear mode (default 10)
    cut : float or bool
        Threshold above which values are shrunk. If True, auto-determines cut value
    skip : float
        Minimum value to plot (values below skip are omitted)
    ylabels : list or None
        Custom y-axis tick labels to also transform
    cut_log : bool
        If True, use logarithmic compression; if False, use linear shrinkage
    verbose : bool
        Whether to print progress messages
    lines_to_plot : array-like
        Additional lines (e.g., significance thresholds) to also transform
    log : Log
        Logging object
    
    Returns
    -------
    tuple
        (transformed_series, maxy, maxticker, cut, cutfactor, ylabels, lines_to_plot)
    """
    log.write(" -Converting data above cut line...",verbose=verbose)
    
    # Step 1: Prepare inputs - convert ylabels to Series if provided
    if ylabels is not None:
        ylabels = pd.Series(ylabels)
    
    # Step 2: Create a copy of the series to avoid modifying the original
    series = series.copy()
    
    # Step 3: Find the maximum value in the series for reporting and calculations
    maxy = series.max()
    if "b" not in mode:
        log.write(" -Maximum -log10(P) value is "+str(maxy) +" .", verbose=verbose)
    elif "b" in mode:
        log.write(" -Maximum DENSITY value is "+str(maxy) +" .", verbose=verbose)
    
    # Step 4: Calculate the maximum ticker value (rounded integer for display purposes)
    maxticker=int(np.round(series.max(skipna=True)))
    
    # Step 5: Process cut transformation if cut is specified
    if cut:
        # Step 5a: Auto mode - automatically determine cut threshold and cutfactor
        if cut==True:
            log.write(" -Cut Auto mode is activated...", verbose=verbose)
            # If maximum value is less than 30, no need to cut (values are not extreme)
            if maxy<30:
                log.write(" - maxy <30 , no need to cut.", verbose=verbose)
                cut=0
            else:
                # Set cut threshold to 20 and calculate cutfactor to fit remaining range
                # Formula: cutfactor = (maxy - cut) / 8 ensures compressed range fits in ~8 units
                cut = 20
                cutfactor = ( maxy - cut )/8
        
        # Step 5b: Apply the cut transformation if cut threshold is still set
        if cut:
            # Step 5b-i: Logarithmic compression mode
            if cut_log==True:
                # Recalculate maxticker for log mode calculations
                maxticker=int(np.round(series.max(skipna=True)))
                
                # Calculate amplitude factor for logarithmic scaling
                # This factor determines how much the log-compressed range will span
                # Formula: amp = (cut - skip) / 2 / log2(maxticker/cut)
                # The division by 2 ensures the compressed range is half the original range
                amp = (cut - skip)/ 2 / np.log2(maxticker/cut)
                
                # Transform data values above cut using log2 compression
                # Formula: new_value = log2(old_value/cut) * amp + cut
                # This compresses values logarithmically while keeping cut as the baseline
                series[series>cut] = (np.log2(series[series>cut]/cut)) * amp + cut
                
                # Transform y-axis labels if provided (same log compression)
                if ylabels is not None:
                    ylabels[ylabels>cut] = (np.log2(ylabels[ylabels>cut]/cut)) * amp +cut 
                
                # Transform additional lines (e.g., significance thresholds) using same compression
                lines_to_plot[lines_to_plot>cut] = (np.log2(lines_to_plot[lines_to_plot>cut]/cut)) * amp +cut 
                
                # Calculate the new maximum y value after log compression
                # This represents the compressed maximum value for setting y-axis limits
                maxy = (np.log2(maxticker) - np.log2(cut)) * amp + cut
            else:
                # Step 5b-ii: Linear shrinkage mode (default)
                if "b" not in mode:
                    log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...", verbose=verbose)
                else:
                    log.write(" -DENSITY values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...", verbose=verbose)

                # Recalculate maxticker for linear mode
                maxticker=int(np.round(series.max(skipna=True)))

                # Transform data values above cut using linear shrinkage
                # Formula: new_value = (old_value - cut) / cutfactor + cut
                # This shrinks values by dividing the excess (above cut) by cutfactor
                # Example: if cut=20, cutfactor=10, value=100 becomes (100-20)/10+20 = 28
                series[series>cut] = (series[series>cut]-cut)/cutfactor+cut
                
                # Transform y-axis labels if provided (same linear shrinkage)
                if ylabels is not None:
                    ylabels[ylabels>cut] = (ylabels[ylabels>cut]-cut)/cutfactor+cut
                
                # Transform additional lines using same linear shrinkage
                lines_to_plot[lines_to_plot>cut] = (lines_to_plot[lines_to_plot>cut]-cut)/cutfactor+cut
                
                # Calculate the new maximum y value after linear shrinkage
                # Formula: maxy = (maxticker - cut) / cutfactor + cut
                maxy = (maxticker-cut)/cutfactor + cut

    # Step 6: Return transformed values and parameters
    return series, maxy, maxticker, cut, cutfactor,ylabels,lines_to_plot

#def _cut_line(level, mode,cutfactor,cut,skip, ylabels, cut_log, verbose, log):




def _set_yticklabels(cut,
                     cutfactor,
                     cut_log,
                     ax1,
                     skip,
                     maxy,
                     maxticker,
                     ystep,
                     sc_linewidth,
                     cut_line_color,
                     fontsize,
                     font_family,
                     ytick3,
                     ylabels,
                     ylabels_converted,
                     log=Log(),
                     verbose=True
                     ):
    """
    Set y-axis tick positions and labels for plots with optional cut transformation.
    
    This function handles the complex task of setting y-axis ticks and labels when a "cut"
    transformation has been applied. It creates separate tick sets for:
    - The region below the cut threshold (normal scale)
    - The region above the cut threshold (compressed scale)
    
    For the compressed region, labels are reverse-transformed to show original values,
    while tick positions use the compressed coordinates.
    
    Parameters
    ----------
    cut : float
        Threshold above which values were compressed (0 means no cut)
    cutfactor : float
        Shrinkage factor used in linear compression mode
    cut_log : bool
        Whether logarithmic compression was used (True) or linear (False)
    ax1 : matplotlib.axes.Axes
        The axes object to modify
    skip : float
        Minimum y-axis value (values below this are not plotted)
    maxy : float
        Maximum y-axis value after transformation
    maxticker : int
        Maximum original value before transformation (for label reverse-transformation)
    ystep : float
        Step size for y-axis ticks (0 means auto-calculate)
    sc_linewidth : float
        Line width for the cut line
    cut_line_color : str
        Color for the cut line
    fontsize : float
        Font size for tick labels
    font_family : str
        Font family for tick labels
    ytick3 : bool
        Whether to show intermediate ticks in the compressed region
    ylabels : list or None
        Custom y-axis labels (if provided, overrides auto-generated labels)
    ylabels_converted : array-like
        Converted positions for custom ylabels (if ylabels is provided)
    log : Log
        Logging object
    verbose : bool
        Whether to print progress messages
    
    Returns
    -------
    matplotlib.axes.Axes
        The modified axes object
    """
    log.write(" -Processing Y tick labels...",verbose=False)
    
    # Step 1: Handle case with no cut transformation
    # If cut==0, no compression was applied, so just set simple y-axis limits
    if cut == 0: 
        # Set y-axis limits from skip to maxy with 20% padding at top
        ax1.set_ylim((skip, ceil(maxy*1.2)) )
    
    # Step 2: Handle case with cut transformation
    if cut!=0:
        # Step 2a: Draw a horizontal line at the cut threshold to visually indicate the compression
        # This line shows where the scale changes from normal to compressed
        cutline = ax1.axhline(y=cut, linewidth = sc_linewidth,linestyle="--",color=cut_line_color,zorder=1)
        
        # Step 2b: Determine step size for tick spacing
        # Default step is 2, but can be auto-calculated or user-specified
        step=2

        if ystep == 0:
            # Auto-calculate step: if there would be more than 8 ticks, increase step size
            # This prevents too many ticks from cluttering the axis
            if (cut - skip ) // step > 8:
                step = (cut - skip ) // 8
        else:
            # Use user-specified step size
            step = ystep

        # Step 2c: Determine upper bound for ticks below cut threshold
        # If (cut-skip) is evenly divisible by step, set upper to cut-1 to avoid overlap
        # Otherwise, set upper to cut to include the cut line
        if (cut-skip)%step==0:
            upper = cut - 1
        else:
            upper = cut
        
        # Step 2d: Generate ticks and labels for the region BELOW the cut threshold (normal scale)
        # ticks1: Regular ticks from skip to upper with step spacing
        ticks1= [x for x in range(skip,upper,step)]
        # ticks2: The cut line itself (always shown)
        ticks2= [cut]
        # Labels are the same as tick positions for the normal scale region
        tickslabel1= [x for x in range(skip,upper,step)]
        tickslabel2= [cut]
        
        # Step 2e: Generate ticks for the region ABOVE the cut threshold (compressed scale)
        if cut_log==True:
            # Logarithmic compression mode: calculate tick spacing based on log scale
            # Step size is determined by log2 of the ratio between max and cut
            # max(..., 1) ensures at least 1 step to avoid division by zero
            ticks3= [x for x in range(cut,int(maxy),max(int(np.log2(maxticker//(cut-skip))),1))] 
            ticks4= [int(maxy)]  # Always include the maximum value
        elif ytick3 == True:
            # Linear compression mode with intermediate ticks enabled
            # Calculate compressed maximum: (maxticker-cut)/cutfactor + cut
            # Generate ticks from cut to compressed max using step spacing
            ticks3= [x for x in range(cut,int((maxticker-cut)/cutfactor + cut),int(step))]
            ticks4= [(maxticker-cut)/cutfactor + cut]  # Compressed maximum
        else:
            # Linear compression mode without intermediate ticks
            # Only show the compressed maximum, no intermediate ticks
            ticks3= []  # No intermediate ticks
            ticks4= [(maxticker-cut)/cutfactor + cut]  # Only the compressed maximum
        
        # Step 2f: Generate LABELS for the compressed region (reverse-transform to show original values)
        # Labels show the original values, not the compressed positions
        if cut_log==True:
            # Logarithmic mode: reverse the log transformation to get original values
            # Formula: original = 2^((compressed_pos - cut) / amp) * cut
            # amp is the same amplitude factor used in _cut function
            amp = (cut - skip)/ 2 / np.log2(maxticker/cut) 
            # For each compressed tick position, calculate the original value it represents
            tickslabel3 = list(map(lambda x: int(2**((x-cut)/amp)* cut) ,ticks3))
        elif ytick3 == True:
            # Linear mode: reverse the linear transformation to get original values
            # Formula: original = (compressed_pos - cut) * cutfactor + cut
            # This is the inverse of the compression formula: (value-cut)/cutfactor+cut
            tickslabel3 = list(map(lambda x: int((x-cut)*cutfactor + cut) ,ticks3))
        else:
            # No intermediate ticks, so no intermediate labels
            tickslabel3=[]
        # Label for the maximum is always the original maxticker value
        tickslabel4= [maxticker]

        # Step 2g: Apply ticks and labels to the axis
        if maxy > cut:
            # If there are values above cut, show all tick sets (below + at + above cut)
            ax1.set_yticks(ticks1+ticks2+ticks3+ticks4)
            ax1.set_yticklabels(tickslabel1+tickslabel2+tickslabel3+tickslabel4,fontsize=fontsize,family=font_family)
        else:
            # If no values above cut, only show ticks for the region below cut
            ax1.set_yticks(ticks1+ticks2)
            ax1.set_yticklabels(tickslabel1+tickslabel2,fontsize=fontsize,family=font_family)
    
    # Step 3: Handle custom y-axis labels if provided
    # Custom labels override the auto-generated labels
    if ylabels is not None:  
        # Use the converted positions (already transformed) for tick positions
        ax1.set_yticks(ylabels_converted)
        # Use the original labels (user-provided) for tick labels
        ax1.set_yticklabels(ylabels,fontsize=fontsize,family=font_family)
    
    # Step 4: Set final y-axis bounds
    # Get the current top limit (may have been set by previous operations)
    ylim_top = ax1.get_ylim()[1]
    # Set lower bound to skip and keep the existing upper bound
    ax1.set_ybound(lower=skip,upper=ylim_top)
    
    # Step 5: Apply font size to y-axis tick labels
    ax1.tick_params(axis='y', labelsize=fontsize)

    return ax1

def _jagged_y(cut,skip,ax1,mode,mqqratio,jagged_len,jagged_wid, log=Log(), verbose=True):
    log.write(" -Processing jagged Y axis...",verbose=verbose)
    tycut = cut +0.3 #(cut - skip)/ (ax1.get_ylim()[1] - skip) + 0.002
    dy= jagged_len * (cut - skip) 
    x0 =  0
    dx= jagged_wid
    if mode>1:
        dx = dx * mqqratio
    kwargs = dict(transform=ax1.get_yaxis_transform(), color='k', clip_on=False,solid_capstyle="round",linewidth=0.8)
    ax1.plot((x0,x0), (tycut,tycut+4*dy), zorder=1000, transform=ax1.get_yaxis_transform(), color="white",clip_on=False,solid_capstyle="butt")
    ax1.plot((x0,-dx), (tycut,tycut+dy), zorder=1001, **kwargs)
    ax1.plot((-dx,+dx), (tycut+dy,tycut+3*dy), zorder=1001, **kwargs)
    ax1.plot((+dx,x0), (tycut+3*dy,tycut+4*dy), zorder=1001,  **kwargs)
    return ax1
