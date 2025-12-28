from gwaslab.viz.viz_aux_quickfix import _get_largenumber
import pandas as pd
from gwaslab.util.util_in_get_density import _get_signal_density2

def _process_density(sumstats, mode, bwindowsizekb, chrom, pos, verbose, log):
    """
    Create Brisbane plot (variant density plot)
    Y axis: number of neiboring variants within a window around each variant)
    X axis: same as Manhattan plot

    Parameters
    ---
    bwindowsizekb : int, default=100
        Window size in kb for Brisbane plot density calculation.
    colors : list, default=['#597FBD', '#74BAD3']
        Color palette used to assign colors to each chromosome in the plot.
        If the number of chromosomes exceeds the number of colors provided,
        the palette will be cycled automatically.
    density_color : bool, default=False
        Whether to color Brisbane plot points by density. If true, `density_palette` override `colors`.
    density_palette : str, default='Reds'
        Palette for density colors. 
    density_range : tuple, default=None.
        Color range for density plot. If None, range will auto-selected as (to_plot["DENSITY"].min(), to_plot["DENSITY"].max()).
    sig_level : int, optional
        Threshold for variants density for annotation. Default is 1.
    windowsizekb : int, default=500
        Window size in kb for determine lead variants for annotation. 
    anno : boolean, str or 'GENENAME', default=None
        Specify which data to use for annotation. Default is None. 
        anno options:
        - None: no annotation
        - True: annotate variants with chromosome and position like chr:pos
        - 'GENENAME': annotate variants with closest gene names
        - str: annotate variants with values in Column with the header
    anno_set : list, optional
        Set of variants IDs to annotate. Default is None. 
    anno_alias : dict, optional
        Dictionary mapping SNP IDs to custom annotation labels. Default is None. 
    anno_d : dict, optional
       Dictionary mapping annotation indices to custom positioning options. 
       For example, {0:"l"} means adjust the direction of the 1st arm end to left. Default is empty dict. Default is None. 
    anno_kwargs : dict, optional
        ictionary of default styling arguments for annotations. Default is None. 
    anno_kwargs_single : dict, optional
        Dictionary mapping SNP IDs to custom styling arguments. Default is None.
    anno_style : 'right' or 'tight' or 'expand', default='right'
        Style of annotation ('right', 'tight', 'expand'). Default is 'right'.
    anno_fixed_arm_length : float, optional
        Fixed arm length for annotations. Default is None. 
    anno_source : str, default='ensembl'
        Source for annotations.
    anno_xshift : float, optional
        X-axis shift for annotations. Default is None. 
    anno_max_iter : int, default=100
        Maximum iterations for text repulsion algorithm. 
    arrow_kwargs : dict, optional
        Arguments for annotation arrows. Default is None.
    arm_offset : float, optional
        Offset for annotation arms. Default is None. 
    arm_scale : float, default=1
        Scaling factor for annotation arm length. 
    anno_height : float, default=1
        Height for annotations. 
    arm_scale_d : dict, optional
        Dictionary mapping annotation indices to custom arm scaling factors.
    pinpoint : list, optional
        List of variants to pinpoint. Default is None. (Used in all modes)
    pinpoint_color : str, default='red'
        Color for pinpointed variants. (Used in all modes)
    verbose : bool, default=True
        Whether to show progress. 
    build : str, optional
        Genomic build version. Default is None. 
    dpi : int, default=200
        Dots per inch for figure resolution. 
    save : str, optional
        File path to save plot. Default is None. 
    save_kwargs : dict, default={"dpi":600,"transparent":True}
        Arguments for saving the plot. 

    
    Less used parameters
    ----------------------
    density_threshold : int, default=5
        Threshold for density coloring. Different coloring makes it clear to see low and high density regions.
        Above density_threshold, coloring using density_range and density_palette. 
        Below density_threshold, coloring using density_trange and density_tpalette. 
    density_trange : tuple, default=(0,10)
        Threshold range for density plot.  For hue_norm.
    density_tpalette : str, default='Blues'
        Palette for thresholded density colors.
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