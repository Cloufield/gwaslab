from gwaslab.viz.viz_aux_quickfix import _get_largenumber
import pandas as pd
from gwaslab.util.util_in_get_density import getsignaldensity2

def _process_density(sumstats, mode, bwindowsizekb, chrom, pos, verbose, log):
    """
    Create Brisbane plot (variant density plot)
    Y axis: number of neiboring variants within a window around each variant)
    X axis: same as Manhattan plot

    Parameters
    ---
    bwindowsizekb : int, default=100
        Window size in kb for Brisbane plot density calculation. (Used in 'b' mode)
    density_color : bool, default=False
        Whether to color Brisbane plot points by density. (Used in 'b' mode)
    density_range : tuple, default=(0,15)
        Color range for density plot. (Used in 'b' mode)
    density_trange : tuple, default=(0,10)
        Threshold range for density plot. (Used in 'b' mode)
    density_threshold : int, default=5
        Threshold for density coloring. (Used in 'b' mode)
    density_tpalette : str, default='Blues'
        Palette for thresholded density colors. (Used in 'b' mode)
    density_palette : str, default='Reds'
        Palette for density colors. (Used in 'b' mode)
    anno : boolean, str or 'GENENAME', default=None.
        Specify which data to use for annotation. Default is None. 
        anno options:
            None: no annotation
            True: annotate variants with chromosome and position like chr:pos
            'GENENAME': annotate variants with closest gene names
            str: annotate variants with values in Column with the header 
    anno_set : list, optional
        Set of variants IDs to annotate. Default is None. 
    anno_alias : dict, optional
        Dictionary mapping SNP IDs to custom annotation labels. Default is None. 
    anno_d : dict, optional
       Dictionary mapping annotation indices to custom positioning options. 
       For example, {0:"l"} means adjust the direction of the 1st arm end to left. Default is empty dict. Default is None. 
    anno_args : dict, optional
        ictionary of default styling arguments for annotations. Default is None. 
    anno_args_single : dict, optional
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
    save_args : dict, default={"dpi":600,"transparent":True}
        Arguments for saving the plot. 
    """
    if "b" in mode and "DENSITY" not in sumstats.columns:

        sumstats["DENSITY"] = getsignaldensity2(insumstats=sumstats,
                                                id="SNPID",
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