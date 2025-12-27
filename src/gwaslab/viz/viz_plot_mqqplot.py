import copy
import gc as garbage_collect
import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy as sp
from math import ceil
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy import stats

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version
from gwaslab.io.io_process_kwargs import _update_arg
from gwaslab.io.io_process_kwargs import _update_kwargs
from gwaslab.qc.qc_build import _process_build
from gwaslab.util.util_in_filter_value import _filter_region
from gwaslab.util.util_in_get_sig import _get_sig, _anno_gene
from gwaslab.viz.viz_aux_annotate_plot import annotate_single
from gwaslab.viz.viz_aux_quickfix import _cut
from gwaslab.viz.viz_aux_quickfix import _jagged_y
from gwaslab.viz.viz_aux_quickfix import _set_yticklabels
from gwaslab.qc.qc_normalize_args import _normalize_group
from gwaslab.qc.qc_normalize_args import _normalize_region
from gwaslab.viz.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz.viz_aux_quickfix import _quick_fix_chr
from gwaslab.viz.viz_aux_quickfix import _quick_fix_eaf
from gwaslab.viz.viz_aux_quickfix import _quick_fix_mlog10p
from gwaslab.viz.viz_aux_quickfix import _quick_fix_p_value
from gwaslab.viz.viz_aux_quickfix import _quick_fix_pos
from gwaslab.viz.viz_aux_save_figure import safefig
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_plot_density import _process_density
from gwaslab.viz.viz_plot_manhattan_like import _configure_fig_save_kwargs, _add_pad_to_x_axis, _configure_cols_to_use, _sanity_check, _process_p_value, _process_highlight, _process_line, _process_cbar, _process_xtick, _process_ytick, _process_xlabel, _process_ylabel, _process_spine, _process_layout
import gwaslab.viz.viz_plot_manhattan_like as manlike
from gwaslab.viz.viz_plot_manhattan_mode import draw_manhattan_panel
from gwaslab.viz.viz_plot_qqplot import _plot_qq
from gwaslab.viz.viz_plot_regional2 import _get_lead_id, _plot_regional, prepare_vcf_context, process_ld, process_vcf
from gwaslab.io.io_process_kwargs import normalize_series_inputs

@normalize_series_inputs(keys=["highlight","pinpoint","anno_set"])
def _setup_and_log_mqq_info(
    insumstats,
    mode,
    sig_level,
    anno_set,
    highlight,
    highlight_chrpos,
    highlight_color,
    highlight_windowkb,
    pinpoint,
    pinpoint_color,
    region,
    build,
    log,
    verbose,
):
    # Derive a concise label for the chosen mode
    if "m" in mode and "qq" in mode:
        plot_label = "Manhattan-QQ plot"
    elif mode == "m":
        plot_label = "Manhattan plot"
    elif mode == "qq":
        plot_label = "QQ plot"
    elif mode == "r":
        plot_label = "Region plot"
    elif mode == "b":
        plot_label = "Density plot"
    else:
        plot_label = "MQQ plot"

    # Header and build info
    log.write("Starting {} creation (Version {})".format(plot_label, _get_version()), verbose=verbose)
    build = _update_arg(build, "19")
    build = _process_build(build, log=log, verbose=verbose)
    log.write(" - Genomic coordinates version: {} ...".format(build), verbose=verbose)
    if build is None or build == "99":
        log.warning("Genomic coordinates version is unknown.")

    # Basic stats and mode
    log.write(" - Genome-wide significance level to plot is set to " + str(sig_level) + " ...", verbose=verbose)
    log.write(" - Input sumstats contains {} variants...".format(len(insumstats)), verbose=verbose)
    log.write(" - {} layout mode selected: {}".format(plot_label, mode), verbose=verbose)

    # Annotation set (only meaningful when Manhattan panel is present)
    if anno_set is not None and len(anno_set) > 0 and ("m" in mode):
        log.write(" -Variants to annotate : " + ",".join(anno_set), verbose=verbose)

    # Highlight configuration (supports lists of loci or direct chr:pos sets)
    if highlight is None:
        highlight = []
    if len(highlight) > 0 and ("m" in mode):
        highlight, highlight_color = _normalize_group(
            items=highlight,
            colors=highlight_color,
            item_name="highlight loci",
            log=log,
            verbose=verbose,
            is_chrpos_mode=highlight_chrpos
        )
        log.write("  -highlight_windowkb is set to: {} kb".format(highlight_windowkb), verbose=verbose)

    # Pinpoint variants
    if pinpoint is None:
        pinpoint = []
    if len(pinpoint) > 0:
        pinpoint, pinpoint_color = _normalize_group(
            items=pinpoint,
            colors=pinpoint_color,
            item_name="pinpoint variant",
            log=log,
            verbose=verbose,
            is_chrpos_mode=False
        )

    # Region specification
    if region is not None:
        chr_, start_, end_ = region[0], region[1], region[2]
        log.write(" -Region to plot : chr{}:{}-{}.".format(chr_, start_, end_), verbose=verbose)

    return plot_label, build, highlight, highlight_color, pinpoint, pinpoint_color

@normalize_series_inputs(keys=["highlight","pinpoint","anno_set"])
@safefig
def _mqqplot(insumstats,            
          chrom="CHR",
          pos="POS",
          p="P",
          snpid="SNPID",
          eaf="EAF",
          ea="EA",
          nea="NEA",
          check = True,
          chr_dict = None,
          xtick_chr_dict = None,
          vcf_path=None,
          vcf_chr_dict = None,
          gtf_path="default",  # Default path for GTF data
          gtf_chr_dict = None,
          gtf_gene_name=None,
          rr_path="default",  # Default path for recombination rate data
          rr_header_dict=None,
          rr_chr_dict = None,
          rr_lim=(0,100),  # Default recombination rate limits
          rr_ylabel=True,  # Default to showing recombination rate y-label
          mlog10p="MLOG10P",
          scaled=False,
          mode="mqq",
          scatter_kwargs=None,
          qq_scatter_kwargs=None,
          qq_line_color = "grey",
          qq_xlabels = None,
          qq_xlim = None,
          region = None,
          region_title=None,
          region_title_kwargs=None,
          region_ref=None,
          region_ref_second=None,
          region_step = 21,
          region_grid = False,
          region_grid_line = None,
          region_lead_grid = True,
          region_lead_grid_line = None,
          region_hspace=0.02,
          region_ld_threshold = None,
          region_ld_legend = True,
          region_ld_colors = None,
          region_ld_colors_m = None,
          region_recombination = True,
          region_protein_coding = True,
          region_flank_factor = 0.05,
          region_anno_bbox_kwargs = None,
          region_marker_shapes=None,
          region_legend_marker=True,
          region_ref_alias = None,
          ld_path=None,
          ld_map_path=None,
          ld_fmt = "npz",
          ld_if_square =  False,
          ld_if_add_T = False,
          ld_map_rename_dic = None,
          ld_map_kwargs = None,
          cbar_title='LD $\mathregular{r^2}$ with variant',
          cbar_fontsize = None,
          cbar_scale=True,
          cbar_font_family = None,
          cbar_bbox_to_anchor = (0,0,1,1),
          cbar_equal_aspect = True,
          cbar_w_scale=1,
          cbar_h_scale=1,
          cbar_downward_offset =1.3, 
          cbar_borderpad=None,
          track_n=4,
          track_n_offset=0,
          track_fontsize_ratio=0.95,
          track_exon_ratio=1,
          track_text_offset=1,
          track_font_family = None,
          taf = None,
          tabix=None,
          mqqratio=3,
          bwindowsizekb = 100,
          density_color=False,
          density_range=None,
          density_trange=(0,10),
          density_threshold=5,
          density_tpalette="Blues",
          density_palette="Reds",
          windowsizekb=500,
          anno=None,
          anno_set=None,
          anno_alias=None,
          anno_d=None,
          anno_kwargs=None,
          anno_kwargs_single=None,
          anno_style="right",  # Default annotation style
          anno_fixed_arm_length=None,
          anno_source = "ensembl",  # Default annotation source
          anno_gtf_path=None,
          anno_adjust=False,  # Default to no annotation adjustment
          anno_xshift=None,
          anno_max_iter=100,  # Default maximum iterations for annotation adjustment
          anno_max_rows=40,  # Default maximum number of annotation rows to display
          arrow_kwargs=None,
          arm_offset=None,
          arm_scale=1,
          anno_height=1,
          arm_scale_d=None,
          cut=0,  # Default cut value for shrinking variants above threshold
          skip=0,  # Default skip value for -log10(P) < skip
          ystep=0,  # Default step size for y-axis
          ylabels=None,  # Default to no custom y-axis labels
          ytick3=True,  # Default to 3 y-axis ticks
          cutfactor=10,  # Default factor for shrinking cut line
          cut_line_color="#ebebeb",  # Default color for cut line
          cut_log = False,  # Default to linear scale for cut line
          jagged=False,
          jagged_len=0.01,
          jagged_wid=0.01,
          sig_line=True,
          sig_level=5e-8,
          anno_sig_level=5e-8,
          sig_line_color="grey",  # Default color for significance line
          suggestive_sig_line=False,  # Default to no suggestive significance line
          suggestive_sig_level=5e-6,  # Default suggestive significance level
          suggestive_sig_line_color="grey",  # Default color for suggestive significance line
          additional_line = None,  # Default to no additional lines
          additional_line_color = None,  # Default colors for additional lines
          sc_linewidth=2,  # Default line width for significance lines
          highlight = None,
          highlight_chrpos = False,
          highlight_color="#CB132D",  # Default highlight color
          highlight_windowkb = 500,  # Default window size for highlighting (kb)
          highlight_anno_kwargs = None,  # Default to no additional highlight annotation arguments
          highlight_lim = None,  # Default to no custom highlight limits
          highlight_lim_mode = "absolute",  # Default highlight limit mode
          pinpoint= None,
          pinpoint_color ="red",
          stratified=False,
          maf_bins=None,
          maf_bin_colors = None,
          gc=True,
          include_chrXYMT = True,
          ylim=None,
          xpad=None,
          xpadl=None,
          xpadr=None,
          xtight=False,
          chrpad=0.03, 
          drop_chr_start=False,
          title =None,
          mtitle=None,
          mtitle_pad=1.08,
          qtitle=None,
          qtitle_pad=1.08,
          ylabel=None,
          xlabel=None,
          title_pad=1.08, 
          title_fontsize=13,
          title_kwargs=None,
          fontsize = 9,
          font_family=None,
          fontfamily="Arial",
          math_fontfamily="dejavusans",
          anno_fontsize = 9,
          fig_kwargs=None,
          figax=None,
          colors=None,
          marker_size=None,
          use_rank=False,
          verbose=True,
          repel_force=0.03,
          build=None,
          _posdiccul=None,
          dpi=200,
          save=None,
          save_kwargs=None,
          _invert=False,
          _chrom_df_for_i=None,
          _if_quick_qc=True,
          _get_region_lead=False,
          expected_min_mlog10p=0,
          log=Log()
          ):
    """
    Create an MQQ plot with multiple modes (Manhattan, QQ, Regional, Brisbane).
    
    Parameters
    ----------
    insumstats : pandas.DataFrame
        Input summary statistics dataframe containing GWAS results.
    chrom : str, default='CHR'
        Column name for chromosome.
    pos : str, default='POS'
        Column name for genomic position.
    p : str, default='P'
        Column name for p-values.
    snpid : str, default='SNPID'
        Column name for variant identifier. Will auto-detect 'SNPID' or 'rsID' if present.
    eaf : str, default='EAF'
        Column name for effect allele frequency.
    ea : str, default='EA'
        Column name for effect allele.
    nea : str, default='NEA'
        Column name for non-effect allele.
    mlog10p : str, default='MLOG10P'
        Column name for -log10(P) values.
    mode : str, default='mqq'
        Plotting mode. Options:
        - 'm' : Manhattan plot only
        - 'qq' : QQ plot only
        - 'mqq' : Manhattan-QQ layout (default)
        - 'r' : Regional plot (requires region parameter)
        - 'b' : Brisbane/density plot
        - 'mb' : Manhattan-Brisbane layout
    mqqratio : int, default=3
        Ratio for MQQ plot layout (Manhattan to QQ width ratio).
    scaled : bool, default=False
        Whether P-values are already in -log10 scale. Auto-detected if mlog10p column exists.
    check : bool, default=True
        Whether to perform input data quality checks.
    scatter_kwargs : dict, optional
        Arguments for scatter plot styling. Default is None.
    qq_scatter_kwargs : dict, optional
        Arguments for QQ plot scatter styling. Default is None.
    qq_line_color : str, default='grey'
        Color for QQ plot reference line.
    qq_xlabels : list, optional
        Custom x-axis labels for QQ plot. Default is None.
    qq_xlim : tuple or list, optional
        X-axis limits for QQ plot. Default is None.
    colors : list, default=['#597FBD','#74BAD3']
        Color palette for plot. Alternates by chromosome.
    marker_size : tuple, default=(5,20)
        Size range for markers. Passed to sns.scatterplot(). Pass like [20,20] for a fixed size. 
        Overwrites "s" in scatter_kwargs.
    use_rank : bool, default=False
        Whether to use rank for plotting instead of position.
    region : tuple or list, optional
        Genomic region specification (chr, start, end). Required in 'r' mode.
        Region can be determined using `get_region_start_and_end`.
    region_title : str, optional
        Title text for regional plot. Default is None.
    region_title_kwargs : dict, optional
        Regional plot title styling args. Default is {'family': 'Arial', 'size': 12}.
    region_ref : str or list, optional
        Reference variant(s) for LD calculation in regional plot. Default is None.
    region_ref_second : str, optional
        Second reference variant for LD calculation. Default is None.
    region_ref_alias : dict, optional
        Dictionary mapping reference variant IDs to display aliases. Default is None.
    region_step : int, default=21
        Step size for regional plot grid lines.
    region_grid : bool, default=False
        Whether to show grid lines in regional plot.
    region_grid_line : dict, optional
        Styling arguments for region grid lines. Default is {'linewidth': 2, 'linestyle': '--'}.
    region_lead_grid : bool, default=True
        Whether to show grid lines at lead variant position.
    region_lead_grid_line : dict, optional
        Styling arguments for lead variant grid lines. 
        Default is {'alpha': 0.5, 'linewidth': 2, 'linestyle': '--', 'color': '#FF0000'}.
    region_hspace : float, default=0.02
        Horizontal spacing between regional plot panels.
    region_ld_threshold : list, optional
        LD r^2 thresholds list for color stratification. Default is [0.2, 0.4, 0.6, 0.8].
    region_ld_legend : bool, default=True
        Whether to show LD legend in regional plot.
    region_ld_colors : list, optional
        LD colors (list or palette) aligned with thresholds. 
        Default is ['#E4E4E4', '#020080', '#86CEF9', '#24FF02', '#FDA400', '#FF0000', '#FF0000'].
    region_ld_colors_m : list, optional
        LD colors for multi-lead mode. 
        Default is ['#E51819', '#367EB7', 'green', '#F07818', '#AD5691', 'yellow', 'purple'].
    region_recombination : bool, default=True
        Whether to show recombination rate track in regional plot.
    region_protein_coding : bool, default=True
        Whether to show only protein-coding genes in gene track.
    region_flank_factor : float, default=0.05
        Flanking region factor for regional plot (proportion of region size).
    region_anno_bbox_kwargs : dict, optional
        Bounding box arguments for region annotations. Default is {'ec': 'None', 'fc': 'None'}.
    region_marker_shapes : list, optional
        Marker shapes for different LD bins. 
        Default is ['o', '^', 's', 'D', '*', 'P', 'X', 'h', '8'].
    region_legend_marker : bool, default=True
        Whether to show marker legend in regional plot.
    vcf_path : str, optional
        Path to VCF file for LD calculation. Default is None. For regional plot with LD, it is required.
    vcf_chr_dict : dict, optional
        Chromosome name mapping dictionary for VCF file. Default is None.
    ld_path : str, optional
        Path to pre-computed LD file. Default is None.
    ld_map_path : str, optional
        Path to LD map file. Default is None.
    ld_fmt : str, default='npz'
        Format of LD file ('npz', 'txt', etc.). Default is 'npz'.
    ld_if_square : bool, default=False
        Whether LD matrix is square format.
    ld_if_add_T : bool, default=False
        Whether to transpose LD matrix.
    ld_map_rename_dic : dict, optional
        Dictionary for renaming LD map columns. Default is None.
    ld_map_kwargs : dict, optional
        Additional arguments for LD map processing. Default is None.
    tabix : pysam.TabixFile, optional
        Tabix file handle for VCF. Auto-created if vcf_path provided. Default is None.
    gtf_path : str, default='default'
        Path to GTF file for gene annotations. 'default' uses auto-detected path.
    gtf_chr_dict : dict, optional
        Chromosome name mapping dictionary for GTF file. Default is None.
    gtf_gene_name : str, optional
        Column name for gene names in GTF. Default is None.
    rr_path : str, default='default'
        Path to recombination rate file. 'default' uses auto-detected path.
    rr_header_dict : dict, optional
        Header mapping dictionary for recombination rate file. Default is None.
    rr_chr_dict : dict, optional
        Chromosome name mapping dictionary for recombination rate file. Default is None.
    rr_lim : tuple, default=(0,100)
        Limits for recombination rate y-axis (cM/Mb).
    rr_ylabel : bool, default=True
        Whether to show y-axis label for recombination rate track.
    cbar_title : str, default='LD $\\mathregular{r^2}$ with variant'
        Title for LD colorbar.
    cbar_fontsize : int, optional
        Font size for colorbar. Default is None (uses fontsize).
    cbar_font_family : str, optional
        Font family for colorbar. Default is None.
    cbar_scale : bool, default=True
        Whether to scale colorbar.
    cbar_bbox_to_anchor : tuple, default=(0,0,1,1)
        Bounding box anchor for colorbar positioning.
    cbar_equal_aspect : bool, default=True
        Whether colorbar has equal aspect ratio.
    cbar_w_scale : float, default=1
        Width scale factor for colorbar.
    cbar_h_scale : float, default=1
        Height scale factor for colorbar.
    cbar_downward_offset : float, default=1.3
        Downward offset for colorbar positioning.
    cbar_borderpad : float, optional
        Border padding for colorbar. Default is None.
    track_n : int, default=4
        Number of gene tracks to display.
    track_n_offset : int, default=0
        Offset for gene track numbering.
    track_fontsize_ratio : float, default=0.95
        Font size ratio for gene track labels.
    track_exon_ratio : float, default=1
        Ratio for exon display size.
    track_text_offset : int, default=1
        Text offset for gene track labels.
    track_font_family : str, optional
        Font family for gene track. Default is None.
    taf : list, optional
        Track appearance factors [n, n_offset, fontsize_ratio, exon_ratio, text_offset]. 
        Default is None (uses individual parameters).
    bwindowsizekb : int, default=100
        Window size in kb for density calculation.
    density_color : bool, default=False
        Whether to color variants by density.
    density_range : tuple, optional
        Range for density coloring. Default is None.
    density_trange : tuple, default=(0,10)
        Threshold range for density coloring.
    density_threshold : float, default=5
        Density threshold value.
    density_tpalette : str, default='Blues'
        Color palette for density below threshold.
    density_palette : str, default='Reds'
        Color palette for density above threshold.
    anno : {None, bool, str, 'GENENAME'}, default=None
        Specify which data to use for annotation. Options:
        - None: no annotation
        - True: annotate variants with chromosome and position like chr:pos
        - 'GENENAME': annotate variants with closest gene names
        - str: annotate variants with values in column with the header
    anno_set : list, optional
        Set of variant IDs to annotate. Default is None.
    anno_alias : dict, optional
        Dictionary mapping SNP IDs to custom annotation labels. Default is None.
    anno_d : dict, optional
        Dictionary mapping annotation indices to custom positioning options ('l' or 'r').
        For example, {"0":"l"} means adjust the direction of the 1st arm end to left. 
        Default is None.
    anno_kwargs : dict, optional
        Dictionary of default styling arguments for annotations. Default is None.
    anno_kwargs_single : dict, optional
        Dictionary mapping SNP IDs to custom styling arguments. Default is None.
    anno_style : {'right', 'tight', 'expand'}, default='right'
        Style of annotation placement.
    anno_fixed_arm_length : float, optional
        Fixed arm length for annotations. Default is None.
    anno_source : str, default='ensembl'
        Source for gene name annotations.
    anno_gtf_path : str, optional
        Path to GTF file for annotations when anno is 'GENENAME'. 
        If None, auto-detected (preferred). Default is None.
    anno_adjust : bool, default=False
        Whether to automatically adjust text positions to prevent overlap.
    anno_xshift : float, optional
        X-axis shift for annotations. Default is None.
    anno_max_iter : int, default=100
        Maximum iterations for text repulsion algorithm.
    anno_max_rows : int, default=40
        Maximum number of annotation rows to display. If more variants are provided,
        they will be sorted by p-value or -log10(p-value) and only the top ones will be shown.
    arrow_kwargs : dict, optional
        Arguments for annotation arrows. Default is None.
    arm_offset : float, optional
        Offset for annotation arms. Default is None.
    arm_scale : float, default=1
        Scaling factor for annotation arm length.
    arm_scale_d : dict, optional
        Dictionary mapping annotation indices to custom arm scaling factors. Default is None.
    anno_height : float, default=1
        Height for annotations.
    repel_force : float, default=0.03
        Force for repelling overlapping labels.
    sig_line : bool, default=True
        Whether to show significance line.
    sig_level : float, default=5e-8
        Significance level threshold to plot a line on the Manhattan plot and for marker sizing.
        Controls the reference line position and marker sizes for variants above this threshold.
    anno_sig_level : float, default=5e-8
        Significance level for extracting lead variants to annotate on the plot.
        Can be set independently from sig_level to control annotation separately from visualization.
        When too many lead variants, it is helpful to set anno_sig_level to a more strict value.
    sig_line_color : str, default='grey'
        Color for significance line.
    suggestive_sig_line : bool, default=False
        Whether to show suggestive significance line.
    suggestive_sig_level : float, default=5e-6
        Suggestive significance level.
    suggestive_sig_line_color : str, default='grey'
        Color for suggestive significance line.
    additional_line : list, optional
        Additional significance lines to plot (as p-values). Default is None.
    additional_line_color : list, optional
        Colors for additional lines. Default is None.
    sc_linewidth : int, default=2
        Line width for significance lines.
    highlight : list, optional
        A list of variant identifiers (e.g., SNPIDs). Each specified variant will be treated 
        as a focal point, and all variants in surrounding loci will be highlighted in distinct colors.
    highlight_chrpos : bool, default=False
        Whether to interpret highlight as chromosome:position format instead of SNP IDs.
    highlight_color : str or list, default='#CB132D'
        Color(s) for highlighted variants. Can be a single color or list of colors.
    highlight_windowkb : int, default=500
        Window size for highlighting variants in kb.
    highlight_anno_kwargs : dict, optional
        Arguments for highlight annotations. Default is None.
    highlight_lim : list, optional
        Custom limits for highlighting. Default is None.
    highlight_lim_mode : {'absolute', 'relative'}, default='absolute'
        Mode for highlight limits.
    pinpoint : list, optional
        List of variants to pinpoint (mark with special styling). Default is None.
    pinpoint_color : str or list, default='red'
        Color(s) for pinpointed variants. Can be a single color or list of colors.
    cut : float, default=0
        Threshold for shrinking extremely large -log10(P) values. Variants with MLOG10P 
        greater than this threshold will be shrunk by a factor of `cutfactor`.
        Useful when the top signal is very strong and dominates the plot.
    skip : float, default=0
        Minimum -log10(P) value required for a variant to be plotted. Variants with 
        MLOG10P < skip will be omitted. This is helpful for fast plotting or focusing on stronger signals.
    cutfactor : int, default=10
        Factor for shrinking cut line.
    cut_line_color : str, default='#ebebeb'
        Color for cut line.
    cut_log : bool, default=False
        Whether to use log scale for cut line.
    ystep : float, default=0
        Step size for y-axis ticks.
    ylabels : list, optional
        Custom y-axis tick labels. Default is None.
    ytick3 : bool, default=True
        Whether to use 3 y-axis ticks.
    ylim : tuple, optional
        Y-axis limits for plot. Default is None.
    jagged : bool, default=False
        Whether to make y-axis jagged (break axis for cut values).
    jagged_len : float, default=0.01
        Length of jagged line segment.
    jagged_wid : float, default=0.01
        Width of jagged line segment.
    xpad : float, optional
        X-axis padding proportion (applies to both sides). Default is None.
    xpadl : float, optional
        Left X-axis padding proportion. Default is None.
    xpadr : float, optional
        Right X-axis padding proportion. Default is None.
    xtight : bool, default=False
        Whether to use tight X-axis padding.
    chrpad : float, default=0.03
        Chromosome padding factor (spacing between chromosomes).
    drop_chr_start : bool, default=False
        Whether to drop chromosome start position from x-axis.
    xlabel : str, optional
        X-axis label. Default is None (auto-generated).
    title : str, optional
        Overall plot title. Only used when both Manhattan and QQ plot are created. Default is None.
    mtitle : str, optional
        Title for Manhattan plot. Default is None.
    mtitle_pad : float, default=1.08
        Padding for Manhattan title.
    qtitle : str, optional
        Title for QQ plot. Default is None.
    qtitle_pad : float, default=1.08
        Padding for QQ plot title.
    title_pad : float, default=1.08
        Padding for overall plot title.
    title_fontsize : int, default=13
        Font size for titles.
    title_kwargs : dict, optional
        Additional styling arguments for titles. Default is None.
    ylabel : str, optional
        Y-axis label. Default is "$\\mathregular{-log_{10}(P)}$".
        In 'b' mode: default is "Density of GWAS \\n SNPs within "+str(bwindowsizekb)+" kb".
    fontsize : int, default=9
        Base font size for plot elements.
    font_family : str, optional
        Font family for text elements. Default is None.
    fontfamily : str, default='Arial'
        Alternative name for font family (used if font_family is None).
    math_fontfamily : str, default='dejavusans'
        Font family for math text rendering.
    anno_fontsize : int, default=9
        Font size for annotations.
    fig_kwargs : dict, optional
        Figure arguments for subplots (e.g., figsize). Default is None.
    figax : dict, optional
        Existing figure and axes for plot. Default is None.
    stratified : bool, default=False
        Whether to create stratified QQ plots by MAF.
    maf_bins : list of lists, optional
        Bins for MAF stratification [[lower1, upper1], ...]. 
        Default is [[0, 0.01], [0.01, 0.05], [0.05, 0.25], [0.25, 0.5]].
    maf_bin_colors : list, optional
        Colors for MAF bins. Default is None.
    gc : bool, default=True
        Whether to calculate and display genomic control lambda.
    include_chrXYMT : bool, default=True
        Whether to include sex chromosomes (X, Y) and mitochondrial chromosome (MT) in QQ plot.
    expected_min_mlog10p : float, default=0
        Expected minimum -log10(P) value for QQ plot.
    windowsizekb : int, default=500
        Window size in kb for obtaining significant variants (lead variant detection).
    chr_dict : dict, optional
        Mapping dictionary for chromosome names. Only used when chromosomes in reference 
        files are not standardized. Default is None.
    xtick_chr_dict : dict, optional
        Chromosome number to name mapping for x-axis ticks. Only used when chromosomes 
        in sumstats are not standardized. Default is None.
    build : str, optional
        Genomic build version (e.g., '19', '38'). Default is None (auto-detected).
    dpi : int, default=200
        Dots per inch for figure resolution.
    save : str, optional
        File path to save plot. Default is None.
    save_kwargs : dict, default={'dpi':600, 'transparent':True}
        Arguments for saving the plot.
    verbose : bool, default=True
        Whether to show progress messages.
    log : gwaslab.Log, default=Log()
        Logging object for operation records.
    _posdiccul : dict, optional
        Internal parameter for position dictionary. Default is None.
    _chrom_df_for_i : pandas.DataFrame, optional
        Internal parameter for chromosome dataframe. Default is None.
    _if_quick_qc : bool, default=True
        Internal parameter for quick QC flag. Default is True.
    _get_region_lead : bool, default=False
        Whether to return lead SNP information in regional mode. Default is False.
    _invert : bool, default=False
        Internal parameter for axis inversion. Default is False.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib figure object.
    log : gwaslab.Log
        Updated logging object with operation records.
    lead_snp_is, lead_snp_is_color : tuple, optional
        Lead SNP positions and colors returned only when _get_region_lead=True in regional ('r') mode.
    
    Notes
    -----
    - The function supports multiple plotting modes: Manhattan plots, QQ plots, regional plots, 
      and density (Brisbane) plots.
    - Regional plots require LD data (VCF or pre-computed LD files) and a region specification.
    - Annotation can be done using gene names, custom columns, or variant IDs.
    - The function automatically handles chromosome name standardization and coordinate conversion.
    - For large datasets, use `skip` parameter to filter low-significance variants for faster plotting.
    """
    # Extract dataframe if Sumstats object is passed
    # Store reference to original object to check meta
    sumstats_obj = None
    if hasattr(insumstats, 'data') and not isinstance(insumstats, pd.DataFrame):
        sumstats_obj = insumstats
        insumstats = insumstats.data
    
    # scan all local args for series and convert to list for further processing
    for arg in locals().values():
        if isinstance(arg, pd.Series):
            arg = arg.tolist()

    if snpid in insumstats.columns:
        snpid=snpid
    elif "SNPID" in insumstats.columns:
        snpid="SNPID"
    elif "rsID" in insumstats.columns:
        snpid="rsID"
    
    if "EAF" not in insumstats.columns:
        eaf=None
    
    chr_dict = _update_kwargs(chr_dict, get_chr_to_number())
    xtick_chr_dict = _update_kwargs(xtick_chr_dict, get_number_to_chr())
    gtf_chr_dict = _update_kwargs(gtf_chr_dict, get_number_to_chr())
    rr_chr_dict = _update_kwargs(rr_chr_dict, get_number_to_chr())

    style = set_plot_style(
        plot="plot_mqq",
        mode=mode,
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs,
        qq_scatter_kwargs=qq_scatter_kwargs,
        title_kwargs=title_kwargs,
        fontsize=fontsize,
        fontfamily=fontfamily,
        colors=colors,
        marker_size=marker_size,
        dpi=dpi,
        verbose=verbose,
        log=log,
    )
    fig_kwargs =style["fig_kwargs"]
    save_kwargs =style["save_kwargs"]
    scatter_kwargs =style["scatter_kwargs"]
    qq_scatter_kwargs = style["qq_scatter_kwargs"]
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    colors = style["colors"]
    marker_size = style["marker_size"]
    dpi = style["dpi"]
    title_kwargs = style["title_kwargs"]

    anno_set = _update_arg(anno_set, [])
    # Auto-detect scaled input when -log10(P) column is present with valid values
    try:
        if (scaled is False) and (mlog10p in insumstats.columns):
            if _mlog10p_series.notna().sum() > 0:
                scaled = True
                log.write(f"Auto-detected scaled input using column '{mlog10p}'.", verbose=verbose)
    except Exception:
        pass

    taf = _update_arg(taf, [track_n,track_n_offset,track_fontsize_ratio,track_exon_ratio,track_text_offset])
    font_family = _update_arg(font_family, fontfamily)

    
    # Step 1: Set arguments and normalize inputs
    if region is not None:
        if region[1] == region[2]:
            raise ValueError("Region should be [chr, start, end] with start < end")
    
    # make region_ref a list of ref variants
    if pd.api.types.is_list_like(region_ref):
        if len(region_ref) == 0 :
            region_ref.append(None)
            if region_ref_second is not None:
                region_ref.append(region_ref_second)
    else:
        region_ref = [region_ref]
        if region_ref_second is not None:
            region_ref.append(region_ref_second)
    region_ref_index_dic = {value: index for index,value in enumerate(region_ref)}

    # scatter_kwargs and qq_scatter_kwargs merged via set_plot_style
        
    # Check if QC is finished in Sumstats meta, skip quick QC if so
    if sumstats_obj is not None:
        qc_performed = sumstats_obj.meta.get("gwaslab", {}).get("basic_check", {}).get("performed", False)
        if qc_performed:
            _if_quick_qc = False
            log.write(" -QC is finished according to meta. Quick QC will be skipped.", verbose=verbose)
        elif check==True and _if_quick_qc==True:
            _if_quick_qc = True
        else:
            _if_quick_qc = False
    else:
        # If not a Sumstats object, use original logic
        if check==True and _if_quick_qc==True:
            _if_quick_qc = True
        else:
            _if_quick_qc = False

    # configure dpi if saving the plot
    fig_kwargs, scatter_kwargs, qq_scatter_kwargs, save_kwargs = _configure_fig_save_kwargs(mode=mode,
                                                                                    save = save, 
                                                                                    fig_kwargs = fig_kwargs, 
                                                                                    scatter_kwargs = scatter_kwargs, 
                                                                                    qq_scatter_kwargs = qq_scatter_kwargs, 
                                                                                    save_kwargs = save_kwargs,
                                                                                    log=log,
                                                                                    verbose=verbose)


    # arm_offset default calculation moved into annotate_single

    plot_label, build, highlight, highlight_color, pinpoint, pinpoint_color = _setup_and_log_mqq_info(
        insumstats=insumstats,
        mode=mode,
        sig_level=sig_level,
        anno_set=anno_set,
        highlight=highlight,
        highlight_chrpos=highlight_chrpos,
        highlight_color=highlight_color,
        highlight_windowkb=highlight_windowkb,
        pinpoint=pinpoint,
        pinpoint_color=pinpoint_color,
        region=region,
        build=build,
        log=log,
        verbose=verbose,
    )
    
    # Build significance threshold series and convert to -log10 scale
    # Add anno_sig_level as last item so it gets transformed by _cut, then remove it after
    if additional_line is None:
        lines_to_plot = pd.Series([sig_level, suggestive_sig_level, anno_sig_level] )
    else:
        lines_to_plot = pd.Series([sig_level, suggestive_sig_level, anno_sig_level] + additional_line ) 
        if additional_line_color is None:
            additional_line_color = ["grey"]
    lines_to_plot = -np.log10(lines_to_plot)

# Step 3: Select layout ####################################################################
    # ax1 : manhattanplot / brisbane plot
    # ax2 : qq plot 
    # ax3 : gene track
    # ax4 : recombination rate
    # cbar : color bar
    # ax5 : miami plot lower panel

    # "m" : Manhattan plot
    # "qq": QQ plot
    # "r" : regional plot
    
    fig, ax1, ax2, ax3, ax4, cbar = _process_layout(mode=mode, 
                                         figax=figax, 
                                         fig_kwargs=fig_kwargs, 
                                         mqqratio=mqqratio, 
                                         region_hspace=region_hspace)
    highlight_i = []
    
# mode specific settings ####################################################################
    if mode=="b":
        sig_level=1,
        sig_line=False,
        #windowsizekb = 100000000   
        mode="mb"
        scatter_kwargs={"marker":"s"}
        marker_size= (marker_size[1],marker_size[1])

# Step 4: Load sumstats #################################################################################################

    usecols = _configure_cols_to_use(insumstats=insumstats, 
                                     snpid=snpid,  
                                     chrom=chrom, 
                                     pos=pos, 
                                     ea=ea, 
                                     nea=nea, 
                                     eaf=eaf, 
                                     p=p, 
                                     mlog10p=mlog10p,
                                     scaled=scaled, 
                                     mode=mode,
                                     stratified=stratified,
                                     anno=anno, 
                                     anno_set=anno_set, 
                                     anno_alias=anno_alias,
                                     _chrom_df_for_i=_chrom_df_for_i,
                                     highlight=highlight,
                                     pinpoint=pinpoint,
                                     density_color=density_color)
    
    sumstats = insumstats[usecols].copy()
    
    #################################################################################################
    
    # Step 5: Standardize and QC
    #Standardize
    ## Annotation
    if (anno == "GENENAME"):
        anno_sig=True
    elif isinstance(anno,str):
        sumstats["Annotation"]=sumstats[anno].astype("string")   
      
    ## P value
    ## m, qq, r
    if "b" not in mode:   
        if scaled is True:
            sumstats["raw_P"] = pd.to_numeric(sumstats[mlog10p], errors='coerce')
        else:
            sumstats["raw_P"] = sumstats[p].astype("float64")
    
    ## CHR & POS
    ## m, qq, b
    if "m" in mode or "r" in mode or "b" in mode: 
        # convert CHR to int
        ## CHR X,Y,MT conversion ############################
        sumstats[pos] = _quick_fix_pos(sumstats[pos], verbose=verbose, log=log)
        sumstats[chrom] = _quick_fix_chr(sumstats[chrom], chr_dict=chr_dict, verbose=verbose, log=log)

    ## r
    
    if region is not None:
        region = _normalize_region(region, chr_dict=chr_dict, sumstats=sumstats, snpid=snpid, rsid="rsID", chrom_col=chrom, pos_col=pos, ea=ea, nea=nea, log=log, verbose=verbose)
        sumstats = _filter_region(sumstats, region, log=log, verbose=verbose)
  
        if len(sumstats)==0:
            log.warning("No valid data! Please check the input.")
            return None
    
    ## EAF
    eaf_raw = pd.Series(dtype="float64")
    if stratified is True: 
        sumstats["MAF"] = _quick_fix_eaf(sumstats[eaf], verbose=verbose, log=log)
        # for stratified qq plot
        eaf_raw = sumstats["MAF"].copy()
        
    if len(highlight)>0 and ("m" in mode):
        sumstats["HUE"] = pd.NA
        sumstats["HUE"] = sumstats["HUE"].astype("Int64")

    log.write("Finished loading specified columns from the statistics",verbose=verbose)


#sanity check############################################################################################################
    log.write("Start data conversion and sanity check:",verbose=verbose)
    
    if _if_quick_qc == False:
        log.write(" -Sanity check will be skipped.", verbose=verbose)
    else:
        sumstats = _sanity_check(sumstats=sumstats, 
                                 mode=mode,
                                 chrom =chrom, 
                                 pos=pos, 
                                 stratified=stratified, 
                                 _if_quick_qc=_if_quick_qc, 
                                 log=log, 
                                 verbose=verbose)
            
    # Step 6: Transform data (highlight/density/p/mlog10p)
    ## configure highlight regions
    if len(highlight)>0 and ("m" in mode):
        # add HUE
        sumstats = _process_highlight(sumstats=sumstats, 
                                                    highlight=highlight, 
                                                    highlight_chrpos=highlight_chrpos, 
                                                    highlight_windowkb=highlight_windowkb, 
                                                    highlight_lim = highlight_lim, 
                                                    highlight_lim_mode = highlight_lim_mode,
                                                    snpid=snpid, 
                                                    chrom=chrom, 
                                                    pos=pos)

# Density mode setup #####################################################################################################              
    if "b" in mode:
        from gwaslab.viz.viz_plot_density_mode import b_mode_setup
        sumstats, bmean, bmedian = b_mode_setup(
            sumstats=sumstats,
            mode=mode,
            bwindowsizekb=bwindowsizekb,
            chrom=chrom,
            pos=pos,
            log=log,
            verbose=verbose,
        )
        lines_to_plot = pd.Series(lines_to_plot.to_list() + [bmean, bmedian])
    
    else:
        bmean, bmedian=0,0 
    
# P value conversion #####################################################################################################  
    
    # add raw_P and scaled_P
    sumstats =  _process_p_value(sumstats=sumstats, 
                                 mode=mode,
                                 p=p, 
                                 mlog10p=mlog10p, 
                                 scaled=scaled, 
                                 log=log, 
                                 verbose=verbose )
    
    # raw -log10(P) per chromosome for QQ panel
    p_toplot_raw = sumstats[["CHR","scaled_P"]].copy()
    
    # filter out variants with -log10p < skip
    sumstats = sumstats.loc[sumstats["scaled_P"]>=skip,:]
    garbage_collect.collect()
    
    # Shrink variants above cut line #########################################################################################
    # Initialize transformed_anno_sig_level (will be updated after _cut if successful)
    transformed_anno_sig_level = -np.log10(anno_sig_level)
    try:
        sumstats["scaled_P"], maxy, maxticker, cut, cutfactor,ylabels_converted, lines_to_plot = _cut(series = sumstats["scaled_P"], 
                                                                        mode =mode, 
                                                                        cut=cut,
                                                                        skip=skip,
                                                                        cutfactor = cutfactor,
                                                                        ylabels=ylabels,
                                                                        cut_log = cut_log,
                                                                        verbose =verbose, 
                                                                        lines_to_plot=lines_to_plot,
                                                                        log = log
                                                                        )
        # Extract transformed anno_sig_level from lines_to_plot (it's at index 2)
        # Then remove it so it doesn't get drawn as a line
        transformed_anno_sig_level = lines_to_plot.iloc[2]
        # Remove anno_sig_level (index 2), keep sig_level (0), suggestive_sig_level (1), and additional_line (3+)
        if additional_line is None:
            lines_to_plot = lines_to_plot.iloc[:2]
        else:
            lines_to_plot = pd.Series(list(lines_to_plot.iloc[:2]) + list(lines_to_plot.iloc[3:])).reset_index(drop=True)
    except:
        log.warning("No valid data! Please check the input.")
        return None
    
    log.write("Finished data conversion and sanity check.",verbose=verbose)
    
    # default: no annotations when not creating Manhattan/Regional panel
    
    
    # Step 7: Assign x-axis index (i) and create panels ##########################################################################################################
    log.write("Start to create {} with ".format(plot_label)+str(len(sumstats))+" variants...",verbose=verbose)
    ## regional data sources
    if vcf_path is not None:
        if tabix is None:
            vcf_chr_dict, tabix = prepare_vcf_context(vcf_path=vcf_path, vcf_chr_dict=vcf_chr_dict, log=log, verbose=verbose)
        sumstats = process_vcf(sumstats=sumstats, 
                               vcf_path=vcf_path,
                               region=region, 
                               region_ref=region_ref, 
                               log=log ,
                               pos=pos,
                               ea=ea,
                               nea=nea,
                               region_ld_threshold=region_ld_threshold,
                               verbose=verbose,
                               vcf_chr_dict=vcf_chr_dict,
                               tabix=tabix)
    elif ld_path is not None:
        sumstats = process_ld(
            ld_path=ld_path,
            ld_map_path = ld_map_path, 
            sumstats=sumstats, 
            region=region, 
            region_ref=region_ref, 
            log=log ,
            pos=pos,
            ea=ea,
            nea=nea,
            region_ld_threshold=region_ld_threshold,
            verbose=verbose,
            ld_fmt = ld_fmt,
            ld_if_square =  ld_if_square,
            ld_if_add_T = ld_if_add_T,
            ld_map_rename_dic = ld_map_rename_dic,
            ld_map_kwargs = ld_map_kwargs
        )
    #sort & add id
    ## Manhatann plot ###################################################
    if ("m" in mode) or ("r" in mode): 
        # Assign index i and tick positions for x-axis
        if _chrom_df_for_i is None:
            sumstats,chrom_df=_quick_assign_i_with_rank(sumstats, chrpad=chrpad, use_rank=use_rank, chrom="CHR",pos="POS",drop_chr_start=drop_chr_start,_posdiccul=_posdiccul)
        else:
            chrom_df = _chrom_df_for_i
        ## Assign marker size ##############################################
        sumstats["s"]=1
        if "b" not in mode:
            sumstats.loc[sumstats["scaled_P"]>-np.log10(5e-4),"s"]=2
            sumstats.loc[sumstats["scaled_P"]>-np.log10(suggestive_sig_level),"s"]=3
            sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level),"s"]=4
        sumstats["chr_hue"]=sumstats[chrom].astype("string")

        if "r" in mode:
            if vcf_path is None and ld_path is None:
                sumstats["LD"]=100
                sumstats["SHAPE"]=1
            sumstats["chr_hue"]=sumstats["LD"]
                
        ## default seetings
        # assign to_plot for scatter plot
        to_plot = None
        palette = sns.color_palette(colors,n_colors=sumstats[chrom].nunique())  

        legend = None
        style=None
        linewidth=0
        edgecolor="black"
        # if regional plot assign colors
        if "r" in mode:
            from gwaslab.viz.viz_plot_regional2 import regional_mode_setup
            legend, linewidth, palette, style, markers, to_plot, edgecolor = regional_mode_setup(
                sumstats=sumstats,
                region_ref=region_ref,
                region_ld_colors=region_ld_colors,
                region_marker_shapes=region_marker_shapes,
                region_ld_colors_m=region_ld_colors_m,
                region_ld_threshold=region_ld_threshold,
                chrom=chrom,
                pos=pos,
                vcf_path=vcf_path,
                ld_path=ld_path,
                log=log,
                verbose=verbose,
            )
            scatter_kwargs["markers"] = markers

        # Step 8: Create Manhattan panel
        highlight_i = draw_manhattan_panel(
            ax1=ax1,
            sumstats=sumstats,
            snpid=snpid,
            palette=palette,
            marker_size=marker_size,
            style=style,
            linewidth=linewidth,
            edgecolor=edgecolor,
            legend=legend,
            scatter_kwargs=scatter_kwargs,
            highlight=highlight,
            highlight_chrpos=highlight_chrpos,
            highlight_color=highlight_color,
            density_color=density_color,
            density_range=density_range,
            density_trange=density_trange,
            density_threshold=density_threshold,
            density_tpalette=density_tpalette,
            density_palette=density_palette,
            pinpoint=pinpoint,
            pinpoint_color=pinpoint_color,
            chrom=chrom,
            log=log,
            verbose=verbose,
        )
            
        # Step 9: Add region panel if specified ##################################################
        if (region is not None) and ("r" in mode):
            
            ax1, ax3, ax4, cbar, lead_snp_is, lead_snp_is_color =_plot_regional(
                                sumstats=sumstats,
                                fig=fig,
                                ax1=ax1,
                                ax3=ax3,
                                region=region,
                                vcf_path=vcf_path,
                                ld_path=ld_path,
                                marker_size=marker_size,
                                build=build,
                                cbar_scale=cbar_scale,
                                cbar_fontsize=cbar_fontsize,
                                cbar_bbox_to_anchor=cbar_bbox_to_anchor,
                                cbar_w_scale=cbar_w_scale,
                                cbar_h_scale=cbar_h_scale,
                                cbar_equal_aspect=cbar_equal_aspect,
                                cbar_downward_offset =cbar_downward_offset, 
                                cbar_borderpad=cbar_borderpad,
                                cut_line_color=cut_line_color,
                                gtf_path=gtf_path,
                                gtf_chr_dict = gtf_chr_dict,
                                gtf_gene_name=gtf_gene_name,
                                rr_path=rr_path,
                                rr_header_dict=rr_header_dict,
                                rr_chr_dict = rr_chr_dict,
                                rr_lim=rr_lim,
                                rr_ylabel=rr_ylabel,
                                mode=mode,
                                region_step = region_step,
                                region_ref = region_ref,
                                region_ref_index_dic = region_ref_index_dic,
                                region_ref_alias = region_ref_alias,
                                region_grid = region_grid,
                                region_grid_line = region_grid_line,
                                region_lead_grid = region_lead_grid,
                                region_lead_grid_line = region_lead_grid_line,
                                region_title=region_title,
                                region_title_kwargs=region_title_kwargs,
                                region_ld_legend = region_ld_legend,
                                region_legend_marker=region_legend_marker,
                                region_ld_threshold = region_ld_threshold,
                                palette = palette,
                                region_marker_shapes = region_marker_shapes,
                                region_recombination = region_recombination,
                                region_protein_coding=region_protein_coding,
                                region_flank_factor =region_flank_factor,
                                track_font_family=track_font_family,
                                taf=taf,
                                pos=pos,
                                verbose=verbose,
                                log=log
                            )
            
        else:
            lead_snp_is =[]
            lead_snp_is_color = []
        
        log.write("Finished creating {} successfully".format(plot_label),verbose=verbose)
        
        # Use transformed anno_sig_level threshold (after _cut transformation) for annotation
        if "b" in mode:
            from gwaslab.viz.viz_plot_density_mode import b_scaled_threshold
            # For bubble mode, use original threshold since density uses different scaling
            scaled_threhosld = b_scaled_threshold(anno_sig_level)
        else:
            # Use the transformed threshold that matches the scaled_P scale after _cut
            scaled_threhosld = transformed_anno_sig_level
            anno_sig_level = np.power(10, -scaled_threhosld)
        # Get top variants for annotation #######################################################
        from gwaslab.viz.viz_plot_manhattan_like import _extract_to_annotate
        to_annotate = pd.DataFrame()
        to_annotate = _extract_to_annotate(
            sumstats=sumstats,
            snpid=snpid,
            chrom=chrom,
            pos=pos,
            scaled_threhosld=scaled_threhosld,
            anno_sig_level=anno_sig_level,
            windowsizekb=windowsizekb,
            scaled=scaled,
            anno=anno,
            anno_set=anno_set,
            anno_source=anno_source,
            anno_gtf_path=anno_gtf_path,
            build=build,
            mode=mode,
            log=log,
            verbose=verbose,
        )

        # Manhatann-like plot Finished #####################################################################

    # Step 10: Process labels and figure arts
    if ("m" in mode) or ("r" in mode):
        ax1, ax3 = _process_xtick(
            ax1=ax1,
            ax3=ax3,
            mode=mode,
            chrom_df=chrom_df,
            xtick_chr_dict=xtick_chr_dict,
            fontsize=fontsize,
            font_family=font_family,
            log=log,
            verbose=verbose,
        )

        ax1, ax3 = _process_xlabel(
            region=region,
            xlabel=xlabel,
            ax1=ax1,
            gtf_path=gtf_path,
            mode=mode,
            fontsize=fontsize,
            font_family=font_family,
            ax3=ax3,
            log=log,
            verbose=verbose,
        )

        ax1, ax4 = _process_ylabel(
            ylabel=ylabel,
            ax1=ax1,
            mode=mode,
            bwindowsizekb=bwindowsizekb,
            fontsize=fontsize,
            font_family=font_family,
            math_fontfamily=math_fontfamily,
            ax4=ax4,
            log=log,
            verbose=verbose,
        )

        ax1 = _set_yticklabels(
            cut=cut,
            cutfactor=cutfactor,
            cut_log=cut_log,
            ax1=ax1,
            skip=skip,
            maxy=maxy,
            maxticker=maxticker,
            ystep=ystep,
            sc_linewidth=sc_linewidth,
            cut_line_color=cut_line_color,
            fontsize=fontsize,
            font_family=font_family,
            ytick3=ytick3,
            ylabels=ylabels,
            ylabels_converted=ylabels_converted,
            log=log,
            verbose=verbose,
        )

        ax1, ax4 = _process_ytick(
            ax1=ax1,
            fontsize=fontsize,
            font_family=font_family,
            ax4=ax4,
            log=log,
            verbose=verbose,
        )

        if cbar is not None:
            cbar = _process_cbar(
                cbar,
                cbar_fontsize=cbar_fontsize,
                cbar_font_family=cbar_font_family,
                cbar_title=cbar_title,
                log=log,
                verbose=verbose,
            )

        ax1 = _process_spine(ax1, mode)

        ax1 = _process_line(
            ax1,
            sig_line,
            suggestive_sig_line,
            additional_line,
            lines_to_plot,
            sc_linewidth,
            sig_line_color,
            suggestive_sig_line_color,
            additional_line_color,
            mode,
            bmean,
            bmedian,
            log=log,
            verbose=verbose,
        )

        if mtitle and anno and to_annotate is not None and len(to_annotate) > 0:
            pad = (ax1.transData.transform((skip, mtitle_pad * maxy))[1] - ax1.transData.transform((skip, maxy)))[1]
            _ax_title_kwargs = dict(title_kwargs) if title_kwargs is not None else {}
            if "family" not in _ax_title_kwargs:
                _ax_title_kwargs["family"] = font_family
            if "fontsize" not in _ax_title_kwargs:
                _ax_title_kwargs["fontsize"] = title_fontsize
            ax1.set_title(mtitle, pad=pad, **_ax_title_kwargs)
        elif mtitle:
            _ax_title_kwargs = dict(title_kwargs) if title_kwargs is not None else {}
            if "family" not in _ax_title_kwargs:
                _ax_title_kwargs["family"] = font_family
            if "fontsize" not in _ax_title_kwargs:
                _ax_title_kwargs["fontsize"] = title_fontsize
            ax1.set_title(mtitle, **_ax_title_kwargs)

        # Y axis jagged
        if jagged==True:
            ax1 = _jagged_y(cut=cut,skip=skip,ax1=ax1,mode=1,mqqratio=mqqratio,jagged_len=jagged_len,jagged_wid=jagged_wid,log=log, verbose=verbose)
            if "qq" in mode:
                ax2 = _jagged_y(cut=cut,skip=skip,ax1=ax2,mode=2,mqqratio=mqqratio,jagged_len=jagged_len,jagged_wid=jagged_wid,log=log, verbose=verbose)
        
        # XY lim
        if ylim is not None and  ylim != []:
            ax1.set_ylim(ylim)
            if "qq" in mode:
                ax2.set_ylim(ylim)
        
        ax1 = _add_pad_to_x_axis(ax1, xpad, xpadl, xpadr, sumstats, pos, chrpad, xtight, log = log, verbose=verbose)

    # Titles 
    has_annotations = ("m" in mode or "r" in mode) and anno and to_annotate is not None and len(to_annotate) > 0
    if title and has_annotations:
        _fig_title_kwargs = dict(title_kwargs) if title_kwargs is not None else {}
        if "family" not in _fig_title_kwargs:
            _fig_title_kwargs["family"] = font_family
        if "fontsize" not in _fig_title_kwargs:
            _fig_title_kwargs["fontsize"] = title_fontsize
        fig.suptitle(title, x=0.5, y=title_pad, **_fig_title_kwargs)
    else:
        title_pad = title_pad -0.05
        _fig_title_kwargs = dict(title_kwargs) if title_kwargs is not None else {}
        if "family" not in _fig_title_kwargs:
            _fig_title_kwargs["family"] = font_family
        if "fontsize" not in _fig_title_kwargs:
            _fig_title_kwargs["fontsize"] = title_fontsize
        fig.suptitle(title, x=0.5, y=title_pad, **_fig_title_kwargs)

    # Step 11: Add QQ panel if specified #########################################################################################################
    if "qq" in mode:
        # ax2 qqplot
        ax2 =_plot_qq(
                    sumstats=sumstats,
                    p_toplot_raw=p_toplot_raw,
                    ax2=ax2,
                    maxticker=maxticker,
                    marker_size=marker_size,
                    gc=gc,
                    cut=cut,
                    cutfactor=cutfactor,
                    cut_log=cut_log,
                    skip=skip,
                    maxy=maxy,
                    ystep=ystep,
                    colors=colors,
                    qq_line_color=qq_line_color,
                    stratified=stratified,
                    eaf_raw=eaf_raw,
                    maf_bins=maf_bins,
                    maf_bin_colors=maf_bin_colors,
                    fontsize=fontsize,
                    font_family=font_family,
                    qtitle=qtitle,
                    qtitle_pad=qtitle_pad,
                    title_fontsize=title_fontsize,
                    include_chrXYMT=include_chrXYMT,
                    cut_line_color=cut_line_color,
                    linewidth=sc_linewidth,
                    ytick3 = ytick3,
                    ylabels = ylabels,
                    xlabels = qq_xlabels,
                    xlim = qq_xlim,
                    ylabels_converted = ylabels_converted,
                    verbose=verbose,
                    qq_scatter_kwargs=qq_scatter_kwargs,
                    expected_min_mlog10p=expected_min_mlog10p,
                    log=log
                )
    
    # Step 12: Annotation
    if (anno is not None):
        ax1 = annotate_single(
            sumstats=sumstats,
            anno=anno,
            mode=mode,
            ax1=ax1,
            highlight_i=highlight_i,
            highlight_chrpos=highlight_chrpos,
            highlight_anno_kwargs=highlight_anno_kwargs,
            to_annotate=to_annotate,
            anno_d=anno_d,
            anno_alias=anno_alias,
            anno_style=anno_style,
            anno_kwargs=anno_kwargs,
            anno_kwargs_single=anno_kwargs_single,
            arm_scale=arm_scale,
            anno_max_iter=anno_max_iter,
            anno_max_rows=anno_max_rows,
            arm_scale_d=arm_scale_d,
            arm_offset=arm_offset,
            anno_adjust=anno_adjust,
            anno_xshift=anno_xshift,
            anno_fixed_arm_length=anno_fixed_arm_length,
            maxy=maxy,
            anno_fontsize=anno_fontsize,
            font_family=font_family,
            region=region,
            region_anno_bbox_kwargs=region_anno_bbox_kwargs,
            skip=skip,
            anno_height=anno_height,
            arrow_kwargs=arrow_kwargs,
            snpid=snpid,
            chrom=chrom,
            pos=pos,
            repel_force=repel_force,
            verbose=verbose,
            log=log,
            _invert=_invert,
        )

    
    # Step 13: Save figure
    save_figure(fig = fig, save = save, keyword=mode, save_kwargs=save_kwargs, log = log, verbose=verbose)

    garbage_collect.collect()
    # Return matplotlib figure object #######################################################################################
    if _get_region_lead==True:
        return fig, log, lead_snp_is, lead_snp_is_color
    
    log.write("Finished creating plot successfully",verbose=verbose)
    return fig, log

##############################################################################################################################################################################
from gwaslab.viz.viz_plot_manhattan_like import _configure_fig_save_kwargs as _configure_fig_save_kwargs
from gwaslab.viz.viz_plot_manhattan_like import _add_pad_to_x_axis as _add_pad_to_x_axis
from gwaslab.viz.viz_plot_manhattan_like import _configure_cols_to_use as _configure_cols_to_use
from gwaslab.viz.viz_plot_manhattan_like import _sanity_check as _sanity_check
from gwaslab.viz.viz_plot_manhattan_like import _process_p_value as _process_p_value
from gwaslab.viz.viz_plot_manhattan_like import _process_highlight as _process_highlight
from gwaslab.viz.viz_plot_manhattan_like import _process_line as _process_line
from gwaslab.viz.viz_plot_manhattan_like import _process_cbar as _process_cbar
from gwaslab.viz.viz_plot_manhattan_like import _process_xtick as _process_xtick
from gwaslab.viz.viz_plot_manhattan_like import _process_ytick as _process_ytick
from gwaslab.viz.viz_plot_manhattan_like import _process_xlabel as _process_xlabel
from gwaslab.viz.viz_plot_manhattan_like import _process_ylabel as _process_ylabel
from gwaslab.viz.viz_plot_manhattan_like import _process_spine as _process_spine
from gwaslab.viz.viz_plot_manhattan_like import _process_layout as _process_layout
