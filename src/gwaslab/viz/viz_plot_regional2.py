import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_hex
import matplotlib.colors
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import seaborn as sns
import numpy as np
import copy
import re
import scipy as sp
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from adjustText import adjust_text
from gwaslab.io.io_gtf import read_gtf

from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_recombination_rate
from gwaslab.io.io_gtf import get_gtf
from gwaslab.viz.viz_aux_save_figure import safefig
from gwaslab.viz.viz_aux_style_options import set_plot_style
@safefig
def _plot_regional(
    sumstats,
    fig,
    ax1,
    ax3,
    region,
    vcf_path,
    marker_size,
    build,
    cut_line_color,
    ld_path=None,
    gtf_path="default",
    gtf_chr_dict = get_number_to_chr(),
    gtf_gene_name=None,
    rr_path="default",
    rr_header_dict=None,
    rr_chr_dict = get_number_to_chr(),
    rr_lim = (0,100),
    rr_ylabel = True,
    region_ld_legend=True,
    region_title=None,
    mode="mqq",
    region_step = 21,
    region_ref=None,
    region_ref_index_dic = None,
    region_ref_alias = None,
    region_grid = False,
    region_grid_line = {"linewidth": 2,"linestyle":"--"},
    region_lead_grid = True,
    region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"},
    region_title_kwargs = None,
    region_ld_threshold = [0.2,0.4,0.6,0.8],
    region_ld_colors = None,
    region_marker_shapes=None,
    cbar_fontsize=None,
    cbar_scale=False,
    cbar_bbox_to_anchor=None,
    cbar_w_scale=1,
    cbar_h_scale=1,
    cbar_downward_offset =1.3, 
    cbar_borderpad=None,
    cbar_equal_aspect=False,
    palette=None,
    region_recombination = True,
    region_protein_coding=True,
    region_legend_marker=True,
    region_flank_factor = 0.05,
    track_font_family="Arial",
    taf=[4,0,0.95,1,1],
    pos="POS",
    ld_link=False,
    ld_link_thr=0,
    ld_link_color="blue",
    ld_link_alpha_scale=0.2,
    ld_link_linewidth=1.0,
    ld_link_sig_level=None,
    show_ld_score=False,
    sig_level=5e-8,
    verbose=True,
    log=Log()
):
    """
    Create an Region plot. By default lead variant will be selected as reference variant for LD coloring.
    Parameters
    ----------
    vcf_path : str, optional, default=None
        Path to VCF file for LD reference in regional plot analysis. Default is None.

    gtf_path : str, optional, default='ensembl'
        Path to GTF file for gene annotation in regional plots. 
        gtf_path options:
        - 'ensembl' : GTF from ensembl. `build` should be specified. 
        - 'refseq' : GTF from refseq. `build` should be specified.
        - str : path for user provided gtf

    scaled : bool, optional, default=False
        Whether P-values are already scaled.
    scatter_kwargs : dict, optional
        Arguments for scatter plot styling. Default is None.
    scatterargs : dict, optional
        Alternative name for scatter plot styling. Default is None.
    region : tuple or list, required
        Genomic region specification (chr, start, end) for regional plots. Note start < end. No numeric expression is allowed.
        Region can be determined using `get_region_start_and_end`.
    region_title : str, optional
        Title for regional plot. Default is None.
    region_title_kwargs : dict, optional
        Arguments for regional plot title styling. Default is None.
    region_ref : str or list, optional
        Reference variant(s) for regional plot. If not provided, the lead variant
        in the region (highest `scaled_P`) is auto-selected as the default reference.
    region_ref2 : str, optional
        Second reference variant for regional plot. Default is None.
    region_ref_second : str, optional,optional
        Alternative name for second reference variant. Default is None.
    region_step : int,optional, default=21
        Step size for regional plot grid.
    region_grid : bool, optional,default=False
        Whether to show grid lines in regional plots.
    region_grid_line : dict,optional, default={"linewidth": 2,"linestyle":"--"}
        Styling arguments for regional plot grid lines.
    region_lead_grid : bool, optional,default=True
        Whether to show lead variant grid in regional plots.
    region_lead_grid_line : dict,optional, default={"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}
        Styling arguments for lead variant grid lines.
    region_hspace : float,optional, default=0.02
        Horizontal space between regional plot panels.
    region_ld_threshold : list,optional, default=[0.2,0.4,0.6,0.8]
        LD thresholds for regional plot coloring.
    region_ld_legend : bool, optional,default=True
        Whether to show LD legend in regional plot.
    region_ld_colors : list, optional,default=["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"]
        Colors for LD levels in regional plot. Maps to r² categories: missing data, low LD (r²<0.2),
        moderate LD (0.2-0.4), medium LD (0.4-0.6), high LD (0.6-0.8), very high LD (≥0.8),
        and reference variant highlighting.
    region_ld_colors_m : list, optional, default=["#E51819","#367EB7","green","#F07818","#AD5691","yellow","purple"]
        Colors for multi-lead variant regional plots.
    region_recombination : bool,optional, default=True
        Whether to plot recombination rate in regional plot.
    region_flank_factor : float, optional,default=0.05
        Flanking factor for regional plot.
    region_anno_bbox_kwargs : dict, optional,default={"ec":"None","fc":"None"}
        Bounding box arguments for regional plot annotations.
    region_marker_shapes : list, optional,default=['o', '^','s','D','*','P','X','h','8']
        Shapes for markers in regional plot. First one, the shape for non-reference variants, the second one is for the first reference variants and so on.
    region_legend_marker : bool, optional,default=True
        Whether to show marker legend in regional plot.
    region_ref_alias : dict, optional
        Alias mapping for reference variants in regional plot. Default is None.
    cbar_title : str, default='LD $\\mathregular{r^2}$ with variant'
        Title for color bar.
    cbar_fontsize : int, optional
        Font size for color bar labels. Default is None.
    cbar_scale : bool, default=True
        Whether to scale color bar.
    cbar_font_family : str, optional
        Font family for color bar. Default is None.
    cbar_bbox_to_anchor : tuple,optional, default=(0,0,1,1)
        Bounding box anchor for color bar.
    cbar_equal_aspect : bool, optional,default=True
        Whether to keep equal aspect in color bar.
    cbar_w_scale : float,optional, default=1
        Width scale for color bar.
    cbar_h_scale : float,optional, default=1
        Height scale for color bar.
    cbar_downward_offset : float, optional,default=1.3
        Downward offset for color bar.
    cbar_borderpad : float, optional
        Border padding for color bar. Default is None.
    track_n : int,optional, default=4
        Number of tracks in regional plot.
    track_n_offset : int,optional,default=0
        Track offset in regional plot.
    track_fontsize_ratio : float,optional, default=0.95
        Font size ratio for tracks.
    track_exon_ratio : int,optional, default=1
        Exon ratio for tracks.
    track_text_offset : int, optional,default=1
        Text offset for tracks.
    track_font_family : str, optional
        Font family for tracks. Default is None.
    windowsizekb : int, optional,default=500
        Window size in kb for obtainning significant variant.
    anno : str, optional,default=None.
        Annotation source (e.g., 'GENENAME'). Default is None.
    anno_set : list, optional
        Set of variants to annotate. Default is None.
    anno_alias : dict, optional
        Alias mapping for annotations. Default is None.
    anno_d : dict, optional
        Additional annotation data. Default is None.
    anno_kwargs : dict, optional
        Arguments for annotation styling. Default is None.
    anno_kwargs_single : dict, optional
        Arguments for single annotation styling. Default is None.
    anno_style : str, optional,default='right'
        Annotation style.
    anno_fixed_arm_length : float, optional
        Fixed arm length for annotations. Default is None.
    anno_adjust : bool, default=False
        Whether to adjust annotation positions for tight style.
    anno_xshift : float, optional
        X-axis shift for annotations. Default is None.
    anno_max_iter : int, optional, default=100
        Maximum iterations for annotation adjustment.
    arrow_kwargs : dict, optional
        Arguments for annotation arrows. Default is None.
    arm_offset : float, optional
        Offset for annotation arms. Default is None.
    arm_scale : float, optional,default=1
        Scale for annotation arms.
    anno_height : float, optional,default=1
        Height for annotations.
    arm_scale_d : float, optional
        Scale for annotation arms in density plots. Default is None.
    cut : float, optional,default=0
        Cut value for shrinking variants above threshold.
    skip : float,optional, default=0
        Skip variants with -log10(P) < skip.
    ystep : float, optional,default=0
        Step size for y-axis.
    ylabels : list, optional
        Custom y-axis labels. Default is None.
    ytick3 : bool, optional,default=True
        Whether to use 3 y-axis ticks.
    cutfactor : int, optional,default=10
        Factor for shrinking cut line.
    cut_line_color : str, optional,default='#ebebeb'
        Color for cut line.
    cut_log : bool,optional, default=False
        Whether to use log scale for cut line.
    sig_line : bool, optional,default=True
        Whether to show significance line.
    sig_level : float, optional, default=5e-8
        Significance level threshold to plot a line on the plot.
    anno_sig_level : float,optional, default=5e-8
        Significance level for lead variants.
    sig_line_color : str, optional,default='grey'
        Color for significance line.
    suggestive_sig_line : bool,optional, default=False
        Whether to show suggestive significance line.
    suggestive_sig_level : float,optional, default=5e-6
        Suggestive significance level.
    suggestive_sig_line_color : str,optional, default='grey'
        Color for suggestive significance line.
    additional_line : list, optional
        Additional lines to plot. Default is None.
    additional_line_color : list, optional
        Colors for additional lines. Default is None.
    sc_linewidth : int,optional, default=2
        Line width for significance lines.
    pinpoint : list, optional
        List of variants to pinpoint. Default is None.
    pinpoint_color : str, optional,default='red'
        Color for pinpointed variants.
    ylim : tuple, optional
        Y-axis limits for plot. Default is None.
    xpad : float, optional
        X-axis padding proportion. Default is None.
    xpadl : float, optional
        Left X-axis padding proportion. Default is None.
    xpadr : float, optional
        Right X-axis padding proportion. Default is None.
    xtight : bool, optional,default=False
        Whether to use tight X-axis padding.
    chrpad : float, optional,default=0.03
        Chromosome padding factor.
    title : str, optional
        Plot title. Default is None.
    xlabel : str, optional
        X-axis label. Default is None.
    title_pad : float,optional,default=1.08
        Padding for plot title.
    title_fontsize : int,optional, default=13
        Font size for title.
    fontsize : int, optional,default=9
        Base font size for plot elements.
    font_family : str, optional
        Font family for text elements. Default is None.
    fontfamily : str, optional,default='Arial'
        Alternative name for font family.
    math_fontfamily : str, optional,default='dejavusans'
        Font family for math text.
    anno_fontsize : int, optional,default=9
        Font size for annotations.
    figargs : dict, optional
        Figure arguments for subplots. Default is None.
    fig_kwargs : dict, optional
        Alternative name for figure arguments. Default is None.
    colors : list,optional, default=['#597FBD','#74BAD3']
        Color palette for plot.
    marker_size : tuple, optional,default=(45,65)
        Size range for markers.
    verbose : bool, optional,default=True
        Whether to show progress.
    repel_force : float, optional,default=0.03
        Force for repelling overlapping labels.
    build : str, optional
        Genomic build version. Default is None.
    dpi : int, optional, default=200
        Dots per inch for figure resolution.
    save : str, optional
        File path to save plot. Default is None.
    save_kwargs : dict, optional,default={"dpi":600,"transparent":True}
        Arguments for saving the plot.
    saveargs : dict, optional
        Alternative name for save arguments. Default is None.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib figure object.
    log : gwaslab.Log
        Updated logging object with operation records.
    lead_snp_is, lead_snp_is_color : tuple (optional)
        Lead SNP positions and colors returned only in regional ('r') mode.

    Less used parameters
    -------
    taf : dict, optional
        Track annotation file. Default is None.
    vcf_chr_dict : dict, optional
        Sumstats CHR (int) to VCF contig (str) dict for loading VCF data. 
        If None, it will be auto-detected (prefered). 
    gtf_chr_dict : dict, optional
        Sumstats CHR (int) to GTF chromosome (str) dict. Default is None. 
        Only used when user provided recombination rate data was used. Default is None. 
    gtf_gene_name : str, optional
        Gene name column in GTF data. Default is None.
        Only used when user provided recombination rate data was used.
    rr_path : str, optional
        Path to recombination rate data. Default is 'default'.
        rr_path options:
        - 'default' : data from hapmap3. `build` should be specified. (preferred)
        - str : path for user provided recombination rate data
    rr_header_dict : dict, optional
        Header mapping dictionary for recombination data. Default is None. 
        Only used when user provided recombination rate data was used.
    rr_chr_dict : dict, optional
        Sumstats CHR (int) to recombination data chromosome (string) dict. Default is None. 
        Only used when user provided recombination rate data was used.
    rr_lim : tuple, optional, default=(0,100)
        Recombination rate limits.
    rr_ylabel : bool,optional, default=True
        Whether to show recombination rate y-label.
    ld_path : str, optional,default=None
        Path to precomputed LD matrix for regional plot. Default is None.
    ld_map_path : str, optional,default=None
        Path to LD map data. Default is None.
    ld_fmt : str, default='npz'
        Format of LD data file.
    ld_if_square : bool, optional,default=False
        Whether LD matrix is square.
    ld_if_add_T : bool,optional, default=False
        Whether to add transpose for LD matrix.
    ld_map_rename_dic : dict, optional
        Dictionary to rename LD map columns. Default is None.
    ld_map_kwargs : dict, optional
        Additional arguments for LD map. Default is None.
    jagged : bool,optional, default=False
        Whether to make y-axis jagged.
    jagged_len : float, optional,default=0.01
        Length of jagged line.
    jagged_wid : float, optional,default=0.01
        Width of jagged line.
    anno_source : str,optional, default='ensembl'
        Source for annotations.
    anno_gtf_path : str, optional
        Path to GTF file for annotations. Default is None.
    ylabel : str, optional
        Y-axis label. Default is None.
    tabix : str, optional
        Path to tabix executable. Default is None.
    """
    style = set_plot_style(
        plot="plot_region",
        mode="r",
        fontfamily=track_font_family,
        verbose=verbose,
        log=log,
    )
    base_region_title = style.get("region_title_kwargs", {})
    explicit_region_title = region_title_kwargs if region_title_kwargs is not None else region_title_kwargs
    region_title_kwargs = dict(base_region_title)
    if explicit_region_title is not None:
        region_title_kwargs.update(explicit_region_title)
    region_title_kwargs = region_title_kwargs

    # Initialize lead_ids for LD score annotation
    lead_ids = []

    if (region is not None) :
        # track_n, track_n_offset,font_ratio,exon_ratio,text_offset
    # x axix: use i to plot (there is a gap between i and pos) 
    # if regional plot : pinpoint lead , add color bar ##################################################
        # pinpoint lead
        lead_ids = []  # Reset for this region
        
        for index, region_ref_single in enumerate(region_ref):
            ax1, lead_id_single = _pinpoint_lead(sumstats = sumstats,
                                        ax1 = ax1, 
                                        region_ref=region_ref_single,
                                        region_ref_total_n = len(region_ref), 
                                        lead_color = palette[(index+1)*100 + len(region_ld_threshold)+2], 
                                        marker_size= marker_size,
                                        region_marker_shapes=region_marker_shapes,
                                        log=log,verbose=verbose)
            #if lead_id_single is not None:
            lead_ids.append(lead_id_single)       
        
        # update region_ref to variant rsID or variantID / skip NAs
        new_region_ref = []
        for i in range(len(lead_ids)):
            if lead_ids[i] is None:
                new_region_ref.append(region_ref[i])
                continue
            if region_ref[i] is None:
                if "SNPID" in sumstats.columns:
                    new_name = sumstats.loc[lead_ids[i],"SNPID"]
                elif "rsID" in sumstats.columns:
                    new_name = sumstats.loc[lead_ids[i],"rsID"]
                else:
                    new_name = "chr{}:{}".format(sumstats.loc[lead_ids[i],"CHR"] , sumstats.loc[lead_ids[i],"POS"])
                new_region_ref.append(new_name)
                region_ref_index_dic[new_name] = region_ref_index_dic[region_ref[i]]
                continue
            else:
                new_region_ref.append(region_ref[i])
        region_ref = new_region_ref
        ##########################################################################################################

        ##########################################################################################################

        if ((vcf_path is not None) or (ld_path is not None)) and region_ld_legend:
            ## plot cbar
            ax1, cbar = _add_ld_legend(sumstats=sumstats, 
                            ax1=ax1, 
                            region_ref=region_ref,
                            region_ld_threshold=region_ld_threshold, 
                            region_ref_index_dic=region_ref_index_dic,
                            region_ref_alias=region_ref_alias, 
                            region_marker_shapes=region_marker_shapes,
                            cbar_fontsize= cbar_fontsize,
                            cbar_scale=cbar_scale,
                            cbar_equal_aspect=cbar_equal_aspect,
                            cbar_bbox_to_anchor=cbar_bbox_to_anchor,
                            cbar_w_scale=cbar_w_scale,
                            cbar_h_scale=cbar_h_scale,
                            cbar_downward_offset =cbar_downward_offset, 
                            cbar_borderpad=cbar_borderpad,
                            palette=palette,
                            region_legend_marker=region_legend_marker,
                            fig=fig)
        else:
            cbar=None

        if region_title is not None:
                ax1 = _add_region_title(region_title, ax1=ax1,region_title_kwargs=region_title_kwargs )
    
    ## recombinnation rate ##################################################       
    if (region is not None) and (region_recombination is True):
        ax4 = _plot_recombination_rate(sumstats = sumstats,
                                        pos =pos,
                                        region= region, 
                                        ax1 = ax1, 
                                        rr_path =rr_path, 
                                        rr_chr_dict = rr_chr_dict, 
                                        rr_header_dict =rr_header_dict, 
                                        build= build,
                                        rr_lim=rr_lim,
                                        rr_ylabel=rr_ylabel)
    else:
        ax4 = None
    
    ## LD link plot ##################################################       
    if (region is not None) and ld_link and (vcf_path is not None):
        from gwaslab.viz.viz_plot_ld_link import _plot_ld_link
        # Use sig_level as default for ld_link_sig_level if not specified
        if ld_link_sig_level is None:
            ld_link_sig_level = sig_level
        _plot_ld_link(
            ax=ax1,
            vcf_path=vcf_path,
            region=region,
            sumstats=sumstats,
            pos_col=pos,
            region_ld_threshold=region_ld_threshold,
            region_ld_colors=region_ld_colors,
            palette=palette,
            link_alpha_scale=ld_link_alpha_scale,
            link_linewidth=ld_link_linewidth,
            sig_level=ld_link_sig_level,
            log=log,
            verbose=verbose
        )
    
    ## LD score annotation and links ##################################################       
    if (region is not None) and show_ld_score and ((vcf_path is not None) or (ld_path is not None)):
        # Get reference panel sample size from VCF if available
        vcf_n_samples = None
        if vcf_path is not None:
            try:
                # Load VCF to get sample count (only need header/sample info)
                # Auto-detect vcf_chr_dict if not provided
                vcf_chr_dict_local, tabix_local = prepare_vcf_context(vcf_path, None, log, verbose)
                # Read a minimal VCF to get sample count
                ref_genotype_sample = read_vcf(vcf_path, 
                                             region=vcf_chr_dict_local[region[0]]+":"+str(region[1])+"-"+str(region[1]+1),
                                             tabix=tabix_local)
                if ref_genotype_sample is not None and "calldata/GT" in ref_genotype_sample:
                    # Get number of samples from genotype array shape: (n_variants, n_samples, ploidy)
                    vcf_n_samples = ref_genotype_sample["calldata/GT"].shape[1]
                    log.write(f" -Reference panel sample size (from VCF): {vcf_n_samples}", verbose=verbose)
            except Exception as e:
                log.warning(f"Could not determine VCF sample size: {e}. Using raw r² values.", verbose=verbose)
        
        ax1 = _plot_ld_score_annotation(
            ax=ax1,
            sumstats=sumstats,
            region_ref=region_ref,
            lead_ids=lead_ids,
            region_ld_threshold=region_ld_threshold,
            region_ld_colors=region_ld_colors,
            palette=palette,
            link_alpha_scale=ld_link_alpha_scale,
            link_linewidth=ld_link_linewidth,
            pos=pos,
            vcf_n_samples=vcf_n_samples,
            log=log,
            verbose=verbose
        )
     
    ## regional plot : gene track ######################################################################
                # calculate offset
    if (region is not None):
        most_left_snp      = sumstats["i"].idxmin()
        
        # distance between leftmost variant position to region left bound
        gene_track_offset  = sumstats.loc[most_left_snp,pos] - region[1]
        
        # rebase i to region[1] : the i value when POS=0
        gene_track_start_i = sumstats.loc[most_left_snp,"i"] - gene_track_offset - region[1]

        lead_snp_ys = []
        lead_snp_is = []
        lead_snp_is_colors = []
        for i,lead_id_single in enumerate(lead_ids):
            if lead_id_single is not None:
                lead_snp_ys.append(sumstats.loc[lead_id_single,"scaled_P"] )
                lead_snp_is.append(sumstats.loc[lead_id_single,"i"])
                lead_color = palette[(region_ref_index_dic[region_ref[i]]+1)*100 + len(region_ld_threshold) +1] # consistent color
                lead_snp_is_colors.append(lead_color)

    if gtf_path is not None:
        # For gene coloring, only use the first reference variant (region_ref[0]), not region_ref_second
        # lead_snp_is_for_gene_coloring will only contain the first lead SNP
        lead_snp_is_for_gene_coloring = []
        if len(lead_ids) > 0 and lead_ids[0] is not None:
            lead_snp_is_for_gene_coloring.append(sumstats.loc[lead_ids[0],"i"])
        # load gtf
        ax3, texts_to_adjust_middle =_plot_gene_track(
                        ax3=ax3,
                        fig=fig,
                        gtf_path=gtf_path,
                        region=region,
                        region_flank_factor=region_flank_factor,
                        region_protein_coding=region_protein_coding,
                        region_lead_grid=region_lead_grid,
                        region_lead_grid_line=region_lead_grid_line,
                        lead_snp_is=lead_snp_is_for_gene_coloring,
                        gene_track_start_i=gene_track_start_i,
                        gtf_chr_dict=gtf_chr_dict,
                        gtf_gene_name=gtf_gene_name, 
                        track_font_family=track_font_family,
                        taf=taf,
                        build=build, 
                        verbose=verbose, 
                        log=log)

    ## regional plot - set X tick
    if region is not None:
        region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 
        
        explicit = {"x","color","zorder"}
        region_grid_line = {k: v for k, v in region_grid_line.items() if k not in explicit}
        explicit = {"x","zorder"}
        region_lead_grid_line = {k: v for k, v in region_lead_grid_line.items() if k not in explicit}

        # set x ticks for gene track
        if "r" in mode:
            if gtf_path is not None: 
                ax3.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
                ax1.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
                ax3.set_xticklabels(region_ticks,rotation=45)
                ax1.set_xticklabels([],rotation=45)
            
            if region_grid==True:
                for i in np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step):
                    ax1.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
                    ax3.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
            
            if region_lead_grid==True:
                for lead_snp_i,  lead_snp_y, lead_snp_is_color in zip(lead_snp_is, lead_snp_ys , lead_snp_is_colors):
                    region_lead_grid_line["color"] = lead_snp_is_color
                    ax1.plot([lead_snp_i,lead_snp_i],[0,lead_snp_y], zorder=1,**region_lead_grid_line)
                    ax3.axvline(x=lead_snp_i, zorder=2,**region_lead_grid_line)

            ax1.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])
            ax3.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])
        else:
            # set x ticks m plot
            ax1.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
            ax1.set_xticklabels(region_ticks,rotation=45)
    
        ax1.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])

    # gene track (ax3) text adjustment
    if (gtf_path is not None ) and ("r" in mode):
        if len(texts_to_adjust_middle)>0:
            adjust_text(texts_to_adjust_middle,
                    autoalign=False, 
                    only_move={'points':'x', 'text':'x', 'objects':'x'},
                    ax=ax3,
                    precision=0,
                    force_text=(0.1,0),
                    expand_text=(1, 1),
                    expand_objects=(1,1),
                    expand_points=(1,1),
                    va="center",
                    ha='center',
                    avoid_points=False,
                    lim =1000)     
    
    return ax1, ax3, ax4, cbar, lead_snp_is, lead_snp_is_colors

def regional_mode_setup(
    sumstats,
    region_ref,
    region_ld_colors,
    region_marker_shapes,
    region_ld_colors_m,
    region_ld_threshold,
    vcf_path,
    ld_path,
    log=Log(),
    verbose=True,
):
    """
    Set up regional plotting parameters and filter data for regional association plots.

    This function configures the color palette, markers, and data filtering for regional
    plots. It handles LD-based coloring and ensures the reference variant is properly
    excluded from the main scatter plot when in single-reference mode.

    For single-reference regional plots, the reference variant is removed from the main
    Manhattan plot data (to_plot) so it can be highlighted separately as a pinpoint marker.
    When region_ref[0] is None, the lead variant (highest -log10(P)) is automatically
    selected and hidden. When region_ref contains a specific variant, that variant is
    found and hidden from the main plot.

    LD Color and Marker System:
    ---------------------------
    The function creates a mapping system where variants are colored and shaped based
    on their LD relationship to the reference variant(s).

    Single-Reference Mode (len(region_ref) == 1):
    - region_ld_colors: List of colors for LD categories (e.g., ['red', 'orange', 'yellow', 'blue'])
    - region_marker_shapes: List of marker shapes (e.g., ['o', '^', 's', 'D'] for circle, triangle, square, diamond)
    - palette: Maps LD values to colors (100+i -> region_ld_colors[i])
    - markers: Maps SHAPE values to marker shapes (1 -> 'o', 2 -> '^')
    - Variants get SHAPE=1 (circle) by default, reference variant gets SHAPE=2 (triangle)

    Multi-Reference Mode (len(region_ref) > 1):
    - region_ld_colors_m: List of colors, one per reference variant
    - region_marker_shapes: Extended to provide unique shapes for each reference
    - palette: Complex mapping with (ref_index+1)*100 + ld_level -> color
    - markers: Extended markers dict with shapes for each reference variant
    - Each reference variant gets a distinct color gradient and marker shape

    LD Thresholds:
    - region_ld_threshold: LD (r²) cutoffs for color categorization
    - Example: [0.2, 0.4, 0.6, 0.8] creates 5 LD categories
    - Variants are assigned to categories based on their r² with reference variant(s)

    Data Filtering:
    - to_plot: Filtered DataFrame with reference variant removed (single-reference mode only)
    - Reference variants are plotted separately as pinpoint markers with distinct styling
    - This prevents visual overlap and ensures clear highlighting of reference variants

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics DataFrame containing variant data with processed P-values.
    region_ref : list
        List of reference variant identifiers. For single-reference mode, the first
        element determines which variant to highlight. Can contain None for automatic
        lead variant detection.
    region_ld_colors : list
        Color palette for LD categories in single-reference regional plots.
    region_marker_shapes : list
        Marker shapes for different LD categories and reference variants.
    region_ld_colors_m : list
        Color palettes for multi-reference regional plots.
    region_ld_threshold : list
        LD (r²) thresholds for color categorization.
    vcf_path : str or None
        Path to VCF file for LD calculations.
    ld_path : str or None
        Path to pre-computed LD file.
    log : gwaslab.Log
        Logging object for recording messages.
    verbose : bool
        Whether to print progress messages.

    Returns
    -------
    legend : dict or None
        Legend configuration for the plot.
    linewidth : int
        Line width for plot elements.
    palette : dict
        Color mapping for LD categories.
    style : str
        Column name for marker style mapping.
    markers : dict
        Marker shape mapping for LD categories.
    to_plot : pandas.DataFrame
        Filtered DataFrame with reference variants removed from main scatter plot.
    edgecolor : str
        Edge color for markers.
    """
    legend = None
    linewidth = 1
    style = None
    edgecolor = "black"
    to_plot = None
    if vcf_path is None and ld_path is None:
        sumstats["LD"] = 100
        sumstats["SHAPE"] = 1
    sumstats["chr_hue"] = sumstats["LD"]
    if len(region_ref) == 1:
        palette = {100 + i: region_ld_colors[i] for i in range(len(region_ld_colors))}
        markers = {(i + 1): m for i, m in enumerate(region_marker_shapes[:2])}
        if region_ref[0] is None:
            id_to_hide = sumstats["scaled_P"].idxmax()
            to_plot = sumstats.drop(id_to_hide, axis=0)
        else:
            id_to_hide = _get_lead_id(sumstats, region_ref, log=log, verbose=verbose)
            if id_to_hide is not None:
                to_plot = sumstats.drop(id_to_hide, axis=0)
        style = "SHAPE"
    else:
        palette = {}
        region_color_maps = []
        for colorgroup in region_ld_colors_m:
            color_map_len = len(region_ld_threshold) + 2
            rgba = LinearSegmentedColormap.from_list("custom", ["white", colorgroup], color_map_len)(range(1, color_map_len))
            output_hex_colors = []
            for i in range(len(rgba)):
                output_hex_colors.append(to_hex(rgba[i]))
                region_ld_colors_single = [region_ld_colors[0]] + output_hex_colors + [output_hex_colors[-1]]
            region_color_maps.append(region_ld_colors_single)
        for i, hex_colors in enumerate(region_color_maps):
            for j, hex_color in enumerate(hex_colors):
                palette[(i + 1) * 100 + j] = hex_color
        # Add key 0 for variants with missing LD for all references
        #palette[0] = region_ld_colors[0]  # Use missing LD color
        edgecolor = "none"
        markers = {(i + 1): m for i, m in enumerate(region_marker_shapes[: len(region_ref)])}
        style = "SHAPE"
        # For multi-reference mode, hide all reference variants from main plot
        ids_to_hide = []
        for ref in region_ref:
            if ref is not None:
                id_to_hide = _get_lead_id(sumstats, [ref], log=log, verbose=verbose)
                if id_to_hide is not None:
                    ids_to_hide.append(id_to_hide)
        if len(ids_to_hide) > 0:
            to_plot = sumstats.drop(ids_to_hide, axis=0)
        else:
            to_plot = sumstats
    return legend, linewidth, palette, style, markers, to_plot, edgecolor

# + ###########################################################################################################################################################################
def _get_lead_id(sumstats=None, region_ref=None, log=None, verbose=True):
    """
    Retrieve the lead variant index from sumstats based on rsID/SNPID/CHR:POS:NEA:EA in region_ref.
    
    Parameters
    ----------
    sumstats : pandas.DataFrame, optional
        Summary statistics DataFrame containing variant data.
    region_ref : str or list, optional
        Reference variant(s) to search for. If None, uses the variant with maximum scaled P-value.
    log : gwaslab.Log, optional
        Logging object for recording messages. Default is None.
    verbose : bool, optional
        Whether to show progress. Default is True.
    
    Returns
    -------
    lead_id : int or None
        Index of the lead variant if found, otherwise None.
    """
    region_ref_to_check = copy.copy(region_ref)
    # region_ref_single (not none) -> specified variant ID
    # convert region_ref_single -> lead_id(index)
    
    #
    try: 
        if len(region_ref_to_check)>0 and type(region_ref_to_check) is not str:
            region_ref_to_check = region_ref_to_check[0]
    except:
        pass
    
    # index of lead variant
    lead_id=None
    
    # match by rsID 
    if "rsID" in sumstats.columns:
        lead_id = sumstats.index[sumstats["rsID"] == region_ref_to_check].to_list()
    # match by SNPID
    if lead_id is None and "SNPID" in sumstats.columns:
        lead_id = sumstats.index[sumstats["SNPID"] == region_ref_to_check].to_list()

    # if duplicated, select the first one
    if type(lead_id) is list:
        if len(lead_id)>0:
            lead_id = int(lead_id[0])


    if region_ref_to_check is not None:
        if type(lead_id) is list:
            if len(lead_id)==0 :
                #try:
                # if region_ref_to_check is in CHR:POS:NEA:EA format
                matched_snpid = re.match("(chr)?[0-9]+:[0-9]+:[ATCG]+:[ATCG]+", region_ref_to_check,  re.IGNORECASE)    
                if matched_snpid is None:
                    # if not, pass
                    pass
                else:
                    # if region_ref_to_check is in CHR:POS:NEA:EA format, match by CHR:POS:NEA:EA
                    lead_snpid = matched_snpid.group(0).split(":")
                    if len(lead_snpid)==4:
                        lead_chr= int(lead_snpid[0])
                        lead_pos= int(lead_snpid[1])
                        lead_ea= lead_snpid[2]
                        lead_nea= lead_snpid[3]
                        chrpos_match = (sumstats["CHR"] == lead_chr) & (sumstats["POS"] == lead_pos)
                        eanea_match = ((sumstats["EA"] == lead_ea) & (sumstats["NEA"] == lead_nea)) | ((sumstats["EA"] == lead_nea) & (sumstats["NEA"] == lead_ea)) 
                        if "rsID" in sumstats.columns:
                            lead_id = sumstats.index[chrpos_match&eanea_match].to_list()
                        if "SNPID" in sumstats.columns:
                            lead_id = sumstats.index[chrpos_match&eanea_match].to_list()     
                if type(lead_id) is list:
                    if len(lead_id)>0:
                        lead_id = int(lead_id[0])   
                        log.warning("Trying matching variant {} using CHR:POS:EA:NEA to {}... ".format(region_ref_to_check,lead_id))

        if type(lead_id) is list:
            if len(lead_id)==0 :
                log.warning("Extracting variant: {} not found in sumstats.. Skipping..".format(region_ref_to_check))
                #lead_id = sumstats["scaled_P"].idxmax()
                lead_id = None
                return lead_id
        else:
            log.write(" -Reference variant ID: {} with Index {}".format(region_ref_to_check, lead_id), verbose=verbose)

    if lead_id is None:
        log.write(" -Extracting lead variant...", verbose=verbose)
        lead_id = sumstats["scaled_P"].idxmax()

    return lead_id

def _pinpoint_lead(sumstats,ax1,region_ref, region_ref_total_n, lead_color, marker_size, log, verbose, region_marker_shapes):
    """
    Extract lead_id
    Draw a single marker for the lead_id
    """
    if region_ref is None:
        log.write(" -Extracting lead variant..." , verbose=verbose)
        lead_id = sumstats["scaled_P"].idxmax()
    else:
        lead_id = _get_lead_id(sumstats, region_ref, log, verbose)
    
    if lead_id is not None:
        if region_ref_total_n <2:
            # single-ref mode: just use SHAPE
            marker_shape = region_marker_shapes[sumstats.loc[lead_id,"SHAPE"]]
        else:
            # multi-ref mode: just use SHAPE - 1
            marker_shape = region_marker_shapes[sumstats.loc[lead_id,"SHAPE"]-1]

    if lead_id is not None:
        ax1.scatter(sumstats.loc[lead_id,"i"],sumstats.loc[lead_id,"scaled_P"],
                color=lead_color,
                zorder=3,
                marker= marker_shape,
                s=marker_size[1]*1.5,
                edgecolor="black")

    return ax1, lead_id
# -############################################################################################################################################################################
def _add_region_title(region_title, ax1,region_title_kwargs):
    explicit = {"region_title","transform","va","ha","fontsize"}
    region_title_kwargs = {k: v for k, v in region_title_kwargs.items() if k not in explicit}
    ax1.text(0.015,0.97, region_title, transform=ax1.transAxes, va="top", ha="left", **region_title_kwargs )
    return ax1

def _add_ld_legend(sumstats, ax1, region_ld_threshold, region_ref,region_ref_index_dic,region_marker_shapes,fig, region_legend_marker=True,
                   cbar_fontsize= None,cbar_scale=False,cbar_equal_aspect=True,cbar_w_scale=1,cbar_h_scale=1,palette =None, 
                   cbar_downward_offset =1.2, cbar_borderpad=None,
                   cbar_bbox_to_anchor=(0, 0, 1, 1),region_ref_alias=None):

    scale = 1
    if cbar_scale:
        base_fontsize = 9
        scale = cbar_fontsize / base_fontsize
        scale = max(1,scale)
    else:
        scale = 1
    
    width_raw= 11 * (scale)*cbar_w_scale
    height_raw=(7 + 7 * len(region_ref))*(scale)*cbar_h_scale

    width_pct = "{}%".format(width_raw)
    height_pct = "{}%".format( height_raw)

    total_y_pixels =(ax1.bbox.get_points()[1][1]-ax1.bbox.get_points()[0][1]) 
    downwards_offset = cbar_fontsize / (total_y_pixels/ fig.dpi * 72) * cbar_downward_offset
    bbox_to_anchor = (cbar_bbox_to_anchor[0],cbar_bbox_to_anchor[1]-downwards_offset,cbar_bbox_to_anchor[2],cbar_bbox_to_anchor[3])
    
    if cbar_borderpad is None:
        borderpad=0.5*(scale)
    else:
        borderpad=cbar_borderpad

    axins1 = inset_axes(ax1,
            width=width_pct,  # width = 50% of parent_bbox width
            height=height_pct,  # height : 5%
            bbox_to_anchor=bbox_to_anchor,
            bbox_transform=ax1.transAxes,
            borderpad=borderpad,
            loc='upper right',
            axes_kwargs={"frameon":True,"facecolor":"white","zorder":999999,"anchor":"NE"})

    ld_ticks = [0]+region_ld_threshold+[1]

    for index, ld_threshold in enumerate(ld_ticks):
        for group_index in range(len(region_ref)):
            if index < len(ld_ticks)-1:            
                x=ld_threshold
                y=0.2*group_index 
                width=0.2
                height=(ld_ticks[index+1]-ld_ticks[index]) 
                hex_color = palette[(region_ref_index_dic[region_ref[group_index]]+1)*100 + index+1] # consistent color
                
                a = Rectangle((x,y),width, height, fill = True, color = hex_color , linewidth = 2)
                #patches.append(a)
                axins1.add_patch(a)
    
    # y snpid
    if region_ref_alias is None:
        region_ref_name = region_ref
    else:
        region_ref_name = [region_ref_alias[i] for i in region_ref]

    yticks_position = (0.1 + 0.2 *np.arange(0,len(region_ref_name)))
    axins1.set_yticks(yticks_position, ["{}".format(x) for x in region_ref_name])
    axins1.set_ylim(0,0.2*len(region_ref_name))    
    ymin, ymax=0,0.2*len(region_ref_name)
    # x ld thresholds
    
    axins1.set_xticks(ticks=ld_ticks) 
    axins1.set_xticklabels([str(i) for i in ld_ticks]) 
    xmin, xmax = 0, 1
    axins1.set_xlim(xmin,xmax)       

    if cbar_equal_aspect==True:
        axins1.set_aspect('equal', adjustable='box',anchor="NE")

    ############### ##############plot marker ############## ##############
    if region_legend_marker==True:
        for group_index, ref in enumerate(region_ref):

            data_to_point_y =((axins1.bbox.get_points()[1][1]-axins1.bbox.get_points()[0][1])*height_raw/(ymax -ymin))
            data_to_point_x =((axins1.bbox.get_points()[1][0]-axins1.bbox.get_points()[0][0])*width_raw/(xmax -xmin))
            y_to_x = data_to_point_y/data_to_point_x
            x_to_y = 1/y_to_x
            xyratio = min(y_to_x, x_to_y)
            
            marker_side_in_data = 0.075
            if cbar_equal_aspect==True:
                xyratio=1
            
            ## change markersize

            if xyratio <1 :
                x = 0 - (marker_side_in_data +0.03) * xyratio
            else:
                x = 0 - (marker_side_in_data +0.03)
            y= (0.1 + 0.2 * group_index)
            
            if len(region_ref) <2:
                # single-ref mode
                marker = region_marker_shapes[group_index+1]
                c =  palette[(region_ref_index_dic[region_ref[group_index]]+1)*100 + len(ld_ticks)]
            else:
                # multi-ref mode
                marker = region_marker_shapes[group_index]
                c =  palette[(region_ref_index_dic[region_ref[group_index]]+1)*100 + len(ld_ticks)-1]
            
            # ([x0,y0][x1,y1])
            #  y pixels / per data 1
            
            data_to_point_y =((axins1.bbox.get_points()[1][1]-axins1.bbox.get_points()[0][1])*height_raw/(ymax -ymin))
            data_to_point_x =((axins1.bbox.get_points()[1][0]-axins1.bbox.get_points()[0][0])*width_raw/(xmax -xmin))
            
            if data_to_point_y < data_to_point_x:
                length_raw = 1 #height_raw
                data_to_point = data_to_point_y
            else:
                length_raw = 1 #width_raw
                data_to_point = data_to_point_x

            # pixels/data 1 -> font points/data 1  
            #  (dpi / 72) = point_per_pixel
            # y pixels / per data 1 / (dpi / 72)  -> y font points/data 1  
            
            font_points_per_data_1 = data_to_point/(fig.dpi/72)
            s =  ((marker_side_in_data*2)* font_points_per_data_1 * length_raw/100 )**2     
            
            axins1.scatter(x, y, s=s, marker=marker,c=c, edgecolors="black", linewidths = 1,  clip_on=False, zorder=100)

            pad = ((marker_side_in_data*2+0.02)* font_points_per_data_1 * length_raw/100)
            tick_length=(abs(x)* font_points_per_data_1 * length_raw/100)
            axins1.tick_params(axis="y", pad=pad-0.5*tick_length, length=tick_length) 

    cbar = axins1
    return ax1, cbar

# -############################################################################################################################################################################
def _plot_ld_score_annotation(
    ax,
    sumstats,
    region_ref,
    lead_ids,
    region_ld_threshold,
    region_ld_colors,
    palette,
    link_alpha_scale,
    link_linewidth,
    pos,
    vcf_n_samples=None,
    log=Log(),
    verbose=True
):
    """
    Annotate ref variants with LD scores and draw links from each variant to ref variants.
    
    LD score is calculated using the LDSC unbiased estimator: 
    l_j = Σ[r²_{j,k} - (1-r²_{j,k})/(n-2)] for all k ≠ j, where r²_{j,k} is the squared 
    correlation (linkage disequilibrium) between ref variant j and variant k, and n is the 
    reference panel sample size (from VCF). If VCF sample size is not available, uses raw 
    r² values: l_j = Σ(r²_{j,k}).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    sumstats : pandas.DataFrame
        Summary statistics with LD information (RSQ_0, RSQ_1, etc.)
        RSQ columns contain r² values between each variant and the corresponding ref variant
    region_ref : list
        List of reference variant identifiers
    lead_ids : list
        List of lead variant indices in sumstats
    region_ld_threshold : list
        LD r² thresholds for color categories
    region_ld_colors : list
        Colors for each LD category
    palette : dict
        Color palette dictionary
    link_alpha_scale : float
        Scaling factor for line transparency
    link_linewidth : float
        Line width for links
    pos : str
        Column name for position
    vcf_n_samples : int, optional
        Reference panel sample size from VCF. If provided and > 2, applies LDSC bias correction.
    log : gwaslab.Log
        Logging object
    verbose : bool
        Whether to show progress
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        Modified axes
    """
    log.write("Adding LD score annotations and links...", verbose=verbose)
    
    # Use region_ld_colors if provided, otherwise extract from palette, otherwise use default
    if region_ld_colors is not None and isinstance(region_ld_colors, list):
        region_ld_colors_for_link = region_ld_colors
    elif palette is not None and isinstance(palette, dict):
        # Extract colors from palette dictionary
        sorted_keys = sorted([k for k in palette.keys() if isinstance(k, (int, float)) and k >= 100 and k < 200])
        if sorted_keys and len(sorted_keys) >= len(region_ld_threshold) + 3:
            region_ld_colors_for_link = [palette[k] for k in sorted_keys]
        else:
            region_ld_colors_for_link = ["#E4E4E4", "#020080", "#86CEF9", "#24FF02", "#FDA400", "#FF0000", "#FF0000"]
    else:
        region_ld_colors_for_link = ["#E4E4E4", "#020080", "#86CEF9", "#24FF02", "#FDA400", "#FF0000", "#FF0000"]
    
    # Ensure we have enough colors
    min_colors_needed = len(region_ld_threshold) + 3
    if len(region_ld_colors_for_link) < min_colors_needed:
        region_ld_colors_for_link = region_ld_colors_for_link + [region_ld_colors_for_link[-1]] * (min_colors_needed - len(region_ld_colors_for_link))
    
    region_ld_colors = region_ld_colors_for_link
    
    # Function to get color index based on LD value
    def get_ld_color_index(ld_value):
        """Get color index for LD value based on thresholds."""
        if pd.isna(ld_value) or ld_value <= 0:
            return 0
        for idx in range(len(region_ld_threshold) - 1, -1, -1):
            threshold = region_ld_threshold[idx]
            if ld_value > threshold:
                return idx + 2
        return 1
    
    # Process each reference variant
    for ref_n, (region_ref_single, lead_id) in enumerate(zip(region_ref, lead_ids)):
        if lead_id is None:
            continue
        
        rsq_col = "RSQ_{}".format(ref_n)
        
        # Check if RSQ column exists
        if rsq_col not in sumstats.columns:
            log.warning(f"RSQ column {rsq_col} not found. Skipping LD score annotation for ref variant {ref_n}.", verbose=verbose)
            continue
        
        # Get ref variant position and y-coordinate
        ref_x = sumstats.loc[lead_id, "i"]
        ref_y = sumstats.loc[lead_id, "scaled_P"]
        
        # Calculate LD score for ref variant: sum of bias-corrected r² values with all other variants
        # LD score l_j = Σ[r²_{j,k} - (1-r²_{j,k})/(n-2)] for all k ≠ j (LDSC unbiased estimator)
        # Exclude the ref variant itself (which would be r²=1.0) and NaN/zero values
        rsq_values = sumstats[rsq_col].copy()
        # Exclude the ref variant itself and invalid values
        rsq_values_excluding_ref = rsq_values.drop(lead_id)
        rsq_values_excluding_ref = rsq_values_excluding_ref.dropna()
        rsq_values_excluding_ref = rsq_values_excluding_ref[rsq_values_excluding_ref > 0]
        
        # Apply LDSC bias correction using reference panel sample size from VCF
        if vcf_n_samples is not None and vcf_n_samples > 2:
            # LDSC unbiased estimator: r²_unbiased = r² - (1-r²)/(n-2)
            # where n is the reference panel sample size
            denom = vcf_n_samples - 2
            # Apply bias correction to each r² value
            rsq_unbiased = rsq_values_excluding_ref - (1 - rsq_values_excluding_ref) / denom
            # Calculate LD score as sum of bias-corrected r² values
            ld_score_value = rsq_unbiased.sum()
        else:
            # No VCF sample size available, use raw r² values
            if vcf_n_samples is None:
                log.warning("Reference panel sample size not available. Using raw r² values without bias correction.", verbose=verbose)
            else:
                log.warning(f"Reference panel sample size ({vcf_n_samples}) <= 2. Using raw r² values without bias correction.", verbose=verbose)
            ld_score_value = rsq_values_excluding_ref.sum()
        
        # Get ref variant name for annotation
        if "SNPID" in sumstats.columns:
            ref_name = sumstats.loc[lead_id, "SNPID"]
        elif "rsID" in sumstats.columns:
            ref_name = sumstats.loc[lead_id, "rsID"]
        else:
            ref_name = f"chr{sumstats.loc[lead_id, 'CHR']}:{sumstats.loc[lead_id, pos]}"
        
        # Format LD score for display with academic notation and SNP ID
        # LD score formula: l_j = Σ(r²_{j,k}) for all k ≠ j
        ld_score_text = f"{ref_name}\n$l_j = {ld_score_value:.2f}$"
        
        # Calculate offset to position annotation above the marker
        # Get the y-axis range to determine appropriate offset
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        # Use 2% of y-axis range as vertical offset
        y_offset = y_range * 0.02
        
        # Add text annotation showing LD score above the marker
        ax.text(ref_x, ref_y + y_offset, ld_score_text, 
                ha='center', va='bottom', 
                fontsize=11, 
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='black', linewidth=1),
                zorder=10)
        
        # Draw links from each variant to the ref variant
        line_count = 0
        for idx, row in sumstats.iterrows():
            # Skip the ref variant itself
            if idx == lead_id:
                continue
            
            # Get LD score for this variant with ref variant
            ld_score = row[rsq_col]
            
            # Skip if LD score is NaN or 0
            if pd.isna(ld_score) or ld_score <= 0:
                continue
            
            # Get variant position and y-coordinate
            var_x = row["i"]
            var_y = row["scaled_P"]
            
            # Get color based on LD category
            color_idx = get_ld_color_index(ld_score)
            if color_idx < len(region_ld_colors):
                line_color = region_ld_colors[color_idx]
            else:
                line_color = region_ld_colors[-1]
            
            # Alpha based on LD value
            alpha = min(ld_score * link_alpha_scale, 1.0)
            
            # Draw line from variant to ref variant
            ax.plot([var_x, ref_x], [var_y, ref_y], 
                   color=line_color, 
                   alpha=alpha, 
                   linewidth=link_linewidth, 
                   zorder=1)
            line_count += 1
        
        log.write(f"Plotted {line_count} links to ref variant {ref_n} ({ref_name}), LD score = {ld_score_value:.2f}", verbose=verbose)
    
    return ax

# -############################################################################################################################################################################
def  _plot_recombination_rate(sumstats,pos, region, ax1, rr_path, rr_chr_dict, rr_header_dict, build,rr_lim, rr_ylabel=True):
    ax4 = ax1.twinx()
    most_left_snp = sumstats["i"].idxmin()
    
    # the i value when pos=0
    rc_track_offset = sumstats.loc[most_left_snp,"i"]-sumstats.loc[most_left_snp,pos]

    if rr_path=="default":
        if rr_chr_dict is not None:
            rr_chr = rr_chr_dict[region[0]]
        else:
            rr_chr = str(region[0])
        rc = get_recombination_rate(chrom=rr_chr,build=build)
    else:
        rc = pd.read_csv(rr_path,sep="\t")
        if rr_header_dict is not None:
            rc = rc.rename(columns=rr_header_dict)

    rc = rc.loc[(rc["Position(bp)"]<region[2]) & (rc["Position(bp)"]>region[1]),:]
    ax4.plot(rc_track_offset+rc["Position(bp)"],rc["Rate(cM/Mb)"],color="#5858FF",zorder=1)
    
    ax1.set_zorder(ax4.get_zorder()+1)
    ax1.patch.set_visible(False)
    
    if rr_ylabel:
        ax4.set_ylabel("Recombination rate(cM/Mb)")
    if rr_lim!="max":
        ax4.set_ylim(rr_lim[0],rr_lim[1])
    else:
        ax4.set_ylim(0, 1.05 * rc["Rate(cM/Mb)"].max())
    ax4.spines["top"].set_visible(False)
    ax4.spines["top"].set(zorder=1) 
    return ax4

# -############################################################################################################################################################################
def _plot_gene_track(
    ax3,
    fig,
    gtf_path,
    region,
    region_flank_factor,
    region_protein_coding,
    region_lead_grid,
    region_lead_grid_line,
    lead_snp_is,
    gene_track_start_i,
    gtf_chr_dict,gtf_gene_name, 
    track_font_family,
    taf,
    build, 
    verbose=True, 
    log=Log()):
    
    # load gtf
    log.write(" -Loading gtf files from:" + gtf_path, verbose=verbose)
    uniq_gene_region,exons = process_gtf(   gtf_path = gtf_path ,
                                            region = region,
                                            region_flank_factor = region_flank_factor,
                                            build=build,
                                            region_protein_coding=region_protein_coding,
                                            gtf_chr_dict=gtf_chr_dict,
                                            gtf_gene_name=gtf_gene_name,
                                            verbose=verbose,
                                            log=log)

    n_uniq_stack = uniq_gene_region["stack"].nunique()
    stack_num_to_plot = max(taf[0],n_uniq_stack)
    ax3.set_ylim((-stack_num_to_plot*2-taf[1]*2,2+taf[1]*2))
    ax3.set_yticks([])
    point_per_pixels = 72/fig.dpi
    pixels_per_point = fig.dpi/72

    pixels_per_track = np.abs(ax3.transData.transform([0,0])[1] - ax3.transData.transform([0,1])[1])                   
    font_size_in_pixels= taf[2] * pixels_per_track
    font_size_in_points =  font_size_in_pixels * point_per_pixels

    linewidth_in_points_per_track=   pixels_per_track * point_per_pixels
    
    log.write(" -plotting gene track..", verbose=verbose)
    
    sig_gene_name = "Undefined"
    sig_gene_name2 = "Undefined"
    texts_to_adjust_left = []
    texts_to_adjust_middle = []
    texts_to_adjust_right = []

    
    sig_gene_names=[]
    sig_gene_lefts=[]
    sig_gene_rights=[]
    log.write(" -plotting genes: {}..".format(len(uniq_gene_region)), verbose=verbose)
    for index,row in uniq_gene_region.iterrows():

        gene_color="#020080"
        #if row[6][0]=="+":
        if row["strand"][0]=="+":
            gene_anno = row["name"] + "->"
        else:
            gene_anno = "<-" + row["name"] 
        


        for lead_snp_i in lead_snp_is:
            if region_lead_grid is True and lead_snp_i > gene_track_start_i+row["start"] and lead_snp_i < gene_track_start_i+row["end"] :
                gene_color=region_lead_grid_line["color"]
                sig_gene_names.append(row["name"])
                sig_gene_lefts.append(gene_track_start_i+row["start"])
                sig_gene_rights.append(gene_track_start_i+row["end"])

        # plot gene line
        ## minimum width = 2 pixel 
        gene_line_width = max(linewidth_in_points_per_track/10, 2/pixels_per_point)
        
        ax3.plot((gene_track_start_i+row["start"],gene_track_start_i+row["end"]),
                    (row["stack"]*2,row["stack"]*2),color=gene_color,linewidth=gene_line_width,solid_capstyle="butt")

        # plot gene name
        if row["end"] >= region[2]:
            #right side
            texts_to_adjust_right.append(ax3.text(x=gene_track_start_i+region[2],
                    y=row["stack"]*2+taf[4],s=gene_anno,ha="right",va="center",color="black",style='italic', size=font_size_in_points,family=track_font_family))

        elif row["start"] <= region[1] :
            #left side
            texts_to_adjust_left.append(ax3.text(x=gene_track_start_i+region[1],
                    y=row["stack"]*2+taf[4],s=gene_anno,ha="left",va="center",color="black",style='italic', size=font_size_in_points,family=track_font_family))
        else:
            texts_to_adjust_middle.append(ax3.text(x=(gene_track_start_i+row["start"]+gene_track_start_i+row["end"])/2,
                    y=row["stack"]*2+taf[4],s=gene_anno,ha="center",va="center",color="black",style='italic',size=font_size_in_points,family=track_font_family))
    
    # plot exons
    log.write(" -plotting exons: {}..".format(len(exons)), verbose=verbose)
    for index,row in exons.iterrows():
        exon_color="#020080" 
        for sig_gene_name, sig_gene_left, sig_gene_right in zip(sig_gene_names,sig_gene_lefts,sig_gene_rights):
            
            if not pd.isnull(row["name"]):
                if (region_lead_grid is True) and row["name"]==sig_gene_name:
                    exon_color = region_lead_grid_line["color"]  
                else:
                    exon_color="#020080" 
            elif gene_track_start_i+row["start"] > sig_gene_left and gene_track_start_i+row["end"] < sig_gene_right:
                exon_color = region_lead_grid_line["color"]  
            else:
                exon_color="#020080"

        ## minimum width = 8 pixel 
        exon_line_width = max(linewidth_in_points_per_track * taf[3],  8/pixels_per_point) 

        ax3.plot((gene_track_start_i+row["start"],gene_track_start_i+row["end"]),
                    (row["stack"]*2,row["stack"]*2),linewidth=exon_line_width,color=exon_color,solid_capstyle="butt")

    log.write(" -Finished plotting gene track..", verbose=verbose)

    return ax3,texts_to_adjust_middle

# -############################################################################################################################################################################
# Helpers    
# -############################################################################################################################################################################
def process_vcf(sumstats, 
                vcf_path, 
                region,
                region_ref, 
                #region_ref_second, 
                log, 
                verbose, 
                pos ,
                nea,
                ea, 
                region_ld_threshold, 
                vcf_chr_dict,
                tabix):
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)

    # load genotype data of the targeted region
    ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
    if ref_genotype is None:
        log.warning("No data was retrieved. Skipping ...")
        ref_genotype=dict()
        ref_genotype["variants/POS"]=np.array([],dtype="int64")
    log.write(" -Retrieving index...", verbose=verbose)
    log.write(" -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
    # match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
     #######################################################################################
    def match_varaint(x):
        # x: "POS,NEA,EA"
        if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
            # position match
            if len(np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0])>1:  
            # multiple position matches
                for j in np.where(ref_genotype["variants/POS"] == x.iloc[0])[0]:
                # for each possible match, compare ref and alt
                    if x.iloc[1] == ref_genotype["variants/REF"][j]:
                        if x.iloc[2] in ref_genotype["variants/ALT"][j]:
                            return j
                    elif x.iloc[1] in ref_genotype["variants/ALT"][j]:
                        if x.iloc[2] == ref_genotype["variants/REF"][j]:
                            return j    
                return None
            else: 
                # single match
                return np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0][0] 
        else:
            # no position match
            return None
    log.write(" -Matching variants using POS, NEA, EA ...", verbose=verbose)
    #############################################################################################
    sumstats["REFINDEX"] = sumstats[[pos,nea,ea]].apply(lambda x: match_varaint(x),axis=1)
    #############################################################################################

    #for loop to add LD information
    #############################################################################################
    for ref_n, region_ref_single in enumerate(region_ref):

        rsq = "RSQ_{}".format(ref_n)
        ld_single = "LD_{}".format(ref_n)
        lead = "LEAD_{}".format(ref_n)
        sumstats[lead]= 0

        # get lead variant id and pos
        if region_ref_single is None:
            # if not specified, use lead variant
            lead_id = sumstats["scaled_P"].idxmax()
        else:
            # figure out lead variant
            lead_id = _get_lead_id(sumstats, region_ref_single, log, verbose)
        
        lead_series = None
        if lead_id is None:
            
            matched_snpid = re.match("(chr)?[0-9]+:[0-9]+:[ATCG]+:[ATCG]+",region_ref_single,  re.IGNORECASE)
            
            if matched_snpid is None:
                sumstats[rsq] = None
                sumstats[rsq] = sumstats[rsq].astype("float")
                sumstats[ld_single] = 0    
                continue    
            else:
                
                lead_snpid = matched_snpid.group(0).split(":")[1:]
                lead_pos = int(lead_snpid[0])
                lead_snpid[0]= int(lead_snpid[0])
                lead_series = pd.Series(lead_snpid)
        else:
            lead_pos = sumstats.loc[lead_id,pos]

        
        # if lead pos is available: 
        if lead_pos in ref_genotype["variants/POS"]:
            
            # get ref index for lead snp
            if lead_series is None:
                lead_snp_ref_index = match_varaint(sumstats.loc[lead_id,[pos,nea,ea]])
                #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]
            else:
                log.warning("Computing LD: {} not found in sumstats but found in reference...Still Computing...".format(region_ref_single))
                lead_snp_ref_index = match_varaint(lead_series)

            # non-na other snp index
            other_snps_ref_index = sumstats["REFINDEX"].dropna().astype("int").values
            # get genotype 
            
            lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
            try:
                if len(set(lead_snp_genotype[0]))==1:
                    log.warning("The variant is mono-allelic in reference VCF. LD can not be calculated.")
            except:
                pass
            other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
            
            log.write(" -Calculating Rsq...", verbose=verbose)
            
            if len(other_snp_genotype)>1:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
            else:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)
            sumstats.loc[~sumstats["REFINDEX"].isna(),rsq] = valid_r2
        else:
            log.write(" -Lead SNP not found in reference...", verbose=verbose)
            sumstats[rsq]=None
            
            # 
            try:
                sumstats.loc[lead_id,rsq]=1
            except KeyError:
                pass
        
        sumstats[rsq] = sumstats[rsq].astype("float")
        sumstats[ld_single] = 0
        
        for index,ld_threshold in enumerate(region_ld_threshold):
            # No data,LD = 0
            # 0, 0.2  LD = 1
            # 1, 0.4  LD = 2
            # 2, 0.6  LD = 3
            # 3, 0.8  LD = 4
            # 4, 1.0  LD = 5
            # lead    LD = 6

            if index==0:
                to_change_color = sumstats[rsq]>-1
                sumstats.loc[to_change_color,ld_single] = 1
            to_change_color = sumstats[rsq]>ld_threshold
            sumstats.loc[to_change_color,ld_single] = index+2
        
        if lead_series is None:
            sumstats.loc[lead_id,ld_single] = len(region_ld_threshold)+2
            sumstats.loc[lead_id,lead] = 1

    ####################################################################################################    
    final_shape_col = "SHAPE"
    final_ld_col = "LD"
    final_rsq_col = "RSQ"

    sumstats[final_ld_col]  = 0
    sumstats[final_shape_col] = 1
    sumstats[final_rsq_col] = 0.0

    if len(region_ref)==1:
        if lead_id is not None:
            sumstats.loc[lead_id, final_shape_col] +=1 

    for i in range(len(region_ref)):
        ld_single = "LD_{}".format(i)
        current_rsq = "RSQ_{}".format(i)
        a_ngt_b = sumstats[final_rsq_col] < sumstats[current_rsq]
        #set levels with interval=100
        sumstats.loc[a_ngt_b, final_ld_col] = 100 * (i+1) + sumstats.loc[a_ngt_b, ld_single]
        sumstats.loc[a_ngt_b, final_rsq_col] = sumstats.loc[a_ngt_b, current_rsq]
        sumstats.loc[a_ngt_b, final_shape_col] = i + 1
    
    sumstats = sumstats.dropna(subset=[pos,nea,ea])
    ####################################################################################################
    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return sumstats

# -############################################################################################################################################################################

def prepare_vcf_context(vcf_path, vcf_chr_dict=None, log=Log(), verbose=True):
    from shutil import which
    from gwaslab.io.io_vcf import auto_check_vcf_chr_dict
    tabix = which("tabix")
    log.write(" -tabix will be used: {}".format(tabix), verbose=verbose)
    vcf_chr_dict = auto_check_vcf_chr_dict(vcf_path, vcf_chr_dict, verbose, log)
    return vcf_chr_dict, tabix

def process_ld(sumstats, 
               ld_path, 
               ld_map_path,
               region,
               region_ref, 
               log, 
               verbose, 
               pos ,
               nea,
               ea, 
               region_ld_threshold, 
               ld_fmt = "npz",
               ld_if_square =  False,
               ld_if_add_T = False,
               ld_map_rename_dic = None,
               ld_map_kwargs = None):
    from gwaslab.io.io_load_ld import process_ld as _io_process_ld
    return _io_process_ld(
        sumstats=sumstats,
        ld_path=ld_path,
        ld_map_path=ld_map_path,
        region=region,
        region_ref=region_ref,
        log=log,
        verbose=verbose,
        pos=pos,
        nea=nea,
        ea=ea,
        region_ld_threshold=region_ld_threshold,
        ld_fmt=ld_fmt,
        ld_if_square=ld_if_square,
        ld_if_add_T=ld_if_add_T,
        ld_map_rename_dic=ld_map_rename_dic,
        ld_map_kwargs=ld_map_kwargs,
    )

def process_gtf(gtf_path,
                region,
                region_flank_factor,
                build,
                region_protein_coding,
                gtf_chr_dict,
                gtf_gene_name,
                verbose=True,
                log=Log()):
    #loading
    log.write(f" -Processing GTF from: {gtf_path}", verbose=verbose)
    
    # chr to string datatype using gtf_chr_dict
    to_query_chrom = gtf_chr_dict[region[0]]

    # loading gtf data
    if gtf_path =="default" or gtf_path =="ensembl":
        log.write(f"  -Loading Ensembl GTF for chromosome {to_query_chrom}, build {build}", verbose=verbose)
        gtf, actual_gtf_path = get_gtf(chrom=to_query_chrom, build=build, source="ensembl", if_return_path=True)
        if actual_gtf_path:
            log.write(f"  -Resolved GTF file path: {actual_gtf_path}", verbose=verbose)
    
    elif gtf_path =="refseq":
        log.write(f"  -Loading RefSeq GTF for chromosome {to_query_chrom}, build {build}", verbose=verbose)
        gtf, actual_gtf_path = get_gtf(chrom=to_query_chrom, build=build, source="refseq", if_return_path=True)
        if actual_gtf_path:
            log.write(f"  -Resolved GTF file path: {actual_gtf_path}", verbose=verbose)
    
    else:
        # if user-provided gtf
        log.write(f"  -Loading user-provided GTF file: {gtf_path}", verbose=verbose)
        #gtf = pd.read_csv(gtf_path,sep="\t",header=None, comment="#", low_memory=False,dtype={0:"string"})
        # Filter by chromosome early for speed
        gtf = read_gtf(gtf_path, chrom=gtf_chr_dict[region[0]])

    # filter in region
    genes_1mb = gtf.loc[(gtf["seqname"]==to_query_chrom)&(gtf["start"]<region[2])&(gtf["end"]>region[1]),:].copy()
    
    # extract biotype
    #genes_1mb.loc[:,"gene_biotype"] = genes_1mb[8].str.extract(r'gene_biotype "([\w\.\_-]+)"')
    
    # extract gene name
    if gtf_gene_name is None:
        if gtf_path=="refseq":
            #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"').astype("string")
            genes_1mb["name"] = genes_1mb["gene_id"]
        elif gtf_path =="default" or gtf_path =="ensembl":
            #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_name "([\w\.-]+)"').astype("string")
            genes_1mb["name"] = genes_1mb["gene_name"]
        else:
            if "gene_name" in gtf.columns:
                genes_1mb["name"] = genes_1mb["gene_name"]
            else:
                #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"').astype("string")
                genes_1mb["name"] = genes_1mb["gene_id"]
    else:
        #pattern = r'{} "([\w\.-]+)"'.format(gtf_gene_name)
        #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(pattern).astype("string")
        genes_1mb["name"] = genes_1mb[gtf_gene_name]

    # extract gene
    #genes_1mb.loc[:,"gene"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"')
    genes_1mb["gene"] = genes_1mb["gene_id"]

    # filter out genes with no name
    n_genes_before = len(genes_1mb)
    genes_1mb = genes_1mb[genes_1mb["name"].notna()]
    n_genes_after = len(genes_1mb)
    if n_genes_before != n_genes_after:
        log.write(f"  -Filtered out {n_genes_before - n_genes_after} genes with missing names ({n_genes_before} -> {n_genes_after})", verbose=verbose)

    # extract protein coding gene
    if region_protein_coding is True:
        #genes_1mb  =  genes_1mb.loc[genes_1mb["gene_biotype"]=="protein_coding",:]
        if "gene_biotype" in genes_1mb.columns:
            pc_genes_1mb_list = genes_1mb.loc[(genes_1mb["feature"]=="gene")& (genes_1mb["gene_biotype"]=="protein_coding") & (genes_1mb["name"]!=""),"name"].values
            genes_1mb = genes_1mb.loc[(genes_1mb["feature"].isin(["exon","gene"])) & (genes_1mb["name"].isin(pc_genes_1mb_list)),:]
        else:
            log.write("  -Warning: 'gene_biotype' column not found in GTF file. Cannot filter for protein coding genes. Showing all genes.", verbose=verbose)
    # extract exon
    exons = genes_1mb.loc[genes_1mb["feature"]=="exon",:].copy()
    
    #uniq genes
    ## get all record with 2nd column == gene
    #uniq_gene_region = genes_1mb.loc[genes_1mb[2]=="gene",:].copy()
    uniq_gene_region = genes_1mb.loc[genes_1mb["feature"]=="gene",:].copy()

    ## extract region + flank
    flank = region_flank_factor * (region[2] - region[1])
    
    ## get left and right boundary
    #uniq_gene_region["left"] = uniq_gene_region[3]-flank
    #uniq_gene_region["right"] = uniq_gene_region[4]+flank
    #
    uniq_gene_region["left"] = uniq_gene_region["start"]-flank
    uniq_gene_region["right"] = uniq_gene_region["end"]+flank

    # arrange gene track
    stack_dic = assign_stack(uniq_gene_region.sort_values(["start"]).loc[:,["name","left","right"]])  

    # map gene to stack and add stack column : minus stack
    uniq_gene_region["stack"] = -uniq_gene_region["name"].map(stack_dic)
    exons.loc[:,"stack"] = -exons.loc[:,"name"].map(stack_dic)

    # return uniq_gene_region (gene records with left and right boundary)
    # return exon records with stack number
    return uniq_gene_region, exons


# -############################################################################################################################################################################
def assign_stack(uniq_gene_region):

    stacks=[] ## stack : gene track
    stack_dic={} # mapping gene name to stack
    
    for index,row in uniq_gene_region.iterrows():
        if len(stacks)==0:
            # add first entry 
            stacks.append([(row["left"],row["right"])])
            stack_dic[row["name"]] = 0
        else:
            for i in range(len(stacks)):
                for j in range(len(stacks[i])):
                    # if overlap
                    if (row["left"]>stacks[i][j][0] and row["left"]<stacks[i][j][1]) or (row["right"]>stacks[i][j][0] and row["right"]<stacks[i][j][1]):
                        # if not last stack : break
                        if i<len(stacks)-1:
                            break
                        # if last stack : add a new stack
                        else:
                            stacks.append([(row["left"],row["right"])])
                            stack_dic[row["name"]] = i+1
                            break
                    # if no overlap       
                    else:
                        # not last in a stack
                        if j<len(stacks[i])-1:
                            #if in the middle
                            if row["left"]>stacks[i][j][1] and row["right"]<stacks[i][j+1][0]:
                                stacks[i].insert(j+1,(row["left"],row["right"]))
                                stack_dic[row["name"]] = i
                                break
                        # last one in a stack
                        elif row["left"]>stacks[i][j][1]:
                            stacks[i].append((row["left"],row["right"]))
                            stack_dic[row["name"]] = i
                            break
                if row["name"] in stack_dic.keys():
                    break         
    return stack_dic

def closest_gene(x,data,chrom="CHR",pos="POS",maxiter=20000,step=50):
    gene = data.gene_names_at_locus(contig=x[chrom], position=x[pos])
    if len(gene)==0:
        i=0
        while i<maxiter:
            distance = i*step
            gene_u = data.gene_names_at_locus(contig=x[chrom], position=x[pos]-distance)
            gene_d = data.gene_names_at_locus(contig=x[chrom], position=x[pos]+distance)
            if len(gene_u)>0 and len(gene_d)>0:
                for j in range(0,step,1):
                    distance = (i-1)*step
                    gene_u = data.gene_names_at_locus(contig=x[chrom], position=x[pos]-distance-j)
                    gene_d = data.gene_names_at_locus(contig=x[chrom], position=x[pos]+distance+j)
                    if len(gene_u)>0:
                        return -distance,",".join(gene_u)
                    else:
                        return distance,",".join(gene_d)
            elif len(gene_u)>0:
                return -distance,",".join(gene_u)
            elif len(gene_d)>0:
                return +distance,",".join(gene_d)
            else:
                i+=1
        return distance,"intergenic"
    else:
        return 0,",".join(gene)
