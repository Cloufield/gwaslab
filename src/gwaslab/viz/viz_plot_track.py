"""
General track plotting functions for GTF, BED, bigWig, and bigBed files.

This module provides functions for plotting genomic tracks from GTF (Gene Transfer Format),
BED (Browser Extensible Data), bigWig (continuous signal), or bigBed (binary BED) files,
with support for stacking overlapping features and customizable styling.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict, Union, Any
from adjustText import adjust_text
from gwaslab.io.io_gtf import read_gtf
from gwaslab.io.io_ucsc_bed import read_bed, bed_to_sumstats_coordinates
from gwaslab.io.io_bigwig_bigbed import (
    read_bigwig, read_bigwig_intervals, read_bigbed,
    PYBigWig_AVAILABLE, BIGWIG_SUFFIXES, BIGBED_SUFFIXES
)
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.info.g_Log import Log


def plot_track(
    track_path: str,
    region: Tuple[int, int, int],
    ax: Optional[plt.Axes] = None,
    fig: Optional[plt.Figure] = None,
    track_start_i: float = 0.0,
    region_flank_factor: float = 0.05,
    file_format: Optional[str] = None,
    chr_dict: Optional[Dict] = None,
    name_column: Optional[str] = None,
    feature_type: Optional[str] = None,
    color: str = "#020080",
    highlight_positions: Optional[List[float]] = None,
    highlight_color: str = "#FF0000",
    track_font_family: str = "Arial",
    taf: List[float] = [4, 0, 0.95, 1, 1],
    min_track_height: int = 4,
    figsize: Tuple[float, float] = (10, 2),
    region_step: int = 21,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    verbose: bool = True,
    log: Log = Log()
) -> Tuple[plt.Axes, List]:
    """
    Plot a genomic track from GTF, BED, bigWig, or bigBed file.
    
    This function creates a track visualization showing genomic features (genes, exons,
    or custom features) from GTF or BED files, or continuous signal from bigWig files.
    Features are automatically stacked to avoid overlaps, and labels are positioned intelligently.
    bigWig tracks are plotted as continuous line/area plots.
    
    Parameters
    ----------
    track_path : str
        Path to GTF, BED, bigWig, or bigBed file. Supports compressed files (.gtf.gz, .bed.gz).
        File format is auto-detected from extension if file_format is None.
    region : tuple of (int, int, int)
        Genomic region as (chromosome, start, end) in 1-based coordinates.
        Example: (1, 1000000, 2000000) for chr1:1000000-2000000
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. If None, a new figure and axes will be created.
    fig : matplotlib.figure.Figure, optional
        Figure object for font size calculations. If None, a new figure will be created.
        If ax is provided but fig is None, fig will be extracted from ax.
    track_start_i : float, default=0.0
        X-axis offset for track positioning. Used to align track with other plots.
    region_flank_factor : float, default=0.05
        Flanking factor for extending region boundaries when loading features.
        Features within flanked region are shown.
    file_format : str, optional
        File format: 'gtf', 'gff', 'bed', 'bigwig', 'bw', or 'bigbed', 'bb'.
        If None, auto-detected from file extension.
    chr_dict : dict, optional
        Dictionary mapping sumstats chromosome (int) to file chromosome (str).
        If None, uses get_number_to_chr().
    name_column : str, optional
        Column name to use for feature labels. 
        For GTF: 'gene_name', 'gene_id', or custom attribute name.
        For BED/bigBed: 'name' (if available in BED4+ format).
        If None, auto-detects: 'gene_name' for GTF, 'name' for BED/bigBed.
        Ignored for bigWig files.
    feature_type : str, optional
        For GTF files: filter by feature type (e.g., 'gene', 'exon', 'transcript').
        If None, uses 'gene' for GTF files. Ignored for BED/bigBed/bigWig files.
    color : str, default="#020080"
        Color for track features (hex color code). For bigWig, this is the line/area color.
    highlight_positions : list of float, optional
        List of x-axis positions to highlight (e.g., lead SNP positions).
        Features overlapping these positions will be colored with highlight_color.
        For bigWig, highlights are shown as vertical lines.
    highlight_color : str, default="#FF0000"
        Color for highlighted features (hex color code).
    track_font_family : str, default="Arial"
        Font family for track labels.
    taf : list of float, default=[4, 0, 0.95, 1, 1]
        Track annotation formatting parameters:
        [min_stacks, y_offset, font_ratio, linewidth_ratio, text_y_offset]
    min_track_height : int, default=4
        Minimum number of track stacks to reserve. For bigWig, this sets the y-axis range.
    figsize : tuple of (float, float), default=(10, 2)
        Figure size (width, height) in inches. Only used when creating a new figure.
        This will be merged into fig_kwargs if fig_kwargs is provided.
    region_step : int, default=21
        Number of x-axis ticks to display. Used for formatting position labels in MB.
    fig_kwargs : dict, optional
        Additional keyword arguments for matplotlib figure creation (e.g., {'dpi': 200, 'facecolor': 'white'}).
        If None, defaults to empty dict. figsize will be merged into this if provided.
    verbose : bool, default=True
        Whether to show progress messages.
    log : Log, default=Log()
        Logger instance for messages.
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object with track plotted.
    texts_to_adjust : list
        List of text objects that may need adjustment (for use with adjustText).
        Empty list for bigWig tracks.
    
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from gwaslab.viz.viz_plot_track import plot_track
    >>> 
    >>> # Simple usage - creates figure and axes automatically
    >>> region = (1, 1000000, 2000000)  # chr1:1Mb-2Mb
    >>> ax, texts = plot_track("genes.gtf", region)
    >>> plt.show()
    >>> 
    >>> # Plot bigWig signal track
    >>> ax, texts = plot_track("signal.bw", region, file_format="bigwig")
    >>> plt.show()
    >>> 
    >>> # Or provide your own figure/axes
    >>> fig, ax = plt.subplots(figsize=(10, 2))
    >>> ax, texts = plot_track("genes.gtf", region, ax=ax, fig=fig)
    >>> plt.show()
    
    Notes
    -----
    - GTF files use 1-based, inclusive coordinates [start, end]
    - BED/bigBed files use 0-based, half-open coordinates [chromStart, chromEnd)
    - bigWig files use 0-based, half-open coordinates [start, end)
    - The function automatically converts BED/bigBed/bigWig coordinates to 1-based for plotting
    - Features are automatically stacked to avoid overlaps (GTF/BED/bigBed)
    - bigWig tracks are plotted as continuous line/area plots
    - Labels are positioned on left, middle, or right based on feature position
    """
    # Initialize fig_kwargs
    if fig_kwargs is None:
        fig_kwargs = {}
    
    # Merge figsize into fig_kwargs if not already specified
    if "figsize" not in fig_kwargs:
        fig_kwargs["figsize"] = figsize
    
    # Create figure and axes if not provided
    if fig is None:
        if ax is None:
            fig, ax = plt.subplots(**fig_kwargs)
        else:
            fig = ax.figure
    elif ax is None:
        ax = fig.gca()
        if ax is None:
            ax = fig.add_subplot(111)
    
    # Auto-detect file format if not specified
    if file_format is None:
        track_path_lower = track_path.lower()
        if track_path_lower.endswith(('.gtf', '.gtf.gz', '.gff', '.gff.gz')):
            file_format = 'gtf'
        elif track_path_lower.endswith(('.bed', '.bed.gz')):
            file_format = 'bed'
        elif track_path_lower.endswith(BIGWIG_SUFFIXES):
            file_format = 'bigwig'
        elif track_path_lower.endswith(BIGBED_SUFFIXES):
            file_format = 'bigbed'
        else:
            log.warning(f"Could not auto-detect file format from {track_path}. Assuming GTF.")
            file_format = 'gtf'
    
    # Set default chromosome dictionary
    if chr_dict is None:
        chr_dict = get_number_to_chr()
    
    # Process track data based on file format
    if file_format.lower() in ['gtf', 'gff']:
        features_df, subfeatures_df = _process_gtf_track(
            gtf_path=track_path,
            region=region,
            region_flank_factor=region_flank_factor,
            chr_dict=chr_dict,
            name_column=name_column,
            feature_type=feature_type,
            verbose=verbose,
            log=log
        )
    elif file_format.lower() == 'bed':
        features_df, subfeatures_df = _process_bed_track(
            bed_path=track_path,
            region=region,
            region_flank_factor=region_flank_factor,
            chr_dict=chr_dict,
            name_column=name_column,
            verbose=verbose,
            log=log
        )
    elif file_format.lower() in ['bigwig', 'bw']:
        # bigWig is handled differently - it's a continuous signal
        return _plot_bigwig_track(
            bw_path=track_path,
            region=region,
            ax=ax,
            fig=fig,
            track_start_i=track_start_i,
            region_flank_factor=region_flank_factor,
            chr_dict=chr_dict,
            color=color,
            highlight_positions=highlight_positions,
            highlight_color=highlight_color,
            min_track_height=min_track_height,
            figsize=figsize,
            region_step=region_step,
            fig_kwargs=fig_kwargs,
            track_font_family=track_font_family,
            verbose=verbose,
            log=log
        )
    elif file_format.lower() in ['bigbed', 'bb']:
        features_df, subfeatures_df = _process_bigbed_track(
            bb_path=track_path,
            region=region,
            region_flank_factor=region_flank_factor,
            chr_dict=chr_dict,
            name_column=name_column,
            verbose=verbose,
            log=log
        )
    else:
        raise ValueError(f"Unsupported file format: {file_format}. Use 'gtf', 'bed', 'bigwig', or 'bigbed'.")
    
    if len(features_df) == 0:
        log.write(" -No features found in specified region.", verbose=verbose)
        ax.set_ylim((-min_track_height*2-taf[1]*2, 2+taf[1]*2))
        ax.set_yticks([])
        return ax, []
    
    # Calculate stacking
    stack_dict = _assign_stack(features_df.sort_values(["start"]).loc[:, ["name", "left", "right"]])
    features_df["stack"] = -features_df["name"].map(stack_dict)
    if len(subfeatures_df) > 0:
        subfeatures_df["stack"] = -subfeatures_df["name"].map(stack_dict)
    
    # Set up track layout
    n_stacks = features_df["stack"].nunique()
    stack_num_to_plot = max(taf[0], n_stacks, min_track_height)
    ax.set_ylim((-stack_num_to_plot*2-taf[1]*2, 2+taf[1]*2))
    ax.set_yticks([])
    
    # Calculate font and line sizes
    point_per_pixels = 72 / fig.dpi
    pixels_per_point = fig.dpi / 72
    pixels_per_track = np.abs(ax.transData.transform([0, 0])[1] - ax.transData.transform([0, 1])[1])
    font_size_in_pixels = taf[2] * pixels_per_track
    font_size_in_points = font_size_in_pixels * point_per_pixels
    linewidth_in_points_per_track = pixels_per_track * point_per_pixels
    
    log.write(f" -Plotting {len(features_df)} features..", verbose=verbose)
    
    # Track highlighted features
    highlighted_names = []
    highlighted_lefts = []
    highlighted_rights = []
    texts_to_adjust_middle = []
    
    # Plot main features
    for index, row in features_df.iterrows():
        feature_color = color
        
        # Check if feature should be highlighted
        if highlight_positions is not None:
            for highlight_pos in highlight_positions:
                if highlight_pos > track_start_i + row["start"] and highlight_pos < track_start_i + row["end"]:
                    feature_color = highlight_color
                    highlighted_names.append(row["name"])
                    highlighted_lefts.append(track_start_i + row["start"])
                    highlighted_rights.append(track_start_i + row["end"])
                    break
        
        # Plot feature line
        gene_line_width = max(linewidth_in_points_per_track / 10, 2 / pixels_per_point)
        
        # Determine strand direction for label
        strand = row.get("strand", "+")
        if isinstance(strand, str) and len(strand) > 0:
            if strand[0] == "+":
                feature_anno = row["name"] + "->"
            elif strand[0] == "-":
                feature_anno = "<-" + row["name"]
            else:
                feature_anno = row["name"]
        else:
            feature_anno = row["name"]
        
        ax.plot(
            (track_start_i + row["start"], track_start_i + row["end"]),
            (row["stack"] * 2, row["stack"] * 2),
            color=feature_color,
            linewidth=gene_line_width,
            solid_capstyle="butt"
        )
        
        # Plot feature name
        if row["end"] >= region[2]:
            # Right side
            ax.text(
                x=track_start_i + region[2],
                y=row["stack"] * 2 + taf[4],
                s=feature_anno,
                ha="right",
                va="center",
                color="black",
                style='italic',
                size=font_size_in_points,
                family=track_font_family
            )
        elif row["start"] <= region[1]:
            # Left side
            ax.text(
                x=track_start_i + region[1],
                y=row["stack"] * 2 + taf[4],
                s=feature_anno,
                ha="left",
                va="center",
                color="black",
                style='italic',
                size=font_size_in_points,
                family=track_font_family
            )
        else:
            # Middle
            text_obj = ax.text(
                x=(track_start_i + row["start"] + track_start_i + row["end"]) / 2,
                y=row["stack"] * 2 + taf[4],
                s=feature_anno,
                ha="center",
                va="center",
                color="black",
                style='italic',
                size=font_size_in_points,
                family=track_font_family
            )
            texts_to_adjust_middle.append(text_obj)
    
    # Plot subfeatures (e.g., exons for GTF)
    if len(subfeatures_df) > 0:
        log.write(f" -Plotting {len(subfeatures_df)} subfeatures..", verbose=verbose)
        for index, row in subfeatures_df.iterrows():
            subfeature_color = color
            
            # Check if subfeature should be highlighted
            for highlighted_name, highlighted_left, highlighted_right in zip(
                highlighted_names, highlighted_lefts, highlighted_rights
            ):
                if not pd.isnull(row.get("name")):
                    if row["name"] == highlighted_name:
                        subfeature_color = highlight_color
                        break
                elif track_start_i + row["start"] > highlighted_left and track_start_i + row["end"] < highlighted_right:
                    subfeature_color = highlight_color
                    break
            
            # Plot subfeature (thicker line)
            exon_line_width = max(linewidth_in_points_per_track * taf[3], 8 / pixels_per_point)
            ax.plot(
                (track_start_i + row["start"], track_start_i + row["end"]),
                (row["stack"] * 2, row["stack"] * 2),
                linewidth=exon_line_width,
                color=subfeature_color,
                solid_capstyle="butt"
            )
    
    log.write(" -Finished plotting track..", verbose=verbose)
    
    # Set x-axis ticks and labels (similar to regional plot)
    if region is not None:
        region_ticks = list(map('{:.3f}'.format, np.linspace(region[1], region[2], num=region_step).astype("int")/1000000))
        ax.set_xticks(np.linspace(track_start_i + region[1], track_start_i + region[2], num=region_step))
        ax.set_xticklabels(region_ticks, rotation=45)
        ax.set_xlim([track_start_i + region[1], track_start_i + region[2]])
        # Set x-axis label (similar to manhattan plot)
        xlabel = "Chromosome " + str(region[0]) + " (MB)"
        # Use font_size_in_points for consistency with other text elements
        ax.set_xlabel(xlabel, fontsize=font_size_in_points, family=track_font_family)
    
    # Adjust text positions to prevent overlapping (similar to regional plot)
    if len(texts_to_adjust_middle) > 0:
        adjust_text(texts_to_adjust_middle,
                    autoalign=False,
                    only_move={'points': 'x', 'text': 'x', 'objects': 'x'},
                    ax=ax,
                    precision=0,
                    force_text=(0.1, 0),
                    expand_text=(1, 1),
                    expand_objects=(1, 1),
                    expand_points=(1, 1),
                    va="center",
                    ha='center',
                    avoid_points=False,
                    lim=1000)
    
    return ax, texts_to_adjust_middle


def _process_gtf_track(
    gtf_path: str,
    region: Tuple[int, int, int],
    region_flank_factor: float,
    chr_dict: Dict,
    name_column: Optional[str],
    feature_type: Optional[str],
    verbose: bool,
    log: Log
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process GTF file for track plotting."""
    log.write(f" -Processing GTF from: {gtf_path}", verbose=verbose)
    
    # Get chromosome string for querying
    to_query_chrom = chr_dict[region[0]]
    
    # Load GTF data
    log.write(f"  -Loading GTF for chromosome {to_query_chrom}", verbose=verbose)
    gtf = read_gtf(gtf_path, chrom=to_query_chrom)
    
    if len(gtf) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Filter by region (with flanking)
    flank = region_flank_factor * (region[2] - region[1])
    gtf_region = gtf.loc[
        (gtf["seqname"] == to_query_chrom) &
        (gtf["start"] < region[2] + flank) &
        (gtf["end"] > region[1] - flank)
    ].copy()
    
    if len(gtf_region) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Determine feature type
    if feature_type is None:
        feature_type = "gene"
    
    # Determine name column
    if name_column is None:
        if "gene_name" in gtf_region.columns:
            name_column = "gene_name"
        elif "gene_id" in gtf_region.columns:
            name_column = "gene_id"
        else:
            name_column = "gene_id"
            log.warning("  -Neither 'gene_name' nor 'gene_id' found. Using 'gene_id'.")
    
    # Extract name
    if name_column in gtf_region.columns:
        gtf_region["name"] = gtf_region[name_column].fillna("").astype(str)
        # Use index for empty names
        empty_mask = gtf_region["name"] == ""
        gtf_region.loc[empty_mask, "name"] = "feature_" + gtf_region.loc[empty_mask].index.astype(str)
    else:
        gtf_region["name"] = "feature_" + gtf_region.index.astype(str)
        log.warning(f"  -Name column '{name_column}' not found. Using index-based names.")
    
    # Filter by feature type
    if feature_type is not None:
        gtf_region = gtf_region.loc[gtf_region["feature"] == feature_type, :].copy()
    
    if len(gtf_region) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Prepare main features
    features_df = gtf_region.copy()
    features_df["left"] = features_df["start"] - flank
    features_df["right"] = features_df["end"] + flank
    
    # Prepare subfeatures (exons if feature_type is gene)
    subfeatures_df = pd.DataFrame()
    if feature_type == "gene":
        # Get exons for the same genes
        gtf_exons = gtf.loc[
            (gtf["seqname"] == to_query_chrom) &
            (gtf["feature"] == "exon") &
            (gtf["start"] < region[2] + flank) &
            (gtf["end"] > region[1] - flank)
        ].copy()
        
        if len(gtf_exons) > 0:
            # Map exons to gene names using the same logic as main features
            if name_column is None:
                if "gene_name" in gtf_exons.columns:
                    exon_name_col = "gene_name"
                elif "gene_id" in gtf_exons.columns:
                    exon_name_col = "gene_id"
                else:
                    exon_name_col = "gene_id"
            else:
                exon_name_col = name_column
            
            if exon_name_col in gtf_exons.columns:
                gtf_exons["name"] = gtf_exons[exon_name_col].fillna("").astype(str)
                # Use index for empty names
                empty_mask = gtf_exons["name"] == ""
                gtf_exons.loc[empty_mask, "name"] = "feature_" + gtf_exons.loc[empty_mask].index.astype(str)
            elif "gene_id" in gtf_exons.columns:
                gtf_exons["name"] = gtf_exons["gene_id"].fillna("").astype(str)
                empty_mask = gtf_exons["name"] == ""
                gtf_exons.loc[empty_mask, "name"] = "feature_" + gtf_exons.loc[empty_mask].index.astype(str)
            else:
                gtf_exons["name"] = "feature_" + gtf_exons.index.astype(str)
            
            # Only keep exons for genes in features_df
            if len(features_df) > 0:
                valid_names = set(features_df["name"].unique())
                gtf_exons = gtf_exons.loc[gtf_exons["name"].isin(valid_names), :].copy()
            
            subfeatures_df = gtf_exons.copy()
    
    return features_df, subfeatures_df


def _process_bed_track(
    bed_path: str,
    region: Tuple[int, int, int],
    region_flank_factor: float,
    chr_dict: Dict,
    name_column: Optional[str],
    verbose: bool,
    log: Log
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process BED file for track plotting."""
    log.write(f" -Processing BED from: {bed_path}", verbose=verbose)
    
    # Load BED file
    log.write(f"  -Loading BED file", verbose=verbose)
    bed_df = read_bed(bed_path, verbose=verbose, log=log)
    
    if len(bed_df) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Convert BED coordinates to 1-based
    mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    bed_df_converted = bed_to_sumstats_coordinates(
        bed_df,
        mapper=mapper,
        chrom_col='chrom',
        start_col='chromStart',
        end_col='chromEnd',
        verbose=verbose,
        log=log
    )
    
    # Get chromosome string for querying
    to_query_chrom = chr_dict[region[0]]
    
    # Filter by region (with flanking)
    flank = region_flank_factor * (region[2] - region[1])
    bed_region = bed_df_converted.loc[
        (bed_df_converted["CHR"] == region[0]) &
        (bed_df_converted["START"] < region[2] + flank) &
        (bed_df_converted["END"] > region[1] - flank)
    ].copy()
    
    if len(bed_region) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Determine name column
    if name_column is None:
        if "name" in bed_region.columns:
            name_column = "name"
        else:
            name_column = None
            log.warning("  -No 'name' column found in BED file. Using index as names.")
    
    # Extract name
    if name_column and name_column in bed_region.columns:
        bed_region["name"] = bed_region[name_column].fillna("").astype(str)
        # Use index for empty names
        empty_mask = bed_region["name"] == ""
        bed_region.loc[empty_mask, "name"] = "feature_" + bed_region.loc[empty_mask].index.astype(str)
    else:
        bed_region["name"] = "feature_" + bed_region.index.astype(str)
    
    # Prepare features (BED doesn't have subfeatures like exons)
    features_df = bed_region.copy()
    features_df["start"] = features_df["START"]
    features_df["end"] = features_df["END"]
    features_df["left"] = features_df["start"] - flank
    features_df["right"] = features_df["end"] + flank
    features_df["strand"] = bed_region.get("strand", "+")
    
    # No subfeatures for BED
    subfeatures_df = pd.DataFrame()
    
    return features_df, subfeatures_df


def _assign_stack(features_df: pd.DataFrame) -> Dict:
    """
    Assign stacking positions to features to avoid overlaps.
    
    This is a simplified version of the assign_stack function from viz_plot_regional2.py.
    Features are assigned to stacks (tracks) such that overlapping features are on different stacks.
    """
    stacks = []  # stack: list of (left, right) tuples
    stack_dict = {}  # mapping feature name to stack number
    
    for index, row in features_df.iterrows():
        if len(stacks) == 0:
            # Add first entry
            stacks.append([(row["left"], row["right"])])
            stack_dict[row["name"]] = 0
        else:
            assigned = False
            for i in range(len(stacks)):
                for j in range(len(stacks[i])):
                    # Check if overlap with existing feature
                    # Overlap occurs if: 
                    # - new left is within existing range, OR
                    # - new right is within existing range, OR
                    # - new feature completely contains existing feature
                    if (row["left"] > stacks[i][j][0] and row["left"] < stacks[i][j][1]) or \
                       (row["right"] > stacks[i][j][0] and row["right"] < stacks[i][j][1]) or \
                       (row["left"] <= stacks[i][j][0] and row["right"] >= stacks[i][j][1]):
                        # If not last stack, break and try next stack
                        if i < len(stacks) - 1:
                            break
                        # If last stack, add a new stack
                        else:
                            stacks.append([(row["left"], row["right"])])
                            stack_dict[row["name"]] = i + 1
                            assigned = True
                            break
                    # If no overlap with this feature
                    else:
                        # Not last in a stack
                        if j < len(stacks[i]) - 1:
                            # If in the middle (between two non-overlapping features)
                            if row["left"] > stacks[i][j][1] and row["right"] < stacks[i][j+1][0]:
                                stacks[i].insert(j + 1, (row["left"], row["right"]))
                                stack_dict[row["name"]] = i
                                assigned = True
                                break
                        # Last one in a stack
                        elif row["left"] > stacks[i][j][1]:
                            stacks[i].append((row["left"], row["right"]))
                            stack_dict[row["name"]] = i
                            assigned = True
                            break
                if assigned:
                    break
            
            # If not assigned to any existing stack, create a new one
            if not assigned:
                stacks.append([(row["left"], row["right"])])
                stack_dict[row["name"]] = len(stacks) - 1
    
    return stack_dict


def _plot_bigwig_track(
    bw_path: str,
    region: Tuple[int, int, int],
    ax: Optional[plt.Axes] = None,
    fig: Optional[plt.Figure] = None,
    track_start_i: float = 0.0,
    region_flank_factor: float = 0.05,
    chr_dict: Optional[Dict] = None,
    color: str = "#020080",
    highlight_positions: Optional[List[float]] = None,
    highlight_color: str = "#FF0000",
    min_track_height: int = 4,
    figsize: Tuple[float, float] = (10, 2),
    region_step: int = 21,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    track_font_family: str = "Arial",
    verbose: bool = True,
    log: Log = Log()
) -> Tuple[plt.Axes, List]:
    """
    Plot a bigWig track as a continuous signal line/area plot.
    
    Parameters
    ----------
    bw_path : str
        Path to bigWig file
    region : tuple of (int, int, int)
        Genomic region as (chromosome, start, end) in 1-based coordinates
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on
    fig : matplotlib.figure.Figure, optional
        Figure object
    track_start_i : float, default=0.0
        X-axis offset for track positioning
    region_flank_factor : float, default=0.05
        Flanking factor for extending region boundaries
    chr_dict : dict, optional
        Dictionary mapping sumstats chromosome (int) to file chromosome (str)
    color : str, default="#020080"
        Color for the signal line/area
    highlight_positions : list of float, optional
        List of x-axis positions to highlight with vertical lines
    highlight_color : str, default="#FF0000"
        Color for highlight lines
    min_track_height : int, default=4
        Sets the y-axis range (not used for actual height calculation)
    figsize : tuple of (float, float), default=(10, 2)
        Figure size
    verbose : bool, default=True
        Whether to show progress messages
    log : Log, default=Log()
        Logger instance
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object with track plotted
    texts_to_adjust : list
        Empty list (no text labels for bigWig)
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    # Initialize fig_kwargs
    if fig_kwargs is None:
        fig_kwargs = {}
    
    # Merge figsize into fig_kwargs if not already specified
    if "figsize" not in fig_kwargs:
        fig_kwargs["figsize"] = figsize
    
    # Create figure and axes if not provided
    if fig is None:
        if ax is None:
            fig, ax = plt.subplots(**fig_kwargs)
        else:
            fig = ax.figure
    elif ax is None:
        ax = fig.gca()
        if ax is None:
            ax = fig.add_subplot(111)
    
    # Set default chromosome dictionary
    if chr_dict is None:
        chr_dict = get_number_to_chr()
    
    # Get chromosome string for querying
    to_query_chrom = chr_dict[region[0]]
    
    # Calculate flanked region
    flank = region_flank_factor * (region[2] - region[1])
    # bigWig uses 0-based coordinates, so convert from 1-based
    start_0based = max(0, region[1] - 1 - int(flank))
    end_0based = region[2] - 1 + int(flank)
    
    if log and verbose:
        log.write(f" -Reading bigWig from: {bw_path}", verbose=verbose)
        log.write(f"  -Querying {to_query_chrom}:{start_0based}-{end_0based} (0-based)", verbose=verbose)
    
    # Read bigWig values
    try:
        values = read_bigwig(
            bw_path=bw_path,
            chrom=to_query_chrom,
            start=start_0based,
            end=end_0based,
            as_dataframe=False,
            numpy=True,
            verbose=verbose,
            log=log
        )
    except Exception as e:
        log.warning(f"  -Error reading bigWig file: {e}")
        ax.set_ylim((0, 1))
        ax.set_yticks([])
        return ax, []
    
    if values is None or len(values) == 0:
        log.write("  -No data found in specified region.", verbose=verbose)
        ax.set_ylim((0, 1))
        ax.set_yticks([])
        return ax, []
    
    # Convert to numpy array and handle NaN values
    if isinstance(values, list):
        values = np.array(values)
    
    # Create x-axis positions (1-based coordinates for plotting)
    x_positions = np.arange(start_0based + 1, start_0based + 1 + len(values))
    
    # Filter to region of interest (remove flanking for display)
    region_mask = (x_positions >= region[1]) & (x_positions <= region[2])
    x_plot = x_positions[region_mask]
    y_plot = values[region_mask]
    
    # Remove NaN values for plotting
    valid_mask = ~np.isnan(y_plot)
    if not np.any(valid_mask):
        log.write("  -No valid values in specified region.", verbose=verbose)
        ax.set_ylim((0, 1))
        ax.set_yticks([])
        return ax, []
    
    x_plot = x_plot[valid_mask]
    y_plot = y_plot[valid_mask]
    
    # Adjust x positions relative to region start and track_start_i
    # x positions are in 1-based coordinates, convert to relative positions
    x_plot = track_start_i + (x_plot - region[1])
    
    # Plot as filled area
    ax.fill_between(x_plot, 0, y_plot, color=color, alpha=0.3, linewidth=0)
    
    # Plot as line
    ax.plot(x_plot, y_plot, color=color, linewidth=1.5)
    
    # Set y-axis limits
    y_min = np.nanmin(y_plot) if len(y_plot) > 0 else 0
    y_max = np.nanmax(y_plot) if len(y_plot) > 0 else 1
    if y_min == y_max:
        y_max = y_min + 1
    y_padding = (y_max - y_min) * 0.1
    ax.set_ylim((max(0, y_min - y_padding), y_max + y_padding))
    
    # Add highlight lines if specified
    if highlight_positions is not None:
        for highlight_pos in highlight_positions:
            # Convert highlight position to plot coordinates
            plot_pos = highlight_pos
            ax.axvline(plot_pos, color=highlight_color, linestyle='--', linewidth=1, alpha=0.7)
    
    ax.set_yticks([])
    ax.set_ylabel("Signal", rotation=90, va='center')
    
    # Set x-axis ticks and labels (similar to regional plot)
    if region is not None:
        region_ticks = list(map('{:.3f}'.format, np.linspace(region[1], region[2], num=region_step).astype("int")/1000000))
        ax.set_xticks(np.linspace(track_start_i + region[1], track_start_i + region[2], num=region_step))
        ax.set_xticklabels(region_ticks, rotation=45)
        ax.set_xlim([track_start_i + region[1], track_start_i + region[2]])
        # Set x-axis label (similar to manhattan plot)
        xlabel = "Chromosome " + str(region[0]) + " (MB)"
        # Use a reasonable fontsize for xlabel (similar to other plots, default 10)
        ax.set_xlabel(xlabel, fontsize=10, family=track_font_family)
    
    if log and verbose:
        log.write(f"  -Plotted {len(x_plot)} data points", verbose=verbose)
    
    return ax, []


def _process_bigbed_track(
    bb_path: str,
    region: Tuple[int, int, int],
    region_flank_factor: float,
    chr_dict: Dict,
    name_column: Optional[str],
    verbose: bool,
    log: Log
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process bigBed file for track plotting."""
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    log.write(f" -Processing bigBed from: {bb_path}", verbose=verbose)
    
    # Get chromosome string for querying
    to_query_chrom = chr_dict[region[0]]
    
    # Calculate flanked region
    flank = region_flank_factor * (region[2] - region[1])
    # bigBed uses 0-based coordinates, so convert from 1-based
    start_0based = max(0, region[1] - 1 - int(flank))
    end_0based = region[2] - 1 + int(flank)
    
    # Ensure start < end
    if start_0based >= end_0based:
        log.warning(f"  -Invalid coordinate range after conversion: start={start_0based}, end={end_0based}")
        return pd.DataFrame(), pd.DataFrame()
    
    # Read bigBed entries
    log.write(f"  -Loading bigBed for chromosome {to_query_chrom} (0-based: {start_0based}-{end_0based})", verbose=verbose)
    try:
        bb_df = read_bigbed(
            bb_path=bb_path,
            chrom=to_query_chrom,
            start=start_0based,
            end=end_0based,
            with_string=True,
            verbose=verbose,
            log=log
        )
    except Exception as e:
        log.warning(f"  -Error reading bigBed file: {e}")
        # Try to get more information about the chromosome
        if PYBigWig_AVAILABLE:
            try:
                # Use the already-imported read_bigbed function
                header_df = read_bigbed(bb_path, chrom=None, verbose=False, log=log)
                if len(header_df) > 0:
                    chrom_info = header_df[header_df["chrom"] == to_query_chrom]
                    if len(chrom_info) > 0:
                        chrom_size = chrom_info.iloc[0]["size"]
                        log.warning(f"  -Chromosome {to_query_chrom} size: {chrom_size}, requested range: {start_0based}-{end_0based}")
            except Exception as e2:
                # Silently ignore errors when trying to get header info
                pass
        return pd.DataFrame(), pd.DataFrame()
    
    if len(bb_df) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Convert 0-based coordinates to 1-based for plotting
    bb_df["START"] = bb_df["start"] + 1  # Convert to 1-based
    bb_df["END"] = bb_df["end"]  # End is already 1-based (exclusive in 0-based = inclusive in 1-based)
    
    # Filter by region (with flanking)
    bb_region = bb_df.loc[
        (bb_df["START"] < region[2] + flank) &
        (bb_df["END"] > region[1] - flank)
    ].copy()
    
    if len(bb_region) == 0:
        return pd.DataFrame(), pd.DataFrame()
    
    # Determine name column
    if name_column is None:
        if "rest" in bb_region.columns:
            # Try to extract name from rest column (BED format: name is often first field in rest)
            # For now, use index as name
            name_column = None
        elif "name" in bb_region.columns:
            name_column = "name"
        else:
            name_column = None
            log.warning("  -No 'name' column found in bigBed file. Using index as names.")
    
    # Extract name
    if name_column and name_column in bb_region.columns:
        bb_region["name"] = bb_region[name_column].fillna("").astype(str)
        empty_mask = bb_region["name"] == ""
        bb_region.loc[empty_mask, "name"] = "feature_" + bb_region.loc[empty_mask].index.astype(str)
    else:
        # Try to extract from rest column if available
        if "rest" in bb_region.columns:
            # Assume first field in rest is the name (BED format)
            bb_region["name"] = bb_region["rest"].apply(
                lambda x: str(x).split("\t")[0] if pd.notna(x) and str(x).strip() else ""
            )
            empty_mask = bb_region["name"] == ""
            bb_region.loc[empty_mask, "name"] = "feature_" + bb_region.loc[empty_mask].index.astype(str)
        else:
            bb_region["name"] = "feature_" + bb_region.index.astype(str)
    
    # Prepare features
    features_df = bb_region.copy()
    features_df["start"] = features_df["START"]
    features_df["end"] = features_df["END"]
    features_df["left"] = features_df["start"] - flank
    features_df["right"] = features_df["end"] + flank
    features_df["strand"] = "+"  # bigBed doesn't always have strand info
    
    # No subfeatures for bigBed
    subfeatures_df = pd.DataFrame()
    
    return features_df, subfeatures_df
