"""
Plot arcs for BEDPE contact pairs (e.g., Hi-C read-pairs).

This module provides functions for visualizing genomic contact pairs from BEDPE files
as arcs connecting the two ends of each pair.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.colors import Normalize
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from typing import Optional, Tuple, List, Dict, Union, Any, Sequence
from gwaslab.io.io_bedpe import read_bedpe
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.info.g_Log import Log


def _arc_cbar_locator(cax: plt.Axes):
    """Return the inset locator underlying a colorbar axes, if present."""
    locator = getattr(cax, "_axes_locator", None)
    if locator is not None and hasattr(locator, "_orig_locator"):
        return locator._orig_locator
    return None


def fit_arc_colorbar_inside_panel(
    ax: plt.Axes,
    cax: plt.Axes,
    pad_frac: float = 0.03,
    max_iter: int = 12,
    reset_anchor: bool = True,
) -> None:
    """
    Shrink the inset anchor box so the horizontal colorbar and its labels
    stay inside ``ax`` (works after plot_panels repositions the parent).
    """
    orig = _arc_cbar_locator(cax)
    if orig is None:
        return

    if reset_anchor:
        orig._bbox_to_anchor = Bbox.from_extents(0.02, 0.02, 0.98, 0.98)

    fig = ax.figure
    for _ in range(max_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        panel = ax.get_position()
        pad_x = pad_frac * panel.width
        pad_y = pad_frac * panel.height
        tight = cax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())

        top_ov = max(0.0, tight.y1 - (panel.y1 - pad_y))
        right_ov = max(0.0, tight.x1 - (panel.x1 - pad_x))
        left_ov = max(0.0, (panel.x0 + pad_x) - tight.x0)
        bottom_ov = max(0.0, (panel.y0 + pad_y) - tight.y0)
        if max(top_ov, right_ov, left_ov, bottom_ov) < 1e-5:
            return

        anchor = orig._bbox_to_anchor
        x0, y0, x1, y1 = anchor.x0, anchor.y0, anchor.x1, anchor.y1
        if top_ov:
            y1 -= top_ov / panel.height
        if bottom_ov:
            y0 += bottom_ov / panel.height
        if right_ov:
            x1 -= right_ov / panel.width
        if left_ov:
            x0 += left_ov / panel.width
        orig._bbox_to_anchor = Bbox.from_extents(x0, y0, x1, y1)


def refit_arc_colorbars(axes: Sequence[plt.Axes]) -> None:
    """Re-fit score colorbars after stacked-panel layout / x-axis alignment."""
    for ax in axes:
        cax = getattr(ax, "_gwaslab_arc_cbar_axes", None)
        if cax is not None:
            fit_arc_colorbar_inside_panel(ax, cax)


def plot_arc(
    bedpe_path: str,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    ax: Optional[plt.Axes] = None,
    fig: Optional[plt.Figure] = None,
    track_start_i: float = 0.0,
    region_flank_factor: float = 0.05,
    chr_dict: Optional[Dict] = None,
    color: str = "#020080",
    alpha: float = 0.3,
    linewidth: float = 0.5,
    arc_height: float = 0.1,
    max_arcs: Optional[int] = None,
    score_col: Optional[str] = None,
    score_threshold: Optional[float] = None,
    color_by_score: bool = False,
    cmap: str = "viridis",
    figsize: Tuple[float, float] = (10, 2),
    region_step: int = 21,
    track_font_family: str = "Arial",
    fig_kwargs: Optional[Dict[str, Any]] = None,
    verbose: bool = True,
    log: Log = Log()
) -> Tuple[plt.Axes, pd.DataFrame]:
    """
    Plot arcs connecting BEDPE contact pairs (cis contacts only).
    
    This function reads a BEDPE file and visualizes cis contact pairs (same chromosome)
    as arcs connecting the two genomic regions. Arcs are drawn above the chromosome
    using Bezier curves for smooth appearance.
    
    Parameters
    ----------
    bedpe_path : str
        Path to BEDPE file. Supports uncompressed (.bedpe) and gzipped (.bedpe.gz) files.
    region : tuple of (chromosome, start, end), optional
        Genomic region to plot. If None, plots all pairs.
        Example: (1, 1000000, 2000000) for chr1:1000000-2000000
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. If None, a new figure and axes will be created.
    fig : matplotlib.figure.Figure, optional
        Figure object. If None, a new figure will be created.
        If ax is provided but fig is None, fig will be extracted from ax.
    track_start_i : float, default=0.0
        X-axis offset for track positioning. Used to align arcs with other plots.
    region_flank_factor : float, default=0.05
        Flanking factor for extending region boundaries when filtering pairs.
    chr_dict : dict, optional
        Dictionary mapping sumstats chromosome (int) to file chromosome (str).
        If None, uses get_number_to_chr().
    color : str, default="#020080"
        Color for arcs (hex color code).
    alpha : float, default=0.3
        Transparency of arcs (0.0 to 1.0).
    linewidth : float, default=0.5
        Line width for arcs.
    arc_height : float, default=0.1
        Height of arcs as fraction of y-axis range. Larger values create taller arcs.
    max_arcs : int, optional
        Maximum number of arcs to plot. If None, plots all pairs.
        Useful for large files to avoid overcrowding.
    score_col : str, optional
        Column name for score/weight. If provided, can be used for filtering or coloring.
    score_threshold : float, optional
        Minimum score to plot. Only pairs with score >= threshold are shown.
    color_by_score : bool, default=False
        If True, color arcs by score using cmap. Requires score_col to be present.
    cmap : str, default="viridis"
        Colormap for score-based coloring. Only used if color_by_score=True.
    figsize : tuple of (float, float), default=(10, 2)
        Figure size (width, height) in inches. Only used when creating a new figure.
    region_step : int, default=21
        Number of x-axis ticks to display. Used for formatting position labels in MB.
    track_font_family : str, default="Arial"
        Font family for axis labels.
    fig_kwargs : dict, optional
        Additional keyword arguments for matplotlib figure creation.
    verbose : bool, default=True
        Whether to show progress messages.
    log : Log, default=Log()
        Logger instance for messages.
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object with arcs plotted.
    bedpe_df : pd.DataFrame
        DataFrame with the BEDPE pairs that were plotted.
    
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from gwaslab.viz.viz_plot_arc import plot_arc
    >>> 
    >>> # Simple usage - creates figure and axes automatically
    >>> region = (1, 1000000, 2000000)  # chr1:1Mb-2Mb
    >>> ax, df = plot_arc("contacts.bedpe.gz", region)
    >>> plt.show()
    >>> 
    >>> # Plot with custom styling
    >>> ax, df = plot_arc(
    ...     "contacts.bedpe.gz",
    ...     region,
    ...     color="#FF0000",
    ...     alpha=0.5,
    ...     linewidth=1.0,
    ...     arc_height=0.2
    ... )
    >>> plt.show()
    >>> 
    >>> # Color by score
    >>> ax, df = plot_arc(
    ...     "contacts.bedpe.gz",
    ...     region,
    ...     score_col="score",
    ...     color_by_score=True,
    ...     cmap="plasma"
    ... )
    >>> plt.show()
    
    Notes
    -----
    - BEDPE files use 0-based, half-open intervals [chromStart, chromEnd)
    - The function converts BEDPE coordinates to 1-based for plotting
    - Only cis contacts (same chromosome) are plotted; trans contacts are filtered out
    - Arcs are drawn using Bezier curves for smooth appearance
    - X-axis formatting matches plot_track: positions in MB with rotated labels
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
    
    # Set default chromosome dictionary
    if chr_dict is None:
        chr_dict = get_number_to_chr()
    
    # Read BEDPE file
    if log and verbose:
        log.write(f"Reading BEDPE file: {bedpe_path}", verbose=verbose)
    
    bedpe_df = read_bedpe(bedpe_path, verbose=verbose, log=log)
    
    if len(bedpe_df) == 0:
        if log:
            log.write("Warning: BEDPE file is empty.", verbose=verbose)
        return ax, bedpe_df
    
    # Filter by region if specified
    if region is not None:
        chrom, start, end = region
        # Convert chromosome to string format used in BEDPE
        # Try multiple formats for matching
        if isinstance(chrom, int):
            chrom_str = chr_dict.get(chrom, str(chrom))
        else:
            chrom_str = str(chrom)
        
        # Normalize chromosome names for matching (handle both "chr10" and "10")
        def normalize_chr(chr_val):
            """Normalize chromosome name for comparison."""
            if pd.isna(chr_val):
                return None
            chr_str = str(chr_val).strip()
            # Remove "chr" prefix if present
            if chr_str.lower().startswith('chr'):
                chr_str = chr_str[3:]
            return chr_str
        
        bedpe_chr1_norm = bedpe_df['chrom1'].apply(normalize_chr)
        bedpe_chr2_norm = bedpe_df['chrom2'].apply(normalize_chr)
        chrom_norm = normalize_chr(chrom_str)
        
        # Extend region by flank factor
        region_size = end - start
        flank_size = int(region_size * region_flank_factor)
        flanked_start = max(0, start - flank_size)
        flanked_end = end + flank_size
        
        # Filter pairs that overlap with the region (either end)
        mask = (
            ((bedpe_chr1_norm == chrom_norm) & 
             (bedpe_df['chromEnd1'] >= flanked_start) & 
             (bedpe_df['chromStart1'] <= flanked_end)) |
            ((bedpe_chr2_norm == chrom_norm) & 
             (bedpe_df['chromEnd2'] >= flanked_start) & 
             (bedpe_df['chromStart2'] <= flanked_end))
        )
        bedpe_df = bedpe_df[mask].copy()
        
        if len(bedpe_df) == 0:
            if log:
                log.write(f"No pairs found in region {chrom}:{start}-{end}", verbose=verbose)
            return ax, bedpe_df
    
    # Filter by score threshold if specified
    if score_threshold is not None and score_col and score_col in bedpe_df.columns:
        bedpe_df = bedpe_df[bedpe_df[score_col] >= score_threshold].copy()
    
    # Filter to only cis contacts (same chromosome)
    bedpe_df = bedpe_df[bedpe_df['chrom1'] == bedpe_df['chrom2']].copy()
    
    if len(bedpe_df) == 0:
        if log:
            log.write("Warning: No cis contacts found in BEDPE file.", verbose=verbose)
        return ax, bedpe_df
    
    # Calculate distances for normalization BEFORE limiting by max_arcs
    # This ensures normalization is based on ALL records in the region, not just the limited subset
    bedpe_df = bedpe_df.copy()
    bedpe_df['mid1'] = (bedpe_df['chromStart1'] + bedpe_df['chromEnd1']) / 2.0
    bedpe_df['mid2'] = (bedpe_df['chromStart2'] + bedpe_df['chromEnd2']) / 2.0
    bedpe_df['distance'] = np.abs(bedpe_df['mid2'] - bedpe_df['mid1'])
    
    # Calculate min/max distances for normalization (based on ALL records in region)
    n_total_in_region = len(bedpe_df)
    dist_min = bedpe_df['distance'].min()
    dist_max = bedpe_df['distance'].max()
    
    # Normalize distances to [0, 1] for scaling arc heights
    # This normalization is based on ALL records in the region
    if dist_max > dist_min:
        bedpe_df['dist_norm'] = (bedpe_df['distance'] - dist_min) / (dist_max - dist_min)
    else:
        # If all distances are the same (or only one record), treat as maximum (height = 1.0)
        bedpe_df['dist_norm'] = 1.0
    
    # Limit number of arcs if specified (AFTER normalization)
    if max_arcs is not None and len(bedpe_df) > max_arcs:
        bedpe_df = bedpe_df.head(max_arcs).copy()
        if log and verbose:
            log.write(f"Limiting to {max_arcs} arcs (normalization based on all {n_total_in_region} records in region)", verbose=verbose)
    
    if len(bedpe_df) == 0:
        return ax, bedpe_df
    
    # Determine if we should use score for coloring
    # Use score if available: either explicitly requested (color_by_score=True) 
    # or if score_col is provided, or if 'score' column exists (default behavior)
    use_score_coloring = False
    if score_col and score_col in bedpe_df.columns:
        if color_by_score or score_col == 'score':  # Default to using 'score' column if it exists
            use_score_coloring = True
    elif 'score' in bedpe_df.columns and color_by_score:  # Use 'score' column if available and requested
        score_col = 'score'
        use_score_coloring = True
    elif 'score' in bedpe_df.columns and not score_col:  # Auto-use 'score' if no other score_col specified
        score_col = 'score'
        use_score_coloring = True
    
    # Get colormap if coloring by score
    if use_score_coloring:
        import matplotlib.cm as cm
        cmap_obj = cm.get_cmap(cmap)
        scores = bedpe_df[score_col].values
        score_min, score_max = scores.min(), scores.max()
        if score_max > score_min:
            norm_scores = (scores - score_min) / (score_max - score_min)
            score_norm = Normalize(vmin=score_min, vmax=score_max)
        else:
            norm_scores = np.ones(len(scores))
            score_norm = Normalize(vmin=score_min, vmax=score_max + 1)
    else:
        cmap_obj = None
        norm_scores = None
        score_norm = None
        scores = None
    
    # Set y-axis limits to 0-1.05 with padding
    # Note: Control points will be set above 1.0 to make curves reach 1.0
    # (Bezier curves reach ~75% of control point height)
    y_base = 0.0  # Base line for arcs
    y_top = 1.05  # Top of plot with padding (curve maximum should reach 1.0)
    ax.set_ylim((y_base, y_top))
    ax.set_yticks([])  # Hide y-axis ticks
    
    # Full panel spines (consistent with track / regional stacked panels)
    ax.spines['left'].set_visible(True)
    ax.spines['top'].set_visible(True)
    ax.spines['top'].set_zorder(1)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    
    # Draw a base line to show the chromosome track
    if region is not None:
        ax.axhline(y=y_base, color='gray', linewidth=0.5, alpha=0.3, linestyle='-', zorder=0)
    
    # Convert BEDPE coordinates for plotting
    # BEDPE uses 0-based, half-open [start, end)
    # We'll plot arcs connecting starts and ends separately
    # Note: distances and normalization already calculated above (before max_arcs limiting)
    
    # Adjust for track_start_i offset
    bedpe_df['x_start1'] = bedpe_df['chromStart1'] + track_start_i
    bedpe_df['x_end1'] = bedpe_df['chromEnd1'] + track_start_i
    bedpe_df['x_start2'] = bedpe_df['chromStart2'] + track_start_i
    bedpe_df['x_end2'] = bedpe_df['chromEnd2'] + track_start_i
    
    # Plot arcs and filled areas (all are cis contacts now)
    for idx, row in bedpe_df.iterrows():
        x_start1 = row['x_start1']
        x_end1 = row['x_end1']
        x_start2 = row['x_start2']
        x_end2 = row['x_end2']
        dist_norm = row['dist_norm']
        
        # Determine fill color based on score
        if cmap_obj is not None and norm_scores is not None:
            fill_color = cmap_obj(norm_scores[bedpe_df.index.get_loc(idx)])
        else:
            fill_color = color
        
        # Calculate arc height proportional to distance, normalized to [0, 1]
        # Maximum distance pair should have height = 1.0 (reaches top)
        # All other pairs scaled proportionally
        # For cubic Bezier curves with symmetric control points:
        # - Curve maximum occurs at t=0.5  
        # - Maximum height = 0.75 * control_point_height
        # To make curve reach dist_norm, we need: control_height = dist_norm / 0.75
        # For dist_norm=1.0, this gives control_height = 4/3 ≈ 1.333...
        # Use a slightly higher value to ensure we reach exactly 1.0 (account for floating point)
        desired_height = dist_norm  # Desired curve maximum height (0 to 1)
        if desired_height > 0:
            # Scale by 4/3 (1/0.75) and add small epsilon to ensure we reach target
            control_height = (desired_height / 0.75) * 1.001  # Small factor to ensure we reach exactly 1.0
        else:
            control_height = 0.0
        
        # Create two arcs: one for starts, one for ends
        # Arc 1: connecting start positions
        verts_start = [
            (x_start1, y_base),  # Start point 1 at base
            (x_start1, y_base + control_height),  # Control point 1 (left side, at peak)
            (x_start2, y_base + control_height),  # Control point 2 (right side, at peak)
            (x_start2, y_base),  # Start point 2 at base
        ]
        codes_start = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path_start = Path(verts_start, codes_start)
        
        # Arc 2: connecting end positions
        verts_end = [
            (x_end1, y_base),  # End point 1 at base
            (x_end1, y_base + control_height),  # Control point 1 (left side, at peak)
            (x_end2, y_base + control_height),  # Control point 2 (right side, at peak)
            (x_end2, y_base),  # End point 2 at base
        ]
        codes_end = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path_end = Path(verts_end, codes_end)
        
        # Create filled area between the two arcs
        # Trace: start1 -> start2 (arc1) -> end2 (base) -> end1 (arc2 reversed) -> start1 (base)
        verts_fill = []
        codes_fill = []
        
        # Start at chromStart1 base
        verts_fill.append((x_start1, y_base))
        codes_fill.append(Path.MOVETO)
        
        # Follow arc 1: start1 -> start2 (using Bezier curve)
        verts_fill.append((x_start1, y_base + control_height))  # Control point 1
        verts_fill.append((x_start2, y_base + control_height))  # Control point 2
        verts_fill.append((x_start2, y_base))  # End point
        codes_fill.extend([Path.CURVE4, Path.CURVE4, Path.CURVE4])
        
        # Draw line along base from start2 to end2
        verts_fill.append((x_end2, y_base))
        codes_fill.append(Path.LINETO)
        
        # Follow arc 2 backwards: end2 -> end1 (reverse Bezier curve)
        verts_fill.append((x_end2, y_base + control_height))  # Control point 2 (reversed)
        verts_fill.append((x_end1, y_base + control_height))  # Control point 1 (reversed)
        verts_fill.append((x_end1, y_base))  # End point
        codes_fill.extend([Path.CURVE4, Path.CURVE4, Path.CURVE4])
        
        # Close the polygon: draw line along base from end1 back to start1
        codes_fill.append(Path.CLOSEPOLY)
        verts_fill.append((0, 0))  # Dummy point for CLOSEPOLY
        
        path_fill = Path(verts_fill, codes_fill)
        
        # Draw filled area
        patch_fill = PathPatch(path_fill, facecolor=fill_color, edgecolor='none',
                              alpha=alpha, zorder=1)
        ax.add_patch(patch_fill)
        
        # Draw arc outlines (optional, for better visibility)
        patch_start = PathPatch(path_start, edgecolor=fill_color, facecolor='none',
                               alpha=min(1.0, alpha * 1.5), linewidth=linewidth * 0.5, zorder=2)
        patch_end = PathPatch(path_end, edgecolor=fill_color, facecolor='none',
                             alpha=min(1.0, alpha * 1.5), linewidth=linewidth * 0.5, zorder=2)
        ax.add_patch(patch_start)
        ax.add_patch(patch_end)
    
    # Set x-axis ticks and labels (same as plot_track)
    if region is not None:
        region_ticks = list(map('{:.3f}'.format, np.linspace(region[1], region[2], num=region_step).astype("int")/1000000))
        ax.set_xticks(np.linspace(track_start_i + region[1], track_start_i + region[2], num=region_step))
        ax.set_xticklabels(region_ticks, rotation=45)
        ax.set_xlim([track_start_i + region[1], track_start_i + region[2]])
        # Set x-axis label (same as plot_track)
        xlabel = "Chromosome " + str(region[0]) + " (MB)"
        ax.set_xlabel(xlabel, fontsize=10, family=track_font_family)
    
    # Add colorbar if using score-based coloring
    if use_score_coloring and score_norm is not None and scores is not None:
        # Create a ScalarMappable for the colorbar
        import matplotlib.cm as cm
        sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=score_norm)
        sm.set_array([])  # Empty array, we just need the colormap

        # inset_axes keeps the colorbar anchored when plot_panels / XAxisManager
        # repositions the parent axis (fig.add_axes figure coords go stale).
        cax = inset_axes(
            ax,
            width="20%",
            height="3%",
            loc="upper right",
            bbox_to_anchor=(0.02, 0.02, 0.98, 0.98),
            bbox_transform=ax.transAxes,
            borderpad=0.5,
        )
        cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
        cbar.set_label(score_col if score_col else 'Score', fontsize=8, family=track_font_family, labelpad=2)
        cbar.ax.tick_params(labelsize=7)
        # Position label and ticks on top of horizontal colorbar (same as plot_ld_block)
        cbar.ax.xaxis.set_label_position('top')
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.tick_params(labeltop=True, labelbottom=False)
        ax._gwaslab_arc_cbar_axes = cax
        fit_arc_colorbar_inside_panel(ax, cax, reset_anchor=True)
    
    if log and verbose:
        log.write(f"Plotted {len(bedpe_df)} cis arcs", verbose=verbose)
        if use_score_coloring:
            log.write(f"  -Colored by {score_col} (range: {score_min:.2f} - {score_max:.2f})", verbose=verbose)
    
    return ax, bedpe_df
