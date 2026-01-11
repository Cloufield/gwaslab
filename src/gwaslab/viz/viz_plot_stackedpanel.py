"""
Plot stacked panels using Panel objects.

This module provides the plot_panels function for creating stacked multi-panel
figures from Panel objects, supporting different panel types like tracks and arcs.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional, Dict, Any, Tuple, Union
from gwaslab.viz.viz_aux_panel import Panel
from gwaslab.viz.viz_plot_track import plot_track
from gwaslab.viz.viz_plot_arc import plot_arc
from gwaslab.viz.viz_plot_ld_block import plot_ld_block
from gwaslab.viz.viz_aux_chromatin import _plot_chromatin_state
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_aux_xaxis_manager import XAxisManager
from gwaslab.info.g_Log import Log


def _add_panel_title(ax, title: str, title_pos: Optional[Union[str, Tuple[float, float]]], title_kwargs: Dict[str, Any]):
    """Helper function to add title to a panel."""
    if title_pos is None:
        title_pos_str = "left"
    elif isinstance(title_pos, str):
        title_pos_str = title_pos
    else:
        # Tuple position - use ax.text for custom positioning
        ax.text(
            title_pos[0], title_pos[1], title,
            transform=ax.transAxes, **title_kwargs
        )
        return  # Early return for tuple position
    
    # String position - use ax.set_title
    if title_pos_str == "left":
        ax.set_title(title, loc="left", **title_kwargs)
    elif title_pos_str == "right":
        ax.set_title(title, loc="right", **title_kwargs)
    elif title_pos_str == "center":
        ax.set_title(title, loc="center", **title_kwargs)
    else:
        ax.set_title(title, **title_kwargs)


def _get_default_height_ratio(panel_type: str) -> float:
    """
    Get default height ratio for a panel type.
    
    Different panel types have different default heights based on their content:
    - region: Large (4.0) - main regional plot with scatter and gene track
    - ld_block: Medium-large (4.0) - LD block visualization
    - track: Small (2.0) - simple genomic track
    - arc: Small (2.0) - contact arcs
    - chromatin: Small (1.0) - chromatin state bars
    - pipcs: Medium (3.0) - PIP and Credible Sets plots
    
    Parameters
    ----------
    panel_type : str
        Type of panel (e.g., "track", "arc", "region", "ld_block", "chromatin", "pipcs")
    
    Returns
    -------
    float
        Default height ratio for the panel type
    """
    default_ratios = {
        "region": 4,      # Regional plots need more space
        "ld_block": 4,    # LD blocks need moderate space
        "track": 2,       # Tracks are compact
        "arc": 2,         # Arcs are compact
        "chromatin": 1,   # Chromatin states are very compact
        "pipcs": 3,       # PIPCS plots need moderate space
    }
    return default_ratios.get(panel_type.lower(), 1.0)  # Default to 1.0 for unknown types


def _auto_adjust_fig_kwargs(
    fig_kwargs: Dict[str, Any],
    panels: List[Panel],
    total_subplot_count: int,
    height_ratios: List[float],
    subplot_height: float,
    log: Log,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Automatically adjust fig_kwargs height based on panel number and types.
    
    Only adjusts the height of the figure, preserving the width from user input
    or using default width if not provided.
    
    Parameters
    ----------
    fig_kwargs : dict
        Current fig_kwargs dictionary
    panels : list of Panel
        List of Panel objects
    total_subplot_count : int
        Total number of subplots (including multi-axis panels)
    height_ratios : list of float
        Height ratios for each subplot
    subplot_height : float
        Base height per subplot in inches
    log : Log
        Logger instance
    verbose : bool
        Whether to show verbose messages
    
    Returns
    -------
    dict
        Adjusted fig_kwargs dictionary (height adjusted, width preserved)
    """
    # Make a copy to avoid modifying the original
    adjusted_kwargs = fig_kwargs.copy()
    
    # Analyze panel types for DPI adjustment
    panel_types = [panel.get_type() for panel in panels]
    n_panels = len(panels)
    n_region = panel_types.count("region")
    n_ld_block = panel_types.count("ld_block")
    
    # Calculate height based on subplot_height and height_ratios
    calculated_height = subplot_height * sum(height_ratios)
    
    # Get or set width (cap at reasonable maximum to avoid too-wide figures)
    max_reasonable_width = 12.0
    default_width = 10.0
    
    if "figsize" not in adjusted_kwargs:
        # No figsize provided: use default width, calculated height
        adjusted_kwargs["figsize"] = (default_width, calculated_height)
        log.write(
            f" -Auto-adjusted figsize: ({default_width:.1f}, {calculated_height:.1f}) "
            f"based on {n_panels} panels ({total_subplot_count} subplots)",
            verbose=verbose
        )
    else:
        # User provided figsize: adjust width if too wide, respect user height if provided
        existing_figsize = adjusted_kwargs["figsize"]
        if isinstance(existing_figsize, (tuple, list)) and len(existing_figsize) >= 2:
            existing_width = existing_figsize[0]
            existing_height = existing_figsize[1]
            
            # Cap width at reasonable maximum (set_plot_style might set 15, which is too wide)
            if existing_width > max_reasonable_width:
                adjusted_width = max_reasonable_width
                log.write(
                    f" -Capped width from {existing_width:.1f} to {adjusted_width:.1f} "
                    f"(too wide, using max {max_reasonable_width:.1f})",
                    verbose=verbose
                )
            else:
                adjusted_width = existing_width
            
            # Respect user-provided height - if they explicitly set it, use it
            # Only calculate height if user didn't provide figsize at all
            final_height = existing_height
            log.write(
                f" -Using user-provided figsize: ({adjusted_width:.1f}, {final_height:.1f}) "
                f"(calculated height would be {calculated_height:.1f} from {total_subplot_count} subplots, "
                f"sum of height_ratios: {sum(height_ratios):.2f})",
                verbose=verbose
            )
            
            adjusted_kwargs["figsize"] = (adjusted_width, final_height)
        else:
            # Invalid figsize format, use default
            adjusted_kwargs["figsize"] = (default_width, calculated_height)
            log.write(
                f" -Auto-adjusted figsize: ({default_width:.1f}, {calculated_height:.1f}) "
                f"(invalid figsize format, using defaults)",
                verbose=verbose
            )
    
    # Adjust DPI based on complexity (only if not explicitly provided)
    if "dpi" not in adjusted_kwargs:
        # More panels or complex panels (region, ld_block) benefit from higher DPI
        if n_region > 0 or n_ld_block > 0 or total_subplot_count >= 5:
            adjusted_kwargs["dpi"] = 200
        elif total_subplot_count >= 3:
            adjusted_kwargs["dpi"] = 150
        else:
            adjusted_kwargs["dpi"] = 100
        log.write(
            f" -Auto-set DPI to {adjusted_kwargs['dpi']} based on plot complexity "
            f"({total_subplot_count} subplots, {n_region} region, {n_ld_block} ld_block panels)",
            verbose=verbose
        )
    
    return adjusted_kwargs


def plot_panels(
    panels: List[Panel],
    region: Optional[Tuple[int, int, int]] = None,
    height_ratios: Optional[List[float]] = None,
    hspace: float = 0.07,
    subplot_height: float = 1.0,
    titles: Optional[List[str]] = None,
    title_pos: Optional[Union[str, Tuple[float, float]]] = None,
    title_kwargs: Optional[Dict[str, Any]] = None,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    save: Optional[Union[str, bool]] = None,
    save_kwargs: Optional[Dict[str, Any]] = None,
    fontsize: int = 9,
    font_family: str = "Arial",
    align_xaxis: bool = True,
    region_step: int = 21,
    track_start_i: float = 0.0,
    verbose: bool = True,
    log: Log = Log(),
    **kwargs
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Create a stacked figure from a list of Panel objects.
    
    This function takes a list of Panel objects (e.g., track panels, arc panels)
    and creates a stacked figure with shared x-axis. Each panel is plotted
    using its stored parameters.
    
    Parameters
    ----------
    panels : list of Panel
        List of Panel objects to plot. Each panel should have a panel_type
        ("track", "arc", "ld_block", "region", "chromatin", "pipcs", etc.) and stored kwargs.
        Note: "region" and "ld_block" panel types require 2 axes each.
    region : tuple of (int, int, int), optional
        Genomic region as (chromosome, start, end) in 1-based coordinates.
        If None, will try to extract from panels. All panels should share
        the same region for proper x-axis alignment.
        Example: (1, 1000000, 2000000) for chr1:1000000-2000000
    height_ratios : list of float, optional
        Height ratios for each panel. If None, default ratios are used based on panel type:
        - region: 4.0 (larger for regional plots)
        - ld_block: 4.0 (medium-large for LD blocks)
        - track: 2.0 (compact for tracks)
        - arc: 2.0 (compact for arcs)
        - chromatin: 1.0 (very compact for chromatin states)
        - pipcs: 3.0 (medium for PIP and Credible Sets plots)
        Length should match number of panels.
    hspace : float, default=0.07
        Space between subplots (vertical spacing).
    subplot_height : float, default=2.0
        Height of each subplot in inches. Used to calculate total figure height
        if figsize is not provided in fig_kwargs.
    titles : list of str, optional
        List of titles for each panel. Length should match number of panels.
    title_pos : str or tuple, optional
        Position of titles. Can be "left", "right", "center", or (x, y) tuple.
    title_kwargs : dict, optional
        Keyword arguments for title styling (e.g., fontsize, fontfamily).
    fig_kwargs : dict, optional
        Additional keyword arguments for matplotlib figure creation
        (e.g., {'figsize': (10, 8), 'dpi': 200}).
    save : str, bool, or None, optional
        If str: file path to save figure.
        If True: save to default path.
        If None/False: skip saving.
    save_kwargs : dict, optional
        Additional arguments for saving (e.g., {'dpi': 300, 'bbox_inches': 'tight'}).
    fontsize : int, default=9
        Default font size for labels.
    font_family : str, default="Arial"
        Font family for text.
    align_xaxis : bool, default=True
        Whether to align x-axes across all panels. If True, applies consistent
        xlim and tick positions to all main axes (excludes position bars and
        gene tracks that use different coordinate systems).
    region_step : int, default=21
        Number of ticks to generate for x-axis when align_xaxis=True.
        Used to calculate tick positions using np.linspace.
    track_start_i : float, default=0.0
        X-axis offset for track-based plots. Used when calculating xlim
        for alignment: xlim = (track_start_i + region[1], track_start_i + region[2]).
    verbose : bool, default=True
        Whether to show progress messages.
    log : Log, default=Log()
        Logger instance for messages.
    **kwargs
        Additional keyword arguments passed to set_plot_style.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib figure object.
    axes : list of matplotlib.axes.Axes
        List of axes objects (one per panel).
    
    Examples
    --------
    >>> import gwaslab as gl
    >>> 
    >>> # Create panels
    >>> panel1 = gl.Panel(
    ...     "track",
    ...     track_path="genes.gtf",
    ...     region=(1, 1000000, 2000000),
    ...     color="#020080"
    ... )
    >>> panel2 = gl.Panel(
    ...     "arc",
    ...     bedpe_path="contacts.bedpe.gz",
    ...     region=(1, 1000000, 2000000),
    ...     color="#FF0000",
    ...     alpha=0.3
    ... )
    >>> 
    >>> # Create a region panel
    >>> panel3 = gl.Panel(
    ...     "region",
    ...     insumstats=sumstats,  # Use insumstats to preserve original data
    ...     region=(1, 1000000, 2000000),
    ...     vcf_path="ld_ref.vcf.gz",
    ...     build="38"
    ... )
    >>> 
    >>> # Plot panels together
    >>> fig, axes = gl.plot_panels(
    ...     [panel1, panel2, panel3],
    ...     region=(1, 1000000, 2000000),
    ...     titles=["Genes", "Contacts", "Regional Plot"],
    ...     save="stacked_panels.png"
    ... )
    """
    
    log.write("Start to create stacked panels plot...", verbose=verbose)
    
    # Validate panels
    if not panels or len(panels) == 0:
        raise ValueError("panels list cannot be empty")
    
    if not all(isinstance(p, Panel) for p in panels):
        raise TypeError("All items in panels must be Panel objects")
    
    n_panels = len(panels)
    log.write(f" -Number of panels: {n_panels}", verbose=verbose)
    
    # Check which panels need multiple axes (e.g., ld_block needs ax_pos and ax, region needs ax1 and ax3)
    panel_needs_multiple_axes = []
    total_subplot_count = 0
    for panel in panels:
        panel_type = panel.get_type()
        if panel_type == "ld_block":
            # ld_block needs 2 axes: ax_pos (position bar) and ax (main LD block)
            panel_needs_multiple_axes.append(True)
            total_subplot_count += 2
        elif panel_type == "region":
            # region needs 2 axes: ax1 (main scatter plot) and ax3 (gene track)
            panel_needs_multiple_axes.append(True)
            total_subplot_count += 2
        else:
            panel_needs_multiple_axes.append(False)
            total_subplot_count += 1
    
    log.write(f" -Total subplot count: {total_subplot_count} (including multi-axis panels)", verbose=verbose)
    
    # Extract region from panels if not provided
    if region is None:
        for panel in panels:
            panel_region = panel.get_kwarg("region")
            if panel_region is not None:
                region = panel_region
                log.write(f" -Extracted region from panels: {region}", verbose=verbose)
                break
        
        if region is None:
            raise ValueError(
                "region must be provided either as parameter or in panel kwargs. "
                "All panels should share the same region for x-axis alignment."
            )
    
    # Validate that all panels have the same region (or set it)
    for i, panel in enumerate(panels):
        panel_region = panel.get_kwarg("region")
        if panel_region is not None and panel_region != region:
            log.warning(
                f"Panel {i} has different region {panel_region} than specified "
                f"region {region}. Using specified region."
            )
        # Set region in panel kwargs if not present
        if panel_region is None:
            panel.set_kwarg("region", region)
    
    # Set up style
    # Remove figsize from fig_kwargs before passing to set_plot_style
    # so that _auto_adjust_fig_kwargs can calculate it properly
    fig_kwargs_for_style = fig_kwargs.copy() if fig_kwargs else {}
    if "figsize" in fig_kwargs_for_style:
        # Temporarily remove figsize to prevent set_plot_style from interfering
        figsize_backup = fig_kwargs_for_style.pop("figsize")
    else:
        figsize_backup = None
    
    style = set_plot_style(
        plot="plot_panels",
        fig_kwargs=fig_kwargs_for_style,
        save_kwargs=save_kwargs,
        save=save,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
        **kwargs
    )
    fig_kwargs = style.get("fig_kwargs", {})
    
    # Remove figsize from set_plot_style output (it may have set a default)
    # We'll calculate it properly in _auto_adjust_fig_kwargs
    if "figsize" in fig_kwargs:
        fig_kwargs.pop("figsize")
    
    # Restore user-provided figsize if it was provided
    if figsize_backup is not None:
        fig_kwargs["figsize"] = figsize_backup
    
    save_kwargs = style.get("save_kwargs", {})
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    
    # Set up title kwargs
    if title_kwargs is None:
        title_kwargs = {}
    if "family" not in title_kwargs:
        title_kwargs["family"] = font_family
    if "fontsize" not in title_kwargs:
        title_kwargs["fontsize"] = fontsize
    
    # Set up height ratios - expand for panels that need multiple axes
    if height_ratios is None:
        expanded_height_ratios = []
        for i, panel in enumerate(panels):
            panel_type = panel.get_type()
            base_ratio = _get_default_height_ratio(panel_type)
            
            if panel_needs_multiple_axes[i]:
                if panel_type == "ld_block":
                    # For ld_block: small ratio for ax_pos, main ratio for ax
                    expanded_height_ratios.extend([base_ratio * 0.1, base_ratio])  # ax_pos gets 10%, ax gets full
                elif panel_type == "region":
                    # For region: main ratio for ax1 (scatter plot), smaller ratio for ax3 (gene track)
                    expanded_height_ratios.extend([base_ratio, base_ratio * 0.5])  # ax1 gets full, ax3 gets 50%
                else:
                    expanded_height_ratios.extend([base_ratio, base_ratio])  # Default for unknown multi-axis types
            else:
                expanded_height_ratios.append(base_ratio)
        
        log.write(
            f" -Using default height ratios: {[f'{r:.2f}' for r in expanded_height_ratios]}",
            verbose=verbose
        )
        height_ratios = expanded_height_ratios
    else:
        # User provided height_ratios - expand them for multi-axis panels
        if len(height_ratios) != n_panels:
            raise ValueError(
                f"height_ratios length ({len(height_ratios)}) must match "
                f"number of panels ({n_panels})"
            )
        expanded_height_ratios = []
        for i, ratio in enumerate(height_ratios):
            if panel_needs_multiple_axes[i]:
                panel_type = panels[i].get_type()
                if panel_type == "ld_block":
                    # For ld_block: small ratio for ax_pos, main ratio for ax
                    expanded_height_ratios.extend([ratio * 0.1, ratio])  # ax_pos gets 10% of ratio, ax gets full ratio
                elif panel_type == "region":
                    # For region: main ratio for ax1, smaller ratio for ax3
                    expanded_height_ratios.extend([ratio, ratio * 0.5])  # ax1 gets full ratio, ax3 gets 50%
                else:
                    expanded_height_ratios.extend([ratio, ratio])  # Default
            else:
                expanded_height_ratios.append(ratio)
        height_ratios = expanded_height_ratios
    
    # Auto-adjust fig_kwargs based on panel number and types
    fig_kwargs = _auto_adjust_fig_kwargs(
        fig_kwargs=fig_kwargs,
        panels=panels,
        total_subplot_count=total_subplot_count,
        height_ratios=height_ratios,
        subplot_height=subplot_height,
        log=log,
        verbose=verbose
    )
    
    # Create figure and subplots
    log.write(f" -Creating figure with {total_subplot_count} subplots...", verbose=verbose)
    fig, axes = plt.subplots(
        total_subplot_count, 1,
        gridspec_kw={'height_ratios': height_ratios},
        **fig_kwargs
    )
    
    # Handle single subplot case (axes is not a list)
    if total_subplot_count == 1:
        axes = [axes]
    
    # Adjust spacing
    #plt.subplots_adjust(hspace=hspace)
    
    # Plot each panel
    log.write(" -Plotting panels...", verbose=verbose)
    axes_index = 0  # Track current position in axes array
    main_axes = []  # Track main axes for alignment (exclude position bars, gene tracks)
    exclude_axes = []  # Track axes to exclude from alignment
    # Start with axes from plt.subplots (preserves order from top to bottom)
    # all_axes[0] is top panel, all_axes[-1] is bottom panel
    all_axes = list(axes)  # Base axes from plt.subplots
    
    for i, panel in enumerate(panels):
        panel_type = panel.get_type()
        panel_kwargs = panel.get_kwargs()
        
        log.write(f"  -Panel {i+1}/{n_panels}: type='{panel_type}'", verbose=verbose)
        
        # Add panel kwargs that are common to all panel types
        panel_kwargs["fig"] = fig
        panel_kwargs["verbose"] = verbose
        panel_kwargs["log"] = log
        
        # Call appropriate plotting function based on panel type
        if panel_type == "track":
            # Required: track_path, region
            if "track_path" not in panel_kwargs:
                raise ValueError(f"Panel {i+1} (type='track') missing required parameter 'track_path'")
            
            # Get axes for this panel
            ax = axes[axes_index]
            panel_kwargs["ax"] = ax
            
            ax, texts = plot_track(**panel_kwargs)
            axes[axes_index] = ax
            # ax is already in all_axes (from plt.subplots), no need to append
            main_axes.append(ax)  # Track panel uses main axis
            
            # Add title if provided
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 1
            
        elif panel_type == "arc":
            # Required: bedpe_path, region (optional but recommended)
            if "bedpe_path" not in panel_kwargs:
                raise ValueError(f"Panel {i+1} (type='arc') missing required parameter 'bedpe_path'")
            
            # Get axes for this panel
            ax = axes[axes_index]
            panel_kwargs["ax"] = ax
            
            ax, bedpe_df = plot_arc(**panel_kwargs)
            axes[axes_index] = ax
            # ax is already in all_axes (from plt.subplots), no need to append
            main_axes.append(ax)  # Arc panel uses main axis
            
            # Add title if provided
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 1
            
        elif panel_type == "ld_block":
            # ld_block needs two axes: ax_pos (position bar) and ax (main LD block)
            # IMPORTANT: axes[axes_index] is TOP (from plt.subplots), axes[axes_index + 1] is BOTTOM
            if axes_index + 1 >= len(axes):
                raise ValueError(f"Panel {i+1} (type='ld_block') needs 2 axes but only {len(axes) - axes_index} available")
            
            ax_pos = axes[axes_index]  # Position bar (TOP, first axis, smaller ratio 0.1)
            ax = axes[axes_index + 1]  # Main LD block (BOTTOM, second axis, larger ratio 1.0)
            
            # Create a copy of panel_kwargs and remove fig (plot_ld_block doesn't accept it)
            ld_block_kwargs = panel_kwargs.copy()
            ld_block_kwargs.pop("fig", None)  # Remove fig if present
            
            ld_block_kwargs["ax"] = ax
            ld_block_kwargs["ax_pos"] = ax_pos
            
            # When both ax and ax_pos are provided, mode doesn't matter much,
            # but we set it to 'standalone' to ensure position bar is shown
            # (in regional mode, position bar behavior might differ)
            if "mode" not in ld_block_kwargs:
                ld_block_kwargs["mode"] = "standalone"
            
            # Save figure size before plot_ld_block (it may modify it in standalone mode)
            fig_size_before = fig.get_size_inches().copy()
            
            # Plot LD block - it will use the provided axes
            # Note: plot_ld_block returns (fig, ax), but we already have the axes
            axes_before = set(fig.axes)  # Track axes before plotting
            plot_ld_block(**ld_block_kwargs)
            
            # Restore figure size after plot_ld_block (it may have modified it)
            fig_size_after = fig.get_size_inches()
            if not np.allclose(fig_size_before, fig_size_after):
                fig.set_size_inches(fig_size_before[0], fig_size_before[1])
        
            
            # For ld_block: 
            # - ax_pos (position bar): uses genomic positions, SHOULD be aligned (main panel)
            # - ax (main LD block): uses rank-based coordinates, NOT genomic positions, should be excluded
            main_axes.append(ax_pos)  # Position bar uses genomic positions, should be aligned
            exclude_axes.append(ax)  # Main LD block uses rank-based coordinates, exclude from alignment
            
            # Add title if provided (on the main ax, not ax_pos)
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 2  # Skip both axes
            
        elif panel_type == "region":
            # region needs two axes: ax1 (main scatter plot) and ax3 (gene track)
            # IMPORTANT: axes[axes_index] is TOP (from plt.subplots), axes[axes_index + 1] is BOTTOM
            if axes_index + 1 >= len(axes):
                raise ValueError(f"Panel {i+1} (type='region') needs 2 axes but only {len(axes) - axes_index} available")
            
            # Required: sumstats or insumstats, region
            # Accept both for backward compatibility, but standardize to insumstats
            if "insumstats" in panel_kwargs:
                insumstats = panel_kwargs["insumstats"]
            elif "sumstats" in panel_kwargs:
                insumstats = panel_kwargs["sumstats"]
            else:
                raise ValueError(f"Panel {i+1} (type='region') missing required parameter 'sumstats' or 'insumstats'")
            
            ax1 = axes[axes_index]  # Main scatter plot (TOP, first axis, larger ratio 1.0)
            ax3 = axes[axes_index + 1]  # Gene track (BOTTOM, second axis, smaller ratio 0.5)
            
            # Use _mqqplot with mode="r" to plot regional plot
            # This will handle all the setup and call _plot_regional internally
            from gwaslab.viz.viz_plot_mqqplot import _mqqplot
            
            # Prepare kwargs for _mqqplot
            mqq_kwargs = panel_kwargs.copy()
            # Standardize to insumstats (plotting functions should not modify input)
            mqq_kwargs["insumstats"] = insumstats
            # Remove sumstats if present to avoid confusion
            mqq_kwargs.pop("sumstats", None)
            # Remove fig (not used by _mqqplot, it uses figax instead)
            mqq_kwargs.pop("fig", None)
            mqq_kwargs["mode"] = "r"
            # For mode="r", figax format is (fig, ax1, ax3) based on _process_layout
            mqq_kwargs["figax"] = (fig, ax1, ax3)
            
            # Call _mqqplot which will handle all setup and plotting
            axes_before = set(fig.axes)  # Track axes before plotting
            _mqqplot(**mqq_kwargs)
            
            # For region: both ax1 (main scatter) and ax3 (gene track) should be aligned
            # Both use genomic positions and should be aligned with other panels
            main_axes.append(ax1)  # Main scatter plot should be aligned
            main_axes.append(ax3)  # Gene track should also be aligned
            
            # Add title if provided (on the main ax1, not ax3)
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax1, titles[i], title_pos, title_kwargs)
            
            axes_index += 2  # Skip both axes
            
        elif panel_type == "chromatin":
            # Required: region_chromatin_files, region_chromatin_labels, region
            if "region_chromatin_files" not in panel_kwargs:
                raise ValueError(f"Panel {i+1} (type='chromatin') missing required parameter 'region_chromatin_files'")
            if "region_chromatin_labels" not in panel_kwargs:
                raise ValueError(f"Panel {i+1} (type='chromatin') missing required parameter 'region_chromatin_labels'")
            
            # Get axes for this panel
            ax = axes[axes_index]
            panel_kwargs["ax"] = ax
            
            # Calculate xlim_i from track_start_i if needed
            if "xlim_i" not in panel_kwargs:
                if region is not None:
                    panel_kwargs["xlim_i"] = [track_start_i + region[1]]
            
            # Call _plot_chromatin_state
            _plot_chromatin_state(**panel_kwargs)
            axes[axes_index] = ax
            main_axes.append(ax)  # Chromatin panel uses main axis
            
            # Add title if provided
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 1
            
        elif panel_type == "pipcs":
            # Required: pipcs_raw, region
            if "pipcs_raw" not in panel_kwargs:
                raise ValueError(f"Panel {i+1} (type='pipcs') missing required parameter 'pipcs_raw'")
            
            # Get axes for this panel
            ax = axes[axes_index]
            # _plot_cs uses figax parameter instead of ax
            panel_kwargs["figax"] = (fig, ax)
            
            # Call _plot_cs (returns fig, log)
            # Note: pipcs data should already be cleaned in _read_pipcs, but we validate here as well
            import pandas as pd
            pipcs_data = panel_kwargs["pipcs_raw"]
            if hasattr(pipcs_data, 'data') and not isinstance(pipcs_data, pd.DataFrame):
                pipcs_data = pipcs_data.data
            
            if not isinstance(pipcs_data, pd.DataFrame):
                raise ValueError(f"Panel {i+1} (type='pipcs'): pipcs_raw must be a DataFrame or Sumstats object")
            
            if "CHR" not in pipcs_data.columns:
                raise ValueError(
                    f"Panel {i+1} (type='pipcs'): pipcs_raw must have a 'CHR' column. "
                    f"Available columns: {list(pipcs_data.columns)}"
                )
            
            fig, log = _plot_cs(**panel_kwargs)
            axes[axes_index] = ax
            main_axes.append(ax)  # PIPCS panel uses main axis
            
            # Add title if provided
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 1
            
        else:
            raise ValueError(
                f"Unsupported panel type '{panel_type}' for panel {i+1}. "
                f"Supported types: 'track', 'arc', 'ld_block', 'region', 'chromatin', 'pipcs'"
            )
    print(all_axes)
    
    # Save figure size before alignment (in case anything modifies it)
    fig_size_before_align = fig.get_size_inches().copy()
    
    # Align x-axes if requested
    if align_xaxis and region is not None and len(main_axes) > 0:
        log.write(" -Aligning x-axes across panels...", verbose=verbose)
        # all_axes starts with axes from plt.subplots (preserves order from top to bottom)
        # all_axes[0] is top panel, all_axes[-1] is bottom panel
        # Any additional axes (like colorbars) are appended at the end
        xm = XAxisManager(
            all_axes=all_axes,  # All axes for spacing adjustment (from plt.subplots + any extras)
            region=region,
            region_step=region_step,
            track_start_i=track_start_i,
            fontsize=fontsize,
            font_family=font_family,
            fig=fig,
            adjust_spacing=True,
            verbose=verbose,
            log=log
        )
        xm.register_align_many(main_axes)  # Register axes for x-tick/label alignment
        xm.align()
        
        # Restore figure size if XAxisManager modified it
        fig_size_after_align = fig.get_size_inches()
        if not np.allclose(fig_size_before_align, fig_size_after_align):
            fig.set_size_inches(fig_size_before_align[0], fig_size_before_align[1])
    
    # Save figure
    save_figure(fig, save, keyword="panels", save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    log.write("Finished creating stacked panels plot.", verbose=verbose)
    return fig, axes
