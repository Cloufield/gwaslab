"""
Plot stacked panels using Panel objects.

This module provides the plot_panels function for creating stacked multi-panel
figures from Panel objects, supporting different panel types like tracks and arcs.
"""

import matplotlib.pyplot as plt
from typing import List, Optional, Dict, Any, Tuple, Union
from gwaslab.viz.viz_aux_panel import Panel
from gwaslab.viz.viz_plot_track import plot_track
from gwaslab.viz.viz_plot_arc import plot_arc
from gwaslab.viz.viz_plot_ld_block import plot_ld_block
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
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


def plot_panels(
    panels: List[Panel],
    region: Optional[Tuple[int, int, int]] = None,
    height_ratios: Optional[List[float]] = None,
    hspace: float = 0.07,
    subplot_height: float = 2.0,
    titles: Optional[List[str]] = None,
    title_pos: Optional[Union[str, Tuple[float, float]]] = None,
    title_kwargs: Optional[Dict[str, Any]] = None,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    save: Optional[Union[str, bool]] = None,
    save_kwargs: Optional[Dict[str, Any]] = None,
    fontsize: int = 9,
    font_family: str = "Arial",
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
        ("track", "arc", etc.) and stored kwargs.
    region : tuple of (int, int, int), optional
        Genomic region as (chromosome, start, end) in 1-based coordinates.
        If None, will try to extract from panels. All panels should share
        the same region for proper x-axis alignment.
        Example: (1, 1000000, 2000000) for chr1:1000000-2000000
    height_ratios : list of float, optional
        Height ratios for each panel. If None, all panels get equal height.
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
    >>> # Plot panels together
    >>> fig, axes = gl.plot_panels(
    ...     [panel1, panel2],
    ...     region=(1, 1000000, 2000000),
    ...     titles=["Genes", "Contacts"],
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
    
    # Check which panels need multiple axes (e.g., ld_block needs ax_pos and ax)
    panel_needs_multiple_axes = []
    total_subplot_count = 0
    for panel in panels:
        panel_type = panel.get_type()
        if panel_type == "ld_block":
            # ld_block needs 2 axes: ax_pos (position bar) and ax (main LD block)
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
    style = set_plot_style(
        plot="plot_panels",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
        **kwargs
    )
    fig_kwargs = style.get("fig_kwargs", {})
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
            if panel_needs_multiple_axes[i]:
                # For ld_block: small ratio for ax_pos, main ratio for ax
                expanded_height_ratios.extend([0.1, 1.0])  # ax_pos gets 0.1, ax gets 1.0
            else:
                expanded_height_ratios.append(1.0)
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
                # For ld_block: small ratio for ax_pos, main ratio for ax
                expanded_height_ratios.extend([ratio * 0.1, ratio])  # ax_pos gets 10% of ratio, ax gets full ratio
            else:
                expanded_height_ratios.append(ratio)
        height_ratios = expanded_height_ratios
    
    # Set up figure size
    if "figsize" not in fig_kwargs:
        # Calculate width (default 10 inches)
        width = fig_kwargs.get("figsize", (10, 1))[0] if isinstance(
            fig_kwargs.get("figsize"), (tuple, list)
        ) else 10
        # Calculate height based on subplot_height and height_ratios
        total_height = subplot_height * sum(height_ratios)
        fig_kwargs["figsize"] = (width, total_height)
    
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
    plt.subplots_adjust(hspace=hspace)
    
    # Plot each panel
    log.write(" -Plotting panels...", verbose=verbose)
    axes_index = 0  # Track current position in axes array
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
            
            # Add title if provided
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 1
            
        elif panel_type == "ld_block":
            # ld_block needs two axes: ax_pos (position bar) and ax (main LD block)
            if axes_index + 1 >= len(axes):
                raise ValueError(f"Panel {i+1} (type='ld_block') needs 2 axes but only {len(axes) - axes_index} available")
            
            ax_pos = axes[axes_index]  # Position bar (first)
            ax = axes[axes_index + 1]  # Main LD block (second)
            
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
            
            # Plot LD block - it will use the provided axes
            # Note: plot_ld_block returns (fig, ax), but we already have the axes
            plot_ld_block(**ld_block_kwargs)
            # Keep the axes references (they're already in the axes array)
            
            # Add title if provided (on the main ax, not ax_pos)
            if titles is not None and i < len(titles) and titles[i] is not None:
                _add_panel_title(ax, titles[i], title_pos, title_kwargs)
            
            axes_index += 2  # Skip both axes
            
        else:
            raise ValueError(
                f"Unsupported panel type '{panel_type}' for panel {i+1}. "
                f"Supported types: 'track', 'arc', 'ld_block'"
            )
    
    # Save figure
    save_figure(fig, save, keyword="panels", save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    log.write("Finished creating stacked panels plot.", verbose=verbose)
    return fig, axes
