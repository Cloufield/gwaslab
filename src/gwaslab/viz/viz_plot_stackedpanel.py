"""
Plot stacked panels using Panel objects.

This module provides the plot_panels function for creating stacked multi-panel
figures from Panel objects, supporting different panel types like tracks and arcs.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional, Dict, Any, Tuple, Union, Sequence
from gwaslab.viz.viz_aux_panel import Panel
from gwaslab.viz.viz_plot_track import plot_track
from gwaslab.viz.viz_plot_arc import plot_arc, refit_arc_colorbars
from gwaslab.viz.viz_plot_ld_block import plot_ld_block
from gwaslab.viz.viz_aux_chromatin import _plot_chromatin_state
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_aux_xaxis_manager import XAxisManager
from gwaslab.viz.viz_aux_materialize import materialize_ag_panels
from gwaslab.viz.viz_plot_alphagenome import (
    finalize_ag_contact_axes,
    finalize_ag_panel_spines,
    plot_ag_tracks,
    plot_ag_overlay,
    plot_ag_contact,
    plot_ag_sashimi,
)
from gwaslab.info.g_Log import Log

_AG_PANEL_TYPES = frozenset({"ag_tracks", "ag_overlay", "ag_contact", "ag_sashimi"})


def _apply_shared_font_defaults(
    panel_kwargs: Dict[str, Any],
    panel_type: str,
    fontsize: int,
    font_family: str,
) -> None:
    """Inject plot_panels-level font settings into per-panel kwargs when unset."""
    if panel_type in ("track", "arc", "region"):
        panel_kwargs.setdefault("track_font_family", font_family)
    if panel_type in ("chromatin", "pipcs", "ld_block", "region"):
        panel_kwargs.setdefault("fontsize", fontsize)
        panel_kwargs.setdefault("font_family", font_family)
    if panel_type == "region":
        panel_kwargs.setdefault("anno_fontsize", fontsize)
        panel_kwargs.setdefault("cbar_fontsize", fontsize)
        panel_kwargs.setdefault("title_fontsize", fontsize)


def _as_axes_list(ax_or_axes: Union[plt.Axes, Sequence[plt.Axes]]) -> List[plt.Axes]:
    if isinstance(ax_or_axes, plt.Axes):
        return [ax_or_axes]
    return list(ax_or_axes)


def _left_label_extent_fig_x(ax: plt.Axes, fig: plt.Figure, renderer) -> float:
    """Leftmost figure-coordinate x of y-axis label and tick labels."""
    extents = []
    ylab = ax.yaxis.get_label()
    if ylab.get_text().strip():
        bb = ylab.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        extents.append(bb.x0)
    for tick in ax.get_yticklabels():
        if tick.get_text().strip():
            bb = tick.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            extents.append(bb.x0)
    if extents:
        return min(extents)
    return ax.get_position().x0


def _panel_vertical_center_fig_y(axes: Sequence[plt.Axes]) -> float:
    y0 = min(ax.get_position().y0 for ax in axes)
    y1 = max(ax.get_position().y1 for ax in axes)
    return y0 + (y1 - y0) / 2.0


def _text_width_fig(fig: plt.Figure, text: str, fontsize: float, family: str, renderer) -> float:
    probe = fig.text(0.0, 0.0, text, fontsize=fontsize, family=family, alpha=0.0)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width = probe.get_window_extent(renderer).transformed(fig.transFigure.inverted()).width
    probe.remove()
    return width


def _draw_margin_title(
    fig: plt.Figure,
    title: str,
    title_left_x: float,
    y: float,
    title_kwargs: Dict[str, Any],
) -> plt.Text:
    """Draw one panel title in the left margin, left-aligned at ``title_left_x``."""
    kw = dict(title_kwargs)
    fontsize = kw.pop("fontsize", 9)
    family = kw.pop("family", kw.pop("fontfamily", "Arial"))
    kw.pop("pad", None)
    kw.pop("xpad", None)
    kw.pop("margin_pad", None)
    return fig.text(
        title_left_x, y, title,
        ha="left",
        va="center",
        fontsize=fontsize,
        family=family,
        clip_on=False,
        **kw,
    )


def _add_panel_title(
    fig: plt.Figure,
    ax_or_axes: Union[plt.Axes, Sequence[plt.Axes]],
    title: str,
    title_pos: Optional[Union[str, Tuple[float, float]]],
    title_kwargs: Dict[str, Any],
    renderer=None,
    title_left_x: Optional[float] = None,
) -> Optional[plt.Text]:
    """Place a panel title outside the axes, left of y-axis labels."""
    axes = _as_axes_list(ax_or_axes)
    ax = axes[0]
    kw = dict(title_kwargs)
    fontsize = kw.pop("fontsize", 9)
    family = kw.pop("family", kw.pop("fontfamily", "Arial"))
    pad = kw.pop("pad", 6)
    xpad = kw.pop("xpad", 4)  # points between title column and y labels

    if isinstance(title_pos, tuple):
        return ax.text(
            title_pos[0], title_pos[1], title,
            transform=ax.transAxes,
            fontsize=fontsize,
            family=family,
            clip_on=False,
            **kw,
        )

    pos = title_pos if title_pos is not None else "left"

    if pos in ("left", "margin"):
        y = _panel_vertical_center_fig_y(axes)
        if title_left_x is not None:
            return _draw_margin_title(fig, title, title_left_x, y, {
                "fontsize": fontsize, "family": family, **kw,
            })
        if renderer is None:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
        ylab_x = min(_left_label_extent_fig_x(a, fig, renderer) for a in axes)
        xpad_frac = xpad / 72.0 / max(fig.get_figwidth(), 1e-6)
        text_w = _text_width_fig(fig, title, fontsize, family, renderer)
        return _draw_margin_title(
            fig, title, ylab_x - xpad_frac - text_w, y,
            {"fontsize": fontsize, "family": family, **kw},
        )

    subplot_h_in = max(ax.get_position().height * fig.get_figheight(), 1e-6)
    y_offset = pad / 72.0 / subplot_h_in

    if pos == "above-left":
        ha, x, va = "left", 0.0, "bottom"
    elif pos in ("above-right", "right"):
        ha, x, va = "right", 1.0, "bottom"
    elif pos in ("above-center", "center"):
        ha, x, va = "center", 0.5, "bottom"
    else:
        return ax.set_title(title, loc=pos, fontsize=fontsize, family=family, pad=pad, **kw)

    return ax.text(
        x, 1.0 + y_offset, title,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=fontsize,
        family=family,
        clip_on=False,
        **kw,
    )


def _apply_panel_titles(
    fig: plt.Figure,
    pending_titles: List[Tuple[Union[plt.Axes, Sequence[plt.Axes]], str]],
    title_pos: Optional[Union[str, Tuple[float, float]]],
    title_kwargs: Dict[str, Any],
) -> None:
    """Draw panel titles after all panels are plotted and layout is final."""
    if not pending_titles:
        return

    if isinstance(title_pos, tuple) or title_pos not in (None, "left", "margin"):
        for ax_or_axes, title in pending_titles:
            _add_panel_title(fig, ax_or_axes, title, title_pos, title_kwargs)
        return

    kw = dict(title_kwargs)
    fontsize = kw.get("fontsize", 9)
    family = kw.get("family", kw.get("fontfamily", "Arial"))
    xpad = kw.get("xpad", 6)
    margin_pad = kw.get("margin_pad", 0.01)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    xpad_frac = xpad / 72.0 / max(fig.get_figwidth(), 1e-6)

    min_ylab_x = 1.0
    max_title_w = 0.0
    for ax_or_axes, title in pending_titles:
        axes = _as_axes_list(ax_or_axes)
        min_ylab_x = min(
            min_ylab_x,
            min(_left_label_extent_fig_x(ax, fig, renderer) for ax in axes),
        )
        max_title_w = max(
            max_title_w,
            _text_width_fig(fig, title, fontsize, family, renderer),
        )

    title_left_x = min_ylab_x - xpad_frac - max_title_w
    required_left = max(fig.subplotpars.left, title_left_x - margin_pad)

    if required_left > fig.subplotpars.left + 1e-4:
        fig.subplots_adjust(left=required_left)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        min_ylab_x = 1.0
        for ax_or_axes, title in pending_titles:
            axes = _as_axes_list(ax_or_axes)
            min_ylab_x = min(
                min_ylab_x,
                min(_left_label_extent_fig_x(ax, fig, renderer) for ax in axes),
            )
            max_title_w = max(
                max_title_w,
                _text_width_fig(fig, title, fontsize, family, renderer),
            )
        title_left_x = min_ylab_x - xpad_frac - max_title_w

    for ax_or_axes, title in pending_titles:
        axes = _as_axes_list(ax_or_axes)
        y = _panel_vertical_center_fig_y(axes)
        _draw_margin_title(fig, title, title_left_x, y, title_kwargs)


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
        "region": 4,
        "ld_block": 4,
        "track": 2,
        "arc": 2,
        "chromatin": 1,
        "pipcs": 3,
        "ag_tracks": 2,
        "ag_overlay": 2,
        "ag_contact": 4,
        "ag_sashimi": 2,
    }
    return default_ratios.get(panel_type.lower(), 1.0)  # Default to 1.0 for unknown types


def _expand_height_ratio_for_panel(panel: Panel, ratio: float) -> List[float]:
    """Expand one panel-level height ratio into per-subplot ratios."""
    panel_type = panel.get_type()
    if panel_type == "ld_block":
        return [ratio * 0.1, ratio]
    if panel_type == "region":
        return [ratio, ratio * 0.5]
    n = panel.get_subplot_count()
    if panel_type == "ag_overlay" and n > 1:
        return [ratio * 0.5] * n
    if panel_type in _AG_PANEL_TYPES and n > 1:
        return [ratio] * n
    return [ratio]


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
            adjusted_width = existing_width
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
    variant_positions: Optional[List[float]] = None,
    variant_line_kwargs: Optional[Dict[str, Any]] = None,
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
    title_pos : str or tuple, default="left"
        Title placement. ``"left"`` (alias ``"margin"``) draws the title outside
        the panel, vertically centered and to the left of y-axis labels.
        ``"above-left"``, ``"above-right"``, and ``"above-center"`` place titles
        above the panel. A ``(x, y)`` tuple uses ``ax.text`` in transAxes.
    title_kwargs : dict, optional
        Title styling. ``fontsize``, ``family``, ``xpad`` (points left of y labels),
        and ``margin_pad`` (figure fraction) are supported for margin titles.
        ``pad`` applies to above-panel titles.
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
        Shared font size for panel labels, legends, and titles. Propagated to
        chromatin, pipcs, ld_block, and region panels (and ag_* panels) unless
        overridden on individual ``Panel`` kwargs. Track/arc gene labels scale
        via ``taf`` on the panel.
    font_family : str, default="Arial"
        Shared font family; also sets ``track_font_family`` on track/arc/region.
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

    # Extract region from panels if not provided (before materialize)
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

    panels = materialize_ag_panels(panels, log=log, verbose=verbose)

    panel_subplot_counts = [p.get_subplot_count() for p in panels]
    total_subplot_count = sum(panel_subplot_counts)
    log.write(
        f" -Total subplot count: {total_subplot_count} (including multi-axis panels)",
        verbose=verbose,
    )

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
    if "pad" not in title_kwargs:
        title_kwargs["pad"] = 6

    if title_pos is None:
        title_pos = "left"

    user_specified_height_ratios = height_ratios is not None

    if height_ratios is None:
        expanded_height_ratios = []
        for panel in panels:
            base_ratio = _get_default_height_ratio(panel.get_type())
            expanded_height_ratios.extend(_expand_height_ratio_for_panel(panel, base_ratio))
        log.write(
            f" -Using default height ratios: {[f'{r:.2f}' for r in expanded_height_ratios]}",
            verbose=verbose,
        )
        height_ratios = expanded_height_ratios
    else:
        if len(height_ratios) != n_panels:
            raise ValueError(
                f"height_ratios length ({len(height_ratios)}) must match "
                f"number of panels ({n_panels})"
            )
        expanded_height_ratios = []
        for panel, ratio in zip(panels, height_ratios):
            expanded_height_ratios.extend(_expand_height_ratio_for_panel(panel, ratio))
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
        gridspec_kw={'height_ratios': height_ratios, 'hspace': hspace},
        **fig_kwargs
    )
    
    # Handle single subplot case (axes is not a list)
    if total_subplot_count == 1:
        axes = [axes]
    
    # Plot each panel
    log.write(" -Plotting panels...", verbose=verbose)
    axes_index = 0  # Track current position in axes array
    main_axes = []  # Track main axes for alignment (exclude position bars, gene tracks)
    exclude_axes = []  # Track axes to exclude from alignment
    ag_axes: List[plt.Axes] = []
    ag_contact_axes: List[plt.Axes] = []
    pending_titles: List[Tuple[Union[plt.Axes, Sequence[plt.Axes]], str]] = []
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
        _apply_shared_font_defaults(panel_kwargs, panel_type, fontsize, font_family)
        
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
            
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append((ax, titles[i]))
            
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
            
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append((ax, titles[i]))
            
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

            # plot_ld_block expects sumstats; accept insumstats from Panel / Sumstats.Panel
            if "sumstats" not in ld_block_kwargs:
                if "insumstats" in ld_block_kwargs:
                    ld_block_kwargs["sumstats"] = ld_block_kwargs.pop("insumstats")
                else:
                    raise ValueError(
                        f"Panel {i+1} (type='ld_block') missing required parameter "
                        "'sumstats' or 'insumstats'"
                    )
            else:
                ld_block_kwargs.pop("insumstats", None)
            
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
            
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append(([ax_pos, ax], titles[i]))
            
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
            
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append(([ax1, ax3], titles[i]))
            
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
            
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append((ax, titles[i]))
            
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
            
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append((ax, titles[i]))
            
            axes_index += 1

        elif panel_type in _AG_PANEL_TYPES:
            if "bundle" not in panel_kwargs:
                raise ValueError(
                    f"Panel {i+1} (type='{panel_type}') missing required parameter 'bundle' "
                    "(or 'ag_spec' for deferred extraction)"
                )
            n_axes = panel.get_subplot_count()
            if axes_index + n_axes > len(axes):
                raise ValueError(
                    f"Panel {i+1} (type='{panel_type}') needs {n_axes} axes but only "
                    f"{len(axes) - axes_index} available"
                )
            panel_axes = axes[axes_index: axes_index + n_axes]
            plot_region = panel_kwargs.get("region", region)
            ag_kw = {
                k: v for k, v in panel_kwargs.items()
                if k not in (
                    "bundle", "ag_spec", "variant_context", "build",
                    "fig", "ax", "axes", "region", "verbose", "log",
                )
            }
            if "fontsize" not in ag_kw:
                ag_kw["fontsize"] = fontsize
            if panel_type == "ag_tracks":
                plot_ag_tracks(
                    bundle=panel_kwargs["bundle"],
                    axes=panel_axes,
                    region=plot_region,
                    verbose=verbose,
                    log=log,
                    **ag_kw,
                )
            elif panel_type == "ag_overlay":
                plot_ag_overlay(
                    bundle=panel_kwargs["bundle"],
                    axes=panel_axes,
                    region=plot_region,
                    verbose=verbose,
                    log=log,
                    **ag_kw,
                )
            elif panel_type == "ag_contact":
                plot_ag_contact(
                    bundle=panel_kwargs["bundle"],
                    axes=panel_axes,
                    region=plot_region,
                    track_start_i=track_start_i,
                    verbose=verbose,
                    log=log,
                    **ag_kw,
                )
                ag_contact_axes.extend(panel_axes)
            elif panel_type == "ag_sashimi":
                plot_ag_sashimi(
                    bundle=panel_kwargs["bundle"],
                    axes=panel_axes,
                    region=plot_region,
                    verbose=verbose,
                    log=log,
                    **ag_kw,
                )
            for ax in panel_axes:
                main_axes.append(ax)
            ag_axes.extend(panel_axes)
            if titles is not None and i < len(titles) and titles[i] is not None:
                pending_titles.append((panel_axes, titles[i]))
            axes_index += n_axes

        else:
            raise ValueError(
                f"Unsupported panel type '{panel_type}' for panel {i+1}. "
                f"Supported types: 'track', 'arc', 'ld_block', 'region', 'chromatin', 'pipcs', "
                f"'ag_tracks', 'ag_overlay', 'ag_contact', 'ag_sashimi'"
            )
    
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

        if ag_axes:
            finalize_ag_panel_spines(ag_axes)

        if ag_contact_axes:
            finalize_ag_contact_axes(ag_contact_axes, region, track_start_i)

        if variant_positions:
            vline_kw = variant_line_kwargs or {
                "color": "#FF0000",
                "linestyle": "--",
                "alpha": 0.7,
                "linewidth": 1.0,
            }
            for ax in main_axes:
                if ax in exclude_axes:
                    continue
                for pos in variant_positions:
                    ax.axvline(pos, **vline_kw)

        # Restore figure size if XAxisManager modified it
        fig_size_after_align = fig.get_size_inches()
        if not np.allclose(fig_size_before_align, fig_size_after_align):
            fig.set_size_inches(fig_size_before_align[0], fig_size_before_align[1])

    refit_arc_colorbars(main_axes)

    _apply_panel_titles(fig, pending_titles, title_pos, title_kwargs)
    
    # Save figure
    save_figure(fig, save, keyword="panels", save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    log.write("Finished creating stacked panels plot.", verbose=verbose)
    return fig, axes
