import matplotlib
from gwaslab.io.io_process_args import _update_args, _update_arg
from gwaslab.g_Log import Log
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

def set_plot_style(
    plot="plot_mqq",
    mode=None,
    pm: VizParamsManager | None = None,
    fig_args=None,
    save_args=None,
    save=None,
    scatter_args=None,
    legend_args=None,
    anno_args=None,
    highlight_anno_args=None,
    anno_args_single=None,
    fig_kwargs=None,
    save_kwargs=None,
    scatter_kwargs=None,
    legend_kwargs=None,
    line_kwargs=None,
    anno_kwargs=None,
    highlight_anno_kwargs=None,
    anno_single_kwargs=None,
    anno_style=None,
    anno_fontsize=None,
    arrow_kwargs=None,
    arm_scale=None,
    anno_height=None,
    anno_xshift=None,
    anno_fixed_arm_length=None,
    repel_force=None,
    region_anno_bbox_args=None,
    line_args=None,
    sig_line_color=None,
    suggestive_sig_line_color=None,
    sc_linewidth=None,
    qq_line_color=None,
    markeredgecolor=None,
    markeredgewidth=None,
    markerfacecolor=None,
    marker=None,
    font=None,
    fontsize=None,
    fontfamily=None,
    font_family="Arial",
    colors=None,
    marker_size=None,
    dpi=None,
    verbose=True,
    log=Log(),
):
    """
    Configure and return a consolidated plotting style dictionary for gwaslab visualizations.

    This helper centralizes common styling knobs and retrieves per-plot/mode defaults
    from `viz_aux_params.txt` via `VizParamsManager`. User-specified values override
    contextual defaults. The result is a style dictionary ready to be forwarded to
    plotting functions (e.g., `fig_args` to `plt.subplots`, `save_args` to `savefig`,
    `scatter_args` to scatter/Seaborn calls, and annotation/line options to helpers).

    Parameters
    ----------
    plot : str, default "plot_mqq"
        Plot identifier (e.g., "plot_mqq", "plot_region", "plot_qq").
    mode : str or None, optional
        Sub-mode when a plot has variants (e.g., "r" for regional, "qq" for QQ).
    pm : VizParamsManager or None, optional
        Parameter manager. If None, a manager is created and populated from
        `viz_aux_params.txt`.
    fig_args, save_args, scatter_args, legend_args : dict or None
        Style dictionaries for figure creation, saving, scatter rendering, and legend.
    anno_args, highlight_anno_args, anno_args_single : dict or None
        Annotation styling for all, highlighted, and per-SNP overrides.
    anno_style : str or None
        Annotation style (e.g., "right", "tight", "expand").
    anno_fontsize : int or None
        Font size for annotations.
    arrow_kwargs : dict or None
        Arrow style kwargs for annotation arms.
    arm_scale, anno_height : float or None
        Scaling factors for annotation arm length and vertical text placement.
    anno_xshift, anno_fixed_arm_length : float or None
        Horizontal shift applied to annotations, and a fixed arm length, respectively.
    repel_force : float or None
        Repulsion strength used in text adjustment routines.
    region_anno_bbox_args : dict or None
        Bounding box kwargs for region annotations.
    line_args : dict or None
        Line styling dict (e.g., for reference lines).
    sig_line_color, suggestive_sig_line_color : str or None
        Colors for significance and suggestive reference lines.
    sc_linewidth : int or None
        Line width for reference lines.
    qq_line_color : str or None
        Reference line color for QQ plots.
    markeredgecolor, markeredgewidth, markerfacecolor, marker : various or None
        Marker-level styling injected into `scatter_args`.
    font, fontsize, fontfamily, font_family : various
        Font configuration; `fontfamily` overrides `font_family`. `font` dict may
        contain `family` and `size`.
    colors : list or None
        Color palette for categorical hues.
    marker_size : int or None
        Default point size; mapped to `scatter_args['s']` if not provided.
    dpi : int or None
        Resolution to apply to both `fig_args['dpi']` and `save_args['dpi']`.
    verbose : bool, default True
        Toggle for progress logging.
    log : gwaslab.g_Log.Log
        Logger instance.

    Returns
    -------
    dict
        Consolidated style configuration containing:
        - `plot`, `mode`
        - `fig_args`, `save_args`, `scatter_args`, `legend_args`
        - `anno_args`, `highlight_anno_args`, `anno_args_single`, `anno_style`, `anno_fontsize`
        - `arrow_kwargs`, `arm_scale`, `anno_height`, `anno_xshift`, `anno_fixed_arm_length`, `repel_force`, `region_anno_bbox_args`
        - `line_args`, `sig_line_color`, `suggestive_sig_line_color`, `sc_linewidth`, `qq_line_color`
        - `font_family`, `fontsize`, `colors`, `marker_size`, `dpi`, `verbose`
    """
    if pm is None:
        pm = VizParamsManager()
        load_viz_config(pm)

    defaults_map = pm.defaults(plot, mode)

    default_fig_args = defaults_map.get("fig_args", {"figsize": (15, 5), "dpi": 200})
    default_save_args = defaults_map.get("save_args", {"dpi": 300, "facecolor": "white"})
    default_save = defaults_map.get("save", None)
    default_scatter_args = defaults_map.get("scatter_args", {})
    default_legend_args = {}
    default_anno_args = defaults_map.get("anno_args", {})
    default_highlight_anno_args = defaults_map.get("highlight_anno_args", {})
    default_anno_args_single = {}
    default_colors = defaults_map.get("colors", ["#597FBD", "#74BAD3"])
    default_fontsize = defaults_map.get("fontsize", 12)
    default_marker_size = defaults_map.get("marker_size", 40)
    default_anno_fontsize = defaults_map.get("anno_fontsize", 9)
    default_anno_style = defaults_map.get("anno_style", "right")
    default_arrow_kwargs = {}
    default_arm_scale = defaults_map.get("arm_scale", 1)
    default_anno_height = defaults_map.get("anno_height", 1)
    default_repel_force = defaults_map.get("repel_force", 0.03)
    default_region_anno_bbox_args = defaults_map.get("region_anno_bbox_args", None)
    default_line_args = {}
    default_sig_line_color = defaults_map.get("sig_line_color", "grey")
    default_suggestive_sig_line_color = defaults_map.get("suggestive_sig_line_color", "grey")
    default_sc_linewidth = defaults_map.get("sc_linewidth", 2)
    default_qq_line_color = defaults_map.get("qq_line_color", "grey")
    default_markeredgecolor = None
    default_markeredgewidth = None
    default_markerfacecolor = None
    default_marker = None

    if fontfamily is not None:
        font_family = fontfamily
    if font is not None:
        if isinstance(font, dict):
            font_family = font.get("family", font_family)
            fontsize = font.get("size", fontsize)

    fontsize = _update_arg(fontsize, default_fontsize)
    colors = _update_arg(colors, default_colors)
    marker_size = _update_arg(marker_size, default_marker_size)
    anno_fontsize = _update_arg(anno_fontsize, default_anno_fontsize)
    anno_style = _update_arg(anno_style, default_anno_style)
    arrow_kwargs = _update_args(arrow_kwargs, default_arrow_kwargs)
    arm_scale = _update_arg(arm_scale, default_arm_scale)
    anno_height = _update_arg(anno_height, default_anno_height)
    anno_xshift = _update_arg(anno_xshift, None)
    anno_fixed_arm_length = _update_arg(anno_fixed_arm_length, None)
    repel_force = _update_arg(repel_force, default_repel_force)
    region_anno_bbox_args = _update_arg(region_anno_bbox_args, default_region_anno_bbox_args)
    line_args = _update_args(line_args, default_line_args)
    sig_line_color = _update_arg(sig_line_color, default_sig_line_color)
    suggestive_sig_line_color = _update_arg(suggestive_sig_line_color, default_suggestive_sig_line_color)
    sc_linewidth = _update_arg(sc_linewidth, default_sc_linewidth)
    qq_line_color = _update_arg(qq_line_color, default_qq_line_color)
    markeredgecolor = _update_arg(markeredgecolor, default_markeredgecolor)
    markeredgewidth = _update_arg(markeredgewidth, default_markeredgewidth)
    markerfacecolor = _update_arg(markerfacecolor, default_markerfacecolor)
    marker = _update_arg(marker, default_marker)

    fig_args = _update_args(fig_args, default_fig_args)
    fig_kwargs = _update_args(fig_kwargs, fig_args)
    save_args = _update_args(save_args, default_save_args)
    save_kwargs = _update_args(save_kwargs, save_args)
    save = _update_arg(save, default_save)
    scatter_args = _update_args(scatter_args, default_scatter_args)
    scatter_kwargs = _update_args(scatter_kwargs, scatter_args)
    legend_args = _update_args(legend_args, default_legend_args)
    legend_kwargs = _update_args(legend_kwargs, legend_args)
    anno_args = _update_args(anno_args, default_anno_args)
    anno_kwargs = _update_args(anno_kwargs, anno_args)
    highlight_anno_args = _update_args(highlight_anno_args, default_highlight_anno_args)
    highlight_anno_kwargs = _update_args(highlight_anno_kwargs, highlight_anno_args)
    anno_args_single = _update_args(anno_args_single, default_anno_args_single)
    anno_single_kwargs = _update_args(anno_single_kwargs, anno_args_single)

    if dpi is None:
        dpi = defaults_map.get("dpi", default_fig_args.get("dpi", 200))
    fig_kwargs["dpi"] = dpi
    save_kwargs["dpi"] = dpi

    if "s" not in scatter_kwargs and marker_size is not None:
        scatter_kwargs["s"] = marker_size
    if marker is not None:
        scatter_kwargs.setdefault("marker", marker)
    if markeredgecolor is not None:
        scatter_kwargs.setdefault("edgecolor", markeredgecolor)
    if markeredgewidth is not None:
        scatter_kwargs.setdefault("lw", markeredgewidth)
        scatter_kwargs.setdefault("linewidth", markeredgewidth)
    if markerfacecolor is not None:
        scatter_kwargs.setdefault("facecolor", markerfacecolor)

    matplotlib.rc("font", family=font_family)
    try:
        import matplotlib.pyplot as plt
        plt.rcParams["font.size"] = fontsize
    except Exception:
        pass

    log.write(f"Configured plot style for {plot}:{mode}", verbose=verbose)
    return {
        "plot": plot,
        "mode": mode,
        "fig_args": fig_kwargs,
        "fig_kwargs": fig_kwargs,
        "save": save,
        "save_args": save_kwargs,
        "save_kwargs": save_kwargs,
        "scatter_args": scatter_kwargs,
        "scatter_kwargs": scatter_kwargs,
        "legend_args": legend_kwargs,
        "legend_kwargs": legend_kwargs,
        "anno_args": anno_kwargs,
        "anno_kwargs": anno_kwargs,
        "highlight_anno_args": highlight_anno_kwargs,
        "highlight_anno_kwargs": highlight_anno_kwargs,
        "anno_args_single": anno_single_kwargs,
        "anno_single_kwargs": anno_single_kwargs,
        "anno_style": anno_style,
        "anno_fontsize": anno_fontsize,
        "arrow_kwargs": arrow_kwargs,
        "arm_scale": arm_scale,
        "anno_height": anno_height,
        "anno_xshift": anno_xshift,
        "anno_fixed_arm_length": anno_fixed_arm_length,
        "repel_force": repel_force,
        "region_anno_bbox_args": region_anno_bbox_args,
        "line_args": line_args,
        "sig_line_color": sig_line_color,
        "suggestive_sig_line_color": suggestive_sig_line_color,
        "sc_linewidth": sc_linewidth,
        "qq_line_color": qq_line_color,
        "font_family": font_family,
        "fontsize": fontsize,
        "colors": colors,
        "marker_size": marker_size,
        "dpi": fig_kwargs.get("dpi"),
        "verbose": verbose,
    }
