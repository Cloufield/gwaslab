import matplotlib
from gwaslab.io.io_process_kwargs import _update_kwargs, _update_arg
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

def set_plot_style(
    plot="plot_mqq",
    mode=None,
    pm: VizParamsManager | None = None,
    fig_kwargs=None,
    save_kwargs=None,
    save=None,
    scatter_kwargs=None,
    qq_scatter_kwargs=None,
    scatter_kwargs_outlier=None,
    legend_kwargs=None,
    title_kwargs=None,
    region_title_kwargs=None,
    title_pos=None,
    anno_kwargs=None,
    anno_set=None,
    anno_alias=None,
    anno_d=None,
    highlight_anno_kwargs=None,
    anno_kwargs_single=None,
    err_kwargs=None,
    eaf_kwargs=None,
    snpr2_kwargs=None,
    font_kwargs=None,
    qqscatterargs=None,
    line_kwargs=None,
    anno_single_kwargs=None,
    anno_style=None,
    anno_fontsize=None,
    arrow_kwargs=None,
    arm_scale=None,
    anno_height=None,
    anno_xshift=None,
    anno_fixed_arm_length=None,
    repel_force=None,
    helper_line_kwargs=None,
    threshold_line_kwargs=None,
    reg_line_kwargs=None,
    histplot_kwargs=None,
    r2_kwargs=None,
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
    cbar_fontsize=None,
    cbar_font_family=None,
    track_font_family=None,
    colors=None,
    marker_size=None,
    dpi=None,
    verbose=True,
    log=Log(),
    **kwargs,
):
    """
    Configure and return a consolidated plotting style dictionary for gwaslab visualizations.

    preset common args for all plots.

    This helper centralizes common styling knobs and retrieves per-plot/mode defaults
    from `viz_aux_params.txt` via `VizParamsManager`. User-specified values override
    contextual defaults. The result is a style dictionary ready to be forwarded to
    plotting functions (e.g., `fig_kwargs` to `plt.subplots`, `save_kwargs` to `savefig`,
    `scatter_kwargs` to scatter/Seaborn calls, and annotation/line options to helpers).

    Parameters
    ----------
    plot : str, default "plot_mqq"
        Plot identifier (e.g., "plot_mqq", "plot_region", "plot_qq").
    mode : str or None, optional
        Sub-mode when a plot has variants (e.g., "r" for regional, "qq" for QQ).
    pm : VizParamsManager or None, optional
        Parameter manager. If None, a manager is created and populated from
        `viz_aux_params.txt`.
    fig_kwargs, save_kwargs, scatter_kwargs, legend_kwargs : dict or None
        Style dictionaries for figure creation, saving, scatter rendering, and legend.
    anno_kwargs, highlight_anno_kwargs, anno_kwargs_single : dict or None
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
    line_kwargs : dict or None
        Line styling dict (e.g., for reference lines).
    sig_line_color, suggestive_sig_line_color : str or None
        Colors for significance and suggestive reference lines.
    sc_linewidth : int or None
        Line width for reference lines.
    qq_line_color : str or None
        Reference line color for QQ plots.
    markeredgecolor, markeredgewidth, markerfacecolor, marker : various or None
        Marker-level styling injected into `scatter_kwargs`.
    font, fontsize, fontfamily, font_family : various
        Font configuration; `fontfamily` overrides `font_family`. `font` dict may
        contain `family` and `size`.
    colors : list or None
        Color palette for categorical hues.
    marker_size : int or None
        Default point size; mapped to `scatter_kwargs['s']` if not provided.
    dpi : int or None
        Resolution to apply to both `fig_kwargs['dpi']` and `save_kwargs['dpi']`.
    verbose : bool, default True
        Toggle for progress logging.
    log : gwaslab.g_Log.Log
        Logger instance.

    Returns
    -------
    dict
        Consolidated style configuration containing:
        - `plot`, `mode`
        - `fig_kwargs`, `save_kwargs`, `scatter_kwargs`, `legend_kwargs`
        - `anno_kwargs`, `highlight_anno_kwargs`, `anno_kwargs_single`, `anno_style`, `anno_fontsize`
        - `arrow_kwargs`, `arm_scale`, `anno_height`, `anno_xshift`, `anno_fixed_arm_length`, `repel_force`,
        - `line_kwargs`, `sig_line_color`, `suggestive_sig_line_color`, `sc_linewidth`, `qq_line_color`
        - `font_family`, `fontsize`, `colors`, `marker_size`, `dpi`, `verbose`
    """

    if pm is None:
        pm = VizParamsManager()
        load_viz_config(pm)

    defaults_map = pm.defaults(plot, mode)
    default_fig_kwargs = defaults_map.get("fig_kwargs", {"figsize": (15, 5), "dpi": 200})
    default_save_kwargs = defaults_map.get("save_kwargs", {"dpi": 300, "facecolor": "white"})
    default_save = defaults_map.get("save", None)
    default_scatter_kwargs = defaults_map.get("scatter_kwargs", {})
    default_qq_scatter_kwargs = defaults_map.get("qq_scatter_kwargs", {})
    default_scatter_kwargs_outlier = defaults_map.get("scatter_kwargs_outlier", {})
    default_legend_kwargs = {}
    default_anno_kwargs = defaults_map.get("anno_kwargs", {})
    default_anno_set = defaults_map.get("anno_set", [])
    default_anno_alias = defaults_map.get("anno_alias", {})
    default_anno_d = defaults_map.get("anno_d", {})
    default_highlight_anno_kwargs = defaults_map.get("highlight_anno_kwargs", {})
    default_anno_kwargs_single = {}
    default_colors = defaults_map.get("colors", ["#597FBD", "#74BAD3"])
    default_fontsize = defaults_map.get("fontsize", 12)
    default_marker_size = defaults_map.get("marker_size", [5,20])
    default_anno_fontsize = defaults_map.get("anno_fontsize", 9)
    default_anno_style = defaults_map.get("anno_style", "right")
    default_arrow_kwargs = {}
    default_arm_scale = defaults_map.get("arm_scale", 1)
    default_anno_height = defaults_map.get("anno_height", 1)
    default_repel_force = defaults_map.get("repel_force", 0.03)
    default_line_kwargs = {}
    default_helper_line_kwargs = defaults_map.get("helper_line_kwargs", {"color":'black', "linestyle":'-',"lw":1})
    default_threshold_line_kwargs = defaults_map.get("threshold_line_kwargs", {})
    default_reg_line_kwargs = defaults_map.get("reg_line_kwargs", {})
    default_histplot_kwargs = defaults_map.get("histplot_kwargs", {})
    default_r2_kwargs = defaults_map.get("r2_kwargs", {})
    default_err_kwargs = defaults_map.get("err_kwargs", {"ecolor":"#cccccc","elinewidth":1})
    default_eaf_kwargs = defaults_map.get("eaf_kwargs", {})
    default_snpr2_kwargs = defaults_map.get("snpr2_kwargs", {})
    default_font_kwargs = defaults_map.get("font_kwargs", {
        "fontsize": default_fontsize,
        "family": "sans",
        "fontname": "Arial",
    })
    default_sig_line_color = defaults_map.get("sig_line_color", "grey")
    default_suggestive_sig_line_color = defaults_map.get("suggestive_sig_line_color", "grey")
    default_sc_linewidth = defaults_map.get("sc_linewidth", 2)
    default_qq_line_color = defaults_map.get("qq_line_color", "grey")
    default_markeredgecolor = None
    default_markeredgewidth = None
    default_markerfacecolor = None
    default_marker = None
    default_title_kwargs = defaults_map.get("title_kwargs", {})
    default_region_title_kwargs = defaults_map.get("region_title_kwargs", {})
    default_title_pos = defaults_map.get("title_pos", None)

    if fontfamily is not None:
        font_family = fontfamily
    if font is not None:
        if isinstance(font, dict):
            font_family = font.get("family", font_family)
            fontsize = font.get("size", fontsize)

    fontsize = _update_arg(fontsize, default_fontsize)
    font_family = _update_arg(font_family, fontfamily)
    cbar_fontsize = _update_arg(cbar_fontsize, fontsize)
    cbar_font_family = _update_arg(cbar_font_family, font_family)
    track_font_family = _update_arg(track_font_family, font_family)
    colors = _update_arg(colors, default_colors)
    marker_size = _update_arg(marker_size, default_marker_size)
    anno_fontsize = _update_arg(anno_fontsize, default_anno_fontsize)
    anno_style = _update_arg(anno_style, default_anno_style)
    arrow_kwargs = _update_kwargs(arrow_kwargs, default_arrow_kwargs)
    arm_scale = _update_arg(arm_scale, default_arm_scale)
    anno_height = _update_arg(anno_height, default_anno_height)
    anno_xshift = _update_arg(anno_xshift, None)
    anno_fixed_arm_length = _update_arg(anno_fixed_arm_length, None)
    repel_force = _update_arg(repel_force, default_repel_force)
    line_kwargs = _update_kwargs(line_kwargs, default_line_kwargs)
    helper_line_kwargs = _update_kwargs(helper_line_kwargs, default_helper_line_kwargs)
    threshold_line_kwargs = _update_kwargs(threshold_line_kwargs, default_threshold_line_kwargs)
    reg_line_kwargs = _update_kwargs(reg_line_kwargs, default_reg_line_kwargs)
    histplot_kwargs = _update_kwargs(histplot_kwargs, default_histplot_kwargs)
    r2_kwargs = _update_kwargs(r2_kwargs, default_r2_kwargs)
    sig_line_color = _update_arg(sig_line_color, default_sig_line_color)
    suggestive_sig_line_color = _update_arg(suggestive_sig_line_color, default_suggestive_sig_line_color)
    sc_linewidth = _update_arg(sc_linewidth, default_sc_linewidth)
    qq_line_color = _update_arg(qq_line_color, default_qq_line_color)
    markeredgecolor = _update_arg(markeredgecolor, default_markeredgecolor)
    markeredgewidth = _update_arg(markeredgewidth, default_markeredgewidth)
    markerfacecolor = _update_arg(markerfacecolor, default_markerfacecolor)
    marker = _update_arg(marker, default_marker)
    title_kwargs = _update_kwargs(title_kwargs, default_title_kwargs)
    region_title_kwargs = _update_kwargs(region_title_kwargs, default_region_title_kwargs)
    title_pos = _update_arg(title_pos, default_title_pos)

    fig_kwargs = _update_kwargs(fig_kwargs, default_fig_kwargs)
    fig_kwargs = _update_kwargs(fig_kwargs, fig_kwargs)
    save_kwargs = _update_kwargs(save_kwargs, default_save_kwargs)
    save_kwargs = _update_kwargs(save_kwargs, save_kwargs)
    save = _update_arg(save, default_save)
    scatter_kwargs = _update_kwargs(scatter_kwargs, default_scatter_kwargs)
    scatter_kwargs = _update_kwargs(scatter_kwargs, scatter_kwargs)
    qq_scatter_kwargs = _update_kwargs(qq_scatter_kwargs, default_qq_scatter_kwargs)
    qq_scatter_kwargs = _update_kwargs(qq_scatter_kwargs, qq_scatter_kwargs)
    qq_scatter_kwargs = _update_kwargs(qq_scatter_kwargs, qqscatterargs)
    scatter_kwargs_outlier = _update_kwargs(scatter_kwargs_outlier, default_scatter_kwargs_outlier)
    legend_kwargs = _update_kwargs(legend_kwargs, default_legend_kwargs)
    legend_kwargs = _update_kwargs(legend_kwargs, legend_kwargs)
    anno_kwargs = _update_kwargs(anno_kwargs, default_anno_kwargs)
    anno_kwargs = _update_kwargs(anno_kwargs, anno_kwargs)
    anno_set = _update_arg(anno_set, default_anno_set)
 
    anno_alias = _update_kwargs(anno_alias, default_anno_alias)
    anno_d = _update_kwargs(anno_d, default_anno_d)
    highlight_anno_kwargs = _update_kwargs(highlight_anno_kwargs, default_highlight_anno_kwargs)
    highlight_anno_kwargs = _update_kwargs(highlight_anno_kwargs, highlight_anno_kwargs)
    anno_kwargs_single = _update_kwargs(anno_kwargs_single, default_anno_kwargs_single)
    anno_single_kwargs = _update_kwargs(anno_single_kwargs, anno_kwargs_single)
    err_kwargs = _update_kwargs(err_kwargs, default_err_kwargs)
    eaf_kwargs = _update_kwargs(eaf_kwargs, default_eaf_kwargs)
    snpr2_kwargs = _update_kwargs(snpr2_kwargs, default_snpr2_kwargs)
    font_kwargs = _update_kwargs(font_kwargs, default_font_kwargs)
    # Remove alias to avoid Matplotlib error: both 'family' and 'fontfamily'
    if "fontfamily" in font_kwargs:
        if "family" not in font_kwargs:
            font_kwargs["family"] = font_kwargs["fontfamily"]
        del font_kwargs["fontfamily"]

    if dpi is None:
        dpi = defaults_map.get("dpi", default_fig_kwargs.get("dpi", 200))
    fig_kwargs["dpi"] = dpi
    save_kwargs["dpi"] = dpi

    if "s" not in scatter_kwargs and marker_size is not None:
        scatter_kwargs["s"] = marker_size[1]
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
        "fig_kwargs": fig_kwargs,
        "save": save,
        "save_kwargs": save_kwargs,
        "scatter_kwargs": scatter_kwargs,
        "scatter_kwargs_outlier": scatter_kwargs_outlier,
        "qq_scatter_kwargs": qq_scatter_kwargs,
        "legend_kwargs": legend_kwargs,
        "title_kwargs": title_kwargs,
        "title_pos": title_pos,
        "anno_kwargs": anno_kwargs,
        "anno_set": anno_set,
        "anno_alias": anno_alias,
        "anno_d": anno_d,
        "region_title_kwargs": region_title_kwargs,
        "highlight_anno_kwargs": highlight_anno_kwargs,
        "anno_kwargs_single": anno_single_kwargs,
        "anno_single_kwargs": anno_single_kwargs,
        "err_kwargs": err_kwargs,
        "eaf_kwargs": eaf_kwargs,
        "snpr2_kwargs": snpr2_kwargs,
        "font_kwargs": font_kwargs,
        "anno_style": anno_style,
        "anno_fontsize": anno_fontsize,
        "arrow_kwargs": arrow_kwargs,
        "arm_scale": arm_scale,
        "anno_height": anno_height,
        "anno_xshift": anno_xshift,
        "anno_fixed_arm_length": anno_fixed_arm_length,
        "repel_force": repel_force,
        "line_kwargs": line_kwargs,
        "helper_line_kwargs": helper_line_kwargs,
        "threshold_line_kwargs": threshold_line_kwargs,
        "reg_line_kwargs": reg_line_kwargs,
        "histplot_kwargs": histplot_kwargs,
        "r2_kwargs": r2_kwargs,
        "sig_line_color": sig_line_color,
        "suggestive_sig_line_color": suggestive_sig_line_color,
        "sc_linewidth": sc_linewidth,
        "qq_line_color": qq_line_color,
        "font_family": font_family,
        "fontsize": fontsize,
        "cbar_fontsize": cbar_fontsize,
        "cbar_font_family": cbar_font_family,
        "track_font_family": track_font_family,
        "colors": colors,
        "marker_size": marker_size,
        "dpi": fig_kwargs.get("dpi"),
        "verbose": verbose,
    }

def figure_kwargs_for_vector_plot(save, fig_kwargs, scatter_kwargs):
    if save is not None:
        if type(save) is not bool:
            sv = str(save)
            if len(sv) > 3:
                suf = sv[-3:].lower()
                if suf == "pdf" or suf == "svg":
                    fig_kwargs["dpi"] = 72
                    if scatter_kwargs is not None:
                        scatter_kwargs["rasterized"] = True
    return fig_kwargs, scatter_kwargs
