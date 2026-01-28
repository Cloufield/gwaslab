import pandas as pd
import matplotlib.pyplot as plt
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.info.g_Log import Log
from gwaslab.io.io_process_kwargs import _extract_kwargs
import seaborn as sns
from gwaslab.viz.viz_aux_style_options import set_plot_style

def _plot_effect(insumstats, 
                 y=None, 
                 y_sort=None, 
                 group=None, 
                 x="BETA", 
                 se="SE", 
                 eaf="EAF", 
                 snpr2="SNPR2", 
                 ylabel="Variant",
                 eaf_panel=True, 
                 snpvar_panel=True, 
                 rename_dic=None, 
                 err_kwargs=None,
                 font_kwargs=None,
                 save=None,
                 title=None,
                 save_kwargs=None,
                 eaf_kwargs=None,
                 snpr2_kwargs=None,
                 fig_kwargs=None,
                 scatter_kwargs=None,
                 effect_label=None,
                 eaf_label=None,
                 snpr2_label=None,
                 xlim_eaf=None,
                 xlim_snpr2 = None,
                 log=Log(),
                 verbose=True,
                 legend_mode=1,
                 ncol=2,
                 gap=1,
                 fontsize=12,
                 font_family="Arial",
                 size=None,
                 hue=None,
                 style_col=None,
                 sort_kwargs=None,
                 **args):
    """Plot effect sizes with optional panels for EAF and SNPR2.

    Parameters
    ----------
    to_plot : pandas.DataFrame
        DataFrame containing the data to plot.
    y : list of str, optional
        Columns to use for y-axis labeling. Default is None.
    y_sort : list of str, optional
        Columns to sort by for y-axis order. Default is None.
    group : list of str, optional
        Columns to group by. Default is None.
    x : str, default="BETA"
        Column name for x-axis (effect size).
    se : str, default="SE"
        Column name for standard error.
    eaf : str, default="EAF"
        Column name for effect allele frequency.
    snpr2 : str, default="SNPR2"
        Column name for SNP R-squared.
    ylabel : str, default="Variant"
        Label for y-axis.
    eaf_panel : bool, default=True
        Whether to include EAF panel.
    snpvar_panel : bool, default=True
        Whether to include SNP variance panel.
    rename_dic : dict, optional
        Dictionary for renaming labels. Default is None.
    err_kwargs : dict, optional
        Arguments for error bars. Default is None.
    font_kwargs : dict, optional
        Font-related arguments. Default is None.
    save : str, optional
        File path to save the figure. Default is None.
    title : str, optional
        Title of the plot. Default is None.
    save_kwargs : dict, optional
        Additional arguments for saving. Default is None.
    eaf_kwargs : dict, optional
        Arguments for EAF panel. Default is None.
    snpr2_kwargs : dict, optional
        Arguments for SNPR2 panel. Default is None.
    fig_kwargs : dict, optional
        Figure size and DPI settings. Default is None.
    scatter_kwargs : dict, optional
        Scatter plot arguments. Default is None.
    effect_label : str, optional
        Custom label for effect size axis. Default is None.
    eaf_label : str, optional
        Custom label for EAF axis. Default is None.
    snpr2_label : str, optional
        Custom label for SNPR2 axis. Default is None.
    xlim_eaf : tuple, optional
        X-axis limits for EAF panel. Default is None.
    xlim_snpr2 : tuple, optional
        X-axis limits for SNPR2 panel. Default is None.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object.
    verbose : bool, default=True
        Whether to print progress messages.
    legend_mode : int, default=1
        Legend positioning mode.
    ncol : int, default=2
        Number of legend columns.
    gap : int, default=1
        Gap between groups.
    fontsize : int, default=12
        Font size for labels.
    font_family : str, default="Arial"
        Font family for text.
    size : str, optional
        Column name for marker size. Default is None.
    hue : str, optional
        Column name for color encoding. Default is None.
    style_col : str, optional
        Column name for marker style. Default is None.
    sort_kwargs : dict, optional
        Additional sorting arguments. Default is None.
    **args : dict
        Additional arguments passed to scatterplot.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object.

    Logic overview: prepare data and style (1–3), build y-axis label column and
    enforce row limit (4), sort and group (5), assign y positions per group (6),
    create 1–3 panels (7), draw effect panel (8), optional EAF/SNPR2 panels (9),
    legend (10), axis labels and save (11).

    How the plot is created
    -----------------------
    The figure is a horizontal forest-style layout:
    - Y-axis: one row per variant (or per study/variant combo). Rows are sorted
      by `group` and `y_sort`; within each group, a contiguous block of y-positions
      is assigned, with `gap` between groups. Tick labels come from the `y` column(s).
    - X-axis: effect size (e.g. BETA). The main panel is a scatter of point
      estimates with horizontal error bars (xerr=SE) and a vertical zero line.
    - Optional side panels (left to right): (1) main effect panel; (2) EAF as
      horizontal bars, if `eaf_panel` and column `eaf` exist; (3) SNPR2 as
      horizontal bars, if `snpvar_panel` and column `snpr2` exist. All panels
      share the same y-axis (same row order). Style (fig_kwargs, scatter_kwargs,
      err_kwargs, etc.) is resolved via set_plot_style and _extract_kwargs before
      drawing.
    """
    # --- 1. Resolve input to a DataFrame and take a working copy ---
    if hasattr(insumstats, 'data') and not isinstance(insumstats, pd.DataFrame):
        insumstats = insumstats.data
    to_plot = insumstats.copy()
    log.write("Starting effect-size plot creation...", verbose=verbose)
    log.write(" -Input contains {} variants.".format(len(to_plot)), verbose=verbose)

    # --- 2. Resolve style and per-plot keyword arguments ---
    style = set_plot_style(
        plot="plot_effect",
        fig_kwargs=fig_kwargs if fig_kwargs is not None else fig_kwargs,
        save_kwargs=save_kwargs if save_kwargs is not None else save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs if scatter_kwargs is not None else scatter_kwargs,
        err_kwargs=err_kwargs,
        eaf_kwargs=eaf_kwargs,
        snpr2_kwargs=snpr2_kwargs,
        font_kwargs=font_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", style.get("fig_kwargs", {}))
    save_kwargs = style.get("save_kwargs", style.get("save_kwargs", {}))
    scatter_kwargs = style.get("scatter_kwargs", style.get("scatter_kwargs", {}))
    err_kwargs = style["err_kwargs"]
    eaf_kwargs = style["eaf_kwargs"]
    snpr2_kwargs = style["snpr2_kwargs"]
    font_kwargs = style["font_kwargs"]
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    if sort_kwargs is None:
        sort_kwargs = {}
    log.write(" -Effect column: {}, SE column: {}.".format(x, se), verbose=verbose)

    # --- 3. Forward hue/size/style to seaborn and track legend titles ---
    legend_titles = []
    if hue is not None:
        args["hue"] = hue
        legend_titles.append(hue)

    if size is not None:
        args["size"] = size
        legend_titles.append(size)

    if style_col is not None:
        args["style_col"] = style_col
        legend_titles.append(style_col)

    save_kwargs =      _extract_kwargs("save", save_kwargs, locals())
    err_kwargs =       _extract_kwargs("err", err_kwargs, locals())
    scatter_kwargs =   _extract_kwargs("scatter", scatter_kwargs, locals())
    font_kwargs =      _extract_kwargs("font", font_kwargs, locals())
    
    def concat_cols(cols):
        """Join multiple column values into one label string (e.g. 'A-B-C')."""
        return "-".join(map(str, cols))

    # --- 4. Cap row count and build one “label per row” for the y-axis ---
    if len(to_plot) > 100:
        log.warning("Too many variants to plot ({}). Maximum is 100. Skipping.".format(len(to_plot)), verbose=verbose)
        return "Too many variants to plot"

    # Build y_name = column name holding the tick label; ensure that column exists.
    if hasattr(y, "__iter__") and not isinstance(y, str):
        y_cols = list(y)
        missing = [c for c in y_cols if c not in to_plot.columns]
        if missing:
            raise KeyError("y column(s) {} not found in sumstats. Available columns: {}".format(
                missing, list(to_plot.columns)))
        y_name = "-".join(str(c) for c in y_cols)
        to_plot[y_name] = to_plot[y_cols].apply(concat_cols, axis=1)
    else:
        if y is not None and y not in to_plot.columns:
            raise KeyError("y column '{}' not found in sumstats. Available columns: {}".format(
                y, list(to_plot.columns)))
        y_name = y

    y_label_display = y_name if isinstance(y_name, str) else ",".join(str(c) for c in (y_name or []))
    log.write(" -Y-axis label column: {}.".format(y_label_display), verbose=verbose)

    # --- 5. Sort and group: order rows for the forest layout ---
    if y_sort is None:
        y_sort = [c for c in ["CHR", "POS", "STUDY"] if c in to_plot.columns]
    if group is None:
        group = ["CHR", "POS"] + y_sort
    sort_columns = group + y_sort  # when group is default, this repeats y_sort (harmless)
    to_plot = to_plot.sort_values(by=sort_columns, **sort_kwargs)
    log.write(" -Grouping by: {}.".format(group), verbose=verbose)

    # --- 6. Assign y positions: one integer per row, groups stacked with gaps ---
    # to_plot is already sorted (section 5); row order = display order. We assign _VAR_INDEX
    # so that row i has the y-coordinate for the i-th line. Later, set_yticks uses
    # to_plot[y].values and to_plot[y_name].tolist() in that same row order, so the
    # i-th tick position and the i-th tick label always refer to the same row.
    #
    # Layout: larger _VAR_INDEX = higher on y-axis (top). Within each group, the first
    # row (in current order) gets the top of that group; groups are stacked with `gap`
    # between them.
    group_sizes = to_plot.groupby(group).size()
    group_end_after_gap = (group_sizes + gap).cumsum()
    to_plot["_GROUP_CUMSUM"] = to_plot.set_index(group).index.map(group_end_after_gap)
    row_rank_in_group = to_plot.groupby(group).cumcount()
    to_plot["_VAR_INDEX"] = to_plot["_GROUP_CUMSUM"] - row_rank_in_group
    y = "_VAR_INDEX"

    if rename_dic is None:
        rename_dic = {
            "BETA": "Per-allele effect size",
            "STUDY": "Study",
        }

    # --- 7. Create figure with 1–3 columns: effect (required), optional EAF, optional SNPR2 ---
    ncols = 1
    panel_names = ["effect"]
    if eaf_panel and eaf in to_plot.columns:
        ncols += 1
        panel_names.append("EAF")
    if snpvar_panel and snpr2 in to_plot.columns:
        ncols += 1
        panel_names.append("SNPR2")
    log.write(" -Creating figure with {} panel(s): {}.".format(ncols, "+".join(panel_names)), verbose=verbose)
    ax2, ax3 = None, None
    if ncols == 1:
        fig, ax1 = plt.subplots(ncols=ncols, **fig_kwargs)
        axes = [ax1]
    elif ncols == 2:
        fig, axes = plt.subplots(ncols=ncols, sharey=True, **fig_kwargs)
        ax1 = axes[0]
        if eaf_panel and eaf in to_plot.columns:
            ax2 = axes[1]
            ax3 = None
        else:
            ax2 = None
            ax3 = axes[1]
    else:
        fig, axes = plt.subplots(ncols=ncols, sharey=True, **fig_kwargs)
        ax1, ax2, ax3 = axes[0], axes[1], axes[2]

    # --- 8. Main panel: effect-size scatter, error bars, zero line, y-tick labels ---
    log.write(" -Plotting effect panel with {} variants...".format(len(to_plot)), verbose=verbose)
    sns.scatterplot(data=to_plot, x=x, y=y, ax=ax1, zorder=100, **{**args, **scatter_kwargs})
    ax1.errorbar(y=to_plot[y], x=to_plot[x], xerr=to_plot[se], **err_kwargs)
    ax1.axvline(x=0, linestyle="dashed", c="grey")
    y_ticks = to_plot[y].values
    y_labels = to_plot[y_name].astype(str).tolist()
    if len(y_ticks) != len(y_labels):
        raise ValueError(
            "Tick/label length mismatch ({} ticks vs {} labels). "
            "Ensure y_column(s) {} have one value per row.".format(
                len(y_ticks), len(y_labels), y_name if isinstance(y_name, str) else list(y_name)
            )
        )
    ax1.set_yticks(y_ticks, labels=y_labels, fontsize=fontsize, family=font_family)
    ax1.set_ylabel(ylabel, fontsize=fontsize, family=font_family)
    ax1.set_xlabel(x, fontsize=fontsize, family=font_family)
    if title is not None:
        ax1.set_title(title, fontsize=fontsize, family=font_family)

    # --- 9. Optional side panels: EAF (ax2) and SNPR2 (ax3) ---
    if ax2 is not None and eaf_panel and eaf in to_plot.columns:
        log.write(" -Adding EAF panel.", verbose=verbose)
        ax2.barh(y=to_plot[y], width=to_plot[eaf], zorder=100, **eaf_kwargs)
        ax2.set_xlabel(eaf, fontsize=fontsize, family=font_family)
        if xlim_eaf is not None:
            ax2.set_xlim(xlim_eaf)
    if ax3 is not None and snpvar_panel and snpr2 in to_plot.columns:
        log.write(" -Adding SNPR2 panel.", verbose=verbose)
        ax3.barh(y=to_plot[y], width=to_plot[snpr2], zorder=100, **snpr2_kwargs)
        ax3.set_xlabel(snpr2, fontsize=fontsize, family=font_family)
        if xlim_snpr2 is not None:
            ax3.set_xlim(xlim_snpr2)
    # --- 10. Legend (when hue/size/style_col are used) ---
    if legend_mode == 1 and len(legend_titles) > 0:
        sns.move_legend(
                ax1, "upper left",
                bbox_to_anchor=(1, 1), title=None, frameon=False, bbox_transform = axes[-1].transAxes, 
                title_fontproperties={"size":fontsize,"family":font_family},
                prop={"size": fontsize, "family": font_family}
        )
    # --- 11. Axis tick params, optional axis-label overrides, then save ---
    #handles, labels = ax1.get_legend_handles_labels()
    #if len(labels)>0:
    #    #new_labels = []
    #    #ncol = len(labels)
    #    max_col=0
    #    new_labels=[]
    #    new_labels_i = []
    #    previous_i = 0
    #    max_string_len=0
    #    for i in range(len(labels)):
    #        if len(labels[i]) > max_string_len:
    #            max_string_len = len(labels[i])
    #        if labels[i] in legend_titles:
    #            new_labels_i.append(i)
    #            col_number = i - previous_i
    #            if col_number > max_col:
    #                max_col = col_number
    #            previous_i = i
    #    for i in labels:
    #        new_labels.append(str(i).ljust(max_string_len))
    #    print(new_labels)
    #    new_labels_i.append(len(labels))
#
    #    legend_rows = []
    #    #new_labels_i[index+1] - i
    #    for index, i in enumerate(new_labels_i):
    #        if index<len(new_labels_i)-1:
    #            legend_row = ax1.legend(labels = new_labels[i:new_labels_i[index+1]],  
    #                                    handles= handles[i:new_labels_i[index+1]],
    #                                    loc="lower left", 
    #                                    bbox_to_anchor=(-0.2, 1.02 + 0.05*index), 
    #                                    ncol=max_col, 
    #                                    scatterpoints=1, 
    #                                    title=None, 
    #                                    borderpad=0,
    #                                    handletextpad=0.1,
    #                                    handlelength=0.7,
    #                                    borderaxespad =0,
    #                                    alignment = "left",
    #                                    fontsize=8,
    #                                    frameon=False)
    #            legend_rows.append(legend_row)
    #    for legend_row in legend_rows[:-1]:
    #        ax1.add_artist(legend_row)



    ax1.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family) 
    
    if effect_label is not None:
        ax1.set_xlabel(effect_label, fontsize=fontsize, family=font_family)
        ax1.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family)
    if eaf_label is not None and eaf in to_plot.columns and ax2 is not None:
        ax2.set_xlabel(eaf_label, fontsize=fontsize, family=font_family)
        ax2.tick_params(axis='x', labelsize=fontsize, labelfontfamily=font_family)
    if snpr2_label is not None and snpr2 in to_plot.columns and ax3 is not None:
        ax3.set_xlabel(snpr2_label, fontsize=fontsize, family=font_family)
        ax3.tick_params(axis='x', labelsize=fontsize, labelfontfamily=font_family)
    save_figure(fig, save, keyword="forest", save_kwargs=save_kwargs, log=log, verbose=verbose)
    log.write("Finished effect-size plot successfully.", verbose=verbose)

    return fig
