import pandas as pd
import seaborn as sns

from gwaslab.info.g_Log import Log


def draw_manhattan_panel(
    ax1,
    sumstats,
    snpid,
    palette,
    marker_size,
    style,
    linewidth,
    edgecolor,
    legend,
    scatter_kwargs,
    highlight,
    highlight_chrpos,
    highlight_color,
    density_color,
    density_range,
    density_trange,
    density_threshold,
    density_tpalette,
    density_palette,
    pinpoint,
    pinpoint_color,
    chrom,
    log=Log(),
    verbose=True,
):
    explicit = {"edgecolor", "edgecolors", "linewidth", "ax", "palette", "hue", "data", "legend", "style", "size", "sizes", "zorder", "s"}
    scatter_kwargs = {k: v for k, v in scatter_kwargs.items() if k not in explicit}

    marker_size = (marker_size[0], marker_size[1])
    highlight_i = pd.DataFrame()
    to_plot = None

    if len(highlight) > 0:
        to_plot = sumstats
        log.write(" -Creating background plot...", verbose=verbose)
        sns.scatterplot(
            data=to_plot,
            x='i',
            y='scaled_P',
            hue='chr_hue',
            palette=palette,
            legend=legend,
            style=style,
            size="s",
            sizes=marker_size,
            linewidth=linewidth,
            zorder=2,
            ax=ax1,
            edgecolor=edgecolor,
            **scatter_kwargs,
        )

        if pd.api.types.is_list_like(highlight[0]) and highlight_chrpos is False:
            for i, highlight_set in enumerate(highlight):
                log.write(" -Highlighting set {} target loci...".format(i + 1), verbose=verbose)
                sns.scatterplot(
                    data=to_plot.loc[to_plot["HUE"] == i],
                    x='i',
                    y='scaled_P',
                    hue="HUE",
                    palette={i: highlight_color[i % len(highlight_color)]},
                    legend=legend,
                    style=style,
                    size="s",
                    sizes=(marker_size[0] + 1, marker_size[1] + 1),
                    linewidth=linewidth,
                    zorder=3 + i,
                    ax=ax1,
                    edgecolor=edgecolor,
                    **scatter_kwargs,
                )
            highlight_i = to_plot.loc[~to_plot["HUE"].isna(), "i"].values
        else:
            log.write(" -Highlighting target loci...", verbose=verbose)
            sns.scatterplot(
                data=to_plot.loc[to_plot["HUE"] == 0],
                x='i',
                y='scaled_P',
                hue="HUE",
                palette={0: highlight_color},
                legend=legend,
                style=style,
                size="s",
                sizes=(marker_size[0] + 1, marker_size[1] + 1),
                linewidth=linewidth,
                zorder=3,
                ax=ax1,
                edgecolor=edgecolor,
                **scatter_kwargs,
            )
            if highlight_chrpos is False:
                highlight_i = to_plot.loc[to_plot[snpid].isin(highlight), "i"].values
            else:
                highlight_i = []
    else:
        if density_color is True:
            hue = "DENSITY_hue"
            s = "DENSITY"
            to_plot = sumstats.sort_values("DENSITY")
            to_plot["DENSITY_hue"] = to_plot["DENSITY"].astype("float")
            if density_range is None:
                density_range = (to_plot["DENSITY"].min(), to_plot["DENSITY"].max())
            if type(density_trange) is list:
                density_trange = tuple(density_trange)
            if type(density_range) is list:
                density_range = tuple(density_range)
            sns.scatterplot(
                data=to_plot.loc[to_plot["DENSITY"] <= density_threshold, :],
                x='i',
                y='scaled_P',
                hue=hue,
                palette=density_tpalette,
                legend=legend,
                style=style,
                size=s,
                sizes=(marker_size[0] + 1, marker_size[0] + 1),
                linewidth=linewidth,
                hue_norm=density_trange,
                zorder=2,
                ax=ax1,
                edgecolor=edgecolor,
                **scatter_kwargs,
            )
            sns.scatterplot(
                data=to_plot.loc[to_plot["DENSITY"] > density_threshold, :],
                x='i',
                y='scaled_P',
                hue=hue,
                palette=density_palette,
                legend=legend,
                style=style,
                size=s,
                sizes=marker_size,
                hue_norm=density_range,
                linewidth=linewidth,
                zorder=2,
                ax=ax1,
                edgecolor=edgecolor,
                **scatter_kwargs,
            )
        else:
            s = "s"
            hue = 'chr_hue'
            hue_norm = None
            to_plot = sumstats
            log.write(" -Creating background plot...", verbose=verbose)
            sns.scatterplot(
                data=to_plot,
                x='i',
                y='scaled_P',
                hue=hue,
                palette=palette,
                legend=legend,
                style=style,
                size=s,
                sizes=marker_size,
                hue_norm=hue_norm,
                linewidth=linewidth,
                edgecolor=edgecolor,
                zorder=2,
                ax=ax1,
                **scatter_kwargs,
            )

    ax1.set_rasterization_zorder(0)

    if len(pinpoint) > 0:
        if isinstance(pinpoint, str):
            pinpoint = [pinpoint]
        if pd.api.types.is_list_like(pinpoint[0]):
            for i, pinpoint_set in enumerate(pinpoint):
                if sum(sumstats[snpid].isin(pinpoint_set)) > 0:
                    to_pinpoint = sumstats.loc[sumstats[snpid].isin(pinpoint_set), :]
                    log.write(" -Pinpointing set {} target variants...".format(i + 1), verbose=verbose)
                    ax1.scatter(to_pinpoint["i"], to_pinpoint["scaled_P"], color=pinpoint_color[i % len(pinpoint_color)], zorder=100, s=marker_size[1] + 1)
        else:
            if sum(sumstats[snpid].isin(pinpoint)) > 0:
                to_pinpoint = sumstats.loc[sumstats[snpid].isin(pinpoint), :]
                log.write(" -Pinpointing target variants...", verbose=verbose)
                ax1.scatter(to_pinpoint["i"], to_pinpoint["scaled_P"], color=pinpoint_color, zorder=100, s=marker_size[1] + 1)

    return highlight_i
