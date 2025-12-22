import matplotlib.pyplot as plt
import seaborn as sns
from gwaslab.info.g_Log import Log
from gwaslab.io.io_process_kwargs import _update_arg
from gwaslab.io.io_process_kwargs import _update_kwargs
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

def _plot( associations, 
            values="Beta",
            sort="P-value",
            mode="gwaslab",
            fontsize=12,
            font_family="Arial",
            cmap=None,
            ylabel="Y",
            xlabel="X",
            xlim=None,
            ylim=None,
            yticks=None,
            xticks=None,
            ytick_labels=None,
            xtick_labels=None,
            title="Title",
            title_pad=1.08, 
            title_fontsize=13,
            title_kwargs=None,
            linewidth=None,
            linestyle=None,
            linecolor=None,
            dpi=200,
            anno_kwargs=None,
            err_kwargs=None,
            fig_kwargs= None,
            scatter_kwargs=None,
            save=None,
            save_kwargs=None,
            log=Log(),
            verbose=True,
            **args               
            ):
    
    # update args
    
    cmap = _update_arg(cmap, "RdBu")

    style = set_plot_style(
        plot="plot_template",
        fig_kwargs=fig_kwargs or {"figsize":(10,10)},
        save_kwargs=save_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig,ax = plt.subplots(**(style.get("fig_kwargs", {})))

    # draw lines
    #horizontal_line = ax.axhline(y=1, linewidth = linewidth,
    #                        linestyle="--",
    #                        color=linecolor,zorder=1000)
    #vertical_line   = ax.axvline(x=1, linewidth = linewidth,
    #                        linestyle="--",
    #                        color=linecolor,zorder=1000)

    # ticks
    if xticks is not None:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks,
                            fontsize=fontsize,
                            family=font_family)    
    
    if yticks is not None:
        ax.set_xticks(yticks)
        ax.set_xticklabels(yticks,
                            fontsize=fontsize,
                            family=font_family)    

    # labels
    ax.set_xlabel(xlabel,
                  fontsize=fontsize, 
                  fontfamily=font_family)
    ax.set_ylabel(ylabel,
                  fontsize=fontsize, 
                  fontfamily=font_family)

    ax.tick_params(axis='x', 
                    labelsize=fontsize,
                    labelfontfamily=font_family)
    ax.tick_params(axis='y', 
                    labelsize=fontsize,
                    labelfontfamily=font_family)
    
    # title
    title_pad = title_pad -0.05
    fig.suptitle(title , 
                 fontsize = title_fontsize, 
                 x=0.5,
                 y=title_pad)
    
    # spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    save_kwargs = style.get("save_kwargs", style.get("save_kwargs", {}))
    save_figure(fig = fig, save = save, keyword=mode, save_kwargs=save_kwargs, log = log, verbose=verbose)

    log.write("Finished creating plots.", verbose=verbose)
    return fig, ax
