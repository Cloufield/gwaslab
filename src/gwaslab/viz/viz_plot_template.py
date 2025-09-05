import seaborn as sns
import matplotlib.pyplot as plt
from gwaslab.io.io_process_args import _update_args
from gwaslab.io.io_process_args import _update_arg
from gwaslab.g_Log import Log
from gwaslab.viz.viz_aux_save_figure import save_figure

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
            title="Title",
            title_pad=1.08, 
            title_fontsize=13,
            linewidth=None,
            linestyle=None,
            linecolor=None,
            dpi=200,
            fig_args= None,
            scatter_args=None,
            save=None,
            save_args=None,
            log=Log(),
            verbose=True,
            **args               
            ):
    
    # update args
    
    cmap = _update_arg(cmap, "RdBu")

    fig_args = _update_args(fig_args, dict(figsize=(10,10)))

    # create fig and ax    
    fig,ax = plt.subplots(**fig_args)

    # draw lines
    horizontal_line = ax.axhline(y=1, linewidth = linewidth,
                            linestyle="--",
                            color=linecolor,zorder=1000)
    vertical_line   = ax.axvline(x=1, linewidth = linewidth,
                            linestyle="--",
                            color=linecolor,zorder=1000)

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

    save_figure(fig = fig, save = save, keyword=mode, save_args=save_args, log = log, verbose=verbose)

    log.write("Finished creating plots.", verbose=verbose)
    return fig, ax