from typing import TYPE_CHECKING, Optional, Dict, Any, Union, List, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
from gwaslab.info.g_Log import Log
from gwaslab.io.io_process_kwargs import _update_arg
from gwaslab.io.io_process_kwargs import _update_kwargs
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes

def _plot(
    associations: Any,
    values: str = "Beta",
    sort: str = "P-value",
    mode: str = "gwaslab",
    fontsize: float = 12,
    font_family: str = "Arial",
    cmap: Optional[Any] = None,
    ylabel: str = "Y",
    xlabel: str = "X",
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    yticks: Optional[List[float]] = None,
    xticks: Optional[List[float]] = None,
    ytick_labels: Optional[List[str]] = None,
    xtick_labels: Optional[List[str]] = None,
    title: str = "Title",
    title_pad: float = 1.08,
    title_fontsize: float = 13,
    title_kwargs: Optional[Dict[str, Any]] = None,
    linewidth: Optional[float] = None,
    linestyle: Optional[str] = None,
    linecolor: Optional[str] = None,
    dpi: int = 200,
    anno_kwargs: Optional[Dict[str, Any]] = None,
    err_kwargs: Optional[Dict[str, Any]] = None,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    scatter_kwargs: Optional[Dict[str, Any]] = None,
    save: Optional[Union[bool, str]] = None,
    save_kwargs: Optional[Dict[str, Any]] = None,
    log: Log = Log(),
    verbose: bool = True,
    **args: Any
) -> Tuple['Figure', 'Axes']:
    
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
