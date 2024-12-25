import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from gwaslab.g_Log import Log
from gwaslab.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz_plot_mqqplot import _process_xtick
from gwaslab.viz_plot_mqqplot import _process_xlabel
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.util_in_filter_value import _filter_region
from gwaslab.io_process_args import _extract_kwargs

def _plot_cs(pipcs,
            region, 
            figax=None, 
            _posdiccul=None,
            xtick_chr_dict=None,
            pip="PIP",
            onlycs=False,
            cs="CREDIBLE_SET_INDEX",
            marker_size=(45,85),
            fontsize = 12,
            font_family = "Arial",
            legend_title="Credible sets",
            log=Log(),
            verbose=True,
            **kwargs):
        '''
        pipcs : a DataFrame of finemapping results
        '''   
        ## parameters ############################# 
        if xtick_chr_dict is None:         
                xtick_chr_dict = get_number_to_chr()

        scatter_kwargs =   _extract_kwargs("scatter", dict(), locals())

        region_marker_shapes = ['o', '^','s','D','*','P','X','h','8']
        region_ld_colors_m = ["grey","#E51819","green","#F07818","#AD5691","yellow","purple"]
        

        ## filter data #############################
        pipcs = _filter_region(pipcs, region)
        if onlycs ==True:
                pipcs = pipcs.loc[pipcs[cs]>0,:]

        pipcs[cs] = pipcs[cs].astype("string")

        ## figure and ax ############################# 
        if figax is not None:
                ax=figax[1]
                fig=figax[0]
        else:
                fig, ax = plt.subplots()

        # assign i
        pipcs,chrom_df=_quick_assign_i_with_rank(pipcs,  chrpad=0.00, 
                                                use_rank=False, 
                                                chrom="CHR",pos="POS",
                                                drop_chr_start=False,
                                                _posdiccul=_posdiccul)
        pipcs = pipcs.sort_values(by=cs,ascending=True)

        ## plot ##########################################
        scatter_kwargs["markers"]= {m:region_marker_shapes[i] for i,m in enumerate(pipcs[cs].unique())}
        palette = sns.color_palette(region_ld_colors_m,n_colors=pipcs[cs].nunique()) 
        edgecolor="none"

        plot = sns.scatterplot(data=pipcs,
                        x="i",
                        y=pip,
                        hue=cs,
                        edgecolor=edgecolor,
                        palette=palette, 
                        style=cs,
                        s=marker_size[1],
                        ax=ax,
                        **scatter_kwargs)

        # process legend
        handles, labels = ax.get_legend_handles_labels()
        new_labels = []
        new_handles = []
        ncol = len(labels)

        for i,label in enumerate(labels):
                if label in [str(j) for j in range(1,10)]:
                        new_labels.append(labels[i])
                        new_handles.append(handles[i])
        
        ax.legend(labels =new_labels,  
                  handles=new_handles, 
                  loc="upper right", 
                  bbox_to_anchor=(0.995, 0.995), 
                  ncol=1, 
                  scatterpoints=2, 
                  title=legend_title, 
                  frameon=True)

        return fig, log