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


def _plot_cs(pipcs,
            region, 
            figax=None, 
            _posdiccul=None,
            xtick_chr_dict=None,
            pip="PIP",
            cs="CREDIBLE_SET_INDEX",
            marker_size=(45,65),
            fontsize = 12,
            font_family = "Arial",
            log=Log(),
            verbose=True,
            region_step=10,
            **kwargs):
        '''
        pipcs : a DataFrame of finemapping results
        '''   
        ## parameters ############################# 
        if xtick_chr_dict is None:         
                xtick_chr_dict = get_number_to_chr()
        
        scatter_args=dict()
        region_marker_shapes = ['o', '^','s','D','*','P','X','h','8']
        
        ## filter data #############################
        pipcs = _filter_region(pipcs, region)
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
        scatter_args["markers"]= {m:region_marker_shapes[i] for i,m in enumerate(pipcs[cs].unique())}

        sns.scatterplot(data=pipcs,
                        x="i",
                        y=pip,
                        hue=cs,
                        style=cs,
                        sizes=marker_size,
                        ax=ax,
                        **scatter_args)

        return fig, log