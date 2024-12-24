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
            xlim_i=None,
            fig=None, 
            ax=None,
            xtick_chr_dict=None,
            pip="PIP",
            cs="CREDIBLE_SET_INDEX",
            fontsize = 12,
            font_family = "Arial",
            log=Log(),
            verbose=True,
            region_step=10):
        '''
        pipcs : a DataFrame of finemapping results
        '''    
        if xtick_chr_dict is None:         
                xtick_chr_dict = get_number_to_chr()
        scatter_args=dict()
        region_marker_shapes = ['o', '^','s','D','*','P','X','h','8']

        if xlim_i is None:
                xlim_i=(0,0)
        target_chr = region[0]
        target_start = region[1]
        target_end = region[2]

        ## filter
        pipcs = _filter_region(pipcs, region)
        pipcs[cs] = pipcs[cs].astype("string")

        pipcs,chrom_df=_quick_assign_i_with_rank(pipcs,  chrpad=0, 
                                                use_rank=False, 
                                                chrom="CHR",pos="POS",
                                                drop_chr_start=False,
                                                _posdiccul=None)

        fig, ax = plt.subplots()

        pipcs = pipcs.sort_values(by=cs,ascending=True)

        scatter_args["markers"]= {m:region_marker_shapes[i] for i,m in enumerate(pipcs[cs].unique())}
        print(scatter_args)
        style=cs



        sns.scatterplot(data=pipcs,
                        x="i",
                        y=pip,
                        hue=cs,
                        style=style,
                        ax=ax,
                        **scatter_args)
    
        most_left_snp      = pipcs["i"].idxmin()
        # distance between leftmost variant position to region left bound
        gap_length  = pipcs.loc[most_left_snp,"POS"] - region[1]
        # rebase i to region[1] : the i value when POS=0
        i_to_pos_offset = pipcs.loc[most_left_snp,"i"] - gap_length - region[1]
    
        ax.set_xlim([i_to_pos_offset+target_start,i_to_pos_offset+target_end])
    
        region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 

        ax.set_xticks(np.linspace(i_to_pos_offset+region[1], i_to_pos_offset+region[2], num=region_step))
        ax.set_xticklabels(region_ticks,rotation=45,fontsize=fontsize,family="sans-serif")

        ax3=None
        ax, ax3 = _process_xlabel(region=region, 
                                                xlabel=None, 
                                                ax1=ax, 
                                                gtf_path=None, 
                                                mode="r", 
                                                fontsize=fontsize, 
                                                font_family=font_family,  
                                                ax3=ax3,
                                                log=log, 
                                                verbose=verbose) 

        return ax