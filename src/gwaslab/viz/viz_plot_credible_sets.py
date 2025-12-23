import copy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.info.g_Log import Log
from gwaslab.io.io_process_kwargs import _extract_kwargs
from gwaslab.util.util_in_filter_value import _filter_region
from gwaslab.viz.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_plot_mqqplot import _process_xlabel, _process_xtick

def _plot_cs(pipcs_raw,
            region=None, 
            locus=None,
            figax=None, 
            _posdiccul=None,
            xtick_chr_dict=None,
            pip="PIP",
            onlycs=False,
            pos="POS",
            chrom="CHR",
            cs="CREDIBLE_SET_INDEX",
            cs_category = "CS_CATEGORY",
            marker_size=(45,85),
            fontsize = 12,
            font_family = "Arial",
            legend_title="Credible set",
            fig_kwargs=None,
            log=Log(),
            verbose=True,
            **kwargs):
        '''
        pipcs : a DataFrame of finemapping results
        '''
        # Extract dataframe if Sumstats object is passed
        if hasattr(pipcs_raw, 'data') and not isinstance(pipcs_raw, pd.DataFrame):
            pipcs_raw = pipcs_raw.data
        
        pipcs = pipcs_raw.copy()
        ## parameters ############################# 
        if xtick_chr_dict is None:         
                xtick_chr_dict = get_number_to_chr()
                
        style = set_plot_style(
                plot="plot_pipcs",
                fig_kwargs=fig_kwargs,
                fontsize=fontsize,
                fontfamily=font_family,
                verbose=verbose,
                log=log,
        )
        fig_kwargs = style["fig_kwargs"]
        fontsize = style["fontsize"]
        font_family = style["font_family"]
        scatter_kwargs =   _extract_kwargs("scatter", dict(), locals())

        region_marker_shapes = ['o', '^','s','D','*','P','X','h','8']
        region_ld_colors_m = ["grey","#E51819","green","#F07818","#AD5691","yellow","purple"]
        
        if region is not None:
                ## filter data #############################
                pipcs = _filter_region(pipcs, region)
                log.write(" -Loading PIP and CS for variants in the region :{}".format(region))
        
        if locus is not None:
                pipcs = pipcs.loc[pipcs["LOCUS"] == locus,:].copy()
                log.write(" -Loading PIP and CS for variants in the locus :{}".format(region))
                if region is None:
                        region = (pipcs[chrom].iloc[0],pipcs[pos].min(),pipcs[pos].max())
                        log.write(" -Extracted region:{}".format(region))
        if onlycs ==True:
                pipcs = pipcs.loc[pipcs[cs]>0,:]
                log.write(" -Loading only variants in CS...")

        pipcs[cs] = pipcs[cs].astype("string")

        ## figure and ax ############################# 
        if figax is not None:
                ax=figax[1]
                fig=figax[0]
        else:
                fig, ax = plt.subplots(**fig_kwargs)

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
        
        if cs_category in pipcs.columns:
                cs_category_dic = pipcs.loc[~pipcs[cs_category].isna(), [cs, cs_category]].drop_duplicates().set_index(cs).to_dict()
                

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
        
        region_step=21
        region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 
        
        most_left_snp      = pipcs["i"].idxmin()
        # distance between leftmost variant position to region left bound
        i_pos_offset = pipcs.loc[most_left_snp,"i"] - pipcs.loc[most_left_snp,pos]
        ax.set_xticks(np.linspace(i_pos_offset+region[1], i_pos_offset+region[2], num=region_step))
        ax.set_xticklabels(region_ticks,rotation=45)
        xlabel = "Chromosome "+str(region[0])+" (MB)"
        ax.set_xlabel(xlabel,fontsize=fontsize,family=font_family)

        # process legend
        handles, labels = ax.get_legend_handles_labels()
        new_labels = []
        new_handles = []
        ncol = len(labels)

        ax.tick_params(axis='y', 
                        labelsize=fontsize,
                        labelfontfamily=font_family) 
        
        for i,label in enumerate(labels):
                if label in [str(j) for j in range(1,10)]:
                        if cs_category in pipcs.columns:
                                new_labels.append("#{} - {}".format(labels[i],cs_category_dic[cs_category][label]) ) 
                        else:
                                new_labels.append("#"+labels[i])
                        new_handles.append(handles[i])
        
        
        ax.legend(labels =new_labels,  
                  handles=new_handles, 
                  loc="upper right", 
                  bbox_to_anchor=(0.995, 0.995), 
                  ncol=1, 
                  markerfirst=False,
                  scatterpoints=1, 
                  title=legend_title, 
                  title_fontproperties={"size":fontsize,"family":font_family},
                  prop={"size":fontsize,"family":font_family},
                  frameon=False)

        return fig, log
