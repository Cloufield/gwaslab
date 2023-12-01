import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import gwaslab as gl
from pyensembl import EnsemblRelease
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gc as garbage_collect
from adjustText import adjust_text
from gwaslab.g_Log import Log
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_get_sig import annogene
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_recombination_rate
from gwaslab.bd_common_data import get_gtf
from gwaslab.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz_aux_quickfix import _quick_fix
#from gwaslab.quickfix import _quick_fix_p
#from gwaslab.quickfix import _quick_fix_mlog10p
#from gwaslab.quickfix import _quick_fix_chr
#from gwaslab.quickfix import _quick_fix_pos
#from gwaslab.quickfix import _quick_fix_eaf
from gwaslab.viz_aux_quickfix import _get_largenumber
from gwaslab.viz_aux_quickfix import _quick_add_tchrpos
from gwaslab.viz_aux_quickfix import _quick_merge_sumstats
from gwaslab.viz_aux_quickfix import _quick_assign_i
from gwaslab.viz_aux_quickfix import _quick_extract_snp_in_region
from gwaslab.viz_aux_quickfix import _quick_assign_highlight_hue_pair
from gwaslab.viz_aux_quickfix import _quick_assign_marker_relative_size
from gwaslab.viz_aux_annotate_plot import annotate_pair
from gwaslab.io_to_pickle import load_pickle
from gwaslab.io_to_pickle import load_data_from_pickle
from gwaslab.g_Sumstats import Sumstats
from gwaslab.viz_aux_save_figure import save_figure

def plot_miami( 
          path1,
          path2,
          cols1=None,
          cols2=None,
          id0="TCHR+POS",
          id1=None,
          id2=None,
          sep=None,
          mode="txt",
          chr_dict  = None,
          chr_dict1 = False,
          chr_dict2 = False,
          scaled=False,
          scaled1=False,
          scaled2=False,
          region=None,
          region_step = 21,
          region_grid = False,
          region_grid_line = None,
          region_lead_grid = True,
          region_lead_grid_line = None,
          anno = None,
          anno_set=None,
          anno_set1=None,
          anno_set2=None,
          anno_alias1=None,
          anno_alias2=None,
          anno_d1=None,
          anno_d2=None,
          anno_args=None,
          anno_style="right",
          anno_fixed_arm_length=None,
          anno_source = "ensembl",
          anno_max_iter=100,
          anno_adjust=False,
          arm_offset=50,
          arm_scale=1,
          arm_scale_d=None,
          highlight  = None,
          highlight1 = None,
          highlight2 = None,
          highlight_color="#CB132D",
          highlight_windowkb = 500,
          pinpoint=None,
          pinpoint1=None,
          pinpoint2=None,
          pinpoint_color ="red",
          titles=None,
          titles_pad=None, 
          cut=0,
          skip=0,
          build="19",
          cutfactor=10,
          readcsv_args=None,
          cut_line_color="#ebebeb",  
          sig_line=True,
          sig_level=5e-8,
          sig_level_lead=5e-8,
          sig_line_color="grey",
          suggestive_sig_line=False,
          suggestive_sig_level=5e-6,
          suggestive_sig_line_color="grey",
          sc_linewidth=2,
          additional_line = None,
          additional_line_color = None,
          region_hspace = 0.1,
          windowsizekb=500,
          dpi=100,
          figargs=None,
          fontsize = 10,
          font_family="Arial",
          colors1=None,
          colors2=None,
          scatter_kwargs=None,
          xlabel_coords=(-0.017, -0.027),
          xtick_label_pad=None,
          marker_size=(5,25),
          verbose=True,
          xtickpad=None,
          repel_force=0.03,
          save=None,
          save_args=None,
          log=Log()
          ):
    ## figuring arguments ###########################################################################################################
    if cols1 is None:
        cols1 = ["CHR","POS","P"]
    if cols2 is None:
        cols2 = ["CHR","POS","P"]
    if highlight is None:
        highlight  = list()
    if highlight1 is None:
        highlight1 = list()
    if highlight2 is None:
        highlight2 = list()
    if pinpoint is None:
        pinpoint= list()
    if pinpoint1 is None:
        pinpoint1 = list()
    if pinpoint2 is None:
        pinpoint2 = list()
    if sep is None:
        sep=["\t","\t"]
    if chr_dict is None:
        chr_dict  = get_chr_to_number()
    if anno_set is None:
        anno_set= list()
    if anno_set1 is None:
        anno_set1=list()
    if anno_set2 is None:
        anno_set2=list()
    if anno_alias1 is None:
        anno_alias1=dict()
    if anno_alias2 is None:
        anno_alias2=dict()
    if anno_d1 is None:
        anno_d1= dict()
    if anno_d2 is None:
        anno_d2= dict()
    if anno_args is None:
        anno_args= dict()
    if readcsv_args is None:
        readcsv_args={}
    if scatter_kwargs is None:
        scatter_kwargs={}
    if figargs is None:
        figargs= {"figsize":(15,5),"dpi":100}
    if save_args is None:
        save_args={"dpi":100,"facecolor":"white"}
    if colors1 is None:
        colors1=["#597FBD","#74BAD3"]
    if colors2 is None:
        colors2=["#597FBD","#74BAD3"]
    if region_lead_grid_line is None:
        region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}
    if region_grid_line is None:
        region_grid_line = {"linewidth": 2,"linestyle":"--"}
    if titles is None:
        titles=["",""]
    if titles_pad is None:
        titles_pad=[0.2,0.2]
    if type(mode) is str:
        modes =[ mode, mode]
    else:
        modes = mode
    if save is not None:
        if type(save) is not bool:
            if len(save)>3:
                if save[-3:]=="pdf" or save[-3:]=="svg":
                    figargs["dpi"]=72
                    scatter_kwargs["rasterized"]=True
    # id1 
    if id1 is None:
        id1 = id0
    else:
        cols1.append(id1)

    #id2
    if id2 is None:
        id2 = id0
    else:
        cols2.append(id2)
    
    if verbose: log.write("Start to plot miami plot with the following basic settings:")
    if verbose: log.write(" -Genome-wide significance level is set to "+str(sig_level)+" ...")
    if len(anno_set)>0 :
        if verbose: log.write(" -Variants to annotate : ", anno_set)    
    if len(highlight)>0 or len(highlight1)>0 or len(highlight2)>0 :
        if verbose: log.write(" -Loci to highlight : ", highlight,highlight1,highlight2)    
        if verbose: log.write(" -Highlight_window is set to: ", highlight_windowkb, " kb")  
    if len(pinpoint)>0 or len(pinpoint1)>0 or len(pinpoint2)>0 :
        if verbose: log.write(" -Variants to pinpoint : ",pinpoint,pinpoint1,pinpoint2)      
    if region is not None:
        if verbose: log.write(" -Region to plot : chr"+str(region[0])+":"+str(region[1])+"-"+str(region[2])+".")  
    if dpi!=100:
        figargs["dpi"] = dpi
    if xtickpad is None:
        if "figsize" not in figargs.keys():
            figargs["figsize"] = (15,5)
        xtickpad =   figargs["figsize"][1] * region_hspace *72 / 6
    if xtick_label_pad is None:
        if "figsize" not in figargs.keys():
            figargs["figsize"] = (15,5)
        xtick_label_pad =  72 * figargs["figsize"][1] * region_hspace / 6
    
    if additional_line is None:
        lines_to_plot = pd.Series([sig_level, suggestive_sig_level] )
    else:
        lines_to_plot = pd.Series([sig_level, suggestive_sig_level] + additional_line ) 
        if additional_line_color is None:
            additional_line_color = ["grey"]
    
    lines_to_plot = -np.log10(lines_to_plot)

    if chr_dict1==False:
        chr_dict1 = chr_dict
    if chr_dict2==False:
        chr_dict2 = chr_dict
    if scaled == True:
        scaled1 = True
        scaled2 = True
    chrom = "CHR"
    pos="POS"
    
    ## load sumstats1 ###########################################################################################################
    if type(path1) is str:
        if verbose: log.write(" -Loading sumstats1 ({} mode):".format(modes[0]) + path1)
    if verbose: log.write(" -Sumstats1 CHR,POS,P information will be obtained from:",cols1)
    
    if modes[0]=="pickle":
        sumstats1 = load_data_from_pickle(path1,usecols=cols1)
    else:
        if type(path1) is Sumstats:
            if verbose: log.write(" -Loading sumstats1 from gwaslab.Sumstats Object")
            sumstats1 = path1.data[cols1].copy()
        elif type(path1) is pd.DataFrame:
            if verbose: log.write(" -Loading sumstats1 from pandas.DataFrame Object")
            sumstats1 = path1[cols1].copy()
        else:
            sumstats1=pd.read_table(path1,sep=sep[0],usecols=cols1,dtype={cols1[0]:"string",cols1[1]:"Int64",cols1[2]:"float64"},**readcsv_args)

    ## load sumstats2 ###########################################################################################################
    if type(path2) is str:
        if verbose: log.write(" -Loading sumstats2 ({} mode):".format(modes[1]) + path2)
    if verbose: log.write(" -Sumstats2 CHR,POS,P information will be obtained from:",cols2)
    if modes[1]=="pickle":
        sumstats2 = load_data_from_pickle(path2,usecols=cols2) 
    else:
        if type(path2) is Sumstats:
            if verbose: log.write(" -Loading sumstats2 from gwaslab.Sumstats Object")
            sumstats2 = path2.data[cols2].copy()
        elif type(path2) is pd.DataFrame:
            if verbose: log.write(" -Loading sumstats2 from pandas.DataFrame Object")
            sumstats2 = path2[cols2].copy()
        else:
            sumstats2=pd.read_table(path2,sep=sep[0],usecols=cols2,dtype={cols2[0]:"string",cols2[1]:"Int64",cols2[2]:"float64"},**readcsv_args)

    sumstats1 = sumstats1.rename(columns={cols1[0]:"CHR",cols1[1]:"POS",cols1[2]:"P"})
    sumstats1 = _quick_fix(sumstats1,chr_dict=chr_dict1, scaled=scaled1, verbose=verbose, log=log)

    sumstats2 = sumstats2.rename(columns={cols2[0]:"CHR",cols2[1]:"POS",cols2[2]:"P"})
    sumstats2 = _quick_fix(sumstats2,chr_dict=chr_dict2, scaled=scaled2, verbose=verbose, log=log)

    # get a large number
    large_number = _get_largenumber(sumstats1["POS"].max(), sumstats2["POS"].max(),log=log)
    
    ## create merge index and then merge
    sumstats1 = _quick_add_tchrpos(sumstats1,large_number=large_number, dropchrpos=False, verbose=verbose, log=log)
    sumstats2 = _quick_add_tchrpos(sumstats2,large_number=large_number, dropchrpos=False, verbose=verbose, log=log)
    if verbose: log.write(" - Merging sumstats using chr and pos...")
    merged_sumstats = _quick_merge_sumstats(sumstats1=sumstats1,sumstats2=sumstats2)
    
    del(sumstats1)
    del(sumstats2)
    garbage_collect.collect()
    
    # merged_sumstats
    # P_1  scaled_P_1  P_2  scaled_P_2  TCHR+POS CHR POS  
    ######################################################################################################################
    #process highlight and pinpoint id
    for i in range(len(highlight)):
        highlight1.append(highlight[i])
        highlight2.append(highlight[i])
        highlight[i] = highlight[i][0] * large_number + highlight[i][1]
    for i in range(len(highlight1)):
        highlight1[i] = highlight1[i][0] * large_number + highlight1[i][1]
    for i in range(len(highlight2)):
        highlight2[i] = highlight2[i][0] * large_number + highlight2[i][1]
    
    for i in range(len(pinpoint)):
        pinpoint1.append(pinpoint[i])
        pinpoint2.append(pinpoint[i])
        pinpoint[i] = pinpoint[i][0] * large_number + pinpoint[i][1] 
    for i in range(len(pinpoint1)):
        pinpoint1[i] = pinpoint1[i][0] * large_number + pinpoint1[i][1]
    for i in range(len(pinpoint2)):
        pinpoint2[i] = pinpoint2[i][0] * large_number + pinpoint2[i][1] 
    
    for i in range(len(anno_set)):
        anno_set1.append(anno_set[i])
        anno_set2.append(anno_set[i])
        anno_set[i] = anno_set[i][0] * large_number + anno_set[i][1]
    for i in range(len(anno_set1)):
        anno_set1[i] = anno_set1[i][0] * large_number + anno_set1[i][1]
    for i in range(len(anno_set2)):
        anno_set2[i] = anno_set2[i][0] * large_number + anno_set2[i][1]

    ## merging ###########################################################################################################
    
    if skip >0:
        sumstats = merged_sumstats.loc[(merged_sumstats["scaled_P_1"]>skip) | ( merged_sumstats["scaled_P_2"]>skip),:]  
    else:
        sumstats = merged_sumstats
    
    if region is not None:
        marker_size=(25,45)
        sumstats = _quick_extract_snp_in_region(sumstats,region, verbose=verbose, log=log)

    # assign i and get tick dictionary
    sumstats, chrom_df = _quick_assign_i(sumstats)
    
    ## assign highlight indicator
    if len(highlight1)>0 or len(highlight2)>0:
        sumstats, to_highlight1, to_highlight2 = _quick_assign_highlight_hue_pair(sumstats = sumstats, 
                                                    highlight1 = highlight1, 
                                                    highlight2 = highlight2, 
                                                    highlight_windowkb = highlight_windowkb,
                                                    verbose=verbose, log=log)
    
    if verbose: log.write("Plotting...")
    ## figure ###########################################################################################################
    
    figargs["figsize"] = (15,10)
    fig, (ax1, ax5) = plt.subplots(2, 1, 
        gridspec_kw={'height_ratios': [1, 1]},**figargs)
    plt.subplots_adjust(hspace=region_hspace)
    
    ###########################################################################################################
    sumstats["s1"] = _quick_assign_marker_relative_size(sumstats["scaled_P_1"])
    sumstats["s2"] = _quick_assign_marker_relative_size(sumstats["scaled_P_2"])
    
    sumstats["chr_hue"]=sumstats[chrom].astype("category")
    

    ##########################################################################################################################
    maxy = max(sumstats["scaled_P_1"].max(skipna=True), sumstats["scaled_P_2"].max(skipna=True))
    maxy1 = sumstats["scaled_P_1"].max(skipna=True)
    maxy5 = sumstats["scaled_P_2"].max(skipna=True)
    
    if cut:
        if cut is True:
            if verbose: log.write(" - Cut Auto mode is activated...")
            if maxy<20:
                if verbose: log.write(" - maxy <20 , no need to cut.")
                cut=0
            else:
                cut = 20
                cutfactor = ( maxy - cut )/5
        if cut:
            if verbose: log.write(" - Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")
            maxticker=int(max(np.round(sumstats["scaled_P_1"].max(skipna=True)) , np.round(sumstats["scaled_P_2"].max(skipna=True))))
            sumstats.loc[sumstats["scaled_P_1"]>cut,"scaled_P_1"] = (sumstats.loc[sumstats["scaled_P_1"]>cut,"scaled_P_1"]-cut)/cutfactor +  cut
            sumstats.loc[sumstats["scaled_P_2"]>cut,"scaled_P_2"] = (sumstats.loc[sumstats["scaled_P_2"]>cut,"scaled_P_2"]-cut)/cutfactor +  cut
            maxy = (maxticker-cut)/cutfactor + cut
            maxy1=( (sumstats["scaled_P_1"].max(skipna=True)) -cut)/cutfactor + cut
            maxy5=( (sumstats["scaled_P_2"].max(skipna=True)) -cut)/cutfactor + cut
            lines_to_plot[lines_to_plot>cut] = (lines_to_plot[lines_to_plot>cut]-cut)/cutfactor+cut
    ##########################################################################################################################
    legend=None
    style=None
    linewidth=0
    
    
    palette = sns.color_palette(colors1,n_colors=sumstats[chrom].nunique())  
    plot1 = sns.scatterplot(data=sumstats, x='i', y='scaled_P_1',
        hue='chr_hue',
        palette= palette,
        legend=legend,
        style=style,
        size="s1",
        sizes=marker_size,
        linewidth=linewidth,
        zorder=2,
        ax=ax1,
        edgecolor="black",
        **scatter_kwargs)     
    palette = sns.color_palette(colors2,n_colors=sumstats[chrom].nunique())  
    plot2 = sns.scatterplot(data=sumstats, x='i', y='scaled_P_2',
        hue='chr_hue',
        palette= palette,
        legend=legend,
        style=style,
        size="s2",
        sizes=marker_size,
        linewidth=linewidth,
        zorder=2,
        ax=ax5,
        edgecolor="black",
        **scatter_kwargs) 
    
    highlight_i = pd.DataFrame()
    if len(highlight1)>0 :
        if len(to_highlight1)>0:
            if verbose: log.write(" -Highlighting target loci for sumstats1...")

            sns.scatterplot(data=sumstats.loc[sumstats["HUE1"]=="0"], x='i', y='scaled_P_1',
                   hue="HUE1",
                   palette={"0":highlight_color},
                   legend=legend,
                   style=style,
                   size="s1",
                   sizes=(marker_size[0]+1,marker_size[1]+1),
                   linewidth=linewidth,
                   zorder=3,ax=ax1,edgecolor="black",**scatter_kwargs)  
            highlight_i = sumstats.loc[sumstats["TCHR+POS"].isin(highlight),"i"].values
            
    if len(highlight2)>0 :
        if len(to_highlight2)>0:
            if verbose: log.write(" -Highlighting target loci for sumstats2.")
            sns.scatterplot(data=sumstats.loc[sumstats["HUE2"]=="0"], x='i', y='scaled_P_2',
                   hue="HUE2",
                   palette={"0":highlight_color},
                   legend=legend,
                   style=style,
                   size="s2",
                   sizes=(marker_size[0]+1,marker_size[1]+1),
                   linewidth=linewidth,
                   zorder=3,ax=ax5,edgecolor="black",**scatter_kwargs)  

            highlight_i = np.append( highlight_i, sumstats.loc[sumstats["TCHR+POS"].isin(highlight),"i"].values)
        
    if len(pinpoint1)>0 or len(pinpoint2)>0: 
        if sum(sumstats["TCHR+POS"].isin(pinpoint1))>0:  
            to_pinpoint = sumstats.loc[sumstats["TCHR+POS"].isin(pinpoint1),:]
            if verbose: log.write(" -Pinpointing target vairants...")
            ax1.scatter(to_pinpoint["i"],to_pinpoint["scaled_P_1"],color=pinpoint_color,zorder=3,s=marker_size[1]+1)
        else:
            if verbose: log.write(" -Target vairants to pinpoint were not found. Skip pinpointing process...")
                
        if sum(sumstats["TCHR+POS"].isin(pinpoint2))>0:  
            to_pinpoint = sumstats.loc[sumstats["TCHR+POS"].isin(pinpoint2),:]
            if verbose: log.write(" -Pinpointing target vairants...")
            ax5.scatter(to_pinpoint["i"],to_pinpoint["scaled_P_2"],color=pinpoint_color,zorder=3,s=marker_size[1]+1)
        else:
            if verbose: log.write(" -Target vairants to pinpoint were not found. Skip pinpointing process...")
        
        
    ## X #######################################################################################################################
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    
    ax1.set_xlabel("")
    ax1.set_xticks(chrom_df)
    ax1.set_xticklabels(chrom_df.index.map(get_number_to_chr()),fontsize=fontsize,family=font_family) 
    
    ax5.set_ylim(ax1.get_ylim())
    ax5.set_xlim(ax1.get_xlim())
    ax5.set_xlabel("")
    ax5.set_xticks(chrom_df)
    ax5.set_xticklabels([])
    ax5.xaxis.set_ticks_position("top")
    
    ## Y #######################################################################################################################
    #sigline = ax1.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
    #sigline = ax5.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
     
    for ax in [ax1,ax5]:
        if cut == 0: 
            ax.set_ylim(skip,maxy*1.2)
        if cut!= 0:
            cutline = ax.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker],fontsize=fontsize,family=font_family)
            else:
                ax.set_yticks([x for x in range(skip,cut+1,2)])
                ax.set_yticklabels([x for x in range(skip,cut+1,2)],fontsize=fontsize,family=font_family)
            ax.set_ylim(bottom = skip)
            
     
    ### Spines ####################################################################################################################
    

    if anno is not None:
        if len(anno_set1)>0 or len(anno_set2)>0:
            if len(anno_set1)>0:
                if id1 is None:
                    to_annotate1=sumstats.loc[sumstats[id0].isin(anno_set1),:]
                else:
                    to_annotate1=sumstats.loc[sumstats[id1].isin(anno_set1),:]
            else:
                to_annotate1 = pd.DataFrame()
            if len(anno_set2)>0:
                if id2 is None:
                    to_annotate5=sumstats.loc[sumstats[id0].isin(anno_set2),:]
                else:
                    to_annotate5=sumstats.loc[sumstats[id2].isin(anno_set2),:]
            else:
                to_annotate5 = pd.DataFrame()
        else:
            to_annotate1 = getsig(sumstats.loc[sumstats["scaled_P_1"]> float(-np.log10(sig_level_lead)),:],
                        "TCHR+POS",
                        "CHR",
                        "POS",
                        "P_1",
                        build=build,
                        source=anno_source,
                        windowsizekb=windowsizekb,
                        verbose=False,
                        sig_level=sig_level)
            
            if verbose: log.write(" -Keeping only the most significant variants to annotate for each position for sumstast #1...")
            to_annotate1 = to_annotate1.sort_values("scaled_P_1").drop_duplicates(subset=["TCHR+POS"]).sort_values("i")
            
            to_annotate5 = getsig(sumstats.loc[sumstats["scaled_P_2"]> float(-np.log10(sig_level_lead)),:],
                        "TCHR+POS",
                        "CHR",
                        "POS",
                        "P_2",
                        build=build,
                            source=anno_source,
                        windowsizekb=windowsizekb,
                        verbose=False,
                        sig_level=sig_level)
            

            if verbose: log.write(" -Keeping only the most significant variants to annotate for each position for sumstast #2...")
            to_annotate5 = to_annotate5.sort_values("scaled_P_2").drop_duplicates(subset=["TCHR+POS"]).sort_values("i")
        #######################################################################################
        if type(anno) is str:
            if (to_annotate1.empty is False) and anno=="GENENAME":
                    to_annotate1 = annogene(to_annotate1,
                                        id=id0,
                                        chrom="CHR",
                                        pos="POS",
                                        log=log,
                                        build=build,
                                        source=anno_source,
                                        verbose=verbose).rename(columns={"GENE":"GENENAME"})

            if (to_annotate5.empty is False) and anno=="GENENAME":
                    to_annotate5 = annogene(to_annotate5,
                                        id=id0,
                                        chrom="CHR",
                                        pos="POS",
                                        log=log,
                                        build=build,
                                        source=anno_source,
                                        verbose=verbose).rename(columns={"GENE":"GENENAME"})        
    else:
        to_annotate1 = pd.DataFrame()
        to_annotate5 = pd.DataFrame()
            
####################################################################################################################
# Add Annotation to manhattan plot #######################################################
    ax1, ax5 = annotate_pair(
                                sumstats=sumstats,
                                anno=anno,
                                ax1=ax1,
                                ax5=ax5,
                                highlight_i=highlight,
                                to_annotate1=to_annotate1,
                                to_annotate5=to_annotate5,
                                anno_d1=anno_d1,
                                anno_d2=anno_d2,
                                anno_alias1=anno_alias1,
                                anno_alias2=anno_alias2,
                                anno_style=anno_style,
                                anno_args=anno_args,
                                anno_max_iter=anno_max_iter,
                                anno_adjust=anno_adjust,
                                arm_scale=arm_scale,
                                arm_scale_d=arm_scale_d,
                                anno_fixed_arm_length=anno_fixed_arm_length,
                                arm_offset=arm_offset,
                                maxy1=maxy1,
                                maxy5=maxy5,
                                fontsize=fontsize,
                                font_family=font_family,
                                region=region,
                                skip=skip,
                                chrom="CHR",
                                pos="POS",
                                repel_force=repel_force,
                                verbose=verbose,
                                log=log)

####################################################################################################################    
# Adjust the visibility for spines #######################################################
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(True)
    ax1.spines["bottom"].set_visible(True)
    
    ax5.spines["top"].set_visible(True)
    ax5.spines["right"].set_visible(False)
    ax5.spines["left"].set_visible(True)
    ax5.spines["bottom"].set_visible(False)
####################################################################################################################

    if region is not None:
        most_left_snp = sumstats["i"].idxmin()            
        gene_track_offset = sumstats.loc[most_left_snp,pos]-region[1]
        gene_track_start_i = sumstats.loc[most_left_snp,"i"] - gene_track_offset - region[1]
        lead_id_1=sumstats["scaled_P_1"].idxmax()
        lead_id_2=sumstats["scaled_P_2"].idxmax()
        lead_snp_i_1 = sumstats.loc[lead_id_1,"i"]
        lead_snp_i_2 = sumstats.loc[lead_id_2,"i"]
        
        
        region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 
        
        if region_grid is True:
            for i in np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step):
                ax1.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
                ax5.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)

        if region_lead_grid is True:
            ax1.vlines(x=lead_snp_i_1, ymin=0,  ymax=maxy1, zorder=1000,**region_lead_grid_line)
            ax5.vlines(x=lead_snp_i_2, ymin=0,  ymax=maxy5, zorder=1000,**region_lead_grid_line)
            # set x ticks m plot
        ax1.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
        ax5.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
        ax1.set_xticklabels(region_ticks,rotation=30,fontsize=fontsize,family=font_family)
        fig.tight_layout()

        ax1.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])
        ax5.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])

####################################################################################################################
    # set labels
    ax1.set_xlabel("Chromosome",fontsize=fontsize,family=font_family)
    ax1.xaxis.set_label_coords(xlabel_coords[0],xlabel_coords[1])
    
    ax1.tick_params(axis='x', which='major', pad=xtick_label_pad)
    
    ax1.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family=font_family)
    ax5.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family=font_family)
    
    ax1.set_title(titles[0],y=1+titles_pad[0],family=font_family)
    ax5.set_title(titles[1],y=-titles_pad[1],family=font_family)

    # genomewide significant line
    for ax in [ax1, ax5]:
        if sig_line is True:
            sigline = ax.axhline(y=lines_to_plot[0], 
                                    linewidth = sc_linewidth,
                                    linestyle="--",
                                    color=sig_line_color,
                                    zorder=1)
        if suggestive_sig_line is True:
            suggestive_sigline = ax.axhline(y=lines_to_plot[1], 
                                                linewidth = sc_linewidth, 
                                                linestyle="--", 
                                                color=suggestive_sig_line_color,
                                                zorder=1)
        if additional_line is not None:
            for index, level in enumerate(lines_to_plot[2:].values):
                ax.axhline(y=level, 
                            linewidth = sc_linewidth, 
                            linestyle="--", 
                            color=additional_line_color[index%len(additional_line_color)],
                            zorder=1)

    ax5.invert_yaxis() 
    #return fig
    #if save is not None:
    #    if verbose: log.write("Saving plot:")
    #    if save==True:
    #        fig.savefig("./miami_plot.png",bbox_inches="tight",**save_args)
    #        log.write(" -Saved to "+ "./miami_plot.png" + " successfully!" )
    #    else:
    #        fig.savefig(save,bbox_inches="tight",**save_args)
    #        log.write(" -Saved to "+ save + " successfully!" )
    save_figure(fig, save, keyword="miami",save_args=save_args, log=log, verbose=verbose)

    garbage_collect.collect()
    if verbose: log.write("Finished creating miami plot successfully")
# Return matplotlib figure object #######################################################################################
    return fig, log
    
    
        