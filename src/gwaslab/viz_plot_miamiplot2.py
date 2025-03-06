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
from gwaslab.viz_aux_quickfix import _get_largenumber
from gwaslab.viz_aux_quickfix import _quick_add_tchrpos
from gwaslab.viz_aux_quickfix import _quick_merge_sumstats
from gwaslab.viz_aux_quickfix import _quick_assign_i
from gwaslab.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz_aux_quickfix import _quick_extract_snp_in_region
from gwaslab.viz_aux_quickfix import _quick_assign_highlight_hue_pair
from gwaslab.viz_aux_quickfix import _quick_assign_marker_relative_size
from gwaslab.viz_aux_annotate_plot import annotate_pair
from gwaslab.io_to_pickle import load_pickle
from gwaslab.io_to_pickle import load_data_from_pickle
from gwaslab.g_Sumstats import Sumstats
from gwaslab.viz_aux_save_figure import save_figure
from gwaslab.viz_plot_mqqplot import mqqplot
from gwaslab.g_version import _get_version

def plot_miami2( 
          path1=None,
          path2=None,
          merged_sumstats=None,
          suffixes = None,
          cols=None,
          cols1=None,
          cols2=None,
          id0="TCHR+POS",
          id1=None,
          id2=None,
          sep=None,
          sep1=None,
          sep2=None,
          mode="m",
          loadmode="tabular",
          loadmode1=None,
          loadmode2=None,
          chr_dict=None,
          chr_dict1=None,
          chr_dict2=None,
          scaled=False,
          scaled1=False,
          scaled2=False,
          titles=None,
          titles_pad=None, 
          use_rank=False,
          chrpad=0.03,
          region_hspace = 0.1,
          dpi=100,
          fontsize = 10,
          font_family="Arial",
          xlabel_coords=(-0.01, -0.027),
          xtick_label_pad=None,
          verbose=True,
          xtickpad=None,
          fig_args=None,
          readcsv_args=None,
          readcsv_args1=None,
          readcsv_args2=None,
          scatter_args=None,
          same_ylim=False,
          save=None,
          save_args=None,
          figax=None,
          log=Log(),
          **mqq_args
          ):
    log.write("Start to create miami plot {}:".format(_get_version()), verbose=verbose)
    ## figuring arguments ###########################################################################################################
    # figure columns to use
    if scaled == True:
        scaled1 = True
        scaled2 = True
    
    if cols is None:
        if scaled == True:
            cols = ["CHR","POS","MLOG10P"]
        else:
            cols = ["CHR","POS","P"]

    if cols1 is None:
        cols1 = cols.copy()
    if cols2 is None:
        cols2 = cols.copy()

    if id1 is not None:
        cols1.append(id1)
    if id2 is not None:
        cols2.append(id2)
    
    if sep is None:
        sep="\t"
    if sep1 is None:
        sep1=sep
    if sep2 is None:
        sep2=sep
    
    if chr_dict is None:
        chr_dict  = get_chr_to_number()
    if chr_dict1 is None:
        chr_dict1 = chr_dict
    if chr_dict2 is None:
        chr_dict2 = chr_dict


    
    if readcsv_args is None:
        readcsv_args={}
    if readcsv_args1 is None:
        readcsv_args1=readcsv_args
    if readcsv_args2 is None:
        readcsv_args2=readcsv_args
    
    if loadmode1 is None:
        loadmode1 = loadmode
    if loadmode2 is None:
        loadmode2 = loadmode

    if scatter_args is None:
        scatter_args={}

    if fig_args is None:
        fig_args= {"figsize":(15,5),"dpi":100}
    if save_args is None:
        save_args={"dpi":100,"facecolor":"white"}

    # figure out mqq args
    mqq_args1,mqq_args2 = _sort_args_to_12(mqq_args)
    
    # figure out args for pdf and svg
    fig_args, scatter_args = _figure_args_for_vector_plot(save, fig_args, scatter_args)

    # add suffix if ids are the same
    id1_1, id2_2, mqq_args1, mqq_args2 = _solve_id_contradictory(id0, id1, id2, mqq_args1, mqq_args2)

    if dpi!=100:
        fig_args["dpi"] = dpi
    if xtickpad is None:
        if "figsize" not in fig_args.keys():
            fig_args["figsize"] = (15,5)
        xtickpad =   fig_args["figsize"][1] * region_hspace *72 / 6
    if xtick_label_pad is None:
        if "figsize" not in fig_args.keys():
            fig_args["figsize"] = (15,5)
        xtick_label_pad =  72 * fig_args["figsize"][1] * region_hspace / 6
    
    if titles is None:
        titles=["",""]
    
    
    titles_pad_adjusted=[1,0]
    if titles_pad is None:
        titles_pad=[0.2,0.2]
        if "anno" in mqq_args.keys():
            if mqq_args["anno"] is not None:
                titles_pad_adjusted= [1 + titles_pad[0], -titles_pad[1]]
        if "anno1" in mqq_args.keys():
            if mqq_args["anno1"] is not None:
                titles_pad_adjusted[0]= 1 + titles_pad[0]
        if "anno2" in mqq_args.keys():
            if mqq_args["anno2"] is not None:
                titles_pad_adjusted[1]=  - titles_pad[1]
    else:
        titles_pad_adjusted=[1 + titles_pad[0], -titles_pad[1]]

    if merged_sumstats is None:
    ## load sumstats1 ###########################################################################################################
        sumstats1 = _figure_type_load_sumstats(name="Sumstats1", 
                                            path=path1, 
                                            sep=sep1, 
                                            cols=cols1, 
                                            readcsv_args=readcsv_args2, 
                                            loadmode=loadmode1, 
                                            log=log,
                                            verbose=verbose)
        ## load sumstats2 ###########################################################################################################
        sumstats2 = _figure_type_load_sumstats(name="Sumstats2", 
                                            path=path2, 
                                            sep=sep2, 
                                            cols=cols2, 
                                            readcsv_args=readcsv_args2, 
                                            loadmode=loadmode2, 
                                            log=log,
                                            verbose=verbose)
    else:
        cols1[2] += suffixes[0]
        cols2[2] += suffixes[1]
        sumstats1 = merged_sumstats[cols1].copy()
        sumstats2 = merged_sumstats[cols2].copy()

    ## rename and quick fix ###########################################################################################################
    renaming_dict1 = {cols1[0]:"CHR",cols1[1]:"POS",cols1[2]:"P"}
    renaming_dict2 = {cols2[0]:"CHR",cols2[1]:"POS",cols2[2]:"P"}

    sumstats1 = sumstats1.rename(columns=renaming_dict1)
    sumstats1 = _quick_fix(sumstats1,chr_dict=chr_dict1, scaled=scaled1, verbose=verbose, log=log)

    sumstats2 = sumstats2.rename(columns=renaming_dict2)
    sumstats2 = _quick_fix(sumstats2,chr_dict=chr_dict2, scaled=scaled2, verbose=verbose, log=log)

    # get a large number ###########################################################################################################
    large_number = _get_largenumber(sumstats1["POS"].max(), sumstats2["POS"].max(),log=log)
    
    ## create merge index ###########################################################################################################
    sumstats1 = _quick_add_tchrpos(sumstats1,large_number=large_number, dropchrpos=False, verbose=verbose, log=log)
    sumstats2 = _quick_add_tchrpos(sumstats2,large_number=large_number, dropchrpos=False, verbose=verbose, log=log)
    log.write(" -Merging sumstats using chr and pos...",verbose=verbose)
    
    ###### merge #####################################################################################################
    merged_sumstats = _quick_merge_sumstats(sumstats1=sumstats1,sumstats2=sumstats2)
    
    #assign index for x axis
    merged_sumstats, chrom_df = _quick_assign_i_with_rank(merged_sumstats, 
                                                          chrpad=chrpad, 
                                                          use_rank=use_rank, 
                                                          chrom="CHR",
                                                          pos="POS",
                                                          drop_chr_start=False)
    
    # P_1  scaled_P_1  P_2  scaled_P_2  TCHR+POS CHR POS    
    log.write(" -Columns in merged sumstats: {}".format(",".join(merged_sumstats.columns)), verbose=verbose)


    del(sumstats1)
    del(sumstats2)
    garbage_collect.collect()
    #####################################################################################################################
    ##plotting
    if figax is None:
        #fig_args["figsize"] = (15,10)
        fig, (ax1, ax5) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 1]},**fig_args)
        plt.subplots_adjust(hspace=region_hspace)    
    else:
        fig, ax1, ax5 = figax

    #if same_ylim==True:
        #maxy = merged_sumstats[["scaled_P_1","scaled_P_2"]].max().max()

    log.write("Start to create Manhattan plot for sumstats1...", verbose=verbose)
    fig,log = mqqplot(merged_sumstats,
                      chrom="CHR",
                      pos="POS",
                      p="P_1",
                      mlog10p="scaled_P_1",
                      snpid=id1_1,
                      scaled=scaled1,
                      log=log, 
                      mode=mode,
                      figax=(fig,ax1),
                      scatter_args=scatter_args,
                      _chrom_df_for_i = chrom_df,
                      _invert=False,
                      _if_quick_qc=False,
                      **mqq_args1
                     )
    log.write("Finished creating Manhattan plot for sumstats1".format(_get_version()), verbose=verbose)

    log.write("Start to create Manhattan plot for sumstats2...", verbose=verbose)
    fig,log = mqqplot(merged_sumstats,
                      chrom="CHR",
                      pos="POS",
                      p="P_2",
                      mlog10p="scaled_P_2",
                      scaled=scaled2,
                      snpid=id2_2,
                      log=log, 
                      mode=mode,
                      figax=(fig,ax5),
                      scatter_args=scatter_args,
                      _chrom_df_for_i = chrom_df,
                       _invert=True,
                      _if_quick_qc=False,
                     **mqq_args2)
    log.write("Finished creating Manhattan plot for sumstats2".format(_get_version()), verbose=verbose)
    

    #####################################################################################################################    
    ax1l, ax1r = ax5.get_xlim()
    ax5l, ax5r = ax1.get_xlim()
    ax1.set_xlim([min(ax1l,ax5l), max(ax1r,ax5r)])
    ax5.set_xlim([min(ax1l,ax5l), max(ax1r,ax5r)])
    #####################################################################################################################
    ax5.set_xlabel("")
    #ax5.set_xticks(chrom_df)
    ax5.set_xticklabels([])
    ax5.xaxis.set_ticks_position("top")

    # Ad#just the visibility for spines #######################################################
    ax1, ax5 = _set_spine_visibility(ax1, ax5)
    ######################################################################################################################
#####################################################################################################################
    # set labels
    ax1.set_xlabel("Chromosome",fontsize=fontsize,family=font_family)
    ax1.xaxis.set_label_coords(xlabel_coords[0],xlabel_coords[1])
    
    ax1.tick_params(axis='x', which='major', pad=xtick_label_pad)
    
    ax1.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family=font_family)
    ax5.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family=font_family)
    
    ax1.set_title(titles[0],y=titles_pad_adjusted[0],family=font_family)
    ax5.set_title(titles[1],y=titles_pad_adjusted[1],family=font_family)

    ax5.invert_yaxis() 

    save_figure(fig, save, keyword="miami",save_args=save_args, log=log, verbose=verbose)

    garbage_collect.collect()
    
    log.write("Finished creating miami plot successfully", verbose=verbose)
    #Return matplotlib figure object #######################################################################################
    return fig, log

def _sort_args_to_12(mqq_args):
    mqq_args1={}
    mqq_args2={}
    for key, value in mqq_args.items():
        if key[-1]=="1":
            # for top mqq
            mqq_args1[key[:-1]] = value
        elif key[-1]=="2":
            # for bottom mqq
            mqq_args2[key[:-1]] = value
        else:
            # for both top and bottom
            mqq_args1[key] = value
            mqq_args2[key] = value
    return mqq_args1, mqq_args2

def _solve_id_contradictory(id0, id1, id2, mqq_args1, mqq_args2):
    if (id1 is not None) and (id2 is not None):
        if id1 == id2:
            id1_1 = id1 + "_1"
            id2_2 = id2 + "_2"
            if "anno" in mqq_args1.keys():
                if mqq_args1["anno"] == id1:
                    mqq_args1["anno"] = id1_1
            if "anno" in mqq_args2.keys():
                if mqq_args2["anno"] == id2:
                    mqq_args2["anno"] = id2_2
        else:
            id1_1 = id1
            id2_2 = id2
    
    if id1 is None:
        id1_1 = id0

    if id2 is None:
        id2_2 = id0

    return (id1_1, id2_2, mqq_args1, mqq_args2)

def _figure_args_for_vector_plot(save, fig_args, scatter_kwargs ):
    if save is not None:
        if type(save) is not bool:
            if len(save)>3:
                if save[-3:]=="pdf" or save[-3:]=="svg":
                    fig_args["dpi"]=72
                    scatter_kwargs["rasterized"]=True
    return fig_args, scatter_kwargs

def _set_spine_visibility(ax1,ax5):
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(True)
    ax1.spines["bottom"].set_visible(True)
    
    ax5.spines["top"].set_visible(True)
    ax5.spines["right"].set_visible(False)
    ax5.spines["left"].set_visible(True)
    ax5.spines["bottom"].set_visible(False)
    return ax1,ax5

def _figure_type_load_sumstats(name, path, sep, cols, readcsv_args, loadmode, log, verbose):
    if type(path) is str:
        log.write(" -Loading {} ({} mode): {}".format(name, loadmode, path), verbose=verbose)
    log.write(" -Obtaining {} CHR, POS, P and annotation from: {}".format(name, cols), verbose=verbose)
    
    if loadmode=="pickle":
        sumstats = load_data_from_pickle(path,usecols=cols) 
    else:
        if type(path) is Sumstats:
            log.write(" -Loading {} from gwaslab.Sumstats Object".format(name), verbose=verbose)
            sumstats = path.data[cols].copy()
        elif type(path) is pd.DataFrame:
            log.write(" -Loading {} from pandas.DataFrame Object".format(name), verbose=verbose)
            sumstats = path[cols].copy()
        else:
            log.write(" -Loading {} from tabular files".format(name), verbose=verbose)
            sumstats=pd.read_table(path,sep=sep,usecols=cols,dtype={cols[0]:"string",cols[1]:"Int64",cols[2]:"float64"},**readcsv_args)
    return sumstats