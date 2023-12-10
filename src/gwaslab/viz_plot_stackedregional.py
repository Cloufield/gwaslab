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

def plot_stacked_region(paths,
                        vcfs,
                        mode="r",
                        region=None,
                        gtf=None,
                        ratio=1,
                        fig_args=None,
                        region_hspace=0.04
                        ):
    

    if fig_args is None:
        fig_args = {}
    if len(vcfs)==1:
        vcfs = vcfs *len(paths)
    # create figure and axes
    number_of_plot = len(paths)

    number_of_plot +=1
    fig_args["figsize"] = [7*number_of_plot , 10]
    fig, axes = plt.subplots(number_of_plot, 1, sharex=True, 
                             gridspec_kw={'height_ratios': [1 for i in range(number_of_plot-1)]+[ratio]},
                             **fig_args)
    
    #plt.subplots_adjust(hspace=region_hspace)
    
    # load sumstats
    sumstats_list = paths
    
    # plot manhattan plot
    for index,sumstats in enumerate(sumstats_list):
        fig,log = mqqplot(sumstats,
                        chrom="CHR",
                        pos="POS",
                        p="P",
                        region=region,
                        mlog10p="MLOG10P",
                        snpid="SNPID",
                        vcf_path=vcfs[index],
                        mode=mode,
                        gtf_path=None,
                        figax=(fig,axes[index],axes[-1]),
                        _invert=False,
                        _if_quick_qc=False
                        )
        if index==len(sumstats_list)-1:
            # plot last m and gene track 
            fig,log = mqqplot(sumstats,
                            chrom="CHR",
                            pos="POS",
                            p="P",
                            region=region,
                            mlog10p="MLOG10P",
                            snpid="SNPID",
                            vcf_path=vcfs[index],
                            mode=mode,
                            gtf_path="default",
                            figax=(fig,axes[index],axes[-1]),
                            _invert=False,
                            _if_quick_qc=False
                            )            
    
    # adjust labels

    # adjust lines
    