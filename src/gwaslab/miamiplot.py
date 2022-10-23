import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import gwaslab as gl
from gwaslab.Log import Log
from gwaslab.getsig import getsig
from gwaslab.getsig import annogene
from gwaslab.CommonData import get_chr_to_number
from gwaslab.CommonData import get_number_to_chr
from gwaslab.CommonData import get_recombination_rate
from gwaslab.CommonData import get_gtf
from pyensembl import EnsemblRelease
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gc

def plot_miami( path1,
          path2,
          cols1=["CHR","POS","P"],
          cols2=["CHR","POS","P"],
          sep=["\t","\t"],
          large_number=10000000000,
          cut=0,
          skip=0,
          cutfactor=10,
          readcsv_args={},
          cut_line_color="#ebebeb",  
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_level=5e-6,
          region_hspace = 0.1,
          figargs= {"figsize":(15,5),"dpi":100},
          fontsize = 10,
          colors=["#597FBD","#74BAD3"],
          marker_size=(5,25),
          verbose=True,
          log=Log()
          ):
    ## load sumstats1 ###########################################################################################################
    if verbose: log.write("Start to plot miami plot with the following basic settings:")
    if verbose: log.write("Loading sumstats1:" + path1)
    if verbose: log.write("Sumstats1 CHR,POS,P information will be obtained from:",cols1)
    sumstats1 = pd.read_csv(path1,sep=sep[0],usecols=cols1,dtype={cols1[0]:"string",cols1[1]:"Int64",cols1[2]:"float64"},**readcsv_args)
    sumstats1 = sumstats1.rename(columns={cols1[0]:"CHR",cols1[1]:"POS",cols1[2]:"P"})
    sumstats1["CHR"] = np.floor(pd.to_numeric(sumstats1["CHR"].map(get_chr_to_number(),na_action="ignore"), errors='coerce')).astype('Int64')
    sumstats1["POS"] = np.floor(pd.to_numeric(sumstats1["POS"], errors='coerce')).astype('Int64')   
    
    ## load sumstats2 ###########################################################################################################
    if verbose: log.write("Loading sumstats2:" + path2)
    if verbose: log.write("Sumstats2 CHR,POS,P information will be obtained from:",cols2)
    sumstats2 = pd.read_csv(path2,sep=sep[1],usecols=cols2,dtype={cols1[0]:"string",cols1[1]:"Int64",cols1[2]:"float64"},**readcsv_args)
    sumstats2 = sumstats2.rename(columns={cols2[0]:"CHR",cols2[1]:"POS",cols2[2]:"P"})
    sumstats2["CHR"] = np.floor(pd.to_numeric(sumstats2["CHR"].map(get_chr_to_number(),na_action="ignore"), errors='coerce')).astype('Int64')
    sumstats2["POS"] = np.floor(pd.to_numeric(sumstats2["POS"], errors='coerce')).astype('Int64')
    
    ## create merge index
    sumstats1["TCHR+POS"]=sumstats1["CHR"]*large_number + sumstats1["POS"]
    sumstats2["TCHR+POS"]=sumstats2["CHR"]*large_number + sumstats2["POS"]
    sumstats1["TCHR+POS"] = sumstats1["TCHR+POS"].astype('Int64')
    sumstats2["TCHR+POS"] = sumstats2["TCHR+POS"].astype('Int64')
    sumstats1.drop(labels=["CHR","POS"],axis=1,inplace=True)
    sumstats2.drop(labels=["CHR","POS"],axis=1,inplace=True)
    sumstats1.dropna(inplace=True)
    sumstats2.dropna(inplace=True)
    ## merging ###########################################################################################################
    if verbose: log.write("Merging sumstats using chr and pos...")
    merged_sumstats = pd.merge(sumstats1,sumstats2,on="TCHR+POS",how="outer",suffixes=('_1', '_2'))
    merged_sumstats["CHR"] = np.divmod(merged_sumstats["TCHR+POS"],large_number)[0]
    merged_sumstats["POS"] = np.divmod(merged_sumstats["TCHR+POS"],large_number*merged_sumstats["CHR"])[1]
    
    del(sumstats1)
    del(sumstats2)
    gc.collect()
    
    chrom = "CHR"
    pos="POS"
    merged_sumstats["scaled_P_1"] = -np.log10(merged_sumstats["P_1"])
    merged_sumstats["scaled_P_2"] = -np.log10(merged_sumstats["P_2"])
    if skip >0:
        sumstats = merged_sumstats.loc[(merged_sumstats["scaled_P_1"]>skip) | ( merged_sumstats["scaled_P_2"]>skip),:]  
    else:
        sumstats = merged_sumstats
        
    sumstats = sumstats.sort_values([chrom,pos])
    sumstats["id"]=range(len(sumstats))
    sumstats=sumstats.set_index("id")

    #create a df , groupby by chromosomes , and get the maximum position
    posdic = sumstats.groupby(chrom)[pos].max()
    # convert to dictionary
    posdiccul = dict(posdic)

    # fill empty chr with 0
    for i in range(0,26):
        if i in posdiccul: continue
        else: posdiccul[i]=0

    # cumulative sum dictionary
    for i in range(2,sumstats[chrom].max()+1):
        posdiccul[i]= posdiccul[i-1] + posdiccul[i] + sumstats[pos].max()*0.05

    # convert base pair postion to x axis position using the cumulative sum dictionary
    sumstats["add"]=sumstats[chrom].apply(lambda x : posdiccul[int(x)-1])
    sumstats["i"]=sumstats[pos]+sumstats["add"]

    #for plot, get the chr text tick position      
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')

    if verbose: log.write("Plotting...")
    ## figure ###########################################################################################################
    
    figargs["figsize"] = (15,10)
    fig, (ax1, ax5) = plt.subplots(2, 1, 
        gridspec_kw={'height_ratios': [1, 1]},**figargs)
    plt.subplots_adjust(hspace=region_hspace)
    
    ###########################################################################################################
    
    sumstats["s1"]=1
    sumstats["s2"]=1
    sumstats.loc[sumstats["scaled_P_1"]>-np.log10(5e-4),"s1"]=2
    sumstats.loc[sumstats["scaled_P_1"]>-np.log10(suggestive_sig_level),"s1"]=3
    sumstats.loc[sumstats["scaled_P_1"]>-np.log10(sig_level),"s1"]=4
    sumstats.loc[sumstats["scaled_P_2"]>-np.log10(5e-4),"s2"]=2
    sumstats.loc[sumstats["scaled_P_2"]>-np.log10(suggestive_sig_level),"s2"]=3
    sumstats.loc[sumstats["scaled_P_2"]>-np.log10(sig_level),"s2"]=4
    
    
    sumstats["chr_hue"]=sumstats[chrom].astype("category")
    palette = sns.color_palette(colors,n_colors=sumstats[chrom].nunique())  
    legend = None
    style=None
    linewidth=0
    ##########################################################################################################################
    maxy = max(sumstats["scaled_P_1"].max(), sumstats["scaled_P_2"].max())
    if cut:
        if cut is True:
            if verbose: log.write(" -Cut Auto mode is activated...")
            if maxy<20:
                if verbose: log.write(" - maxy <20 , no need to cut.")
                cut=0
            else:
                cut = 20
                cutfactor = ( maxy - cut )/5
        if cut:
            if verbose: log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")

            maxticker=int(max(np.round(sumstats["scaled_P_1"].max(skipna=True)) , np.round(sumstats["scaled_P_2"].max(skipna=True))))
            
            sumstats.loc[sumstats["scaled_P_1"]>cut,"scaled_P_1"] = (sumstats.loc[sumstats["scaled_P_1"]>cut,"scaled_P_1"]-cut)/cutfactor +  cut
            sumstats.loc[sumstats["scaled_P_2"]>cut,"scaled_P_2"] = (sumstats.loc[sumstats["scaled_P_2"]>cut,"scaled_P_2"]-cut)/cutfactor +  cut
            maxy = (maxticker-cut)/cutfactor + cut
    ##########################################################################################################################
    plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P_1',
        hue='chr_hue',
        palette= palette,
        legend=legend,
        style=style,
        size="s1",
        sizes=marker_size,
        linewidth=linewidth,
        zorder=2,ax=ax1,edgecolor="black")     
    
    plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P_2',
        hue='chr_hue',
        palette= palette,
        legend=legend,
        style=style,
        size="s2",
        sizes=marker_size,
        linewidth=linewidth,
        zorder=2,ax=ax5,edgecolor="black") 
    
    ## X #######################################################################################################################
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    
    ax1.set_xlabel("")
    ax1.set_xticks(chrom_df)
    ax1.set_xticklabels(chrom_df.index.map(get_number_to_chr()),fontsize=fontsize,family="sans-serif") 
    
    ax5.set_ylim(ax1.get_ylim())
    ax5.set_xlim(ax1.get_xlim())
    ax5.set_xlabel("")
    ax5.set_xticks(chrom_df)
    ax5.set_xticklabels([])
    ax5.xaxis.set_ticks_position("top")
    
    ## Y #######################################################################################################################
    sigline = ax1.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
    sigline = ax5.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
     
    for ax in [ax1,ax5]:
        if cut == 0: 
            ax.set_ylim(skip,maxy*1.2)
        if cut:
            cutline = ax.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker],fontsize=fontsize,family="sans-serif")
            else:
                ax.set_yticks([x for x in range(skip,cut+1,2)])
                ax.set_yticklabels([x for x in range(skip,cut+1,2)],fontsize=fontsize,family="sans-serif")
            ax.set_ylim(bottom = skip)
            
    ax5.invert_yaxis() 
    ### Spines ####################################################################################################################
    
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(True)
    ax1.spines["bottom"].set_visible(True)
    ax5.spines["top"].set_visible(True)
    ax5.spines["right"].set_visible(False)
    ax5.spines["left"].set_visible(True)
    ax5.spines["bottom"].set_visible(False)
    ##########################################################################################################################
    
    #return fig
    
    
    
        