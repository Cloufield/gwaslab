import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
from gwaslab import getsig

def mplot(insumstats,chrom,pos,p,
          scaled=False,
          cut=0,
          cutfactor=10,
          windowsizekb=500,
          anno=None,
          sig_level=5e-8,
          title=None,
          figsize =(5,15),
          colors = ["#000042", "#7878BA"],
          suggestive_sig_level=5e-6,
          verbose=True):
    
# printing meta info ###############################
    if verbose: print("Basic settings:")
    if verbose: print("  - Genome-wide significance level is set to "+str(sig_level)+" ...")
    if verbose: print("  - Raw input contains "+str(len(insumstats))+" variants...")
    if verbose: print("Strat conversion and QC:")
    
# read sumstats ###############################
    sumstats=pd.DataFrame()
    sumstats["P-value_association"]=insumstats[p].astype(np.float64)
    sumstats["CHROM"]=insumstats[chrom].astype(int)
    sumstats["POS"]=insumstats[pos].astype(int)
    if anno and anno!=True:
            sumstats["Annotation"]=insumstats[anno].astype(object)
            
# P value conversion ###############################  
    if scaled:
        if verbose:print("  - P values are already -log10 scaled!")
        sumstats["scaled_P"] = sumstats["P-value_association"]
    else:
        if verbose:print("  - P values are being -log10 scaled...")
        sumstats["scaled_P"] = -np.log10(sumstats["P-value_association"])

    if verbose: print("  - "+str(len(sumstats[sumstats["scaled_P"].isin([np.nan, np.inf, -np.inf, float('inf'), None])])) + " na/inf/-inf variants will be removed..." )
    sumstats = sumstats[~sumstats["scaled_P"].isin([np.nan, np.inf, -np.inf, float('inf'), None])]
# shrink at a certain value #############################################
    if cut:
        if verbose: print("  - Minus log10(P) values will be shrunk at " + str(cut)+" with a shrinkage factor of " + str(cutfactor)+"...")
        maxticker=int(np.round(sumstats["scaled_P"].max()))
        sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
        maxy = (maxticker-cut)/cutfactor + cut

# if metal format#############################################
    #sumstats["CHROM"]=sumstats["MarkerName"].str.split(":").apply(lambda x: x[0]).astype(int)
    #sumstats["POS"]=sumstats["MarkerName"].str.split(":").apply(lambda x: x[1]).astype(int)
# create index ##########################

    if verbose:print("Plotting "+str(len(sumstats))+" variants...")
    #sort & add id
    sumstats = sumstats.sort_values(["CHROM","POS"])
    sumstats["id"]=range(len(sumstats))
    sumstats=sumstats.set_index("id")

    #create a position dictionary
    posdic=sumstats.groupby("CHROM")["POS"].max()
    posdiccul=dict(posdic)
    posdiccul[0]=0
    for i in range(2,sumstats["CHROM"].nunique()+1):
        posdiccul[i]=posdiccul[i-1]+posdiccul[i]

    #convert base pair postion to x axis position
    sumstats["add"]=sumstats["CHROM"].apply(lambda x : posdiccul[int(x)-1])
    sumstats["i"]=sumstats["POS"]+sumstats["add"]

    #for plot
    chrom_df=sumstats.groupby('CHROM')['i'].median()
    chrom_df=chrom_df+((chrom_df.index.astype(int))-1)*len(sumstats)*0.02

    sumstats["i"]=sumstats["i"]+((sumstats["CHROM"].astype(int))-1)*len(sumstats)*0.01

# assign marker size ##############################################
    sumstats["s"]=1
    sumstats.loc[sumstats["scaled_P"]>-np.log10(5e-4),"s"]=2
    sumstats.loc[sumstats["scaled_P"]>-np.log10(suggestive_sig_level),"s"]=3
    sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level),"s"]=10

# Plotting ########################################################  
    plt.figure(figsize=figsize,constrained_layout=True)
    plot = sns.relplot(data=sumstats, x='i', y='scaled_P', aspect=2.3,
                       hue='CHROM',
                       palette=sns.color_palette(colors,n_colors=sumstats["CHROM"].nunique()),
                       legend=None,
                       size="s",
                       sizes=(10,40),
                       linewidth=0,
                       zorder=2)   
        
    plot.ax.set_xlabel('CHROM'); plot.ax.set_xticks(chrom_df);
    plot.ax.set_xticklabels(chrom_df.index)
    
    sigline = plot.ax.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color="grey",zorder=1)
    if cut:
        cutline=plot.ax.axhline(y=cut, linewidth = 2,linestyle="--",color="#ebebeb",zorder=1)
        plot.ax.set_yticks([x for x in range(0,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
        plot.ax.set_yticklabels([x for x in range(0,cut+1,2)]+[maxticker])

# Annotation #######################################################
    if anno and anno!=True:
        to_annotate=getsig(sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level)],"Annotation",'CHROM',"POS","P-value_association",sig_level=sig_level,windowsizekb=windowsizekb,verbose=False)
    else:
        to_annotate=getsig(sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level)],"i",'CHROM',"POS","P-value_association",windowsizekb=windowsizekb,verbose=False,sig_level=sig_level)
    if verbose:print("  - Found "+str(len(to_annotate))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")

# Annotation #######################################################
    if anno and len(to_annotate)>0:
        text = []
        last_pos=0
        
        for rowi,row in to_annotate.iterrows():
            if row["i"]>last_pos+0.03*sumstats["i"].max():
                last_pos=row["i"]
            else:
                last_pos+=0.03*sumstats["i"].max()
                
            armB_length_in_point=plot.ax.transData.transform((0, 0.95*maxy))[1]-plot.ax.transData.transform((0, row["scaled_P"]+1))[1]
            armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else plot.ax.transData.transform((0, maxy+2))[1]-plot.ax.transData.transform((0,  row["scaled_P"]+1))[1] 
            
            if anno==True:
                annotation_text=str(int(row["CHROM"])) +":"+ str(int(row["POS"]))
                annotation_col="CHR:POS"
            elif anno:
                annotation_text=row["Annotation"]
                annotation_col=anno
                
            text.append(plot.ax.annotate(annotation_text,
                             xy=(row["i"],row["scaled_P"]+1),
                             xytext=(last_pos,1.15*maxy),rotation=40,
                                        ha="left",va="bottom",
                             arrowprops=dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                             connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                                        )
            
                       )
        #plt.subplots_adjust(top=0.7, wspace=None, hspace=None)
        
        if verbose: print("  - Annotating using column "+annotation_col+"...")
    else:
        if verbose: print("  - Skip annotating")

    plot.ax.set_ylabel("-log10(P)")
    plot.ax.set_xlabel("Chromosomes")
    plot.ax.spines["top"].set_visible(False)
    plot.ax.spines["right"].set_visible(False)
    plot.ax.spines["left"].set_visible(True)
    
    if verbose: print("  - Created Manhattan plot successfully!")
    if title and anno and len(to_annotate)>0: 
        pad=(plot.ax.transData.transform((0, 1.35*maxy))[1]-plot.ax.transData.transform((0, maxy)))[1]
        plt.title(title,pad=pad)
    elif title:
        plt.title(title)
    return plot