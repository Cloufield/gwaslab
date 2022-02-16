import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp

def mplot(insumstats,chrom,pos,p,title="Manhatann-plot",figsize=(8,8),scaled=False,cut=0,cutfactor=10):
    colors = ["#000042", "#7878BA"]
    print("Processing raw input of "+str(len(insumstats))+" variants...")
    #read sumstats
    sumstats=pd.DataFrame()
    sumstats["P-value_association"]=insumstats[p].astype(np.float64)
    sumstats["CHROM"]=insumstats[chrom].astype(int)
    sumstats["POS"]=insumstats[pos].astype(int)
    
    if scaled:
        print("P values are already -log10 scaled!")
        sumstats["scaled_P"] = sumstats["P-value_association"]
    else:
        print("P values are being -log10 scaled...")
        sumstats["scaled_P"] = -np.log10(sumstats["P-value_association"])
    
    print(str(len(sumstats[sumstats["scaled_P"].isin([np.nan, np.inf, -np.inf, float('inf'), None])])) + " na/inf/-inf variants will be removed..." )
    sumstats = sumstats[~sumstats["scaled_P"].isin([np.nan, np.inf, -np.inf, float('inf'), None])]
    
    if cut:  
        print("-log10(P) values will be cut at " + str(cut)+" with cut factor of " + str(cutfactor)+"...")
        maxticker=int(np.round(sumstats["scaled_P"].max()))
        sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
    
    
    #sumstats["CHROM"]=sumstats["MarkerName"].str.split(":").apply(lambda x: x[0]).astype(int)
    #sumstats["POS"]=sumstats["MarkerName"].str.split(":").apply(lambda x: x[1]).astype(int)
    print("Plotting "+str(len(sumstats))+" variants...")
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
    

    
    plt.figure(figsize=figsize)
    plot = sns.relplot(data=sumstats, x='i', y='scaled_P', aspect=2, 
                       hue='CHROM',
                       palette=sns.color_palette(colors,n_colors=sumstats["CHROM"].nunique()), 
                       legend=None,
                       s=10,
                       linewidth=0) 
    
    #chrom_df=sumstats.groupby('CHROM')['i'].median()
    plot.ax.set_xlabel('CHROM'); plot.ax.set_xticks(chrom_df);
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.axhline(y=-np.log10(5e-8), linewidth = 2,linestyle="--",color="grey")
    if cut:
        plot.ax.axhline(y=cut, linewidth = 2,linestyle="--",color="grey")
        plot.ax.set_yticks([x for x in range(0,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
        plot.ax.set_yticklabels([x for x in range(0,cut+1,2)]+[maxticker])
    plot.ax.set_ylabel("-log10(P)")
    plot.ax.set_xlabel("Chromosomes")
    plot.fig.suptitle(title)
    return plot