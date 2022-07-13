import sys
import gwaslab as gl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib
from statsmodels.stats.multitest import fdrcorrection
#################################################################################################
def convert_p_to_width(p,sig_level):
    width_factor= -np.log10(sig_level)
    # if significant, full square
    if p<sig_level: 
        return 1
    else:
        #scaled using mlog10(p)
        return max(-np.log10(p)/width_factor,0.1)

def conver_rg_to_color(rg,cmap):
    #(1,1)
    if rg>1: rg=1
    if rg<-1: rg=-1
    return cmap((rg+1)/2)
####################################################################################################

def plot_rg(ldscrg,
        p1="p1",p2="p2",rg="rg",p="p",
        sig_level=0.05,
        rganno=False,
        correction="",
        cmap = matplotlib.cm.get_cmap('RdBu'),
        log=gl.Log(),
        verbose=True,
        sort_key=None,
        square=False,
        colorbarargs={"shrink":0.82},
        **args):
    
    if verbose: log.write("Total non-NA records:",len(ldscrg.dropna(subset=[p])))
    df=ldscrg.dropna(subset=[p]).copy()
    df["p1p2"]=df.apply(lambda x:"_".join(sorted([x[p1],x[p2]])),axis=1)
    
    dfp=ldscrg.loc[ldscrg[p1]!=ldscrg[p2],:].dropna(subset=[p]).copy()
    dfp["p1p2"]=dfp.apply(lambda x:"_".join(sorted([x[p1],x[p2]])),axis=1)
    dfp = dfp.drop_duplicates(subset=["p1p2"]).copy()
    
    if verbose: log.write("Valid unique records:",len(dfp))
    if verbose: log.write("Significant correlations after Bonferroni correction:",sum(dfp[p]<0.05/len(dfp)))
    
    if correction=="fdr":
        dfp["fdr_p"]=fdrcorrection(dfp[p],alpha=1)[1]
        print(dfp["fdr_p"])
        dfp["fdr"]=fdrcorrection(dfp[p],alpha=sig_level)[0]
        if verbose: log.write("Significant correlations after FDR correction:",sum(dfp["fdr"]))
        dfp=dfp.set_index("p1p2").loc[:,"fdr_p"].to_dict()
        print(dfp)
    else:
        dfp=dfp.set_index("p1p2").loc[:,p].to_dict()
        print(dfp)
    #########ticks dict###########################################   
    dic_p1={}
    dic_p2={}
    dic_p1_r={}
    dic_p2_r={}
    
    if sort_key is None:
        for i,p1_name in enumerate(df[p1].sort_values(ascending=False).drop_duplicates()):
            dic_p1[p1_name]  = i
            dic_p1_r[i] = p1_name
        for i,p2_name in enumerate(df[p2].sort_values().drop_duplicates()):
            dic_p2[p2_name]  = i
            dic_p2_r[i] = p2_name
    else:
        for i,p1_name in enumerate(df[p1].sort_values(ascending=False,key=sort_key).drop_duplicates()):
            dic_p1[p1_name]  = i
            dic_p1_r[i] = p1_name
        for i,p2_name in enumerate(df[p2].sort_values(key=sort_key).drop_duplicates()):
            dic_p2[p2_name]  = i
            dic_p2_r[i] = p2_name

    df["y"]=df[p1].map(dic_p1)
    df["y_x"]=df[p1].map(dic_p2)
    df["x"]=df[p2].map(dic_p2)
    df["x_y"]=df[p2].map(dic_p1)
    
    if verbose: log.write("Plotting...")
    ########ticks###############################################
    fig,ax = plt.subplots(dpi=300,**args)
    xticks=df["x"].sort_values().drop_duplicates().astype(int)
    yticks=df["y"].sort_values().drop_duplicates().astype(int)  
    ax.xaxis.tick_top()
    ax.set_xticks(xticks,minor=False)
    ax.set_xticks(xticks+0.5,minor=True)
    ax.set_yticks(yticks,minor=False)
    ax.set_yticks(yticks+0.5,minor=True)
    ax.grid(visible=True,which="minor")
    ax.set_xlim( -0.5 ,df["x"].max()+0.5)
    ax.set_ylim( -0.5 ,df["y"].max()+0.5)
    ax.tick_params('both', length=5, width=2, which='major')
    ax.tick_params('both', length=0, width=0, which='minor')
    
    #labels
    ax.set_yticklabels(yticks.map(dic_p1_r),fontsize=15)
    ax.set_xticklabels(xticks.map(dic_p2_r),rotation=45,horizontalalignment="left", verticalalignment="bottom",fontsize=15)
    
    width_max=1

    
    #########patches###########################################
    
    squares=[]
    panno=[]
    rgtoanno=[]
    maxsigp=sig_level
    
    #if correction=="fdr":
    #    if len(df.loc[df["fdr"]==True,p])>=1:
    #        maxsigp = df.loc[df["fdr"]==True,p].max()*1.0001
    #        
    #    else:
    #        maxsigp = sig_level/len(df.dropna(subset=[p]))
    if correction=="fdr":
        p="fdr_p"
    
        
    for i,row in df.iterrows():
        xcenter=row["x"]
        ycenter=row["y"]
        if row[p1]==row[p2]:
            width=1
            x=xcenter-width/2
            y=ycenter-width/2
            rgba = conver_rg_to_color(1,cmap)
            
        else:    
            adjusted_p = dfp["_".join(sorted([row[p1],row[p2]]))]
            if adjusted_p<0.05 and square is True:
                if xcenter + ycenter < len(df[p1].unique()):
                    panno.append([xcenter,ycenter,adjusted_p])
            elif adjusted_p<0.05:
                panno.append([xcenter,ycenter,adjusted_p])
            
            width= convert_p_to_width(adjusted_p,sig_level)    
            x=xcenter-width/2
            y=ycenter-width/2           
            rgba = conver_rg_to_color(row[rg],cmap)
            
            if xcenter + ycenter > len(df[p1].unique())-1 and (square is True) and (rganno is True):
                rgtoanno.append([xcenter,ycenter,row[rg],rgba])  
        
        if xcenter + ycenter < len(df[p1].unique()) and (square is True) and (rganno is True):
            squares.append(patches.Rectangle((x,y),width=width,height=width,fc=rgba,ec="white",lw=0))  
        elif (square is not True):
            squares.append(patches.Rectangle((x,y),width=width,height=width,fc=rgba,ec="white",lw=0))  
    
    
    
    
    squares_collection = matplotlib.collections.PatchCollection(squares,match_original=True)
    ax.add_collection(squares_collection)       
    
    if rganno is not False:
        for i in rgtoanno:
            if i[2]>1: i[2]=1
            if i[2]<-1: i[2]=-1
            ax.text(i[0],i[1],"{:.3f}".format(i[2]),color=i[3],weight="bold",ha="center", va="center",font="Arial")
    
    for i in panno:
        if i[2]<sig_level/len(dfp):
            ax.text(i[0],i[1],"**",color="white",weight="bold",ha="center", va="center",font="Arial")
        else:
            ax.text(i[0],i[1],"*",color="white",weight="bold",ha="center", va="center",font="Arial")

            
    ## color bar ###############################################
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, **colorbarargs)

    return fig,ax,log
    