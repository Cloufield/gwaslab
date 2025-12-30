from typing import TYPE_CHECKING, Optional, List, Dict, Any, Union, Tuple, Callable
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib
from gwaslab.info.g_Log import Log
import scipy.stats as ss
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes

#################################################################################################
def convert_p_to_width(p: float, sig_level: float) -> float:
    width_factor= -np.log10(sig_level)
    # if significant, full square
    if p<sig_level: 
        return 1
    else:
        #scaled using mlog10(p)
        return max(-np.log10(p)/width_factor,0.1)

def convert_rg_to_color(rg: float, cmap: Any) -> Any:
    #(1,1)
    if rg>1: rg=1
    if rg<-1: rg=-1
    return cmap((rg+1)/2)
####################################################################################################

def plot_rg(ldscrg,
        p1="p1",
        p2="p2",
        rg="rg",
        p="p",
        sig_levels=None,
        rganno="non",
        panno=True,
        corrections=None,
        panno_texts=None,
        equal_aspect=True,
        cmap = None,
        full_cell =None,
        log=Log(),
        panno_kwargs=None,
        rganno_kwargs=None,
        verbose=True,
        asize=10,
        sort_key=None,
        square=False,
        colorbar_kwargs=None,
        fig_kwargs=None,
        xticklabel_kwargs=None,
        yticklabel_kwargs=None,
        fdr_method="bh",
        fontsize=10,
        save=None,
        save_kwargs=None):
    
    log.write("Start to create ldsc genetic correlation heatmap..." ,verbose=verbose)
    # configure arguments
    if cmap is None:
        try: #matplotlib >=3.7
            cmap = matplotlib.colormaps.get_cmap('RdBu')
        except:
            # Fallback for older matplotlib versions
            cmap = matplotlib.cm.get_cmap('RdBu')
    style = set_plot_style(
        plot="plot_rg",
        fig_kwargs=fig_kwargs or {"dpi":300},
        save_kwargs=save_kwargs,
        fontsize=fontsize,
        fontfamily="Arial",
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    fontsize = style["fontsize"]
    if colorbar_kwargs is None:
        colorbar_kwargs={"shrink":0.82}
    if yticklabel_kwargs is None:
        yticklabel_kwargs={"fontsize":fontsize, "fontfamily":style["font_family"]}
    if xticklabel_kwargs is None:    
        xticklabel_kwargs={"rotation":45,"horizontalalignment":"left", "verticalalignment":"bottom","fontsize":fontsize, "fontfamily":style["font_family"]}
    if sig_levels is None:
        sig_levels = [0.05]
    if corrections is None:
        corrections = ["non", "fdr","bon"]
    if panno_texts is None:
        panno_texts = ["*"*(i+1) for i in range(len(sig_levels)*len(corrections))]
    if full_cell is None:
        full_cell = ("fdr",0.05)
    if rganno_kwargs is None:
        rganno_kwargs ={}
    if save_kwargs is None:
        save_kwargs = {}
    
    #drop na records in P column 
    log.write("Raw dataset records:",len(ldscrg) ,verbose=verbose)
    df=ldscrg.dropna(subset=[p]).copy()
    
    log.write(" -Raw dataset non-NA records:",len(df) ,verbose=verbose)
    # create unique pair column
    df["p1p2"]=df.apply(lambda x:"_".join(sorted([x[p1],x[p2]])),axis=1)
    
    log.write("Filling diagonal line and duplicated pair for plotting..." ,verbose=verbose)
    # fill na
    df_fill_reverse = df.loc[(df[p2].isin(df[p1].values)) & (df[p1].isin(df[p2].values)),:].copy()
    df_fill_reverse = df_fill_reverse.rename(columns={p1:p2,p2:p1})
    
    # fill dia
    df_fill_dia = pd.DataFrame(columns=df.columns)
    p1_dup_list = list(df.loc[(df[p2].isin(df[p1].values)),"p2"].values)
    p2_dup_list = list(df.loc[(df[p1].isin(df[p2].values)),"p1"].values)
    p_dup_list = p2_dup_list + p1_dup_list
    if len(set(p_dup_list)) > 0:
        log.write(" -Diagonal records:", len(set(p_dup_list)) ,verbose=verbose)
    df_fill_dia["p1"] = p_dup_list
    df_fill_dia["p2"] = df_fill_dia["p1"] 
    df_fill_dia["rg"] = 1

    df_fill_na = pd.DataFrame(columns=df.columns)
    df_fill_na[[p1,p2]] = [(i,j) for i in df[p1].sort_values(ascending=False).drop_duplicates() for j in df[p2].sort_values(ascending=False).drop_duplicates()]
    
    to_concate=[]
    for i in [df,df_fill_reverse,df_fill_dia,df_fill_na]:
        if len(i)>0:
            to_concate.append(i.dropna(axis=1))
    
    # fill diagonal
    df = pd.concat(to_concate,ignore_index=True).sort_values(by=p).drop_duplicates(subset=[p1,p2])
    
    #log.write(" -Dataset shape match:", len(df)==)
    #
    ## remove record with p1 = p2, dropna in P column
    dfp=ldscrg.loc[ldscrg[p1]!=ldscrg[p2],:].dropna(subset=[p]).copy()
    
    ## create pair column
    dfp["p1p2"]=dfp.apply(lambda x:"_".join(sorted([x[p1],x[p2]])),axis=1)
    
    ## drop duplicate and keep only unique pairs 
    dfp = dfp.drop_duplicates(subset=["p1p2"]).copy()
    
    log.write("Valid unique trait pairs:",len(dfp) ,verbose=verbose)
    log.write(" -Valid unique trait1:",dfp["p1"].nunique() ,verbose=verbose)
    log.write(" -Valid unique trait2:",dfp["p2"].nunique() ,verbose=verbose)
    log.write(" -Significant correlations with P < 0.05:",sum(dfp[p]<0.05) ,verbose=verbose)
    log.write(" -Significant correlations after Bonferroni correction:",sum(dfp[p]<(0.05/len(dfp))) ,verbose=verbose)
    
    #if correction=="fdr":
        # fdr corrected p
    #dfp["fdr_p"]=fdrcorrection(dfp[p],alpha=1,method=fdr_method)[1]
        # is fdr < sig_level
    #dfp["fdr"]=fdrcorrection(dfp[p],alpha=0.05,method=fdr_method)[0]
    
    dfp["fdr_p"]=ss.false_discovery_control(dfp[p],method=fdr_method)
    dfp["fdr"]  =ss.false_discovery_control(dfp[p],method=fdr_method) < 0.05

    log.write(" -Significant correlations with FDR <0.05:",sum(dfp["fdr"]) ,verbose=verbose)
        # convert to dict for annotation and plotting
    df_rawp = dfp.set_index("p1p2").loc[:,p].to_dict()
    dfp = dfp.set_index("p1p2").loc[:,"fdr_p"].to_dict()

    #########ticks dict###########################################   
    dic_p1={}
    dic_p2={}
    
    dic_p1_r={}
    dic_p2_r={}
    
    ## sort position 
    if sort_key is None:
        # alphabetic order
        for i,p1_name in enumerate(df[p1].sort_values(ascending=False).drop_duplicates()):
            dic_p1[p1_name]  = i
            dic_p1_r[i] = p1_name
        for i,p2_name in enumerate(df[p2].sort_values().drop_duplicates()):
            dic_p2[p2_name]  = i
            dic_p2_r[i] = p2_name
    else:
        # user-provided order
        for i,p1_name in enumerate(df[p1].sort_values(ascending=False,key=sort_key).drop_duplicates()):
            dic_p1[p1_name]  = i
            dic_p1_r[i] = p1_name
        for i,p2_name in enumerate(df[p2].sort_values(key=sort_key).drop_duplicates()):
            dic_p2[p2_name]  = i
            dic_p2_r[i] = p2_name

    # assign coordinate
    df["y"]=df[p1].map(dic_p1)
    df["y_x"]=df[p1].map(dic_p2)
    df["x"]=df[p2].map(dic_p2)
    df["x_y"]=df[p2].map(dic_p1)
    
    log.write("Plotting heatmap..." ,verbose=verbose)
    ########ticks###############################################
    fig,ax = plt.subplots(**fig_kwargs)
    
    # configure x/y ticks
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
    ax.set_yticklabels(yticks.map(dic_p1_r),**yticklabel_kwargs)

    ax.set_xticklabels(xticks.map(dic_p2_r),**xticklabel_kwargs)
    
    #########patches###########################################
    
    squares=[]
    panno_list={1:{},2:{}}
    rgtoanno=[]
    
    log.write("Full cell : {}-corrected P == {}".format(full_cell[0],full_cell[1]) ,verbose=verbose)

    for i,row in df.iterrows():
        xcenter=row["x"]
        ycenter=row["y"]

        if np.isnan(row[rg]):
            width=1
            x=xcenter-width/2
            y=ycenter-width/2
            ax.plot([x,x+width],[y,y+width],c="grey")
            ax.plot([x,x+width],[y+width,y],c="grey")
        
        else:
            if row[p1]==row[p2]:
                # diagonal line
                width=1
                x=xcenter-width/2
                y=ycenter-width/2
                rgba = convert_rg_to_color(1,cmap)    
            else:    
                # get the adjusted p value from dict
                if  xcenter + ycenter < len(df[p1].unique()): 
                    panno_set=1
                else:
                    panno_set=2
                for i,correction in enumerate(corrections):
                    for j,sig_level in enumerate(sig_levels):
                        
                        index = len(sig_levels)*i + j
                        
                        p1p2="_".join(sorted([row[p1],row[p2]]))
                        
                        raw_p = df_rawp[p1p2]
                    
                        if correction in ["B","bonferroni ","bon","Bon","b"]:
                            current_threhold = sig_level/len(dfp)
                            if raw_p < current_threhold:
                                panno_list[panno_set][p1p2] = [xcenter,ycenter,raw_p,"bon",panno_texts[index]]

                        elif correction in ["fdr","FDR","F","f"]:
                            adjusted_p = dfp[p1p2]
                            if adjusted_p < sig_level and square is True:
                                #if square is True, only annotate half 
                                if xcenter + ycenter < len(df[p1].unique()):
                                    panno_list[panno_set][p1p2]=[xcenter,ycenter,adjusted_p,"fdr",panno_texts[index]]
                            elif adjusted_p < sig_level:
                                    panno_list[panno_set][p1p2]=[xcenter,ycenter,adjusted_p,"fdr",panno_texts[index]]

                        elif correction == "non":
                            if raw_p < sig_level:
                                panno_list[panno_set][p1p2]=[xcenter,ycenter,"raw",raw_p,panno_texts[index]]
                
                # configuring the square 
                if full_cell[0] == "fdr":
                    width= convert_p_to_width(adjusted_p,full_cell[1])
                elif full_cell[0] == "bon":
                    width= convert_p_to_width(raw_p*len(dfp),full_cell[1])
                else:
                    width= convert_p_to_width(raw_p,full_cell[1])

                x=xcenter-width/2
                y=ycenter-width/2           
                rgba = convert_rg_to_color(row[rg],cmap)
                if xcenter + ycenter > len(df[p1].unique())-1 and (square is True) and (rganno == "half"):
                    rgtoanno.append([xcenter,ycenter,row[rg],rgba])  
                elif "full" in rganno:
                    rgtoanno.append([xcenter,ycenter,row[rg],rgba])  
            
            #if xcenter + ycenter < len(df[p1].unique()) and (square is True) and (rganno == "half"):
            #    squares.append(patches.Rectangle((x,y),width=width,height=width,fc=rgba,ec="white",lw=0))  
            #elif (square is not True):
            if ("nb" not in rganno):
                if rganno == "half":
                    if xcenter + ycenter < len(df[p1].unique()) and (square is True):
                        squares.append(patches.Rectangle((x,y),width=width,height=width,fc=rgba,ec="white",lw=0))  
                else:
                    squares.append(patches.Rectangle((x,y),width=width,height=width,fc=rgba,ec="white",lw=0))
    
    squares_collection = matplotlib.collections.PatchCollection(squares,match_original=True)
    ax.add_collection(squares_collection)       
    
    if rganno is not False:
        rganno_default_kwargs = {"weight":"bold","ha":"center", "va":"center", "fontfamily":"Arial","fontsize":fontsize}
        for key, value in rganno_kwargs.items():
            rganno_default_kwargs[key] = value
        for i in rgtoanno:
            if i[2]>1: i[2]=1
            if i[2]<-1: i[2]=-1
            if "color" in rganno_default_kwargs.keys() or "c" in rganno_default_kwargs.keys():
                ax.text(i[0],i[1],"{:.3f}".format(i[2]),**rganno_default_kwargs)
            else:
                ax.text(i[0],i[1],"{:.3f}".format(i[2]),color=i[3],**rganno_default_kwargs)
    
    # configure args for p annotation
    panno_default_kwargs={"size":asize,"color":"white","weight":"bold","horizontalalignment":"center","verticalalignment":"center_baseline","font":"Arial"}
    if panno_kwargs is not None:
        for key, value in panno_kwargs.items():
            panno_default_kwargs[key] = value

    # annotate p
    if panno is True:
        log.write("P value annotation text (Order: Bon -> FDR -> Pnom): " ,verbose=verbose)
        for i,correction in enumerate(corrections):
            for j,sig_level in enumerate(sig_levels):
                index = len(sig_levels)*i + j
                log.write(" -{} : {}-corrected P < {} ".format(panno_texts[index], correction, sig_level) ,verbose=verbose)
        for panno_set_number in panno_list.keys():
            for key, i in panno_list[panno_set_number].items():
                if panno_set_number == 1:
                    ax.text(i[0],i[1]-0.1,i[4], **panno_default_kwargs)
                else:
                    ax.text(i[0],i[1]-0.1,i[4], **panno_default_kwargs)
            
    ## color bar ###############################################
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, **colorbar_kwargs)

    if equal_aspect is True:
        ax.set_aspect('equal', adjustable='box')
    
    save_figure(fig, save, keyword="ldscrg",save_kwargs=save_kwargs, log=log, verbose=verbose)

    log.write("Finished creating ldsc genetic correlation heatmap!" ,verbose=verbose)

    return fig,ax,log,df
    
