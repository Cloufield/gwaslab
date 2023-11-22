import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy import stats, optimize
from statsmodels.stats.meta_analysis import combine_effects
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
# plot_forest
# plot row

def plot_forest(data, 
                study_col, 
                group_col=False,
                compact_factor=1,
                width_ratios = [2,6,2],
                sharex="col",
                meta=False,
                combine_effects_args={}):
    if data is pd.DataFrame():
        meta_df_all = data
    else:
        meta_df_all = pd.read_csv(data,sep="\s+")
    if group_col == False:
        group_col = "_GROUP"
        meta_df_all[group_col] = "Group1"

    df_array = []
    df_array_name = []
    df_array_het=[]
    
    # rename 


    for i in meta_df_all[group_col].unique():
        meta_df = meta_df_all.loc[meta_df_all[group_col]==i,:]
        combined_results = combine_effects(meta_df["beta"], 
                                           meta_df["se"]**2, 
                                           row_names=meta_df[study_col],
                                           **combine_effects_args)
        if meta==True:
            Q=combined_results.test_homogeneity().statistic
            Qdf=combined_results.test_homogeneity().df
            Qp=combined_results.test_homogeneity().pvalue
            I2 = (Q-Qdf)/Q *100
       
        df_to_plot = combined_results.summary_frame().iloc[:-3,:]
        df_to_plot["YORDER"]=range(len(df_to_plot),0,-1)
        df_to_plot["OR"]= np.exp(df_to_plot["eff"])
        df_to_plot["OR_95L"]= np.exp(df_to_plot["ci_low"])
        df_to_plot["OR_95U"]= np.exp(df_to_plot["ci_upp"])
        df_to_plot["OR_95L-OR"] = df_to_plot["OR"] - df_to_plot["OR_95L"]
        df_to_plot["OR_95U-OR"] = df_to_plot["OR_95U"] - df_to_plot["OR"]
        df_to_plot["MARKERSIZE"] = 1 / (df_to_plot["sd_eff"]**2)
        df_array.append(df_to_plot)
        df_array_name.append(i)
        if meta==True:
            df_array_het.append((Q,Qdf,Qp,I2))
        else:
            df_array_het.append(list())

    nrows=len(df_array)
    fig, axes = plt.subplots(nrows=nrows,ncols=3, sharex=sharex,gridspec_kw={'width_ratios': width_ratios},figsize=(15,len(df_to_plot)*nrows/compact_factor),dpi=300)
    if nrows==1:
        plot_row(axes,
                 df_to_plot = df_array[0], 
                 title=0,
                 group_name=df_array_name[0],
                 het_q_qdf_qp_i=df_array_het[0],
                 meta=meta)
    else:
        for i in range(nrows):
            plot_row(axes[i],
                     df_to_plot = df_array[i], 
                     title=i, 
                     group_name=df_array_name[i],
                     het_q_qdf_qp_i=df_array_het[i],
                     meta=meta)
    fig.tight_layout()




def plot_row(axes,group_name,df_to_plot,meta,het_q_qdf_qp_i=None,title=0):
    studies = df_to_plot.loc[df_to_plot.index!="fixed effect",:].copy()
    studies["MARKERSIZE_N"] = studies["MARKERSIZE"]/np.nanmax(studies["MARKERSIZE"]) * 0.8
    pixels = axes[1].transData.transform((0,len(df_to_plot)+0.5))[1] - axes[1].transData.transform((0,0))[1]
    studies["MARKERSIZE_P"] =  studies["MARKERSIZE_N"] * (pixels/(len(df_to_plot)+0.5)) 
    axes[1].scatter(x=studies["OR"],y=studies["YORDER"],s=studies["MARKERSIZE_P"]*2 ,marker="s",color="grey")
    axes[1].errorbar(x=studies["OR"],y=studies["YORDER"],xerr=studies[["OR_95L-OR","OR_95U-OR"]].T ,elinewidth=2,lw=0,marker="s",color="black")

    # diamond patch coordinates
    xm=df_to_plot.loc[df_to_plot.index=="fixed effect","OR"]
    xl=df_to_plot.loc[df_to_plot.index=="fixed effect","OR_95L"]
    xr=df_to_plot.loc[df_to_plot.index=="fixed effect","OR_95U"]
    ym=df_to_plot.loc[df_to_plot.index=="fixed effect","YORDER"]
    yu=df_to_plot.loc[df_to_plot.index=="fixed effect","YORDER"]+0.25
    yd=df_to_plot.loc[df_to_plot.index=="fixed effect","YORDER"]-0.25

    xy1=(xl.values[0],ym.values[0])
    xy2=(xm.values[0],yu.values[0])
    xy3=(xr.values[0],ym.values[0])
    xy4=(xm.values[0],yd.values[0])

    xy=[xy1,xy2,xy3,xy4]
    patches = []
    polygon = Polygon(xy, closed=True)
    patches.append(polygon)
    p = PatchCollection(patches,color="grey",edgecolor="black",alpha=0.8)
    axes[1].add_collection(p)

    ## add lines
    axes[1].axvline(x=1,color="black")
    axes[1].axvline(x=df_to_plot.loc[df_to_plot.index=="fixed effect","OR"].values[0],color="grey",ls="--")
    axes[1].set_ylim(0.25,len(df_to_plot)+0.75)
    axes[2].set_ylim(0.25,len(df_to_plot)+0.75)
    axes[2].set_xlim([0,1])
    ## set spine invisible
    axes[1].spines["left"].set(visible=False)
    axes[1].spines["right"].set(visible=False)

    axes[2].spines["top"].set(visible=False)
    axes[2].spines["bottom"].set(visible=False)
    axes[2].spines["left"].set(visible=False)
    axes[2].spines["right"].set(visible=False)

    axes[0].spines["top"].set(visible=False)
    axes[0].spines["bottom"].set(visible=False)
    axes[0].spines["left"].set(visible=False)
    axes[0].spines["right"].set(visible=False)

    ## set tickers invisible 
    axes[0].set_xticklabels([])
    axes[0].set_yticklabels([])
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[2].set_xticklabels([])
    axes[2].set_yticklabels([])
    axes[2].set_xticks([])
    axes[2].set_yticks([])
    

    # title
    if title ==0:
        #axes[0].set_title("Group",size=15)
        axes[0].set_title("Groups",size=15,loc="left",fontweight ="bold")
        axes[1].set_title("Odds Ratio",size=15,fontweight ="bold")
        axes[1].set_title("Studies",size=15,loc="left",horizontalalignment="right",fontweight ="bold")
        #axes[1].set_ylabel("Studies",size=15,loc="top")
        axes[2].set_title("Odds Ratio, [95% CI]",size=15,loc="right",fontweight ="bold")
    if meta==True:
        het_string1 = r'$Q={:.2f},Qdf={}$'.format(het_q_qdf_qp_i[0],het_q_qdf_qp_i[1])
        het_string2 = r'$p={:.1e},I^2={:.1f}%$'.format(het_q_qdf_qp_i[2],het_q_qdf_qp_i[3])
        axes[0].text(0,0,s=het_string1+"\n"+het_string2,size=15) 
    
    axes[1].set_yticks(df_to_plot["YORDER"])
    axes[1].set_yticklabels(df_to_plot.index,size=15)
    axes[1].get_yticklabels()[-1].set(fontweight = "bold")


    axes[0].get_shared_y_axes().join(axes[0], axes[1])
    axes[0].text(x=0,y=0.5,s=group_name,size=20)

    for index, row in df_to_plot.iterrows():
        fontweight ="normal"
        if index=="fixed effect":
            fontweight ="bold"
        axes[2].text(x=1,y=row["YORDER"],
                     s='{:.2f} [ {:.2f} , {:.2f} ]'.format(row["OR"],row["OR_95L"],row["OR_95U"]),
                     size=15,fontweight = fontweight,ha="right",va="center")
