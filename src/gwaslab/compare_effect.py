import sys
import gwaslab as gl
import os, psutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import scipy.stats as ss
import seaborn as sns

#20220310
def compare_effect(path1,
                   cols_name_list_1,
                   path2,
                   cols_name_list_2,
                   label1="Sumstats_1",
                   label2="Sumstats_2",
                   label3="Both",
                   label0="Not in any",
                   snplist=None,
                   anno=False,
                   verbose=False,
                   sig_level=5e-8,
                   scatterargs={"s":30},
                   reg=True,
                   legend_pos='upper left',
                   fontargs={'family':'sans','fontname':'Arial','fontsize':12},
                   errargs={"ecolor":"#cccccc","elinewidth":1}):
    
    #[snpid,chr,pos,p,ea,nea,effect,se]
    
    
    #snplist2
    if verbose: print("Loading "+label2+" SNP list...")
    sumstats = pd.read_table(path2,sep="\s+",usecols=[cols_name_list_2[0]])
    common_snp_set=set(sumstats[cols_name_list_2[0]])
    #snplist1

    #load sumstats1
    if verbose: print("Loading "+label1+"...")
    sumstats = pd.read_table(path1,sep="\s+",usecols=[cols_name_list_1[0],
                                                       cols_name_list_1[1],
                                                       cols_name_list_1[2],
                                                       cols_name_list_1[3]])
    #extract the common set
    common_snp_set = common_snp_set.intersection(sumstats[cols_name_list_1[0]])
    if verbose: print("Counting common snps : ",len(common_snp_set)," variants...")
    
    ### rename    
    sumstats.rename(columns={
                        cols_name_list_1[0]:"SNPID",
                        cols_name_list_1[1]:"CHR",
                        cols_name_list_1[2]:"POS",
                        cols_name_list_1[3]:"P"
                              },
                     inplace=True)
    
    ## 
    sumstats=sumstats.loc[sumstats["SNPID"].isin(common_snp_set),:]
    if verbose: print("Using only common snps : ",len(common_snp_set)," variants...")
    
    #extact SNPs for comparison 
    if verbose: print("Extract top snps from "+label1+"...")
    
    #if a snplist is provided, use the snp list
    if snplist: 
        if verbose: print("Extract snps in the given list from "+label1+"...")
        sig_list_1 = sumstats.loc[sumstats["SNPID"].isin(snplist),:]
    else:
        # otherwise use the sutomatically detected lead SNPs
        sig_list_1 = gl.getsig(sumstats,"SNPID","CHR","POS","P",
                               verbose=verbose,sig_level=sig_level)

    # do the same for sumstats 2
    if verbose: print("Loading "+label2+"...")
    sumstats = pd.read_table(path2,sep="\s+",usecols=[cols_name_list_2[0],
                                                       cols_name_list_2[1],
                                                       cols_name_list_2[2],
                                                       cols_name_list_2[3]])
    sumstats.rename(columns={
                        cols_name_list_2[0]:"SNPID",
                        cols_name_list_2[1]:"CHR",
                        cols_name_list_2[2]:"POS",
                        cols_name_list_2[3]:"P"
                              },
                     inplace=True)   
    ##
    
    sumstats=sumstats.loc[sumstats["SNPID"].isin(common_snp_set),:]
    
    if verbose: print("Using only common snps : ",len(common_snp_set)," variants...")
    
    if verbose: print("Extract top snps from "+label2+"...")
    if snplist: 
        if verbose: print("Extract snps in the given list from "+label2+"...")
        sig_list_2 = sumstats.loc[sumstats["SNPID"].isin(snplist),:]
    else: 
        sig_list_2 = gl.getsig(sumstats,"SNPID","CHR","POS","P",
                                 verbose=verbose,sig_level=sig_level)
    
    
    #Merge two list using SNPID
    ##############################################################################
    if verbose: print("Merging snps from "+label1+" and "+label2+"...")
    sig_list_merged = pd.merge(sig_list_1,sig_list_2,left_on="SNPID",right_on="SNPID",how="outer",suffixes=('_1', '_2'))
    
    ###############################################################################
    
    #### merging sumstats1
    if verbose: print("Extract the EFFECT_ALLELE, NON_EFFECT_ALLELE, EFFECT_SIZE,and SE of selected snps from "+label1+"...")
    sumstats = pd.read_table(path1,sep="\s+",
                              usecols=[cols_name_list_1[0],
                                       cols_name_list_1[4],
                                       cols_name_list_1[5],
                                       cols_name_list_1[6],
                                       cols_name_list_1[7]])
    sumstats.rename(columns={
                        cols_name_list_1[0]:"SNPID",
                        cols_name_list_1[4]:"EA_1",
                        cols_name_list_1[5]:"NEA_1",
                        cols_name_list_1[6]:"EFFECT_1",
                        cols_name_list_1[7]:"SE_1",
                              },
                     inplace=True)
    
    if verbose: print("Merging "+label1+" effect information...")
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")
    
    #### merging sumstats2
    if verbose: print("Extract the  EFFECT_ALLELE, NON_EFFECT_ALLELE, EFFECT_SIZE,and SE of selected snps from "+label2+"...")
    sumstats = pd.read_table(path2,sep="\s+",
                              usecols=[cols_name_list_2[0],
                                       cols_name_list_2[4],
                                       cols_name_list_2[5],
                                       cols_name_list_2[6],
                                       cols_name_list_2[7]])
    
    sumstats.rename(columns={
                        cols_name_list_2[0]:"SNPID",
                        cols_name_list_2[4]:"EA_2",
                        cols_name_list_2[5]:"NEA_2",
                        cols_name_list_2[6]:"EFFECT_2",
                        cols_name_list_2[7]:"SE_2"
                              },
                     inplace=True)
    if verbose: print("Merging "+label2+" effect information...")
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")
    
    sig_list_merged.set_index("SNPID",inplace=True)
    
    ###update sumstats1
    sumstats = pd.read_table(path1,sep="\s+",usecols=[cols_name_list_1[0],cols_name_list_1[3]])
    sumstats.rename(columns={
                        cols_name_list_1[0]:"SNPID",
                        cols_name_list_1[3]:"P_1"
                              },
                     inplace=True)
    sumstats.set_index("SNPID",inplace=True)
    sig_list_merged.update(sumstats)
    
    ###update sumstats2
    sumstats = pd.read_table(path2,sep="\s+",usecols=[cols_name_list_2[0],cols_name_list_2[3]])
    sumstats.rename(columns={
                        cols_name_list_2[0]:"SNPID",
                        cols_name_list_2[3]:"P_2"
                              },
                     inplace=True)
    sumstats.set_index("SNPID",inplace=True)
    sig_list_merged.update(sumstats)
    
    ####
    
    sig_list_merged["indicator"] = 0
    sig_list_merged.loc[sig_list_merged["P_1"]<sig_level,"indicator"]=1+sig_list_merged.loc[sig_list_merged["P_1"]<sig_level,"indicator"]
    sig_list_merged.loc[sig_list_merged["P_2"]<sig_level,"indicator"]=2+sig_list_merged.loc[sig_list_merged["P_2"]<sig_level,"indicator"]
    
    
    sig_list_merged["CHR"]=np.max(sig_list_merged[["CHR_1","CHR_2"]], axis=1).astype(int)
    sig_list_merged["POS"]=np.max(sig_list_merged[["POS_1","POS_2"]], axis=1).astype(int)
    sig_list_merged.drop(labels=['CHR_1', 'CHR_2','POS_1', 'POS_2'], axis=1,inplace=True)
    
    #### align allele effect with sumstats 1
    
    sig_list_merged["EA_2_aligned"]=sig_list_merged["EA_2"]
    sig_list_merged["NEA_2_aligned"]=sig_list_merged["NEA_2"]
    sig_list_merged["EFFECT_2_aligned"]=sig_list_merged["EFFECT_2"]

    sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2"]
    sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2"]
    sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EFFECT_2_aligned"]= -sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EFFECT_2"]
    
    save_path = label1+"_"+label2+"_beta_sig_list_merged.tsv"
    sig_list_merged.to_csv(save_path,"\t")
    
    sum0 = sig_list_merged.loc[sig_list_merged["indicator"]==0,:].dropna(axis=0)
    sum1only = sig_list_merged.loc[sig_list_merged["indicator"]==1,:].dropna(axis=0)
    sum2only = sig_list_merged.loc[sig_list_merged["indicator"]==2,:].dropna(axis=0)
    both     = sig_list_merged.loc[sig_list_merged["indicator"]==3,:].dropna(axis=0)
    
    if verbose: print("Identified "+str(len(sum0)) + " variants which are not significant in " + label3+".")
    if verbose: print("Identified "+str(len(sum1only)) + " variants which are only significant in " + label1+".")
    if verbose: print("Identified "+str(len(sum2only)) + " variants which are only significant in " + label2+".")
    if verbose: print("Identified "+str(len(both)) + " variants which are significant in " + label3 + ".")
    
    ##plot########################################################################################
    if verbose: print("Plotting the scatter plot for effect size comparison...")
    #plt.style.use("ggplot")
    sns.set_style("ticks")
    fig,ax = plt.subplots(figsize=(8,8)) 
    if len(sum0)>0:
        ax.errorbar(sum0["EFFECT_1"],sum0["EFFECT_2_aligned"], xerr=sum0["SE_1"],yerr=sum0["SE_2"],
                    linewidth=0,zorder=1,**errargs)
        ax.scatter(sum0["EFFECT_1"],sum0["EFFECT_2_aligned"],label=label0,zorder=2,color="#cccccc",marker="^",**scatterargs)
    if len(sum1only)>0:
        ax.errorbar(sum1only["EFFECT_1"],sum1only["EFFECT_2_aligned"], xerr=sum1only["SE_1"],yerr=sum1only["SE_2"],
                    linewidth=0,zorder=1,**errargs)
        ax.scatter(sum1only["EFFECT_1"],sum1only["EFFECT_2_aligned"],label=label1,zorder=2,color="#3eb489",marker="^",**scatterargs)
        
    if len(sum2only)>0:
        ax.errorbar(sum2only["EFFECT_1"],sum2only["EFFECT_2_aligned"], xerr=sum2only["SE_1"],yerr=sum2only["SE_2"],
                    linewidth=0,zorder=1,**errargs)
        ax.scatter(sum2only["EFFECT_1"],sum2only["EFFECT_2_aligned"],label=label2,zorder=2,color="#3264a8",marker="o",**scatterargs)

    if len(both)>0:
        ax.errorbar(both["EFFECT_1"],both["EFFECT_2_aligned"], xerr=both["SE_1"],yerr=both["SE_2"],
                    linewidth=0,zorder=1,**errargs)
        ax.scatter(both["EFFECT_1"],both["EFFECT_2_aligned"],label=label3,zorder=2,color="#fc8d14",marker="s",**scatterargs)
    
    ## annotation
    if anno==True:
        sig_list_toanno = sig_list_merged.dropna(axis=0)
        texts=[]
        for index, row in sig_list_toanno.iterrows():
            texts.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],index, ha='center', va='center')) 
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='grey'))
    
    elif type(anno) is dict:
        # if input is a dict
        sig_list_toanno = sig_list_toanno.loc[sig_list_toanno.index.isin(list(anno.keys())),:]
        texts=[]
        for index, row in sig_list_toanno.iterrows():
            texts.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],anno[index], ha='right', va='top')) 
        adjust_text(texts,ha='right', va='top',arrowprops=dict(arrowstyle='->', color='grey'))
        
    ax.set_xlabel("Per-allele effect size in "+label1,**fontargs)
    ax.set_ylabel("Per-allele effect size in "+label2,**fontargs)
    
    # plot x=0,y=0, and a 45 degree line
    xl,xh=ax.get_xlim()
    yl,yh=ax.get_ylim()
    ax.axline([min(xl,yl),min(xl,yl)], [max(xh, yh),max(xh, yh)],color='#cccccc', linestyle='--',zorder=1)
    ax.axhline(y=0, color='#cccccc', linestyle='--',zorder=1)
    ax.axvline(x=0, color='#cccccc', linestyle='--',zorder=1)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    if reg:
        reg = ss.linregress(sig_list_merged["EFFECT_1"],sig_list_merged["EFFECT_2_aligned"])
        if reg[0] > 0:
            ax.text(0.98,0.02,"y = "+"{:.2f}".format(reg[1]) +" + "+ "{:.2f}".format(reg[0])+"x ,p ="+ "{:.3e}".format(reg[3])+ ",$r^{2}$ =" +"{:.2f}".format(reg[2]),va="bottom",ha="right",transform=ax.transAxes,**fontargs)
        else:
            ax.text(0.98,0.02,"y = "+"{:.2f}".format(reg[1]) +" - "+ "{:.2f}".format(abs(reg[0]))+"x ,p ="+"{:.3e}".format(reg[3])+ ",$r^{2}$ =" +"{:.2f}".format(reg[2]),va="bottom",ha="right",transform=ax.transAxes,**fontargs)
            
        ax.axline(xy1=(0,reg[1]),slope=reg[0],color="#cccccc",zorder=1)
    
    L = ax.legend(title=r'$ P < 5 x 10^{-8}$ in:',loc=legend_pos)
    plt.setp(L.texts,**fontargs)
    plt.setp(L.get_title(),**fontargs)
    ##plot finished########################################################################################
    
    return [sig_list_merged, fig]