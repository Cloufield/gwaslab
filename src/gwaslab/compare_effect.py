import sys
import os, psutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns
import matplotlib.patches as mpatches
from gwaslab.getsig import getsig
from gwaslab.Log import Log
from gwaslab.winnerscurse import wc_correct
import gc

#20220422
def compare_effect(path1,
                   cols_name_list_1, effect_cols_list_1,
                   path2,
                   cols_name_list_2, effect_cols_list_2,
                   eaf=[],
                   maf_level=None,
                   label=["Sumstats_1","Sumstats_2","Both","None"],
                   snplist=None,
                   mode="beta",
                   anno=False,
                   wc_correction=False, 
                   null_beta=0,
                   is_q=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   wc_sig_level=5e-8,
                   # reg
                   reg_box=dict(boxstyle='round', facecolor='white', alpha=1,edgecolor="grey"),
                   is_reg=True,
                   allele_match=False,
                   r2_se=False,
                   is_45_helper_line=True,
                   legend_title=r'$ P < 5 x 10^{-8}$ in:',
                   legend_pos='upper left',
                   scatterargs={"s":20},
                   plt_args={"figsize":(8,8),"dpi":300},
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_args={"color":'black', "linestyle":'-',"lw":1},
                   fontargs={'fontsize':12},
                   # 'family':'sans','fontname':'Arial',
                   errargs={"ecolor":"#cccccc","elinewidth":1},
                   sep=["\t","\t"],
                   log = Log(),
                   verbose=False):
    
    #[snpid,p,ea,nea]      ,[effect,se]
    #[snpid,p,ea,nea,chr,pos],[effect,se]
    #[snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]
    
    if verbose: log.write("Start to process the raw sumstats for plotting...")
    
    ######### 1 check the value used to plot
    if mode not in ["Beta","beta","BETA","OR","or"]:
        raise ValueError("Please input Beta or OR")
    
    ######### 2 extract snplist2
    if verbose: log.write(" -Loading "+label[1]+" SNP list in memory...")    
    sumstats=pd.read_table(path2,sep=sep[1],usecols=[cols_name_list_2[0]])
    common_snp_set=set(sumstats[cols_name_list_2[0]].values)
    
    ######### 3 extract snplist1
    if snplist is not None:
        cols_to_extract = [cols_name_list_1[0],cols_name_list_1[1]]
    else:
        cols_to_extract = [cols_name_list_1[0],cols_name_list_1[1],cols_name_list_1[4],cols_name_list_1[5]]
 
    ######### 4 load sumstats1
    if verbose: log.write(" -Loading sumstats for "+label[0]+":",",".join(cols_to_extract))
    sumstats = pd.read_table(path1,sep=sep[0],usecols=cols_to_extract)
    
    ######### 5 extract the common set
    common_snp_set = common_snp_set.intersection(sumstats[cols_name_list_1[0]].values)
    if verbose: log.write(" -Counting  variants available for both datasets:",len(common_snp_set)," variants...")
    
    ######### 6 rename the sumstats
    rename_dict = { cols_name_list_1[0]:"SNPID",
               cols_name_list_1[1]:"P",
               }
    
    if snplist is None: 
        rename_dict[cols_name_list_1[4]]="CHR"
        rename_dict[cols_name_list_1[5]]="POS"
    
    sumstats.rename(columns=rename_dict,inplace=True)
    
    ######### 7 exctract only available variants from sumstats1 
    sumstats = sumstats.loc[sumstats["SNPID"].isin(common_snp_set),:]
    
    if verbose: log.write(" -Using only variants available for both datasets...")
    ######### 8 extact SNPs for comparison 
    
    if snplist is not None: 
        ######### 8.1 if a snplist is provided, use the snp list
        if verbose: log.write(" -Extract variants in the given list from "+label[0]+"...")
        sig_list_1 = sumstats.loc[sumstats["SNPID"].isin(snplist),:].copy()
    else:
        ######### 8,2 otherwise use the sutomatically detected lead SNPs
        if verbose: log.write(" -Extract lead variants from "+label[0]+"...")
        sig_list_1 = getsig(sumstats,"SNPID","CHR","POS","P", verbose=verbose,sig_level=sig_level)
    
    sig_list_1 = drop_duplicate_and_na(sig_list_1, sort_by="P", log=log ,verbose=verbose)

    ######### 9 extract snplist2
    if snplist is not None:
        cols_to_extract = [cols_name_list_2[0],cols_name_list_2[1]]
    else:
        cols_to_extract = [cols_name_list_2[0],cols_name_list_2[1],cols_name_list_2[4],cols_name_list_2[5]]
    
    if verbose: log.write(" -Loading sumstats for "+label[1]+":",",".join(cols_to_extract))
    sumstats = pd.read_table(path2,sep=sep[1],usecols=cols_to_extract)
    
    ######### 10 rename sumstats2
    rename_dict = { cols_name_list_2[0]:"SNPID",
                    cols_name_list_2[1]:"P",
                }
    if snplist is None: 
        rename_dict[cols_name_list_2[4]]="CHR"
        rename_dict[cols_name_list_2[5]]="POS"
    sumstats.rename(columns=rename_dict,inplace=True)
    
    ######### 11 exctract only overlapping variants from sumstats2
    sumstats = sumstats.loc[sumstats["SNPID"].isin(common_snp_set),:]
    
    ######## 12 extact SNPs for comparison 
    if snplist: 
        ######### 12.1 if a snplist is provided, use the snp list
        if verbose: log.write(" -Extract snps in the given list from "+label[1]+"...")
        sig_list_2 = sumstats.loc[sumstats["SNPID"].isin(snplist),:].copy()
    else: 
        if verbose: log.write(" -Extract lead snps from "+label[1]+"...")
        ######### 12.2 otherwise use the sutomatically detected lead SNPs
        sig_list_2 = getsig(sumstats,"SNPID","CHR","POS","P",
                                 verbose=verbose,sig_level=sig_level)

    sig_list_2 = drop_duplicate_and_na(sig_list_2, sort_by="P", log=log ,verbose=verbose)

    ######### 13 Merge two list using SNPID
    ##############################################################################
    if verbose: log.write("Merging snps from "+label[0]+" and "+label[1]+"...")
    
    sig_list_merged = pd.merge(sig_list_1,sig_list_2,left_on="SNPID",right_on="SNPID",how="outer",suffixes=('_1', '_2'))
    #     SNPID       P_1       P_2
    #0   rs117986209  0.142569  0.394455
    #1     rs6704312  0.652104  0.143750

    ###############################################################################

    ########## 14 Merging sumstats1
    
    if mode=="beta" or mode=="BETA" or mode=="Beta":
         #[snpid,p,ea,nea]      ,[effect,se]
        #[snpid,p,ea,nea,chr,pos],[effect,se]
        #[snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]
        cols_to_extract = [cols_name_list_1[0],cols_name_list_1[1], cols_name_list_1[2],cols_name_list_1[3], effect_cols_list_1[0], effect_cols_list_1[1]]
    else:
        cols_to_extract = [cols_name_list_1[0],cols_name_list_1[1], cols_name_list_1[2],cols_name_list_1[3], effect_cols_list_1[0], effect_cols_list_1[1], effect_cols_list_1[2]]
    
    if len(eaf)>0: cols_to_extract.append(eaf[0])   
    if verbose: log.write(" -Extract statistics of selected variants from "+label[0]+" : ",",".join(cols_to_extract) )
    sumstats = pd.read_table(path1,sep=sep[0],usecols=cols_to_extract)


    if mode=="beta" or mode=="BETA" or mode=="Beta":
        rename_dict = { cols_name_list_1[0]:"SNPID",
                        cols_name_list_1[1]:"P_1",
                        cols_name_list_1[2]:"EA_1",
                        cols_name_list_1[3]:"NEA_1",
                        effect_cols_list_1[0]:"EFFECT_1",
                        effect_cols_list_1[1]:"SE_1",
    }
        
    else:
        # if or
        rename_dict = { cols_name_list_1[0]:"SNPID",
                        cols_name_list_1[1]:"P_1",
                        cols_name_list_1[2]:"EA_1",
                        cols_name_list_1[3]:"NEA_1",
                        effect_cols_list_1[0]:"OR_1",
                        effect_cols_list_1[1]:"OR_L_1",
                        effect_cols_list_1[2]:"OR_H_1"
    }
    ## check if eaf column is provided.
    if len(eaf)>0: rename_dict[eaf[0]]="EAF_1"
    sumstats.rename(columns=rename_dict, inplace=True)
    
    # drop na and duplicate
    sumstats = drop_duplicate_and_na(sumstats,  sort_by="P_1", log=log , verbose=verbose)
    sumstats.drop("P_1",axis=1,inplace=True)

    if verbose: log.write(" -Merging "+label[0]+" effect information...")
    
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")

    ############ 15 merging sumstats2
    
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        cols_to_extract = [cols_name_list_2[0],cols_name_list_2[1],cols_name_list_2[2],cols_name_list_2[3], effect_cols_list_2[0], effect_cols_list_2[1]]
    else:
        # if or
        cols_to_extract = [cols_name_list_2[0],cols_name_list_2[1],cols_name_list_2[2],cols_name_list_2[3], effect_cols_list_2[0], effect_cols_list_2[1], effect_cols_list_2[2]]
    ## check if eaf column is provided.
    if len(eaf)>0: cols_to_extract.append(eaf[1])
    
    if verbose: log.write(" -Extract statistics of selected variants from "+label[1]+" : ",",".join(cols_to_extract) )
    sumstats = pd.read_table(path2,sep=sep[1],usecols=cols_to_extract)
    gc.collect()
    if mode=="beta" or mode=="BETA" or mode=="Beta":
          rename_dict = { cols_name_list_2[0]:"SNPID",
                        cols_name_list_2[1]:"P_2",
                        cols_name_list_2[2]:"EA_2",
                        cols_name_list_2[3]:"NEA_2",
                        effect_cols_list_2[0]:"EFFECT_2",
                        effect_cols_list_2[1]:"SE_2",
    }
    else:
                    rename_dict = { cols_name_list_2[0]:"SNPID",
                        cols_name_list_2[1]:"P_2",
                        cols_name_list_2[2]:"EA_2",
                        cols_name_list_2[3]:"NEA_2",
                        effect_cols_list_2[0]:"OR_2",
                        effect_cols_list_2[1]:"OR_L_2",
                        effect_cols_list_2[2]:"OR_H_2"
    }
    if len(eaf)>0: rename_dict[eaf[1]]="EAF_2"
    sumstats.rename(columns=rename_dict, inplace=True)         
    # drop na and duplicate
    sumstats = drop_duplicate_and_na(sumstats, sort_by="P_2", log=log, verbose=verbose)
    sumstats.drop("P_2",axis=1,inplace=True)

    if verbose: log.write(" -Merging "+label[1]+" effect information...")
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")
    
    sig_list_merged.set_index("SNPID",inplace=True)

    ################ 16 update sumstats1
    if verbose: log.write(" -Updating missing information for "+label[0]+" ...")
    sumstats = pd.read_table(path1,sep=sep[0],usecols=[cols_name_list_1[0],cols_name_list_1[1]])
    sumstats.rename(columns={
                        cols_name_list_1[0]:"SNPID",
                        cols_name_list_1[1]:"P_1"
                              },
                     inplace=True)
    # drop na and duplicate
    sumstats = drop_duplicate_and_na(sumstats, sort_by="P_1", log=log, verbose=verbose)
    
    sumstats.set_index("SNPID",inplace=True)
    sig_list_merged.update(sumstats)
    
    ################# 17 update sumstats2
    if verbose: log.write(" -Updating missing information for "+label[1]+" ...")
    sumstats = pd.read_table(path2,sep=sep[1],usecols=[cols_name_list_2[0],cols_name_list_2[1]])
    sumstats.rename(columns={
                        cols_name_list_2[0]:"SNPID",
                        cols_name_list_2[1]:"P_2"
                              },
                     inplace=True)
    # drop na and duplicate
    sumstats = drop_duplicate_and_na(sumstats, sort_by="P_2", log=log, verbose=verbose)              
    
    sumstats.set_index("SNPID",inplace=True)
    sig_list_merged.update(sumstats)

    ####
#################################################################################
    ############## 18 init indicator
    if verbose: log.write(" -Assigning indicator  ...")
    # 0-> 0
    # 1 -> sig in sumstats1
    # 2 -> sig in sumsatts2
    # 3->  sig in both sumstats1 + sumstats2
    sig_list_merged["indicator"] = 0
    sig_list_merged.loc[sig_list_merged["P_1"]<sig_level,"indicator"]=1+sig_list_merged.loc[sig_list_merged["P_1"]<sig_level,"indicator"]
    sig_list_merged.loc[sig_list_merged["P_2"]<sig_level,"indicator"]=2+sig_list_merged.loc[sig_list_merged["P_2"]<sig_level,"indicator"]
    
    if snplist is None:
        sig_list_merged["CHR"]=np.max(sig_list_merged[["CHR_1","CHR_2"]], axis=1).astype(int)
        sig_list_merged["POS"]=np.max(sig_list_merged[["POS_1","POS_2"]], axis=1).astype(int)
        sig_list_merged.drop(labels=['CHR_1', 'CHR_2','POS_1', 'POS_2'], axis=1,inplace=True)
    
    if verbose: log.write(" -Aligning "+label[1]+" EA with "+label[0]+" EA ...")
    ############### 19 align allele effect with sumstats 1
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        # copy raw
        sig_list_merged["EA_2_aligned"]=sig_list_merged["EA_2"]
        sig_list_merged["NEA_2_aligned"]=sig_list_merged["NEA_2"]
        sig_list_merged["EFFECT_2_aligned"]=sig_list_merged["EFFECT_2"]
        
        #filp ea/nea and beta for sumstats2
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EFFECT_2_aligned"]= -sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EFFECT_2"]
    else:
        #flip for OR or - +
        sig_list_merged["OR_L_1"]=np.abs(sig_list_merged["OR_L_1"]-sig_list_merged["OR_1"])
        sig_list_merged["OR_H_1"]=np.abs(sig_list_merged["OR_H_1"]-sig_list_merged["OR_1"])
        sig_list_merged["OR_L_2"]=np.abs(sig_list_merged["OR_L_2"]-sig_list_merged["OR_2"])
        sig_list_merged["OR_H_2"]=np.abs(sig_list_merged["OR_H_2"]-sig_list_merged["OR_2"])
        
        sig_list_merged["EA_2_aligned"]=sig_list_merged["EA_2"]
        sig_list_merged["NEA_2_aligned"]=sig_list_merged["NEA_2"]
        sig_list_merged["OR_2_aligned"]=sig_list_merged["OR_2"]
        sig_list_merged["OR_L_2_aligned"]=sig_list_merged["OR_L_2"]
        sig_list_merged["OR_H_2_aligned"]=sig_list_merged["OR_H_2"]

        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_2_aligned"]= 1/sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_L_2_aligned"]= 1/sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_L_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_H_2_aligned"]= 1/sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_H_2"]
    
    if len(eaf)>0:
        # flip eaf
        sig_list_merged["EAF_2_aligned"]=sig_list_merged["EAF_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EAF_2_aligned"]= 1 -sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EAF_2"]
    
    # checking effect allele matching
    nonmatch = np.nansum(sig_list_merged["EA_1"] != sig_list_merged["EA_2_aligned"])
    if verbose: log.write(" -Aligned all EAs in {} with EAs in {} ...".format(label[1],label[0]))
    if nonmatch>0:
        if verbose: log.write(" -Warning: Alleles for {} variants do not match...".format(nonmatch))
    if allele_match==True:
        if nonmatch>0:
            sig_list_merged = sig_list_merged.loc[sig_list_merged["EA_1"] == sig_list_merged["EA_2_aligned"]]
        else:
            if verbose: log.write(" -No variants with EA not matching...")
    ####################################################################################################################################
    ## winner's curse correction using aligned beta
    if wc_correction == True:
        if verbose: log.write(" -Correcting BETA for winner's curse with threshold at {}...".format(sig_level))
        sig_list_merged["EFFECT_1_RAW"] = sig_list_merged["EFFECT_1"].copy()
        sig_list_merged["EFFECT_2_aligned_RAW"] = sig_list_merged["EFFECT_2_aligned"].copy()
        sig_list_merged["EFFECT_1"] = sig_list_merged[["EFFECT_1_RAW","SE_1"]].apply(lambda x: wc_correct(x[0],x[1],sig_level),axis=1)
        sig_list_merged["EFFECT_2_aligned"] = sig_list_merged[["EFFECT_2_aligned_RAW","SE_2"]].apply(lambda x: wc_correct(x[0],x[1],sig_level),axis=1)
    ########################## Het test############################################################
    ## heterogeneity test
    if (is_q is True) and (mode=="beta" or mode=="BETA" or mode=="Beta"):
        if verbose: log.write(" -Calculating Cochran's Q statistics and peform chisq test...")
        sig_list_merged = test_q(sig_list_merged,"EFFECT_1","SE_1","EFFECT_2_aligned","SE_2")
        
    ######################### save ###############################################################
    ## save the merged data
    save_path = label[0]+"_"+label[1]+"_beta_sig_list_merged.tsv"
    if verbose: log.write(" -Saving the merged data to:",save_path)
    sig_list_merged.to_csv(save_path,"\t")
    
    ########################## maf_threshold#############################################################
    if (len(eaf)>0) and (maf_level is not None):
        both_eaf_clear =  (sig_list_merged["EAF_1"]>maf_level)&(sig_list_merged["EAF_1"]<1-maf_level)&(sig_list_merged["EAF_2"]>maf_level)&(sig_list_merged["EAF_2"]<1-maf_level)
        if verbose: log.write(" -Exclude "+str(len(sig_list_merged) -sum(both_eaf_clear))+ " variants with maf <",maf_level)
        sig_list_merged = sig_list_merged.loc[both_eaf_clear,:]
    # heterogeneity summary
    if (is_q is True) and (mode=="beta" or mode=="BETA" or mode=="Beta"):
        if verbose: log.write(" -Significant het:" ,len(sig_list_merged.loc[sig_list_merged["HetP"]<0.05,:]))
        if verbose: log.write(" -All sig:" ,len(sig_list_merged))
        if verbose: log.write(" -Het rate:" ,len(sig_list_merged.loc[sig_list_merged["HetP"]<0.05,:])/len(sig_list_merged))   
    
    # extract group
    sum0 = sig_list_merged.loc[sig_list_merged["indicator"]==0,:].dropna(axis=0)
    sum1only = sig_list_merged.loc[sig_list_merged["indicator"]==1,:].dropna(axis=0)
    sum2only = sig_list_merged.loc[sig_list_merged["indicator"]==2,:].dropna(axis=0)
    both     = sig_list_merged.loc[sig_list_merged["indicator"]==3,:].dropna(axis=0)
    
    if verbose: log.write(" -Identified "+str(len(sum0)) + " variants which are not significant in " + label[3]+".")
    if verbose: log.write(" -Identified "+str(len(sum1only)) + " variants which are only significant in " + label[0]+".")
    if verbose: log.write(" -Identified "+str(len(sum2only)) + " variants which are only significant in " + label[1]+".")
    if verbose: log.write(" -Identified "+str(len(both)) + " variants which are significant in " + label[2] + ".")
    
    ##plot########################################################################################
    if verbose: log.write("Creating the scatter plot for effect sizes comparison...")
    #plt.style.use("ggplot")
    sns.set_style("ticks")
    fig,ax = plt.subplots(**plt_args) 
    legend_elements=[]
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        if len(sum0)>0:
            ax.errorbar(sum0["EFFECT_1"],sum0["EFFECT_2_aligned"], xerr=sum0["SE_1"],yerr=sum0["SE_2"],
                        linewidth=0,zorder=1,**errargs)
            
            ax.scatter(sum0["EFFECT_1"],sum0["EFFECT_2_aligned"],label=label[3],zorder=2,color="#cccccc",edgecolors=sum0["Edge_color"],marker=".",**scatterargs)
            #legend_elements.append(mpatches.Circle(facecolor='#cccccc', edgecolor='white', label=label[3]))
        if len(sum1only)>0:
            ax.errorbar(sum1only["EFFECT_1"],sum1only["EFFECT_2_aligned"], xerr=sum1only["SE_1"],yerr=sum1only["SE_2"],
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(sum1only["EFFECT_1"],sum1only["EFFECT_2_aligned"],label=label[0],zorder=2,color="#e6320e",edgecolors=sum1only["Edge_color"],marker="^",**scatterargs)
            #legend_elements.append(mpatches.Patch(facecolor='#e6320e', edgecolor='white', label=label[0]))

        if len(sum2only)>0:
            ax.errorbar(sum2only["EFFECT_1"],sum2only["EFFECT_2_aligned"], xerr=sum2only["SE_1"],yerr=sum2only["SE_2"],
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(sum2only["EFFECT_1"],sum2only["EFFECT_2_aligned"],label=label[1],zorder=2,color="#41e620",edgecolors=sum2only["Edge_color"],marker="o",**scatterargs)
            #legend_elements.append(mpatches.Circle(facecolor='#41e620', edgecolor='white', label=label[1]))

        if len(both)>0:
            ax.errorbar(both["EFFECT_1"],both["EFFECT_2_aligned"], xerr=both["SE_1"],yerr=both["SE_2"],
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(both["EFFECT_1"],both["EFFECT_2_aligned"],label=label[2],zorder=2,color="#205be6",edgecolors=both["Edge_color"],marker="s",**scatterargs)  
            #legend_elements.append(mpatches.Patch(facecolor='#205be6', edgecolor='white', label=label[2]))
    else:
        ## if OR
        if len(sum0)>0:
            ax.errorbar(sum0["OR_1"],sum0["OR_2_aligned"], xerr=sum0[["OR_L_1","OR_H_1"]].T,yerr=sum0[["OR_L_2_aligned","OR_H_2_aligned"]].T,
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(sum0["OR_1"],sum0["OR_2_aligned"],label=label[3],zorder=2,color="#cccccc",marker=".",**scatterargs)
        if len(sum1only)>0:
            ax.errorbar(sum1only["OR_1"],sum1only["OR_2_aligned"], xerr=sum1only[["OR_L_1","OR_H_1"]].T,yerr=sum1only[["OR_L_2_aligned","OR_H_2_aligned"]].T,
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(sum1only["OR_1"],sum1only["OR_2_aligned"],label=label[0],zorder=2,color="#205be6",marker="^",**scatterargs)

        if len(sum2only)>0:
            ax.errorbar(sum2only["OR_1"],sum2only["OR_2_aligned"], xerr=sum2only[["OR_L_1","OR_H_1"]].T,yerr=sum2only[["OR_L_2_aligned","OR_H_2_aligned"]].T,
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(sum2only["OR_1"],sum2only["OR_2_aligned"],label=label[1],zorder=2,color="#41e620",marker="o",**scatterargs)

        if len(both)>0:
            ax.errorbar(both["OR_1"],both["OR_2_aligned"], xerr=both[["OR_L_1","OR_H_1"]].T,yerr=both[["OR_L_2_aligned","OR_H_2_aligned"]].T,
                        linewidth=0,zorder=1,**errargs)
            ax.scatter(both["OR_1"],both["OR_2_aligned"],label=label[2],zorder=2,color="#e6320e",marker="s",**scatterargs)
    
    ## annotation #################################################################################################################
    if anno==True:
        from adjustText import adjust_text
        sig_list_toanno = sig_list_merged.dropna(axis=0)
        texts=[]
        for index, row in sig_list_toanno.iterrows():
            if mode=="beta" or mode=="BETA" or mode=="Beta":
                texts.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],index, ha='center', va='center')) 
            else:
                texts.append(plt.text(row["OR_1"], row["OR_2_aligned"],index, ha='center', va='center')) 
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='grey'))
    
    elif type(anno) is dict:
        # if input is a dict
        sig_list_toanno = sig_list_toanno.loc[sig_list_toanno.index.isin(list(anno.keys())),:]
        texts=[]
        for index, row in sig_list_toanno.iterrows():
            texts.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],anno[index], ha='right', va='top')) 
        adjust_text(texts,ha='right', va='top',arrowprops=dict(arrowstyle='->', color='grey'))
    #################################################################################################################################
    
    # plot x=0,y=0, and a 45 degree line
    xl,xh=ax.get_xlim()
    yl,yh=ax.get_ylim()
    
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        #if using beta
        ax.axhline(y=0, zorder=1,**helper_line_args)
        ax.axvline(x=0, zorder=1,**helper_line_args)
    else:
        #if using OR
        ax.axhline(y=1, zorder=1,**helper_line_args)
        ax.axvline(x=1, zorder=1,**helper_line_args)
    
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    

    ###regression line##############################################################################################################################
    if len(sig_list_merged)<3: is_reg=False
    if is_reg is True:
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            reg = ss.linregress(sig_list_merged["EFFECT_1"],sig_list_merged["EFFECT_2_aligned"])
            
            # estimate se for r2
            if r2_se==True:
                if verbose:log.write(" -Estimating SE for Rsq using Jackknife method.")
                r2_se_jackknife = jackknife_r2(sig_list_merged)
                r2_se_jackknife_string = " ({:.2f})".format(r2_se_jackknife)
            else:
                r2_se_jackknife_string= ""
        else:
            reg = ss.linregress(sig_list_merged["OR_1"],sig_list_merged["OR_2_aligned"])
        #### calculate p values based on selected value , default = 0 
        if verbose:log.write(" -Calculating p values based on given null slope :",null_beta)
        t_score = (reg[0]-null_beta) / reg[4]
        degree = len(sig_list_merged.dropna())-2
        p = ss.t.sf(abs(t_score), df=degree)*2
        if verbose:log.write(" -Beta_se = ", reg[4])
        if verbose:log.write(" -H0 beta = ", null_beta, ", recalculated p = ", "{:.2e}".format(p))
        if verbose:log.write(" -H0 beta =  0",", default p = ", "{:.2e}".format(reg[3]))
        if r2_se==True:
            if verbose:log.write(" -R2 se (jackknife) = {:.2e}".format(r2_se_jackknife))

        if reg[0] > 0:
            #if regression coeeficient >0 : auxiliary line slope = 1
            if is_45_helper_line is True:
                ax.axline([min(xl,yl),min(xl,yl)], [max(xh, yh),max(xh, yh)],zorder=1,**helper_line_args)

            #add text
            p12=str("{:.2e}".format(p)).split("e")[0]
            pe =str(int("{:.2e}".format(p).split("e")[1]))
            p_text="$p = " + p12 + " \\times  10^{"+pe+"}$"
            p_latex= f'{p_text}'
            ax.text(0.98,0.02,"$y =$ "+"{:.2f}".format(reg[1]) +" $+$ "+ "{:.2f}".format(reg[0])+" $x$, "+ p_latex + ", $r^{2} =$" +"{:.2f}".format(reg[2])+r2_se_jackknife_string, va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **fontargs)
        else:
            #if regression coeeficient <0 : auxiliary line slope = -1
            if is_45_helper_line is True:
                if mode=="beta" or mode=="BETA" or mode=="Beta": 
                    ax.axline([min(xl,yl),-min(xl,yl)], [max(xh, yh),-max(xh, yh)],zorder=1,**helper_line_args)
                else:
                    ax.axline([min(xl,yl),-min(xl,yl)], [max(xh, yh),-max(xh, yh)],zorder=1,**helper_line_args)
            #add text
            p12=str("{:.2e}".format(p)).split("e")[0]
            pe =str(int("{:.2e}".format(p).split("e")[1]))
            p_text="$p = " + p12 + " \\times  10^{"+pe+"}$"
            p_latex= f'{p_text}'
            ax.text(0.98,0.02,"$y =$ "+"{:.2f}".format(reg[1]) +" $-$ "+ "{:.2f}".format(abs(reg[0]))+" $x$, "+ p_latex + ", $r^{2} =$" +"{:.2f}".format(reg[2])+r2_se_jackknife_string, va="bottom",ha="right",transform=ax.transAxes,bbox=reg_box,**fontargs)
            
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            middle = sig_list_merged["EFFECT_1"].mean()
        else:
            middle = sig_list_merged["OR_1"].mean()
        
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            ax.axline(xy1=(0,reg[1]),slope=reg[0],color="#cccccc",linestyle='--',zorder=1)
        else:
            ax.axline(xy1=(1,reg[0]+reg[1]),slope=reg[0],color="#cccccc",linestyle='--',zorder=1)
        
    
    ax.set_xlabel(xylabel_prefix+label[0],**fontargs)
    ax.set_ylabel(xylabel_prefix+label[1],**fontargs)
    
    L = ax.legend(title=legend_title,loc=legend_pos,framealpha=1,edgecolor="grey")
    
    for i, handle in enumerate(L.legendHandles):
        handle.set_edgecolor("white")
    
    plt.setp(L.texts,**fontargs)
    plt.setp(L.get_title(),**fontargs)
    ##plot finished########################################################################################
    gc.collect()
    return [sig_list_merged, fig,log]



################################################################################################################################
def plotdaf(sumstats,
             eaf="EAF",
             daf="DAF",
             scatter_args={"s":1},
             threshold=0.16,
             is_reg=True,
             is_45_helper_line=True,
             is_threshold=True,
             helper_line_args={"color":'black', "linestyle":'-',"lw":1},
             threshold_line_args={"color":'#cccccc', "linestyle":'dotted'},
             reg_line_args={"color":'#cccccc', "linestyle":'--'},
             plt_args={"figsize":(8,4),"dpi":300},
            histplot_args={"log_scale":(False,True)},
            fontargs={'family':'sans','fontname':'Arial','fontsize':8},
            verbose=True,
            log=Log()
           ):
    
    
    if verbose: log.write("Start to plot Reference frequency vs Effect allele frequency plot...")
    if not ((eaf in sumstats.columns) and (daf in sumstats.columns)):
        raise ValueError("EAF and/or DAF columns were not detected.")
    
    
    sumstats = sumstats.loc[(~sumstats[eaf].isna())&(~sumstats[daf].isna()),[eaf,daf]].copy()
    sumstats.loc[:,daf] = sumstats.loc[:,daf].astype("float")
    sumstats.loc[:,eaf] = sumstats.loc[:,eaf].astype("float")
    if verbose: log.write(" -Plotting valriants:" + str(len(sumstats)))
    
    sumstats.loc[:,"RAF"]=sumstats[eaf] - sumstats[daf]
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(1, 2,**plt_args)
    ax1.scatter(sumstats["RAF"],sumstats[eaf],**scatter_args)
    
    if is_reg is True:
        if verbose: log.write(" -Plotting regression line...")
        reg = ss.linregress(sumstats["RAF"],sumstats[eaf])
        if verbose:log.write(" -Beta = ", reg[0])
        if verbose:log.write(" -Intercept = ", reg[1])
        if verbose:log.write(" -R2 = ", reg[2])
        ax1.axline(xy1=(0,reg[1]),slope=reg[0],zorder=1,**reg_line_args)
    if is_threshold is True:
        if verbose: log.write(" -Threshold : " + str(threshold))
        num = sum(np.abs(sumstats[daf])>threshold )
        if verbose: log.write(" -Variants with relatively large DAF : ",num )
        if verbose: log.write(" -Percentage for variants with relatively large DAF : ",num/len(sumstats) )
        ax1.axline(xy1=(0,threshold),slope=1,zorder=1,**threshold_line_args)
        ax1.axline(xy1=(threshold,0),slope=1,zorder=1,**threshold_line_args)
    xl,xh=ax1.get_xlim()
    yl,yh=ax1.get_ylim()
    if is_45_helper_line is True:
        ax1.axline([0,0], [1,1],zorder=1, **helper_line_args)
    ax1.set_xlabel("Alternative Allele Frequency in Reference Population (RAF)",**fontargs)
    ax1.set_ylabel("Effect Allele Frequency in Sumstats (EAF)",**fontargs)
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    
    sumstats.loc[:,"ID"] = sumstats.index
    to_plot = pd.melt(sumstats,id_vars=['ID'], value_vars=['EAF',"RAF"], var_name='Types', value_name='Allele Frequency')
    
    sns.histplot(data=to_plot, x="Allele Frequency", hue="Types", fill=True, ax=ax2,**histplot_args)
    ax2.set_xlabel("Allele Frequency",**fontargs)
    plt.tight_layout()
    return fig
    
def test_q(df,beta1,se1,beta2,se2):
    w1="Weight_1"
    w2="Weight_2"
    beta="BETA_FE"
    q="Q"
    pq="HetP"
    i2="I2"
    df[w1]=1/(df[se1])**2
    df[w2]=1/(df[se2])**2
    df[beta] =(df[w1]*df[beta1] + df[w2]*df[beta2])/(df[w1]+df[w2])
    
    # Cochran(1954)
    df[q] = df[w1]*(df[beta1]-df[beta])**2 + df[w2]*(df[beta2]-df[beta])**2
    df[pq] = ss.chi2.sf(df[q], 1)
    df["Edge_color"]="white"
    df.loc[df[pq]<0.05,"Edge_color"]="black"
    df.drop(columns=["Weight_1","Weight_2","BETA_FE"],inplace=True)
    # Huedo-Medina, T. B., Sánchez-Meca, J., Marín-Martínez, F., & Botella, J. (2006). Assessing heterogeneity in meta-analysis: Q statistic or I² index?. Psychological methods, 11(2), 193.
    
    # calculate I2
    df[i2] = (df[q] - 1)/df[q]
    df.loc[df[i2]<0,i2] = 0 
    
    return df

def jackknife_r2(df,x="EFFECT_1",y="EFFECT_2_aligned"):
    """Jackknife estimation of se for rsq

    """

    # dropna
    df_nona = df.loc[:,[x,y]].dropna()
    
    # non-empty entries
    n=len(df)
    
    # assign row number
    df_nona["nrow"] = range(n)
    
    # a list to store r2
    r2_list=[]
    
    # estimate r2
    for i in range(n):
        # exclude 1 record
        records_to_use = df_nona["nrow"]!=i
        # estimate r2
        reg_jackknife = ss.linregress(df_nona.loc[records_to_use, x],df_nona.loc[records_to_use,y])
        # add r2_i to list
        r2_list.append(reg_jackknife[2])

    # convert list to array
    r2s = np.array(r2_list)
    # https://en.wikipedia.org/wiki/Jackknife_resampling
    r2_se = np.sqrt( (n-1)/n * np.sum((r2s - np.mean(r2s))**2) )
    return r2_se

def drop_duplicate_and_na(df,snpid="SNPID",sort_by=False,log=Log(),verbose=True):
    length_before = len(df)
    if sort_by!=False:
        df.sort_values(by = sort_by, inplace=True)
    df.dropna(axis="index",subset=[snpid],inplace=True)
    df.drop_duplicates(subset=[snpid], keep='first', inplace=True) 
    length_after= len(df)
    if length_before !=  length_after:
        if verbose: log.write(" -Dropped {} duplicates or NAs...".format(length_before - length_after))
    return df