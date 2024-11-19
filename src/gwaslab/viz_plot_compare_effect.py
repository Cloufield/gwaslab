import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns
import gc
import math
import scipy.stats as ss
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from gwaslab.viz_aux_save_figure import save_figure
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_get_sig import annogene
from gwaslab.g_Log import Log
from gwaslab.util_in_correct_winnerscurse import wc_correct
from gwaslab.util_in_correct_winnerscurse import wc_correct_test
from gwaslab.g_Sumstats import Sumstats
from gwaslab.io_process_args import _merge_and_sync_dic
from gwaslab.io_process_args import _extract_kwargs
#20220422
def compare_effect(path1,
                   path2,
                   cols_name_list_1=None, effect_cols_list_1=None,
                   cols_name_list_2=None, effect_cols_list_2=None,
                   eaf=[],
                   maf_level=None,
                   label=None,
                   snplist=None,
                   mode="beta",
                   anno=False,
                   anno_het=False,
                   anno_min=0,
                   anno_min1=0,
                   anno_min2=0,
                   anno_diff=0,
                   anno_args=None,
                   scaled=False,
                   scaled1=False,
                   scaled2=False,
                   wc_correction=False, 
                   null_beta=0,
                   is_q=False,
                   is_q_mc = False,
                   include_all=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   get_lead_args=None,
                   drop=False,
                   wc_sig_level=5e-8,
                   # reg
                   reg_box=None,
                   is_reg=True,
                   fdr=False,
                   allele_match=False,
                   r_se=False,
                   is_45_helper_line=True,
                   legend_mode="full",
                   legend_title=r'$ P < 5 x 10^{-8}$ in:',
                   legend_title2=r'Heterogeneity test:',
                   legend_pos='upper left',
                   scatterargs=None,
                   plt_args=None,
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_args=None,
                   adjust_text_kwargs = None,
                   adjust_text_kwargs_l = None,
                   adjust_text_kwargs_r = None,
                   font_args=None,
                   fontargs=None,
                   build="19",
                   r_or_r2="r",
                   errargs=None,
                   legend_args=None,
                   sep=["\t","\t"],
                   log = Log(),
                   save=False,
                   save_args=None,
                   verbose=False,
                   **kwargs):

    #[snpid,p,ea,nea]      ,[effect,se]
    #[snpid,p,ea,nea,chr,pos],[effect,se]
    #[snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]
    if scaled == True:
        scaled1 = True
        scaled2 = True

    if legend_title== r'$ P < 5 x 10^{-8}$ in:' and sig_level!=5e-8:
        
        exponent = math.floor(math.log10(sig_level))
        mantissa = sig_level / 10**exponent

        legend_title =  '$ P < {} x 10^{{{}}}$ in:'.format(mantissa, exponent)

    if is_q_mc=="fdr" or is_q_mc=="bon":
        is_q = True
    if is_q == True:
        if is_q_mc not in [False,"fdr","bon","non"]:
            raise ValueError('Please select either "fdr" or "bon" or "non"/False for is_q_mc.')
    if save_args is None:
        save_args = {"dpi":300,"facecolor":"white"}
    if reg_box is None:
        reg_box = dict(boxstyle='round', facecolor='white', alpha=1,edgecolor="grey")
    if sep is None:
        sep = ["\t","\t"]
    if get_lead_args is None:
        get_lead_args = {}
    if anno=="GENENAME":
        get_lead_args["anno"]=True
    if anno_args is None:
        anno_args = {}
    if errargs is None:
        errargs={"ecolor":"#cccccc","elinewidth":1}
    if fontargs is None:
        fontargs={'fontsize':12,'family':'sans','fontname':'Arial'}
    if helper_line_args is None:
        helper_line_args={"color":'black', "linestyle":'-',"lw":1}
    if plt_args is None:
        plt_args={"figsize":(8,8),"dpi":300}
    if scatterargs is None:
        scatterargs={"s":20}
    if label is None:
        label = ["Sumstats_1","Sumstats_2","Both","None"]
    if anno_het ==True:
        is_q=True
    
    adjust_text_kwargs_r_default = {"autoalign":False,"precision":0.001,"lim":1000,"ha":"left","va":"top","expand_text":(1,1.8),"expand_objects":(0.1,0.1),"expand_points":(1.8,1.8),"force_objects":(0.8,0.8),"arrowprops":dict(arrowstyle='-|>', color='grey')}
    adjust_text_kwargs_l_default = {"autoalign":False,"precision":0.001,"lim":1000,"ha":"right","va":"bottom","expand_text":(1,1.8),"expand_objects":(0.1,0.1),"expand_points":(1.8,1.8),"force_objects":(0.8,0.8),"arrowprops":dict(arrowstyle='-|>', color='grey')}

    if adjust_text_kwargs_l is None:
        adjust_text_kwargs_l = adjust_text_kwargs_l_default
    else:
        for key, value in adjust_text_kwargs_l_default.items():
            if key not in adjust_text_kwargs_l:
                adjust_text_kwargs_l[key] = value

    if adjust_text_kwargs_r is None:
        adjust_text_kwargs_r = adjust_text_kwargs_r_default
    else:
        for key, value in adjust_text_kwargs_r_default.items():
            if key not in adjust_text_kwargs_r:
                adjust_text_kwargs_r[key] = value
    
    if adjust_text_kwargs is not None:
        for key, value in adjust_text_kwargs.items():
            adjust_text_kwargs_l[key] = value
            adjust_text_kwargs_r[key] = value
    else:
        adjust_text_kwargs = {}


    save_kwargs =      _extract_kwargs("save", save_args, locals())
    anno_kwargs =      _extract_kwargs("anno", anno_args, locals())
    err_kwargs =       _extract_kwargs("err", errargs, locals())
    plt_kwargs =       _extract_kwargs("plt", plt_args,  locals())
    scatter_kwargs =   _extract_kwargs("scatter", scatterargs, locals())
    font_kwargs =      _extract_kwargs("font",fontargs, locals())

    log.write("Start to process the raw sumstats for plotting...", verbose=verbose)
    
    # configure headers
    cols_name_list_1,cols_name_list_2, effect_cols_list_1, effect_cols_list_2 = configure_headers(mode, 
                                                                                                  path1, 
                                                                                                  path2, 
                                                                                                  cols_name_list_1,
                                                                                                  cols_name_list_2,
                                                                                                  effect_cols_list_1,
                                                                                                  effect_cols_list_2,
                                                                                                  scaled1,
                                                                                                  scaled2,
                                                                                                  log, 
                                                                                                  verbose)
    
    # extract common variants / load sumstats 1
    sumstats, common_snp_set = configure_common_snp_set(path1,path2,
                               snplist,
                               label,
                               cols_name_list_1,
                               cols_name_list_2,
                               sep,
                               scaled1,
                               scaled2,
                               log,verbose)
    
    # rename sumstats headers -> keywords in gwaslab
    sumstats = rename_sumtats(sumstats=sumstats,
                              cols_name_list = cols_name_list_1, 
                              scaled=scaled1,
                              snplist=snplist)
    
    # exctract only available variants from sumstats1 
    sumstats = sumstats.loc[sumstats["SNPID"].isin(common_snp_set),:]
    log.write(" -Using only variants available for both datasets...", verbose=verbose)
    
    ######### 8 extact SNPs for comparison 
    sig_list_1 = extract_snp_for_comparison(sumstats, 
                                            snplist, 
                                            label=label[0], 
                                            get_lead_args=get_lead_args, 
                                            build=build,
                                            drop=drop, 
                                            anno=anno, 
                                            sig_level=sig_level, 
                                            scaled = scaled1,
                                            log = log, 
                                            verbose = verbose)

 
    ######### load sumstats1

    ######### 9 extract snplist2
    if snplist is not None:
        cols_to_extract = [cols_name_list_2[0],cols_name_list_2[1]]
    else:
        cols_to_extract = [cols_name_list_2[0],cols_name_list_2[1],cols_name_list_2[4],cols_name_list_2[5]]
    
    sumstats = load_sumstats(path=path2, 
                            usecols=cols_to_extract, 
                            label=label[1], 
                            log=log, 
                            verbose= verbose, 
                            sep=sep[1])
    gc.collect()
    
    #if scaled2==True:
    #    sumstats[cols_name_list_2[1]] = np.power(10,-sumstats[cols_name_list_2[1]])

    sumstats = rename_sumtats(sumstats=sumstats,
                              cols_name_list = cols_name_list_2, 
                              scaled=scaled2,
                              snplist=snplist)
    ######### 11 exctract only overlapping variants from sumstats2
    sumstats = sumstats.loc[sumstats["SNPID"].isin(common_snp_set),:]
    sig_list_2 = extract_snp_for_comparison(sumstats, 
                                        snplist, 
                                        label=label[1], 
                                        get_lead_args=get_lead_args, 
                                        build=build,
                                        drop=drop, 
                                        anno=anno, 
                                        sig_level=sig_level, 
                                        scaled = scaled2,
                                        log = log, 
                                        verbose = verbose)

    ######### 13 Merge two list using SNPID
    sig_list_merged = merge_list(sig_list_1, 
                                 sig_list_2, 
                                 anno = anno, 
                                 labels=label,
                                 log=log, 
                                 verbose=verbose)

    ###############################################################################
    cols_to_extract = configure_cols_to_extract(mode=mode, 
                                                cols_name_list = cols_name_list_1, 
                                                effect_cols_list= effect_cols_list_1, 
                                                eaf = eaf)
    sumstats = load_sumstats(path=path1, 
                            usecols=cols_to_extract, 
                            label=label[0], 
                            log=log, 
                            verbose= verbose, 
                            sep=sep[0])
    
    #if scaled1==True:
    #    sumstats[cols_name_list_1[1]] = np.power(10,-sumstats[cols_name_list_1[1]])
    sumstats = rename_sumstats_full(mode, sumstats,
                                    index=1, 
                                    cols_name_list = cols_name_list_1, 
                                    effect_cols_list = effect_cols_list_1, 
                                    eaf = eaf, 
                                    drop = drop, 
                                    scaled=scaled1,
                                    log=log, verbose=verbose)

    log.write(" -Merging "+label[0]+" effect information...", verbose=verbose)
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")

    ############ 15 merging sumstats2
    cols_to_extract = configure_cols_to_extract(mode=mode, 
                                                cols_name_list = cols_name_list_2, 
                                                effect_cols_list= effect_cols_list_2, 
                                                eaf = eaf)
    
    sumstats = load_sumstats(path=path2, 
                            usecols=cols_to_extract, 
                            label=label[1], 
                            log=log, 
                            verbose= verbose, 
                            sep=sep[1])
    
    #if scaled2==True:
    #    sumstats[cols_name_list_2[1]] = np.power(10,-sumstats[cols_name_list_2[1]])
    
    gc.collect()
    
    sumstats = rename_sumstats_full(mode, sumstats,
                                    index=2, 
                                    cols_name_list = cols_name_list_2, 
                                    effect_cols_list = effect_cols_list_2, 
                                    eaf = eaf, 
                                    drop = drop, 
                                    scaled=scaled2,
                                    log=log, verbose=verbose)

    log.write(" -Merging "+label[1]+" effect information...", verbose=verbose)
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")
    
    sig_list_merged.set_index("SNPID",inplace=True)

    ################ 16 update sumstats1

    sig_list_merged = update_stats(sig_list_merged = sig_list_merged, 
                                   path = path1, 
                                   cols_name_list = cols_name_list_1,
                                   index=1, 
                                   sep=sep[0], 
                                   snplist = snplist, 
                                   label=label[0], 
                                   drop = drop, 
                                   scaled=scaled1,
                                   log=log, 
                                   verbose = verbose)
    
    ################# 17 update sumstats2
    sig_list_merged = update_stats(sig_list_merged = sig_list_merged, 
                                   path = path2, 
                                   cols_name_list = cols_name_list_2,
                                   index=2, 
                                   sep=sep[1], 
                                   snplist = snplist, 
                                   label=label[1], 
                                   drop = drop, 
                                   scaled=scaled2,
                                   log=log, 
                                   verbose = verbose)

    #if scaled1 ==True :
    #    log.write(" -Sumstats -log10(P) values are being converted to P...", verbose=verbose)
    #    sig_list_merged["P_1"] = np.power(10,-sig_list_merged["P_1"])
    #if scaled2 ==True :
    #    log.write(" -Sumstats -log10(P) values are being converted to P...", verbose=verbose)
    #    sig_list_merged["P_2"] = np.power(10,-sig_list_merged["P_2"])
    
#################################################################################
    sig_list_merged = assign_indicator(sig_list_merged, snplist, sig_level, scaled1, scaled2, log, verbose)
    
    sig_list_merged = align_alleles(sig_list_merged, label, mode, eaf, log, verbose)

    sig_list_merged = check_allele_match(sig_list_merged, allele_match, label, log,verbose)

    sig_list_merged = filter_by_maf(sig_list_merged, eaf, maf_level, log, verbose)

    if fdr==True and scaled==False:
        log.write(" -Using FDR...", verbose=verbose)
        #sig_list_merged["P_1"] = fdrcorrection(sig_list_merged["P_1"])[1]
        #sig_list_merged["P_2"] = fdrcorrection(sig_list_merged["P_2"])[1]
        sig_list_merged["P_1"] =ss.false_discovery_control(sig_list_merged["P_1"])
        sig_list_merged["P_2"] =ss.false_discovery_control(sig_list_merged["P_2"])

    ####################################################################################################################################
    ## winner's curse correction using aligned beta
    sig_list_merged = winnerscurse_correction(sig_list_merged, mode, wc_correction, sig_level,scaled1, scaled2, log, verbose)

    ########################## Het test############################################################
    ## heterogeneity test
    if (is_q == True):
        log.write(" -Calculating Cochran's Q statistics and peform chisq test...", verbose=verbose)
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            sig_list_merged = test_q(sig_list_merged,"EFFECT_1","SE_1","EFFECT_2_aligned","SE_2",q_level=q_level,is_q_mc=is_q_mc, log=log, verbose=verbose)
        else:
            sig_list_merged = test_q(sig_list_merged,"BETA_1","SE_1","BETA_2_aligned","SE_2",q_level=q_level,is_q_mc=is_q_mc, log=log, verbose=verbose)

        # heterogeneity summary
        log.write(" -Significant het:" ,len(sig_list_merged.loc[sig_list_merged["HetP"]<0.05,:]), verbose=verbose)
        log.write(" -All sig:" ,len(sig_list_merged), verbose=verbose)
        log.write(" -Het rate:" ,len(sig_list_merged.loc[sig_list_merged["HetP"]<0.05,:])/len(sig_list_merged), verbose=verbose)   
    
    ######################### save ###############################################################
    ## save the merged data
    save_path = label[0]+"_"+label[1]+"_beta_sig_list_merged.tsv"
    log.write(" -Saving the merged data to:",save_path, verbose=verbose)
    sig_list_merged = reorder_columns(sig_list_merged)
    sig_list_merged.to_csv(save_path,sep="\t")
    
    # extract group
    if include_all==True:
        sum0 = sig_list_merged.loc[sig_list_merged["indicator"]==0,:].dropna(axis=0)
    else:
        sum0 = pd.DataFrame()

    sum1only = sig_list_merged.loc[sig_list_merged["indicator"]==1,:].copy()
    sum2only = sig_list_merged.loc[sig_list_merged["indicator"]==2,:].copy()
    both     = sig_list_merged.loc[sig_list_merged["indicator"]==3,:].copy()
    
    if is_q==False:
        sum0["Edge_color"]="none"
        sum1only["Edge_color"]="none"
        sum2only["Edge_color"]="none"
        both["Edge_color"]="none"

    log.write(" -Identified "+str(len(sum0)) + " variants which are not significant in " + label[3]+".", verbose=verbose)
    log.write(" -Identified "+str(len(sum1only)) + " variants which are only significant in " + label[0]+".", verbose=verbose)
    log.write(" -Identified "+str(len(sum2only)) + " variants which are only significant in " + label[1]+".", verbose=verbose)
    log.write(" -Identified "+str(len(both)) + " variants which are significant in " + label[2] + ".", verbose=verbose)
    
    ##plot########################################################################################
    log.write("Creating the scatter plot for effect sizes comparison...", verbose=verbose)
    #plt.style.use("ggplot")
    sns.set_style("ticks")
    fig,ax = plt.subplots(**plt_kwargs) 
    legend_elements=[]
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        if len(sum0)>0:
            ax.errorbar(sum0["EFFECT_1"],sum0["EFFECT_2_aligned"], xerr=sum0["SE_1"],yerr=sum0["SE_2"],
                        linewidth=0,zorder=1,**err_kwargs)
            
            ax.scatter(sum0["EFFECT_1"],sum0["EFFECT_2_aligned"],label=label[3],zorder=2,color="#cccccc",edgecolors=sum0["Edge_color"],marker=".",**scatter_kwargs)
            #legend_elements.append(mpatches.Circle(facecolor='#cccccc', edgecolor='white', label=label[3]))
            legend_elements.append(label[3])
        if len(sum1only)>0:
            ax.errorbar(sum1only["EFFECT_1"],sum1only["EFFECT_2_aligned"], xerr=sum1only["SE_1"],yerr=sum1only["SE_2"],
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(sum1only["EFFECT_1"],sum1only["EFFECT_2_aligned"],label=label[0],zorder=2,color="#e6320e",edgecolors=sum1only["Edge_color"],marker="^",**scatter_kwargs)
            #legend_elements.append(mpatches.Patch(facecolor='#e6320e', edgecolor='white', label=label[0]))
            legend_elements.append(label[0])
        if len(sum2only)>0:
            ax.errorbar(sum2only["EFFECT_1"],sum2only["EFFECT_2_aligned"], xerr=sum2only["SE_1"],yerr=sum2only["SE_2"],
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(sum2only["EFFECT_1"],sum2only["EFFECT_2_aligned"],label=label[1],zorder=2,color="#41e620",edgecolors=sum2only["Edge_color"],marker="o",**scatter_kwargs)
            #legend_elements.append(mpatches.Circle(facecolor='#41e620', edgecolor='white', label=label[1]))
            legend_elements.append(label[1])
        if len(both)>0:
            ax.errorbar(both["EFFECT_1"],both["EFFECT_2_aligned"], xerr=both["SE_1"],yerr=both["SE_2"],
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(both["EFFECT_1"],both["EFFECT_2_aligned"],label=label[2],zorder=2,color="#205be6",edgecolors=both["Edge_color"],marker="s",**scatter_kwargs)  
            #legend_elements.append(mpatches.Patch(facecolor='#205be6', edgecolor='white', label=label[2]))
            legend_elements.append(label[2])
    else:
        ## if OR
        if len(sum0)>0:
            ax.errorbar(sum0["OR_1"],sum0["OR_2_aligned"], xerr=sum0[["OR_L_1_err","OR_H_1_err"]].T,yerr=sum0[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T,
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(sum0["OR_1"],sum0["OR_2_aligned"],label=label[3],zorder=2,color="#cccccc",edgecolors=sum0["Edge_color"],marker=".",**scatter_kwargs)
            legend_elements.append(label[3])
        if len(sum1only)>0:
            ax.errorbar(sum1only["OR_1"],sum1only["OR_2_aligned"], xerr=sum1only[["OR_L_1_err","OR_H_1_err"]].T,yerr=sum1only[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T,
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(sum1only["OR_1"],sum1only["OR_2_aligned"],label=label[0],zorder=2,color="#e6320e",edgecolors=sum1only["Edge_color"],marker="^",**scatter_kwargs)
            legend_elements.append(label[0])
        if len(sum2only)>0:
            ax.errorbar(sum2only["OR_1"],sum2only["OR_2_aligned"], xerr=sum2only[["OR_L_1_err","OR_H_1_err"]].T,yerr=sum2only[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T,
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(sum2only["OR_1"],sum2only["OR_2_aligned"],label=label[1],zorder=2,color="#41e620",edgecolors=sum2only["Edge_color"],marker="o",**scatter_kwargs)
            legend_elements.append(label[1])
        if len(both)>0:
            ax.errorbar(both["OR_1"],both["OR_2_aligned"], xerr=both[["OR_L_1_err","OR_H_1_err"]].T,yerr=both[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T,
                        linewidth=0,zorder=1,**err_kwargs)
            ax.scatter(both["OR_1"],both["OR_2_aligned"],label=label[2],zorder=2,color="#205be6",edgecolors=both["Edge_color"],marker="s",**scatter_kwargs)
            legend_elements.append(label[2])
    ## annotation #################################################################################################################
    ax = scatter_annotation(ax, sig_list_merged,anno, anno_het, is_q, mode, 
                       anno_min,anno_min1,anno_min2,anno_diff,anno_kwargs,adjust_text_kwargs_l,adjust_text_kwargs_r,
                       log,verbose
                       )
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
    ax = confire_regression_line(is_reg,reg_box, sig_list_merged, ax, mode,xl,yl,xh,yh, null_beta, r_se, 
                            is_45_helper_line,helper_line_args, font_kwargs,
                            log, verbose)
        
    
    ax.set_xlabel(xylabel_prefix+label[0],**font_kwargs)
    ax.set_ylabel(xylabel_prefix+label[1],**font_kwargs)
    
    ax = configure_legend(fig, ax, legend_mode, is_q, is_q_mc, legend_elements, legend_pos, q_level, 
                        font_kwargs,scatterargs,legend_args,
                        legend_title, legend_title2 )
    ##plot finished########################################################################################
    gc.collect()

    save_figure(fig, save, keyword="esc",save_args=save_kwargs, log=log, verbose=verbose)
    
    sig_list_merged = reorder_columns(sig_list_merged)
    
    return [sig_list_merged, fig,log]

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

def load_sumstats(path, usecols, label, log, verbose, sep):
    if type(usecols) is not list:
        usecols = [usecols]
    
    log.write(" -Loading sumstats for {} : {}".format(label,",".join(usecols)), verbose=verbose)
    #log.write(" -Loading {} SNP list in memory...".format(label), verbose=verbose)   
    
    if type(path) is Sumstats:
        sumstats = path.data.loc[:,usecols].copy()
    elif type(path) is pd.DataFrame:
        sumstats = path.loc[:,usecols].copy()
    else:
        sumstats=pd.read_table(path,sep=sep,usecols=usecols)
    return sumstats

def configure_headers(mode, 
                    path1, 
                    path2, 
                    cols_name_list_1,
                    cols_name_list_2,
                    effect_cols_list_1,
                    effect_cols_list_2,
                    scaled1,
                    scaled2,
                    log, 
                    verbose):
    
    if mode not in ["Beta","beta","BETA","OR","or"]:
        raise ValueError("Please input Beta or OR")
    
    if type(path1) is Sumstats:
        log.write("Path1 is gwaslab Sumstats object...", verbose=verbose)
        if cols_name_list_1 is None:
            cols_name_list_1 = ["SNPID","P","EA","NEA","CHR","POS"]
            if scaled1==True:
                cols_name_list_1 = ["SNPID","MLOG10P","EA","NEA","CHR","POS"]
        if effect_cols_list_1 is None:
            if mode=="beta":
                effect_cols_list_1 = ["BETA","SE"]
            else:
                effect_cols_list_1 = ["OR","OR_95L","OR_95U"]
    elif type(path1) is pd.DataFrame:
        log.write("Path1 is pandas DataFrame object...", verbose=verbose)

    if type(path2) is Sumstats:
        log.write("Path2 is gwaslab Sumstats object...", verbose=verbose)
        if cols_name_list_2 is None:
            cols_name_list_2 = ["SNPID","P","EA","NEA","CHR","POS"]
            if scaled2==True:
                cols_name_list_2 = ["SNPID","MLOG10P","EA","NEA","CHR","POS"]
        if effect_cols_list_2 is None:
            if mode=="beta":
                effect_cols_list_2 = ["BETA","SE"]
            else:
                effect_cols_list_2 = ["OR","OR_95L","OR_95U"]
    elif type(path2) is pd.DataFrame:
        log.write("Path2 is pandas DataFrame object...", verbose=verbose)

    return cols_name_list_1,cols_name_list_2, effect_cols_list_1, effect_cols_list_2

def configure_common_snp_set(path1,path2,
                             snplist,
                             label,
                             cols_name_list_1,cols_name_list_2,
                             sep,
                             scaled1,
                             scaled2,
                             log,verbose):
    
    ######### load sumstats2
    sumstats = load_sumstats(path=path2, 
                             usecols=cols_name_list_2[0], 
                             label=label[1], 
                             log=log, 
                             verbose= verbose, 
                             sep=sep[1])
    
    common_snp_set=set(sumstats[cols_name_list_2[0]].values)
    
    ######### extract snplist1
    if snplist is not None:
        #use only SNPID, P
        cols_to_extract = [cols_name_list_1[0],cols_name_list_1[1]]
    else:
        # use SNPID, P, chr pos
        cols_to_extract = [cols_name_list_1[0],cols_name_list_1[1],cols_name_list_1[4],cols_name_list_1[5]]
 
    ######### load sumstats1
    sumstats = load_sumstats(path=path1, 
                             usecols=cols_to_extract, 
                             label=label[0], 
                             log=log, 
                             verbose= verbose, 
                             sep=sep[0])
    
    gc.collect()

    #if scaled1==True:
    #    sumstats[cols_name_list_1[1]] = np.power(10,-sumstats[cols_name_list_1[1]])
    ######### 5 extract the common set
    
    common_snp_set = common_snp_set.intersection(sumstats[cols_name_list_1[0]].values)
    
    log.write(" -Counting  variants available for both datasets:",len(common_snp_set)," variants...", verbose=verbose)
    
    return sumstats, common_snp_set

def rename_sumtats(sumstats, cols_name_list, snplist, scaled,suffix=""):
    ######### 6 rename the sumstats
    rename_dict = { cols_name_list[0]:"SNPID",
                   cols_name_list[1]:"P{}".format(suffix),
               }
    if scaled==True:
        rename_dict[cols_name_list[1]] = "MLOG10P{}".format(suffix)
    
    if snplist is None: 
        rename_dict[cols_name_list[4]]="CHR"
        rename_dict[cols_name_list[5]]="POS"
    
    sumstats = sumstats.rename(columns=rename_dict)
    return sumstats


def extract_snp_for_comparison(sumstats, snplist, label, 
                               get_lead_args, build, drop, anno, 
                               sig_level,scaled, log, verbose):
    ######### 8 extact SNPs for comparison 
    if snplist is not None: 
        ######### 8.1 if a snplist is provided, use the snp list
        log.write(" -Extract variants in the given list from "+label+"...")
        sig_list = sumstats.loc[sumstats["SNPID"].isin(snplist),:].copy()
        if anno=="GENENAME":
            sig_list = annogene(sig_list,"SNPID","CHR","POS", build=build, verbose=verbose, **get_lead_args)
    else:
        ######### 8,2 otherwise use the automatically detected lead SNPs
        log.write(" -Extract lead variants from "+label +"...", verbose=verbose)
        sig_list = getsig(sumstats,"SNPID","CHR","POS","P","MLOG10P", build=build, verbose=verbose,sig_level=sig_level,**get_lead_args)
    
    if drop==True:
        if scaled==True:
            sig_list = drop_duplicate_and_na(sig_list,  sort_by="MLOG10P",ascending=False, log=log , verbose=verbose)
        else:
            sig_list = drop_duplicate_and_na(sig_list,  sort_by="P", ascending=True, log=log , verbose=verbose)

    return sig_list

def merge_list(sig_list_1, sig_list_2, anno,labels,log, verbose):
    
    log.write("Merging snps from "+labels[0]+" and "+labels[1]+"...", verbose=verbose)
    
    if anno == "GENENAME":
        if "GENE" not in sig_list_1.columns:
            sig_list_1["GENE"]=pd.NA
            sig_list_1["LOCATION"]=pd.NA
        if "GENE" not in sig_list_2.columns:
            sig_list_2["GENE"]=pd.NA
            sig_list_2["LOCATION"]=pd.NA

    sig_list_merged = pd.merge(sig_list_1,sig_list_2,left_on="SNPID",right_on="SNPID",how="outer",suffixes=('_1', '_2'))
    
    if anno == "GENENAME":
        sig_list_merged.loc[sig_list_merged["SNPID"].isin((sig_list_1["SNPID"])),"GENENAME"] = sig_list_merged.loc[sig_list_merged["SNPID"].isin((sig_list_1["SNPID"])),"GENE_1"]
        sig_list_merged.loc[~sig_list_merged["SNPID"].isin((sig_list_1["SNPID"])),"GENENAME"] = sig_list_merged.loc[~sig_list_merged["SNPID"].isin((sig_list_1["SNPID"])),"GENE_2"]
        sig_list_merged = sig_list_merged.drop(columns=["GENE_1","GENE_2","LOCATION_1","LOCATION_2"])
    #     SNPID       P_1       P_2
    #0   rs117986209  0.142569  0.394455
    #1     rs6704312  0.652104  0.143750
    return sig_list_merged

def configure_cols_to_extract(mode, 
                              cols_name_list, 
                              effect_cols_list, 
                              eaf):
    
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        #[snpid,p,ea,nea]      ,[effect,se]
        #[snpid,p,ea,nea,chr,pos],[effect,se]
        #[snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]
        cols_to_extract = [cols_name_list[0],cols_name_list[1], cols_name_list[2],cols_name_list[3], effect_cols_list[0], effect_cols_list[1]]
    else:
        cols_to_extract = [cols_name_list[0],cols_name_list[1], cols_name_list[2],cols_name_list[3], effect_cols_list[0], effect_cols_list[1], effect_cols_list[2]]
    
    if len(eaf)>0: 
        cols_to_extract.append(eaf[0])   
    
    return cols_to_extract

def rename_sumstats_full(mode, sumstats, cols_name_list, effect_cols_list, eaf, drop, index, scaled, log, verbose):
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        rename_dict = { cols_name_list[0]:"SNPID",
                        cols_name_list[1]:"P_{}".format(index),
                        cols_name_list[2]:"EA_{}".format(index),
                        cols_name_list[3]:"NEA_{}".format(index),
                        effect_cols_list[0]:"EFFECT_{}".format(index),
                        effect_cols_list[1]:"SE_{}".format(index)}

        
    else:
        # if or
        rename_dict = { cols_name_list[0]:"SNPID",
                        cols_name_list[1]:"P_{}".format(index),
                        cols_name_list[2]:"EA_{}".format(index),
                        cols_name_list[3]:"NEA_{}".format(index),
                        effect_cols_list[0]:"OR_{}".format(index),
                        effect_cols_list[1]:"OR_L_{}".format(index),
                        effect_cols_list[2]:"OR_H_{}".format(index)}
    if scaled==True:
            rename_dict[cols_name_list[1]]="MLOG10P_{}".format(index)
    ## check if eaf column is provided.
    if len(eaf)>0: 
        rename_dict[eaf[index-1]]="EAF_{}".format(index)
    sumstats = sumstats.rename(columns=rename_dict)
    
    # drop na and duplicate
    if drop==True:
        if scaled==True:
            sumstats = drop_duplicate_and_na(sumstats,  sort_by="MLOG10P_{}".format(index),ascending=False, log=log , verbose=verbose)
        else:
            sumstats = drop_duplicate_and_na(sumstats,  sort_by="P_{}".format(index), ascending=True, log=log , verbose=verbose)

    if scaled==True:
        sumstats.drop("MLOG10P_{}".format(index),axis=1,inplace=True)
    else:
        sumstats.drop("P_{}".format(index),axis=1,inplace=True)
    return sumstats
    
def update_stats(sig_list_merged, 
                 path, 
                 cols_name_list,
                 sep, 
                 snplist, 
                 label, 
                 drop, 
                 index, 
                 scaled,
                 log, 
                 verbose):
    
    log.write(" -Updating missing information for "+label+" ...", verbose=verbose)
    cols_to_extract = [cols_name_list[0], cols_name_list[1]]
    
    sumstats = load_sumstats(path=path, 
                             usecols=cols_to_extract, 
                             label=label, 
                             log=log, 
                             verbose= verbose, 
                             sep=sep)
    #if scaled1==True:
    #    sumstats[cols_name_list_1[1]] = np.power(10,-sumstats[cols_name_list_1[1]])
    
    sumstats = rename_sumtats(sumstats = sumstats, 
                              cols_name_list = cols_name_list, 
                              snplist = snplist,
                              scaled=scaled, 
                              suffix="_{}".format(index))
    # drop na and duplicate
    if drop==True:
        if scaled==True:
            sumstats = drop_duplicate_and_na(sumstats,  sort_by="MLOG10P_{}".format(index),ascending=False, log=log , verbose=verbose)
        else:
            sumstats = drop_duplicate_and_na(sumstats,  sort_by="P_{}".format(index), ascending=True, log=log , verbose=verbose)

    
    sumstats = sumstats.set_index("SNPID")
    sig_list_merged.update(sumstats)

    return sig_list_merged


def assign_indicator(sig_list_merged, snplist, sig_level, scaled1, scaled2, log, verbose):
    ############## 18 init indicator
    log.write(" -Assigning indicator  ...", verbose=verbose)
    # 0-> 0
    # 1 -> sig in sumstats1
    # 2 -> sig in sumsatts2
    # 3->  sig in both sumstats1 + sumstats2
    sig_list_merged["indicator"] = 0
    
    if scaled1==True:
        sig_list_merged.loc[sig_list_merged["MLOG10P_1"]>-np.log10(sig_level),"indicator"]=1+sig_list_merged.loc[sig_list_merged["MLOG10P_1"]>-np.log10(sig_level),"indicator"]
    else:
        sig_list_merged.loc[sig_list_merged["P_1"]<sig_level,"indicator"]=1+sig_list_merged.loc[sig_list_merged["P_1"]<sig_level,"indicator"]
    
    if scaled2==True:
        sig_list_merged.loc[sig_list_merged["MLOG10P_2"]>-np.log10(sig_level),"indicator"]=2+sig_list_merged.loc[sig_list_merged["MLOG10P_2"]>-np.log10(sig_level),"indicator"]
    else:
        sig_list_merged.loc[sig_list_merged["P_2"]<sig_level,"indicator"]=2+sig_list_merged.loc[sig_list_merged["P_2"]<sig_level,"indicator"]
    
    if snplist is None:
        sig_list_merged["CHR"]=np.max(sig_list_merged[["CHR_1","CHR_2"]], axis=1).astype(int)
        sig_list_merged["POS"]=np.max(sig_list_merged[["POS_1","POS_2"]], axis=1).astype(int)
        sig_list_merged.drop(labels=['CHR_1', 'CHR_2','POS_1', 'POS_2'], axis=1,inplace=True)
    return sig_list_merged

def align_alleles(sig_list_merged, label,mode,eaf, log, verbose):
    log.write(" -Aligning "+label[1]+" EA with "+label[0]+" EA ...", verbose=verbose)
    ############### 19 align allele effect with sumstats 1
    sig_list_merged["EA_1"]=sig_list_merged["EA_1"].astype("string")
    sig_list_merged["EA_2"]=sig_list_merged["EA_2"].astype("string")
    sig_list_merged["NEA_1"]=sig_list_merged["NEA_1"].astype("string")
    sig_list_merged["NEA_2"]=sig_list_merged["NEA_2"].astype("string")
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

        sig_list_merged["EA_2_aligned"]=sig_list_merged["EA_2"]
        sig_list_merged["NEA_2_aligned"]=sig_list_merged["NEA_2"]
        sig_list_merged["OR_2_aligned"]=sig_list_merged["OR_2"]
        sig_list_merged["OR_L_2_aligned"]=sig_list_merged["OR_L_2"]
        sig_list_merged["OR_H_2_aligned"]=sig_list_merged["OR_H_2"]

        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"NEA_2_aligned"]= sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EA_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_2_aligned"]= 1/sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_H_2_aligned"]= 1/sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_L_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_L_2_aligned"]= 1/sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"OR_H_2"]
        
        sig_list_merged["BETA_1"]=np.log(sig_list_merged["OR_1"])
        sig_list_merged["BETA_2_aligned"]=np.log(sig_list_merged["OR_2_aligned"])
        sig_list_merged["SE_1"]=(np.log(sig_list_merged["OR_H_1"]) - np.log(sig_list_merged["OR_1"]))/ss.norm.ppf(0.975)
        sig_list_merged["SE_2"]=(np.log(sig_list_merged["OR_H_2_aligned"]) - np.log(sig_list_merged["OR_2_aligned"]))/ss.norm.ppf(0.975)
        
        sig_list_merged["OR_L_1_err"]=np.abs(sig_list_merged["OR_L_1"]-sig_list_merged["OR_1"])
        sig_list_merged["OR_H_1_err"]=np.abs(sig_list_merged["OR_H_1"]-sig_list_merged["OR_1"])
        sig_list_merged["OR_L_2_aligned_err"]=np.abs(sig_list_merged["OR_L_2_aligned"]-sig_list_merged["OR_2_aligned"])
        sig_list_merged["OR_H_2_aligned_err"]=np.abs(sig_list_merged["OR_H_2_aligned"]-sig_list_merged["OR_2_aligned"])
        
    if len(eaf)>0:
        # flip eaf
        sig_list_merged["EAF_2_aligned"]=sig_list_merged["EAF_2"]
        sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EAF_2_aligned"]= 1 -sig_list_merged.loc[sig_list_merged["EA_1"]!=sig_list_merged["EA_2"],"EAF_2"]
    return sig_list_merged

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

def check_allele_match(sig_list_merged, allele_match, label, log,verbose):
    # checking effect allele matching
    nonmatch = np.nansum(sig_list_merged["EA_1"] != sig_list_merged["EA_2_aligned"])
    log.write(" -Aligned all EAs in {} with EAs in {} ...".format(label[1],label[0]), verbose=verbose)
    if nonmatch>0:
        log.warning("Alleles for {} variants do not match...".format(nonmatch))
    if allele_match==True:
        if nonmatch>0:
            sig_list_merged = sig_list_merged.loc[sig_list_merged["EA_1"] == sig_list_merged["EA_2_aligned"]]
        else:
            log.write(" -No variants with EA not matching...", verbose=verbose)
    return sig_list_merged

def winnerscurse_correction(sig_list_merged, mode, wc_correction, sig_level, scaled1, scaled2, log, verbose):
    if mode=="beta":
        if scaled1==True:            
            match1=  sig_list_merged["MLOG10P_1"]>-np.log10(sig_level)
        else:
            match1 = sig_list_merged["P_1"]<sig_level
        if scaled2==True:            
            match2=  sig_list_merged["MLOG10P_2"]>-np.log10(sig_level)
        else:
            match2 = sig_list_merged["P_2"]<sig_level

        if wc_correction == "all":
            log.write(" -Correcting BETA for winner's curse with threshold at {} for all variants...".format(sig_level), verbose=verbose)
            sig_list_merged["EFFECT_1_RAW"] = sig_list_merged["EFFECT_1"].copy()
            sig_list_merged["EFFECT_2_aligned_RAW"] = sig_list_merged["EFFECT_2_aligned"].copy()
            
            log.write("  -Correcting BETA for {} variants in sumstats1...".format(sum(~sig_list_merged["EFFECT_1"].isna())), verbose=verbose)
            sig_list_merged["EFFECT_1"] = sig_list_merged[["EFFECT_1_RAW","SE_1"]].apply(lambda x: wc_correct(x[0],x[1],sig_level),axis=1)

            log.write("  -Correcting BETA for {} variants in sumstats2...".format(sum(~sig_list_merged["EFFECT_2_aligned"].isna())), verbose=verbose)
            sig_list_merged["EFFECT_2_aligned"] = sig_list_merged[["EFFECT_2_aligned_RAW","SE_2"]].apply(lambda x: wc_correct(x[0],x[1],sig_level),axis=1)
        
        elif wc_correction == "sig" :

            log.write(" - Correcting BETA for winner's curse with threshold at {} for significant variants...".format(sig_level), verbose=verbose)
            sig_list_merged["EFFECT_1_RAW"] = sig_list_merged["EFFECT_1"].copy()
            sig_list_merged["EFFECT_2_aligned_RAW"] = sig_list_merged["EFFECT_2_aligned"].copy()
            log.write("  -Correcting BETA for {} variants in sumstats1...".format(sum(match1)), verbose=verbose)
            sig_list_merged.loc[match1, "EFFECT_1"]         = sig_list_merged.loc[match1, ["EFFECT_1_RAW","SE_1"]].apply(lambda x: wc_correct_test(x[0],x[1],sig_level),axis=1)
            log.write("  -Correcting BETA for {} variants in sumstats2...".format(sum(match2)), verbose=verbose)
            sig_list_merged.loc[match2, "EFFECT_2_aligned"] = sig_list_merged.loc[match2, ["EFFECT_2_aligned_RAW","SE_2"]].apply(lambda x: wc_correct_test(x[0],x[1],sig_level),axis=1)
        
        elif wc_correction == "sumstats1" :
            log.write(" - Correcting BETA for winner's curse with threshold at {} for significant variants in sumstats1...".format(sig_level), verbose=verbose)
            sig_list_merged["EFFECT_1_RAW"] = sig_list_merged["EFFECT_1"].copy()
            log.write("  -Correcting BETA for {} variants in sumstats1...".format(sum(match1)), verbose=verbose)
            sig_list_merged.loc[match1, "EFFECT_1"]         = sig_list_merged.loc[match1, ["EFFECT_1_RAW","SE_1"]].apply(lambda x: wc_correct_test(x[0],x[1],sig_level),axis=1)
            
        elif wc_correction == "sumstats2" :
            log.write(" - Correcting BETA for winner's curse with threshold at {} for significant variants in sumstats2...".format(sig_level), verbose=verbose)
            sig_list_merged["EFFECT_2_aligned_RAW"] = sig_list_merged["EFFECT_2_aligned"].copy()
            log.write("  -Correcting BETA for {} variants in sumstats2...".format(sum(match2)), verbose=verbose)
            sig_list_merged.loc[match2, "EFFECT_2_aligned"] = sig_list_merged.loc[match2, ["EFFECT_2_aligned_RAW","SE_2"]].apply(lambda x: wc_correct_test(x[0],x[1],sig_level),axis=1)
    return sig_list_merged

def filter_by_maf(sig_list_merged, eaf, maf_level, log, verbose):
    if (len(eaf)>0) and (maf_level is not None):
        both_eaf_clear =  (sig_list_merged["EAF_1"]>maf_level)&(sig_list_merged["EAF_1"]<1-maf_level)&(sig_list_merged["EAF_2"]>maf_level)&(sig_list_merged["EAF_2"]<1-maf_level)
        log.write(" -Exclude "+str(len(sig_list_merged) -sum(both_eaf_clear))+ " variants with maf <",maf_level, verbose=verbose)
        sig_list_merged = sig_list_merged.loc[both_eaf_clear,:]
    return sig_list_merged

    



def test_q(df,beta1,se1,beta2,se2,q_level=0.05,is_q_mc=False, log=Log(), verbose=False):
    w1="Weight_1"
    w2="Weight_2"
    beta="BETA_FE"
    q="Q"
    pq="HetP"
    rawpq="RAW_HetP"
    i2="I2"
    df[w1]=1/(df[se1])**2
    df[w2]=1/(df[se2])**2
    df[beta] =(df[w1]*df[beta1] + df[w2]*df[beta2])/(df[w1]+df[w2])
    
    # Cochran(1954)
    df[q] = df[w1]*(df[beta1]-df[beta])**2 + df[w2]*(df[beta2]-df[beta])**2
    df[pq] = ss.chi2.sf(df[q], 1)
    df["Edge_color"]="white"
    
    if is_q_mc=="fdr":
        log.write(" -FDR correction applied...", verbose=verbose)
        df[rawpq] = df[pq] 
        df[pq] = ss.false_discovery_control(df[pq])
        
    elif is_q_mc=="bon":
        log.write(" -Bonferroni correction applied...", verbose=verbose)
        df[rawpq] = df[pq] 
        df[pq] = df[pq] * len(df[pq])

    df.loc[df[pq]<q_level,"Edge_color"]="black"
    df.drop(columns=["Weight_1","Weight_2","BETA_FE"],inplace=True)
    # Huedo-Medina, T. B., Sánchez-Meca, J., Marín-Martínez, F., & Botella, J. (2006). Assessing heterogeneity in meta-analysis: Q statistic or I² index?. Psychological methods, 11(2), 193.
    
    # calculate I2
    df[i2] = (df[q] - 1)/df[q]
    df.loc[df[i2]<0,i2] = 0 
    
    return df

def jackknife_r(df,x="EFFECT_1",y="EFFECT_2_aligned"):
    """Jackknife estimation of se for rsq

    """

    # dropna
    df_nona = df.loc[:,[x,y]].dropna()
    
    # non-empty entries
    n=len(df)
    
    # assign row number
    df_nona["nrow"] = range(n)
    
    # a list to store r2
    r_list=[]
    
    # estimate r
    for i in range(n):
        # exclude 1 record
        records_to_use = df_nona["nrow"]!=i
        # estimate r
        reg_jackknife = ss.linregress(df_nona.loc[records_to_use, x],df_nona.loc[records_to_use,y])
        # add r_i to list
        r_list.append(reg_jackknife[2])

    # convert list to array
    rs = np.array(r_list)
    # https://en.wikipedia.org/wiki/Jackknife_resampling
    r_se = np.sqrt( (n-1)/n * np.sum((rs - np.mean(rs))**2) )
    return r_se

def drop_duplicate_and_na(df,snpid="SNPID",sort_by=False,log=Log(),ascending=True,verbose=True):
    
    length_before = len(df)
    
    if sort_by!=False:
        df.sort_values(by = sort_by, ascending=ascending, inplace=True)
    
    df.dropna(axis="index",subset=[snpid],inplace=True)
    df.drop_duplicates(subset=[snpid], keep='first', inplace=True) 
    
    length_after= len(df)
    if length_before !=  length_after:
        log.write(" -Dropped {} duplicates or NAs...".format(length_before - length_after), verbose=verbose)
    return df



#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

def scatter_annotation(ax, sig_list_merged,anno, anno_het, is_q, mode, 
                       anno_min,anno_min1,anno_min2,anno_diff,anno_kwargs,adjust_text_kwargs_l,adjust_text_kwargs_r,
                       log,verbose
                       ):
    if anno==True or anno=="GENENAME":
        sig_list_toanno = sig_list_merged.dropna(axis=0)
        if is_q==True and anno_het == True:
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["Edge_color"]=="black",:]

        if mode=="beta":
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["EFFECT_1"].abs() >=anno_min1 ,:]
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["EFFECT_2_aligned"].abs() >=anno_min2 ,:]
            sig_list_toanno = sig_list_toanno.loc[(sig_list_toanno["EFFECT_1"].abs() >=anno_min) & (sig_list_toanno["EFFECT_2_aligned"].abs() >=anno_min) ,:]
            sig_list_toanno = sig_list_toanno.loc[np.abs(sig_list_toanno["EFFECT_1"] - sig_list_toanno["EFFECT_2_aligned"]) >=anno_diff,:]
        else:            
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["OR_1"].abs() >=anno_min1 ,:]
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["OR_2_aligned"].abs() >=anno_min2 ,:]
            sig_list_toanno = sig_list_toanno.loc[(sig_list_toanno["OR_1"].abs() >=anno_min) & (sig_list_toanno["OR_2_aligned"].abs() >=anno_min) ,:]
            sig_list_toanno = sig_list_toanno.loc[np.abs(sig_list_toanno["OR_1"] - sig_list_toanno["OR_2_aligned"]) >=anno_diff,:]

        texts_l=[]
        texts_r=[]
        
        if anno==True:
            log.write("Annotating variants using {}".format("SNPID"), verbose=verbose)
        elif anno=="GENENAME":
            log.write("Annotating variants using {}".format("GENENAME"), verbose=verbose)
        
        for index, row in sig_list_toanno.iterrows():
            #log.write("Annotating {}...".format(row), verbose=verbose)
            if anno==True:
                to_anno_text = index
            elif type(anno) is str:
                if not pd.isna(row[anno]):
                    to_anno_text = row[anno]
                else:
                    to_anno_text = index

            if mode=="beta" or mode=="BETA" or mode=="Beta":
                if row["EFFECT_1"] <  row["EFFECT_2_aligned"]:
                    texts_l.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],to_anno_text,ha="right",va="bottom", **anno_kwargs))
                else:
                    texts_r.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],to_anno_text,ha="left",va="top", **anno_kwargs))
            else:
                if row["OR_1"] <  row["OR_2_aligned"]:
                    texts_l.append(plt.text(row["OR_1"], row["OR_2_aligned"],to_anno_text, ha='right', va='bottom', **anno_kwargs)) 
                else:
                    texts_r.append(plt.text(row["OR_1"], row["OR_2_aligned"],to_anno_text, ha='left', va='top', **anno_kwargs)) 
        if len(texts_l)>0:
            adjust_text(texts_l,ax=ax,**adjust_text_kwargs_l)
        if len(texts_r)>0:
            adjust_text(texts_r,ax=ax,**adjust_text_kwargs_r)
    elif type(anno) is dict:
        sig_list_toanno = sig_list_merged.dropna(axis=0)
        # if input is a dict
        sig_list_toanno = sig_list_toanno.loc[sig_list_toanno.index.isin(list(anno.keys())),:]
        if is_q==True and anno_het == True:
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["Edge_color"]=="black",:]

        if mode=="beta":
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["EFFECT_1"].abs() >=anno_min1 ,:]
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["EFFECT_2_aligned"].abs() >=anno_min2 ,:]
            sig_list_toanno = sig_list_toanno.loc[(sig_list_toanno["EFFECT_1"].abs() >=anno_min) & (sig_list_toanno["EFFECT_2_aligned"].abs() >=anno_min) ,:]
            sig_list_toanno = sig_list_toanno.loc[np.abs(sig_list_toanno["EFFECT_1"] - sig_list_toanno["EFFECT_2_aligned"]) >=anno_diff,:]
        else:            
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["OR_1"].abs() >=anno_min1 ,:]
            sig_list_toanno = sig_list_toanno.loc[sig_list_toanno["OR_2_aligned"].abs() >=anno_min2 ,:]
            sig_list_toanno = sig_list_toanno.loc[(sig_list_toanno["OR_1"].abs() >=anno_min) & (sig_list_toanno["OR_2_aligned"].abs() >=anno_min) ,:]
            sig_list_toanno = sig_list_toanno.loc[np.abs(sig_list_toanno["OR_1"] - sig_list_toanno["OR_2_aligned"]) >=anno_diff,:]
        
        texts_l=[]
        texts_r=[]
        for index, row in sig_list_toanno.iterrows():
            if mode=="beta" or mode=="BETA" or mode=="Beta":
                if row["EFFECT_1"] <  row["EFFECT_2_aligned"]:
                    texts_l.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],anno[index],ha="right",va="bottom", **anno_kwargs))
                else:
                    texts_r.append(plt.text(row["EFFECT_1"], row["EFFECT_2_aligned"],anno[index],ha="left",va="top", **anno_kwargs))
            else:
                if row["OR_1"] <  row["OR_2_aligned"]:
                    texts_l.append(plt.text(row["OR_1"], row["OR_2_aligned"],anno[index], ha='right', va='bottom', **anno_kwargs)) 
                else:
                    texts_r.append(plt.text(row["OR_1"], row["OR_2_aligned"],anno[index], ha='left', va='top', **anno_kwargs)) 
        if len(texts_l)>0:
            adjust_text(texts_l,ax=ax,**adjust_text_kwargs_l)
        if len(texts_r)>0:
            adjust_text(texts_r,ax=ax,**adjust_text_kwargs_r)
    return ax


def confire_regression_line(is_reg, reg_box, sig_list_merged,  ax, mode,xl,yl,xh,yh, null_beta, r_se, 
                            is_45_helper_line,helper_line_args, font_kwargs,
                            log, verbose):
    if len(sig_list_merged)<3: is_reg=False
    if is_reg is True:
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            reg = ss.linregress(sig_list_merged["EFFECT_1"],sig_list_merged["EFFECT_2_aligned"])
            
            # estimate se for r
            if r_se==True:
                log.write(" -Estimating SE for rsq using Jackknife method.", verbose=verbose)
                r_se_jackknife = jackknife_r(sig_list_merged)
                r_se_jackknife_string = " ({:.2f})".format(r_se_jackknife)
            else:
                r_se_jackknife_string= ""
        else:
            reg = ss.linregress(sig_list_merged["OR_1"],sig_list_merged["OR_2_aligned"])
            r_se_jackknife_string= ""

        #### calculate p values based on selected value , default = 0 
        log.write(" -Calculating p values based on given null slope :",null_beta, verbose=verbose)
        t_score = (reg[0]-null_beta) / reg[4]
        degree = len(sig_list_merged.dropna())-2
        p =  reg[3]
        #ss.t.sf(abs(t_score), df=degree)*2
        log.write(" -Beta = ", reg[0], verbose=verbose)
        log.write(" -Beta_se = ", reg[4], verbose=verbose)
        #log.write(" -H0 beta = ", null_beta, ", recalculated p = ", "{:.2e}".format(p), verbose=verbose)
        log.write(" -H0 beta =  0",", default p = ", "{:.2e}".format(reg[3]), verbose=verbose)
        log.write(" -Peason correlation coefficient =  ", "{:.2f}".format(reg[2]), verbose=verbose)
        log.write(" -r2 =  ", "{:.2f}".format(reg[2]**2), verbose=verbose)
        if r_se==True:
            log.write(" -R se (jackknife) = {:.2e}".format(r_se_jackknife), verbose=verbose)

        if reg[0] > 0:
            #if regression coeeficient >0 : auxiliary line slope = 1
            if is_45_helper_line is True:
                ax.axline([min(xl,yl),min(xl,yl)], [max(xh, yh),max(xh, yh)],zorder=1,**helper_line_args)

            #add text
            try:
                p12=str("{:.2e}".format(p)).split("e")[0]
                pe =str(int("{:.2e}".format(p).split("e")[1]))
            except:
                p12="0"
                pe="0"
            p_text="$p = " + p12 + " \\times  10^{"+pe+"}$"
            p_latex= f'{p_text}'
            ax.text(0.98,0.02,"$y =$ "+"{:.2f}".format(reg[1]) +" $+$ "+ "{:.2f}".format(reg[0])+" $x$, "+ p_latex + ", $r =$" +"{:.2f}".format(reg[2])+r_se_jackknife_string, va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **font_kwargs)
        else:
            #if regression coeeficient <0 : auxiliary line slope = -1
            if is_45_helper_line is True:
                if mode=="beta" or mode=="BETA" or mode=="Beta": 
                    ax.axline([min(xl,yl),-min(xl,yl)], [max(xh, yh),-max(xh, yh)],zorder=1,**helper_line_args)
                else:
                    ax.axline([min(xl,yl),-min(xl,yl)], [max(xh, yh),-max(xh, yh)],zorder=1,**helper_line_args)
            #add text
            try:
                p12=str("{:.2e}".format(p)).split("e")[0]
                pe =str(int("{:.2e}".format(p).split("e")[1]))
            except:
                p12="0"
                pe="0"
            p_text="$p = " + p12 + " \\times  10^{"+pe+"}$"
            p_latex= f'{p_text}'
            ax.text(0.98,0.02,"$y =$ "+"{:.2f}".format(reg[1]) +" $-$ "+ "{:.2f}".format(abs(reg[0]))+" $x$, "+ p_latex + ", $r =$" +"{:.2f}".format(reg[2])+r_se_jackknife_string, va="bottom",ha="right",transform=ax.transAxes,bbox=reg_box,**font_kwargs)
            
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            middle = sig_list_merged["EFFECT_1"].mean()
        else:
            middle = sig_list_merged["OR_1"].mean()
        
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            ax.axline(xy1=(0,reg[1]),slope=reg[0],color="#cccccc",linestyle='--',zorder=1)
        else:
            ax.axline(xy1=(1,reg[0]+reg[1]),slope=reg[0],color="#cccccc",linestyle='--',zorder=1)
    return ax


def configure_legend(fig, ax, legend_mode, is_q, is_q_mc, legend_elements, legend_pos, q_level, 
                     font_kwargs,scatterargs,legend_args,
                     legend_title, legend_title2  ):
    legend_args_to_use ={
            "framealpha":1,
            "handlelength":0.7,
            "handletextpad":0.8,
            "edgecolor":"grey",
            "borderpad":0.3,
            "alignment":"left"
        }

    if legend_args is not None:
        for key,value in legend_args.items():
            legend_args_to_use[key] = value

    if legend_mode == "full" and is_q==True :
        title_proxy = Rectangle((0,0), 0, 0, color='w',label=legend_title)
        title_proxy2 = Rectangle((0,0), 0, 0, color='w',label=legend_title2)
        if is_q_mc=="fdr":
            het_label_sig = r"$FDR_{het} < $" + "${}$".format(q_level)
            het_label_sig2 = r"$FDR_{het} > $" + "${}$".format(q_level)
        elif is_q_mc=="bon":
            het_label_sig = r"$P_{het,bon} < $" + "${}$".format(q_level)
            het_label_sig2 = r"$P_{het,bon} > $" + "${}$".format(q_level)
        else:
            het_label_sig = r"$P_{het} < $" + "${}$".format(q_level)
            het_label_sig2 = r"$P_{het} > $" + "${}$".format(q_level)
        het_sig = Rectangle((0,0), 0, 0, facecolor='#cccccc',edgecolor="black", linewidth=1, label=het_label_sig)
        het_nonsig = Rectangle((0,0), 0, 0, facecolor='#cccccc',edgecolor="white",linewidth=1, label=het_label_sig2)
        
        ax.add_patch(title_proxy)
        ax.add_patch(title_proxy2)
        ax.add_patch(het_sig)
        ax.add_patch(het_nonsig)

        legend_order = [legend_title] + legend_elements + [legend_title2] +[het_label_sig, het_label_sig2]
        handles, labels = reorderLegend(ax=ax, order=legend_order)
        
        #handles.append([het_sig,het_nonsig])
        #labels.append([het_label_sig,het_label_sig2])
        L = ax.legend(
            handles=handles, 
            labels=labels,
            #title=legend_title,
            loc=legend_pos,
            **legend_args_to_use
            )
    else:
        L = ax.legend(
            title=legend_title,
            loc=legend_pos,
            **legend_args_to_use
            )
    
    #for i, handle in enumerate(L.legendHandles):
    #    handle.set_edgecolor("white")

    ## Move titles to the left 
    for item, label in zip(L.legendHandles, L.texts):
        if label._text  in legend_elements:
            item.set_edgecolor("white")
            #item._legmarker.set_markersize(scatterargs["s"]*1.5)
            item._sizes = [scatterargs["s"]*2]
        if legend_mode == "full":
            if label._text  in [legend_title, legend_title2]:
                width=item.get_window_extent(fig.canvas.get_renderer()).width
                label.set_ha('left')
                label.set_position((-8*width,0))

    ax.tick_params(axis='both', labelsize=font_kwargs["fontsize"])
    plt.setp(L.texts,**font_kwargs)
    plt.setp(L.get_title(),**font_kwargs)
    return ax

def reorderLegend(ax=None, order=None, add=None):
    handles, labels = ax.get_legend_handles_labels()
    info = dict(zip(labels, handles))

    new_handles = [info[l] for l in order]
    return new_handles, order

def reorder_columns(sig_list_merged):
    order=[ 'CHR', 'POS', 'GENENAME', 
            'EA_1', 'NEA_1', 'EFFECT_1', 'SE_1', 'P_1', 'MLOG10P_1', 
            'EA_2_aligned','NEA_2_aligned', 'EFFECT_2_aligned', 'SE_2','P_2','MLOG10P_2',  'EA_2', 'NEA_2', 'EFFECT_2', 
            'indicator' ]
    
    new_order=[]
    for i in order:
        if i in sig_list_merged.columns:
            new_order.append(i)
    for i in sig_list_merged.columns:
        if i not in new_order:
            new_order.append(i)
    
    return sig_list_merged[new_order]