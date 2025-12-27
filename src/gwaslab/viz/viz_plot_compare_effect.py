import gc
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

from gwaslab.util.util_in_get_sig import _get_sig, _anno_gene
from gwaslab.util.util_in_correct_winnerscurse import wc_correct
from gwaslab.util.util_in_correct_winnerscurse import wc_correct_test
from gwaslab.viz.viz_plot_scatter_with_reg import jackknife_r as scatter_jackknife_r
from gwaslab.viz.viz_plot_scatter_with_reg import create_helper_line as _create_helper_line
from gwaslab.viz.viz_plot_scatter_with_reg import create_reg_line as _create_reg_line
from gwaslab.viz.viz_plot_scatter_with_reg import plot_scatter_with_err

#20220422
def compare_effect(path1,
                   path2,
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
                   anno_kwargs=None,
                   wc_correction=False, 
                   null_beta=0,
                   is_q=False,
                   is_q_mc = False,
                   include_all=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   get_lead_kwargs=None,
                   drop=False,
                   wc_sig_level=5e-8,
                   # reg
                   reg_box=None,
                   is_reg=True,
                   fdr=False,
                   reg_text="full",
                   allele_match=False,
                   r_se=False,
                   is_45_helper_line=True,
                   legend_mode="full",
                   legend_title=r'$\mathregular{ P < 5 x 10^{-8}}$ in:',
                   legend_title2=r'Heterogeneity test:',
                   legend_pos='upper left',
                   scatter_kwargs=None,
                   fig_kwargs=None,
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_kwargs=None,
                   adjust_text_kwargs = None,
                   adjust_text_kwargs_l = None,
                   adjust_text_kwargs_r = None,
                   font_kwargs=None,
                   build="19",
                   r_or_r2="r",
                   err_kwargs=None,
                   legend_kwargs=None,
                   clean_output=False,
                   log = Log(),
                   save=False,
                   save_kwargs=None,
                   verbose=False,
                   **kwargs):

    # Input column examples:
    # [snpid,p,ea,nea] / [effect,se]
    # [snpid,p,ea,nea,chr,pos] / [effect,se]
    # [snpid,p,ea,nea,chr,pos] / [OR,OR_l,OR_h]

    # Auto-format legend title based on significance threshold
    if legend_title == r'$\mathregular{ P < 5 x 10^{-8}}$ in:' and sig_level != 5e-8:
        exponent = math.floor(math.log10(sig_level))
        mantissa = sig_level / 10 ** exponent
        legend_title = '$\mathregular{{ P < {} x 10^{{{}}} }}$ in:'.format(mantissa, exponent)

    # Enable heterogeneity flag when using multiple correction
    if is_q_mc in {"fdr", "bon"}:
        is_q = True

    # Validate heterogeneity options
    if is_q:
        if is_q_mc not in [False, "fdr", "bon", "non"]:
            raise ValueError('Please select either "fdr" or "bon" or "non"/False for is_q_mc.')

    # Merge style
    style = set_plot_style(
        plot="compare_effect",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs,
        err_kwargs=err_kwargs,
        font_kwargs=font_kwargs,
        helper_line_kwargs=helper_line_kwargs,
        anno_kwargs=anno_kwargs,
        fontsize=font_kwargs.get('fontsize') if font_kwargs else None,
        fontfamily=None,
        verbose=verbose,
        log=log,
    )

    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    scatter_kwargs = style.get("scatter_kwargs", {})
    err_kwargs = style.get("err_kwargs", {})
    font_kwargs = style.get("font_kwargs", {})
    helper_line_kwargs = style.get("helper_line_kwargs", {})
    anno_kwargs = style.get("anno_kwargs", {})

    if reg_box is None:
        reg_box= None
    elif reg_box==True:
        reg_box = dict(boxstyle='round', facecolor='white', alpha=1,edgecolor="grey")
    if get_lead_kwargs is None:
        get_lead_kwargs = {}
    if anno=="GENENAME":
        get_lead_kwargs["anno"]=True
    if label is None:
        label = ["Sumstats_1","Sumstats_2","Both","None"]
    if anno_het:
        is_q = True
    
    adjust_text_kwargs_r_default = {"autoalign":False,"precision":0.001,"lim":1000,"ha":"left","va":"top","expand_text":(1,1.8),"expand_objects":(0.1,0.1),"expand_points":(1.8,1.8),"force_objects":(0.8,0.8),"arrowprops":dict(arrowstyle='-|>', color='grey')}
    adjust_text_kwargs_l_default = {"autoalign":False,"precision":0.001,"lim":1000,"ha":"right","va":"bottom","expand_text":(1,1.8),"expand_objects":(0.1,0.1),"expand_points":(1.8,1.8),"force_objects":(0.8,0.8),"arrowprops":dict(arrowstyle='-|>', color='grey')}

    base_l = dict(adjust_text_kwargs_l_default)
    base_r = dict(adjust_text_kwargs_r_default)
    
    if isinstance(adjust_text_kwargs_l, dict):
        base_l.update({k: v for k, v in adjust_text_kwargs_l.items() if v is not None})
    if isinstance(adjust_text_kwargs_r, dict):
        base_r.update({k: v for k, v in adjust_text_kwargs_r.items() if v is not None})
    if isinstance(adjust_text_kwargs, dict):
        base_l.update({k: v for k, v in adjust_text_kwargs.items() if v is not None})
        base_r.update({k: v for k, v in adjust_text_kwargs.items() if v is not None})
    adjust_text_kwargs_l = base_l
    adjust_text_kwargs_r = base_r

    # 0. Start processing and validate inputs
    log.write("Start to process the raw sumstats for plotting...", verbose=verbose)

    if not isinstance(path1, Sumstats) or not isinstance(path2, Sumstats):
        raise ValueError("Please provide GWASLab Sumstats objects for `path1` and `path2`.")
    scaled1 = "MLOG10P" in path1.data.columns
    scaled2 = "MLOG10P" in path2.data.columns
    has_eaf1 = "EAF" in path1.data.columns
    has_eaf2 = "EAF" in path2.data.columns
    # 1. Load sumstats1 minimal columns
    sumstats, common_snp_set = configure_common_snp_set(path1, path2,
                               snplist,
                               label,
                               log, verbose)
    
    # 2. Rename sumstats1 headers to GWASLab keywords
    sumstats = rename_sumtats(sumstats=sumstats,
                              snplist=snplist,
                              scaled=scaled1)
    
    # 3. Extract SNPs for comparison (provided list or auto lead SNPs)
    sig_list_1 = extract_snp_for_comparison(sumstats, 
                                            snplist, 
                                            label=label[0], 
                                            get_lead_kwargs=get_lead_kwargs, 
                                            build=build,
                                            drop=drop, 
                                            anno=anno, 
                                            sig_level=sig_level, 
                                            scaled = scaled1,
                                            log = log, 
                                            verbose = verbose)

 
    ######### load sumstats1

    # 4. Load sumstats2 minimal columns
    if snplist is not None:
        cols_to_extract = ["SNPID", "MLOG10P" if scaled2 else "P"]
    else:
        cols_to_extract = ["SNPID", "MLOG10P" if scaled2 else "P", "CHR", "POS"]
    
    sumstats = load_sumstats(path=path2, 
                            usecols=cols_to_extract, 
                            label=label[1], 
                            log=log, 
                            verbose= verbose)
    gc.collect()
    
    #if scaled2==True:
    #    sumstats[cols_name_list_2[1]] = np.power(10,-sumstats[cols_name_list_2[1]])

    # 5. Rename sumstats2 headers to GWASLab keywords
    sumstats = rename_sumtats(sumstats=sumstats,
                              snplist=snplist,
                              scaled=scaled2)
    # 6. Extract SNPs from sumstats2
    sig_list_2 = extract_snp_for_comparison(sumstats, 
                                        snplist, 
                                        label=label[1], 
                                        get_lead_kwargs=get_lead_kwargs, 
                                        build=build,
                                        drop=drop, 
                                        anno=anno, 
                                        sig_level=sig_level, 
                                        scaled = scaled2,
                                        log = log, 
                                        verbose = verbose)

    # 7. Merge SNP lists by SNPID
    sig_list_merged = merge_list(sig_list_1, 
                                 sig_list_2, 
                                 anno = anno, 
                                 labels=label,
                                 log=log, 
                                 verbose=verbose)

    # 8. Load effect columns for sumstats1
    cols_to_extract = configure_cols_to_extract(mode=mode, has_eaf=has_eaf1, scaled=scaled1)
    sumstats = load_sumstats(path=path1, 
                            usecols=cols_to_extract, 
                            label=label[0], 
                            log=log, 
                            verbose= verbose)
    
    #if scaled1==True:
    #    sumstats[cols_name_list_1[1]] = np.power(10,-sumstats[cols_name_list_1[1]])
    # 9. Rename and clean sumstats1 effects
    sumstats = rename_sumstats_full(mode, sumstats,
                                    index=1, 
                                    drop = drop, 
                                    scaled=scaled1,
                                    has_eaf=has_eaf1,
                                    log=log, verbose=verbose)

    # 10. Merge effects from sumstats1
    log.write(" -Merging "+label[0]+" effect information...", verbose=verbose)
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")

    # 11. Load effect columns for sumstats2
    cols_to_extract = configure_cols_to_extract(mode=mode, has_eaf=has_eaf2, scaled=scaled2)
    
    sumstats = load_sumstats(path=path2, 
                            usecols=cols_to_extract, 
                            label=label[1], 
                            log=log, 
                            verbose= verbose)
    
    #if scaled2==True:
    #    sumstats[cols_name_list_2[1]] = np.power(10,-sumstats[cols_name_list_2[1]])
    
    gc.collect()
    
    # 12. Rename and clean sumstats2 effects
    sumstats = rename_sumstats_full(mode, sumstats,
                                    index=2, 
                                    drop = drop, 
                                    scaled=scaled2,
                                    has_eaf=has_eaf2,
                                    log=log, verbose=verbose)

    # 13. Merge effects from sumstats2
    log.write(" -Merging "+label[1]+" effect information...", verbose=verbose)
    sig_list_merged = pd.merge(sig_list_merged,sumstats,
                               left_on="SNPID",right_on="SNPID",
                               how="left")
    
    sig_list_merged.set_index("SNPID",inplace=True)

    # 14. Update p-values for sumstats1

    sig_list_merged = update_stats(sig_list_merged = sig_list_merged, 
                                   path = path1, 
                                   index=1, 
                                   snplist = snplist, 
                                   label=label[0], 
                                   drop = drop, 
                                   scaled=scaled1,
                                   log=log, 
                                   verbose = verbose)
    
    # 15. Update p-values for sumstats2
    sig_list_merged = update_stats(sig_list_merged = sig_list_merged, 
                                   path = path2, 
                                   index=2, 
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
    
    # 16. Assign significance indicators
    sig_list_merged = assign_indicator(sig_list_merged, snplist, sig_level, scaled1, scaled2, log, verbose)
    # 17. Align alleles across datasets
    sig_list_merged = align_alleles(sig_list_merged, label, mode, log, verbose)
    # 18. Check allele match consistency
    sig_list_merged = check_allele_match(sig_list_merged, allele_match, label, log,verbose)
    # 19. Filter by MAF if EAF present
    sig_list_merged = filter_by_maf(sig_list_merged, maf_level, log, verbose)

    if fdr is True and (not scaled1) and (not scaled2):
        log.write(" -Using FDR...", verbose=verbose)
        #sig_list_merged["P_1"] = fdrcorrection(sig_list_merged["P_1"])[1]
        #sig_list_merged["P_2"] = fdrcorrection(sig_list_merged["P_2"])[1]
        sig_list_merged["P_1"] =ss.false_discovery_control(sig_list_merged["P_1"])
        sig_list_merged["P_2"] =ss.false_discovery_control(sig_list_merged["P_2"])

    # 20. Winner's curse correction (optional)
    sig_list_merged = winnerscurse_correction(sig_list_merged, mode, wc_correction, sig_level,scaled1, scaled2, log, verbose)

    # 21. Heterogeneity test (optional)
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
    fig,ax = plt.subplots(**fig_kwargs)
    legend_elements=[]
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        if len(sum0)>0:
            plot_scatter_with_err(ax=ax, df=sum0, x="EFFECT_1", y="EFFECT_2_aligned",
                                  xerr=sum0["SE_1"], yerr=sum0["SE_2"], engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[3], "zorder": 2, "color": "#cccccc", "edgecolors": sum0["Edge_color"], "marker": "."},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[3])
        if len(sum1only)>0:
            plot_scatter_with_err(ax=ax, df=sum1only, x="EFFECT_1", y="EFFECT_2_aligned",
                                  xerr=sum1only["SE_1"], yerr=sum1only["SE_2"], engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[0], "zorder": 2, "color": "#e6320e", "edgecolors": sum1only["Edge_color"], "marker": "^"},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[0])
        if len(sum2only)>0:
            plot_scatter_with_err(ax=ax, df=sum2only, x="EFFECT_1", y="EFFECT_2_aligned",
                                  xerr=sum2only["SE_1"], yerr=sum2only["SE_2"], engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[1], "zorder": 2, "color": "#41e620", "edgecolors": sum2only["Edge_color"], "marker": "o"},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[1])
        if len(both)>0:
            plot_scatter_with_err(ax=ax, df=both, x="EFFECT_1", y="EFFECT_2_aligned",
                                  xerr=both["SE_1"], yerr=both["SE_2"], engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[2], "zorder": 2, "color": "#205be6", "edgecolors": both["Edge_color"], "marker": "s"},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[2])
    else:
        ## if OR
        if len(sum0)>0:
            plot_scatter_with_err(ax=ax, df=sum0, x="OR_1", y="OR_2_aligned",
                                  xerr=sum0[["OR_L_1_err","OR_H_1_err"]].T, yerr=sum0[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T, engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[3], "zorder": 2, "color": "#cccccc", "edgecolors": sum0["Edge_color"], "marker": "."},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[3])
        if len(sum1only)>0:
            plot_scatter_with_err(ax=ax, df=sum1only, x="OR_1", y="OR_2_aligned",
                                  xerr=sum1only[["OR_L_1_err","OR_H_1_err"]].T, yerr=sum1only[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T, engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[0], "zorder": 2, "color": "#e6320e", "edgecolors": sum1only["Edge_color"], "marker": "^"},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[0])
        if len(sum2only)>0:
            plot_scatter_with_err(ax=ax, df=sum2only, x="OR_1", y="OR_2_aligned",
                                  xerr=sum2only[["OR_L_1_err","OR_H_1_err"]].T, yerr=sum2only[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T, engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[1], "zorder": 2, "color": "#41e620", "edgecolors": sum2only["Edge_color"], "marker": "o"},
                                  err_kwargs=err_kwargs)
            legend_elements.append(label[1])
        if len(both)>0:
            plot_scatter_with_err(ax=ax, df=both, x="OR_1", y="OR_2_aligned",
                                  xerr=both[["OR_L_1_err","OR_H_1_err"]].T, yerr=both[["OR_L_2_aligned_err","OR_H_2_aligned_err"]].T, engine="plt",
                                  scatter_kwargs={**scatter_kwargs, "label": label[2], "zorder": 2, "color": "#205be6", "edgecolors": both["Edge_color"], "marker": "s"},
                                  err_kwargs=err_kwargs)
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
        ax.axhline(y=0, zorder=1,**helper_line_kwargs)
        ax.axvline(x=0, zorder=1,**helper_line_kwargs)
    else:
        #if using OR
        ax.axhline(y=1, zorder=1,**helper_line_kwargs)
        ax.axvline(x=1, zorder=1,**helper_line_kwargs)
    
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    ###regression line##############################################################################################################################
    ax = configure_regression_line(is_reg=is_reg,
                                 reg_text=reg_text,
                                 reg_box=reg_box, 
                                 sig_list_merged=sig_list_merged, 
                                 ax=ax, 
                                 mode=mode,
                                 xl=xl,
                                 yl=yl,
                                 xh=xh,
                                 yh=yh, 
                                 null_beta=null_beta, 
                                 r_se=r_se, 
                                 is_45_helper_line=is_45_helper_line,
                                 helper_line_kwargs=helper_line_kwargs, 
                                 font_kwargs=font_kwargs,
                                 log=log, 
                                 verbose=verbose)
        
    
    ax.set_xlabel(xylabel_prefix+label[0],**font_kwargs)
    ax.set_ylabel(xylabel_prefix+label[1],**font_kwargs)
    
    ax = configure_legend(fig, ax, legend_mode, is_q, is_q_mc, legend_elements, legend_pos, q_level, 
                        font_kwargs,scatter_kwargs,legend_kwargs,
                        legend_title, legend_title2 )
    ##plot finished########################################################################################
    gc.collect()

    save_figure(fig, save, keyword="esc",save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    # Final reorder and save merged data at the end
    if clean_output:
        sig_list_merged = reorder_columns_clean(sig_list_merged, mode)
    else:
        sig_list_merged = reorder_columns(sig_list_merged)

    save_path = label[0]+"_"+label[1]+"_beta_sig_list_merged.tsv"
    log.write(" -Saving the merged data to:",save_path, verbose=verbose)
    sig_list_merged.to_csv(save_path,sep="\t")

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

def load_sumstats(path, usecols, label, log, verbose):
    if type(usecols) is not list:
        usecols = [usecols]
    
    log.write(" -Loading sumstats for {} : {}".format(label,",".join(usecols)), verbose=verbose)
    if isinstance(path, Sumstats):
        sumstats = path.data.loc[:,usecols].copy()
    else:
        raise ValueError("Please provide GWASLab Sumstats objects; file paths or DataFrames are not supported.")
    return sumstats

def configure_headers(mode, scaled1, scaled2, log, verbose):
    if mode not in ["Beta","beta","BETA","OR","or"]:
        raise ValueError("Please input Beta or OR")
    log.write("Using fixed headers from Sumstats object.", verbose=verbose)
    cols_name_list_1 = ["SNPID", "MLOG10P" if scaled1 else "P", "EA", "NEA", "CHR", "POS"]
    cols_name_list_2 = ["SNPID", "MLOG10P" if scaled2 else "P", "EA", "NEA", "CHR", "POS"]
    if mode == "beta" or mode == "Beta" or mode == "BETA":
        effect_cols_list_1 = ["BETA", "SE"]
        effect_cols_list_2 = ["BETA", "SE"]
    else:
        effect_cols_list_1 = ["OR", "OR_95L", "OR_95U"]
        effect_cols_list_2 = ["OR", "OR_95L", "OR_95U"]
    return cols_name_list_1, cols_name_list_2, effect_cols_list_1, effect_cols_list_2

def configure_common_snp_set(path1, path2,
                             snplist,
                             label,
                             log, verbose):
    # favor snplist if provided
    if snplist is not None:
        log.write(" -Using provided snplist for intersection...", verbose=verbose)
        # ensure set for fast lookup
        snp_set = set(snplist)
        # load SNPID only to reduce memory IO
        snp2 = load_sumstats(path=path2, usecols="SNPID", label=label[1], log=log, verbose=verbose)
        common_snp_set = snp_set.intersection(snp2["SNPID"].values)
        # prepare columns from path1 according to downstream needs
        use_pcol = "MLOG10P" if ("MLOG10P" in path1.data.columns) else "P"
        cols_to_extract = ["SNPID", use_pcol]
    else:
        snp2 = load_sumstats(path=path2, usecols="SNPID", label=label[1], log=log, verbose=verbose)
        common_snp_set = set(snp2["SNPID"].values)
        use_pcol = "MLOG10P" if ("MLOG10P" in path1.data.columns) else "P"
        cols_to_extract = ["SNPID", use_pcol, "CHR", "POS"]

    sumstats = load_sumstats(path=path1, usecols=cols_to_extract, label=label[0], log=log, verbose=verbose)
    gc.collect()

    # intersect with path1 SNPs
    common_snp_set = common_snp_set.intersection(sumstats["SNPID"].values)

    log.write(" -Counting  variants available for both datasets:", len(common_snp_set), " variants...", verbose=verbose)

    return sumstats, common_snp_set

def rename_sumtats(sumstats, snplist, scaled, suffix=""):
    rename_dict = {"SNPID": "SNPID"}
    if scaled:
        rename_dict["MLOG10P"] = "MLOG10P{}".format(suffix)
    else:
        rename_dict["P"] = "P{}".format(suffix)
    if snplist is None:
        rename_dict["CHR"] = "CHR"
        rename_dict["POS"] = "POS"
    return sumstats.rename(columns=rename_dict)


def extract_snp_for_comparison(sumstats, snplist, label, 
                               get_lead_kwargs, build, drop, anno, 
                               sig_level,scaled, log, verbose):
    ######### 8 extact SNPs for comparison 
    if snplist is not None: 
        ######### 8.1 if a snplist is provided, use the snp list
        log.write(" -Extract variants in the given list from "+label+"...")
        sig_list = sumstats.loc[sumstats["SNPID"].isin(snplist),:].copy()
        if anno=="GENENAME":
            sig_list = _anno_gene(sig_list,"SNPID","CHR","POS", build=build, verbose=verbose, **get_lead_kwargs)
    else:
        ######### 8,2 otherwise use the automatically detected lead SNPs
        log.write(" -Extract lead variants from "+label +"...", verbose=verbose)
        sig_list = _get_sig(sumstats,"SNPID","CHR","POS","P","MLOG10P", build=build, verbose=verbose,sig_level=sig_level,**get_lead_kwargs)
    
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

def configure_cols_to_extract(mode, has_eaf, scaled):
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        cols = ["SNPID", "MLOG10P" if scaled else "P", "EA", "NEA", "BETA", "SE"]
    else:
        cols = ["SNPID", "MLOG10P" if scaled else "P", "EA", "NEA", "OR", "OR_95L", "OR_95U"]
    if has_eaf:
        cols.append("EAF")
    return cols

def rename_sumstats_full(mode, sumstats, drop, index, scaled, has_eaf, log, verbose):
    if mode=="beta" or mode=="BETA" or mode=="Beta":
        rename_dict = {
            "SNPID": "SNPID",
            ("MLOG10P" if scaled else "P"): ("MLOG10P_{}".format(index) if scaled else "P_{}".format(index)),
            "EA": "EA_{}".format(index),
            "NEA": "NEA_{}".format(index),
            "BETA": "EFFECT_{}".format(index),
            "SE": "SE_{}".format(index),
        }
    else:
        rename_dict = {
            "SNPID": "SNPID",
            ("MLOG10P" if scaled else "P"): ("MLOG10P_{}".format(index) if scaled else "P_{}".format(index)),
            "EA": "EA_{}".format(index),
            "NEA": "NEA_{}".format(index),
            "OR": "OR_{}".format(index),
            "OR_95L": "OR_L_{}".format(index),
            "OR_95U": "OR_H_{}".format(index),
        }
    if has_eaf:
        rename_dict["EAF"] = "EAF_{}".format(index)
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
                 snplist, 
                 label, 
                 drop, 
                 index, 
                 scaled,
                 log, 
                 verbose):
    
    log.write(" -Updating missing information for "+label+" ...", verbose=verbose)
    cols_to_extract = ["SNPID", "MLOG10P" if scaled else "P"]
    
    sumstats = load_sumstats(path=path, 
                             usecols=cols_to_extract, 
                             label=label, 
                             log=log, 
                             verbose= verbose)
    #if scaled1==True:
    #    sumstats[cols_name_list_1[1]] = np.power(10,-sumstats[cols_name_list_1[1]])
    
    sumstats = rename_sumtats(sumstats = sumstats, 
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

def align_alleles(sig_list_merged, label, mode, log, verbose):
    log.write(" -Aligning "+label[1]+" EA with "+label[0]+" EA ...", verbose=verbose)
    cols = ["EA_1", "EA_2", "NEA_1", "NEA_2"]
    sig_list_merged[cols] = sig_list_merged[cols].astype("string")
    mismatch = sig_list_merged["EA_1"] != sig_list_merged["EA_2"]
    if mode in ("beta", "BETA", "Beta"):
        sig_list_merged["EA_2_aligned"] = np.where(mismatch, sig_list_merged["NEA_2"], sig_list_merged["EA_2"])
        sig_list_merged["NEA_2_aligned"] = np.where(mismatch, sig_list_merged["EA_2"], sig_list_merged["NEA_2"])
        sig_list_merged["EFFECT_2_aligned"] = np.where(mismatch, -sig_list_merged["EFFECT_2"], sig_list_merged["EFFECT_2"])
    else:
        sig_list_merged["EA_2_aligned"] = np.where(mismatch, sig_list_merged["NEA_2"], sig_list_merged["EA_2"])
        sig_list_merged["NEA_2_aligned"] = np.where(mismatch, sig_list_merged["EA_2"], sig_list_merged["NEA_2"])
        sig_list_merged["OR_2_aligned"] = np.where(mismatch, 1 / sig_list_merged["OR_2"], sig_list_merged["OR_2"])
        sig_list_merged["OR_L_2_aligned"] = np.where(mismatch, 1 / sig_list_merged["OR_H_2"], sig_list_merged["OR_L_2"])
        sig_list_merged["OR_H_2_aligned"] = np.where(mismatch, 1 / sig_list_merged["OR_L_2"], sig_list_merged["OR_H_2"])
        sig_list_merged["BETA_1"] = np.log(sig_list_merged["OR_1"])
        sig_list_merged["BETA_2_aligned"] = np.log(sig_list_merged["OR_2_aligned"])
        sig_list_merged["SE_1"] = (np.log(sig_list_merged["OR_H_1"]) - np.log(sig_list_merged["OR_1"])) / ss.norm.ppf(0.975)
        sig_list_merged["SE_2"] = (np.log(sig_list_merged["OR_H_2_aligned"]) - np.log(sig_list_merged["OR_2_aligned"])) / ss.norm.ppf(0.975)
        sig_list_merged["OR_L_1_err"] = np.abs(sig_list_merged["OR_L_1"] - sig_list_merged["OR_1"])
        sig_list_merged["OR_H_1_err"] = np.abs(sig_list_merged["OR_H_1"] - sig_list_merged["OR_1"])
        sig_list_merged["OR_L_2_aligned_err"] = np.abs(sig_list_merged["OR_L_2_aligned"] - sig_list_merged["OR_2_aligned"])
        sig_list_merged["OR_H_2_aligned_err"] = np.abs(sig_list_merged["OR_H_2_aligned"] - sig_list_merged["OR_2_aligned"])
    if "EAF_2" in sig_list_merged.columns:
        sig_list_merged["EAF_2_aligned"] = np.where(mismatch, 1 - sig_list_merged["EAF_2"], sig_list_merged["EAF_2"])
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

def filter_by_maf(sig_list_merged, maf_level, log, verbose):
    if maf_level is not None and "EAF_1" in sig_list_merged.columns and "EAF_2" in sig_list_merged.columns:
        both_eaf_clear =  (sig_list_merged["EAF_1"]>maf_level)&(sig_list_merged["EAF_1"]<1-maf_level)&(sig_list_merged["EAF_2"]>maf_level)&(sig_list_merged["EAF_2"]<1-maf_level)
        log.write(" -Exclude "+str(len(sig_list_merged) -sum(both_eaf_clear))+ " variants with maf <"+str(maf_level), verbose=verbose)
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
        # P value upper bound -> 1
        df.loc[df[pq]>1,pq] = 1

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
    xcol = "EFFECT_1" if (mode=="beta" or mode=="BETA" or mode=="Beta") else "OR_1"
    ycol = "EFFECT_2_aligned" if (mode=="beta" or mode=="BETA" or mode=="Beta") else "OR_2_aligned"
    df = sig_list_merged.dropna(axis=0)
    if is_q==True and anno_het == True:
        df = df.loc[df["Edge_color"]=="black",:]
    m = (df[xcol].abs()>=anno_min1) & (df[ycol].abs()>=anno_min2) & (df[xcol].abs()>=anno_min) & (df[ycol].abs()>=anno_min) & (np.abs(df[xcol]-df[ycol])>=anno_diff)
    df = df.loc[m,:]
    texts_l=[]
    texts_r=[]
    if anno==True:
        log.write("Annotating variants using {}".format("SNPID"), verbose=verbose)
    elif anno=="GENENAME":
        log.write("Annotating variants using {}".format("GENENAME"), verbose=verbose)
    if anno==False:
        return ax

    if isinstance(anno, dict):
        df = df.loc[df.index.isin(list(anno.keys())),:]
        iterator = ((index, row, anno[index]) for index, row in df.iterrows())
    else:
        def pick_text(index, row):
            if anno==True:
                return index
            elif type(anno) is str and not pd.isna(row[anno]):
                return row[anno]
            else:
                return index
        iterator = ((index, row, pick_text(index, row)) for index, row in df.iterrows())
    for index, row, t in iterator:
        if row[xcol] < row[ycol]:
            texts_l.append(plt.text(row[xcol], row[ycol], t, ha="right", va="bottom", **anno_kwargs))
        else:
            texts_r.append(plt.text(row[xcol], row[ycol], t, ha="left", va="top", **anno_kwargs))
    if len(texts_l)>0:
        adjust_text(texts_l, ax=ax, **adjust_text_kwargs_l)
    if len(texts_r)>0:
        adjust_text(texts_r, ax=ax, **adjust_text_kwargs_r)
    return ax


def configure_regression_line(is_reg, 
                              reg_box, 
                              reg_text,
                              sig_list_merged,  
                              ax, 
                              mode,
                              xl,
                              yl,
                              xh,
                              yh, 
                              null_beta, 
                              r_se, 
                            is_45_helper_line,
                            helper_line_kwargs, 
                            font_kwargs,
                            log, 
                            verbose):
    if len(sig_list_merged)<3: is_reg=False
    if is_reg is True:
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            reg = ss.linregress(sig_list_merged["EFFECT_1"],sig_list_merged["EFFECT_2_aligned"])
            
            # estimate se for r
            if r_se==True:
                log.write(" -Estimating SE for rsq using Jackknife method.", verbose=verbose)
                r_se_jackknife = scatter_jackknife_r(sig_list_merged, x="EFFECT_1", y="EFFECT_2_aligned", log=log, verbose=verbose)
                r_se_jackknife_string = " ({:.2f})".format(r_se_jackknife)
            else:
                r_se_jackknife_string= ""
        else:
            reg = ss.linregress(sig_list_merged["OR_1"],sig_list_merged["OR_2_aligned"])
            r_se_jackknife_string= ""

        #### calculate p values based on selected value , default = 0 
        log.write(" -Calculating p values based on given null slope :",null_beta, verbose=verbose)
        #t_score = (reg[0]-null_beta) / reg[4]
        #degree = len(sig_list_merged.dropna())-2
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
                _create_helper_line(ax, reg[0], is_45_helper_line, helper_line_kwargs)

            #add text
            try:
                p12=str("{:.2e}".format(p)).split("e")[0]
                pe =str(int("{:.2e}".format(p).split("e")[1]))
            except:
                p12="0"
                pe="0"
            if p > 1e-300:
                p_text="$\mathregular{p = " + p12 + " \\times  10^{"+pe+"}}$"
            else:
                p_text="$\mathregular{p < 1 \\times 10^{-300}}$"
            p_latex= f'{p_text}'

            if reg_text=="full":
                reg_string = "y = "+"{:.2f}".format(reg[1]) +" + "+ "{:.2f}".format(reg[0])+" x, "+ p_latex + ", r = " +"{:.2f}".format(reg[2])+r_se_jackknife_string
                ax.text(0.98,0.02,
                        reg_string,
                        va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **font_kwargs)
            elif reg_text=="r":
                reg_string ="r = " +"{:.2f}".format(reg[2])+r_se_jackknife_string
                ax.text(0.98,0.02,
                        reg_string, va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **font_kwargs)
            elif reg_text=="r2":
                reg_string = "$\mathregular{r^{2}} = " +"{:.2f}".format(reg[2]**2)
                ax.text(0.98,0.02,
                        reg_string, va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **font_kwargs)
        else:
            #if regression coeeficient <0 : auxiliary line slope = -1
            if is_45_helper_line is True:
                _create_helper_line(ax, reg[0], is_45_helper_line, helper_line_kwargs)
            #add text
            try:
                p12=str("{:.2e}".format(p)).split("e")[0]
                pe =str(int("{:.2e}".format(p).split("e")[1]))
            except:
                p12="0"
                pe="0"

            if p > 1e-300:
                p_text="$\mathregular{p = " + p12 + " \\times  10^{"+pe+"}}$"
            else:
                p_text="$\mathregular{p < 1 \\times 10^{-300}}$"
            p_latex= f'{p_text}'

            if reg_text=="full":
                ax.text(0.98,0.02,
                    "y = "+"{:.2f}".format(reg[1]) +" - "+ "{:.2f}".format(abs(reg[0]))+" x, "+ p_latex + ", r = " +"{:.2f}".format(reg[2])+r_se_jackknife_string, 
                    va="bottom",ha="right",transform=ax.transAxes,bbox=reg_box,**font_kwargs)
            elif reg_text=="r":
                ax.text(0.98,0.02,
                        "r = " +"{:.2f}".format(reg[2])+r_se_jackknife_string, 
                        va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **font_kwargs)
            elif reg_text=="r2":
                ax.text(0.98,0.02,
                        "$\mathregular{r^{2}} = " +"{:.2f}".format(reg[2]**2), 
                        va="bottom",ha="right",transform=ax.transAxes, bbox=reg_box, **font_kwargs)
            
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            middle = sig_list_merged["EFFECT_1"].mean()
        else:
            middle = sig_list_merged["OR_1"].mean()
        
        if mode=="beta" or mode=="BETA" or mode=="Beta":
            _create_reg_line(ax, reg, reg_xmin=0)
        else:
            _create_reg_line(ax, reg, reg_xmin=1)
    return ax


def configure_legend(fig, ax, legend_mode, is_q, is_q_mc, legend_elements, legend_pos, q_level, 
                     font_kwargs, scatter_kwargs, legend_kwargs,
                     legend_title, legend_title2  ):
    legend_kwargs_to_use ={
            "framealpha":1,
            "handlelength":0.7,
            "handletextpad":0.8,
            "edgecolor":"grey",
            "borderpad":0.3,
            "alignment":"left",
            "frameon":False
        }

    if legend_kwargs is not None:
        for key,value in legend_kwargs.items():
            legend_kwargs_to_use[key] = value

    if legend_mode == "full" and is_q==True :
        title_proxy = Rectangle((0,0), 0, 0, color='w',label=legend_title)
        title_proxy2 = Rectangle((0,0), 0, 0, color='w',label=legend_title2)
        if is_q_mc=="fdr":
            het_label_sig = r"$\mathregular{FDR_{het} < }$" + "{}".format(q_level)
            het_label_sig2 = r"$\mathregular{FDR_{het} > }$" + "{}".format(q_level)
        elif is_q_mc=="bon":
            het_label_sig = r"$\mathregular{P_{het,bon} < }$" + "{}".format(q_level)
            het_label_sig2 = r"$\mathregular{P_{het,bon} > }$" + "{}".format(q_level)
        else:
            het_label_sig = r"$\mathregular{P_{het} < }$" + "{}".format(q_level)
            het_label_sig2 = r"$\mathregular{P_{het} > }$" + "{}".format(q_level)
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
            **legend_kwargs_to_use
            )
    else:
        L = ax.legend(
            title=legend_title,
            loc=legend_pos,
            **legend_kwargs_to_use
            )
    
    #for i, handle in enumerate(L.legend_handles):
    #    handle.set_edgecolor("white")

    ## Move titles to the left 
    try:
        for item, label in zip(L.legend_handles, L.texts):
            if label._text  in legend_elements:
                item.set_edgecolor("white")
                #item._legmarker.set_markersize(scatter_kwargs["s"]*1.5)
                item._sizes = [scatter_kwargs["s"]*2]
            if legend_mode == "full":
                if label._text  in [legend_title, legend_title2]:
                    width=item.get_window_extent(fig.canvas.get_renderer()).width
                    label.set_ha('left')
                    label.set_position((-8*width,0))
    except:
        pass
    
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

def reorder_columns_clean(sig_list_merged, mode):
    base = sig_list_merged.copy()
    if 'EA_1' in base.columns and 'EA' not in base.columns:
        base['EA'] = base['EA_1']
    if 'NEA_1' in base.columns and 'NEA' not in base.columns:
        base['NEA'] = base['NEA_1']
    if mode in ("beta", "BETA", "Beta"):
        if 'EFFECT_2_aligned' in base.columns:
            base['BETA_2'] = base['EFFECT_2_aligned']
        if 'BETA_1' not in base.columns and 'EFFECT_1' in base.columns:
            base['BETA_1'] = base['EFFECT_1']
        clean_cols = ['CHR','POS','GENENAME','EA','NEA','BETA_1','SE_1','P_1','MLOG10P_1','BETA_2','SE_2','P_2','MLOG10P_2','indicator']
    else:
        if 'OR_2_aligned' in base.columns:
            base['OR_2'] = base['OR_2_aligned']
        clean_cols = ['CHR','POS','GENENAME','EA','NEA','OR_1','SE_1','P_1','MLOG10P_1','OR_2','SE_2','P_2','MLOG10P_2','indicator']

    # Return only the clean columns that exist
    return base[[c for c in clean_cols if c in base.columns]].copy()
