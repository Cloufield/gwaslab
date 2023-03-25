import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
from math import ceil
from shutil import which
from gwaslab.Log import Log
from gwaslab.calculate_gc import lambdaGC
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
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import MaxNLocator
import gc as garbage_collect
from adjustText import adjust_text
from gwaslab.textreposition import adjust_text_position
from gwaslab.annotateplot import annotate_single
from gwaslab.qqplot import _plot_qq
from gwaslab.regionalplot import _plot_regional
from gwaslab.regionalplot import process_vcf
from gwaslab.quickfix import _get_largenumber
from gwaslab.quickfix import _quick_fix_p_value
from gwaslab.quickfix import _quick_fix_pos
from gwaslab.quickfix import _quick_fix_chr
from gwaslab.quickfix import _quick_fix_eaf
from gwaslab.quickfix import _quick_fix_mlog10p
from gwaslab.quickfix import _quick_add_tchrpos
from gwaslab.quickfix import _quick_merge_sumstats
from gwaslab.quickfix import _quick_assign_i
from gwaslab.quickfix import _quick_assign_i_with_rank
from gwaslab.quickfix import _quick_extract_snp_in_region
from gwaslab.quickfix import _quick_assign_highlight_hue_pair
from gwaslab.quickfix import _quick_assign_marker_relative_size
from gwaslab.quickfix import _cut
# 20230202 ######################################################################################################

def mqqplot(insumstats,            
          chrom=None,
          pos=None,
          p=None,
          snpid=None,
          eaf=None,
          chr_dict = get_chr_to_number(),
          xtick_chr_dict = get_number_to_chr(),
          vcf_path=None,
          vcf_chr_dict = get_number_to_chr(),
          gtf_path="default",
          gtf_chr_dict = get_number_to_chr(),
          gtf_gene_name=None,
          rr_path="default",
          rr_header_dict=None,
          rr_chr_dict = get_number_to_chr(),
          mlog10p="MLOG10P",
          scaled=False,
          mode="mqq",
          scatter_kwargs=dict(),
          # region
          region = None,
          region_step = 21,
          region_grid = False,
          region_grid_line = {"linewidth": 2,"linestyle":"--"},
          region_lead_grid = True,
          region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"},
          region_hspace=0.02,
          region_ld_threshold = [0.2,0.4,0.6,0.8],
          region_ld_colors = ["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"],
          region_recombination = True,
          region_protein_coding=True,
          region_flank_factor = 0.05,
          region_anno_bbox_args=dict(),
          taf=[4,0,0.95,1,1],
          # track_n, track_n_offset,font_ratio,exon_ratio,text_offset
          tabix=None,
          mqqratio=3,
          bwindowsizekb = 100,
          density_color=False,
          density_range=(0,15),
          density_trange=(0,10),
          density_threshold=5,
          density_tpalette="Blues",
          density_palette="Reds",
          windowsizekb=500,
          anno=None,
          anno_set=list(),
          anno_alias=dict(),
          anno_d=dict(),
          anno_args=dict(),
          anno_style="right",
          anno_fixed_arm_length=None,
          anno_source = "ensembl",
          anno_adjust=False,
          anno_max_iter=100,
          arm_offset=50,
          arm_scale=1,
          arm_scale_d=None,
          cut=0,
          skip=0,
          ystep=0,
          cutfactor=10,
          cut_line_color="#ebebeb",  
          sig_line=True,
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_line=False,
          suggestive_sig_level=5e-6,
          suggestive_sig_line_color="grey",
          highlight = list(),
          highlight_color="#CB132D",
          highlight_windowkb = 500,
          pinpoint=list(),
          pinpoint_color ="red",
          stratified=False,
          maf_bins=[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)],
          maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"],
          gc=True,
          include_chrXYMT = True,
          ylim=None,
          title =None,
          mtitle=None,
          qtitle=None,
          title_pad=1.08, 
          title_fontsize=13,
          fontsize = 10,
          anno_fontsize = 10,
          figargs= dict(figsize=(15,5)),
          colors=["#597FBD","#74BAD3"],
          marker_size=(5,25),
          use_rank=False,
          verbose=True,
          repel_force=0.03,
          build="19",
          dpi=200,
          save=None,
          saveargs={"dpi":300,"facecolor":"white"},
          _invert=False,
          log=Log()
          ):

# log.writeing meta info #######################################################################################

    if verbose: log.write("Start to plot manhattan/qq plot with the following basic settings:")
    if verbose: log.write(" -Genome-wide significance level is set to "+str(sig_level)+" ...")
    if verbose: log.write(" -Raw input contains "+str(len(insumstats))+" variants...")
    if verbose: log.write(" -Plot layout mode is : "+mode)
    if len(anno_set)>0 and ("m" in mode):
        if verbose: log.write(" -Variants to annotate : "+",".join(anno_set))    
    if len(highlight)>0 and ("m" in mode):
        if verbose: log.write(" -Loci to highlight : "+",".join(highlight))    
        if verbose: log.write("  -Highlight_window is set to: ", highlight_windowkb, " kb")  
    if len(pinpoint)>0 :
        if verbose: log.write(" -Variants to pinpoint : "+",".join(pinpoint))  
    if region is not None:
        if verbose: log.write(" -Region to plot : chr"+str(region[0])+":"+str(region[1])+"-"+str(region[2])+".")  
    if "dpi" not in figargs.keys():
        figargs["dpi"] = dpi

# Plotting mode selection : layout ####################################################################
    # ax1 : manhattanplot / brisbane plot
    # ax2 : qq plot 
    # ax3 : gene track
    # ax4 : recombination rate
    # ax5 : miami plot lower panel

    # "m" : Manhattan plot
    # "qq": QQ plot
    # "r" : regional plot
    
    if  mode=="qqm": 
        fig, (ax2, ax1) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [1, mqqratio]},**figargs)
    elif mode=="mqq":
        fig, (ax1, ax2) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [mqqratio, 1]},**figargs)
    elif mode=="m":
        fig, ax1 = plt.subplots(1, 1,**figargs)
    elif mode=="qq": 
        fig, ax2 = plt.subplots(1, 1,**figargs) 
    elif mode=="r": 
        figargs["figsize"] = (15,10)
        fig, (ax1, ax3) = plt.subplots(2, 1, sharex=True, 
                                       gridspec_kw={'height_ratios': [mqqratio, 1]},**figargs)
        plt.subplots_adjust(hspace=region_hspace)
    elif mode =="b" :
        figargs["figsize"] = (15,5)
        fig, ax1 = plt.subplots(1, 1,**figargs)
        sig_level=1,
        sig_line=False,
        windowsizekb = 100000000   
        mode="mb"
        scatter_kwargs={"marker":"s"}
        marker_size= (marker_size[1],marker_size[1])
    else:
        raise ValueError("Please select one from the 5 modes: mqq/qqm/m/qq/r/b")


# Read sumstats #################################################################################################

    usecols=[]
    #  P ###############################################################################
    if mlog10p in insumstats.columns and scaled is True:
        usecols.append(mlog10p)
    elif p in insumstats.columns:
        # p value is necessary for all modes.
        usecols.append(p)  
    elif "b" in mode:
        pass
    else :
        raise ValueError("Please make sure "+p+" column is in input sumstats.")
    
    # CHR and POS ########################################################################
    # chrom and pos exists && (m || r mode)
    if (chrom is not None) and (pos is not None) and (("m" in mode) or ("r" in mode)):
        # when manhattan plot, chrom and pos is needed.
        if chrom in insumstats.columns:
            usecols.append(chrom)
        else:
            raise ValueError("Please make sure "+chrom+" column is in input sumstats.")
        if pos in insumstats.columns:
            usecols.append(pos)
        else:
            raise ValueError("Please make sure "+pos+" column is in input sumstats.")
            
    # SNPID ###############################################################################
    if len(highlight)>0 or len(pinpoint)>0 or (snpid is not None):
        # read snpid when highlight/pinpoint is needed.
        if snpid in insumstats.columns:
            usecols.append(snpid)
        else:
            raise ValueError("Please make sure "+snpid+" column is in input sumstats.")
    
    # EAF #################################################################################
    if (stratified is True) and (eaf is not None):
        # read eaf when stratified qq plot is needed.
        if eaf in insumstats.columns:
            usecols.append(eaf)
        else:
            raise ValueError("Please make sure "+eaf+" column is in input sumstats.")
    
    # ANNOTATion ##########################################################################
    #
    if len(anno_set)>0 or len(anno_alias)>0:
        if anno is None:
            anno=True
    if (anno is not None) and (anno is not True):
        if anno=="GENENAME":
            pass
        elif (anno in insumstats.columns):
            if (anno not in usecols):
                usecols.append(anno)
        else:
            raise ValueError("Please make sure "+anno+" column is in input sumstats.")
    
    if (density_color==True) or ("b" in mode and "DENSITY" in insumstats.columns):
        usecols.append("DENSITY")

    #################################################################################################
    sumstats = insumstats.loc[:,usecols].copy()


    #Standardize
    ## Annotation
    if (anno == "GENENAME"):
        anno_sig=True
    elif (anno is not None) and (anno is not True):
        sumstats["Annotation"]=sumstats.loc[:,anno].astype("string")   
      
    ## P value
    ## m, qq, r
    if "b" not in mode:   
        if scaled is True:
            sumstats["raw_P"] = pd.to_numeric(sumstats[mlog10p], errors='coerce')
        else:
            sumstats["raw_P"] = sumstats[p].astype("float64")
    
    ## CHR & POS
    ## m, qq, b
    if "m" in mode or "r" in mode or "b" in mode: 
        # convert CHR to int
        ## CHR X,Y,MT conversion ############################
        sumstats[pos] = _quick_fix_pos(sumstats[pos])
        sumstats[chrom] = _quick_fix_chr(sumstats[chrom], chr_dict=chr_dict)
    ## r
    if region is not None:
        region_chr = region[0]
        region_start = region[1]
        region_end = region[2]
        marker_size=(25,45)
        if verbose:log.write(" -Extract SNPs in region : chr"+str(region_chr)+":"+str(region[1])+"-"+str(region[2])+ "...")
        
        in_region_snp = (sumstats[chrom]==region_chr) &(sumstats[pos]<region_end) &(sumstats[pos]>region_start)
        if verbose:log.write(" -Extract SNPs in specified regions: "+str(sum(in_region_snp)))
        sumstats = sumstats.loc[in_region_snp,:]
        if len(sumstats)==0:
            log.write(" -Warning : No valid data! Please check the input.")
            return None
    
    ## EAF
    eaf_raw = pd.Series()
    if stratified is True: 
        sumstats["MAF"] = _quick_fix_eaf(sumstats[eaf])
        # for stratified qq plot
        eaf_raw = sumstats["MAF"].copy()
        
    if len(highlight)>0 and ("m" in mode):
        sumstats["HUE"] = sumstats[chrom].astype("string")
    
    if verbose: log.write("Finished loading specified columns from the sumstats.")


#sanity check############################################################################################################
    if verbose: log.write("Start conversion and sanity check:")
    if ("m" in mode or "r" in mode): 
        pre_number=len(sumstats)
        #sanity check : drop variants with na values in chr and pos df
        sumstats = sumstats.dropna(subset=[chrom,pos])
        after_number=len(sumstats)
        if verbose:log.write(" -Removed "+ str(pre_number-after_number) +" variants with nan in CHR or POS column ...")
    
    if stratified is True: 
        pre_number=len(sumstats)
        sumstats = sumstats.dropna(subset=["MAF"])
        after_number=len(sumstats)
        if verbose:log.write(" -Removed "+ str(pre_number-after_number) +" variants with nan in EAF column ...")
            
        ## Highlight
    if len(highlight)>0 and ("m" in mode or "r" in mode):
        to_highlight = sumstats.loc[sumstats[snpid].isin(highlight),:]
        #assign colors: 0 is hightlight color
        for i,row in to_highlight.iterrows():
            target_chr = int(row[chrom])
            target_pos = int(row[pos])
            right_chr=sumstats[chrom]==target_chr
            up_pos=sumstats[pos]>target_pos-highlight_windowkb*1000
            low_pos=sumstats[pos]<target_pos+highlight_windowkb*1000
            sumstats.loc[right_chr&up_pos&low_pos,"HUE"]="0"

# Density #####################################################################################################              
    if "b" in mode and "DENSITY" not in sumstats.columns:
        if verbose:log.write(" -Calculating DENSITY with windowsize of ",bwindowsizekb ," kb")
        large_number = _get_largenumber(sumstats[pos].max(),log=log)

        stack=[]
        sumstats["TCHR+POS"] = sumstats[chrom]*large_number +  sumstats[pos]
        sumstats = sumstats.sort_values(by="TCHR+POS")
        for index,row in sumstats.iterrows():
            stack.append([row["SNPID"],row["TCHR+POS"],0])  
            for i in range(2,len(stack)+1):
                if stack[-i][1]>= (row["TCHR+POS"]- 1000*bwindowsizekb):
                    stack[-i][2]+=1
                    stack[-1][2]+=1
                else:
                    break
        df = pd.DataFrame(stack,columns=["SNPID","TCHR+POS","DENSITY"])
        sumstats["DENSITY"] = df["DENSITY"].values
        bmean=sumstats["DENSITY"].mean()
        bmedian=sumstats["DENSITY"].median()
    elif "b" in mode and "DENSITY" in sumstats.columns:
        bmean=sumstats["DENSITY"].mean()
        bmedian=sumstats["DENSITY"].median()
        if verbose:log.write(" -DENSITY column exists. Skipping calculation...")
     
    #############################
# P value conversion #####################################################################################################  
    ## m,qq,r -> dropna
    if "b" not in mode:
        pre_number=len(sumstats)
        sumstats = sumstats.dropna(subset=["raw_P"])
        after_number=len(sumstats)
        if verbose:log.write(" -Removed "+ str(pre_number-after_number) +" variants with nan in P column ...")
    
    ## b: value to plot is density
    if "b" in mode:
        sumstats["scaled_P"] = sumstats["DENSITY"].copy()
        sumstats["raw_P"] = -np.log10(sumstats["DENSITY"].copy()+2)
    elif scaled is True:
        if verbose:log.write(" -P values are already converted to -log10(P)!")
        sumstats["scaled_P"] = sumstats["raw_P"].copy()
        sumstats["raw_P"] = np.power(10,-sumstats["scaled_P"].astype("float64"))
    else:
        if not scaled:
            # quick fix p
            sumstats = _quick_fix_p_value(sumstats, verbose=verbose, log=log)

        # quick fix mlog10p
        sumstats = _quick_fix_mlog10p(sumstats, scaled=scaled, verbose=verbose, log=log)

    # raw p for calculate lambda
    p_toplot_raw = sumstats[["CHR","scaled_P"]].copy()
    
    # filter out variants with -log10p < skip
    sumstats = sumstats.loc[sumstats["scaled_P"]>=skip,:]
    garbage_collect.collect()
    
    # shrink variants above cut line #########################################################################################
    try:
        sumstats["scaled_P"], maxy, maxticker, cut, cutfactor = _cut(series = sumstats["scaled_P"], 
                                                                        mode =mode, 
                                                                        cut=cut,
                                                                        cutfactor = cutfactor,
                                                                        verbose =verbose, 
                                                                        log = log)
    except:
        log.write(" -Warning : No valid data! Please check the input.")
        return None
    
    #maxy = sumstats["scaled_P"].max()
    #if "b" not in mode:
    #    if verbose: log.write(" -Maximum -log10(P) values is "+str(maxy) +" .")
    #elif "b" in mode:
    #    if verbose: log.write(" -Maximum DENSITY values is "+str(maxy) +" .")
#
    #maxticker=int(np.round(sumstats["scaled_P"].max(skipna=True)))
    #if cut:
    #    if cut is True:
    #        if verbose: log.write(" -Cut Auto mode is activated...")
    #        if maxy<20:
    #            if verbose: log.write(" - maxy <20 , no need to cut.")
    #            cut=0
    #        else:
    #            cut = 20
    #            cutfactor = ( maxy - cut )/5
    #    if cut:
    #        if "b" not in mode:
    #            if verbose: log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")
    #        else:
    #            if verbose: log.write(" -Minus DENSITY values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")
#
    #        maxticker=int(np.round(sumstats["scaled_P"].max(skipna=True)))
#
    #        sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
    #        maxy = (maxticker-cut)/cutfactor + cut
    #if verbose: log.write("Finished data conversion and sanity check.")


    # Manhattan plot ##########################################################################################################
    ## regional plot ->rsq
        #calculate rsq]
    if vcf_path is not None:
        if tabix is None:
            tabix = which("tabix")
        sumstats = process_vcf(sumstats=sumstats, vcf_path=vcf_path,region=region, 
                               log=log ,pos=pos,region_ld_threshold=region_ld_threshold,verbose=verbose,vcf_chr_dict=vcf_chr_dict,tabix=tabix)

    #sort & add id
    ## Manhatann plot ###################################################
    if ("m" in mode) or ("r" in mode): 
        # assign index i and tick position
        sumstats,chrom_df=_quick_assign_i_with_rank(sumstats, use_rank=use_rank, chrom="CHR",pos="POS")
        
        ## Assign marker size ##############################################
        sumstats["s"]=1
        if "b" not in mode:
            sumstats.loc[sumstats["scaled_P"]>-np.log10(5e-4),"s"]=2
            sumstats.loc[sumstats["scaled_P"]>-np.log10(suggestive_sig_level),"s"]=3
            sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level),"s"]=4
        sumstats["chr_hue"]=sumstats[chrom].astype("string")

        if vcf_path is not None:
            sumstats["chr_hue"]=sumstats["LD"]

        if verbose:log.write("Start to create manhattan plot with "+str(len(sumstats))+" variants:")
        ## default seetings
        
        palette = sns.color_palette(colors,n_colors=sumstats[chrom].nunique())  
        

        legend = None
        style=None
        linewidth=0
        # if regional plot assign colors
        if vcf_path is not None:
            palette = { i:region_ld_colors[i] for i in range(len(region_ld_colors))}
            legend=None
            style=None
            linewidth=1
            
        ## if highlight 
        highlight_i = pd.DataFrame()
        if len(highlight) >0:
            plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P',
                               hue='chr_hue',
                               palette=palette,
                               legend=legend,
                               style=style,
                               size="s",
                               sizes=marker_size,
                               linewidth=linewidth,
                               zorder=2,ax=ax1,edgecolor="black", **scatter_kwargs)   
            
            if verbose: log.write(" -Highlighting target loci...")
            sns.scatterplot(data=sumstats.loc[sumstats["HUE"]=="0"], x='i', y='scaled_P',
                   hue="HUE",
                   palette={"0":highlight_color},
                   legend=legend,
                   style=style,
                   size="s",
                   sizes=(marker_size[0]+1,marker_size[1]+1),
                   linewidth=linewidth,
                   zorder=3,ax=ax1,edgecolor="black",**scatter_kwargs)  
            highlight_i = sumstats.loc[sumstats[snpid].isin(highlight),"i"].values
        
        ## if not highlight    
        else:
            if density_color == True:
                hue = "DENSITY_hue"
                s = "DENSITY"
                to_plot = sumstats.sort_values("DENSITY")
                to_plot["DENSITY_hue"] = to_plot["DENSITY"].astype("float")
                plot = sns.scatterplot(data=to_plot.loc[to_plot["DENSITY"]<=density_threshold,:], x='i', y='scaled_P',
                       hue=hue,
                       palette= density_tpalette,
                       legend=legend,
                       style=style,
                       size=s,
                       sizes=(marker_size[0]+1,marker_size[0]+1),
                       linewidth=linewidth,
                       hue_norm=density_trange,
                       zorder=2,ax=ax1,edgecolor="black",**scatter_kwargs) 

                plot = sns.scatterplot(data=to_plot.loc[to_plot["DENSITY"]>density_threshold,:], x='i', y='scaled_P',
                   hue=hue,
                   palette= density_palette,
                   legend=legend,
                   style=style,
                   size=s,
                   sizes=marker_size,
                   hue_norm=density_range,
                   linewidth=linewidth,
                   zorder=2,ax=ax1,edgecolor="black",**scatter_kwargs)   
            else:
                s = "s"
                hue = 'chr_hue'
                hue_norm=None
                to_plot = sumstats
                plot = sns.scatterplot(data=to_plot, x='i', y='scaled_P',
                       hue=hue,
                       palette= palette,
                       legend=legend,
                       style=style,
                       size=s,
                       sizes=marker_size,
                       hue_norm=hue_norm,
                       linewidth=linewidth,
                       zorder=2,ax=ax1,edgecolor="black",**scatter_kwargs)   
        
        
        ## if pinpoint variants
        if (len(pinpoint)>0):
            if sum(sumstats[snpid].isin(pinpoint))>0:
                to_pinpoint = sumstats.loc[sumstats[snpid].isin(pinpoint),:]
                if verbose: log.write(" -Pinpointing target vairants...")
                ax1.scatter(to_pinpoint["i"],to_pinpoint["scaled_P"],color=pinpoint_color,zorder=3,s=marker_size[1]+1)
            else:
                if verbose: log.write(" -Target vairants to pinpoint were not found. Skip pinpointing process...")
        
        #ax1.set_xticks(chrom_df.astype("float64"))
        #ax1.set_xticklabels(chrom_df.index.astype("Int64").map(xtick_chr_dict),fontsize=fontsize,family="sans-serif")
        
        # if regional plot : pinpoint lead , add color bar ##################################################
        if (region is not None) and ("r" in mode):
            ax1, ax3 =_plot_regional(
                                sumstats=sumstats,
                                fig=fig,
                                ax1=ax1,
                                ax3=ax3,
                                region=region,
                                vcf_path=vcf_path,
                                marker_size=marker_size,
                                fontsize=fontsize,
                                build=build,
                                chrom_df=chrom_df,
                                xtick_chr_dict=xtick_chr_dict,
                                cut_line_color=cut_line_color,
                                vcf_chr_dict =vcf_chr_dict,
                                gtf_path=gtf_path,
                                gtf_chr_dict = gtf_chr_dict,
                                gtf_gene_name=gtf_gene_name,
                                rr_path=rr_path,
                                rr_header_dict=rr_header_dict,
                                rr_chr_dict = rr_chr_dict,
                                mode=mode,
                                region_step = region_step,
                                region_grid = region_grid,
                                region_grid_line = region_grid_line,
                                region_lead_grid = region_lead_grid,
                                region_lead_grid_line = region_lead_grid_line,
                                region_hspace=region_hspace,
                                region_ld_threshold = region_ld_threshold,
                                region_ld_colors = region_ld_colors,
                                region_recombination = region_recombination,
                                region_protein_coding=region_protein_coding,
                                region_flank_factor =region_flank_factor,
                                taf=taf,
                                tabix=tabix,
                                chrom=chrom,
                                pos=pos,
                                verbose=verbose,
                                log=log
                            )

        if region is None:
            #plot.set_xlabel(chrom)
            ax1.set_xticks(chrom_df.astype("float64"))
            ax1.set_xticklabels(chrom_df.index.astype("Int64").map(xtick_chr_dict),fontsize=fontsize,family="sans-serif")
        
        # genomewide significant line
        if sig_line is True:
            sigline = ax1.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
        if suggestive_sig_line is True:
            suggestive_sig_line = ax1.axhline(y=-np.log10(suggestive_sig_level), linewidth = 2,linestyle="--",color=suggestive_sig_line_color,zorder=1)
        
        # for brisbane plot, add median and mean line
        if "b" in mode:    
            meanline = ax1.axhline(y=bmean, linewidth = 2,linestyle="-",color=sig_line_color,zorder=1000)
            medianline = ax1.axhline(y=bmedian, linewidth = 2,linestyle="--",color=sig_line_color,zorder=1000)
        
        # cut line
        if cut == 0: 
            ax1.set_ylim(skip, ceil(maxy*1.2) )
        
        if cut:
            cutline = ax1.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            step=2
            if ((maxticker-cut)/cutfactor + cut) > cut:
                if ystep == 0:
                    if (cut - skip ) // step > 10:
                        step = (cut - skip ) // 10
                else:
                    step = ystep

                ax1.set_yticks([x for x in range(skip,cut-step,step)]+[cut]+[(maxticker-cut)/cutfactor + cut])
                ax1.set_yticklabels([x for x in range(skip,cut-step,step)]+[cut]+[maxticker],fontsize=fontsize,family="sans-serif")
                #ax1.set_yticks([x for x in range(skip,cut+1,step)]+[(maxticker-cut)/cutfactor + cut])
                #ax1.set_yticklabels([x for x in range(skip,cut+1,step)]+[maxticker],fontsize=fontsize,family="sans-serif")
            else:
                if ystep == 0:
                    if (cut - skip ) // step > 10:
                        step = (cut - skip ) // 10
                else:
                    step = ystep

                ax1.set_yticks([x for x in range(skip,cut-step,step)]+[cut])
                ax1.set_yticklabels([x for x in range(skip,cut-step,step)]+[cut],fontsize=fontsize,family="sans-serif")
                #ax1.set_yticks([x for x in range(skip,cut+1,step)])
                #ax1.set_yticklabels([x for x in range(skip,cut+1,step)],fontsize=fontsize,family="sans-serif")
            ax1.set_ylim(bottom = skip)


        # Get top variants for annotation #######################################################
        if (anno and anno!=True) or (len(anno_set)>0):
            if len(anno_set)>0:
                to_annotate=sumstats.loc[sumstats[snpid].isin(anno_set),:]
                if to_annotate.empty is not True:
                    if verbose: log.write(" -Found "+str(len(to_annotate))+" specified variants to annotate...")
            else:
                to_annotate=getsig(sumstats.loc[sumstats["scaled_P"]> float(-np.log10(sig_level)),:],
                               snpid,
                               chrom,
                               pos,
                               "raw_P",
                               sig_level=sig_level,
                               windowsizekb=windowsizekb,
                               verbose=False)
                if (to_annotate.empty is not True) and ("b" not in mode):
                    if verbose: log.write(" -Found "+str(len(to_annotate))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")
        else:
            to_annotate=getsig(sumstats.loc[sumstats["scaled_P"]> float(-np.log10(sig_level)),:],
                               "i",
                               chrom,
                               pos,
                               "raw_P",
                               windowsizekb=windowsizekb,
                               verbose=False,
                               sig_level=sig_level)
            if (to_annotate.empty is not True) and ("b" not in mode):
                if verbose: log.write(" -Found "+str(len(to_annotate))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")
        if (to_annotate.empty is not True) and anno=="GENENAME":
            to_annotate = annogene(to_annotate,
                                   id=snpid,
                                   chrom=chrom,
                                   pos=pos,
                                   log=log,
                                   build=build,
                                   source=anno_source,
                                   verbose=verbose).rename(columns={"GENE":"Annotation"})

        # Add Annotation to manhattan plot #######################################################
        
        if "b" in mode:
            ax1.set_ylabel("Density of GWAS \n SNPs within "+str(bwindowsizekb)+" kb",ha="center",va="bottom",fontsize=fontsize,family="sans-serif")
        else:
            ax1.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family="sans-serif")

        if region is not None:
            if (gtf_path is not None ) and ("r" in mode):
                ax3.set_xlabel("Chromosome "+str(region[0])+" (MB)",fontsize=fontsize,family="sans-serif")
            else:
                ax1.set_xlabel("Chromosome "+str(region[0])+" (MB)",fontsize=fontsize,family="sans-serif")
        else:
            ax1.set_xlabel("Chromosomes",fontsize=fontsize,family="sans-serif")
        ##
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.spines["left"].set_visible(True)
        if mode=="r":
            ax1.spines["top"].set_visible(True)
            ax1.spines["top"].set_zorder(1)    
            ax1.spines["right"].set_visible(True)
        
        if verbose: log.write("Finished creating Manhattan plot successfully!")
        if mtitle and anno and len(to_annotate)>0: 
            pad=(ax1.transData.transform((skip, title_pad*maxy))[1]-ax1.transData.transform((skip, maxy)))[1]
            ax1.set_title(mtitle,pad=pad,fontsize=title_fontsize,family="sans-serif")
        elif mtitle:
            ax1.set_title(mtitle,fontsize=title_fontsize,family="sans-serif")
            
        # add annotation arrows and texts
        ax1 = annotate_single(
                                sumstats=sumstats,
                                anno=anno,
                                mode=mode,
                                ax1=ax1,
                                highlight_i=highlight_i,
                                to_annotate=to_annotate,
                                anno_d=anno_d,
                                anno_alias=anno_alias,
                                anno_style=anno_style,
                                anno_args=anno_args,
                                arm_scale=arm_scale,
                                anno_max_iter=anno_max_iter,
                                arm_scale_d=arm_scale_d,
                                arm_offset=arm_offset,
                                anno_adjust=anno_adjust,
                                anno_fixed_arm_length=anno_fixed_arm_length,
                                maxy=maxy,
                                anno_fontsize= anno_fontsize,
                                region=region,
                                region_anno_bbox_args=region_anno_bbox_args,
                                skip=skip,
                                snpid=snpid,
                                chrom=chrom,
                                pos=pos,
                                repel_force=repel_force,
                                verbose=verbose,
                                log=log
                            )  
    # Creating Manhatann plot Finished #####################################################################

    # QQ plot #########################################################################################################
    if "qq" in mode:
        # ax2 qqplot
        ax2 =_plot_qq(
                    sumstats=sumstats,
                    p_toplot_raw=p_toplot_raw,
                    ax2=ax2,
                    maxticker=maxticker,
                    gc=gc,
                    cut=cut,
                    cutfactor=cutfactor,
                    skip=skip,
                    maxy=maxy,
                    colors=colors,
                    sig_line_color=sig_line_color,
                    stratified=stratified,
                    eaf_raw=eaf_raw,
                    maf_bins=maf_bins,
                    maf_bin_colors=maf_bin_colors,
                    fontsize=fontsize,
                    qtitle=qtitle,
                    title_fontsize=title_fontsize,
                    include_chrXYMT=include_chrXYMT,
                    cut_line_color=cut_line_color,
                    verbose=verbose,
                    log=log
                )
    
    if ylim is not None:
        ax1.set_ylim(ylim)
            
    # Saving plot ##########################################################################################################
    if save:
        if verbose: log.write("Saving plot:")
        if save==True:
            fig.savefig("./"+mode+"_plot.png",bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ "./"+mode+"_plot.png" + " successfully!" )
        else:
            fig.savefig(save,bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ save + " successfully!" )
    
    # add title 
    if title and anno and len(to_annotate)>0:
        # increase height if annotation 
        fig.suptitle(title , fontsize = title_fontsize ,x=0.5, y=1.05)
    else:
        fig.suptitle(title , fontsize = title_fontsize, x=0.5,y=1)
    
    garbage_collect.collect()
    
    # Return matplotlib figure object #######################################################################################
    return fig, log
