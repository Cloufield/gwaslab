import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
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
import gc as garbage_collect
from adjustText import adjust_text
from gwaslab.textreposition import adjust_text_position

# 20220310 ######################################################################################################

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
          scatter_kwargs={},
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
          region_anno_bbox_args={},
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
          anno_set=[],
          anno_alias={},
          anno_d={},
          anno_args={},
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
          cutfactor=10,
          cut_line_color="#ebebeb",  
          sig_line=True,
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_level=5e-6,
          highlight = [],
          highlight_color="#CB132D",
          highlight_windowkb = 500,
          pinpoint=[],
          pinpoint_color ="red",
          stratified=False,
          maf_bins=[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)],
          maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"],
          gc=True,
          include_chrXYMT = True,
          title =None,
          mtitle=None,
          qtitle=None,
          title_pad=1.08, 
          title_fontsize=13,
          fontsize = 10,
          anno_fontsize = 10,
          figargs= {"figsize":(15,5)},
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
        if verbose: log.write(" -Highlight_window is set to: ", highlight_windowkb, " kb")  
    if len(pinpoint)>0 :
        if verbose: log.write(" -Variants to pinpoint : "+",".join(pinpoint))  
    if region is not None:
        if verbose: log.write(" -Region to plot : chr"+str(region[0])+":"+str(region[1])+"-"+str(region[2])+".")  
    if "dpi" not in figargs.keys():
        figargs["dpi"] = dpi


# Plotting mode selection : layout ####################################################################
    # ax1 : manhattanplot : 
    # ax2 : qq plot 
    # ax3 : gene track
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
    if p in insumstats.columns:
        # p value is necessary for all modes.
        usecols.append(p)
    elif mlog10p in insumstats.columns and scaled is True:
        usecols.append(mlog10p)
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
        if pd.api.types.is_string_dtype(sumstats[chrom]) == True:
        ## if chr is string dtype: convert using chr_dict
            sumstats[chrom] = sumstats[chrom].map(chr_dict,na_action="ignore")
        ## CHR
        sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
        ## POS
        sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    
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
    
    ## EAF
    if stratified is True: 
        sumstats["MAF"] = pd.to_numeric(sumstats[eaf], errors='coerce')
        sumstats.loc[sumstats["MAF"]>0.5,"MAF"] = 1 - sumstats.loc[sumstats["MAF"]>0.5,"MAF"]
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
        
        large_number = 1000000000
        for i in range(7):
            if sumstats[pos].max()*10 >  large_number:
                large_number = int(large_number * 10)
            else:
                break
        
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
        if verbose:log.write(" -P values are being converted to -log10(P)...") 
        bad_p_value=len(sumstats.loc[(sumstats["raw_P"]>1)|(sumstats["raw_P"]<=0),:])
        if verbose:log.write(" -Sanity check after conversion: "+ str(bad_p_value) +" variants with P value outside of (0,1] will be removed...")
        sumstats = sumstats.loc[(sumstats["raw_P"]<=1)&(sumstats["raw_P"]>0),:]
        sumstats["scaled_P"] = -np.log10(sumstats["raw_P"])
        is_inf = sumstats["scaled_P"].isin([np.inf, -np.inf, float('inf'),-float('inf')])
        is_na = sumstats["scaled_P"].isna()
        bad_p = sum(is_inf | is_na)
        if verbose: log.write(" -Sanity check: "+str(bad_p) + " na/inf/-inf variants will be removed..." )
          
        sumstats = sumstats.loc[~(is_inf | is_na),:]

    # raw p for calculate lambda
    p_toplot_raw = sumstats[["CHR","scaled_P"]].copy()
    
    # filter out variants with -log10p < skip
    sumstats = sumstats.loc[sumstats["scaled_P"]>=skip,:]
    garbage_collect.collect()
    
    # shrink variants above cut line #########################################################################################
    maxy = sumstats["scaled_P"].max()
    if "b" not in mode:
        if verbose: log.write(" -Maximum -log10(P) values is "+str(maxy) +" .")
    elif "b" in mode:
        if verbose: log.write(" -Maximum DENSITY values is "+str(maxy) +" .")
            
    if cut:
        if cut is True:
            if verbose: log.write(" -Cut Auto mode is activated...")
            if maxy<20:
                if verbose: log.write(" - maxy <20 , no need to cut.")
                cut=0
            else:
                cut = 20
                cutfactor = ( maxy - cut )/5
        if cut:
            if "b" not in mode:
                if verbose: log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")
            else:
                if verbose: log.write(" -Minus DENSITY values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")

            maxticker=int(np.round(sumstats["scaled_P"].max(skipna=True)))

            sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
            maxy = (maxticker-cut)/cutfactor + cut
    if verbose: log.write("Finished data conversion and sanity check.")


# Manhattan plot ##########################################################################################################
## regional plot ->rsq
     #calculate rsq]
    if vcf_path is not None:
        if tabix is None:
            tabix = which("tabix")
        sumstats = process_vcf(sumstats=sumstats, vcf_path=vcf_path,region=region, 
                               log=log ,pos=pos,region_ld_threshold=region_ld_threshold,verbose=verbose,vcf_chr_dict=vcf_chr_dict,tabix=tabix)

# # Create index for plotting ###################################################


    #sort & add id
    if ("m" in mode) or ("r" in mode): 
        sumstats = sumstats.sort_values([chrom,pos])
        if use_rank is True: sumstats["POS_RANK"] = sumstats.groupby(chrom)[pos].rank("dense", ascending=True)
        sumstats["id"]=range(len(sumstats))
        sumstats=sumstats.set_index("id")

        #create a df , groupby by chromosomes , and get the maximum position
        if use_rank is True: 
            posdic = sumstats.groupby(chrom)["POS_RANK"].max()
        else:
            posdic = sumstats.groupby(chrom)[pos].max()
        
        # convert to dictionary
        posdiccul = dict(posdic)
        
        # fill empty chr with 0
        for i in range(0,26):
            if i in posdiccul: continue
            else: posdiccul[i]=0
        
        # cumulative sum dictionary
        for i in range(2,sumstats[chrom].max()+1):
            posdiccul[i]= posdiccul[i-1] + posdiccul[i] + sumstats[pos].max()*0.05

        # convert base pair postion to x axis position using the cumulative sum dictionary
        sumstats["add"]=sumstats[chrom].apply(lambda x : posdiccul[int(x)-1])
        
        if use_rank is True: 
            sumstats["i"]=sumstats["POS_RANK"]+sumstats["add"]
        else:
            sumstats["i"]=sumstats[pos]+sumstats["add"]
        

        #for plot, get the chr text tick position      
        chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
        #sumstats["i"] = sumstats["i"]+((sumstats[chrom].map(dict(chrom_df)).astype("int")))*0.02
        #sumstats["i"] = sumstats["i"].astype("Int64")
        sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
        
        
        ## Assign marker size ##############################################
        
        sumstats["s"]=1
        if "b" not in mode:
            sumstats.loc[sumstats["scaled_P"]>-np.log10(5e-4),"s"]=2
            sumstats.loc[sumstats["scaled_P"]>-np.log10(suggestive_sig_level),"s"]=3
            sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level),"s"]=4
        sumstats["chr_hue"]=sumstats[chrom].astype("string")

        if vcf_path is not None:
            sumstats["chr_hue"]=sumstats["LD"]
## Manhatann plot ###################################################
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
                hue = "DENSITY"
                s = "DENSITY"
                to_plot = sumstats.sort_values("DENSITY")

                plot = sns.scatterplot(data=to_plot.loc[to_plot["DENSITY"]<=density_threshold], x='i', y='scaled_P',
                       hue=hue,
                       palette= density_tpalette,
                       legend=legend,
                       style=style,
                       size=s,
                       sizes=(marker_size[0]+1,marker_size[0]+1),
                       linewidth=linewidth,
                       hue_norm=density_trange,
                       zorder=2,ax=ax1,edgecolor="black",**scatter_kwargs) 
                plot = sns.scatterplot(data=to_plot.loc[to_plot["DENSITY"]>density_threshold], x='i', y='scaled_P',
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
        
        # if regional plot : pinpoint lead , add color bar ##################################################
        if (region is not None) :
            # pinpoint lead
            lead_id = sumstats["scaled_P"].idxmax()
            ax1.scatter(sumstats.loc[lead_id,"i"],sumstats.loc[lead_id,"scaled_P"],
                        color=region_ld_colors[-1],marker="D",zorder=3,s=marker_size[1]+1,edgecolor="black")
            if (vcf_path is not None):
                # add a in-axis axis for colorbar
                axins1 = inset_axes(ax1,
                        width="5%",  # width = 50% of parent_bbox width
                        height="40%",  # height : 5%
                        loc='upper right',axes_kwargs={"frameon":True,"facecolor":"white","zorder":999999})

                cmp= mpl.colors.ListedColormap(region_ld_colors[1:])
                norm= mpl.colors.BoundaryNorm([0]+region_ld_threshold+[1],cmp.N)     
                cbar= plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmp),
                                   cax=axins1,
                                   fraction=1,
                                   orientation="vertical",
                                   ticklocation="left")
                cbar.set_ticks(ticks=region_ld_threshold) 
                cbar.set_ticklabels([str(i) for i in region_ld_threshold]) 
                cbar.ax.set_title('LD $r^{2}$',loc="center",y=-0.2)
                rect = patches.Rectangle((0.91,0.50),
                                        height=0.5,
                                        width=0.09, 
                                        transform=ax1.transAxes,
                                        linewidth=1, 
                                        edgecolor='black', 
                                        facecolor='white',
                                        zorder=999998)
                ax1.add_patch(rect)
        ## recombinnation rate ##################################################       
        if (region is not None) and (region_recombination is True):
            ax4=ax1.twinx()
            most_left_snp = sumstats["i"].idxmin()
            rc_track_offset = sumstats.loc[most_left_snp,"i"]-sumstats.loc[most_left_snp,pos]
            if rr_path=="default":
                if rr_chr_dict is not None:
                    rr_chr = rr_chr_dict[region[0]]
                else:
                    rr_chr = str(region[0])
                rc = get_recombination_rate(chrom=rr_chr,build=build)
            else:
                rc = pd.read_csv(rr_path,sep="\t")
                if rr_header_dict is not None:
                    rc = rc.rename(columns=rr_header_dict)

            rc = rc.loc[(rc["Position(bp)"]<region[2]) & (rc["Position(bp)"]>region[1]),:]
            ax4.plot(rc_track_offset+rc["Position(bp)"],rc["Rate(cM/Mb)"],color="#5858FF",zorder=1)
            ax4.set_ylabel("Recombination rate(cM/Mb)")
            ax4.set_ylim(0,100)
            ax4.spines["top"].set_visible(False)
            ax4.spines["top"].set(zorder=1)    
        
        ## regional plot : gene track ######################################################################
                    # calculate offset
        if (region is not None):
            most_left_snp = sumstats["i"].idxmin()
            gene_track_offset = sumstats.loc[most_left_snp,pos]-region[1]
            gene_track_start_i = sumstats.loc[most_left_snp,"i"] - gene_track_offset - region[1]
            lead_id=sumstats["scaled_P"].idxmax()
            lead_snp_y = sumstats["scaled_P"].max()
            lead_snp_i = sumstats.loc[lead_id,"i"]
            
        if (gtf_path is not None ) and ("r" in mode):
            # load gtf
            if verbose: log.write(" -Loading gtf files from:" + gtf_path)
            uniq_gene_region,exons = process_gtf(gtf_path = gtf_path ,
                                                    region = region, 
                                                    region_flank_factor = region_flank_factor,
                                                    build=build,
                                                    region_protein_coding=region_protein_coding,
                                                    gtf_chr_dict=gtf_chr_dict,
                                                    gtf_gene_name=gtf_gene_name)
            n_uniq_stack = uniq_gene_region["stack"].nunique()
            stack_num_to_plot = max(taf[0],n_uniq_stack)
            ax3.set_ylim((-stack_num_to_plot*2-taf[1]*2,2+taf[1]*2))
            ax3.set_yticks([])
            pixels_per_point = 72/fig.dpi
            pixels_per_track = np.abs(ax3.transData.transform([0,0])[1] - ax3.transData.transform([0,1])[1])                   
            font_size_in_pixels= taf[2] * pixels_per_track
            font_size_in_points =  font_size_in_pixels * pixels_per_point
            linewidth_in_points=   pixels_per_track * pixels_per_point
            if verbose: log.write(" -plotting gene track..")
            
            sig_gene_name = "Undefined"
            texts_to_adjust_left = []
            texts_to_adjust_middle = []
            texts_to_adjust_right = []
            for index,row in uniq_gene_region.iterrows():
                if row[6][0]=="+":
                    gene_anno = row["name"] + "->"
                else:
                    gene_anno = "<-" + row["name"] 
                
                if region_lead_grid is True and lead_snp_i > gene_track_start_i+row[3] and lead_snp_i < gene_track_start_i+row[4] :
                        gene_color="#FF0000"
                        sig_gene_name = row["name"]
                        sig_gene_left = gene_track_start_i+row[3]
                        sig_gene_right= gene_track_start_i+row[4]
                else:
                    gene_color="#020080"
                
                # plot gene line
                ax3.plot((gene_track_start_i+row[3],gene_track_start_i+row[4]),
                         (row["stack"]*2,row["stack"]*2),color=gene_color,linewidth=linewidth_in_points/10)
                
                # plot gene name
                if row[4] >= region[2]:
                    #right side
                    texts_to_adjust_right.append(ax3.text(x=gene_track_start_i+region[2],
                         y=row["stack"]*2+taf[4],s=gene_anno,ha="right",va="center",color="black",style='italic', size=font_size_in_points))
                    
                elif row[3] <= region[1] :
                    #left side
                    texts_to_adjust_left.append(ax3.text(x=gene_track_start_i+region[1],
                         y=row["stack"]*2+taf[4],s=gene_anno,ha="left",va="center",color="black",style='italic', size=font_size_in_points))
                else:
                    texts_to_adjust_middle.append(ax3.text(x=(gene_track_start_i+row[3]+gene_track_start_i+row[4])/2,
                         y=row["stack"]*2+taf[4],s=gene_anno,ha="center",va="center",color="black",style='italic',size=font_size_in_points))
            
            # plot exons
            for index,row in exons.iterrows():
                if not pd.isnull(row["name"]):
                    if (region_lead_grid is True) and row["name"]==sig_gene_name:
                        exon_color = region_lead_grid_line["color"]  
                    else:
                        exon_color="#020080" 
                elif gene_track_start_i+row[3] > sig_gene_left and gene_track_start_i+row[4] < sig_gene_right:
                    exon_color = region_lead_grid_line["color"]  
                else:
                    exon_color="#020080"
                #if row["gene_biotype"]!="protein_coding":
                #    exon_color="grey"
                ax3.plot((gene_track_start_i+row[3],gene_track_start_i+row[4]),
                         (row["stack"]*2,row["stack"]*2),linewidth=linewidth_in_points*taf[3],color=exon_color,solid_capstyle="butt")

            if verbose: log.write(" -Finished plotting gene track..")
                           
        #plot.set_xlabel(chrom); 
        ax1.set_xticks(chrom_df.astype("float64"))

        ax1.set_xticklabels(chrom_df.index.astype("Int64").map(xtick_chr_dict),fontsize=fontsize,family="sans-serif")
        
        ## regional plot - set X tick
        if region is not None:
            region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 
            
            # set x ticks for gene track
            if (gtf_path is not None ) and ("r" in mode):
                ax3.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
                ax3.set_xticklabels(region_ticks,rotation=45,fontsize=fontsize,family="sans-serif")
                if region_grid is True:
                    for i in np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step):
                        ax1.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
                        ax3.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
                if region_lead_grid is True:
                    ax1.axvline(x=lead_snp_i,ymax=1/1.2, zorder=1,**region_lead_grid_line)
                    ax3.axvline(x=lead_snp_i, zorder=1,**region_lead_grid_line)
            else:
                # set x ticks m plot
                ax1.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
                ax1.set_xticklabels(region_ticks,rotation=45,fontsize=fontsize,family="sans-serif")
        
            ax1.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])
        # genomewide significant line
        if sig_line is True:
            sigline = ax1.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
        if "b" in mode:
            meanline = ax1.axhline(y=bmean, linewidth = 2,linestyle="-",color=sig_line_color,zorder=1000)
            medianline = ax1.axhline(y=bmedian, linewidth = 2,linestyle="--",color=sig_line_color,zorder=1000)
        # cut line
        if cut == 0: ax1.set_ylim(skip,maxy*1.2)
        if cut:
            cutline = ax1.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax1.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax1.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker],fontsize=fontsize,family="sans-serif")
            else:
                ax1.set_yticks([x for x in range(skip,cut+1,2)])
                ax1.set_yticklabels([x for x in range(skip,cut+1,2)],fontsize=fontsize,family="sans-serif")
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
            #if build=="19":
            #    data = EnsemblRelease(75)
            #    if verbose:log.write(" -Assigning Gene name using Ensembl Release",75)
            #elif build=="38":
            #    data = EnsemblRelease(77)
            #    if verbose:log.write(" -Assigning Gene name using Ensembl Release",77)
            #to_annotate = to_annotate.copy()
            #to_annotate.loc[:,["LOCATION","Annotation"]] = pd.DataFrame(to_annotate.apply(lambda x:closest_gene(x,data=data), axis=1).tolist(), index=to_annotate.index).values

# Add Annotation to manhattan plot #######################################################
           ## final
        
        ax1.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
        if "b" in mode:
            ax1.set_ylabel("Density of GWAS \n SNPs within "+str(bwindowsizekb)+" kb",ha="center",va="bottom",fontsize=fontsize,family="sans-serif")
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
        if anno and (to_annotate.empty is not True):
            #initiate a list for text and a starting position
            text = []
            last_pos=0
            anno_count=0
            to_annotate = to_annotate.sort_values(by=[chrom,pos])
            ## log : annotation column
            if anno==True:
                    annotation_col="CHR:POS"
            elif anno:
                    annotation_col=anno
            if verbose: log.write(" -Annotating using column "+annotation_col+"...")
                            ## calculate y span
            if region is not None:
                y_span = region[2] - region[1]
            else:
                y_span = sumstats["i"].max()-sumstats["i"].min()
            
            if verbose: log.write(" -Adjusting text positions with repel_force={}...".format(repel_force))
            if anno_style == "expand" :
                to_annotate.loc[:, "ADJUSTED_i"] = adjust_text_position(to_annotate["i"].values.copy(), y_span, repel_force,max_iter=anno_max_iter,log=log,verbose=verbose)
            ##  iterate through variants to be annotated
            anno_to_adjust_list = list()
            for rowi,row in to_annotate.iterrows():
                
                # avoid text overlapping
                ## adjust x to avoid overlapping################################################################
                if anno_style == "right" :
                    #right style
                    if row["i"]>last_pos+repel_force*y_span:
                        last_pos=row["i"]
                    else:
                        last_pos+=repel_force*y_span
                elif anno_style == "expand" :
                    #expand style
                    last_pos = row["ADJUSTED_i"]
                    anno_args["rotation"] = 90
                elif anno_style == "tight" :
                    #tight style
                    anno_fixed_arm_length = 1
                    anno_adjust = True
                    anno_args["rotation"] = 90
                else:
                    pass
                ################################################################
                #shrink or increase the arm
                if arm_scale_d is not None:
                    if anno_count not in arm_scale_d.keys():
                        arm_scale =1
                    else:
                        arm_scale = arm_scale_d[anno_count]
                ################################################################

                # vertical arm length in pixels
                armB_length_in_point = ax1.transData.transform((skip,1.15*maxy))[1]-ax1.transData.transform((skip, row["scaled_P"]+1))[1]-arm_offset/2
                # scale if needed
                armB_length_in_point = armB_length_in_point*arm_scale
                ################################################################
                if arm_scale>=1:
                    armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else ax1.transData.transform((skip, maxy+2))[1]-plot.transData.transform((skip,  row["scaled_P"]+1))[1] 
                ###if anno_fixed_arm_length #############################################################
                if anno_fixed_arm_length is not None:
                    anno_fixed_arm_length_factor = ax1.transData.transform((skip,anno_fixed_arm_length))[1]-ax1.transData.transform((skip,0))[1] 
                    armB_length_in_point = anno_fixed_arm_length_factor
                ################################################################################################################################
                # annotation alias
                if anno==True:
                    if row[snpid] in anno_alias.keys():
                        annotation_text = anno_alias[row[snpid]]
                    else:
                        annotation_text="Chr"+ str(row[chrom]) +":"+ str(int(row[pos]))
                elif anno:
                    annotation_text=row["Annotation"]
                
                #
                fontweight = "normal"
                if len(highlight) >0:
                    if row["i"] in highlight_i:
                        fontweight = "bold"
                #

                
                xy=(row["i"],row["scaled_P"]+0.2)
                xytext=(last_pos,1.15*maxy*arm_scale)
                
                if anno_fixed_arm_length is not None:
                    armB_length_in_point = anno_fixed_arm_length
                    xytext=(row["i"],row["scaled_P"]+0.2+anno_fixed_arm_length)
                
                if anno_count not in anno_d.keys():
                    #arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                    #                         connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                    arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                             connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                else:
                    # adjuest direction
                    xy=(row["i"],row["scaled_P"])
                    if anno_d[anno_count] in ["right","left","l","r"]:
                        if anno_d[anno_count]=="right" or anno_d[anno_count]=="r": 
                            armoffsetall = (ax1.transData.transform(xytext)[0]-ax1.transData.transform(xy)[0])*np.sqrt(2)
                            armoffsetb = arm_offset 
                            armoffseta = armoffsetall - armoffsetb   
                            arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                 connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                        elif anno_d[anno_count]=="left" or anno_d[anno_count]=="l":
                            arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                 connectionstyle="arc,angleA=-135,armA="+str(arm_offset)+",angleB=135,armB="+str(arm_offset)+",rad=0")
                    else:
                        if anno_d[anno_count][0]=="right" or anno_d[anno_count][0]=="r": 
                            armoffsetall = (ax1.transData.transform(xytext)[0]-ax1.transData.transform(xy)[0])*np.sqrt(2)
                            armoffsetb = anno_d[anno_count][1] 
                            armoffseta = armoffsetall - armoffsetb   
                            arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                 connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                        elif anno_d[anno_count]=="left" or anno_d[anno_count]=="l":
                            arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                 connectionstyle="arc,angleA=-135,armA="+str( anno_d[anno_count][1])+",angleB=135,armB="+str( anno_d[anno_count][1])+",rad=0")
                
                
                if "r" in mode:
                    arrowargs["color"] = "black" 
                    bbox_para=dict(boxstyle="round", fc="white",zorder=3)
                    for key,value in region_anno_bbox_args.items():
                        bbox_para[key]=value
                else:
                    bbox_para=None
                
                anno_default = {"rotation":40,"style":"italic","ha":"left","va":"bottom","fontsize":anno_fontsize,"fontweight":fontweight}
                for key,value in anno_args.items():
                    anno_default[key]=value

                if anno_adjust==True:
                    arrowargs=dict(arrowstyle='-|>', color='grey', shrinkA=10, linewidth=0.1, relpos=(0,0.5))
                anno_to_adjust = ax1.annotate(annotation_text,
                         xy=xy,
                         xytext=xytext,
                         bbox=bbox_para,
                         arrowprops=arrowargs,
                         zorder=100,
                         **anno_default
                         )
                anno_to_adjust_list.append(anno_to_adjust) 
                anno_count +=1
            #anno_adjust_keyargs = {"arrowprops":dict(arrowstyle='->', color='grey', linewidth=0.1,relpos=(0.5,0.5))}
            if anno_adjust==True:
                if verbose: log.write(" -Auto-adjusting text positions...")
                adjust_text(texts = anno_to_adjust_list,
                            autoalign=False, 
                            only_move={'points':'x', 'text':'x', 'objects':'x'},
                            ax=ax1,
                            precision=0.02,
                            force_text=(repel_force,repel_force),
                            expand_text=(1,1),
                            expand_objects=(0,0),
                            expand_points=(0,0),
                            va="bottom",
                            ha='left',
                            avoid_points=False,
                            lim =100
                            #kwargs = anno_adjust_keyargs
                            )
            
        else:
            if verbose: log.write(" -Skip annotating")
        
        # gene track (ax3) text adjustment
        if (gtf_path is not None ) and ("r" in mode):
            if len(texts_to_adjust_middle)>0:
                adjust_text(texts_to_adjust_middle,
                        autoalign=False, 
                        only_move={'points':'x', 'text':'x', 'objects':'x'},
                        ax=ax3,
                        precision=0,
                        force_text=(0.1,0),
                        expand_text=(1, 1),
                        expand_objects=(1,1),
                        expand_points=(1,1),
                        va="center",
                        ha='center',
                        avoid_points=False,
                        lim =1000)

            
# Creating Manhatann plot Finished #####################################################################

# QQ plot #########################################################################################################
    # ax2 qqplot
    if "qq" in mode:
        if verbose:log.write("Start to create QQ plot with "+str(len(sumstats))+" variants:")
        # plotting qq plots using processed data after cut and skip
        # select -log10 scaled p to plot
        p_toplot = sumstats["scaled_P"]
        # min p value for uniform distribution 
        minit=1/len(p_toplot)
        if stratified is False:
            # sort x,y for qq plot
            # high to low
            observed = p_toplot.sort_values(ascending=False)
            # uniform distribution using raw number -> -log10 -> observed number (omit variants with low -log10p)
            expected = -np.log10(np.linspace(minit,1,len(p_toplot_raw)))[:len(observed)]
            #p_toplot = sumstats["scaled_P"]
            ax2.scatter(expected,observed,s=8,color=colors[0])
        else:
            # stratified qq plot
            for i,(lower, upper) in enumerate(maf_bins):
                # extract data for a maf_bin
                databin = sumstats.loc[(sumstats["MAF"]>lower) &( sumstats["MAF"]<=upper),["MAF","scaled_P"]]
                # raw data : varaints with maf(eaf_raw) in maf_bin
                databin_raw = eaf_raw[(eaf_raw>lower) & (eaf_raw<=upper)]
                # sort x,y for qq plot
                # high to low
                observed = databin["scaled_P"].sort_values(ascending=False)
                # uniform distribution using raw number -> -log10 -> observed number (omit variants with low -log10p)
                expected = -np.log10(np.linspace(minit,1,max(len(databin_raw),len(databin))))[:len(observed)]
                label ="("+str(lower)+","+str(upper) +"]"
                ax2.scatter(expected,observed,s=8,color=maf_bin_colors[i],label=label)
                ax2_legend= ax2.legend(loc="best",fontsize=fontsize,markerscale=3,frameon=False)
                plt.setp(ax2_legend.texts, family="sans-serif")
        
        ax2.plot([skip,-np.log10(minit)],[skip,-np.log10(minit)],linestyle="--",color=sig_line_color)
        ax2.set_xlabel("Expected $-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
        ax2.set_ylabel("Observed $-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(True)
        
        # calculate genomic inflation factor and add annotation
        if gc:
            #gc was calculated using raw data (before cut and skip)
            p_toplot_raw = p_toplot_raw.rename(columns={"scaled_P":"MLOG10P"})
            if verbose and not include_chrXYMT : log.write(" -Excluding chrX,Y, MT from calculation of lambda GC.")
            lambdagc = lambdaGC(p_toplot_raw, mode="MLOG10P", include_chrXYMT=include_chrXYMT,log=log,verbose=False)
            #observedMedianChi2 = sp.stats.chi2.isf( np.median(np.power(10,-p_toplot_raw)) ,1)
            #expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
            #lambdagc=observedMedianChi2/expectedMedianChi2
            if verbose: log.write(" -Calculating lambda GC:",lambdagc)
            ax2.text(0.10, 1.03,"$\\lambda_{GC}$ = "+"{:.4f}".format(lambdagc),
                     horizontalalignment='left',
                     verticalalignment='top',
                     transform=ax2.transAxes,
                     fontsize=fontsize,family="sans-serif")
        
        #
        if cut == 0: ax2.set_ylim(skip,maxy*1.2)
        if cut:
            qcutline=ax2.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax2.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax2.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker],fontsize=fontsize,family="sans-serif")
            else:
                ax2.set_yticks([x for x in range(skip,cut+1,2)])
                ax2.set_yticklabels([x for x in range(skip,cut+1,2)],fontsize=fontsize,family="sans-serif")
        
        ax2.tick_params(axis='both', which='both', labelsize=fontsize)
        #
        if qtitle:
            ax2.set_title(qtitle,fontsize=title_fontsize,pad=10,family="sans-serif")

        if verbose: log.write("Finished creating QQ plot successfully!")
# Creating QQ plot Finished #############################################################################################
        

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


# +
#Helpers

def process_vcf(sumstats, vcf_path, region, log, verbose, pos , region_ld_threshold, vcf_chr_dict,tabix):
    if verbose: log.write("Start to load reference genotype...")
    if verbose: log.write(" -reference vcf path : "+ vcf_path)
    # load genotype data of the targeted region
    ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
    
    if verbose: log.write(" -Retrieving index...")
    if verbose: log.write(" -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])))
    # match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
    sumstats["REFINDEX"] = sumstats[pos].apply(lambda x: np.where(ref_genotype["variants/POS"] == x )[0][0] if np.any(ref_genotype["variants/POS"] == x) else None)
    # get lead variant id and pos
    lead_id = sumstats["scaled_P"].idxmax()
    lead_pos = sumstats.loc[lead_id,pos]
    
    # if lead pos is available: 
    if lead_pos in ref_genotype["variants/POS"]:
        
        # get ref index for lead snp
        lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

        # non-na other snp index
        other_snps_ref_index = sumstats["REFINDEX"].dropna().astype("int").values
        # get genotype 
        lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
        other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
        
        if verbose: log.write(" -Calculating Rsq...")
        
        if len(other_snp_genotype)>1:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
        else:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)
        sumstats.loc[~sumstats["REFINDEX"].isna(),"RSQ"] = valid_r2
    else:
        if verbose: log.write(" -Lead SNP not found in reference...")
        sumstats["RSQ"]=None
        
    sumstats["RSQ"] = sumstats["RSQ"].astype("float")
    sumstats["LD"] = 0
    for index,ld_threshold in enumerate(region_ld_threshold):
        if index==0:
            to_change_color = sumstats["RSQ"]>-1
            sumstats.loc[to_change_color,"LD"] = 1
        else:
            to_change_color = sumstats["RSQ"]>ld_threshold
            sumstats.loc[to_change_color,"LD"] = index+1
    
    sumstats.loc[lead_id,"LD"] = len(region_ld_threshold)+2
    sumstats["LEAD"]="Other variants"
    sumstats.loc[lead_id,"LEAD"] = "Lead variants"
    if verbose: log.write("Finished loading reference genotype successfully!")
    return sumstats


# -

def process_gtf(gtf_path,region,region_flank_factor,build,region_protein_coding,gtf_chr_dict,gtf_gene_name):
    #loading
    
    # chr to string datatype using gtf_chr_dict
    to_query_chrom = gtf_chr_dict[region[0]]

    # loading gtf data
    if gtf_path =="default" or gtf_path =="ensembl":
        
        gtf = get_gtf(chrom=to_query_chrom, build=build, source="ensembl")
    
    elif gtf_path =="refseq":

        gtf = get_gtf(chrom=to_query_chrom, build=build, source="refseq")
    
    else:
        # if user-provided gtf
        gtf = pd.read_csv(gtf_path,sep="\t",header=None, comment="#", low_memory=False,dtype={0:"string"})
        gtf = gtf.loc[gtf[0]==gtf_chr_dict[region[0]],:]

    # filter in region
    genes_1mb = gtf.loc[(gtf[0]==to_query_chrom)&(gtf[3]<region[2])&(gtf[4]>region[1]),:].copy()
    
    # extract biotype
    genes_1mb.loc[:,"gene_biotype"] = genes_1mb[8].str.extract(r'gene_biotype "([\w\.\_-]+)"')
    
    # extract gene name
    if gtf_gene_name is None:
        if gtf_path=="refseq":
            genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"').astype("string")
        elif gtf_path =="default" or gtf_path =="ensembl":
            genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_name "([\w\.-]+)"').astype("string")
        else:
            genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"').astype("string")
    else:
        pattern = r'{} "([\w\.-]+)"'.format(gtf_gene_name)
        genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(pattern).astype("string")

    # extract gene
    genes_1mb.loc[:,"gene"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"')
    
    # extract exon
    exons = genes_1mb.loc[genes_1mb[2]=="exon",:].copy()

    # extract protein coding gene
    if region_protein_coding is True:
        genes_1mb  =  genes_1mb.loc[genes_1mb["gene_biotype"]=="protein_coding",:].copy()
    
    #uniq genes
    ## get all record with 2nd column == gene
    uniq_gene_region = genes_1mb.loc[genes_1mb[2]=="gene",:].copy()
    
    ## extract region + flank
    flank = region_flank_factor * (region[2] - region[1])
    
    ## get left and right boundary
    uniq_gene_region["left"] = uniq_gene_region[3]-flank
    uniq_gene_region["right"] = uniq_gene_region[4]+flank

    # arrange gene track
    stack_dic = assign_stack(uniq_gene_region.sort_values([3]).loc[:,["name","left","right"]])  

    # map gene to stack and add stack column : minus stack
    uniq_gene_region["stack"] = -uniq_gene_region["name"].map(stack_dic)
    exons.loc[:,"stack"] = -exons.loc[:,"name"].map(stack_dic)

    # return uniq_gene_region (gene records with left and right boundary)
    # return exon records with stack number
    return uniq_gene_region, exons

def assign_stack(uniq_gene_region):

    stacks=[] ## stack : gene track
    stack_dic={} # mapping gene name to stack
    
    for index,row in uniq_gene_region.iterrows():
        if len(stacks)==0:
            # add first entry 
            stacks.append([(row["left"],row["right"])])
            stack_dic[row["name"]] = 0
        else:
            for i in range(len(stacks)):
                for j in range(len(stacks[i])):
                    # if overlap
                    if (row["left"]>stacks[i][j][0] and row["left"]<stacks[i][j][1]) or (row["right"]>stacks[i][j][0] and row["right"]<stacks[i][j][1]):
                        # if not last stack : break
                        if i<len(stacks)-1:
                            break
                        # if last stack : add a new stack
                        else:
                            stacks.append([(row["left"],row["right"])])
                            stack_dic[row["name"]] = i+1
                            break
                    # if no overlap       
                    else:
                        # not last in a stack
                        if j<len(stacks[i])-1:
                            #if in the middle
                            if row["left"]>stacks[i][j][1] and row["right"]<stacks[i][j+1][0]:
                                stacks[i].insert(j+1,(row["left"],row["right"]))
                                stack_dic[row["name"]] = i
                                break
                        #last one in a stack
                        elif row["left"]>stacks[i][j][1]:
                            stacks[i].append((row["left"],row["right"]))
                            stack_dic[row["name"]] = i
                            break
                if row["name"] in stack_dic.keys():
                    break         
    return stack_dic

def closest_gene(x,data,chrom="CHR",pos="POS",maxiter=20000,step=50):
        gene = data.gene_names_at_locus(contig=x[chrom], position=x[pos])
        if len(gene)==0:
            i=0
            while i<maxiter:
                distance = i*step
                gene_u = data.gene_names_at_locus(contig=x[chrom], position=x[pos]-distance)
                gene_d = data.gene_names_at_locus(contig=x[chrom], position=x[pos]+distance)
                if len(gene_u)>0 and len(gene_d)>0:
                    for j in range(0,step,1):
                        distance = (i-1)*step
                        gene_u = data.gene_names_at_locus(contig=x[chrom], position=x[pos]-distance-j)
                        gene_d = data.gene_names_at_locus(contig=x[chrom], position=x[pos]+distance+j)
                        if len(gene_u)>0:
                            return -distance,",".join(gene_u)
                        else:
                            return distance,",".join(gene_d)
                elif len(gene_u)>0:
                    return -distance,",".join(gene_u)
                elif len(gene_d)>0:
                    return +distance,",".join(gene_d)
                else:
                    i+=1
            return distance,"intergenic"
        else:
            return 0,",".join(gene)

