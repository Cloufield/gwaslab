import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import gwaslab as gl

# 20220310 ######################################################################################################

def mqqplot(insumstats,            
          chrom=None,
          pos=None,
          p=None,
          snpid=None,
          eaf=None,
          mlog10p="MLOG10P",
          scaled=False,
          mode="mqq",
          mqqratio=3,
          windowsizekb=500,
          anno=None,
          anno_set=[],
          anno_d={},
          arm_offset=50,
          arm_scale=1,
          cut=0,
          skip=0,
          cutfactor=10,
          cut_line_color="#ebebeb",  
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
          title =None,
          mtitle=None,
          qtitle=None,
          figargs= {"figsize":(15,5),"dpi":100},
          fontsize = 10,
          colors=["#597FBD","#74BAD3"],
          marker_size=(5,25),
          use_rank=False,
          verbose=True,
          repel_force=0.03,
          title_pad=1.08, 
          save=None,
          saveargs={"dpi":400,"facecolor":"white"},
          log=gl.Log()
          ):
    
# log.writeing meta info #######################################################################################
    
    if verbose: log.write("Start to plot manhattan/qq plot with the following basic settings:")
    if verbose: log.write(" -Genome-wide significance level is set to "+str(sig_level)+" ...")
    if verbose: log.write(" -Raw input contains "+str(len(insumstats))+" variants...")
    if verbose: log.write(" -Plot layout mode is : "+mode)
    if len(anno_set)>0 and ("m" in mode):
        if verbose: log.write(" -Variants to pinpoint : "+",".join(anno_set))    
    if len(highlight)>0 and ("m" in mode):
        if verbose: log.write(" -Loci to highlight : "+",".join(highlight))    
        if verbose: log.write(" -Highlight_window is set to: ", highlight_windowkb, " kb")  
    if len(pinpoint)>0 :
        if verbose: log.write(" -Variants to pinpoint : "+",".join(highlight))    
     
    

# Plotting mode selection ######################################################################################
    # ax1 : manhattanplot : 
    # ax2 : qq plot 
    if  mode=="qqm":   
        fig, (ax2, ax1) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [1, mqqratio]},**figargs)
    elif mode=="mqq":
        fig, (ax1, ax2) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [mqqratio, 1]},**figargs)
    elif mode=="m":
        fig, ax1 = plt.subplots(1, 1,**figargs)
    elif mode=="qq":
        fig, ax2 = plt.subplots(1, 1,**figargs) 
    else:
        raise ValueError("Please select one from the 4 modes: mqq/qqm/m/qq/")
        
# Read sumstats #################################################################################################
    usecols=[]
    if p in insumstats.columns:
        # p value is necessary for all modes.
            usecols.append(p)
    elif scaled is True:
            usecols.append(mlog10p)
    else:
        raise ValueError("Please make sure "+p+" column is in input sumstats.")
    
    if (chrom is not None) and (pos is not None) and ("m" in mode):
        # when manhattan plot, chrom and pos is needed.
        if chrom in insumstats.columns:
            usecols.append(chrom)
        else:
            raise ValueError("Please make sure "+chrom+" column is in input sumstats.")
        if pos in insumstats.columns:
            usecols.append(pos)
        else:
            raise ValueError("Please make sure "+pos+" column is in input sumstats.")
    
    if len(highlight)>0 or len(pinpoint)>0 or (snpid is not None):
        # read snpid when highlight/pinpoint is needed.
        if snpid in insumstats.columns:
            usecols.append(snpid)
        else:
            raise ValueError("Please make sure "+snpid+" column is in input sumstats.")
    
    if (stratified is True) and (eaf is not None):
        # read eaf when stratified qq plot is needed.
        if eaf in insumstats.columns:
            usecols.append(eaf)
        else:
            raise ValueError("Please make sure "+eaf+" column is in input sumstats.")
         
    if (anno is not None) and (anno != True):
        if (anno in insumstats.columns):
            if (anno not in usecols):
                usecols.append(anno)
        else:
            raise ValueError("Please make sure "+anno+" column is in input sumstats.")
    
    sumstats = insumstats.loc[:,usecols].copy()

#####################################################################################################################   
    #Standardize
    ## Annotation
    if len(anno_set)>0 and (anno is False):
        anno=True
    if (anno is not None) and (anno != True):
        sumstats["Annotation"]=sumstats.loc[:,anno].astype("string")   
    
    ## P value
    if scaled is True:
        sumstats["raw_P"] = pd.to_numeric(sumstats[mlog10p], errors='coerce')
    else:
        sumstats["raw_P"] = sumstats[p].astype("float64")
    
    ## CHR & POS
    if "m" in mode: 
        # CHR X,Y,MT conversion ############################
        if sumstats[chrom].dtype =="string":
            sumstats[chrom] = sumstats[chrom].map(gl.get_chr_to_number(out_chr=True),na_action="ignore")
        ## CHR
        sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
        ## POS
        sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')

      
    ## EAF
    if stratified is True: 
        sumstats["MAF"] = pd.to_numeric(insumstats[eaf], errors='coerce')
        sumstats.loc[sumstats["MAF"]>0.5,"MAF"] = 1 - sumstats.loc[sumstats["MAF"]>0.5,"MAF"]
        eaf_raw = sumstats["MAF"].copy()
    if len(highlight)>0 and ("m" in mode):
        sumstats["HUE"] = sumstats[chrom].astype("string")
            
#sanity check############################################################################################################
    if verbose: log.write("Start conversion and QC:")
    if "m" in mode: 
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
    if len(highlight)>0 and ("m" in mode):
        to_highlight = sumstats.loc[sumstats[snpid].isin(highlight),:]
        #assign colors: 0 is hightlight color
        for i,row in to_highlight.iterrows():
            target_chr = int(row[chrom])
            target_pos = int(row[pos])
            right_chr=sumstats[chrom]==target_chr
            up_pos=sumstats[pos]>target_pos-highlight_windowkb*1000
            low_pos=sumstats[pos]<target_pos+highlight_windowkb*1000
            sumstats.loc[right_chr&up_pos&low_pos,"HUE"]="0"

# P value conversion #####################################################################################################  
    pre_number=len(sumstats)
    sumstats = sumstats.dropna(subset=["raw_P"])
    after_number=len(sumstats)
    if verbose:log.write(" -Removed "+ str(pre_number-after_number) +" variants with nan in P column ...")
    
    if scaled:
        if verbose:log.write(" -P values are already converted to -log10(P)!")
        sumstats["scaled_P"] = sumstats["raw_P"].copy()
        sumstats["raw_P"] = np.power(10,-sumstats["scaled_P"])
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
    
        
    p_toplot_raw = sumstats["scaled_P"].copy()
    
    sumstats = sumstats.loc[sumstats["scaled_P"]>skip,:]
    
    # shrink at a certain value #########################################################################################
    maxy = sumstats["scaled_P"].max()
    if verbose: log.write(" -Maximum -log10(P) values is "+str(maxy) +" .")
    if cut:
        if verbose: log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")

        maxticker=int(np.round(sumstats["scaled_P"].max(skipna=True)))
        
        sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
        maxy = (maxticker-cut)/cutfactor + cut

# Manhattan plot ##########################################################################################################
## Create index for plotting ###################################################

    if verbose:log.write("Plotting "+str(len(sumstats))+" variants:")
    #sort & add id
    if "m" in mode: 
        sumstats = sumstats.sort_values([chrom,pos])
        if use_rank is True: sumstats["POS_RANK"] = sumstats.groupby(chrom)[pos].rank("dense", ascending=True)
        sumstats["id"]=range(len(sumstats))
        sumstats=sumstats.set_index("id")

        #create a position dictionary
        if use_rank is True: 
            posdic = sumstats.groupby(chrom)["POS_RANK"].max()
        else:
            posdic = sumstats.groupby(chrom)[pos].max()
        posdiccul = dict(posdic)
        for i in range(0,26):
            if i in posdiccul: continue
            else: posdiccul[i]=0

        for i in range(2,sumstats[chrom].max()+1):
            posdiccul[i]=posdiccul[i-1]+posdiccul[i]

        #convert base pair postion to x axis position
        sumstats["add"]=sumstats[chrom].apply(lambda x : posdiccul[int(x)-1])
        if use_rank is True: 
            sumstats["i"]=sumstats["POS_RANK"]+sumstats["add"]
        else:
            sumstats["i"]=sumstats[pos]+sumstats["add"]
        

        #for plot
        chrom_df=sumstats.groupby(chrom)['i'].median()
        sumstats["i"]=sumstats["i"]+((sumstats[chrom].map(dict(chrom_df)).astype("int")))*0.02
        sumstats["i"] = sumstats["i"].astype("Int64")
## Assign marker size ##############################################
        
        sumstats["s"]=1
        sumstats.loc[sumstats["scaled_P"]>-np.log10(5e-4),"s"]=2
        sumstats.loc[sumstats["scaled_P"]>-np.log10(suggestive_sig_level),"s"]=3
        sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level),"s"]=4
        sumstats["chr_hue"]=sumstats[chrom].astype("string")
        
## Manhatann plot ###################################################
        if len(highlight) >0:
            plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P',
                               hue='chr_hue',
                               palette=sns.color_palette(colors,n_colors=sumstats[chrom].nunique()),
                               legend=None,
                               size="s",
                               sizes=marker_size,
                               linewidth=0,
                               zorder=2,ax=ax1)   
            if verbose: log.write(" -Highlighting target loci...")
            sns.scatterplot(data=sumstats.loc[sumstats["HUE"]=="0"], x='i', y='scaled_P',
                   hue="HUE",
                   palette={"0":highlight_color},
                   legend=None,
                   size="s",
                   sizes=(marker_size[0]+1,marker_size[1]+1),
                   linewidth=0,
                   zorder=3,ax=ax1)  
            highlight_i = sumstats.loc[sumstats[snpid].isin(highlight),"i"].values
            
        else:
            plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P',
                   hue='chr_hue',
                   palette= sns.color_palette(colors,n_colors=sumstats[chrom].nunique()),
                   legend=None,
                   size="s",
                   sizes=marker_size,
                   linewidth=0,
                   zorder=2,ax=ax1)   
        
        if (len(pinpoint)>0):
            if sum(sumstats[snpid].isin(pinpoint))>0:
                to_pinpoint = sumstats.loc[sumstats[snpid].isin(pinpoint),:]
                if verbose: log.write(" -Pinpointing target vairants...")
                ax1.scatter(to_pinpoint["i"],to_pinpoint["scaled_P"],color=pinpoint_color,zorder=3,s=marker_size[1]+1)
            else:
                if verbose: log.write(" -Target vairants to pinpoint were not found. Skip pinpointing process...")
 
                
        chrom_df=sumstats.groupby(chrom)['i'].median()  
        plot.set_xlabel(chrom); 
        plot.set_xticks(chrom_df);

        plot.set_xticklabels(chrom_df.index.map(gl.get_number_to_chr()))

        sigline = plot.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
        if cut == 0: plot.set_ylim(skip,maxy*1.2)
        if cut:
            cutline=plot.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                plot.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                plot.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker])
            else:
                plot.set_yticks([x for x in range(skip,cut+1,2)])
                plot.set_yticklabels([x for x in range(skip,cut+1,2)])

# Get top variants for annotation #######################################################
        if (anno and anno!=True) or (len(anno_set)>0):
            if len(anno_set)>0:
                to_annotate=sumstats.loc[sumstats[snpid].isin(anno_set),:]
                if to_annotate.empty is not True:
                    if verbose: log.write(" -Found "+str(len(to_annotate))+" specified variants to annotate...")
            else:
                to_annotate=gl.getsig(sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level)],
                               snpid,
                               chrom,
                               pos,
                               "raw_P",
                               sig_level=sig_level,
                               windowsizekb=windowsizekb,
                               verbose=False)
                if to_annotate.empty is not True:
                    if verbose: log.write(" -Found "+str(len(to_annotate))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")
        else:
            to_annotate=gl.getsig(sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level)],
                               "i",
                               chrom,
                               pos,
                               "raw_P",
                               windowsizekb=windowsizekb,
                               verbose=False,
                               sig_level=sig_level)
            if to_annotate.empty is not True:
                if verbose: log.write(" -Found "+str(len(to_annotate))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")
        
        
# Add Annotation to manhattan plot #######################################################
        if anno and (to_annotate.empty is not True):
            #initiate a list for text and a starting position
            text = []
            last_pos=0
            anno_count=0
            if anno==True:
                    annotation_col="CHR:POS"
            elif anno:
                    annotation_col=anno
            if verbose: log.write(" -Annotating using column "+annotation_col+"...")
                
            for rowi,row in to_annotate.iterrows():
                # avoid text overlapping
                if row["i"]>last_pos+repel_force*sumstats["i"].max():
                    last_pos=row["i"]
                else:
                    last_pos+=repel_force*sumstats["i"].max()
                # data to pixels
                armB_length_in_point = plot.transData.transform((skip,1.15*maxy))[1]-plot.transData.transform((skip, row["scaled_P"]+1))[1]-arm_offset/2
                #  
                armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else plot.transData.transform((skip, maxy+2))[1]-plot.transData.transform((skip,  row["scaled_P"]+1))[1] 
                
                if anno==True:
                    annotation_text="Chr"+ str(row[chrom]) +":"+ str(int(row[pos]))
                    annotation_col="CHR:POS"
                elif anno:
                    annotation_text=row["Annotation"]
                    annotation_col=anno
                    
                fontweight = "normal"
                if len(highlight) >0:
                    if row["i"] in highlight_i:
                        fontweight = "bold"
                
                xy=(row["i"],row["scaled_P"]+0.2)
                xytext=(last_pos,1.15*maxy)
                if anno_count not in anno_d.keys():
                    arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                             connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                else:
                    if anno_d[anno_count]=="right": 
                        armoffsetall = (plot.transData.transform(xytext)[0]-plot.transData.transform(xy)[0])*np.sqrt(2)
                        armoffsetb = arm_offset 
                        armoffseta = armoffsetall - armoffsetb
                        
                        arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                             connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                    else:
                        arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                             connectionstyle="arc,angleA=-135,armA="+str(arm_offset)+",angleB=135,armB="+str(arm_offset)+",rad=0")
                
                text.append(plot.annotate(annotation_text,
                                             style='italic',
                                             xy=xy,
                                             xytext=xytext,rotation=40,
                                             ha="left",va="bottom",
                                             fontsize=fontsize,
                                             fontweight=fontweight,
                                             arrowprops=arrowargs
                                            )
                           )
                anno_count +=1

            
        else:
            if verbose: log.write(" -Skip annotating")

        plot.set_ylabel("$-log_{10}(P)$",fontsize=fontsize)
        plot.set_xlabel("Chromosomes",fontsize=fontsize)
        plot.spines["top"].set_visible(False)
        plot.spines["right"].set_visible(False)
        plot.spines["left"].set_visible(True)

        if verbose: log.write(" -Created Manhattan plot successfully!")
        if mtitle and anno and len(to_annotate)>0: 
            pad=(plot.transData.transform((skip, title_pad*maxy))[1]-plot.transData.transform((skip, maxy)))[1]
            plot.set_title(mtitle,pad=pad)
        elif mtitle:
            plot.set_title(mtitle,fontsize=fontsize)
# Creating Manhatann plot Finished #####################################################################

# QQ plot #########################################################################################################
    # ax2 qqplot
    if "qq" in mode:
        
        p_toplot = sumstats["scaled_P"]
            # select -log10 scaled p to plot
            # sort x,y for qq plot
        minit=1/len(p_toplot)
        if stratified is False:
            observed = p_toplot.sort_values(ascending=False)
            expected = -np.log10(np.linspace(minit,1,len(p_toplot_raw)))[:len(observed)]
            #p_toplot = sumstats["scaled_P"]
            ax2.scatter(expected,observed,s=8,color=colors[0])
        else:
            # stratified qq plot
            for i,(lower, upper) in enumerate(maf_bins):
                databin = sumstats.loc[(sumstats["MAF"]>lower) &( sumstats["MAF"]<=upper),["MAF","scaled_P"]]
                databin_raw = eaf_raw[(eaf_raw>lower) & (eaf_raw<=upper)]
                observed = databin["scaled_P"].sort_values(ascending=False)
                expected = -np.log10(np.linspace(minit,1,max(len(databin_raw),len(databin))))[:len(observed)]
                label ="("+str(lower)+","+str(upper) +"]"
                ax2.scatter(expected,observed,s=8,color=maf_bin_colors[i],label=label)
                ax2.legend(loc="best",fontsize=10,markerscale=3,frameon=False)
        
        ax2.plot([skip,-np.log10(minit)],[skip,-np.log10(minit)],linestyle="--",color=sig_line_color)

        ax2.set_xlabel("Expected $-log_{10}(P)$",fontsize=fontsize)
        ax2.set_ylabel("Observed $-log_{10}(P)$",fontsize=fontsize)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(True)
        
        # calculate genomic inflation factor and add annotation
        if gc:
            
            observedMedianChi2 = sp.stats.chi2.isf( np.median(np.power(10,-p_toplot_raw)) ,1)
            expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
            lambdagc=observedMedianChi2/expectedMedianChi2
            if verbose: log.write(" -Calculating GC using P :",lambdagc)
            ax2.text(0.10, 1.03,"$\\lambda_{GC}$ = "+"{:.4f}".format(lambdagc),
                     horizontalalignment='left',
                     verticalalignment='top',
                     transform=ax2.transAxes,
                     fontsize=fontsize)
        
        #
        if cut == 0: ax2.set_ylim(skip,maxy*1.2)
        if cut:
            qcutline=ax2.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax2.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax2.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker])
            else:
                ax2.set_yticks([x for x in range(skip,cut+1,2)])
                ax2.set_yticklabels([x for x in range(skip,cut+1,2)])
        
        #
        if qtitle:
            ax2.set_title(qtitle,fontsize=fontsize,pad=10)

        if verbose: log.write(" -Created QQ plot successfully!")
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
        fig.suptitle(title ,x=0.5, y=1.05)
    else:
        fig.suptitle(title ,x=0.5,y=1)

# Return matplotlib figure object #######################################################################################
    return fig, log