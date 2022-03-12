import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import gwaslab as gl
#20220310
def mqqplot(insumstats,
          chrom,
          pos,
          p,
          snpid=None,
          scaled=False,
          cut=0,
          cutfactor=10,
          mode="mqq",
          mqqratio=3,
          cut_line_color="#ebebeb",
          windowsizekb=500,
          anno=None,
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_level=5e-6,
          stratified=False,
          eaf=None,
          maf_bins=[(0, 0.01), (0.01, 0.05), (0.05, 0.25),(0.25,0.5)],
          maf_bin_colors = ["#f0ad4e","#5cb85c", "#5bc0de","#000042"],
          highlight = [],
          highlight_color ="#33FFA0",
          highlight_windowkb = 500,
          title =None,
          mtitle=None,
          qtitle=None,
          figsize =(15,5),
          fontsize = 10,
          colors = ["#000042", "#7878BA"],
          verbose=True,
          repel_force=0.03,
          gc=True,
          save=None,
          saveargs={"dpi":400,"facecolor":"white"}  
          ):
    
# Printing meta info ##################################################################
    if verbose: print("Basic settings:")
    if verbose: print("  - Genome-wide significance level is set to "+str(sig_level)+" ...")
    if verbose: print("  - Raw input contains "+str(len(insumstats))+" variants...")
    if verbose: print("  - Ploy layout mode is : "+mode)
    if verbose: print("Start conversion and QC:")

# Plotting mode selection ##############################################################
    # ax1 : manhattanplot : 
    # ax2 : qq plot 
    
    if   mode=="qqm":   
        fig, (ax2, ax1) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [1, mqqratio]},figsize=figsize)
    elif mode=="mqq":
        fig, (ax1, ax2) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [mqqratio, 1]},figsize=figsize)
    elif mode=="m":
        fig, ax1 = plt.subplots(1, 1,figsize=figsize)
    elif mode=="qq":
        fig, ax2 = plt.subplots(1, 1,figsize=figsize) 
    else:
        raise ValueError("Please select one from the 4 modes: mqq/qqm/m/qq/")
        
# Read sumstats ###############################
    
    sumstats = pd.DataFrame()
    sumstats["P-value_association"] = insumstats[p].astype(np.float64)
    if "m" in mode: 
        sumstats["CHROM"] = insumstats[chrom].copy()
        # CHR X,Y,MT conversion ############################
        sumstats.loc[sumstats["CHROM"]=="X","CHROM"] = "23"
        sumstats.loc[sumstats["CHROM"]=="Y","CHROM"] = "24"
        sumstats.loc[sumstats["CHROM"]=="MT","CHROM"] = "25"
        sumstats["CHROM"] = sumstats["CHROM"].astype("int")
        
    if "m" in mode: 
        sumstats["POS"] = insumstats[pos].astype("int")
    if anno and anno!=True:
        sumstats["Annotation"]=insumstats[anno].astype("string")
    if stratified is True:
        sumstats["MAF"] = insumstats[eaf].copy()
        sumstats.loc[sumstats["MAF"]>0.5,"MAF"] = 1 - sumstats.loc[sumstats["MAF"]>0.5,"MAF"]
    
    if len(highlight)>0:
        sumstats["MARKERNAME"] = insumstats[snpid]
        sumstats["HUE"] = sumstats["CHROM"]
        to_highlight = sumstats.loc[sumstats["MARKERNAME"].isin(highlight)]
        for i,row in to_highlight.iterrows():
            sumstats.loc[(sumstats["CHROM"]==row["CHROM"])
                         &(sumstats["POS"]>row["POS"]-highlight_windowkb*1000)
                         &(sumstats["POS"]<row["POS"]+highlight_windowkb*1000),"HUE"]="0"

    
            
# P value conversion ###############################  
    if scaled:
        if verbose:print("  - P values are already converted to -log10(P)!")
        sumstats["scaled_P"] = sumstats["P-value_association"]
    else:
        if verbose:print("  - P values are being converted to -log10(P)...")
        #sanity check
        outside_01_p_value_number=len(sumstats.loc[(sumstats["P-value_association"]>1)|(sumstats["P-value_association"]<=0),:])
        if verbose:print("  - Sanity check: "+ str(outside_01_p_value_number) +" variants with P value outside of (0,1] will be removed...")
        sumstats = sumstats.loc[(sumstats["P-value_association"]<=1)&(sumstats["P-value_association"]>0),:]
        sumstats["scaled_P"] = -np.log10(sumstats["P-value_association"])

    if verbose: print("  - Sanity check: "+str(len(sumstats[sumstats["scaled_P"].isin([np.nan, np.inf, -np.inf, float('inf'), None])])) + " na/inf/-inf variants will be removed..." )
    sumstats = sumstats[~sumstats["scaled_P"].isin([np.nan, np.inf, -np.inf, float('inf'), None])]

    # shrink at a certain value #############################################
    maxy = sumstats["scaled_P"].max()
    if verbose: print("  - Maximum -log10(P) values is "+str(maxy) +" .")
    if cut:
        if verbose: print("  - Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")
        maxticker=int(np.round(sumstats["scaled_P"].max()))
        sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
        maxy = (maxticker-cut)/cutfactor + cut

# Manhattan plot #####################################################################################################
## Create index for plotting ###################################################

    if verbose:print("Plotting "+str(len(sumstats))+" variants:")
    #sort & add id
    if "m" in mode: 
        sumstats = sumstats.sort_values(["CHROM","POS"])
        sumstats["id"]=range(len(sumstats))
        sumstats=sumstats.set_index("id")

        #create a position dictionary
        posdic = sumstats.groupby("CHROM")["POS"].max()

        posdiccul = dict(posdic)
        for i in range(0,26):
            if i in posdiccul: continue
            else: posdiccul[i]=0



        for i in range(2,sumstats["CHROM"].nunique()+1):
            posdiccul[i]=posdiccul[i-1]+posdiccul[i]

        #convert base pair postion to x axis position
        sumstats["add"]=sumstats["CHROM"].apply(lambda x : posdiccul[int(x)-1])
        sumstats["i"]=sumstats["POS"]+sumstats["add"]

        #for plot
        chrom_df=sumstats.groupby('CHROM')['i'].median()
        chrom_df=chrom_df+((chrom_df.index.astype(int))-1)*len(sumstats)*0.02

        sumstats["i"]=sumstats["i"]+((sumstats["CHROM"].astype(int))-1)*len(sumstats)*0.01

## Assign marker size ##############################################

        sumstats["s"]=1
        sumstats.loc[sumstats["scaled_P"]>-np.log10(5e-4),"s"]=2
        sumstats.loc[sumstats["scaled_P"]>-np.log10(suggestive_sig_level),"s"]=3
        sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level),"s"]=4

## Manhatann plot ###################################################
        if len(highlight) >0:
            plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P',
                               hue='CHROM',
                               palette=sns.color_palette(colors,n_colors=sumstats["CHROM"].nunique()),
                               legend=None,
                               size="s",
                               sizes=(10,35),
                               linewidth=0,
                               zorder=2,ax=ax1)   
            if verbose: print("  - Loci to highlight :",highlight)
            if verbose: print("  - Highlight range +- :",highlight_windowkb , " kb")
            sns.scatterplot(data=sumstats.loc[sumstats["HUE"]=="0"], x='i', y='scaled_P',
                   hue="HUE",
                   palette={"0":highlight_color},
                   legend=None,
                   size="s",
                   sizes=(45,45),
                   linewidth=0,
                   zorder=3,ax=ax1)   
        else:
            plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P',
                   hue='CHROM',
                   palette= sns.color_palette(colors,n_colors=sumstats["CHROM"].nunique()),
                   legend=None,
                   size="s",
                   sizes=(10,40),
                   linewidth=0,
                   zorder=2,ax=ax1)   

        chrom_df = (sumstats.groupby('CHROM')['i'].max() + sumstats.groupby('CHROM')['i'].min())/2
        plot.set_xlabel('CHROM'); 
        plot.set_xticks(chrom_df);
        chromosome_conversion_dict = {i:str(i) for i in range(1,23)}
        chromosome_conversion_dict[23] = "X"
        chromosome_conversion_dict[24] = "Y"
        chromosome_conversion_dict[25] = "MT"
        plot.set_xticklabels(chrom_df.index.map(chromosome_conversion_dict))

        sigline = plot.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
        if cut:
            cutline=plot.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                plot.set_yticks([x for x in range(0,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                plot.set_yticklabels([x for x in range(0,cut+1,2)]+[maxticker])
            else:
                plot.set_yticks([x for x in range(0,cut+1,2)])
                plot.set_yticklabels([x for x in range(0,cut+1,2)])

# Get top variants for annotation #######################################################
        if anno and anno!=True:
            to_annotate=gl.getsig(sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level)],
                               "Annotation",
                               'CHROM',
                               "POS",
                               "P-value_association",
                               sig_level=sig_level,
                               windowsizekb=windowsizekb,
                               verbose=False)
        else:
            to_annotate=gl.getsig(sumstats.loc[sumstats["scaled_P"]>-np.log10(sig_level)],
                               "i",
                               'CHROM',
                               "POS",
                               "P-value_association",
                               windowsizekb=windowsizekb,
                               verbose=False,
                               sig_level=sig_level)
        
        if to_annotate.empty is not True:
            if verbose: print("  - Found "+str(len(to_annotate))+" significant variants with a sliding window size of "+str(windowsizekb)+" kb...")
# Add Annotation to manhattan plot #######################################################
        if anno and (to_annotate.empty is not True):
            #initiate a list for text and a starting position
            text = []
            last_pos=0

            for rowi,row in to_annotate.iterrows():
                # avoid text overlapping
                if row["i"]>last_pos+repel_force*sumstats["i"].max():
                    last_pos=row["i"]
                else:
                    last_pos+=repel_force*sumstats["i"].max()
                # data to pixels
                armB_length_in_point=plot.transData.transform((0, 0.95*maxy))[1]-plot.transData.transform((0, row["scaled_P"]+1))[1]

                #  
                armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else plot.transData.transform((0, maxy+2))[1]-plot.transData.transform((0,  row["scaled_P"]+1))[1] 

                if anno==True:
                    annotation_text="Chr"+ row["CHROM"] +":"+ str(int(row["POS"]))
                    annotation_col="CHR:POS"
                elif anno:
                    annotation_text=row["Annotation"]
                    annotation_col=anno

                text.append(plot.annotate(annotation_text,
                                             xy=(row["i"],row["scaled_P"]+1),
                                             xytext=(last_pos,1.15*maxy),rotation=40,
                                             ha="left",va="bottom",
                                             fontsize=fontsize,
                                             arrowprops=dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                             connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                                            )
                           )

            if verbose: print("  - Annotating using column "+annotation_col+"...")
        else:
            if verbose: print("  - Skip annotating")

        plot.set_ylabel("$-log_{10}(P)$",fontsize=fontsize)
        plot.set_xlabel("Chromosomes",fontsize=fontsize)
        plot.spines["top"].set_visible(False)
        plot.spines["right"].set_visible(False)
        plot.spines["left"].set_visible(True)

        if verbose: print("  - Created Manhattan plot successfully!")
        if mtitle and anno and len(to_annotate)>0: 
            pad=(plot.transData.transform((0, 1.35*maxy))[1]-plot.transData.transform((0, maxy)))[1]
            plot.set_title(mtitle,pad=pad)
        elif title:
            plot.set_title(mtitle,fontsize=fontsize)
# Creating Manhatann plot Finished #####################################################################

# QQ plot #########################################################################################################
    # ax2 qqplot
    if "qq" in mode:
        
            # select -log10 scaled p to plot
        p_toplot = sumstats["scaled_P"]
            # sort x,y for qq plot
        minit=1/len(p_toplot)
        if stratified is False:
            observed = p_toplot.sort_values(ascending=False)
            expected = -np.log10(np.linspace(minit,1,len(observed)))
            ax2.scatter(expected,observed,s=8,color=colors[0])
        else:
            # stratified qq plot
            for i,(lower, upper) in enumerate(maf_bins):
                databin = sumstats.loc[(sumstats["MAF"]>lower) &( sumstats["MAF"]<=upper),["MAF","scaled_P"]]
                observed = databin["scaled_P"].sort_values(ascending=False)
                expected = -np.log10(np.linspace(minit,1,len(observed)))
                label ="("+str(lower)+","+str(upper) +"]"
                ax2.scatter(expected,observed,s=8,color=maf_bin_colors[i],label=label)
                ax2.legend(loc="center left",fontsize=10,markerscale=3)
        
        ax2.plot([0,-np.log10(minit)],[0,-np.log10(minit)],linestyle="--",color=sig_line_color)

        ax2.set_xlabel("Expected $-log_{10}(P)$",fontsize=fontsize)
        ax2.set_ylabel("Observed $-log_{10}(P)$",fontsize=fontsize)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(True)
        
        # calculate genomic inflation factor and add annotation
        if gc:
            observedMedianChi2 = sp.stats.chi2.isf( np.median(np.power(10,-p_toplot)) ,1)
            expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
            lambdagc=observedMedianChi2/expectedMedianChi2
            ax2.text(0.05, 0.95,"$\\lambda_{GC}$ = "+"{:.4f}".format(lambdagc),
                     horizontalalignment='left',
                     verticalalignment='top',
                     transform=ax2.transAxes,
                     fontsize=fontsize)
        
        #
        if cut:
            qcutline=ax2.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax2.set_yticks([x for x in range(0,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax2.set_yticklabels([x for x in range(0,cut+1,2)]+[maxticker])
            else:
                ax2.set_yticks([x for x in range(0,cut+1,2)])
                ax2.set_yticklabels([x for x in range(0,cut+1,2)])
        
        #
        if qtitle:
            ax2.set_title(qtitle,fontsize=fontsize,pad=10)

        if verbose: print("  - Created QQ plot successfully!")
# Creating QQ plot Finished #############################################################################################


# Saving plot ##########################################################################################################
    if save:
        if verbose: print("Saving plot:")
        if save==True:
            fig.savefig("./"+mode+"_plot.png",bbox_inches="tight",**saveargs)
            print("  - Saved to "+ "./"+mode+"_plot.png" + " successfully!" )
        else:
            fig.savefig(save,bbox_inches="tight",**saveargs)
            print("  - Saved to "+ save + " successfully!" )
    
    # add title 
    if title and anno and len(to_annotate)>0:
        # increase height if annotation 
        fig.suptitle(title ,x=0.5, y=1.2)
    else:
        fig.suptitle(title ,x=0.5,y=1)

# Return matplotlib figure object #######################################################################################
    return fig