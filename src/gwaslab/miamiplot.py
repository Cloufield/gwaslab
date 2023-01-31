import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import gwaslab as gl
from gwaslab.Log import Log
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
import gc as garbage_collect

def plot_miami( 
          path1,
          path2,
          cols1=["CHR","POS","P"],
          cols2=["CHR","POS","P"],
          sep=["\t","\t"],
          region=None,
          region_step = 21,
          region_grid = False,
          region_grid_line = {"linewidth": 2,"linestyle":"--"},
          region_lead_grid = True,
          region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"},
          anno=None,
          anno_set=[],
          anno_set1=[],
          anno_set2=[],
          anno_alias1={},
          anno_alias2={},
          anno_d1={},
          anno_d2={},
          anno_args={},
          anno_fixed_arm_length=None,
          anno_source = "ensembl",
          arm_offset=50,
          arm_scale=1,
          arm_scale_d=None,
          highlight  = [],
          highlight1 = [],
          highlight2 = [],
          highlight_color="#CB132D",
          highlight_windowkb = 500,
          pinpoint=[],
          pinpoint1=[],
          pinpoint2=[],
          pinpoint_color ="red",
          titles=["",""],
          titles_pad=[0.2,0.2], 
          cut=0,
          skip=0,
          build="19",
          cutfactor=10,
          readcsv_args={},
          cut_line_color="#ebebeb",  
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_level=5e-6,
          region_hspace = 0.1,
          windowsizekb=500,
          dpi=100,
          figargs= {"figsize":(15,5),"dpi":100},
          fontsize = 10,
          colors1=["#597FBD","#74BAD3"],
          colors2=["#597FBD","#74BAD3"],
          scatter_kwargs={},
          marker_size=(5,25),
          verbose=True,
          large_number=10000000000,
          repel_force=0.03,
          save=None,
          saveargs={"dpi":100,"facecolor":"white"},
          log=Log()
          ):
    ## load sumstats1 ###########################################################################################################
    if verbose: log.write("Start to plot miami plot with the following basic settings:")
    if verbose: log.write(" -Genome-wide significance level is set to "+str(sig_level)+" ...")
    if len(anno_set)>0 :
        if verbose: log.write(" -Variants to annotate : ", anno_set)    
    if len(highlight)>0 or len(highlight1)>0 or len(highlight2)>0 :
        if verbose: log.write(" -Loci to highlight : ", highlight,highlight1,highlight2)    
        if verbose: log.write(" -Highlight_window is set to: ", highlight_windowkb, " kb")  
    if len(pinpoint)>0 or len(pinpoint1)>0 or len(pinpoint2)>0 :
        if verbose: log.write(" -Variants to pinpoint : ",pinpoint,pinpoint1,pinpoint2)      
    if region is not None:
        if verbose: log.write(" -Region to plot : chr"+str(region[0])+":"+str(region[1])+"-"+str(region[2])+".")  
    if dpi!=100:
        figargs["dpi"] = dpi
        
    for i in range(len(highlight)):
        highlight1.append(highlight[i])
        highlight2.append(highlight[i])
        highlight[i] = highlight[i][0] * large_number + highlight[i][1]
    for i in range(len(highlight1)):
        highlight1[i] = highlight1[i][0] * large_number + highlight1[i][1]
    for i in range(len(highlight2)):
        highlight2[i] = highlight2[i][0] * large_number + highlight2[i][1]
    
    for i in range(len(pinpoint)):
        pinpoint1.append(pinpoint[i])
        pinpoint2.append(pinpoint[i])
        pinpoint[i] = pinpoint[i][0] * large_number + pinpoint[i][1] 
    for i in range(len(pinpoint1)):
        pinpoint1[i] = pinpoint1[i][0] * large_number + pinpoint1[i][1]
    for i in range(len(pinpoint2)):
        pinpoint2[i] = pinpoint2[i][0] * large_number + pinpoint2[i][1] 
    
    for i in range(len(anno_set)):
        anno_set1.append(anno_set[i])
        anno_set2.append(anno_set[i])
        anno_set[i] = anno_set[i][0] * large_number + anno_set[i][1]
    for i in range(len(anno_set1)):
        anno_set1[i] = anno_set1[i][0] * large_number + anno_set1[i][1]
    for i in range(len(anno_set2)):
        anno_set2[i] = anno_set2[i][0] * large_number + anno_set2[i][1]

        
        
        
    if verbose: log.write(" -Loading sumstats1:" + path1)
    if verbose: log.write(" -Sumstats1 CHR,POS,P information will be obtained from:",cols1)
    sumstats1 = pd.read_csv(path1,sep=sep[0],usecols=cols1,dtype={cols1[0]:"string",cols1[1]:"Int64",cols1[2]:"float64"},**readcsv_args)
    sumstats1 = sumstats1.rename(columns={cols1[0]:"CHR",cols1[1]:"POS",cols1[2]:"P"})
    sumstats1["CHR"] = np.floor(pd.to_numeric(sumstats1["CHR"].map(get_chr_to_number(),na_action="ignore"), errors='coerce')).astype('Int64')
    sumstats1["POS"] = np.floor(pd.to_numeric(sumstats1["POS"], errors='coerce')).astype('Int64')   
    
    ## load sumstats2 ###########################################################################################################
    if verbose: log.write(" -Loading sumstats2:" + path2)
    if verbose: log.write(" -Sumstats2 CHR,POS,P information will be obtained from:",cols2)
    sumstats2 = pd.read_csv(path2,sep=sep[1],usecols=cols2,dtype={cols1[0]:"string",cols1[1]:"Int64",cols1[2]:"float64"},**readcsv_args)
    sumstats2 = sumstats2.rename(columns={cols2[0]:"CHR",cols2[1]:"POS",cols2[2]:"P"})
    sumstats2["CHR"] = np.floor(pd.to_numeric(sumstats2["CHR"].map(get_chr_to_number(),na_action="ignore"), errors='coerce')).astype('Int64')
    sumstats2["POS"] = np.floor(pd.to_numeric(sumstats2["POS"], errors='coerce')).astype('Int64')
    
    ## create merge index
    sumstats1["TCHR+POS"]=sumstats1["CHR"]*large_number + sumstats1["POS"]
    sumstats2["TCHR+POS"]=sumstats2["CHR"]*large_number + sumstats2["POS"]
    sumstats1["TCHR+POS"] = sumstats1["TCHR+POS"].astype('Int64')
    sumstats2["TCHR+POS"] = sumstats2["TCHR+POS"].astype('Int64')
    sumstats1.drop(labels=["CHR","POS"],axis=1,inplace=True)
    sumstats2.drop(labels=["CHR","POS"],axis=1,inplace=True)
    sumstats1.dropna(inplace=True)
    sumstats2.dropna(inplace=True)
    ## merging ###########################################################################################################
    if verbose: log.write(" - Merging sumstats using chr and pos...")
    merged_sumstats = pd.merge(sumstats1,sumstats2,on="TCHR+POS",how="outer",suffixes=('_1', '_2'))
    merged_sumstats["CHR"] = np.divmod(merged_sumstats["TCHR+POS"],large_number)[0]
    merged_sumstats["POS"] = np.divmod(merged_sumstats["TCHR+POS"],large_number*merged_sumstats["CHR"])[1]
    
    del(sumstats1)
    del(sumstats2)
    garbage_collect.collect()
    
    chrom = "CHR"
    pos="POS"
    
    merged_sumstats["scaled_P_1"] = -np.log10(merged_sumstats["P_1"])
    merged_sumstats["scaled_P_2"] = -np.log10(merged_sumstats["P_2"])
    
    if skip >0:
        sumstats = merged_sumstats.loc[(merged_sumstats["scaled_P_1"]>skip) | ( merged_sumstats["scaled_P_2"]>skip),:]  
    else:
        sumstats = merged_sumstats
    
    
    
    
    
    if region is not None:
        region_chr = region[0]
        region_start = region[1]
        region_end = region[2]
        marker_size=(25,45)
        if verbose:log.write(" -Extract SNPs in region : chr"+str(region_chr)+":"+str(region[1])+"-"+str(region[2])+ "...")
        in_region_snp = (sumstats[chrom]==region_chr) &(sumstats[pos]<region_end) &(sumstats[pos]>region_start)
        if verbose:log.write(" -Extract SNPs in specified regions: "+str(sum(in_region_snp)))
        sumstats = sumstats.loc[in_region_snp,:]

    
    
    
    
    
    
    
    sumstats = sumstats.sort_values([chrom,pos])
    sumstats["id"]=range(len(sumstats))
    sumstats=sumstats.set_index("id")

    #create a df , groupby by chromosomes , and get the maximum position
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
    sumstats["i"]=sumstats[pos]+sumstats["add"]

    #for plot, get the chr text tick position      
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
    
    
    
    ## highlight
    if len(highlight1)>0 :
        to_highlight1 = sumstats.loc[sumstats["TCHR+POS"].isin(highlight1),:]
        #assign colors: 0 is hightlight color
        if len(to_highlight1)>0:
            for i,row in to_highlight1.iterrows():
                target_chr = int(row[chrom])
                target_pos = int(row[pos])
                right_chr = sumstats["CHR"]==target_chr  
                up_pos = sumstats["POS"]>target_pos-highlight_windowkb*1000
                low_pos = sumstats["POS"]<target_pos+highlight_windowkb*1000
                sumstats.loc[right_chr&up_pos&low_pos,"HUE1"]="0"
    
    if len(highlight2)>0:
        to_highlight2 = sumstats.loc[sumstats["TCHR+POS"].isin(highlight2),:]
        if len(to_highlight2)>0:
            for i,row in to_highlight2.iterrows():
                target_chr = int(row[chrom])
                target_pos = int(row[pos])
                right_chr = sumstats["CHR"]==target_chr  
                up_pos = sumstats["POS"]>target_pos-highlight_windowkb*1000
                low_pos = sumstats["POS"]<target_pos+highlight_windowkb*1000
                sumstats.loc[right_chr&up_pos&low_pos,"HUE2"]="0"
    
    
    
    if verbose: log.write("Plotting...")
    ## figure ###########################################################################################################
    
    figargs["figsize"] = (15,10)
    fig, (ax1, ax5) = plt.subplots(2, 1, 
        gridspec_kw={'height_ratios': [1, 1]},**figargs)
    plt.subplots_adjust(hspace=region_hspace)
    
    ###########################################################################################################
    
    sumstats["s1"]=1
    sumstats["s2"]=1
    sumstats.loc[sumstats["scaled_P_1"]>-np.log10(5e-4),"s1"]=2
    sumstats.loc[sumstats["scaled_P_1"]>-np.log10(suggestive_sig_level),"s1"]=3
    sumstats.loc[sumstats["scaled_P_1"]>-np.log10(sig_level),"s1"]=4
    sumstats.loc[sumstats["scaled_P_2"]>-np.log10(5e-4),"s2"]=2
    sumstats.loc[sumstats["scaled_P_2"]>-np.log10(suggestive_sig_level),"s2"]=3
    sumstats.loc[sumstats["scaled_P_2"]>-np.log10(sig_level),"s2"]=4
    
    
    sumstats["chr_hue"]=sumstats[chrom].astype("category")
    

    ##########################################################################################################################
    maxy = max(sumstats["scaled_P_1"].max(), sumstats["scaled_P_2"].max())
    maxy1 = sumstats["scaled_P_1"].max()
    maxy5 = sumstats["scaled_P_2"].max()
    
    if cut:
        if cut is True:
            if verbose: log.write(" - Cut Auto mode is activated...")
            if maxy<20:
                if verbose: log.write(" - maxy <20 , no need to cut.")
                cut=0
            else:
                cut = 20
                cutfactor = ( maxy - cut )/5
        if cut:
            if verbose: log.write(" - Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")

            maxticker=int(max(np.round(sumstats["scaled_P_1"].max(skipna=True)) , np.round(sumstats["scaled_P_2"].max(skipna=True))))
            
            sumstats.loc[sumstats["scaled_P_1"]>cut,"scaled_P_1"] = (sumstats.loc[sumstats["scaled_P_1"]>cut,"scaled_P_1"]-cut)/cutfactor +  cut
            sumstats.loc[sumstats["scaled_P_2"]>cut,"scaled_P_2"] = (sumstats.loc[sumstats["scaled_P_2"]>cut,"scaled_P_2"]-cut)/cutfactor +  cut
            maxy = (maxticker-cut)/cutfactor + cut
            maxy1=( int(np.round(sumstats["scaled_P_1"].max(skipna=True))) -cut)/cutfactor + cut
            maxy5=( int(np.round(sumstats["scaled_P_2"].max(skipna=True))) -cut)/cutfactor + cut
    ##########################################################################################################################
    legend = None
    style=None
    linewidth=0
    
    
    palette = sns.color_palette(colors1,n_colors=sumstats[chrom].nunique())  
    plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P_1',
        hue='chr_hue',
        palette= palette,
        legend=legend,
        style=style,
        size="s1",
        sizes=marker_size,
        linewidth=linewidth,
        zorder=2,ax=ax1,edgecolor="black")     
    palette = sns.color_palette(colors2,n_colors=sumstats[chrom].nunique())  
    plot = sns.scatterplot(data=sumstats, x='i', y='scaled_P_2',
        hue='chr_hue',
        palette= palette,
        legend=legend,
        style=style,
        size="s2",
        sizes=marker_size,
        linewidth=linewidth,
        zorder=2,ax=ax5,edgecolor="black") 
    
    if len(highlight1)>0 :
        if len(to_highlight1)>0:
            if verbose: log.write(" -Highlighting target loci for sumstats1...")

            sns.scatterplot(data=sumstats.loc[sumstats["HUE1"]=="0"], x='i', y='scaled_P_1',
                   hue="HUE1",
                   palette={"0":highlight_color},
                   legend=legend,
                   style=style,
                   size="s1",
                   sizes=(marker_size[0]+1,marker_size[1]+1),
                   linewidth=linewidth,
                   zorder=3,ax=ax1,edgecolor="black",**scatter_kwargs)  
            highlight_i = sumstats.loc[sumstats["TCHR+POS"].isin(highlight),"i"].values
            
    if len(highlight2)>0 :
        if len(to_highlight2)>0:
            if verbose: log.write(" -Highlighting target loci for sumstats2.")
            sns.scatterplot(data=sumstats.loc[sumstats["HUE2"]=="0"], x='i', y='scaled_P_2',
                   hue="HUE2",
                   palette={"0":highlight_color},
                   legend=legend,
                   style=style,
                   size="s2",
                   sizes=(marker_size[0]+1,marker_size[1]+1),
                   linewidth=linewidth,
                   zorder=3,ax=ax5,edgecolor="black",**scatter_kwargs)  

            highlight_i = np.append( highlight_i, sumstats.loc[sumstats["TCHR+POS"].isin(highlight),"i"].values)
        
    if len(pinpoint1)>0 or len(pinpoint2)>0: 
        if sum(sumstats["TCHR+POS"].isin(pinpoint1))>0:  
            to_pinpoint = sumstats.loc[sumstats["TCHR+POS"].isin(pinpoint1),:]
            if verbose: log.write(" -Pinpointing target vairants...")
            ax1.scatter(to_pinpoint["i"],to_pinpoint["scaled_P_1"],color=pinpoint_color,zorder=3,s=marker_size[1]+1)
        else:
            if verbose: log.write(" -Target vairants to pinpoint were not found. Skip pinpointing process...")
                
        if sum(sumstats["TCHR+POS"].isin(pinpoint2))>0:  
            to_pinpoint = sumstats.loc[sumstats["TCHR+POS"].isin(pinpoint2),:]
            if verbose: log.write(" -Pinpointing target vairants...")
            ax5.scatter(to_pinpoint["i"],to_pinpoint["scaled_P_2"],color=pinpoint_color,zorder=3,s=marker_size[1]+1)
        else:
            if verbose: log.write(" -Target vairants to pinpoint were not found. Skip pinpointing process...")
        
        
    ## X #######################################################################################################################
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    
    ax1.set_xlabel("")
    ax1.set_xticks(chrom_df)
    ax1.set_xticklabels(chrom_df.index.map(get_number_to_chr()),fontsize=fontsize,family="sans-serif") 
    
    ax5.set_ylim(ax1.get_ylim())
    ax5.set_xlim(ax1.get_xlim())
    ax5.set_xlabel("")
    ax5.set_xticks(chrom_df)
    ax5.set_xticklabels([])
    ax5.xaxis.set_ticks_position("top")
    
    ## Y #######################################################################################################################
    sigline = ax1.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
    sigline = ax5.axhline(y=-np.log10(sig_level), linewidth = 2,linestyle="--",color=sig_line_color,zorder=1)
     
    for ax in [ax1,ax5]:
        if cut == 0: 
            ax.set_ylim(skip,maxy*1.2)
        if cut:
            cutline = ax.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
            if ((maxticker-cut)/cutfactor + cut) > cut:
                ax.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
                ax.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker],fontsize=fontsize,family="sans-serif")
            else:
                ax.set_yticks([x for x in range(skip,cut+1,2)])
                ax.set_yticklabels([x for x in range(skip,cut+1,2)],fontsize=fontsize,family="sans-serif")
            ax.set_ylim(bottom = skip)
            
     
    ### Spines ####################################################################################################################
    
    if len(anno_set1)>0 or len(anno_set2)>0:
        if len(anno_set1)>0:
            to_annotate1=sumstats.loc[sumstats["TCHR+POS"].isin(anno_set1),:]
        if len(anno_set2)>0:
            to_annotate2=sumstats.loc[sumstats["TCHR+POS"].isin(anno_set2),:]
    elif anno is not None:
        to_annotate1 = getsig(sumstats.loc[sumstats["scaled_P_1"]> float(-np.log10(sig_level)),:],
                       "TCHR+POS",
                       "CHR",
                       "POS",
                       "P_1",
                        build=build,
                        source=anno_source,
                       windowsizekb=windowsizekb,
                       verbose=False,
                       sig_level=sig_level)


        to_annotate5 = getsig(sumstats.loc[sumstats["scaled_P_2"]> float(-np.log10(sig_level)),:],
                       "TCHR+POS",
                       "CHR",
                       "POS",
                       "P_2",
                       build=build,
                        source=anno_source,
                       windowsizekb=windowsizekb,
                       verbose=False,
                       sig_level=sig_level)

        #######################################################################################
        if (to_annotate1.empty is False) and anno=="GENENAME":
                to_annotate1 = annogene(to_annotate1,
                                       id="TCHR+POS",
                                       chrom="CHR",
                                       pos="POS",
                                       log=log,
                                       build=build,
                                       source=anno_source,
                                       verbose=verbose).rename(columns={"GENE":"GENENAME"})
        if (to_annotate5.empty is False) and anno=="GENENAME":
                to_annotate5 = annogene(to_annotate5,
                                       id="TCHR+POS",
                                       chrom="CHR",
                                       pos="POS",
                                       log=log,
                                       build=build,
                                       source=anno_source,
                                       verbose=verbose).rename(columns={"GENE":"GENENAME"})
            
####################################################################################################################

# Add Annotation to manhattan plot #######################################################
           ## final
    if anno is not None:
        for index,ax,to_annotate,anno_d, anno_alias in [(0,ax1,to_annotate1,anno_d1,anno_alias1),(1,ax5,to_annotate5,anno_d2,anno_alias2)]:
            ###################### annotate() args
            fontweight = "normal"

            
            anno_default = {"rotation":40,"style":"italic","ha":"left","va":"bottom","fontsize":fontsize,"fontweight":fontweight}
            ########################
            
            if index==0:
                to_annotate["scaled_P"] = to_annotate1["scaled_P_1"]
                maxy = maxy1
            else:
                to_annotate["scaled_P"] = to_annotate5["scaled_P_2"]
                anno_default["rotation"] = -40
                anno_default["ha"]="left"
                anno_default["va"]="top"
                maxy = maxy5
            
            snpid = "TCHR+POS"

            if (anno is not None) and (to_annotate.empty!=True):
                #initiate a list for text and a starting position
                text = []
                last_pos=0
                anno_count=0

                ## log : annotation column
                if anno==True:
                        annotation_col="CHR:POS"
                elif anno is not None:
                        annotation_col = anno
                if verbose: log.write(" -Annotating using column "+annotation_col+"...")

                ##   
                for rowi,row in to_annotate.iterrows():
                    # avoid text overlapping
                    ## calculate y span
                    
                    if region is not None:
                        y_span = region[2] - region[1]
                    else:
                        y_span = sumstats["i"].max()-sumstats["i"].min()

                    ## adjust x to avoid overlapping
                    if row["i"]>last_pos+repel_force*y_span:
                        last_pos=row["i"]
                    else:
                        last_pos+=repel_force*y_span

                    if arm_scale_d is not None:
                        if anno_count not in arm_scale_d.keys():
                            arm_scale =1
                        else:
                            arm_scale = arm_scale_d[anno_count]

                    # vertical arm length in pixels
                    armB_length_in_point = plot.transData.transform((skip,1.15*maxy))[1]-plot.transData.transform((skip, row["scaled_P"]+1))[1]-arm_offset/2
                    #
                    armB_length_in_point = armB_length_in_point*arm_scale
                    
                    if arm_scale>=1:
                        armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else plot.transData.transform((skip, maxy+2))[1]-plot.transData.transform((skip,  row["scaled_P"]+1))[1] 

                    if anno_fixed_arm_length is not None:
                        anno_fixed_arm_length_factor = plot.transData.transform((skip,anno_fixed_arm_length))[1]-plot.transData.transform((skip,0))[1] 
                        armB_length_in_point = anno_fixed_arm_length_factor
                    
                    if anno==True:
                        if row[snpid] in anno_alias.keys():
                            annotation_text = anno_alias[row[snpid]]
                        else:
                            annotation_text="Chr"+ str(row[chrom]) +":"+ str(int(row[pos]))
                    elif anno == "GENENAME":
                        annotation_text=row["GENENAME"]

                    # annoatte position
                    xy=(row["i"],row["scaled_P"]+0.2)  
                    
                    xytext=(last_pos,1.15*maxy*arm_scale)
                    
                    if anno_fixed_arm_length is not None:
                        armB_length_in_point = anno_fixed_arm_length
                        xytext=(row["i"],row["scaled_P"]+0.2+anno_fixed_arm_length)
                    
                    if anno_count not in anno_d.keys():
                        if index==0:
                            arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                        else:
                            arrowargs = dict(arrowstyle="-|>",relpos=(0,1),color="#ebebeb",
                                                connectionstyle="arc,angleA=0,armA=0,angleB=-90,armB="+str(armB_length_in_point)+",rad=0")
                    
                    else:
                        xy=(row["i"],row["scaled_P"])
                        if anno_d[anno_count] in ["right","left","l","r"]:
                            if anno_d[anno_count]=="right" or anno_d[anno_count]=="r": 
                                armoffsetall = (plot.transData.transform(xytext)[0]-plot.transData.transform(xy)[0])*np.sqrt(2)
                                armoffsetb = arm_offset 
                                armoffseta = armoffsetall - armoffsetb   
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                            elif anno_d[anno_count]=="left" or anno_d[anno_count]=="l":
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str(arm_offset)+",angleB=135,armB="+str(arm_offset)+",rad=0")
                        else:
                            if anno_d[anno_count][0]=="right" or anno_d[anno_count][0]=="r": 
                                armoffsetall = (plot.transData.transform(xytext)[0]-plot.transData.transform(xy)[0])*np.sqrt(2)
                                armoffsetb = anno_d[anno_count][1] 
                                armoffseta = armoffsetall - armoffsetb   
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                            elif anno_d[anno_count]=="left" or anno_d[anno_count]=="l":
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str( anno_d[anno_count][1])+",angleB=135,armB="+str( anno_d[anno_count][1])+",rad=0")

                    bbox_para=None            
                    
                    if len(highlight) >0:
                        if row["i"] in highlight_i:
                            anno_default["fontweight"] = "bold"
                        else:
                            anno_default["fontweight"] = "normal"
                    
                    for key,value in anno_args.items():
                        anno_default[key]=value

                    ax.annotate(annotation_text,
                            xy=xy,
                            xytext=xytext,
                            bbox=bbox_para,
                            arrowprops=arrowargs,
                            zorder=100,
                            **anno_default
                            )

                    anno_count +=1


            else:
                if verbose: log.write(" -Skip annotating")
    else:
        if verbose: log.write(" -Skip annotating")

####################################################################################################################    
           
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(True)
    ax1.spines["bottom"].set_visible(True)
    
    ax5.spines["top"].set_visible(True)
    ax5.spines["right"].set_visible(False)
    ax5.spines["left"].set_visible(True)
    ax5.spines["bottom"].set_visible(False)
####################################################################################################################


    if region is not None:
        most_left_snp = sumstats["i"].idxmin()            
        gene_track_offset = sumstats.loc[most_left_snp,pos]-region[1]
        gene_track_start_i = sumstats.loc[most_left_snp,"i"] - gene_track_offset - region[1]
        lead_id_1=sumstats["scaled_P_1"].idxmax()
        lead_id_2=sumstats["scaled_P_2"].idxmax()
        lead_snp_i_1 = sumstats.loc[lead_id_1,"i"]
        lead_snp_i_2 = sumstats.loc[lead_id_2,"i"]
        
        
        region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 
        
        if region_grid is True:
            for i in np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step):
                ax1.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
                ax5.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)

        if region_lead_grid is True:
            ax1.vlines(x=lead_snp_i_1, ymin=0,  ymax=maxy1, zorder=1000,**region_lead_grid_line)
            ax5.vlines(x=lead_snp_i_2, ymin=0,  ymax=maxy5, zorder=1000,**region_lead_grid_line)
            # set x ticks m plot
        ax1.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
        ax5.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
        ax1.set_xticklabels(region_ticks,rotation=45,fontsize=fontsize,family="sans-serif")

        ax1.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])
        ax5.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])

####################################################################################################################
    # set labels
    ax1.set_xlabel("CHR",fontsize=fontsize,family="sans-serif")
    ax1.xaxis.set_label_coords(-0.02, -0.02)
    
    ax1.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
    ax5.set_ylabel("$-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
    
    ax1.set_title(titles[0],y=1+titles_pad[0])
    ax5.set_title(titles[1],y=-titles_pad[1])
    ax5.invert_yaxis() 
    
    #return fig
    if save is not None:
        if verbose: log.write("Saving plot:")
        if save==True:
            fig.savefig("./miami_plot.png",bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ "./miami_plot.png" + " successfully!" )
        else:
            fig.savefig(save,bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ save + " successfully!" )
    
    garbage_collect.collect()
    if verbose: log.write("Finished creating miami plot successfully")
# Return matplotlib figure object #######################################################################################
    return fig, log
    
    
        