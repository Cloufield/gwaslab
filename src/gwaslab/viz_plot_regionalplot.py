import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import copy
import scipy as sp
from pyensembl import EnsemblRelease
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from adjustText import adjust_text
from gtfparse import read_gtf
from gwaslab.g_Log import Log
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_recombination_rate
from gwaslab.bd_common_data import get_gtf

def _plot_regional(
    sumstats,
    fig,
    ax1,
    ax3,
    region,
    vcf_path,
    marker_size,
    fontsize,
    build,
    chrom_df,
    xtick_chr_dict,
    cut_line_color,
    vcf_chr_dict = None,
    gtf_path="default",
    gtf_chr_dict = get_number_to_chr(),
    gtf_gene_name=None,
    rr_path="default",
    rr_header_dict=None,
    rr_chr_dict = get_number_to_chr(),
    rr_lim = (0,100),
    rr_ylabel = True,
    rr_title=None,
    region_ld_legend=True,
    region_title=None,
    mode="mqq",
    region_step = 21,
    region_ref=None,
    region_ref_second=None,
    region_grid = False,
    region_grid_line = {"linewidth": 2,"linestyle":"--"},
    region_lead_grid = True,
    region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"},
    region_title_args = None,
    region_hspace=0.02,
    region_ld_threshold = [0.2,0.4,0.6,0.8],
    region_ld_colors = ["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"],
    region_ld_colors1=None,
    region_ld_colors2=None,
    region_recombination = True,
    region_protein_coding=True,
    region_flank_factor = 0.05,
    track_font_family="Arial",
    taf=[4,0,0.95,1,1],
    # track_n, track_n_offset,font_ratio,exon_ratio,text_offset
    tabix=None,
    chrom="CHR",
    pos="POS",
    verbose=True,
    log=Log()
):  

    # x axix: use i to plot (there is a gap between i and pos) 

    # if regional plot : pinpoint lead , add color bar ##################################################
    if (region is not None) :
        # pinpoint lead
        if region_ref_second is None:
            ax1, lead_id = _pinpoint_lead(sumstats = sumstats,
                                        ax1 = ax1, 
                                        region_ref=region_ref,
                                        region_ld_threshold = region_ld_threshold, 
                                        region_ld_colors = region_ld_colors, 
                                        marker_size= marker_size,
                                        log=log,verbose=verbose)
        else:
            ax1, lead_id = _pinpoint_lead(sumstats = sumstats,
                                            ax1 = ax1, 
                                            region_ref=region_ref,
                                            region_ld_threshold = region_ld_threshold, 
                                            region_ld_colors = region_ld_colors1, 
                                            marker_size= marker_size,
                                            log=log,verbose=verbose)
            ax1, lead_id2 = _pinpoint_lead(sumstats = sumstats,
                                            ax1 = ax1, 
                                            region_ref=region_ref_second,
                                            region_ld_threshold = region_ld_threshold, 
                                            region_ld_colors = region_ld_colors2, 
                                            marker_size= marker_size,
                                            log=log,verbose=verbose)

        if (vcf_path is not None) and region_ld_legend:
            if region_ref_second is None:
                ax1, cbar = _add_ld_legend(sumstats=sumstats, 
                                ax1=ax1, 
                                region_ld_threshold=region_ld_threshold, 
                                region_ld_colors=region_ld_colors)
            else:
                
                ax1, cbar1 = _add_ld_legend(sumstats=sumstats, 
                                ax1=ax1, 
                                region_ld_threshold=region_ld_threshold, 
                                region_ld_colors=region_ld_colors1,
                                position=1)
                ax1, cbar2 = _add_ld_legend(sumstats=sumstats, 
                                ax1=ax1, 
                                region_ld_threshold=region_ld_threshold, 
                                region_ld_colors=region_ld_colors2,
                                position=2)
                cbar = [cbar1, cbar2]
        else:
            cbar=None
        if region_title is not None:
                ax1 = _add_region_title(region_title, ax1=ax1,region_title_args=region_title_args )
    ## recombinnation rate ##################################################       
    if (region is not None) and (region_recombination is True):
        ax4 = _plot_recombination_rate(sumstats = sumstats,
                                        pos =pos,
                                        region= region, 
                                        ax1 = ax1, 
                                        rr_path =rr_path, 
                                        rr_chr_dict = rr_chr_dict, 
                                        rr_header_dict =rr_header_dict, 
                                        build= build,
                                        rr_lim=rr_lim,
                                        rr_ylabel=rr_ylabel)
     
    ## regional plot : gene track ######################################################################
                # calculate offset
    if (region is not None):
        most_left_snp      = sumstats["i"].idxmin()
        
        # distance between leftmost variant position to region left bound
        gene_track_offset  = sumstats.loc[most_left_snp,pos] - region[1]
        
        # rebase i to region[1] : the i value when POS=0
        gene_track_start_i = sumstats.loc[most_left_snp,"i"] - gene_track_offset - region[1]
        

        #lead_id            = sumstats["scaled_P"].idxmax()
        lead_snp_y         = sumstats.loc[lead_id,"scaled_P"] 
        #sumstats["scaled_P"].max()
        lead_snp_i         = sumstats.loc[lead_id,"i"]
        
        if region_ref_second is not None:
            lead_snp_y2    = sumstats.loc[lead_id2,"scaled_P"] 
            #sumstats["scaled_P"].max()
            lead_snp_i2    = sumstats.loc[lead_id2,"i"]
        else:
            lead_snp_i2 = None

    if gtf_path is not None:
        # load gtf
        ax3, texts_to_adjust_middle =_plot_gene_track(
                        ax3=ax3,
                        fig=fig,
                        gtf_path=gtf_path,
                        region=region,
                        region_flank_factor=region_flank_factor,
                        region_protein_coding=region_protein_coding,
                        region_lead_grid=region_lead_grid,
                        region_lead_grid_line=region_lead_grid_line,
                        lead_snp_i=lead_snp_i,
                        lead_snp_i2 = lead_snp_i2,
                        region_ld_colors1=region_ld_colors1,
                        region_ld_colors2=region_ld_colors2,
                        gene_track_start_i=gene_track_start_i,
                        gtf_chr_dict=gtf_chr_dict,
                        gtf_gene_name=gtf_gene_name, 
                        track_font_family=track_font_family,
                        taf=taf,
                        build=build, 
                        verbose=verbose, 
                        log=log)

    ## regional plot - set X tick
    if region is not None:
        region_ticks = list(map('{:.3f}'.format,np.linspace(region[1], region[2], num=region_step).astype("int")/1000000)) 
        
        # set x ticks for gene track
        if "r" in mode:
            if gtf_path is not None: 
                ax3.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
                ax3.set_xticklabels(region_ticks,rotation=45,fontsize=fontsize,family="sans-serif")
            if region_grid==True:
                for i in np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step):
                    ax1.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
                    ax3.axvline(x=i, color=cut_line_color,zorder=1,**region_grid_line)
            if region_lead_grid==True:
                if region_ref_second is None:
                    ax1.plot([lead_snp_i,lead_snp_i],[0,lead_snp_y], zorder=1,**region_lead_grid_line)
                    ax3.axvline(x=lead_snp_i, zorder=2,**region_lead_grid_line)
                    #ax3.plot([lead_snp_i,lead_snp_i],[0,1], zorder=1,transform=ax3.transAxes, **region_lead_grid_line)
                else:
                    region_lead_grid_line["color"] = region_ld_colors1[-1]
                    ax1.plot([lead_snp_i,lead_snp_i],[0,lead_snp_y], zorder=1,**region_lead_grid_line)
                    ax3.axvline(x=lead_snp_i, zorder=1,**region_lead_grid_line)
                    
                    region_lead_grid_line["color"] = region_ld_colors2[-1]
                    ax1.plot([lead_snp_i2,lead_snp_i2],[0,lead_snp_y2], zorder=1,**region_lead_grid_line)
                    ax3.axvline(x=lead_snp_i2, zorder=1,**region_lead_grid_line)
        else:
            # set x ticks m plot
            ax1.set_xticks(np.linspace(gene_track_start_i+region[1], gene_track_start_i+region[2], num=region_step))
            ax1.set_xticklabels(region_ticks,rotation=45,fontsize=fontsize,family="sans-serif")
    
        ax1.set_xlim([gene_track_start_i+region[1], gene_track_start_i+region[2]])

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
    
    return ax1, ax3, ax4, cbar, lead_snp_i, lead_snp_i2

# + ###########################################################################################################################################################################
def _get_lead_id(sumstats=None, region_ref=None, log=None, verbose=True):
    region_ref_to_check = copy.copy(region_ref)
    try: 
        if len(region_ref_to_check)>0 and type(region_ref_to_check) is not str:
            region_ref_to_check = region_ref_to_check[0]
    except:
        pass
    lead_id=None
    if "rsID" in sumstats.columns:
        lead_id = sumstats.index[sumstats["rsID"] == region_ref_to_check].to_list()

    if lead_id is None and "SNPID" in sumstats.columns:
        lead_id = sumstats.index[sumstats["SNPID"] == region_ref_to_check].to_list()

    if type(lead_id) is list:
        if len(lead_id)>0:
            lead_id = int(lead_id[0])
    if region_ref_to_check is not None:
        if type(lead_id) is list:
            if len(lead_id)==0 :
                log.warning("{} not found. Roll back to lead variant...".format(region_ref_to_check))
                lead_id = sumstats["scaled_P"].idxmax()
        else:
            log.write(" -Reference variant ID: {} - {}".format(region_ref_to_check, lead_id))

    if lead_id is None:
        log.write(" -Extracting lead variant...", verbose=verbose)
        lead_id = sumstats["scaled_P"].idxmax()

    return lead_id

def _pinpoint_lead(sumstats,ax1,region_ref, region_ld_threshold, region_ld_colors, marker_size, log, verbose):
    if region_ref is None:
        log.write(" -Extracting lead variant..." , verbose=verbose)
        lead_id = sumstats["scaled_P"].idxmax()
    else:
        lead_id = _get_lead_id(sumstats, region_ref, log, verbose)
        
    ax1.scatter(sumstats.loc[lead_id,"i"],sumstats.loc[lead_id,"scaled_P"],
            color=region_ld_colors[-1],
            marker="D",
            zorder=3,
            s=marker_size[1]+2,
            edgecolor="black")

    return ax1, lead_id
# -############################################################################################################################################################################
def _add_region_title(region_title, ax1,region_title_args):
    ax1.text(0.015,0.97, region_title, transform=ax1.transAxes, va="top", ha="left", **region_title_args )
    return ax1

def _add_ld_legend(sumstats, ax1, region_ld_threshold, region_ld_colors,position=1):
    # add a in-axis axis for colorbar
    if position==1:
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
    elif position==2:
        axins1 = inset_axes(ax1,
                width="5%",  # width = 50% of parent_bbox width
                height="40%",  # height : 5%
                loc='upper left',axes_kwargs={"frameon":True,"facecolor":"white","zorder":999999})

        cmp= mpl.colors.ListedColormap(region_ld_colors[1:])
        norm= mpl.colors.BoundaryNorm([0]+region_ld_threshold+[1],cmp.N)     
        cbar= plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmp),
                            cax=axins1,
                            fraction=1,
                            orientation="vertical",
                            ticklocation="right")
        cbar.set_ticks(ticks=region_ld_threshold) 
        cbar.set_ticklabels([str(i) for i in region_ld_threshold]) 
        cbar.ax.set_title('LD $r^{2}$',loc="center",y=-0.2)
        rect = patches.Rectangle((0.0,0.50),
                                height=0.5,
                                width=0.09, 
                                transform=ax1.transAxes,
                                linewidth=1, 
                                edgecolor='black', 
                                facecolor='white',
                                zorder=999998)
    ax1.add_patch(rect)
    return ax1, cbar

# -############################################################################################################################################################################
def  _plot_recombination_rate(sumstats,pos, region, ax1, rr_path, rr_chr_dict, rr_header_dict, build,rr_lim,rr_ylabel=True):
    ax4 = ax1.twinx()
    most_left_snp = sumstats["i"].idxmin()
    
    # the i value when pos=0
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
    
    ax1.set_zorder(ax4.get_zorder()+1)
    ax1.patch.set_visible(False)
    
    if rr_ylabel:
        ax4.set_ylabel("Recombination rate(cM/Mb)")
    if rr_lim!="max":
        ax4.set_ylim(rr_lim[0],rr_lim[1])
    else:
        ax4.set_ylim(0, 1.05 * rc["Rate(cM/Mb)"].max())
    ax4.spines["top"].set_visible(False)
    ax4.spines["top"].set(zorder=1) 
    return ax4

# -############################################################################################################################################################################
def _plot_gene_track(
    ax3,
    fig,
    gtf_path,
    region,
    region_flank_factor,
    region_protein_coding,
    region_lead_grid,
    region_lead_grid_line,
    lead_snp_i,
    lead_snp_i2,
    region_ld_colors1,
    region_ld_colors2,
    gene_track_start_i,
    gtf_chr_dict,gtf_gene_name, 
    track_font_family,
    taf,
    build, 
    verbose=True, 
    log=Log()):
    
    # load gtf
    log.write(" -Loading gtf files from:" + gtf_path, verbose=verbose)
    uniq_gene_region,exons = process_gtf(   gtf_path = gtf_path ,
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
    log.write(" -plotting gene track..", verbose=verbose)
    
    sig_gene_name = "Undefined"
    sig_gene_name2 = "Undefined"
    texts_to_adjust_left = []
    texts_to_adjust_middle = []
    texts_to_adjust_right = []
    for index,row in uniq_gene_region.iterrows():

        gene_color="#020080"
        #if row[6][0]=="+":
        if row["strand"][0]=="+":
            gene_anno = row["name"] + "->"
        else:
            gene_anno = "<-" + row["name"] 
        
        if region_lead_grid is True and lead_snp_i > gene_track_start_i+row["start"] and lead_snp_i < gene_track_start_i+row["end"] :
                gene_color=region_lead_grid_line["color"]
                sig_gene_name = row["name"]
                sig_gene_left = gene_track_start_i+row["start"]
                sig_gene_right= gene_track_start_i+row["end"]

        # plot gene line
        ax3.plot((gene_track_start_i+row["start"],gene_track_start_i+row["end"]),
                    (row["stack"]*2,row["stack"]*2),color=gene_color,linewidth=linewidth_in_points/10)
        
        if lead_snp_i2 is not None:
            if region_lead_grid is True and lead_snp_i2 > gene_track_start_i+row["start"] and lead_snp_i2 < gene_track_start_i+row["end"] :
                    gene_color=region_lead_grid_line["color"]
                    sig_gene_name2 = row["name"]
                    sig_gene_left2 = gene_track_start_i+row["start"]
                    sig_gene_right2= gene_track_start_i+row["end"]

            # plot gene line
            ax3.plot((gene_track_start_i+row["start"],gene_track_start_i+row["end"]),
                    (row["stack"]*2,row["stack"]*2),color=gene_color,linewidth=linewidth_in_points/10)
        
        

        
        # plot gene name
        if row["end"] >= region[2]:
            #right side
            texts_to_adjust_right.append(ax3.text(x=gene_track_start_i+region[2],
                    y=row["stack"]*2+taf[4],s=gene_anno,ha="right",va="center",color="black",style='italic', size=font_size_in_points,family=track_font_family))
            
        elif row["start"] <= region[1] :
            #left side
            texts_to_adjust_left.append(ax3.text(x=gene_track_start_i+region[1],
                    y=row["stack"]*2+taf[4],s=gene_anno,ha="left",va="center",color="black",style='italic', size=font_size_in_points,family=track_font_family))
        else:
            texts_to_adjust_middle.append(ax3.text(x=(gene_track_start_i+row["start"]+gene_track_start_i+row["end"])/2,
                    y=row["stack"]*2+taf[4],s=gene_anno,ha="center",va="center",color="black",style='italic',size=font_size_in_points,family=track_font_family))
    
    # plot exons
    for index,row in exons.iterrows():
        
        exon_color="#020080" 
        if not pd.isnull(row["name"]):
            if (region_lead_grid is True) and row["name"]==sig_gene_name:
                exon_color = region_lead_grid_line["color"]  
            else:
                exon_color="#020080" 
        elif gene_track_start_i+row["starts"] > sig_gene_left and gene_track_start_i+row["end"] < sig_gene_right:
            exon_color = region_lead_grid_line["color"]  
        else:
            exon_color="#020080"

        if lead_snp_i2 is not None:
            if not pd.isnull(row["name"]):
                if (region_lead_grid is True) and row["name"]==sig_gene_name2:
                    exon_color = region_lead_grid_line["color"] 
            elif gene_track_start_i+row["starts"] > sig_gene_left2 and gene_track_start_i+row["end"] < sig_gene_right2:
                exon_color = region_lead_grid_line["color"]  

        
        #if row["gene_biotype"]!="protein_coding":
        #    exon_color="grey"
        ax3.plot((gene_track_start_i+row["start"],gene_track_start_i+row["end"]),
                    (row["stack"]*2,row["stack"]*2),linewidth=linewidth_in_points*taf[3],color=exon_color,solid_capstyle="butt")

    log.write(" -Finished plotting gene track..", verbose=verbose)

    return ax3,texts_to_adjust_middle

# -############################################################################################################################################################################
# Helpers    
# -############################################################################################################################################################################
def process_vcf(sumstats, vcf_path, region,region_ref, region_ref_second, log, verbose, pos ,nea,ea, region_ld_threshold, vcf_chr_dict,tabix):
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)



    # load genotype data of the targeted region
    ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
    if ref_genotype is None:
        log.warning("No data was retrieved. Skipping ...")
        ref_genotype=dict()
        ref_genotype["variants/POS"]=np.array([],dtype="int64")
    log.write(" -Retrieving index...", verbose=verbose)
    log.write(" -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
    # match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
     #######################################################################################
    def match_varaint(x):
        # x: "POS,NEA,EA"
        if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
            # position match
            if len(np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0])>1:  
            # multiple position matches
                for j in np.where(ref_genotype["variants/POS"] == x.iloc[0])[0]:
                # for each possible match, compare ref and alt
                    if x.iloc[1] == ref_genotype["variants/REF"][j]:
                        if x.iloc[2] in ref_genotype["variants/ALT"][j]:
                            return j
                    elif x.iloc[1] in ref_genotype["variants/ALT"][j]:
                        if x.iloc[2] == ref_genotype["variants/REF"][j]:
                            return j    
                return None
            else: 
                # single match
                return np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0][0] 
        else:
            # no position match
            return None
    log.write(" -Matching variants using POS, NEA, EA ...", verbose=verbose)
    sumstats["REFINDEX"] = sumstats[[pos,nea,ea]].apply(lambda x: match_varaint(x),axis=1)
    #############################################################################################
    #sumstats["REFINDEX"] = sumstats[pos].apply(lambda x: np.where(ref_genotype["variants/POS"] == x )[0][0] if np.any(ref_genotype["variants/POS"] == x) else None)
    
    # get lead variant id and pos
    #lead_id = sumstats["scaled_P"].idxmax()
    if region_ref is None:
        lead_id = sumstats["scaled_P"].idxmax()
    else:
        lead_id = _get_lead_id(sumstats, region_ref, log, verbose)
    lead_pos = sumstats.loc[lead_id,pos]
    
    # if lead pos is available: 
    if lead_pos in ref_genotype["variants/POS"]:
        
        # get ref index for lead snp
        lead_snp_ref_index = match_varaint(sumstats.loc[lead_id,[pos,nea,ea]])
        #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

        # non-na other snp index
        other_snps_ref_index = sumstats["REFINDEX"].dropna().astype("int").values
        # get genotype 
        lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
        try:
            if len(set(lead_snp_genotype[0]))==1:
                log.warning("The variant is mono-allelic in reference VCF. LD can not be calculated.")
        except:
            pass
        other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
        
        log.write(" -Calculating Rsq...", verbose=verbose)
        
        if len(other_snp_genotype)>1:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
        else:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)
        sumstats.loc[~sumstats["REFINDEX"].isna(),"RSQ"] = valid_r2
    else:
        log.write(" -Lead SNP not found in reference...", verbose=verbose)
        sumstats["RSQ"]=None
    
    sumstats["RSQ"] = sumstats["RSQ"].astype("float")
    sumstats["LD"] = 0
    
    for index,ld_threshold in enumerate(region_ld_threshold):
        # No data,LD = 0
        # 0, 0.2  LD = 1
        # 1, 0.4  LD = 2
        # 2, 0.6  LD = 3
        # 3, 0.8  LD = 4
        # 4, 1.0  LD = 5
        # lead    LD = 6

        if index==0:
            to_change_color = sumstats["RSQ"]>-1
            sumstats.loc[to_change_color,"LD"] = 1
        to_change_color = sumstats["RSQ"]>ld_threshold
        sumstats.loc[to_change_color,"LD"] = index+2

    sumstats.loc[lead_id,"LD"] = len(region_ld_threshold)+2
    sumstats["LEAD"]="Other variants"
    sumstats.loc[lead_id,"LEAD"] = "Lead variants"
    sumstats["SHAPE"] = 1
    #####################################################################################################
    if region_ref_second is not None:

        lead_id2 = _get_lead_id(sumstats, region_ref_second, log, verbose)

        lead_pos2 = sumstats.loc[lead_id2,pos]
        if lead_pos2 in ref_genotype["variants/POS"]:
            
            # get ref index for lead snp
            lead_snp_ref_index = match_varaint(sumstats.loc[lead_id2,[pos,nea,ea]])
            #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

            # non-na other snp index
            other_snps_ref_index = sumstats["REFINDEX"].dropna().astype("int").values
            # get genotype 
            lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
            try:
                if len(set(lead_snp_genotype[0]))==1:
                    log.warning("The variant is mono-allelic in reference VCF. LD can not be calculated.")
            except:
                pass
            other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
            
            log.write(" -Calculating Rsq...", verbose=verbose)
            
            if len(other_snp_genotype)>1:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
            else:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)
            sumstats.loc[~sumstats["REFINDEX"].isna(),"RSQ2"] = valid_r2
        else:
            log.write(" -Lead SNP not found in reference...", verbose=verbose)
            sumstats["RSQ2"]=None

        sumstats["RSQ2"] = sumstats["RSQ2"].astype("float")
        sumstats["LD2"] = 0

        for index,ld_threshold in enumerate(region_ld_threshold):
            if index==0:
                to_change_color = sumstats["RSQ2"]>-1
                sumstats.loc[to_change_color,"LD2"] = 1
            #else:
            #    to_change_color = sumstats["RSQ2"]>ld_threshold
            #    sumstats.loc[to_change_color,"LD2"] = index+1
            to_change_color = sumstats["RSQ2"]>ld_threshold
            sumstats.loc[to_change_color,"LD2"] = index+2
            
        sumstats.loc[lead_id,"LD2"] = len(region_ld_threshold)+2
        sumstats["LEAD2"]="Other variants"
        sumstats.loc[lead_id,"LEAD2"] = "Lead variants"

        a_ngt_b = sumstats["RSQ"] < sumstats["RSQ2"]
        sumstats.loc[a_ngt_b, "LD"] = 100 + sumstats.loc[a_ngt_b, "LD2"]
        sumstats.loc[a_ngt_b, "SHAPE"] =2 
        #sumstats.loc[lead_id,"LEAD2"]
    ####################################################################################################

    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return sumstats

# -############################################################################################################################################################################

def process_gtf(gtf_path,
                region,
                region_flank_factor,
                build,
                region_protein_coding,
                gtf_chr_dict,
                gtf_gene_name):
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
        #gtf = pd.read_csv(gtf_path,sep="\t",header=None, comment="#", low_memory=False,dtype={0:"string"})
        gtf = read_gtf(gtf_path)
        gtf = gtf.loc[gtf["seqname"]==gtf_chr_dict[region[0]],:]

    # filter in region
    genes_1mb = gtf.loc[(gtf["seqname"]==to_query_chrom)&(gtf["start"]<region[2])&(gtf["end"]>region[1]),:].copy()
    
    # extract biotype
    #genes_1mb.loc[:,"gene_biotype"] = genes_1mb[8].str.extract(r'gene_biotype "([\w\.\_-]+)"')
    
    # extract gene name
    if gtf_gene_name is None:
        if gtf_path=="refseq":
            #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"').astype("string")
            genes_1mb.loc[:,"name"] = genes_1mb["gene_id"]
        elif gtf_path =="default" or gtf_path =="ensembl":
            #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_name "([\w\.-]+)"').astype("string")
            genes_1mb.loc[:,"name"] = genes_1mb["gene_name"]
        else:
            #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"').astype("string")
            genes_1mb.loc[:,"name"] = genes_1mb["gene_id"]
    else:
        #pattern = r'{} "([\w\.-]+)"'.format(gtf_gene_name)
        #genes_1mb.loc[:,"name"] = genes_1mb[8].str.extract(pattern).astype("string")
        genes_1mb.loc[:,"name"] = genes_1mb[gtf_gene_name]

    # extract gene
    #genes_1mb.loc[:,"gene"] = genes_1mb[8].str.extract(r'gene_id "([\w\.-]+)"')
    genes_1mb.loc[:,"gene"] = genes_1mb["gene_id"]
    

    # extract protein coding gene
    if region_protein_coding is True:
        #genes_1mb  =  genes_1mb.loc[genes_1mb["gene_biotype"]=="protein_coding",:].copy()
        pc_genes_1mb_list = genes_1mb.loc[(genes_1mb["feature"]=="gene")& (genes_1mb["gene_biotype"]=="protein_coding") & (genes_1mb["name"]!=""),"name"].values
        genes_1mb = genes_1mb.loc[(genes_1mb["feature"].isin(["exon","gene"])) & (genes_1mb["name"].isin(pc_genes_1mb_list)),:]
    # extract exon
    exons = genes_1mb.loc[genes_1mb["feature"]=="exon",:].copy()
    
    #uniq genes
    ## get all record with 2nd column == gene
    #uniq_gene_region = genes_1mb.loc[genes_1mb[2]=="gene",:].copy()
    uniq_gene_region = genes_1mb.loc[genes_1mb["feature"]=="gene",:].copy()
    
    ## extract region + flank
    flank = region_flank_factor * (region[2] - region[1])
    
    ## get left and right boundary
    #uniq_gene_region["left"] = uniq_gene_region[3]-flank
    #uniq_gene_region["right"] = uniq_gene_region[4]+flank
    #
    uniq_gene_region["left"] = uniq_gene_region["start"]-flank
    uniq_gene_region["right"] = uniq_gene_region["end"]+flank

    # arrange gene track
    stack_dic = assign_stack(uniq_gene_region.sort_values(["start"]).loc[:,["name","left","right"]])  

    # map gene to stack and add stack column : minus stack
    uniq_gene_region["stack"] = -uniq_gene_region["name"].map(stack_dic)
    exons.loc[:,"stack"] = -exons.loc[:,"name"].map(stack_dic)

    # return uniq_gene_region (gene records with left and right boundary)
    # return exon records with stack number
    return uniq_gene_region, exons


# -############################################################################################################################################################################
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
