import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import copy
import re
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
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

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
    region_ref_index_dic = None,
    #region_ref_second=None,
    region_grid = False,
    region_grid_line = {"linewidth": 2,"linestyle":"--"},
    region_lead_grid = True,
    region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"},
    region_title_args = None,
    region_hspace=0.02,
    region_ld_threshold = [0.2,0.4,0.6,0.8],
    region_ld_colors = ["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"],
    region_marker_shapes=None,
    palette=None,
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
        lead_ids = []
        
        for index, region_ref_single in enumerate(region_ref):
            ax1, lead_id_single = _pinpoint_lead(sumstats = sumstats,
                                        ax1 = ax1, 
                                        region_ref=region_ref_single,
                                        lead_color = palette[(index+1)*100 + len(region_ld_threshold)+2], 
                                        marker_size= marker_size,
                                        region_marker_shapes=region_marker_shapes,
                                        log=log,verbose=verbose)
            #if lead_id_single is not None:
            lead_ids.append(lead_id_single)       
        
        # update region_ref to variant rsID or variantID / skip NAs
        new_region_ref = []
        for i in range(len(lead_ids)):
            if lead_ids[i] is None:
                new_region_ref.append(region_ref[i])
                continue
            if region_ref[i] is None:
                if "SNPID" in sumstats.columns:
                    new_name = sumstats.loc[lead_ids[i],"SNPID"]
                elif "rsID" in sumstats.columns:
                    new_name = sumstats.loc[lead_ids[i],"rsID"]
                else:
                    new_name = "chr{}:{}".format(sumstats.loc[lead_ids[i],"CHR"] , sumstats.loc[lead_ids[i],"POS"])
                new_region_ref.append(new_name)
                region_ref_index_dic[new_name] = region_ref_index_dic[region_ref[i]]
                continue
            else:
                new_region_ref.append(region_ref[i])
        region_ref = new_region_ref
        ##########################################################################################################

        ##########################################################################################################

        if (vcf_path is not None) and region_ld_legend:
            ## plot cbar
            ax1, cbar = _add_ld_legend(sumstats=sumstats, 
                            ax1=ax1, 
                            region_ref=region_ref,
                            region_ld_threshold=region_ld_threshold, 
                            region_ref_index_dic=region_ref_index_dic,
                            palette=palette)
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
        
        lead_snp_ys = []
        lead_snp_is = []
        lead_snp_is_colors = []
        for i,lead_id_single in enumerate(lead_ids):
            if lead_id_single is not None:
                lead_snp_ys.append(sumstats.loc[lead_id_single,"scaled_P"] )
                lead_snp_is.append(sumstats.loc[lead_id_single,"i"])
                lead_color = palette[(region_ref_index_dic[region_ref[i]]+1)*100 + len(region_ld_threshold) +1] # consistent color
                lead_snp_is_colors.append(lead_color)

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
                        lead_snp_is=lead_snp_is,
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
                for lead_snp_i,  lead_snp_y, lead_snp_is_color in zip(lead_snp_is, lead_snp_ys , lead_snp_is_colors):
                    region_lead_grid_line["color"] = lead_snp_is_color
                    ax1.plot([lead_snp_i,lead_snp_i],[0,lead_snp_y], zorder=1,**region_lead_grid_line)
                    ax3.axvline(x=lead_snp_i, zorder=2,**region_lead_grid_line)

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
    
    return ax1, ax3, ax4, cbar, lead_snp_is, lead_snp_is_colors

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
                #try:
                matched_snpid = re.match("(chr)?[0-9]+:[0-9]+:[ATCG]+:[ATCG]+", region_ref_to_check,  re.IGNORECASE)    
                if matched_snpid is None:
                    pass
                else:
                    lead_snpid = matched_snpid.group(0).split(":")
                    if len(lead_snpid)==4:
                        lead_chr= int(lead_snpid[0])
                        lead_pos= int(lead_snpid[1])
                        lead_ea= lead_snpid[2]
                        lead_nea= lead_snpid[3]
                        chrpos_match = (sumstats["CHR"] == lead_chr) & (sumstats["POS"] == lead_pos)
                        eanea_match = ((sumstats["EA"] == lead_ea) & (sumstats["NEA"] == lead_nea)) | ((sumstats["EA"] == lead_nea) & (sumstats["NEA"] == lead_ea)) 
                        if "rsID" in sumstats.columns:
                            lead_id = sumstats.index[chrpos_match&eanea_match].to_list()
                        if "SNPID" in sumstats.columns:
                            lead_id = sumstats.index[chrpos_match&eanea_match].to_list()     
                if type(lead_id) is list:
                    if len(lead_id)>0:
                        lead_id = int(lead_id[0])   
                        log.warning("Trying matching variant {} using CHR:POS:EA:NEA to {}... ".format(region_ref_to_check,lead_id))

        if type(lead_id) is list:
            if len(lead_id)==0 :
                log.warning("Extracting variant: {} not found in sumstats.. Skipping..".format(region_ref_to_check))
                #lead_id = sumstats["scaled_P"].idxmax()
                lead_id = None
                return lead_id
        else:
            log.write(" -Reference variant ID: {} - {}".format(region_ref_to_check, lead_id), verbose=verbose)

    if lead_id is None:
        log.write(" -Extracting lead variant...", verbose=verbose)
        lead_id = sumstats["scaled_P"].idxmax()

    return lead_id

def _pinpoint_lead(sumstats,ax1,region_ref, lead_color, marker_size, log, verbose, region_marker_shapes):
    
    if region_ref is None:
        log.write(" -Extracting lead variant..." , verbose=verbose)
        lead_id = sumstats["scaled_P"].idxmax()
    else:
        lead_id = _get_lead_id(sumstats, region_ref, log, verbose)
    
    if lead_id is not None:
        ax1.scatter(sumstats.loc[lead_id,"i"],sumstats.loc[lead_id,"scaled_P"],
                color=lead_color,
                zorder=3,
                marker= region_marker_shapes[sumstats.loc[lead_id,"SHAPE"]-1],
                s=marker_size[1]+2,
                edgecolor="black")

    return ax1, lead_id
# -############################################################################################################################################################################
def _add_region_title(region_title, ax1,region_title_args):
    ax1.text(0.015,0.97, region_title, transform=ax1.transAxes, va="top", ha="left", region_ref=None, **region_title_args )
    return ax1

def _add_ld_legend(sumstats, ax1, region_ld_threshold, region_ref,region_ref_index_dic,palette =None, position=1):
    
    width_pct = "11%"
    height_pct = "{}%".format( 14 + 7 * len(region_ref))
    axins1 = inset_axes(ax1,
            width=width_pct,  # width = 50% of parent_bbox width
            height=height_pct,  # height : 5%
            loc='upper right',axes_kwargs={"frameon":True,"facecolor":"white","zorder":999999})
    
    ld_ticks = [0]+region_ld_threshold+[1]

    for index, ld_threshold in enumerate(ld_ticks):
        for group_index in range(len(region_ref)):
            if index < len(ld_ticks)-1:            
                x=ld_threshold
                y=0.2*group_index
                width=0.2
                height=ld_ticks[index+1]-ld_ticks[index]
                hex_color = palette[(region_ref_index_dic[region_ref[group_index]]+1)*100 + index+1] # consistent color
                
                a = Rectangle((x,y),width, height, fill = True, color = hex_color , linewidth = 2)
                #patches.append(a)
                axins1.add_patch(a)

    # y snpid
    yticks_position = 0.1 + 0.2 *np.arange(0,len(region_ref))
    axins1.set_yticks(yticks_position, ["{}".format(x) for x in region_ref])
    axins1.set_ylim(0,0.2*len(region_ref))    
    
    # x ld thresholds
    axins1.set_xticks(ticks=ld_ticks) 
    axins1.set_xticklabels([str(i) for i in ld_ticks]) 
    axins1.set_xlim(0,1)    

    axins1.set_aspect('equal', adjustable='box')
    axins1.set_title('LD $r^{2}$ with variant',loc="center",y=-0.2)
    cbar = axins1
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
    lead_snp_is,
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

    
    sig_gene_names=[]
    sig_gene_lefts=[]
    sig_gene_rights=[]
    for index,row in uniq_gene_region.iterrows():

        gene_color="#020080"
        #if row[6][0]=="+":
        if row["strand"][0]=="+":
            gene_anno = row["name"] + "->"
        else:
            gene_anno = "<-" + row["name"] 
        


        for lead_snp_i in lead_snp_is:
            if region_lead_grid is True and lead_snp_i > gene_track_start_i+row["start"] and lead_snp_i < gene_track_start_i+row["end"] :
                gene_color=region_lead_grid_line["color"]
                sig_gene_names.append(row["name"])
                sig_gene_lefts.append(gene_track_start_i+row["start"])
                sig_gene_rights.append(gene_track_start_i+row["end"])

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
        for sig_gene_name, sig_gene_left, sig_gene_right in zip(sig_gene_names,sig_gene_lefts,sig_gene_rights):
            
            if not pd.isnull(row["name"]):
                if (region_lead_grid is True) and row["name"]==sig_gene_name:
                    exon_color = region_lead_grid_line["color"]  
                else:
                    exon_color="#020080" 
            elif gene_track_start_i+row["starts"] > sig_gene_left and gene_track_start_i+row["end"] < sig_gene_right:
                exon_color = region_lead_grid_line["color"]  
            else:
                exon_color="#020080"
        
        ax3.plot((gene_track_start_i+row["start"],gene_track_start_i+row["end"]),
                    (row["stack"]*2,row["stack"]*2),linewidth=linewidth_in_points*taf[3],color=exon_color,solid_capstyle="butt")

    log.write(" -Finished plotting gene track..", verbose=verbose)

    return ax3,texts_to_adjust_middle

# -############################################################################################################################################################################
# Helpers    
# -############################################################################################################################################################################
def process_vcf(sumstats, 
                vcf_path, 
                region,
                region_ref, 
                #region_ref_second, 
                log, 
                verbose, 
                pos ,
                nea,
                ea, 
                region_ld_threshold, 
                vcf_chr_dict,
                tabix):

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
    #############################################################################################
    sumstats["REFINDEX"] = sumstats[[pos,nea,ea]].apply(lambda x: match_varaint(x),axis=1)
    #############################################################################################

    #for loop to add LD information
    #############################################################################################
    for ref_n, region_ref_single in enumerate(region_ref):

        rsq = "RSQ_{}".format(ref_n)
        ld_single = "LD_{}".format(ref_n)
        lead = "LEAD_{}".format(ref_n)
        sumstats[lead]= 0

        # get lead variant id and pos
        if region_ref_single is None:
            # if not specified, use lead variant
            lead_id = sumstats["scaled_P"].idxmax()
        else:
            # figure out lead variant
            lead_id = _get_lead_id(sumstats, region_ref_single, log, verbose)
        

        lead_series = None
        if lead_id is None:
            
            matched_snpid = re.match("(chr)?[0-9]+:[0-9]+:[ATCG]+:[ATCG]+",region_ref_single,  re.IGNORECASE)
            
            if matched_snpid is None:
                sumstats[rsq] = None
                sumstats[rsq] = sumstats[rsq].astype("float")
                sumstats[ld_single] = 0    
                continue    
            else:
                
                lead_snpid = matched_snpid.group(0).split(":")[1:]
                lead_pos = int(lead_snpid[0])
                lead_snpid[0]= int(lead_snpid[0])
                lead_series = pd.Series(lead_snpid)
        else:
            lead_pos = sumstats.loc[lead_id,pos]

        
        # if lead pos is available: 
        if lead_pos in ref_genotype["variants/POS"]:
            
            # get ref index for lead snp
            if lead_series is None:
                lead_snp_ref_index = match_varaint(sumstats.loc[lead_id,[pos,nea,ea]])
                #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]
            else:
                log.warning("Computing LD: {} not found in sumstats but found in reference...Still Computing...".format(region_ref_single))
                lead_snp_ref_index = match_varaint(lead_series)

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
            sumstats.loc[~sumstats["REFINDEX"].isna(),rsq] = valid_r2
        else:
            log.write(" -Lead SNP not found in reference...", verbose=verbose)
            sumstats[rsq]=None
        
        sumstats[rsq] = sumstats[rsq].astype("float")
        sumstats[ld_single] = 0
        
        for index,ld_threshold in enumerate(region_ld_threshold):
            # No data,LD = 0
            # 0, 0.2  LD = 1
            # 1, 0.4  LD = 2
            # 2, 0.6  LD = 3
            # 3, 0.8  LD = 4
            # 4, 1.0  LD = 5
            # lead    LD = 6

            if index==0:
                to_change_color = sumstats[rsq]>-1
                sumstats.loc[to_change_color,ld_single] = 1
            to_change_color = sumstats[rsq]>ld_threshold
            sumstats.loc[to_change_color,ld_single] = index+2
        
        if lead_series is None:
            sumstats.loc[lead_id,ld_single] = len(region_ld_threshold)+2
            sumstats.loc[lead_id,lead] = 1

    ####################################################################################################    
    final_shape_col = "SHAPE"
    final_ld_col = "LD"
    final_rsq_col = "RSQ"

    sumstats[final_ld_col]  = 0
    sumstats[final_shape_col] = 1
    sumstats[final_rsq_col] = 0.0
    
    for i in range(len(region_ref)):
        ld_single = "LD_{}".format(i)
        current_rsq = "RSQ_{}".format(i)
        a_ngt_b = sumstats[final_rsq_col] < sumstats[current_rsq]
        #set levels with interval=100
        sumstats.loc[a_ngt_b, final_ld_col] = 100 * (i+1) + sumstats.loc[a_ngt_b, ld_single]
        sumstats.loc[a_ngt_b, final_rsq_col] = sumstats.loc[a_ngt_b, current_rsq]
        sumstats.loc[a_ngt_b, final_shape_col] = i + 1

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
