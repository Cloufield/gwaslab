import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import matplotlib as mpl
from scipy import stats
"""
Annotate significant variants in a single Manhattan/QQ plot with customizable styling and positioning options.
"""

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from adjustText import adjust_text
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_reposition_text import adjust_text_position
from pandas.api.types import is_string_dtype

# single mqqplot
def annotate_single(
    sumstats,
    anno,
    mode,
    ax1,
    highlight_i,
    highlight_chrpos,
    highlight_anno_kwargs,
    to_annotate,
    anno_d,
    anno_alias,
    anno_style,
    anno_kwargs,
    anno_kwargs_single,
    arm_scale,
    anno_max_iter,
    arm_scale_d,
    arm_offset,
    anno_adjust,
    anno_xshift,
    anno_fixed_arm_length,
    maxy,
    anno_fontsize,
    font_family,
    region,
    region_anno_bbox_kwargs,
    skip,
    anno_max_rows=40,
    arrow_kwargs=None,
    anno_height=1,
    amode="int",
    snpid="SNPID",
    chrom="CHR",
    pos="POS",
    repel_force=0.02,
    verbose=True,
    _invert=False,
    log=Log()
):
    """
    Annotate significant variants in a single Manhattan/QQ plot with customizable styling and positioning options.
    
    Parameters
    ----------
    sumstats : pandas.DataFrame
        DataFrame containing summary statistics with columns for chromosome (CHR), position (POS), and p-values (P)
    anno : bool or str, optional
        Column name to use for annotation. If True, uses "CHR:POS". Default is True
    mode : str, optional
        Annotation mode ('r' for right style, 'tight' for compact style, 'expand' for spaced style). Default is 'r'
    ax1 : matplotlib.axes.Axes
        Matplotlib axis object to add annotations to
    highlight_i : list, optional
        List of indices to highlight with special styling. Default is empty list
    highlight_chrpos : dict, optional
        Dictionary mapping highlighted SNPs to custom CHR:POS values. Default is empty dict
    highlight_anno_kwargs : dict, optional
        Dictionary of styling arguments for highlighted variants. Default is empty dict
    to_annotate : pandas.DataFrame
        DataFrame containing variants to annotate with their positions and p-values
    anno_d : dict, optional
        Dictionary mapping annotation indices to custom positioning options. For example, {"0":"l"} means adjust the direction of the 0th to left. Default is empty dict
    anno_alias : dict, optional
        Dictionary mapping SNP IDs to custom annotation labels. Default is empty dict
    anno_style : str, optional
        Style of annotation ('right', 'tight', 'expand'). Default is 'right'
    anno_kwargs : dict, optional
        Dictionary of default styling arguments for annotations. Default is empty dict
    anno_kwargs_single : dict, optional
        Dictionary mapping SNP IDs to custom styling arguments. Default is empty dict
    arm_scale : float, optional
        Scaling factor for annotation arm length. Default is 1.0
    anno_max_iter : int, optional
        Maximum iterations for text repulsion algorithm. Default is 50
    anno_max_rows : int, optional
        Maximum number of annotation rows to display. If more variants are provided, they will be sorted by p-value or -log10(p-value) and only the top ones will be shown. Default is 40
    arm_scale_d : dict, optional
        Dictionary mapping annotation indices to custom arm scaling factors. Default is None
    arm_offset : float, optional
        Offset for custom arm positioning. Default is 0.0
    anno_adjust : bool, optional
        Whether to automatically adjust text positions to prevent overlap. Default is True
    anno_xshift : float, optional
        Horizontal shift for all annotations. Default is None
    anno_fixed_arm_length : float, optional
        Fixed arm length for annotations. Default is None
    maxy : float
        Maximum y-value for the plot
    anno_fontsize : int, optional
        Font size for annotations. Default is 10
    font_family : str, optional
        Font family for annotations. Default is 'sans-serif'
    region : tuple, optional
        Genomic region coordinates (start, end). Default is None
    region_anno_bbox_kwargs : dict, optional
        Dictionary of bounding box arguments for region annotations. Default is empty dict
    skip : float, optional
        Spacing between annotations. Default is 0.1
    arrow_kwargs : dict, optional
        Dictionary of arrow styling arguments. Default is None
    anno_height : float, optional
        Vertical spacing factor for annotations. Default is 1.0
    amode : str, optional
        Algorithm mode for text positioning ('int' or 'float'). Default is 'int'
    snpid : str, optional
        Column name containing SNP identifiers. Default is 'SNPID'
    chrom : str, optional
        Column name containing chromosome numbers. Default is 'CHR'
    pos : str, optional
        Column name containing genomic positions. Default is 'POS'
    repel_force : float, optional
        Force strength for text repulsion algorithm. Default is 0.02
    verbose : bool, optional
        Whether to print progress messages. Default is True
    _invert : bool, optional
        Whether to invert the plot orientation. Default is False
    log : gwaslab.Log, optional
        Log object for tracking progress. Default is new Log instance

    Returns
    -------
    matplotlib.axes.Axes
        Modified axis object with annotations added
    """
    if anno_d is not None:
        anno_d = {int(key): value for key, value in anno_d.items()}

    if arrow_kwargs is None:
        arrow_kwargs=dict()
    
    if (anno_d is not None) and (len(anno_d) > 0) and (arm_offset is None):
        try:
            dpi = ax1.figure.dpi
            width_in = ax1.figure.get_size_inches()[0]
            arm_offset = dpi * repel_force * width_in * 0.5
        except Exception:
            arm_offset = 0.0
        
    if anno and (to_annotate.empty is not True):
        # Limit number of annotations if exceeding anno_max_rows
        if len(to_annotate) > anno_max_rows:
            log.write(" -Found {} variants to annotate, limiting to top {} by significance...".format(len(to_annotate), anno_max_rows), verbose=verbose)
            # Try to sort by p-value or -log10(p-value)
            sort_column = None
            ascending = True
            
            if "scaled_P" in to_annotate.columns:
                sort_column = "scaled_P"
                ascending = False  # Higher -log10(p) is more significant
            elif "MLOG10P" in to_annotate.columns:
                sort_column = "MLOG10P"
                ascending = False  # Higher -log10(p) is more significant
            elif "raw_P" in to_annotate.columns:
                sort_column = "raw_P"
                ascending = True  # Lower p-value is more significant
            elif "P" in to_annotate.columns:
                sort_column = "P"
                ascending = True  # Lower p-value is more significant
            
            if sort_column is not None:
                to_annotate = to_annotate.sort_values(by=sort_column, ascending=ascending).head(anno_max_rows)
                log.write(" -Sorted by {} and selected top {} variants...".format(sort_column, anno_max_rows), verbose=verbose)
            else:
                # If no p-value column found, just take first anno_max_rows
                to_annotate = to_annotate.head(anno_max_rows)
                log.write(" -No p-value column found, selecting first {} variants...".format(anno_max_rows), verbose=verbose)
        
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
        log.write(" -Annotating using column "+annotation_col+"...", verbose=verbose)
        ################################################################################################################################
        ## calculate y span
        if region is not None:
            y_span = region[2] - region[1]
        else:
            y_span = sumstats["i"].max()-sumstats["i"].min()
        log.write(" -Adjusting text positions with repel_force={}...".format(repel_force), verbose=verbose)
        if anno_style == "expand" :
            to_annotate.loc[:, "ADJUSTED_i"] = adjust_text_position(to_annotate["i"].values.copy(), y_span, repel_force,max_iter=anno_max_iter,log=log,amode=amode,verbose=verbose)
        ##  iterate through variants to be annotated
        ################################################################################################################################

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
            elif anno_style == "tight" :
                #tight style
                anno_fixed_arm_length = 1
                anno_adjust = True
            else:
                pass
            ################################################################
            # shrink or increase the arm by a factor (arm_scale)
            if arm_scale_d is not None:
                if anno_count not in arm_scale_d.keys():
                    arm_scale =1
                else:
                    arm_scale = arm_scale_d[anno_count]
            ################################################################

            # vertical arm length in pixels
            # Annotation y : 1.15 * maxy_anno
            # Top dot:  1 * maxy_anno
            # armB_length_in_point_raw = 0.15 * maxy_anno -> gap_pixel
            # Fixed Offset: 0.5 * 0.15 * gap_pixel

            #Calculate armB length in pixels
            # arm_scale: raise up the ceiling
            
            # gap : 0.5* space between top variant and annotation text
            gap_pixel =                (ax1.transData.transform((0,1.15*maxy*arm_scale))[1]-ax1.transData.transform((0, maxy*arm_scale))[1])*0.5

            # armB_length_in_pixel_raw : distance between variant to annotate and annotation text
            armB_length_in_pixel_raw =  ax1.transData.transform((0,1.15*maxy*arm_scale))[1]-ax1.transData.transform((0, row["scaled_P"]+1))[1] 
            
            armB_length_in_pixel = armB_length_in_pixel_raw - gap_pixel

            ################################################################
            # armB_length_in_pixel should not be negative
            if arm_scale>=1:
                armB_length_in_pixel = max(0, armB_length_in_pixel)
            
            ################################################################
            #if setting anno_fixed_arm_length 
            if anno_fixed_arm_length is not None:
                armB_length_in_pixel = ax1.transData.transform((skip,anno_fixed_arm_length))[1]-ax1.transData.transform((skip,0))[1] 
            
            ################################################################################################################################
            # annotation alias
            if anno==True:
                if row[snpid] in anno_alias.keys():
                    annotation_text = anno_alias[row[snpid]]
                else:
                    annotation_text="Chr"+ str(row[chrom]) +":"+ str(int(row[pos]))
            elif anno:
                if row[snpid] in anno_alias.keys():
                    annotation_text = anno_alias[row[snpid]]
                else:
                    annotation_text=row["Annotation"]
            
            ################################################################################################################################
            # setting arrow xy and text xy
            # add a small space between variant and arrow head
            xy=(row["i"],row["scaled_P"]+0.01*maxy)

            # text xy is of the same height
            # anno_height can be used to adjust the height of annotation text
            xytext=(last_pos,1.15*maxy*(arm_scale + anno_height -1))
            
            # for anno_fixed_arm_length
            if anno_fixed_arm_length is not None:
                xytext=(row["i"],row["scaled_P"] + 0.2 + anno_fixed_arm_length)

            if anno_xshift is not None:
                xytext = (xytext[0] +(anno_xshift*y_span), xytext[1])
            ################################################################################################################################
            # if not changing the directions of some annotation arror arms 
            if anno_count not in anno_d.keys():
                if _invert==False:
                    arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                            connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_pixel)+",rad=0")
                else:
                    arrowargs = dict(arrowstyle="-|>",relpos=(0,1),color="#ebebeb",
                                            connectionstyle="arc,angleA=0,armA=0,angleB=-90,armB="+str(armB_length_in_pixel)+",rad=0")
            else:
                # if not changing the directions of some annotation arror arms 
                # adjust horizontal direction
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
            ################################################################################################################################
            
            if "r" in mode:
                arrowargs["color"] = "black" 
                bbox_para=dict(boxstyle="round", fc="white",zorder=3)
                for key,value in region_anno_bbox_kwargs.items():
                    bbox_para[key]=value
            else:
                bbox_para=None

            ################################################################################################################################
            if  _invert==False:
                anno_default = {"rotation":40, "fontstyle":"italic","ha":"left","va":"bottom","fontsize":anno_fontsize,"fontweight":"normal","fontfamily":font_family}
            else:
                anno_default = {"rotation":-40,"fontstyle":"italic","ha":"left","va":"top",   "fontsize":anno_fontsize,"fontweight":"normal","fontfamily":font_family}

            ################################################################################################################################
            if anno_style == "expand" :
                anno_default["rotation"] = 90
            if anno_style == "tight" :
                anno_default["rotation"] = 90
            ################################################################################################################################
            # anno args for all
            for key,value in anno_kwargs.items():
                anno_default[key]=value
            
            # anno args for highlight group
            if len(highlight_i) >0:
                if row["i"] in highlight_i:
                    if highlight_anno_kwargs is None:
                        highlight_anno_kwargs = dict()
                    for key,value in highlight_anno_kwargs.items():
                        anno_default[key]=value
            
            # anno args for specifc
            #try:
            if row[snpid] in anno_kwargs_single.keys():
                for key,value in anno_kwargs_single[row[snpid]].items():
                    anno_default[key]=value
            #except:
            #    pass
            ################################################################################################################################
            if anno_adjust==True:
                if  _invert==False:
                    arrowargs=dict(arrowstyle='-|>', color='grey', shrinkA=10, linewidth=0.1, relpos=(0,0.5))
                else:
                    arrowargs=dict(arrowstyle='-|>', color='grey', shrinkA=10, linewidth=0.1, relpos=(1,0.5))
            ################################################################################################################################
            for key,value in arrow_kwargs.items():
                arrowargs[key]=value

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
            ################################################################################################################################
        
        #anno_adjust_keyargs = {"arrowprops":dict(arrowstyle='->', color='grey', linewidth=0.1,relpos=(0.5,0.5))}
        if anno_adjust==True:
            log.write(" -Auto-adjusting text positions...", verbose=verbose)
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
        log.write(" -Skip annotating", verbose=verbose)
    
    return ax1



# miami plot
def annotate_pair(
    sumstats,
    anno,
    ax1,
    ax5,
    highlight_i,
    to_annotate1,
    to_annotate5,
    anno_d1,
    anno_d2,
    anno_alias1,
    anno_alias2,
    anno_style,
    anno_kwargs,
    arm_scale,
    anno_max_iter,
    arm_scale_d,
    arm_offset,
    anno_adjust,
    anno_fixed_arm_length,
    maxy1,
    maxy5,
    fontsize,
    font_family,
    region,
    skip,
    chrom="CHR",
    pos="POS",
    repel_force=0.02,
    verbose=True,
    arrow_kwargs=None,
    log=Log()
):
    if anno is not None:
        for index,ax,to_annotate_df,anno_d, anno_alias in [(0,ax1,to_annotate1,anno_d1,anno_alias1),(1,ax5,to_annotate5,anno_d2,anno_alias2)]:
            ###################### annotate() args
            if to_annotate_df.empty is True:
                log.write(" -Skipping annotation...", verbose=verbose)
                continue

            fontweight = "normal"

            anno_default =dict({"rotation":40,"fontstyle":"italic","ha":"left","va":"bottom","fontsize":fontsize,"fontweight":fontweight,"fontfamily":font_family})
            ########################
            to_annotate = to_annotate_df.copy()

            if index == 0:
                to_annotate["scaled_P"] = to_annotate1["scaled_P_1"].copy()
                maxy_anno = maxy1
            else:
                to_annotate["scaled_P"] = to_annotate5["scaled_P_2"].copy()
                anno_default["rotation"] = -40
                anno_default["ha"]="left"
                anno_default["va"]="top"
                maxy_anno = maxy5
            
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
                    if type(anno) == list:
                        annotation_col=anno[index]
                    else:
                        if anno =="GENENAME":
                            annotation_col=anno
                        else:
                            annotation_col=anno+"_"+str(index+1)
                log.write(" -Annotating using column "+annotation_col+"...", verbose=verbose)
                
                ## calculate y span
                if region is not None:
                    y_span = region[2] - region[1]
                else:
                    y_span = sumstats["i"].max()-sumstats["i"].min()
                ##   
                if anno_style == "expand" :
                    to_annotate.loc[:, "ADJUSTED_i"] = adjust_text_position(to_annotate["i"].values.copy(), y_span, repel_force, max_iter=anno_max_iter,log=log,verbose=verbose)
                
                anno_to_adjust_list = list()
                for rowi,row in to_annotate.iterrows():
                    
                    # avoid text overlapping
                    if anno_style == "right" :
                    #right style
                        if row["i"]>last_pos+repel_force*y_span:
                            last_pos=row["i"]
                        else:
                            last_pos+=repel_force*y_span
                    elif anno_style == "expand" :
                        #expand style
                        last_pos = row["ADJUSTED_i"]
                        anno_kwargs["rotation"] = 90
                    elif anno_style == "tight" :
                        #tight style
                        anno_fixed_arm_length = 1
                        anno_adjust = True
                        anno_kwargs["rotation"] = 90
                    else:
                        pass

                    if arm_scale_d is not None:
                        if anno_count not in arm_scale_d.keys():
                            arm_scale =1
                        else:
                            arm_scale = arm_scale_d[anno_count]

                    # vertical arm length in pixels
                    # Annotation y : 1.15 * maxy_anno
                    # Top dot:  1 * maxy_anno
                    # armB_length_in_point_raw = 0.15 * maxy_anno -> gap_pixel
                    # Fixed Offset: 0.5 * 0.15 * gap_pixel

                    #Calculate armB length in pixels
                    # arm_scale: raise up the ceiling
                    gap_pixel =                (ax1.transData.transform((0,1.15*maxy_anno*arm_scale))[1]-ax1.transData.transform((0, maxy_anno*arm_scale))[1])*0.5

                    armB_length_in_pixel_raw =  ax1.transData.transform((0,1.15*maxy_anno*arm_scale))[1]-ax1.transData.transform((0, row["scaled_P"]+1))[1] 
                    
                    armB_length_in_pixel = armB_length_in_pixel_raw - gap_pixel
                    
                    if arm_scale>=1:
                        armB_length_in_pixel= armB_length_in_pixel if armB_length_in_pixel>0 else 0
                    
                    if anno_fixed_arm_length is not None:
                        anno_fixed_arm_length_factor = ax.transData.transform((skip,anno_fixed_arm_length))[1]-ax.transData.transform((skip,0))[1] 
                        armB_length_in_pixel = anno_fixed_arm_length_factor
                    
                    if anno==True:
                        if row[snpid] in anno_alias.keys():
                            annotation_text = anno_alias[row[snpid]]
                        else:
                            annotation_text="Chr"+ str(row[chrom]) +":"+ str(int(row[pos]))
                    elif anno == "GENENAME":
                        annotation_text=row["GENENAME"]
                    elif anno is not None:
                        if type(anno) == list:
                            annotation_text=row[anno[index]]
                        else: 
                            annotation_text=row[anno+"_"+str(index+1)]

                    # annoatte position
                    xy=(row["i"],row["scaled_P"]+0.2)  
                    
                    xytext=(last_pos,1.15*maxy_anno*arm_scale)
                    
                    if anno_fixed_arm_length is not None:
                        armB_length_in_pixel = anno_fixed_arm_length
                        xytext=(row["i"],row["scaled_P"]+0.2+anno_fixed_arm_length)
                    
                    if anno_count not in anno_d.keys():
                        if index==0:
                            #upper panel
                            if armB_length_in_pixel <5:
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",connectionstyle="arc,armA=0,armB=0,rad=0.")
                            else:  
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_pixel)+",rad=0")
                        else:
                            #lower panel
                            if armB_length_in_pixel <5:
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",connectionstyle="arc,armA=0,armB=0,rad=0.")
                            else:
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,1),color="#ebebeb",
                                                    connectionstyle="arc,angleA=0,armA=0,angleB=-90,armB="+str(armB_length_in_pixel)+",rad=0")
                    
                    else:
                        xy=(row["i"],row["scaled_P"])
                        if anno_d[anno_count] in ["right","left","l","r"]:
                            if anno_d[anno_count]=="right" or anno_d[anno_count]=="r": 
                                armoffsetall = (ax.transData.transform(xytext)[0]-ax.transData.transform(xy)[0])*np.sqrt(2)
                                armoffsetb = arm_offset 
                                armoffseta = armoffsetall - armoffsetb   
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                            elif anno_d[anno_count]=="left" or anno_d[anno_count]=="l":
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str(arm_offset)+",angleB=135,armB="+str(arm_offset)+",rad=0")
                        else:
                            if anno_d[anno_count][0]=="right" or anno_d[anno_count][0]=="r": 
                                armoffsetall = (ax.transData.transform(xytext)[0]-ax.transData.transform(xy)[0])*np.sqrt(2)
                                armoffsetb = anno_d[anno_count][1] 
                                armoffseta = armoffsetall - armoffsetb   
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str(armoffseta)+",angleB=45,armB="+str(armoffsetb)+",rad=0")
                            elif anno_d[anno_count]=="left" or anno_d[anno_count]=="l":
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                    connectionstyle="arc,angleA=-135,armA="+str( anno_d[anno_count][1])+",angleB=135,armB="+str( anno_d[anno_count][1])+",rad=0")

                    bbox_para=None            
                    
                    if len(highlight_i) >0:
                        if row["i"] in highlight_i:
                            anno_default["fontweight"] = "bold"
                        else:
                            anno_default["fontweight"] = "normal"
                    
                    for key,value in anno_kwargs.items():
                        anno_default[key]=value

                    if anno_adjust==True:
                        if index==0:
                            arrowargs=dict(arrowstyle='-|>', color='grey', shrinkA=10, linewidth=0.1, relpos=(0,0.5))
                        else:
                            arrowargs=dict(arrowstyle='-|>', color='grey', shrinkA=10, linewidth=0.1, relpos=(1,0.5))

                    
                    anno_to_adjust = ax.annotate(annotation_text,
                            xy=xy,
                            xytext=xytext,
                            bbox=bbox_para,
                            arrowprops=arrowargs,
                            zorder=100,
                            **anno_default
                            )
                    anno_to_adjust_list.append(anno_to_adjust) 
                    anno_count +=1
            
            if anno_adjust==True:
                log.write(" -Auto-adjusting text positions for plot {}...".format(index), verbose=verbose)
                if index==0:
                    va="bottom"
                    ha='left'
                else:
                    va="top"
                    ha='left'
                
                adjust_text(texts = anno_to_adjust_list,
                            autoalign=False, 
                            only_move={'points':'x', 'text':'x', 'objects':'x'},
                            ax=ax,
                            precision=0.02,
                            force_text=(repel_force,repel_force),
                            expand_text=(1,1),
                            expand_objects=(0,0),
                            expand_points=(0,0),
                            va=va,
                            ha=ha,
                            avoid_points=False,
                            lim =anno_max_iter
                            )
    else:
        log.write(" -Skip annotating", verbose=verbose)
    return ax1,ax5
