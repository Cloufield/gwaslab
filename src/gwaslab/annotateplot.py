import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from adjustText import adjust_text
from gwaslab.Log import Log
from gwaslab.textreposition import adjust_text_position


def annotate_single(
    sumstats,
    anno,
    mode,
    ax1,
    highlight_i,
    to_annotate,
    anno_d,
    anno_alias,
    anno_style,
    anno_args,
    arm_scale,
    anno_max_iter,
    arm_scale_d,
    arm_offset,
    anno_adjust,
    anno_fixed_arm_length,
    maxy,
    anno_fontsize,
    region,
    region_anno_bbox_args,
    skip,
    snpid="SNPID",
    chrom="CHR",
    pos="POS",
    repel_force=0.02,
    verbose=True,
    log=Log()
):
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
                armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else ax1.transData.transform((skip, maxy+2))[1]-ax1.transData.transform((skip,  row["scaled_P"]+1))[1] 
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
            if len(highlight_i) >0:
                if row["i"] in highlight_i:
                    fontweight = "bold"

            
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
    
    return ax1




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
    anno_args,
    arm_scale,
    anno_max_iter,
    arm_scale_d,
    arm_offset,
    anno_adjust,
    anno_fixed_arm_length,
    maxy1,
    maxy5,
    fontsize,
    region,
    skip,
    chrom="CHR",
    pos="POS",
    repel_force=0.02,
    verbose=True,
    log=Log()
):

    if anno is not None:
        for index,ax,to_annotate_df,anno_d, anno_alias in [(0,ax1,to_annotate1,anno_d1,anno_alias1),(1,ax5,to_annotate5,anno_d2,anno_alias2)]:
            ###################### annotate() args
            fontweight = "normal"

            anno_default =dict({"rotation":40,"style":"italic","ha":"left","va":"bottom","fontsize":fontsize,"fontweight":fontweight})
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
                        annotation_col = anno
                if verbose: log.write(" -Annotating using column "+annotation_col+"...")
                
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
                        anno_args["rotation"] = 90
                    elif anno_style == "tight" :
                        #tight style
                        anno_fixed_arm_length = 1
                        anno_adjust = True
                        anno_args["rotation"] = 90
                    else:
                        pass

                    if arm_scale_d is not None:
                        if anno_count not in arm_scale_d.keys():
                            arm_scale =1
                        else:
                            arm_scale = arm_scale_d[anno_count]

                    # vertical arm length in pixels
                    armB_length_in_point = ax.transData.transform((skip,1.15*maxy_anno))[1]-ax.transData.transform((skip, row["scaled_P"]+1))[1]
                    # times arm_scale to increase or reduce the length
                    armB_length_in_point = armB_length_in_point*arm_scale
                    
                    if arm_scale>=1:
                        armB_length_in_point= armB_length_in_point if armB_length_in_point>0 else 0 #ax.transData.transform((skip, maxy_anno+2))[1]-ax.transData.transform((skip,  row["scaled_P"]+1))[1] 
                    if anno_fixed_arm_length is not None:
                        anno_fixed_arm_length_factor = ax.transData.transform((skip,anno_fixed_arm_length))[1]-ax.transData.transform((skip,0))[1] 
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
                    
                    xytext=(last_pos,1.15*maxy_anno*arm_scale)
                    
                    if anno_fixed_arm_length is not None:
                        armB_length_in_point = anno_fixed_arm_length
                        xytext=(row["i"],row["scaled_P"]+0.2+anno_fixed_arm_length)
                    
                    if anno_count not in anno_d.keys():
                        if index==0:
                            #upper panel
                            if armB_length_in_point <5:
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",connectionstyle="arc,armA=0,armB=0,rad=0.")
                            else:  
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",
                                                connectionstyle="arc,angleA=0,armA=0,angleB=90,armB="+str(armB_length_in_point)+",rad=0")
                        else:
                            #lower panel
                            if armB_length_in_point <5:
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,0),color="#ebebeb",connectionstyle="arc,armA=0,armB=0,rad=0.")
                            else:
                                arrowargs = dict(arrowstyle="-|>",relpos=(0,1),color="#ebebeb",
                                                    connectionstyle="arc,angleA=0,armA=0,angleB=-90,armB="+str(armB_length_in_point)+",rad=0")
                    
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
                    
                    for key,value in anno_args.items():
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
                if verbose: log.write(" -Auto-adjusting text positions for plot {}...".format(index))
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
        if verbose: log.write(" -Skip annotating")
    return ax1,ax5





