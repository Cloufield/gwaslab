import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from matplotlib import ticker
import matplotlib.pyplot as plt
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from math import ceil

def _quick_fix(sumstats, chr_dict=get_chr_to_number(), scaled=False, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",log=Log(), verbose=True):
    '''
    quick sanity check for input sumstats
    '''

    # quick fix chr
    sumstats[chrom] = _quick_fix_chr(sumstats[chrom], chr_dict=chr_dict, verbose=verbose, log=log)

    # quick fix pos
    sumstats[pos] = _quick_fix_pos(sumstats[pos], verbose=verbose, log=log)

    # quick fix p
    sumstats = _quick_fix_p_value(sumstats, p=p, mlog10p=mlog10p, scaled=scaled, verbose=verbose, log=log)

    # quick fix mlog10p
    sumstats = _quick_fix_mlog10p(sumstats, p=p, mlog10p=mlog10p, scaled=scaled, verbose=verbose, log=log)

    return sumstats


def _quick_fix_p_value(sumstats, p="P", mlog10p="MLOG10P", scaled=False,verbose=True, log=Log()):
    '''
    drop variants with bad P values
    '''
    if scaled==True:
        # if scaled, add scaled P and P col 
        log.write(" -P values are already scaled...", verbose=verbose)
        log.write(" -Sumstats -log10(P) values are being converted to P...", verbose=verbose)
        sumstats["scaled_P"] = sumstats[mlog10p].copy()
        sumstats[p]= np.power(10,-sumstats[mlog10p].astype("float64"))
        return sumstats
    # bad p : na and outside (0,1]
    bad_p_value = (sumstats[p].isna()) | (sumstats[p] > 1) | (sumstats[p] <= 0)

    log.write(" -Sanity check after conversion: " + str(sum(bad_p_value)) +
                  " variants with P value outside of (0,1] will be removed...", verbose=verbose)
    sumstats = sumstats.loc[~bad_p_value, :]
    return sumstats


def _quick_fix_mlog10p(insumstats,p="P", mlog10p="MLOG10P", scaled=False, log=Log(), verbose=True):
    '''
    drop variants with bad -log10(P) values
    '''
    sumstats = insumstats.copy()
    if scaled != True:
        log.write(" -Sumstats P values are being converted to -log10(P)...", verbose=verbose)
        sumstats["scaled_P"] = -np.log10(sumstats[p].astype("float64"))
        
    #with pd.option_context('mode.use_inf_as_na', True):
    #    is_na = sumstats["scaled_P"].isna()
    if_inf_na = np.isinf(sumstats["scaled_P"]) | sumstats["scaled_P"].isna()

    log.write(" -Sanity check: "+str(sum(if_inf_na)) +
                  " na/inf/-inf variants will be removed...", verbose=verbose)
    sumstats = sumstats.loc[~if_inf_na, :]
    return sumstats


def _quick_fix_eaf(seires,log=Log(), verbose=True):
    '''
    conversion of eaf to maf
    '''
    seires = pd.to_numeric(seires, errors='coerce')
    flipped = seires > 0.5
    seires[flipped] = 1 - seires[flipped]
    return seires.copy()


def _quick_fix_chr(seires, chr_dict,log=Log(), verbose=True):
    '''
    conversion and check for chr
    '''
    if pd.api.types.is_string_dtype(seires) == True:
        # if chr is string dtype: convert using chr_dict
        seires = seires.astype("string")
        seires = seires.map(chr_dict, na_action="ignore")
    seires = np.floor(pd.to_numeric(seires, errors='coerce')).astype('Int64')
    return seires


def _quick_fix_pos(seires,log=Log(), verbose=True):
    '''
    force conversion for pos
    '''
    seires = np.floor(pd.to_numeric(seires, errors='coerce')).astype('Int64')
    return seires


def _dropna_in_cols(sumstats, cols, log=Log(), verbose=True):
    to_drop = sumstats[cols].isna().any(axis=1)
    log.write(" -Dropping {} variants due to missing values in {}.".format(sum(to_drop),cols))
    return sumstats.loc[~to_drop,:]


def _get_largenumber(*args,log=Log(), verbose=True):
    '''
    get a helper large number, >> max(pos)
    '''
    large_number = 1000000000
    max_number = np.nanmax(args)
    for i in range(7):
        if max_number*10 > large_number:
            large_number = int(large_number * 10)
        else:
            break
        if i == 7:
            log.warning("Max POS is too large!")
    return large_number


def _quick_add_tchrpos(sumstats, chr="chr", pos="POS", large_number=10000000000, dropchrpos=False,log=Log(), verbose=True):
    sumstats["TCHR+POS"] = sumstats["CHR"]*large_number + sumstats["POS"]
    sumstats["TCHR+POS"] = sumstats["TCHR+POS"].astype('Int64')
    if dropchrpos == True:
        sumstats = sumstats.drop(labels=["CHR", "POS"], axis=1)
    sumstats = sumstats.dropna()
    return sumstats


def _quick_merge_sumstats(sumstats1, sumstats2, log=Log(), verbose=True):
    merged_sumstats = pd.merge(sumstats1, sumstats2, on="TCHR+POS", how="outer", suffixes=('_1', '_2'))
    merged_sumstats["CHR"] = merged_sumstats["CHR_1"]
    merged_sumstats["POS"] = merged_sumstats["POS_1"]
    merged_sumstats.loc[merged_sumstats["CHR_1"].isna(), "CHR"] = merged_sumstats.loc[merged_sumstats["CHR_1"].isna(), "CHR_2"]
    merged_sumstats.loc[merged_sumstats["POS_1"].isna(), "POS"] = merged_sumstats.loc[merged_sumstats["POS_1"].isna(), "POS_2"]
    merged_sumstats = merged_sumstats.drop(labels=["CHR_1", "CHR_2", "POS_1", "POS_2"],axis=1)
    return merged_sumstats

def _quick_assign_i(sumstats, chrom="CHR",pos="POS",log=Log(), verbose=True):
    # sort by CHR an POS
    sumstats = sumstats.sort_values([chrom,pos])
    # set new id
    sumstats["_ID"]=range(len(sumstats))
    sumstats = sumstats.set_index("_ID")
    #create a df , groupby by chromosomes , and get the maximum position
    posdic = sumstats.groupby(chrom)[pos].max()
    # convert to dictionary
    posdiccul = dict(posdic)
    # fill empty chr with 0
    for i in range(0,26):
        if i in posdiccul: 
            continue
        else: 
            posdiccul[i]=0
    # get interval
    interval_between_chr = sumstats[pos].max()*0.05
    # cumulative sum dictionary
    for i in range(1,sumstats[chrom].max()+1):
        posdiccul[i] =  posdiccul[i-1] + posdiccul[i] + interval_between_chr
    # convert base pair postion to x axis position using the cumulative sum dictionary
    sumstats["_ADD"] = sumstats[chrom].apply(lambda x : posdiccul[int(x)-1])
    sumstats["i"]   = sumstats[pos] + sumstats["_ADD"]
    # drop add
    sumstats = sumstats.drop(labels=["_ADD"],axis=1)
    #for plot, get the chr text tick position  
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    # fix dtype for i
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
    return sumstats, chrom_df

def _quick_assign_i_with_rank(sumstats, chrpad, use_rank=False, chrom="CHR",pos="POS",drop_chr_start=False,_posdiccul=None,log=Log(), verbose=True):
    # align all variants on a single axis (i)
    sumstats = sumstats.sort_values([chrom,pos])
    if use_rank is True: 
        sumstats["_POS_RANK"] = sumstats.groupby(chrom)[pos].rank("dense", ascending=True)
        pos="_POS_RANK"
    sumstats["_ID"]=range(len(sumstats))
    sumstats=sumstats.set_index("_ID")

    #create a df , groupby by chromosomes , and get the maximum position
    if use_rank is True: 
        posdic = sumstats.groupby(chrom)["_POS_RANK"].max()
    else:
        posdic = sumstats.groupby(chrom)[pos].max()
    
    if _posdiccul is None:
        # convert to dictionary
        posdiccul = dict(posdic)
        
        # fill empty chr with 0
        for i in range(0,sumstats[chrom].max()+1):
            if i in posdiccul: 
                continue
            else: 
                posdiccul[i]=0

        # cumulative sum dictionary
        for i in range(1,sumstats[chrom].max()+1):
            posdiccul[i]= posdiccul[i-1] + posdiccul[i] + sumstats[pos].max()*chrpad
    else:
        posdiccul = _posdiccul

    # convert base pair postion to x axis position using the cumulative sum dictionary
    sumstats["_ADD"]=sumstats[chrom].apply(lambda x : posdiccul[int(x)-1])
    
    if drop_chr_start==True:
            posdic_min =  sumstats.groupby(chrom)[pos].min()
            posdiccul_min= dict(posdic_min)
            for i in range(0,sumstats[chrom].max()+1):
                if i in posdiccul_min: 
                    continue
                else: 
                    posdiccul_min[i]=0
            for i in range(1,sumstats[chrom].max()+1):
                posdiccul_min[i]= posdiccul_min[i-1] + posdiccul_min[i]
            sumstats["_ADD"]=sumstats["_ADD"] - sumstats[chrom].apply(lambda x : posdiccul_min[int(x)])
        
    if use_rank is True: 
        sumstats["i"]=sumstats["_POS_RANK"]+sumstats["_ADD"]
    else:
        sumstats["i"]=sumstats[pos]+sumstats["_ADD"]
    

    #for plot, get the chr text tick position      
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    #sumstats["i"] = sumstats["i"]+((sumstats[chrom].map(dict(chrom_df)).astype("int")))*0.02
    #sumstats["i"] = sumstats["i"].astype("Int64")
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
    return sumstats, chrom_df

def _quick_assign_marker_relative_size(series, sig_level = 5e-8, suggestive_sig_level=5e-6, lower_level=5e-4,log=Log(), verbose=True):
    size_series = series.copy()
    size_series[:] = 1

    is_lower_level          = series > -np.log10(lower_level)
    is_suggestive_sig_level = series > -np.log10(suggestive_sig_level)
    is_sig_level            = series > -np.log10(sig_level)
    
    size_series[is_lower_level] = 2
    size_series[is_suggestive_sig_level] = 3
    size_series[is_sig_level] = 4
    return size_series

def _quick_assign_highlight_hue(sumstats,highlight,highlight_windowkb, snpid="SNPID",chrom="CHR",pos="POS",log=Log(), verbose=True):
    to_highlight = sumstats.loc[sumstats[snpid].isin(highlight),:]
    #assign colors: 0 is hightlight color
    for i,row in to_highlight.iterrows():
        target_chr = int(row[chrom])
        target_pos = int(row[pos])
        right_chr=sumstats[chrom]==target_chr
        up_pos=sumstats[pos]>target_pos-highlight_windowkb*1000
        low_pos=sumstats[pos]<target_pos+highlight_windowkb*1000
        sumstats.loc[right_chr&up_pos&low_pos,"HUE"]="0"
    return sumstats

def _quick_assign_highlight_hue_pair(sumstats, highlight1, highlight2, highlight_windowkb, chrom="CHR",pos="POS",log=Log(), verbose=True):
    #assign colors: 0 is hightlight color
    to_highlight1 = pd.DataFrame()
    to_highlight2 = pd.DataFrame()
    
    if len(highlight1)>0 :
        to_highlight1 = sumstats.loc[sumstats["TCHR+POS"].isin(highlight1),:]
    if len(highlight2)>0:
        to_highlight2 = sumstats.loc[sumstats["TCHR+POS"].isin(highlight2),:]
    
    if len(to_highlight1)>0:
        for i,row in to_highlight1.iterrows():
            target_chr = int(row[chrom])
            target_pos = int(row[pos])
            right_chr = sumstats["CHR"]==target_chr  
            up_pos = sumstats["POS"]>target_pos-highlight_windowkb*1000
            low_pos = sumstats["POS"]<target_pos+highlight_windowkb*1000
            sumstats.loc[right_chr&up_pos&low_pos,"HUE1"]="0"
    if len(to_highlight2)>0:
        for i,row in to_highlight2.iterrows():
            target_chr = int(row[chrom])
            target_pos = int(row[pos])
            right_chr = sumstats["CHR"]==target_chr  
            up_pos = sumstats["POS"]>target_pos-highlight_windowkb*1000
            low_pos = sumstats["POS"]<target_pos+highlight_windowkb*1000
            sumstats.loc[right_chr&up_pos&low_pos,"HUE2"]="0"
    return sumstats, to_highlight1, to_highlight2

def _quick_extract_snp_in_region(sumstats, region, chrom="CHR",pos="POS",log=Log(), verbose=True):
    region_chr = region[0]
    region_start = region[1]
    region_end = region[2]
    log.write(" -Extract SNPs in region : chr"+str(region_chr)+":"+str(region[1])+"-"+str(region[2])+ "...", verbose=verbose)
    is_in_region_snp = (sumstats[chrom]==region_chr) &(sumstats[pos]<region_end) &(sumstats[pos]>region_start)
    log.write(" -Extract SNPs in specified regions: "+str(sum(is_in_region_snp)), verbose=verbose)
    sumstats = sumstats.loc[is_in_region_snp,:]
    return sumstats

def _cut(series, mode,cutfactor,cut,skip, ylabels, cut_log, verbose, lines_to_plot, log):
    log.write(" -Converting data above cut line...",verbose=verbose)
    if ylabels is not None:
        ylabels = pd.Series(ylabels)
    series = series.copy()
    
    maxy = series.max()
    if "b" not in mode:
        log.write(" -Maximum -log10(P) value is "+str(maxy) +" .", verbose=verbose)
    elif "b" in mode:
        log.write(" -Maximum DENSITY value is "+str(maxy) +" .", verbose=verbose)
    
    maxticker=int(np.round(series.max(skipna=True)))
    
    if cut:
        # auto mode : determine curline and cut factor
        if cut==True:
            log.write(" -Cut Auto mode is activated...", verbose=verbose)
            if maxy<30:
                log.write(" - maxy <30 , no need to cut.", verbose=verbose)
                cut=0
            else:
                cut = 20
                cutfactor = ( maxy - cut )/8
        
        if cut:
            #cut log mode
            if cut_log==True:
                maxticker=int(np.round(series.max(skipna=True)))
                
                amp = (cut - skip)/ 2 / np.log2(maxticker/cut)
                
                # scaled_P
                series[series>cut] = (np.log2(series[series>cut]/cut)) * amp + cut
                
                # y labels
                if ylabels is not None:
                    ylabels[ylabels>cut] = (np.log2(ylabels[ylabels>cut]/cut)) * amp +cut 
                
                # lines
                lines_to_plot[lines_to_plot>cut] = (np.log2(lines_to_plot[lines_to_plot>cut]/cut)) * amp +cut 
                
                maxy = (np.log2(maxticker) - np.log2(cut)) * amp + cut
            else:
                # cut linear mode
                if "b" not in mode:
                    log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...", verbose=verbose)
                else:
                    log.write(" -DENSITY values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...", verbose=verbose)

                maxticker=int(np.round(series.max(skipna=True)))

                series[series>cut] = (series[series>cut]-cut)/cutfactor+cut
                
                if ylabels is not None:
                    ylabels[ylabels>cut] = (ylabels[ylabels>cut]-cut)/cutfactor+cut
                
                lines_to_plot[lines_to_plot>cut] = (lines_to_plot[lines_to_plot>cut]-cut)/cutfactor+cut
                #sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
                
                maxy = (maxticker-cut)/cutfactor + cut

    return series, maxy, maxticker, cut, cutfactor,ylabels,lines_to_plot

#def _cut_line(level, mode,cutfactor,cut,skip, ylabels, cut_log, verbose, log):




def _set_yticklabels(cut,
                     cutfactor,
                     cut_log,
                     ax1,
                     skip,
                     maxy,
                     maxticker,
                     ystep,
                     sc_linewidth,
                     cut_line_color,
                     fontsize,
                     font_family,
                     ytick3,
                     ylabels,
                     ylabels_converted,
                     log=Log(),
                     verbose=True
                     ):
    log.write(" -Processing Y tick lables...",verbose=verbose)
    # if no cut
    if cut == 0: 
            ax1.set_ylim((skip, ceil(maxy*1.2)) )
    # if cut     
    if cut!=0:
        # add cut line
        
        cutline = ax1.axhline(y=cut, linewidth = sc_linewidth,linestyle="--",color=cut_line_color,zorder=1)
        
        # default step
        step=2

        if ystep == 0:
            if (cut - skip ) // step > 8:
                step = (cut - skip ) // 8
        else:
            step = ystep

        if (cut-skip)%step==0:
            upper = cut - 1
        else:
            upper = cut
        
        ticks1= [x for x in range(skip,upper,step)]
        ticks2= [cut]
        tickslabel1= [x for x in range(skip,upper,step)]
        tickslabel2= [cut]
        
        if cut_log==True:
            ticks3= [x for x in range(cut,int(maxy),max(int(np.log2(maxticker//(cut-skip))),1))] 
            ticks4= [int(maxy)]
        elif ytick3 == True:
            ticks3= [x for x in range(cut,int((maxticker-cut)/cutfactor + cut),int(step))]
            ticks4= [(maxticker-cut)/cutfactor + cut]
        else:
            ticks3= []
            ticks4= [(maxticker-cut)/cutfactor + cut]
        
        #tickslabel3= [x for x in range(cut,int(maxticker),int(step*cutfactor))]
        if cut_log==True:
            amp = (cut - skip)/ 2 / np.log2(maxticker/cut) 
            tickslabel3 = list(map(lambda x: int(2**((x-cut)/amp)* cut) ,ticks3))
        elif ytick3 == True:
            tickslabel3 = list(map(lambda x: int((x-cut)*cutfactor + cut) ,ticks3))
        else:
            tickslabel3=[]
        tickslabel4= [maxticker]

        if maxy > cut:
            ax1.set_yticks(ticks1+ticks2+ticks3+ticks4)
            ax1.set_yticklabels(tickslabel1+tickslabel2+tickslabel3+tickslabel4,fontsize=fontsize,family=font_family)
        else:
            ax1.set_yticks(ticks1+ticks2)
            ax1.set_yticklabels(tickslabel1+tickslabel2,fontsize=fontsize,family=font_family)
    
    if ylabels is not None:  
        ax1.set_yticks(ylabels_converted)
        ax1.set_yticklabels(ylabels,fontsize=fontsize,family=font_family)
    
    ylim_top = ax1.get_ylim()[1]
    ax1.set_ybound(lower=skip,upper=ylim_top)
    ax1.tick_params(axis='y', labelsize=fontsize)

    return ax1

def _jagged_y(cut,skip,ax1,mode,mqqratio,jagged_len,jagged_wid, log=Log(), verbose=True):
    log.write(" -Processing jagged Y axis...",verbose=verbose)
    tycut = cut +0.3 #(cut - skip)/ (ax1.get_ylim()[1] - skip) + 0.002
    dy= jagged_len * (cut - skip) 
    x0 =  0
    dx= jagged_wid
    if mode>1:
        dx = dx * mqqratio
    kwargs = dict(transform=ax1.get_yaxis_transform(), color='k', clip_on=False,solid_capstyle="round",linewidth=0.8)
    ax1.plot((x0,x0), (tycut,tycut+4*dy), zorder=1000, transform=ax1.get_yaxis_transform(), color="white",clip_on=False,solid_capstyle="butt")
    ax1.plot((x0,-dx), (tycut,tycut+dy), zorder=1001, **kwargs)
    ax1.plot((-dx,+dx), (tycut+dy,tycut+3*dy), zorder=1001, **kwargs)
    ax1.plot((+dx,x0), (tycut+3*dy,tycut+4*dy), zorder=1001,  **kwargs)
    return ax1
