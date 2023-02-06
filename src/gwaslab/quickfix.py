import pandas as pd
import numpy as np
from gwaslab.Log import Log
from gwaslab.CommonData import get_chr_to_number
from gwaslab.CommonData import get_number_to_chr


def _quick_fix(sumstats, chr_dict=get_chr_to_number(), scaled=False, chrom="CHR", pos="POS", p="P", verbose=True, log=Log()):
    '''
    quick sanity check for input sumstats
    '''

    # quick fix chr
    sumstats[chrom] = _quick_fix_chr(sumstats[chrom], chr_dict=chr_dict, verbose=verbose, log=log)

    # quick fix pos
    sumstats[pos] = _quick_fix_pos(sumstats[pos], verbose=verbose, log=log)

    # quick fix p
    sumstats = _quick_fix_p_value(sumstats, scaled=scaled, verbose=verbose, log=log)

    # quick fix mlog10p
    sumstats = _quick_fix_mlog10p(sumstats, scaled=scaled, verbose=verbose, log=log)

    return sumstats


def _quick_fix_p_value(sumstats, scaled=False,verbose=True, log=Log()):
    '''
    drop variants with bad P values
    '''
    if scaled==True:
        # if scaled, add scaled P and P col 
        if verbose:log.write(" -P values are already scaled...")
        if verbose:log.write(" -Sumstats -log10(P) values are being converted to P...")
        sumstats["scaled_P"] = sumstats["P"].copy()
        sumstats["P"]= np.power(10,-sumstats["P"])
    
    # bad p : na and outside (0,1]
    bad_p_value = (sumstats["P"].isna()) | (
        sumstats["P"] > 1) | (sumstats["P"] <= 0)
    if verbose:
        log.write(" -Sanity check after conversion: " + str(sum(bad_p_value)) +
                  " variants with P value outside of (0,1] will be removed...")
    sumstats = sumstats.loc[~bad_p_value, :]
    return sumstats


def _quick_fix_mlog10p(sumstats, scaled=False, verbose=True, log=Log()):
    '''
    drop variants with bad -log10(P) values
    '''
    if scaled != True:
        if verbose:log.write(" -sumstats P values are being converted to -log10(P)...")
        sumstats["scaled_P"] = -np.log10(sumstats["P"])
        
    with pd.option_context('mode.use_inf_as_na', True):
        is_na = sumstats["scaled_P"].isna()
    if verbose:
        log.write(" -Sanity check: "+str(sum(is_na)) +
                  " na/inf/-inf variants will be removed...")
    sumstats = sumstats.loc[~is_na, :]
    return sumstats


def _quick_fix_eaf(seires, verbose=True, log=Log()):
    '''
    conversion of eaf to maf
    '''
    seires = pd.to_numeric(seires, errors='coerce')
    flipped = seires > 0.5
    seires[flipped] = 1 - seires[flipped]
    return seires.copy()


def _quick_fix_chr(seires, chr_dict, verbose=True, log=Log()):
    '''
    conversion and check for chr
    '''
    if pd.api.types.is_string_dtype(seires) == True:
        # if chr is string dtype: convert using chr_dict
        seires = seires.map(chr_dict, na_action="ignore")
    seires = np.floor(pd.to_numeric(seires, errors='coerce')).astype('Int64')
    return seires


def _quick_fix_pos(seires, verbose=True, log=Log()):
    '''
    force conversion for pos
    '''
    seires = np.floor(pd.to_numeric(seires, errors='coerce')).astype('Int64')
    return seires


def _get_largenumber(*args, log=Log()):
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
            log.write(" -Warning: max POS is too large!")
    return large_number


def _quick_add_tchrpos(sumstats, chr="chr", pos="POS", large_number=10000000000, dropchrpos=False, verbose=True, log=Log()):
    sumstats["TCHR+POS"] = sumstats["CHR"]*large_number + sumstats["POS"]
    sumstats["TCHR+POS"] = sumstats["TCHR+POS"].astype('Int64')
    if dropchrpos == True:
        sumstats = sumstats.drop(labels=["CHR", "POS"], axis=1)
    sumstats = sumstats.dropna()
    return sumstats


def _quick_merge_sumstats(sumstats1, sumstats2):
    merged_sumstats = pd.merge(sumstats1, sumstats2, on="TCHR+POS", how="outer", suffixes=('_1', '_2'))
    merged_sumstats["CHR"] = merged_sumstats["CHR_1"]
    merged_sumstats["POS"] = merged_sumstats["POS_1"]
    merged_sumstats.loc[merged_sumstats["CHR_1"].isna(), "CHR"] = merged_sumstats.loc[merged_sumstats["CHR_1"].isna(), "CHR_2"]
    merged_sumstats.loc[merged_sumstats["POS_1"].isna(), "POS"] = merged_sumstats.loc[merged_sumstats["POS_1"].isna(), "POS_2"]
    merged_sumstats = merged_sumstats.drop(labels=["CHR_1", "CHR_2", "POS_1", "POS_2"],axis=1)
    return merged_sumstats

def _quick_assign_i(sumstats, chrom="CHR",pos="POS"):
    # sort by CHR an POS
    sumstats = sumstats.sort_values([chrom,pos])
    # set new id
    sumstats["_id"]=range(len(sumstats))
    sumstats = sumstats.set_index("_id")
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
    for i in range(2,sumstats[chrom].max()+1):
        posdiccul[i] =  posdiccul[i-1] + posdiccul[i] + interval_between_chr
    # convert base pair postion to x axis position using the cumulative sum dictionary
    sumstats["add"] = sumstats[chrom].apply(lambda x : posdiccul[int(x)-1])
    sumstats["i"]   = sumstats[pos] + sumstats["add"]
    # drop add
    sumstats = sumstats.drop(labels=["add"],axis=1)
    #for plot, get the chr text tick position  
    chrom_df=sumstats.groupby(chrom)['i'].agg(lambda x: (x.min()+x.max())/2)
    # fix dtype for i
    sumstats["i"] = np.floor(pd.to_numeric(sumstats["i"], errors='coerce')).astype('Int64')
    return sumstats, chrom_df

def _quick_assign_i_with_rank(sumstats, use_rank=False, chrom="CHR",pos="POS"):
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
        return sumstats, chrom_df

def _quick_assign_marker_relative_size(series, sig_level = 5e-8, suggestive_sig_level=5e-6, lower_level=5e-4):
    size_series = series.copy()
    size_series[:] = 1

    is_lower_level          = series > -np.log10(lower_level)
    is_suggestive_sig_level = series > -np.log10(suggestive_sig_level)
    is_sig_level            = series > -np.log10(sig_level)
    
    size_series[is_lower_level] = 2
    size_series[is_suggestive_sig_level] = 3
    size_series[is_sig_level] = 4
    return size_series

def _quick_assign_highlight_hue(sumstats,highlight,highlight_windowkb, snpid="SNPID",chrom="CHR",pos="POS",verbose=True, log=Log()):
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

def _quick_assign_highlight_hue_pair(sumstats, highlight1, highlight2, highlight_windowkb,chrom="CHR",pos="POS",verbose=True, log=Log()):
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

def _quick_extract_snp_in_region(sumstats, region, chrom="CHR",pos="POS",verbose=True, log=Log()):
    region_chr = region[0]
    region_start = region[1]
    region_end = region[2]
    if verbose:log.write(" -Extract SNPs in region : chr"+str(region_chr)+":"+str(region[1])+"-"+str(region[2])+ "...")
    is_in_region_snp = (sumstats[chrom]==region_chr) &(sumstats[pos]<region_end) &(sumstats[pos]>region_start)
    if verbose:log.write(" -Extract SNPs in specified regions: "+str(sum(is_in_region_snp)))
    sumstats = sumstats.loc[is_in_region_snp,:]
    return sumstats

def _cut(series, mode,cutfactor,cut, verbose, log):
    maxy = series.max()
    if "b" not in mode:
        if verbose: log.write(" -Maximum -log10(P) values is "+str(maxy) +" .")
    elif "b" in mode:
        if verbose: log.write(" -Maximum DENSITY values is "+str(maxy) +" .")
    maxticker=int(np.round(series.max(skipna=True)))
    if cut:
        if cut is True:
            if verbose: log.write(" -Cut Auto mode is activated...")
            if maxy<30:
                if verbose: log.write(" - maxy <20 , no need to cut.")
                cut=0
            else:
                cut = 20
                cutfactor = ( maxy - cut )/8
        if cut:
            if "b" not in mode:
                if verbose: log.write(" -Minus log10(P) values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")
            else:
                if verbose: log.write(" -Minus DENSITY values above " + str(cut)+" will be shrunk with a shrinkage factor of " + str(cutfactor)+"...")

            maxticker=int(np.round(series.max(skipna=True)))
            series[series>cut] = (series[series>cut]-cut)/cutfactor +  cut
            #sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"] = (sumstats.loc[sumstats["scaled_P"]>cut,"scaled_P"]-cut)/cutfactor +  cut
            maxy = (maxticker-cut)/cutfactor + cut
    if verbose: log.write("Finished data conversion and sanity check.")
    return series, maxy, maxticker, cut, cutfactor

def _set_yticklabels():
    pass