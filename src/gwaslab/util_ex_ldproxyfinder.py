from shutil import which
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
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
from gwaslab.util_in_filter_value import _get_flanking
from gwaslab.hm_harmonize_sumstats import auto_check_vcf_chr_dict
# unmatched SNP list 1 

# for each SNP in unmatched SNP list 1:
    # calculate LD for the variants
    
    # if the variant not in reference : return
    
    # if no variants reached threshold : return

    # if there are vailable possible proxies : 
        # from high to low: 
            #check if in outcome and exposure snp list
            #replace

def _extract_with_ld_proxy( snplist=None,
                            common_sumstats=None,
                            sumstats1=None,
                            vcf_path=None, 
                            vcf_chr_dict=None,
                            tabix=None,
                            log=Log(), 
                            verbose=True, 
                            windowsizekb=100,
                            ld_threshold=0.8
                            ):
    ### Load vcf#######################################################################################
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)
    if tabix is None:
        tabix = which("tabix")
    vcf_chr_dict = auto_check_vcf_chr_dict(vcf_path=vcf_path, vcf_chr_dict=vcf_chr_dict, verbose=verbose, log=log)

    is_needed=[]
    no_need  =[]
    
    print(common_sumstats.head())
    for i in snplist:
        if i in common_sumstats["SNPID"].values:
            no_need.append(i)
        else:
            is_needed.append(i)

    extracted_sumstats  = common_sumstats.loc[common_sumstats["SNPID"].isin(no_need),:]
    
    ld_proxies  = pd.DataFrame()
    in_sumstats = sumstats1.loc[sumstats1["SNPID"].isin(is_needed),:]
    
    if len(in_sumstats)==0:
        log.write(" -No available variants for LD proxy checking...Skipping... ", verbose=verbose)
    else:
        log.write(" -{} available variants for LD proxy checking... ".format(len(in_sumstats)), verbose=verbose)

    for index,row in in_sumstats.iterrows():
        # determine SNP and select region
        snpid = row["SNPID"]
        chrom= int(row["CHR"])
        start= int(row["POS"]-windowsizekb*1000)
        end=   int(row["POS"]+windowsizekb*1000)

        region = (chrom, start, end)
        
        ###  #######################################################################################
        #is_flanking = common_sumstats["CHR"] == chrom & common_sumstats["CHR"]>start & common_sumstats["CHR"]<end
        #flanking_sumstats = common_sumstats.loc[is_flanking,:]
        flanking_sumstats = common_sumstats.query('CHR == @chrom and @start < POS < @end',engine='python').copy()
        
        log.write(" -Extract {} variants in flanking region of {} for checking: {}:{}-{}".format(len(flanking_sumstats), snpid, chrom, start, end), verbose=verbose)

        if len(flanking_sumstats)==0:
            log.write("  -No availble variants in the region...Skipping!", verbose=verbose)
            continue
        
        _get_rsq_single(in_sumstats.loc[index,["POS","NEA_1","EA_1"]], 
                        row_pos=row["POS"], 
                        vcf_path=vcf_path, 
                        region=region,
                        log=log, 
                        verbose=verbose, 
                        vcf_chr_dict=vcf_chr_dict, 
                        tabix=tabix)


        flanking_sumstats = _get_rsq(row =in_sumstats.loc[index,["POS","NEA_1","EA_1"]],
                                     sumstats = flanking_sumstats, 
                                     row_pos=row["POS"], 
                                     vcf_path=vcf_path, 
                                     region=region,
                                     log=log, 
                                     verbose=verbose, 
                                     vcf_chr_dict=vcf_chr_dict, 
                                     tabix=tabix)
        if flanking_sumstats is None:
            log.write("  -{} is not found in the vcf...Skipping!".format(snpid))
            continue
        flanking_sumstats = flanking_sumstats.loc[flanking_sumstats["RSQ"]>ld_threshold,:]
        
        log.write("  -Variants in LD with {} (RSQ > {}): {}".format(snpid, ld_threshold,len(flanking_sumstats)), verbose=verbose)
        
        if len(flanking_sumstats)>0:
            flanking_sumstats["LD_REF_VARIANT"]= snpid
            for i,row_with_rsq in flanking_sumstats.iterrows():
                if row_with_rsq["SNPID"] in common_sumstats["SNPID"].values:
                    log.write("  -Proxy for {} is found: {} (LD RSQ= {})".format(snpid, row_with_rsq["SNPID"], row_with_rsq["RSQ"]))
                    row_with_rsq = pd.DataFrame(row_with_rsq)
                    ld_proxies = pd.concat([ld_proxies, row_with_rsq.T], ignore_index=True)
                    break
            
    
    extracted_sumstats = pd.concat([extracted_sumstats, ld_proxies],ignore_index=True)

    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return extracted_sumstats


def _extract_ld_proxy(  snplist=None,
                        common_sumstats=None,
                        vcf_path=None, 
                        vcf_chr_dict=None,
                        tabix=None,
                        log=Log(), 
                        verbose=True, 
                        windowsizekb=100,
                        ld_threshold=0.8
                            ):
    ### Load vcf#######################################################################################
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)
    if tabix is None:
        tabix = which("tabix")
    vcf_chr_dict = auto_check_vcf_chr_dict(vcf_path=vcf_path, vcf_chr_dict=vcf_chr_dict, verbose=verbose, log=log)
    
    ld_proxies  = pd.DataFrame()
    in_sumstats = common_sumstats.loc[common_sumstats["SNPID"].isin(snplist),:]
    
    if len(in_sumstats)==0:
        log.write(" -No available variants for LD proxy checking...Skipping... ", verbose=verbose)
    else:
        log.write(" -{} available variants for LD proxy checking... ".format(len(in_sumstats)), verbose=verbose)

    for index,row in in_sumstats.iterrows():
        # determine SNP and select region
        snpid = row["SNPID"]
        chrom= int(row["CHR"])
        start= int(row["POS"]-windowsizekb*1000)
        end=   int(row["POS"]+windowsizekb*1000)

        region = (chrom, start, end)
        
        ###  #######################################################################################
        #is_flanking = common_sumstats["CHR"] == chrom & common_sumstats["CHR"]>start & common_sumstats["CHR"]<end
        #flanking_sumstats = common_sumstats.loc[is_flanking,:]
        flanking_sumstats = common_sumstats.query('CHR == @chrom and @start < POS < @end',engine='python').copy()
        
        log.write(" -Extract {} variants in flanking region of {} for checking: {}:{}-{}".format(len(flanking_sumstats), snpid, chrom, start, end), verbose=verbose)

        if len(flanking_sumstats)==0:
            log.write("  -No availble variants in the region...Skipping!", verbose=verbose)
            continue

        flanking_sumstats = _get_rsq(row =in_sumstats.loc[index,["POS","NEA","EA"]],
                                     sumstats = flanking_sumstats, 
                                     row_pos=row["POS"], 
                                     vcf_path=vcf_path, 
                                     region=region,
                                     log=log, 
                                     verbose=verbose, 
                                     vcf_chr_dict=vcf_chr_dict, 
                                     tabix=tabix)
        if flanking_sumstats is None:
            log.write("  -{} is not found in the vcf...Skipping!".format(snpid))
            continue
        flanking_sumstats = flanking_sumstats.loc[flanking_sumstats["RSQ"]>ld_threshold,:]
        
        log.write("  -Variants in LD with {} (RSQ > {}): {}".format(snpid, ld_threshold,len(flanking_sumstats)), verbose=verbose)
        
        if len(flanking_sumstats)>0:
            flanking_sumstats["LD_REF_VARIANT"]= snpid
            for i,row_with_rsq in flanking_sumstats.iterrows():
                if row_with_rsq["SNPID"] in common_sumstats["SNPID"].values:
                    log.write("  -Top Proxy for {} is found: {} (LD RSQ= {})".format(snpid, row_with_rsq["SNPID"], row_with_rsq["RSQ"]))
                    break
            #row_with_rsq = pd.DataFrame(row_with_rsq)
            ld_proxies = pd.concat([ld_proxies, flanking_sumstats], ignore_index=True)
                    

    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return ld_proxies.sort_values(by="RSQ",ascending=False)


def _get_rsq(    row,
                 sumstats,
                 row_pos,
                 vcf_path, 
                 region,
                 log, 
                 verbose, 
                 vcf_chr_dict,
                 tabix):
        #load genotype data of the targeted region
        ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
        
        if ref_genotype is None:
            log.warning("No data was retrieved. Skipping ...", verbose=verbose)
            ref_genotype=dict()
            ref_genotype["variants/POS"]=np.array([],dtype="int64")
            return None
        
        log.write("  -Retrieving index...", verbose=verbose)
        log.write("  -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
        #  match sumstats pos and ref pos: 
        # get ref index for its first appearance of sumstats pos
        #######################################################################################
        def match_varaint(x):
            # x: "POS,NEA,EA"
            if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
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
                        else:
                            return None
                else: 
                    # single match
                    return np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0][0] 
            else:
                # no position match
                return None
        log.write("  -Matching variants using POS, NEA, EA ...", verbose=verbose)

        sumstats["REFINDEX"] = sumstats.loc[:,["POS","NEA","EA"]].apply(lambda x: match_varaint(x), axis=1)
        log.write("  -Matched variants in sumstats and vcf:{} ".format(sum(~sumstats["REFINDEX"].isna())))
        #############################################################################################
        lead_pos = row_pos

        # if lead pos is available: 
        if lead_pos in ref_genotype["variants/POS"]:
            
            # get ref index for lead snp
            lead_snp_ref_index = match_varaint(row)
            #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

            # non-na other snp index
            other_snps_ref_index = sumstats["REFINDEX"].dropna().astype("int").values
            # get genotype 
            lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
            other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
            
            log.write("  -Calculating Rsq...", verbose=verbose)
            
            if len(other_snp_genotype)>1:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
            else:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)
            sumstats.loc[~sumstats["REFINDEX"].isna(),"RSQ"] = valid_r2
        else:
            log.write("  -Lead SNP not found in reference...", verbose=verbose)
            sumstats["RSQ"]=None
        
        sumstats["RSQ"] = sumstats["RSQ"].astype("float")
        return sumstats

def _check_if_in_sumstats2(row, sumstast):
    pass

def _get_rsq_single(  row,
                 row_pos,
                 vcf_path, 
                 region,
                 log, 
                 verbose, 
                 vcf_chr_dict,
                 tabix):
    #load genotype data of the targeted region
    ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
    
    if ref_genotype is None:
        log.warning("No data was retrieved. Skipping ...", verbose=verbose)
        ref_genotype=dict()
        ref_genotype["variants/POS"]=np.array([],dtype="int64")
        return None
    
    log.write("  -Retrieving index...", verbose=verbose)
    log.write("  -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
    #  match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
    #######################################################################################
    def match_varaint(x):
        # x: "POS,NEA,EA"
        if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
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
                    else:
                        return None
            else: 
                # single match
                return np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0][0] 
        else:
            # no position match
            return None

    #############################################################################################
    lead_pos = row_pos

    # if lead pos is available: 
    if lead_pos in ref_genotype["variants/POS"]:
        
        # get ref index for lead snp
        lead_snp_ref_index = match_varaint(row)
        #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

        # non-na other snp index
        other_snps_ref_index = list(range(len(ref_genotype["calldata/GT"])))
        other_snps_ref_index.remove(lead_snp_ref_index)

        # get genotype 
        lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
        other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
        
        log.write("  -Calculating Rsq...", verbose=verbose)
        
        if len(other_snp_genotype)>1:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
        else:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)

        ld_proxy = pd.DataFrame( {"SNPID":ref_genotype["variants/ID"][other_snps_ref_index],"RSQ":valid_r2 })

    return  ld_proxy.sort_values(by="RSQ",ascending=False)