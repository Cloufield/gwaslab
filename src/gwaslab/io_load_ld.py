
import scipy.sparse as sparse
import numpy as np
import pandas as pd
from gwaslab.hm_casting import _merge_mold_with_sumstats_by_chrpos
import subprocess
import os
import re
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_ex_process_ref import _process_plink_input_files
from gwaslab.g_version import _checking_plink_version
from gwaslab.util_in_filter_value import _exclude_hla
from gwaslab.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.util_ex_calculate_ldmatrix import _export_snplist_and_locus_sumstats
from gwaslab.viz_plot_regional2 import _get_lead_id
from gwaslab.util_ex_calculate_ldmatrix import _extract_variants_in_locus

def tofinemapping_using_ld(sumstats, 
                  study=None, 
                  ld_map_path=None, 
                  ld_path=None,
                  ld_fmt = "npz",
                  ld_if_square =  False,
                  ld_if_add_T = False,
                  ld_map_rename_dic = None,
                  ld_map_kwargs = None,
                  loci=None,
                  out="./",
                  windowsizekb=1000,
                  n_cores=1, 
                  mode="r",
                  exclude_hla=False, 
                  getlead_args=None, 
                  memory=None, 
                  overwrite=False,
                  log=Log(),
                  suffixes=None,
                  verbose=True,
                  **kwargs):
    ##start function with col checking##########################################################
    _start_line = "calculate LD matrix"
    _end_line = "calculating LD matrix"
    _start_cols =["SNPID","CHR","POS","EA","NEA"]
    _start_function = ".calculate_ld_matrix()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: raise ValueError("Not enough columns for calculating LD matrix")
    ############################################################################################
    if suffixes is None:
        suffixes=[""]
    if getlead_args is None:
        getlead_args={"windowsizekb":1000}
    if ld_map_kwargs is None:
        ld_map_kwargs={}
    
    if loci is None:
        log.write(" -Loci were not provided. All significant loci will be automatically extracted...",verbose=verbose)
        sig_df = getsig(sumstats,id="SNPID",chrom="CHR",pos="POS",p="P"+suffixes[0],**getlead_args)
    else:
        sig_df = sumstats.loc[sumstats["SNPID"].isin(loci),:]

    # Drop duplicate!!!!
    log.write(" -Dropping duplicated SNPIDs...",verbose=verbose)
    sumstats = sumstats.drop_duplicates(subset=["SNPID"]).copy()

    # init Filelist DataFrame
    output_file_list = pd.DataFrame(columns=["SNPID","SNPID_LIST","LD_R_MATRIX","LOCUS_SUMSTATS"])
    
    plink_log=""

    if exclude_hla==True:
        sig_df = _exclude_hla(sig_df, log=log, verbose=verbose)
    
    sig_df = sig_df.reset_index()
    
    ## for each lead variant 
    for index, row in sig_df.iterrows():
        # extract snplist in each locus
        gc.collect()
        log.write(" -Locus #{}---------------------------------------------------------------".format(index+1))
        log.write(" -Processing locus with lead variant {} at CHR {} POS {} ...".format(row["SNPID"],row["CHR"],row["POS"]))
        locus_sumstats = _extract_variants_in_locus(sumstats, windowsizekb, locus = (row["CHR"],row["POS"]))
        
        ld_map = _load_ld_map(ld_map_path, ld_map_rename_dic = ld_map_rename_dic, **ld_map_kwargs )

        ## check available snps with reference file
        matched_sumstats = _merge_ld_map_with_sumstats(row=row, 
                                                    locus_sumstats=locus_sumstats, 
                                                    ld_map=ld_map,
                                                    log=log,suffixes=suffixes)
        if len(matched_sumstats)==0:
            log.write(" -No matching LD information... Skipping...")
            continue
        
        #########################################################################################################
        # create matched snp list
        matched_snp_list_path, matched_sumstats_path=_export_snplist_and_locus_sumstats(matched_sumstats=matched_sumstats, 
                                                                                       out=out, 
                                                                                       study=study, 
                                                                                       row=row, 
                                                                                       windowsizekb=windowsizekb,
                                                                                       log=log,
                                                                                       suffixes=suffixes)
        #########################################################################################################

        ## Calculate ld matrix using PLINK
        r_matrix = _load_ld_matrix(ld_path, fmt=ld_fmt, if_square=ld_if_square, if_add_T=ld_if_add_T, log=log, verbose=verbose)
        
        matched_ld_matrix_path = _extract_variants(matched_sumstats, r_matrix, out, study, row, windowsizekb, log=log, verbose=verbose)

        # print file list
        row_dict={}
        row_dict["SNPID"]=row["SNPID"]
        row_dict["SNPID_LIST"] = matched_snp_list_path
        row_dict["LD_R_MATRIX"] = matched_ld_matrix_path
        row_dict["LOCUS_SUMSTATS"] = matched_sumstats_path
        file_row = pd.Series(row_dict).to_frame().T
        output_file_list = pd.concat([output_file_list, file_row],ignore_index=True)
    
    if len(output_file_list)>0:
        output_file_list["STUDY"] = study
        nloci = len(output_file_list)
        output_file_list_path =  "{}/{}_{}loci_{}kb.filelist".format(out.rstrip("/"), study,nloci, windowsizekb)
        output_file_list.to_csv(output_file_list_path,index=None,sep="\t")
        log.write(" -File list is saved to: {}".format(output_file_list_path),verbose=verbose)
        log.write(" -Finished LD matrix calculation.",verbose=verbose)
    else:
        output_file_list_path=None
        log.write(" -No avaialable lead variants.",verbose=verbose)
        log.write(" -Stopped LD matrix calculation.",verbose=verbose)
    finished(log=log, verbose=verbose, end_line=_end_line)
    return output_file_list_path, output_file_list, plink_log




def process_ld(sumstats, 
                ld_path, 
                ld_map_path,
                region,
                region_ref, 
                log, 
                verbose, 
                pos ,
                nea,
                ea, 
                region_ld_threshold, 
                ld_fmt = "npz",
                ld_if_square =  False,
                ld_if_add_T = False,
                ld_map_rename_dic = None,
                ld_map_kwargs = None):
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference ld matrix path : "+ ld_path, verbose=verbose)

    # load genotype data of the targeted region


    log.write(" -Retrieving index...", verbose=verbose)
    
    # match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
     #######################################################################################
    if ld_map_kwargs is None:
        ld_map_kwargs={}

    ld_map = _load_ld_map(ld_map_path, 
                          ld_map_rename_dic = ld_map_rename_dic, 
                          **ld_map_kwargs )
    
    log.write(" -Ref variants: {}".format( len(ld_map) ), verbose=verbose)

    ## check available snps with reference file
    sumstats = _merge_ld_map_with_sumstats_for_regional(
                                           locus_sumstats=sumstats, 
                                           ld_map=ld_map,
                                           log=log,
                                           suffixes=None,verbose=verbose)
    sumstats["REFINDEX"] = sumstats["_INDEX_BIM"]

    #############################################################################################
    
    r_matrix = _load_ld_matrix(ld_path, 
                               fmt=ld_fmt, 
                               if_square=ld_if_square, 
                               if_add_T=ld_if_add_T, 
                               log=log, 
                               verbose=verbose)

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
                lead_snpid[0]= int(lead_snpid[0])
                lead_series = pd.Series(lead_snpid)

        # if lead pos is available: 
        if sumstats.loc[lead_id, "REFINDEX"] is not None:
            lead_snp_ref_index = sumstats.loc[lead_id, "REFINDEX"]

            is_matched = ~sumstats["REFINDEX"].isna()

            ref_index = sumstats.loc[is_matched,"REFINDEX"].astype("Int64")
            
            sumstats.loc[is_matched, rsq] = r_matrix[int(lead_snp_ref_index), list(ref_index.values)]

        else:
            log.write(" -Lead SNP not found in reference...", verbose=verbose)
            sumstats[rsq]=None
            # 
            try:
                sumstats.loc[lead_id,rsq]=1
            except KeyError:
                pass
        
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

    if len(region_ref)==1:
        if lead_id is not None:
            sumstats.loc[lead_id, final_shape_col] +=1 

    for i in range(len(region_ref)):
        ld_single = "LD_{}".format(i)
        current_rsq = "RSQ_{}".format(i)
        a_ngt_b = sumstats[final_rsq_col] < sumstats[current_rsq]
        #set levels with interval=100
        sumstats.loc[a_ngt_b, final_ld_col] = 100 * (i+1) + sumstats.loc[a_ngt_b, ld_single]
        sumstats.loc[a_ngt_b, final_rsq_col] = sumstats.loc[a_ngt_b, current_rsq]
        sumstats.loc[a_ngt_b, final_shape_col] = i + 1
    
    sumstats = sumstats.dropna(subset=[pos,nea,ea])

    ####################################################################################################
    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return sumstats
####################################################################################################
































####################################################################################################
def _load_ld_matrix(path, 
                    fmt="npz", 
                    if_square=False, 
                    if_add_T=False,
                    log=Log(),
                    verbose=True):
    
    if fmt == "npz":
        log.write("   -Loading LD matrix from npz file...",verbose=verbose)
        r_matrix = sparse.load_npz(path).toarray()
    if fmt == "txt":
        log.write("   -Loading LD matrix from text file...",verbose=verbose)
        r_matrix = np.loadtxt(path)

    if if_add_T==True:
        log.write("   -Transforming LD matrix by adding its transpose...",verbose=verbose)
        r_matrix += r_matrix.T
    if if_square==True:
        log.write("   -Transforming LD matrix by squaring all elements...",verbose=verbose)
        r_matrix = np.power(r_matrix,2)
    return r_matrix
    
def _load_ld_map(path, 
                 snpid="rsid", 
                 chrom="chromosome", 
                 pos="position", 
                 ref="allele1", 
                 alt="allele2",
                 ld_map_rename_dic = None,
                 **ld_map_kwargs):
    
    if ld_map_rename_dic is not None:
        if type(ld_map_rename_dic) is dict:
            ld_map_rename_dic_to_use={ld_map_rename_dic["EA"]:'EA_bim', 
                                        ld_map_rename_dic["NEA"]:'NEA_bim', 
                                        ld_map_rename_dic["POS"]:'POS', 
                                        ld_map_rename_dic["CHR"]:'CHR',
                                        ld_map_rename_dic["SNPID"]:'SNPID_bim'
                                        }
            ld_map_kwargs["usecols"]=list(ld_map_rename_dic.values())
        else:
            ld_map_rename_dic_to_use={ld_map_rename_dic[4]:'EA_bim', 
                                        ld_map_rename_dic[3]:'NEA_bim', 
                                        ld_map_rename_dic[2]:'POS', 
                                        ld_map_rename_dic[1]:'CHR',
                                        ld_map_rename_dic[0]:'SNPID_bim'
                                        }
            ld_map_kwargs["usecols"]=ld_map_rename_dic
    else:
        ld_map_rename_dic_to_use={alt:'EA_bim', 
                                  ref:'NEA_bim', 
                                  pos:'POS', 
                                  chrom:'CHR',
                                  snpid:"SNPID_bim"
                                    }
        ld_map_kwargs["usecols"]=[chrom, pos, ref, alt, snpid]
    #rsid    chromosome      position        allele1 allele2
    if "sep" not in ld_map_kwargs:
        ld_map_kwargs["sep"] = "\s+"
    
    ld_map = pd.read_csv(path,**ld_map_kwargs)
    ld_map = ld_map.rename(columns=ld_map_rename_dic_to_use, errors='ignore')
    # "SNPID",0:"CHR_bim",3:"POS_bim",4:"EA_bim",5:"NEA_bim"
    return ld_map

def _extract_variants(merged_sumstats, r_matrix, out, study, row, windowsizekb, log, verbose):
    
    avaiable_index = merged_sumstats["_INDEX_BIM"].values 

    flipped = merged_sumstats["_FLIPPED"].values 

    reduced_r_matrix = r_matrix[np.ix_(avaiable_index, avaiable_index)]
    
    log.write(" -Flipping LD matrix for {} variants...".format(sum(flipped)),verbose=verbose)
    reduced_r_matrix[flipped,:] = -1 * reduced_r_matrix[flipped,:]
    reduced_r_matrix[:,flipped] = -1 * reduced_r_matrix[:,flipped]

    snplist_path =   "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
    output_prefix =  "{}/{}_{}_{}".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
    output_path = "{}.ld.gz".format(output_prefix)
    
    pd.DataFrame(reduced_r_matrix).to_csv(output_path,sep="\t",index=None,header=None)
    #reduced_r_matrix.to_csv("{}.ld.gz".format(output_prefix),se="\t")
    return output_path

def _merge_ld_map_with_sumstats(row, 
                             locus_sumstats, 
                             ld_map, 
                             log=Log(),
                             suffixes=None):
    '''
    align sumstats with bim
    '''

    index1= "_INDEX_SUMSTATS"
    index2= "_INDEX_BIM"
    locus_sumstats[index1] = locus_sumstats.index
    ld_map[index2] =  ld_map.index
    locus_sumstats["_FLIPPED"] = False
    
    if suffixes is None:
            suffixes=[""]
    
    log.write("   -Variants in locus ({}): {}".format(row["SNPID"],len(locus_sumstats)))
    # convert category to string
    locus_sumstats["EA"] = locus_sumstats["EA"].astype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].astype("string")
    
    # matching by SNPID
    # preserve bim keys (use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.)
    combined_df = pd.merge(ld_map, locus_sumstats, on=["CHR","POS"],how="inner")
    
    # match allele
    perfect_match =  ((combined_df["EA"] == combined_df["EA_bim"]) & (combined_df["NEA"] == combined_df["NEA_bim"]) ) 
    log.write("   -Variants with perfect matched alleles:{}".format(sum(perfect_match)))

    # fliipped allele
    #ea_mis_match = combined_df["EA"] != combined_df["EA_bim"]
    flipped_match = ((combined_df["EA"] == combined_df["NEA_bim"])& (combined_df["NEA"] == combined_df["EA_bim"]))
    log.write("   -Variants with flipped alleles:{}".format(sum(flipped_match)))
    
    allele_match = perfect_match | flipped_match
    log.write("   -Total Variants matched:{}".format(sum(allele_match)))

    if row["SNPID"] not in combined_df.loc[allele_match,"SNPID"].values:
        log.warning("Lead variant was not available in reference!")
    
    # adjust statistics
    output_columns=["SNPID","CHR","POS","EA","NEA","_INDEX_BIM","_FLIPPED"]
    for suffix in suffixes:
        if ("BETA"+suffix in locus_sumstats.columns) and ("SE"+suffix in locus_sumstats.columns):
            #log.write("   -Flipping BETA{} for variants with flipped alleles...".format(suffix))
            #combined_df.loc[flipped_match,"BETA"+suffix] = - combined_df.loc[flipped_match,"BETA"+suffix]
            output_columns.append("BETA"+suffix)
            output_columns.append("SE"+suffix)
        if "Z" in locus_sumstats.columns:
            #log.write("   -Flipping Z{} for variants with flipped alleles...".format(suffix))
            #combined_df.loc[flipped_match,"Z"+suffix] = - combined_df.loc[flipped_match,"Z"+suffix]
            output_columns.append("Z"+suffix)
        if "EAF" in locus_sumstats.columns:
            #log.write("   -Flipping EAF{} for variants with flipped alleles...".format(suffix))
            #combined_df.loc[flipped_match,"EAF"+suffix] = 1 - combined_df.loc[flipped_match,"EAF"+suffix]
            output_columns.append("EAF"+suffix)
        if "N" in locus_sumstats.columns:
            output_columns.append("N"+suffix)
    combined_df.loc[flipped_match,"_FLIPPED"] = True
    return combined_df.loc[allele_match,output_columns]

def _merge_ld_map_with_sumstats_for_regional(
                             locus_sumstats, 
                             ld_map, 
                             log=Log(),
                             suffixes=None,
                             verbose=True):
    '''
    align sumstats with bim
    '''

    index1= "_INDEX_SUMSTATS"
    index2= "_INDEX_BIM"
    locus_sumstats[index1] = locus_sumstats.index
    ld_map[index2] =  ld_map.index

    if suffixes is None:
            suffixes=[""]
    
    # convert category to string
    locus_sumstats["EA"] = locus_sumstats["EA"].astype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].astype("string")
    
    # matching by SNPID
    # preserve bim keys (use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.)
    combined_df = pd.merge(locus_sumstats, ld_map, on=["CHR","POS"],how="left")
    combined_df[["EA_bim","NEA_bim"]] = combined_df[["EA_bim","NEA_bim"]].fillna("N")
    # match allele
    perfect_match =  ((combined_df["EA"] == combined_df["EA_bim"]) & (combined_df["NEA"] == combined_df["NEA_bim"]) ) 

    # fliipped allele
    #ea_mis_match = combined_df["EA"] != combined_df["EA_bim"]
    flipped_match = ((combined_df["EA"] == combined_df["NEA_bim"])& (combined_df["NEA"] == combined_df["EA_bim"]))
    
    not_matched = combined_df[index2].isna() 

    allele_match = perfect_match | flipped_match 

    log.write("   -Total Variants matched:{}".format( sum(allele_match) ),verbose=verbose)
    log.write("   -Total Variants not in reference:{}".format(sum(not_matched)),verbose=verbose)
    
    return combined_df.loc[allele_match | not_matched,:]

############################################################################################################################################################################################################################################################