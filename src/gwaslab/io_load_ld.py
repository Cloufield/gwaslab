
import scipy.sparse as sparse
import numpy as np
import pandas as pd
from gwaslab.hm_casting import _merge_mold_with_sumstats_by_chrpos
import subprocess
import os
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

from gwaslab.util_ex_calculate_ldmatrix import _extract_variants_in_locus


def _load_ld_matrix(path, 
                    fmt="npz", 
                    if_square=False, 
                    if_add_T=False):
    
    if fmt == "npz":
        r_matrix = sparse.load_npz(path).toarray()
    if fmt == "txt":
        r_matrix = np.loadtxt(path,)

    if if_add_T==True:
        r_matrix += r_matrix.T
    if if_square==True:
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
        ld_map_rename_dic_to_use={ld_map_rename_dic["EA"]:'EA_bim', 
                                    ld_map_rename_dic["NEA"]:'NEA_bim', 
                                    ld_map_rename_dic["POS"]:'POS', 
                                    ld_map_rename_dic["CHR"]:'CHR',
                                    ld_map_rename_dic["SNPID"]:'SNPID_bim'
                                    }
        ld_map_kwargs["usecols"]=list(ld_map_rename_dic.values())
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

def _extract_variants(merged_sumstats, r_matrix, out, study, row, windowsizekb):
    
    avaiable_index= merged_sumstats["_INDEX_BIM"].values 
    
    reduced_r_matrix = r_matrix[np.ix_(avaiable_index, avaiable_index)]
    
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
    output_columns=["SNPID","CHR","POS","EA_bim","NEA_bim","_INDEX_BIM"]
    for suffix in suffixes:
        if ("BETA"+suffix in locus_sumstats.columns) and ("SE"+suffix in locus_sumstats.columns):
            log.write("   -Flipping BETA{} for variants with flipped alleles...".format(suffix))
            combined_df.loc[flipped_match,"BETA"+suffix] = - combined_df.loc[flipped_match,"BETA"+suffix]
            output_columns.append("BETA"+suffix)
            output_columns.append("SE"+suffix)
        if "Z" in locus_sumstats.columns:
            log.write("   -Flipping Z{} for variants with flipped alleles...".format(suffix))
            combined_df.loc[flipped_match,"Z"+suffix] = - combined_df.loc[flipped_match,"Z"+suffix]
            output_columns.append("Z"+suffix)
        if "EAF" in locus_sumstats.columns:
            log.write("   -Flipping EAF{} for variants with flipped alleles...".format(suffix))
            combined_df.loc[flipped_match,"EAF"+suffix] = 1 - combined_df.loc[flipped_match,"EAF"+suffix]
            output_columns.append("EAF"+suffix)
        if "N" in locus_sumstats.columns:
            output_columns.append("N"+suffix)
    
    return combined_df.loc[allele_match,output_columns]


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
        r_matrix = _load_ld_matrix(ld_path, fmt=ld_fmt, if_square=ld_if_square, if_add_T=ld_if_add_T)
        
        matched_ld_matrix_path = _extract_variants(matched_sumstats, r_matrix, out, study, row, windowsizekb)

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