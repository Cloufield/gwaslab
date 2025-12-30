import scipy.sparse as sparse
import numpy as np
import pandas as pd

import subprocess
import os
import re
import gc
from typing import TYPE_CHECKING, Optional, List, Dict, Any, Tuple, Union

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_plink_version

from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos

from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.io.io_plink import _process_plink_input_files
from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.util.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.qc.qc_decorator import with_logging

@with_logging(
        start_to_msg="calculate LD matrix",
        finished_msg="calculating LD matrix",
        start_cols=["SNPID","CHR","POS","EA","NEA"],
        start_function=".calculate_ld_matrix()"
)
def tofinemapping_m(
    sumstats: pd.DataFrame, 
    studies: Optional[List[str]] = None, 
    group: Optional[str] = None,
    ld_paths: Optional[List[str]] = None,
    ld_types: Optional[List[str]] = None, 
    ld_maps: Optional[List[str]] = None,
    ld_map_dics: Optional[List[Dict[str, str]]] = None,
    bfile: Optional[str] = None, 
    vcf: Optional[str] = None, 
    locus: Optional[Tuple[int, int]] = None,
    loci: Optional[List[str]] = None,
    loci_chrpos: Optional[List[Tuple[int, int]]] = None,
    out: str = "./",
    plink: str = "plink",
    plink2: str = "plink2",
    windowsizekb: int = 1000,
    threads: int = 1, 
    mode: str = "r",
    exclude_hla: bool = False, 
    getlead_kwargs: Optional[Dict[str, Any]] = None, 
    memory: Optional[int] = None, 
    overwrite: bool = False,
    log: Log = Log(),
    suffixes: Optional[List[str]] = None,
    ld_map_kwargs: Optional[Dict[str, Any]] = None,
    extra_plink_option: str = "",
    verbose: bool = True,
    **kwargs: Any
) -> Tuple[Optional[str], pd.DataFrame, str]:
    
    ############################################################################################
    if suffixes is None:
        suffixes=[""]
    if getlead_kwargs is None:
        getlead_kwargs={"windowsizekb":1000}
    if ld_map_kwargs is None:
        ld_map_kwargs={}
    
    if loci is None:
        log.write(" -Loci were not provided. All significant loci will be automatically extracted...",verbose=verbose)
        sig_df = _get_sig(sumstats,variant_id="SNPID",chrom="CHR",pos="POS",p="P"+suffixes[0],**getlead_kwargs)
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
    row = sig_df.iloc[0,:]

    matched_sumstats = _extract_variants_in_locus(sumstats, windowsizekb, locus = (row["CHR"],row["POS"]))

    for i in range(2):
    # for each study
        gc.collect()
        # get ld path and dic
        ld_map_path = ld_maps[i]
        ld_map_rename_dic =  ld_map_dics[i]
        
        
        log.write(" -Processing locus with lead variant {} at CHR {} POS {} ...".format(row["SNPID"],row["CHR"],row["POS"]))
        ld_map = _load_ld_map(ld_map_path, 
                              ld_map_rename_dic = ld_map_rename_dic, 
                              **ld_map_kwargs )

        ## check available snps with reference file
        matched_sumstats = _merge_ld_map_with_sumstats_m(row=row, 
                                                        locus_sumstats=matched_sumstats, 
                                                        ld_map=ld_map,
                                                        log=log,
                                                        index=i)
        
        if len(matched_sumstats)==0:
            log.write(" -No matching LD information... Skipping...")
            continue
    
    # drop na
    matched_sumstats = matched_sumstats.dropna()

    # export common variants list
    matched_snp_list_path, matched_sumstats_paths=_export_snplist_and_locus_sumstats_m(matched_sumstats=matched_sumstats, 
                                                                                    out=out, 
                                                                                    group=group, 
                                                                                    row=row, 
                                                                                    windowsizekb=windowsizekb,
                                                                                    log=log)  
    
    for i in range(2):
        ld_path = ld_paths[i]

        r_matrix = _load_ld_matrix(ld_path, 
                                   fmt="txt", 
                                   if_square=False, 
                                   if_add_T=False, 
                                   log=log, 
                                   verbose=verbose)
        
        matched_ld_matrix_path = _extract_variants_from_ld_matrix_m(merged_sumstats = matched_sumstats, 
                                                   r_matrix = r_matrix, 
                                                   out = out, 
                                                   group = group, 
                                                   row = row, 
                                                   windowsizekb = windowsizekb, 
                                                   index=i,
                                                   log=log, verbose=verbose)
#     #########################################################################################################

        row_dict={}
        row_dict["SNPID"]=row["SNPID"]
        row_dict["SNPID_LIST"] = matched_snp_list_path
        row_dict["LD_R_MATRIX"] = matched_ld_matrix_path
        row_dict["LOCUS_SUMSTATS"] = matched_sumstats_paths[i] + ".gz"
        row_dict["LOCUS"] = row["SNPID"]
        row_dict["SUBSTUDY"]= i+1
        row_dict["STUDY"] = studies[i]
        file_row = pd.Series(row_dict).to_frame().T
        output_file_list = pd.concat([output_file_list, file_row],ignore_index=True)

    if len(output_file_list)>0:
        
        output_file_list["GROUP"] = group
        nloci = len(output_file_list)
        output_file_list_path =  "{}/{}_{}study_{}_{}kb.filelist".format(out.rstrip("/"), group,nloci, row["SNPID"], windowsizekb)
        output_file_list.to_csv(output_file_list_path,index=None,sep="\t")
        log.write(" -File list is saved to: {}".format(output_file_list_path),verbose=verbose)
        log.write(" -Finished LD matrix calculation.",verbose=verbose)
    else:
        output_file_list_path=None
        log.write(" -No avaialable lead variants.",verbose=verbose)
        log.write(" -Stopped LD matrix calculation.",verbose=verbose)


    return output_file_list_path, output_file_list,  plink_log


def _export_snplist_and_locus_sumstats_m(
    matched_sumstats: pd.DataFrame, 
    out: str, 
    group: str, 
    row: pd.Series, 
    windowsizekb: int, 
    log: Log
) -> Tuple[str, List[str]]:
        # study suffixes starting from 1
        suffixes=["_{}".format(i+1) for i in range(2)]
        
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), group, row["SNPID"] ,windowsizekb)
        
        matched_sumstats["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)
        log.write(" -Exporting SNP list of {}  to: {}...".format(len(matched_sumstats) ,matched_snp_list_path))

        # create locus-sumstats EA, NEA, (BETA, SE), Z 
        matched_sumstats_paths =[]

        
        for i in range(2):
            # export sumstats for each study
            suffix = suffixes[i]
            
            matched_sumstats_path =  "{}/{}_{}_{}_{}.sumstats".format(out.rstrip("/"), group, row["SNPID"] ,windowsizekb, i + 1)
            matched_sumstats_paths.append(matched_sumstats_path)
            to_export_columns=["CHR","POS","EA","NEA"]

            if "Z"+suffix in matched_sumstats.columns :
                to_export_columns.append("Z"+suffix)
            if ("BETA"+suffix in matched_sumstats.columns) and ("SE"+suffix in matched_sumstats.columns):
                to_export_columns.append("BETA"+suffix)
                to_export_columns.append("SE"+suffix)
            if "EAF"+suffix in matched_sumstats.columns :
                to_export_columns.append("EAF"+suffix)
            if "N"+suffix in matched_sumstats.columns:
                to_export_columns.append("N"+suffix)

            log.write(" -Exporting locus sumstats to: {}...".format(matched_sumstats_path))
            log.write(" -Exported columns: {}...".format(["SNPID"]+to_export_columns))
            #matched_sumstats[ ["SNPID"]+to_export_columns].to_csv(matched_sumstats_path, sep="\t",index=None)
            rename_dic={
                "BETA"+suffix:"Beta",
                "SE"+suffix:"Se",
                "SNPID":"SNP"
            }
            matched_sumstats[ ["SNPID"]+to_export_columns].rename(columns=rename_dic).to_csv(matched_sumstats_path, sep="\t",index=None)
            matched_sumstats[ ["SNPID"]+to_export_columns].rename(columns=rename_dic).to_csv(matched_sumstats_path+".gz", sep="\t",index=None)
        
        return matched_snp_list_path, matched_sumstats_paths

###################################################################################################################################################################
####################################################################################################
def _load_ld_matrix(
    path: str, 
    fmt: str = "npz", 
    if_square: bool = False, 
    if_add_T: bool = False,
    log: Log = Log(),
    verbose: bool = True
) -> np.ndarray:
    
    if fmt == "npz":
        log.write("   -Loading LD matrix from npz file...",verbose=verbose)
        r_matrix = sparse.load_npz(path).toarray()
    if fmt == "txt":
        log.write("   -Loading LD matrix from text file...",verbose=verbose)
        r_matrix = np.loadtxt(path,delimiter="\t")
        log.write("   -LD matrix shape : {}".format(r_matrix.shape) ,verbose=verbose)

    if if_add_T==True:
        log.write("   -Transforming LD matrix by adding its transpose...",verbose=verbose)
        r_matrix += r_matrix.T
    if if_square==True:
        log.write("   -Transforming LD matrix by squaring all elements...",verbose=verbose)
        r_matrix = np.power(r_matrix,2)
    return r_matrix
    
def _load_ld_map(
    path: str, 
    snpid: str = "rsid", 
    chrom: str = "chromosome", 
    pos: str = "position", 
    ref: str = "allele1", 
    alt: str = "allele2",
    ld_map_rename_dic: Optional[Union[Dict[str, str], List[str]]] = None,
    **ld_map_kwargs: Any
) -> pd.DataFrame:
    
    if ld_map_rename_dic is not None:
        # ld map format
        # SNPID_bim,CHR，POS, NEA_bim, EA_bim 
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

def _extract_variants_from_ld_matrix_m(
    merged_sumstats: pd.DataFrame, 
    r_matrix: np.ndarray, 
    out: str, 
    group: str, 
    row: pd.Series, 
    windowsizekb: int, 
    log: Log, 
    verbose: bool, 
    index: int
) -> str:
    # study suffixes starting from 1
    index_bim_header = "_INDEX_BIM_{}".format(index + 1) 
    flipped_header = "_FLIPPED_{}".format(index + 1) 
    
    # a series of int to indicate if the variant index in raw ld matrix
    avaiable_index = merged_sumstats[index_bim_header].values 

    # a series of boolean values to indicate if the variants is flipped
    flipped = merged_sumstats[flipped_header].values 

    # extract the sub-matrix
    reduced_r_matrix = r_matrix[np.ix_(avaiable_index, avaiable_index)]
    
    log.write(" -Flipping LD matrix for {} variants...".format(sum(flipped)),verbose=verbose)
    reduced_r_matrix[flipped,:] = -1 * reduced_r_matrix[flipped,:]
    reduced_r_matrix[:,flipped] = -1 * reduced_r_matrix[:,flipped]

    output_prefix =  "{}/{}_{}_{}_{}".format(out.rstrip("/"),group,row["SNPID"],windowsizekb, index + 1)
    output_path = "{}.ld.gz".format(output_prefix)
    
    pd.DataFrame(reduced_r_matrix).to_csv(output_path,sep="\t",index=None,header=None)
    #reduced_r_matrix.to_csv("{}.ld.gz".format(output_prefix),se="\t")
    return output_path

def _merge_ld_map_with_sumstats_m(
    row: pd.Series, 
    locus_sumstats: pd.DataFrame, 
    ld_map: pd.DataFrame, 
    log: Log = Log(),
    index: Optional[int] = None
) -> pd.DataFrame:
    '''
    align sumstats with ld map
    '''
    # study suffixes starting from 1
    index_suffix = "_{}".format(index+1)
    
    index1= "_INDEX_SUMSTATS" 
    index2= "_INDEX_BIM" +index_suffix
    
    # Sumstats index
    locus_sumstats[index1] = locus_sumstats.index
    
    # ld map index
    ld_map[index2] =  ld_map.index
    
    # init a column to show if the variants in LD map are flipped or not
    locus_sumstats["_FLIPPED"+index_suffix] = False
    
    
    log.write("   -Variants in locus ({}): {}".format(row["SNPID"],len(locus_sumstats)))
    # convert category to string
    locus_sumstats["EA"] = locus_sumstats["EA"].astype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].astype("string")
    
    # matching by SNPID
    # preserve bim keys (use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.)
    # variants without a match were removed
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
    
    combined_df.loc[flipped_match,"_FLIPPED"+index_suffix] = True

    if row["SNPID"] not in combined_df.loc[allele_match,"SNPID"].values:
        log.warning("Lead variant was not available in reference!")
    
    # adjust output columns
    output_columns=["SNPID","CHR","POS","EA","NEA"]
    for i in combined_df.columns:
        if "_INDEX_BIM" in i:
            output_columns.append(i)
        if "_FLIPPED" in i:
            output_columns.append(i)
    
    for i in range(2):
        # study suffixes starting from 1
        index_suffix = "_{}".format(i+1)
        if ("BETA"+index_suffix in combined_df.columns) and ("SE"+index_suffix in combined_df.columns):
            output_columns.append("BETA"+index_suffix)
            output_columns.append("SE"+index_suffix)
        if "Z"+index_suffix in combined_df.columns:
            output_columns.append("Z"+index_suffix)
        if "EAF"+index_suffix in combined_df.columns:
            output_columns.append("EAF"+index_suffix)
        if "N"+index_suffix in combined_df.columns:
            output_columns.append("N"+index_suffix)
    
    return combined_df.loc[allele_match,output_columns]
