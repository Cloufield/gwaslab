from typing import TYPE_CHECKING, Optional, List, Dict, Any, Union, Tuple
import scipy.sparse as sparse
import numpy as np
import pandas as pd
import subprocess
import os
import re
import gc
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
from gwaslab.extension import _checking_plink_version
from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.io.io_plink import _process_plink_input_files
from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.util.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.util.util_ex_calculate_ldmatrix import _export_snplist_and_locus_sumstats
from gwaslab.util.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.viz.viz_plot_regional2 import _get_lead_id
from gwaslab.qc.qc_decorator import with_logging
@with_logging(
        start_to_msg="calculate LD matrix",
        finished_msg="calculating LD matrix",
        start_cols=["SNPID","CHR","POS","EA","NEA"],
        start_function="calculate_ld_matrix"
)
def _to_finemapping_using_ld(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    study: Optional[str] = None,
    ld_map_path: Optional[str] = None,
    ld_path: Optional[str] = None,
    ld_fmt: str = "npz",
    ld_if_square: bool = False,
    ld_if_add_T: bool = False,
    ld_map_rename_dic: Optional[Union[Dict[str, str], List[str]]] = None,
    ld_map_kwargs: Optional[Dict[str, Any]] = None,
    loci: Optional[List[str]] = None,
    out: str = "./",
    windowsizekb: int = 1000,
    threads: int = 1,
    mode: str = "r",
    exclude_hla: bool = False,
    getlead_kwargs: Optional[Dict[str, Any]] = None,
    memory: Optional[int] = None,
    overwrite: bool = False,
    log: Log = Log(),
    suffixes: Optional[List[str]] = None,
    verbose: bool = True,
    **kwargs: Any
) -> pd.DataFrame:
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data

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
    log.write(" -Number of loci: {}...".format(len(sig_df)),verbose=verbose)
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
    return output_file_list_path, output_file_list, plink_log




def process_ld(
    sumstats: pd.DataFrame,
    ld_path: str,
    ld_map_path: str,
    region: Tuple[int, int, int],
    region_ref: List[Optional[str]],
    log: Log,
    verbose: bool,
    pos: str,
    nea: str,
    ea: str,
    region_ld_threshold: float,
    ld_fmt: str = "npz",
    ld_if_square: bool = False,
    ld_if_add_T: bool = False,
    ld_map_rename_dic: Optional[Union[Dict[str, str], List[str]]] = None,
    ld_map_kwargs: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
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
        r_matrix = np.loadtxt(path)

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

def _extract_variants(
    merged_sumstats: pd.DataFrame,
    r_matrix: np.ndarray,
    out: str,
    study: str,
    row: pd.Series,
    windowsizekb: int,
    log: Log,
    verbose: bool
) -> str:
    
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

def _merge_ld_map_with_sumstats(
    row: pd.Series,
    locus_sumstats: pd.DataFrame,
    ld_map: pd.DataFrame,
    log: Log = Log(),
    suffixes: Optional[List[str]] = None
) -> pd.DataFrame:
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
    locus_sumstats: pd.DataFrame,
    ld_map: pd.DataFrame,
    log: Log = Log(),
    suffixes: Optional[List[str]] = None,
    verbose: bool = True
) -> pd.DataFrame:
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
