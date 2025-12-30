from typing import TYPE_CHECKING, Optional, Union, Dict, Any, Tuple, List
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.io.io_plink import _process_plink_input_files, _load_bim_single, _load_pvar_single
from gwaslab.extension import _checking_plink_version

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _calculate_prs(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], 
          bfile: Optional[str] = None, 
          vcf: Optional[str] = None, 
          study: Optional[str] = None, 
          out: str = "./",
          id_to_use: Optional[str] = None,
          threads: int = 1, 
          memory: Optional[int] = None, 
          overwrite: bool = False,
          mode: Optional[str] = None,
          delete: bool = True,
          plink: str = "plink",
          plink2: str = "plink2",
          log: Log = Log(),
          **kwargs: Any) -> None:
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    #matching_alleles
        #read_bim
        #match_id
        #merge
        #output
        #calculate using plink2
        chrlist = list(sumstats["CHR"].unique())
        chrlist.sort()
        plink_log = ""
        #process reference fileWWW
        bfile_prefix, plink_log, ref_bim, filetype = _process_plink_input_files(  
                                                                    chrlist=chrlist,
                                                                    bfile=bfile, 
                                                                    vcf=vcf, 
                                                                    plink_log=plink_log,
                                                                    threads=threads, 
                                                                    log=log,
                                                                    load_bim=False,
                                                                    overwrite=overwrite,
                                                                    plink=plink,
                                                                    plink2=plink2,
                                                                    **kwargs)
        score_file_path_list =[]
        for index, chrom in enumerate(chrlist): 
            chr_sumstats = sumstats.loc[sumstats["CHR"]==chrom,:].copy()
            
            chr_sumstats_matched = _match_snpid_with_bim(chrom = chrom, 
                                                              chr_sumstats = chr_sumstats, 
                                                              bfile_prefix= bfile_prefix, 
                                                              id_to_use=id_to_use,
                                                              log=log,filetype=filetype)
            
            model_path =   "{}/{}_{}.model".format(out.rstrip("/"),study, chrom)
            log.write("  -Model file is saved to : {} .".format(model_path))
            chr_sumstats_matched.to_csv(model_path,sep="\t",index=None)

            score_file_path, plink_log = _run_calculate_prs(study=study, 
                               chrom=chrom , 
                               model_path=model_path,
                               bfile_prefix=bfile_prefix, 
                               threads=threads, 
                               out=out, 
                               plink_log=plink_log, 
                               log=log, 
                               memory=memory, 
                               mode=mode,filetype=filetype,plink2=plink2)
            score_file_path_list.append(score_file_path)
            if delete == True:
                os.remove(model_path)
        combined_results_summary = _combine_all_chr_prs(score_file_path_list)
        return combined_results_summary




def _run_calculate_prs(
    study: str,
    chrom: int,
    model_path: str,
    bfile_prefix: str,
    threads: int,
    out: str,
    plink_log: str,
    log: Log,
    memory: Optional[int],
    filetype: str,
    plink2: str,
    mode: Optional[str] = None
) -> Tuple[str, str]:
    
    log.write(" -Start to calculate PRS for Chr {}...".format(chrom))
    _checking_plink_version(plink2=plink2, log=log)
    
    if "@" in bfile_prefix:
        bpfile_to_use = bfile_prefix.replace("@",str(chrom))
    else:
        bpfile_to_use = bfile_prefix
    
    if filetype=="bfile":
        file_flag = "--bfile {}".format(bpfile_to_use) 
    else:
        file_flag = "--pfile {}".format(bpfile_to_use) 

    output_prefix =  model_path.replace(".model","")
    
    if memory is not None:
        memory_flag = "--memory {}".format(memory)
    
    script_vcf_to_bfile = """
    {} \
        {} \
        --score {} 1 2 3 header {} cols=+scoresums,+denom ignore-dup-ids \
        --chr {} \
        --threads {} {}\
        --out {}
    """.format(plink2, file_flag, model_path ,  mode if mode is not None else "", chrom, threads, memory_flag if memory is not None else "", output_prefix)

    try:
        output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
        plink_log+=output + "\n"
        log.write(" -Finished calculating PRS for Chr {} : {}...".format(chrom, output_prefix+".sscore"))
    except subprocess.CalledProcessError as e:
        log.write(e.output)
    gc.collect()
    return output_prefix+".sscore",plink_log




def _match_snpid_with_bim(
    chrom: int,
    chr_sumstats: pd.DataFrame,
    bfile_prefix: str,
    filetype: str,
    id_to_use: Optional[str] = None,
    log: Log = Log()
) -> pd.DataFrame:
    
    #load bim 
    log.write("   -#variants in model for Chr {} : {}".format(chrom, len(chr_sumstats))) 
    if filetype == "bfile":
        single_bim = _load_bim_single(chrom, bfile_prefix, log)
    else:
        single_bim = _load_pvar_single(chrom, bfile_prefix, log)
    if id_to_use is None:
        
        # create temp ID
        
        chr_sumstats["_TEMP_ID"]=chr_sumstats["POS"].astype("Int64")
        
        single_bim["_TEMP_ID"]=single_bim["POS_bim"].astype("Int64")

        # merge by pos
        chr_sumstats_bim = pd.merge(chr_sumstats, single_bim, on="_TEMP_ID", how="inner")
        chr_sumstats_bim["EA"] = chr_sumstats_bim["EA"].astype("string")
        chr_sumstats_bim["NEA"] = chr_sumstats_bim["NEA"].astype("string")
        
        # effect allele in NEA_bim or EA_bim
        if "NEA" in chr_sumstats.columns:
            # match both EA and NEA
            is_allele_perfect_match = (chr_sumstats_bim["EA"] == chr_sumstats_bim["EA_bim"]) & (chr_sumstats_bim["NEA"] == chr_sumstats_bim["NEA_bim"])
            is_allele_flipped_match = (chr_sumstats_bim["EA"] == chr_sumstats_bim["NEA_bim"])& (chr_sumstats_bim["NEA"] == chr_sumstats_bim["EA_bim"])
            is_allele_matched = is_allele_flipped_match | is_allele_perfect_match
        else:
            # match only EA
            is_allele_rough_match = (chr_sumstats_bim["EA"] == chr_sumstats_bim["EA_bim"]) & (chr_sumstats_bim["EA"] == chr_sumstats_bim["NEA_bim"])
            is_allele_matched = is_allele_rough_match   
    else:
        chr_sumstats_bim = pd.merge(chr_sumstats_bim, single_bim, left_on=id_to_use, right_on="SNPID_bim", how="inner")
    
    matched_cols = ["SNPID_bim", "EA", "BETA"]
    
    matched_sumstats = chr_sumstats_bim.loc[is_allele_matched, matched_cols]
    log.write("   -#variants matched for Chr {} : {}".format(chrom, len(matched_sumstats))) 
    return matched_sumstats





def _combine_all_chr_prs(score_file_path_list: List[str]) -> pd.DataFrame:
    combined_results = pd.DataFrame()
    for index, score_file_path in enumerate(score_file_path_list):
        score_file = pd.read_csv(score_file_path,sep="\t")
        combined_results = pd.concat([combined_results, score_file],ignore_index=True)
        gc.collect()
    
    grouby_columns=[]
    if "#FID" in combined_results.columns:
        grouby_columns.append("#FID")
    if "IID" in combined_results.columns:
        grouby_columns.append("IID")
    if "#IID" in combined_results.columns:
        grouby_columns.append("#IID")

    combined_results_summary = combined_results.groupby(grouby_columns)[["ALLELE_CT","DENOM","NAMED_ALLELE_DOSAGE_SUM","SCORE1_SUM"]].sum()
    combined_results_summary["SCORE1_AVG_recal"] = combined_results_summary["SCORE1_SUM"] / combined_results_summary["DENOM"] 
    return combined_results_summary.reset_index()
