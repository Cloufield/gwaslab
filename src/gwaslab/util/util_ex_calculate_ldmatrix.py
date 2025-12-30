import subprocess
import os
import gc
import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Optional, List, Dict, Any, Tuple, Union

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.io.io_plink import _process_plink_input_files
from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.extension import _checking_plink_version

from gwaslab.qc.qc_decorator import with_logging
@with_logging(
        start_to_msg="calculate LD matrix",
        finished_msg="calculating LD matrix",
        start_cols=["SNPID","CHR","POS","EA","NEA"],
        start_function="calculate_ld_matrix"
)
def _to_finemapping(
    gls: 'Sumstats', 
    study: Optional[str] = None, 
    bfile: Optional[str] = None, 
    vcf: Optional[str] = None, 
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
    extra_plink_option: str = "",
    verbose: bool = True,
    **kwargs: Any
) -> Tuple[str, pd.DataFrame, str]:
    
    ##start function with col checking##########################################################
    
    sumstats = gls.data
    gls.offload()

    ############################################################################################
    if suffixes is None:
        suffixes=[""]
    if getlead_kwargs is None:
        getlead_kwargs={"windowsizekb":1000}
    
    if loci_chrpos is None:
        if loci is None:
            log.write(" -Loci were not provided. All significant loci will be automatically extracted...",verbose=verbose)
            sig_df = _get_sig(sumstats,variant_id="SNPID",chrom="CHR",pos="POS",p="P"+suffixes[0],**getlead_kwargs)
        else:
            sig_df = sumstats.loc[sumstats["SNPID"].isin(loci),:]
    else:
        sig_df = pd.DataFrame()
        for chrpos in loci_chrpos:
            chrpos_row_dict={}
            chrpos_row_dict["SNPID"]="{}:{}".format(chrpos[0], chrpos[1])
            chrpos_row_dict["CHR"] = chrpos[0]
            chrpos_row_dict["POS"] = chrpos[1]
            chrpos_row = pd.Series(chrpos_row_dict).to_frame().T
            sig_df = pd.concat([sig_df, chrpos_row],ignore_index=True)        
    
    log.write(" -plink1.9 path: {}".format(plink),verbose=verbose)
    log.write(" -plink2 path: {}".format(plink2),verbose=verbose)

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
        
        #process reference file
        bfile_prefix, plink_log, ref_bim, filetype = _process_plink_input_files(  chrlist=[row["CHR"]],
                                                                    bfile=bfile, 
                                                                    vcf=vcf, 
                                                                    plink_log=plink_log,
                                                                    threads=threads, 
                                                                    log=log,
                                                                    load_bim=True,
                                                                    overwrite=overwrite,
                                                                    plink=plink,
                                                                    plink2=plink2,
                                                                    **kwargs)

        ## check available snps with reference file
        matched_sumstats = _align_sumstats_with_bim(row=row, 
                                                    locus_sumstats=locus_sumstats, 
                                                    ref_bim=ref_bim[0],
                                                    log=log,suffixes=suffixes)
        del locus_sumstats
        gc.collect()
        #########################################################################################################
        # create matched snp list
        matched_snp_list_path,matched_sumstats_path=_export_snplist_and_locus_sumstats(matched_sumstats=matched_sumstats, 
                                                                                       out=out, 
                                                                                       study=study, 
                                                                                       row=row, 
                                                                                       windowsizekb=windowsizekb,
                                                                                       log=log,
                                                                                       suffixes=suffixes)
        #########################################################################################################

        ## Calculate ld matrix using PLINK
        matched_ld_matrix_path,plink_log = _calculate_ld_r(study=study,
                                                           mode=mode,
                                                            memory=memory,
                                                            matched_sumstats_snpid= matched_sumstats["SNPID"],
                                                            row=row, 
                                                            bfile_prefix=bfile_prefix, 
                                                            threads=threads, 
                                                            windowsizekb=windowsizekb,
                                                            out=out,
                                                            plink_log=plink_log,
                                                            log=log,
                                                            filetype=filetype,
                                                            plink=plink,
                                                            plink2=plink2,
                                                            extra_plink_option=extra_plink_option,
                                                            ref_allele_path = matched_sumstats_path,
                                                            verbose=verbose)
        del matched_sumstats
        gc.collect()
    
        # print file list
        row_dict={}
        row_dict["SNPID"]=row["SNPID"]
        row_dict["SNPID_LIST"] = matched_snp_list_path
        row_dict["LD_R_MATRIX"] = matched_ld_matrix_path
        row_dict["LOCUS_SUMSTATS"] = matched_sumstats_path+".gz"
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
    
    del sumstats

    gls.reload()

    return output_file_list_path, output_file_list, plink_log



def _calculate_ld_r(
    study: str, 
    matched_sumstats_snpid: pd.Series, 
    row: pd.Series, 
    bfile_prefix: str, 
    threads: int, 
    windowsizekb: int,
    out: str,
    plink_log: str,
    log: Log,
    memory: Optional[int],
    mode: str,
    filetype: str,
    plink: str,
    plink2: str,
    ref_allele_path: str, 
    extra_plink_option: str = "",
    verbose: bool = True
) -> str:
    '''
    Calculate LD r matrix by calling PLINK; return file name and log
    '''
    log.write(" -Start to calculate LD r matrix...",verbose=verbose)
    log = _checking_plink_version(plink=plink, log=log)
    if "@" in bfile_prefix:
        bfile_to_use = bfile_prefix.replace("@",str(row["CHR"]))
    else:
        bfile_to_use = bfile_prefix
    
    if os.path.exists(bfile_to_use+".bed"):
        snplist_path =   "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
        output_prefix =  "{}/{}_{}_{}".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
        
        if memory is not None:
            memory_flag = "--memory {}".format(memory)
        
        if filetype=="pfile":
            raise ValueError("Please use bfile instead of pfile for PLINK1.")
        
        #log.write(" -Flipping plink file ref allele to match...",verbose=verbose)
        #script_vcf_to_bfile = """
        #{} \
        #    --bfile {} \
        #    --extract {} \
        #    --chr {} \
        #    --ref-allele 'force' {} 4 1 \
        #    --threads {} {} \
        #    --make-bed \
        #    --out {}

        #""".format(plink2, bfile_to_use, snplist_path,  row["CHR"],ref_allele_path, threads, memory_flag if memory is not None else "", output_prefix+"_gwaslab_tmp")

        log.write(" -Calculating r matrix...",verbose=verbose)
        script_vcf_to_bfile = """
        {} \
            --bfile {} \
            --a2-allele {} 4 1 \
            --extract {} \
            --chr {} \
            --{} square gz \
            --allow-no-sex \
            --threads {} {}\
            --write-snplist \
            --out {} {}
        """.format(plink, bfile_to_use, ref_allele_path,  snplist_path , row["CHR"], mode, threads, memory_flag if memory is not None else "", output_prefix, extra_plink_option)

        try:
            output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
            #plink_process = subprocess.Popen("exec "+script_vcf_to_bfile, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
            #output1,output2 = plink_process.communicate()
            plink_log+=output + "\n"
            #plink_process.kill()
            log.write(" -Finished calculating LD r for locus with lead variant {} at CHR {} POS {}...".format(row["SNPID"],row["CHR"],row["POS"]))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
        
        _check_snpid_order(snplist_path.replace(".raw",""), matched_sumstats_snpid,log)
        gc.collect()
        return output_prefix+".ld.gz",plink_log

def _align_sumstats_with_bim(
    row: pd.Series, 
    locus_sumstats: pd.DataFrame, 
    ref_bim: pd.DataFrame, 
    log: Log = Log(),
    suffixes: Optional[List[str]] = None
) -> pd.DataFrame:
    '''
    align sumstats with bim
    '''
    if suffixes is None:
            suffixes=[""]
    
    log.write("   -Variants in locus ({}): {}".format(row["SNPID"],len(locus_sumstats)))
    # convert category to string
    locus_sumstats["EA"] = locus_sumstats["EA"].astype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].astype("string")

    # matching by SNPID
    # preserve bim keys (use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.)
    combined_df = pd.merge(ref_bim, locus_sumstats, on="SNPID",how="inner")
    
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
    output_columns=["SNPID","CHR","POS","EA","NEA"]
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
    
    return combined_df.loc[allele_match,output_columns]


def _export_snplist_and_locus_sumstats(
    matched_sumstats: pd.DataFrame, 
    out: str, 
    study: str, 
    row: pd.Series, 
    windowsizekb: int,
    log: Log,
    suffixes: Optional[List[str]] = None
) -> Tuple[str, str]:
        if suffixes is None:
            suffixes=[""]
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        
        matched_sumstats["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)
        log.write(" -Exporting SNP list of {}  to: {}...".format(len(matched_sumstats) ,matched_snp_list_path))

        # create locus-sumstats EA, NEA, (BETA, SE), Z 
        matched_sumstats_path =  "{}/{}_{}_{}.sumstats".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        
        to_export_columns=["CHR","POS","EA","NEA"]
        for suffix in suffixes:
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
        matched_sumstats[ ["SNPID"]+to_export_columns].to_csv(matched_sumstats_path, sep="\t",index=None)
        matched_sumstats[ ["SNPID"]+to_export_columns].to_csv(matched_sumstats_path+".gz", sep="\t",index=None)
        return matched_snp_list_path, matched_sumstats_path

def _check_snpid_order(
    snplist_path: str, 
    matched_sumstats_snpid: pd.Series, 
    log: Log
) -> None:
    snpid_list = pd.read_csv(snplist_path,dtype="string",header=None)[0]
    if list(matched_sumstats_snpid) == list(snpid_list):
        log.write(" -Sumstats SNPID order and LD matrix SNPID order are matched.")
    else:
        log.warning("Sumstats SNPID order and LD matrix SNPID order are not matched!")

def _extract_variants_in_locus(
    sumstats: pd.DataFrame, 
    windowsizekb: int, 
    locus: Tuple[int, int], 
    chrom: str = "CHR", 
    pos: str = "POS"
) -> pd.DataFrame:
    
    is_in_locus = (sumstats["CHR"] == locus[0]) & (sumstats["POS"] >= locus[1] - windowsizekb*1000) & (sumstats["POS"] < locus[1] + windowsizekb*1000)
    ## extract snp list from sumstats
    locus_sumstats = sumstats.loc[is_in_locus,:].copy()
    return locus_sumstats
