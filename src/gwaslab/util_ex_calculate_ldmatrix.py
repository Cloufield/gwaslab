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

def tofinemapping(sumstats, 
                  study=None, 
                  bfile=None, 
                  vcf=None, 
                  loci=None,
                  out="./",
                  plink="plink",
                  plink2="plink2",
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
                                                                    n_cores=n_cores, 
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
                                                            n_cores=n_cores, 
                                                            windowsizekb=windowsizekb,
                                                            out=out,
                                                            plink_log=plink_log,
                                                            log=log,
                                                            filetype=filetype,
                                                            plink=plink,
                                                            plink2=plink2,
                                                            verbose=verbose)
    
    
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



def _calculate_ld_r(study, matched_sumstats_snpid, row, bfile_prefix, n_cores, windowsizekb,out,plink_log,log,memory,mode,filetype,plink,plink2,verbose=True):
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
        
        script_vcf_to_bfile = """
        {} \
            --bfile {} \
            --keep-allele-order \
            --extract {} \
            --chr {} \
            --{} square gz \
            --allow-no-sex \
            --threads {} {}\
            --write-snplist \
            --out {}
        """.format(plink, bfile_to_use, snplist_path , row["CHR"], mode, n_cores, memory_flag if memory is not None else "", output_prefix)

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

def _align_sumstats_with_bim(row, locus_sumstats, ref_bim, log=Log(),suffixes=None):
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
    output_columns=["SNPID","CHR","POS","EA_bim","NEA_bim"]
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


def _export_snplist_and_locus_sumstats(matched_sumstats, out, study, row, windowsizekb,log,suffixes=None):
        if suffixes is None:
            suffixes=[""]
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        
        matched_sumstats["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)
        log.write(" -Exporting SNP list of {}  to: {}...".format(len(matched_sumstats) ,matched_snp_list_path))

        # create locus-sumstats EA, NEA, (BETA, SE), Z 
        matched_sumstats_path =  "{}/{}_{}_{}.sumstats.gz".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        
        to_export_columns=["CHR","POS","EA_bim","NEA_bim"]
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
        matched_sumstats[ ["SNPID"]+to_export_columns].to_csv(matched_sumstats_path, index=None)
        return matched_snp_list_path, matched_sumstats_path

def _check_snpid_order(snplist_path, matched_sumstats_snpid,log):
    snpid_list = pd.read_csv(snplist_path,dtype="string",header=None)[0]
    if list(matched_sumstats_snpid) == list(snpid_list):
        log.write(" -Sumstats SNPID order and LD matrix SNPID order are matched.")
    else:
        log.warning("Sumstats SNPID order and LD matrix SNPID order are not matched!")

def _extract_variants_in_locus(sumstats, windowsizekb, locus, chrom = "CHR", pos="POS"):
    
    is_in_locus = (sumstats["CHR"] == locus[0]) & (sumstats["POS"] >= locus[1] - windowsizekb*1000) & (sumstats["POS"] < locus[1] + windowsizekb*1000)
    ## extract snp list from sumstats
    locus_sumstats = sumstats.loc[is_in_locus,:].copy()
    return locus_sumstats