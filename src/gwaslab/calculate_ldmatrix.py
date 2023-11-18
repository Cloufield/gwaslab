import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.Log import Log
from gwaslab.getsig import getsig
from gwaslab.processreference import _process_vcf_and_bfile

def tofinemapping(sumstats, study=None, bfile=None, vcf=None, out="./",windowsizekb=1000,n_cores=2, overwrite=False,log=Log()):
## for each lead variant 
    ## extract snp list from sumstats
    sig_df = getsig(sumstats,id="SNPID",chrom="CHR",pos="POS",p="P")

    # drop duplicate!!!!
    log.write(" -Dropping duplicated SNPIDs...")
    sumstats = sumstats.drop_duplicates(subset=["SNPID"]).copy()


    output_file_list = pd.DataFrame(columns=["SNPID","SNPID_List","LD_r_matrix","Locus_sumstats"])
    
    plink_log=""

    for index, row in sig_df.iterrows():
        # extract snplist in each locus
        gc.collect()

        log.write(" -Processing locus with lead variant {} at CHR {} POS {} #########################...".format(row["SNPID"],row["CHR"],row["POS"]))
        
        is_in_locus = (sumstats["CHR"] == row["CHR"]) & (sumstats["POS"] >= row["POS"] - windowsizekb*1000) & (sumstats["POS"] < row["POS"] + windowsizekb*1000)
        
        locus_sumstats = sumstats.loc[is_in_locus,:].copy()
        
        
        #process reference file
        bfile_prefix, plink_log, ref_bim = _process_vcf_and_bfile(  chrlist=[row["CHR"]],
                                                                    bfile=bfile, 
                                                                    vcf=vcf, 
                                                                    plink_log=plink_log,
                                                                    n_cores=n_cores, 
                                                                    log=log,
                                                                    load_bim=True,
                                                                    overwrite=overwrite)

        ## check available snps with reference file
        matched_sumstats = _align_sumstats_with_bim(row=row, 
                                                    locus_sumstats=locus_sumstats, 
                                                    ref_bim=ref_bim[0],
                                                    log=log)
        
        #########################################################################################################
        # create matched snp list
        matched_snp_list_path,matched_sumstats_path=_export_snplist_and_locus_sumstats(matched_sumstats=matched_sumstats, 
                                                                                       out=out, 
                                                                                       study=study, 
                                                                                       row=row, 
                                                                                       windowsizekb=windowsizekb,
                                                                                       log=log)
        #########################################################################################################

        ## Calculate ld matrix using PLINK
        matched_ld_matrix_path = _calculate_ld_r(study=study,
                                                 matched_sumstats_snpid= matched_sumstats["SNPID"],
                                                row=row, 
                                                bfile_prefix=bfile_prefix, 
                                                n_cores=n_cores, 
                                                windowsizekb=windowsizekb,
                                                out=out,
                                                plink_log=plink_log,
                                                log=log)
    
    
    # print file list
        row_dict={}
        row_dict["SNPID"]=row["SNPID"]
        row_dict["SNPID_List"] = matched_snp_list_path
        row_dict["LD_r_matrix"] = matched_ld_matrix_path
        row_dict["Locus_sumstats"] = matched_sumstats_path
        file_row = pd.Series(row_dict).to_frame().T
        output_file_list = pd.concat([output_file_list, file_row],ignore_index=True)
        output_file_list["study"] = study
        output_file_list_path =  "{}/{}_{}.filelist".format(out.rstrip("/"), study,windowsizekb)
        output_file_list.to_csv(output_file_list_path,index=None,sep="\t")
        log.write(" -File list is saved to: {}".format(output_file_list_path))

    log.write(" -Finished LD matrix calculation.")
    return output_file_list_path



def _calculate_ld_r(study, matched_sumstats_snpid, row, bfile_prefix, n_cores, windowsizekb,out,plink_log,log):
    log.write(" -#Calculating LD r...")
    
    if "@" in bfile_prefix:
        bfile_to_use = bfile_prefix.replace("@",str(row["CHR"]))
    else:
        bfile_to_use = bfile_prefix
    
    if os.path.exists(bfile_to_use+".bed"):
        snplist_path =   "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
        output_prefix =  "{}/{}_{}_{}".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
        script_vcf_to_bfile = """
        plink \
            --bfile {} \
            --keep-allele-order \
            --extract {} \
            --chr {} \
            --r square gz \
            --allow-no-sex \
            --threads {} \
            --write-snplist \
            --out {}
        """.format(bfile_to_use, snplist_path , row["CHR"], n_cores, output_prefix)
        
        try:
            output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
            plink_log+=output + "\n"
            log.write(" -Finished calculating LD r for locus with lead variant {} at CHR {} POS {} ################...".format(row["SNPID"],row["CHR"],row["POS"]))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
        
        _check_snpid_order(snplist_path.replace(".raw",""), matched_sumstats_snpid,log)

        return output_prefix+".ld.gz"

def _align_sumstats_with_bim(row, locus_sumstats, ref_bim, log=Log()):
    
    log.write(" -#variants in locus ({}): {}".format(row["SNPID"],len(locus_sumstats)))
    # convert category to string
    locus_sumstats["EA"] = locus_sumstats["EA"].astype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].astype("string")

    # matching by SNPID
    # preserve bim keys (use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.)
    combined_df = pd.merge(ref_bim, locus_sumstats, on="SNPID",how="inner")
    
    # match allele
    allele_match =  ((combined_df["EA"] == combined_df["EA_bim"]) & (combined_df["NEA"] == combined_df["NEA_bim"]) ) | ((combined_df["EA"] == combined_df["NEA_bim"])& (combined_df["NEA"] == combined_df["EA_bim"]))
    log.write(" -#Variants with matched alleles:{}".format(sum(allele_match)))

    # fliipped allele
    ea_mis_match = combined_df["EA"] != combined_df["EA_bim"]
    log.write(" -#Variants with flipped alleles:{}".format(sum(ea_mis_match)))
    
    # adjust statistics
    output_columns=["SNPID","CHR","POS","EA_bim","NEA_bim"]

    if ("BETA" in locus_sumstats.columns) and ("SE" in locus_sumstats.columns):
        combined_df.loc[ea_mis_match,"BETA"] = - combined_df.loc[ea_mis_match,"BETA"]
        output_columns.append("BETA")
        output_columns.append("SE")
    if "Z" in locus_sumstats.columns:
        combined_df.loc[ea_mis_match,"Z"] = - combined_df.loc[ea_mis_match,"Z"]
        output_columns.append("Z")
    if "EAF" in locus_sumstats.columns:
        combined_df.loc[ea_mis_match,"EAF"] = 1 - combined_df.loc[ea_mis_match,"EAF"]
        output_columns.append("EAF")
    if "N" in locus_sumstats.columns:
        output_columns.append("N")
    
    return combined_df.loc[allele_match,output_columns]


def _export_snplist_and_locus_sumstats(matched_sumstats, out, study, row, windowsizekb,log):
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        matched_sumstats["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)

        # create locus-sumstats EA, NEA, (BETA, SE), Z 
        matched_sumstats_path =  "{}/{}_{}_{}.sumstats.gz".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        
        to_export_columns=["CHR","POS","EA_bim","NEA_bim"]
        if "Z" in matched_sumstats.columns :
            to_export_columns.append("Z")
        if ("BETA" in matched_sumstats.columns) and ("SE" in matched_sumstats.columns):
            to_export_columns.append("BETA")
            to_export_columns.append("SE")
        if "EAF" in matched_sumstats.columns :
            to_export_columns.append("EAF")
        if "N" in matched_sumstats.columns:
            to_export_columns.append("N")
        matched_sumstats.loc[:, ["SNPID"]+to_export_columns].to_csv(matched_sumstats_path, index=None)
        return matched_snp_list_path, matched_sumstats_path

def _check_snpid_order(snplist_path, matched_sumstats_snpid,log):
    snpid_list = pd.read_csv(snplist_path,dtype="string",header=None)[0]
    if list(matched_sumstats_snpid) == list(snpid_list):
        log.write(" -Order matched.")
    else:
        log.write(" -Warning: Order not matched...")