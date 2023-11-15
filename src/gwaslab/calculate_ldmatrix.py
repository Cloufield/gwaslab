import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.Log import Log
from gwaslab.getsig import getsig

def tofinemapping(sumstats, study="sumstats1", bfile=None, vcf=None, out="./",windowsizekb=1000,n_cores=2, log=Log()):
## for each lead variant 
    ## extract snp list from sumstats
    sig_df = getsig(sumstats,id="SNPID",chrom="CHR",pos="POS",p="P")

    # drop duplicate!!!!
    log.write(" -Dropping duplicated SNPIDs...")
    sumstats = sumstats.drop_duplicates(subset=["SNPID"]).copy()


    output_file_list = pd.DataFrame(columns=["SNPID","SNPID_List","LD_r_matrix"])
    
    for index, row in sig_df.iterrows():
        # extract snplist in each locus
        gc.collect()

        log.write(" -Processing locus with lead variant {} at CHR {} POS {} #########################...".format(row["SNPID"],row["CHR"],row["POS"]))
        
        is_in_locus = (sumstats["CHR"] == row["CHR"]) & (sumstats["POS"] >= row["POS"] - windowsizekb*1000) & (sumstats["POS"] < row["POS"] + windowsizekb*1000)
        
        locus_sumstats = sumstats.loc[is_in_locus,:].copy()
        
        
        #process reference file
        ref_bim , plink_log, bfile_gwaslab = _process_vcf_and_bfile(
                                                                    bfile=bfile, 
                                                                    vcf=vcf, 
                                                                    row=row, 
                                                                    n_cores=n_cores, 
                                                                    log=log)

        ## check available snps with reference file
        matched_sumstats = _align_sumstats_with_bim(row=row, 
                                                    locus_sumstats=locus_sumstats, 
                                                    ref_bim=ref_bim,
                                                    log=log)
        
        #avaiable_snplist = locus_snplist[locus_snplist.isin(ref_snplist[1])]
        log.write(" -#variants available in sumstats and LD panel: {}".format(len(matched_sumstats)))
        
        # create matched snp list
        matched_snp_list_path = "{}/{}_{}_{}.snplist".format(out.rstrip("/"), study, row["SNPID"] ,windowsizekb)
        matched_sumstats["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)

        ## Calculate ld matrix using PLINK
        matched_ld_matrix_path = _calculate_ld_r(study=study,
                                                row=row, 
                                                bfile_gwaslab=bfile_gwaslab, 
                                                n_cores=n_cores, 
                                                windowsizekb=windowsizekb,
                                                out=out,
                                                plink_log=plink_log,
                                                log=log)
    
    
    # print file list
        row_dict={}
        row_dict["SNPID"]=row["SNPID"]
        row_dict["SNPID_List"]=matched_snp_list_path
        row_dict["LD_r_matrix"]=matched_ld_matrix_path
        file_row = pd.Series(row_dict).to_frame().T
        output_file_list = pd.concat([output_file_list, file_row],ignore_index=True)
        output_file_list_path =  "{}/{}_{}.filelist".format(out.rstrip("/"), study,windowsizekb)
        output_file_list.to_csv(output_file_list_path,index=None,sep="\t")
        log.write(" -File list is saved to: {}".format(output_file_list_path))

    log.write(" -Finished LD matrix calculation.")


def _process_vcf_and_bfile(bfile, vcf, row, n_cores, log):
    plink_log =""
    if bfile is None:
        if vcf is not None:
            log.write(" -Processing VCF : {}...".format(vcf))
            bfile_gwaslab = vcf.replace(".vcf.gz","") + ".{}".format(row["CHR"])
            bfile_prefix= vcf.split("/")[-1].replace(".vcf.gz","") + ".{}".format(row["CHR"])
            log.write("  -Processing VCF for CHR {}...".format(row["CHR"]))
            
            if not os.path.exists(bfile_gwaslab+".bed"):
                script_vcf_to_bfile = """
                plink2 \
                    --vcf {} \
                    --chr {} \
                    --make-bed \
                    --rm-dup force-first \
                    --threads {}\
                    --out {}
                """.format(vcf, row["CHR"], n_cores, bfile_gwaslab)
                
                try:
                    output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                    plink_log+=output + "\n"
                except subprocess.CalledProcessError as e:
                    log.write(e.output)
                    
            else:
                log.write("  -Plink bfile for CHR {} exists. Skipping...".format(row["CHR"]))
            ## load bim file
            bim_path =bfile_gwaslab+".bim"
            ref_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
            log.write("#variants in ref file: {}".format(len(ref_bim)))
        else:
            log.write("  -Please provide PLINK bfile or VCF as reference!")
    else:
        log.write(" -PLINK bfile as LD reference panel: {}".format(bfile))  
        ## load bim file
        bim_path =bfile_gwaslab+".bim"
        ref_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
        log.write(" -#variants in ref file: {}".format(len(ref_bim)))
    
    return ref_bim , plink_log, bfile_gwaslab

def _calculate_ld_r(study, row, bfile_gwaslab, n_cores, windowsizekb,out,plink_log,log):
    log.write(" -#Calculating LD r...")
    if os.path.exists(bfile_gwaslab+".bed"):
        snplist_path =   "{}/{}_{}_{}.snplist".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
        output_prefix =  "{}/{}_{}_{}".format(out.rstrip("/"),study,row["SNPID"],windowsizekb)
        script_vcf_to_bfile = """
        plink \
            --bfile {} \
            --extract {} \
            --chr {} \
            --r square gz \
            --allow-no-sex \
            --threads {} \
            --out {}
        """.format(bfile_gwaslab, snplist_path , row["CHR"], n_cores, output_prefix)
        
        try:
            output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
            plink_log+=output + "\n"
            log.write(" -Finished calculating LD r for locus with lead variant {} at CHR {} POS {} ################...".format(row["SNPID"],row["CHR"],row["POS"]))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
        return output_prefix+".ld.gz"

def _align_sumstats_with_bim(row, locus_sumstats, ref_bim, log=Log()):
    
    log.write(" -#variants in locus ({}): {}".format(row["SNPID"],len(locus_sumstats)))
    # convert category to string
    locus_sumstats["EA"] = locus_sumstats["EA"].astype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].astype("string")

    # matching by SNPID
    combined_df = pd.merge(locus_sumstats, ref_bim, on="SNPID",how="inner")
    
    # match allele
    allele_match =  ((combined_df["EA"] == combined_df["EA_bim"]) & (combined_df["NEA"] == combined_df["NEA_bim"]) ) | ((combined_df["EA"] == combined_df["NEA_bim"])& (combined_df["NEA"] == combined_df["EA_bim"]))
    log.write(" -#Variants with matched alleles:{}".format(sum(allele_match)))

    # fliipped allele
    ea_mis_match = combined_df["EA"] != combined_df["EA_bim"]
    log.write(" -#Variants with flipped alleles:{}".format(sum(ea_mis_match)))
    
    # adjust statistics
    output_columns=["SNPID","EA_bim","NEA_bim"]

    if "BETA" in locus_sumstats.columns:
        combined_df.loc[ea_mis_match,"BETA"] = - combined_df.loc[ea_mis_match,"BETA"]
        output_columns.append("BETA")
    if "Z" in locus_sumstats.columns:
        combined_df.loc[ea_mis_match,"Z"] = - combined_df.loc[ea_mis_match,"Z"]
        output_columns.append("Z")
    if "EAF" in locus_sumstats.columns:
        combined_df.loc[ea_mis_match,"EAF"] = 1 - combined_df.loc[ea_mis_match,"EAF"]
        output_columns.append("EAF")
    
    return combined_df.loc[allele_match,output_columns]