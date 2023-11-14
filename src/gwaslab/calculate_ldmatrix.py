import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.Log import Log
from gwaslab.getsig import getsig

def tofinemapping(sumstats, bfile=None, vcf=None, out="./",windowsizekb=1000,n_cores=2, log=Log()):
## for each lead variant 
    ## extract snp list from sumstats
    sig_df = getsig(sumstats,id="SNPID",chrom="CHR",pos="POS",p="P")
    plink_log=""

    # drop duplicate!!!!

    for index, row in sig_df.iterrows():
        # extract snplist in each locus
        gc.collect()
        log.write(" -Processing locus with lead variant {} at CHR {} POS {} #########################...".format(row["SNPID"],row["CHR"],row["POS"]))
        
        is_in_locus = (sumstats["CHR"] == row["CHR"]) & (sumstats["POS"] >= row["POS"] - windowsizekb*1000) & (sumstats["POS"] < row["POS"] + windowsizekb*1000)
        
        locus_sumstats = sumstats.loc[is_in_locus,:].copy()
        
        log.write(" -#variants in locus ({}): {}".format(row["SNPID"],len(locus_sumstats)))
        #process reference file
        if bfile is None:
            if vcf is not None:
                log.write(" -Processing VCF : {}...".format(vcf))
                i = row["CHR"]
                bfile_gwaslab = vcf.replace(".vcf.gz","") + ".{}".format(i)
                bfile_prefix= vcf.split("/")[-1].replace(".vcf.gz","") + ".{}".format(i)
                log.write("  -Processing VCF for CHR {}...".format(i))
                
                if not os.path.exists(bfile_gwaslab+".bed"):
                    script_vcf_to_bfile = """
                    plink2 \
                        --vcf {} \
                        --chr {} \
                        --make-bed \
                        --rm-dup force-first \
                        --threads {}\
                        --out {}
                    """.format(vcf, i, n_cores, bfile_gwaslab)
                    
                    try:
                        output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                        plink_log+=output + "\n"
                    except subprocess.CalledProcessError as e:
                        log.write(e.output)
                        
                else:
                    log.write("  -Plink bfile for CHR {} exists. Skipping...".format(i))
                bim_path =bfile_gwaslab+".bim"
                ref_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
                log.write("#variants in ref file: {}".format(len(ref_bim)))
            else:
                log.write("  -Please provide PLINK bfile or VCF as reference!")
        
        else:
            log.write(" -PLINK bfile as LD reference panel: {}".format(bfile))  
            bim_path =bfile_gwaslab+".bim"
            ref_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
            log.write(" -#variants in ref file: {}".format(len(ref_bim)))

        align_sumstats_with_bim(locus_sumstats, ref_bim):
        #avaiable_snplist = locus_snplist[locus_snplist.isin(ref_snplist[1])]
        log.write(" -#variants available in sumstats and LD panel: {}".format(len(avaiable_snplist)))
        
        avaiable_snplist.to_csv("{}/available_{}_{}.snplist".format(out.rstrip("/"),row["SNPID"],windowsizekb),index=None,header=None)

        log.write(" -#Calculating LD r...")
        if os.path.exists(bfile_gwaslab+".bed"):
            script_vcf_to_bfile = """
            plink \
                --bfile {} \
                --extract {} \
                --chr {} \
                --r square gz \
                --allow-no-sex \
                --threads {} \
                --out {}
            """.format(bfile_gwaslab, "{}/available_{}_{}.snplist".format(out.rstrip("/"),row["SNPID"],windowsizekb), i, n_cores, "{}/{}_{}".format(out.rstrip("/"),row["SNPID"],windowsizekb))
            
            try:
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
                log.write(" -Finished calculating LD r for locus with lead variant {} at CHR {} POS {} ################...".format(row["SNPID"],row["CHR"],row["POS"]))
            except subprocess.CalledProcessError as e:
                log.write(e.output)

                

    ## check available snps with reference file

    ## matching alleles

    ## extract available SNPs and calculate ld matrix using PLINK
def calculateldr():
    pass

def align_sumstats_with_bim(locus_sumstats, bim_path, log=Log()):
    
    locus_sumstats["EA"] = locus_sumstats["EA"].atype("string")
    locus_sumstats["NEA"] = locus_sumstats["NEA"].atype("string")

    
    
    combined_df = pd.merge(locus_sumstats, bim, on="SNPID",how="inner")
    
    allele_match =  ((combined_df["EA"] == combined_df["EA_bim"]) & (combined_df["NEA"] == combined_df["NEA_bim"]) ) | ((combined_df["EA"] == combined_df["NEA_bim"])& (combined_df["NEA"] == combined_df["EA_bim"]))
    log.write("#Variants with matched alleles:{}".format(sum(allele_match)))
    ea_mis_match = combined_df["EA"] != combined_df["EA_bim"]
    log.write("#Variants with flipped alleles:{}".format(sum(allele_match)))
    
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
#plink \
#  --bfile ${plinkFile} \
#  --keep-allele-order \
#  --r square \
#  --extract sig_locus.snplist \
#  --out sig_locus_mt
#
#plink \
#  --bfile ${plinkFile} \
#  --keep-allele-order \
#  --r2 square \
#  --extract sig_locus.snplist \
#  --out sig_locus_mt_r2
