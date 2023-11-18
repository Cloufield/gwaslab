import subprocess
import numpy as np
import os
import pandas as pd
from gwaslab.Log import Log
from gwaslab.processreference import _process_vcf_and_bfile
from gwaslab.version import _checking_plink_version

def _clump(insumstats, vcf=None, scaled=False, out="clumping_plink2", overwrite=False, study=None, bfile=None, n_cores=2, 
          chrom=None, clump_p1=5e-8, clump_p2=5e-8, clump_r2=0.2, clump_kb=250,log=Log()):
    ## process reference
    log.write("Start to perform clumping...")
    log.write(" -Clumping parameters for PLINK2:")
    log.write("  -clump_p1 : {}...".format(clump_p1))
    log.write("  -clump_p2 : {}...".format(clump_p2))
    log.write("  -clump_kb : {}...".format(clump_kb))
    log.write("  -clump_r2 : {}...".format(clump_r2))
    if scaled == True:
        log.write(" -Clumping will be performed using MLOG10P")
        clump_log10_p1=-np.log10(clump_p1)
        clump_log10_p2=-np.log10(clump_p2)
        log.write("  -clump_log10_p1 : {}...".format(clump_log10_p1))
        log.write("  -clump_log10_p2 : {}...".format(clump_log10_p2))
        sumstats = insumstats.loc[insumstats["MLOG10P"]>min(clump_log10_p1,clump_log10_p2),:].copy()
    # extract lead variants
    else:
        log.write(" -Clumping will be performed using P")
        sumstats = insumstats.loc[insumstats["P"]<max(clump_p1,clump_p2),:].copy()
    log.write(" -Significant variants on CHR: ",list(sumstats["CHR"].unique()))
    
    plink_log=""

    # process reference file
    bfile, plink_log, ref_bim = _process_vcf_and_bfile(chrlist=sumstats["CHR"].unique(), 
                                                       bfile=bfile, 
                                                       vcf=vcf, 
                                                       n_cores=n_cores,
                                                       plink_log=plink_log, 
                                                       log=log,
                                                       overwrite=overwrite)           
    
    ## process sumstats by CHR
    for i in sumstats["CHR"].unique():
        log.write(" -Processing sumstats for CHR {}...".format(i))
        
        if "@" in bfile:
            bfile_to_use = bfile.replace("@",str(i))
        else:
            bfile_to_use = bfile
        
        # checking # variants
        try:
            bim = pd.read_csv(bfile_to_use + ".bim",usecols=[1],header=None,sep="\s+")[1]
            snplist = sumstats.loc[sumstats["CHR"]==i,"SNPID"]
            is_on_both = snplist.isin(bim)
            log.write(" -#variants in reference file: {}...".format(len(bim)))
            log.write(" -#variants in sumstats: {}...".format(len(snplist)))
            log.write(" -#variants available in both reference and sumstats: {}...".format(sum(is_on_both)))

            if scaled == True:
                sumstats.loc[(sumstats["CHR"]==i) &(is_on_both),["SNPID","MLOG10P"]].to_csv("_gwaslab_tmp.{}.SNPIDP".format(i),index=None,sep="\t")
            else:
                sumstats.loc[(sumstats["CHR"]==i) & (is_on_both),["SNPID","P"]].to_csv("_gwaslab_tmp.{}.SNPIDP".format(i),index=None,sep="\t")
        except:
            log.write(" -Not available for: {}...".format(i))
        
    # create a empty dataframe for combining results from each CHR 
    results = pd.DataFrame()

    
    # clumping using plink
    for i in sumstats["CHR"].unique():
        chrom = i
        # temp file  
        clump = "_gwaslab_tmp.{}.SNPIDP".format(chrom)
        # output prefix
        out_single_chr= out + ".{}".format(chrom)
        
        if "@" in bfile:
            bfile_to_use = bfile.replace("@",str(i))
        else:
            bfile_to_use = bfile

        log.write(" -Performing clumping for CHR {}...".format(i))
        log = _checking_plink_version(v=2, log=log)
        if scaled == True:
            # clumping using LOG10P
            script = """
            plink2 \
                --bfile {}\
                --chr {} \
                --clump {} \
                --clump-log10 \
                --clump-field MLOG10P \
                --clump-snp-field SNPID \
                --clump-log10-p1 {} \
                --clump-log10-p2 {} \
                --clump-r2 {} \
                --clump-kb {} \
                --threads {} \
                --out {}
            """.format(bfile_to_use, chrom, clump, clump_log10_p1, clump_log10_p2, clump_r2, clump_kb, n_cores, out_single_chr)    
        else:
            # clumping using P
            script = """
            plink2 \
                --bfile {}\
                --chr {} \
                --clump {} \
                --clump-field P \
                --clump-snp-field SNPID \
                --clump-p1 {} \
                --clump-p2 {} \
                --clump-r2 {} \
                --clump-kb {} \
                --threads {} \
                --out {}
            """.format(bfile_to_use, chrom, clump, clump_p1, clump_p2, clump_r2, clump_kb, n_cores, out_single_chr)
        
        try:
            output = subprocess.check_output(script, stderr=subprocess.STDOUT, shell=True,text=True)
            log.write(" -Saved results for CHR {} to : {}".format(i,"{}.clumps".format(out_single_chr)))
            plink_log +=output + "\n"
        except subprocess.CalledProcessError as e:
            log.write(e.output)
        #os.system(script)
        
        clumped = pd.read_csv("{}.clumps".format(out_single_chr),usecols=[2,0,1,3],sep="\s+")
        results = pd.concat([results,clumped],ignore_index=True)
        
        # remove temp SNPIDP file 
        os.remove(clump)
    
    results = results.sort_values(by=["#CHROM","POS"]).rename(columns={"#CHROM":"CHR","ID":"SNPID"})
    log.write("Finished clumping.")
    results_sumstats = insumstats.loc[insumstats["SNPID"].isin(results["SNPID"]),:].copy()
    return results_sumstats, plink_log


       