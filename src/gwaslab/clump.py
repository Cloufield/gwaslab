import subprocess
import numpy as np
import os
import pandas as pd
from gwaslab.Log import Log

def _clump(insumstats, vcf=None, scaled=False, out="clumping_plink2", overwrite=False, bfile=None, n_cores=2, 
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

    chrlist = sumstats["CHR"].unique()
    
    # process reference file
    bfile,plink_log = _process_reference(chrlist, bfile, vcf, n_cores,plink_log, log,overwrite)           
    
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


def _process_reference(chrlist, bfile, vcf, n_cores, plink_log, log, overwrite):
    
    # return bfile_path and plink_log

    if bfile is None:
        if vcf is not None:
            log.write(" -Processing VCF : {}...".format(vcf))
            
            for i in chrlist:
                log.write("  -Processing VCF for CHR {}...".format(i))
                if "@" in vcf:
                    vcf_to_load = vcf.replace("@",str(i))
                else:
                    vcf_to_load = vcf

                bfile_root = vcf.replace(".vcf.gz","")
                bfile_prefix = bfile_root + ".{}".format(i)
                
                if (not os.path.exists(bfile_prefix+".bed")) or overwrite:
                    script_vcf_to_bfile = """
                    plink2 \
                        --vcf {} \
                        --chr {} \
                        --make-bed \
                        --rm-dup force-first \
                        --threads {}\
                        --out {}
                    """.format(vcf_to_load, i, n_cores, bfile_prefix)
                    
                    try:
                        log.write("  -Converting VCF to bed: {}.bim/bed/fam...".format(bfile_prefix))
                        output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                        plink_log+=output + "\n"
                    except subprocess.CalledProcessError as e:
                        log.write(e.output)
                        
                else:
                    log.write("  -Plink bfile for CHR {} exists. Skipping...".format(i))
                
            return bfile_root+".@", plink_log
        else:
            log.write("  -Please provide PLINK bfile or VCF as reference!")
    else:
        log.write(" -PLINK bfile as LD reference panel: {}".format(bfile))      
        return bfile, plink_log          