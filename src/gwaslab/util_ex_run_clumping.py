import subprocess
import numpy as np
import os
import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.util_ex_process_ref import _process_plink_input_files
from gwaslab.g_version import _checking_plink_version

def _clump(insumstats, vcf=None, scaled=False, out="clumping_plink2", 
           p="P",mlog10p="MLOG10P", overwrite=False, study=None, bfile=None, 
           n_cores=1, memory=None, chrom=None, clump_p1=5e-8, clump_p2=5e-8, clump_r2=0.01, clump_kb=250,
           log=Log(),verbose=True,plink="plink",plink2="plink2"):
    ##start function with col checking##########################################################
    _start_line = "perfrom clumping"
    _end_line = "clumping"
    _start_cols =["SNPID","CHR","POS"]
    _start_function = ".clump()"
    _must_args ={}

    is_enough_info = start_to(sumstats=insumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: raise ValueError("Not enough columns for clumping")
    ############################################################################################
    ## process reference
    log.write("Start to perform clumping...",verbose=verbose)
    log.write(" -Clumping parameters for PLINK2:",verbose=verbose)
    log.write("  -clump_p1 : {}...".format(clump_p1),verbose=verbose)
    log.write("  -clump_p2 : {}...".format(clump_p2),verbose=verbose)
    log.write("  -clump_kb : {}...".format(clump_kb),verbose=verbose)
    log.write("  -clump_r2 : {}...".format(clump_r2),verbose=verbose)
    if scaled == True:
        log.write(" -Clumping will be performed using {}".format(mlog10p),verbose=verbose)
        clump_log10_p1=-np.log10(clump_p1)
        clump_log10_p2=-np.log10(clump_p2)
        log.write("  -clump_log10_p1 : {}...".format(clump_log10_p1),verbose=verbose)
        log.write("  -clump_log10_p2 : {}...".format(clump_log10_p2),verbose=verbose)
        sumstats = insumstats.loc[insumstats[mlog10p]>min(clump_log10_p1,clump_log10_p2),:].copy()
    # extract lead variants
    else:
        log.write(" -Clumping will be performed using {}".format(p),verbose=verbose)
        sumstats = insumstats.loc[insumstats[p]<max(clump_p1,clump_p2),:].copy()
    log.write(" -Significant variants on CHR: ",list(sumstats["CHR"].unique()),verbose=verbose)
    
    plink_log=""

    # process reference file
    bfile, plink_log, ref_bim,filetype = _process_plink_input_files(chrlist=sumstats["CHR"].unique(), 
                                                       bfile=bfile, 
                                                       vcf=vcf, 
                                                       n_cores=n_cores,
                                                       plink_log=plink_log, 
                                                       log=log,
                                                       overwrite=overwrite)           
    
    ## process sumstats by CHR
    for i in sumstats["CHR"].unique():
        log.write(" -Processing sumstats for CHR {}...".format(i),verbose=verbose)
        
        if "@" in bfile:
            bfile_to_use = bfile.replace("@",str(i))
        else:
            bfile_to_use = bfile
        
        # checking # variants
        try:
            if filetype=="bfile":
                bim = pd.read_csv(bfile_to_use + ".bim",usecols=[1],header=None,sep="\s+")[1]
            else:
                bim = pd.read_csv(bfile_to_use + ".pvar",usecols=[2],header=None,comment="#",sep="\s+")[2]
            
            snplist = sumstats.loc[sumstats["CHR"]==i,"SNPID"]
            
            is_on_both = sumstats["SNPID"].isin(bim)

            log.write(" -Variants in reference file: {}...".format(len(bim)),verbose=verbose)
            log.write(" -Variants in sumstats: {}...".format(len(snplist)),verbose=verbose)
            log.write(" -Variants available in both reference and sumstats: {}...".format(sum(is_on_both)),verbose=verbose)

            is_avaialable_variant = (sumstats["CHR"]==i) & (is_on_both)

            if scaled == True:
                sumstats.loc[is_avaialable_variant,["SNPID",mlog10p]].to_csv("_gwaslab_tmp.{}.SNPIDP".format(i),index=False,sep="\t")
            else:
                sumstats.loc[is_avaialable_variant,["SNPID",p]].to_csv("_gwaslab_tmp.{}.SNPIDP".format(i),index=False,sep="\t")
        except:
            log.write(" -Not available for: {}...".format(i),verbose=verbose)
        
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

        log.write(" -Performing clumping for CHR {}...".format(i),verbose=verbose)
        log = _checking_plink_version(plink2=plink2, log=log)
        if memory is not None:
            memory_flag = "--memory {}".format(memory)
        
        if filetype=="bfile":
            file_flag = "--bfile {}".format(bfile_to_use) 
        else:
            file_flag = "--pfile {}".format(bfile_to_use) 
    
        if scaled == True:
            # clumping using LOG10P
            script = """
            {} \
                {}\
                --chr {} \
                --clump {} \
                --clump-log10 \
                --clump-field {} \
                --clump-snp-field SNPID \
                --clump-log10-p1 {} \
                --clump-log10-p2 {} \
                --clump-r2 {} \
                --clump-kb {} \
                --threads {} {}\
                --out {}
            """.format(plink2, file_flag, chrom, clump, mlog10p,clump_log10_p1, clump_log10_p2, clump_r2, clump_kb, n_cores, memory_flag if memory is not None else "", out_single_chr)    
        else:
            # clumping using P
            script = """
            {} \
                {}\
                --chr {} \
                --clump {} \
                --clump-field {} \
                --clump-snp-field SNPID \
                --clump-p1 {} \
                --clump-p2 {} \
                --clump-r2 {} \
                --clump-kb {} \
                --threads {} {}\
                --out {}
            """.format(plink2,file_flag, chrom, clump, p, clump_p1, clump_p2, clump_r2, clump_kb, n_cores,memory_flag if memory is not None else "", out_single_chr)
        
        try:
            output = subprocess.check_output(script, stderr=subprocess.STDOUT, shell=True,text=True)
            log.write(" -Saved results for CHR {} to : {}".format(i,"{}.clumps".format(out_single_chr)),verbose=verbose)
            plink_log +=output + "\n"
        except subprocess.CalledProcessError as e:
            log.write(e.output)
        #os.system(script)
        
        clumped = pd.read_csv("{}.clumps".format(out_single_chr),sep="\s+")
        results = pd.concat([results,clumped],ignore_index=True)
        
        # remove temp SNPIDP file 
        os.remove(clump)
    
    results = results.sort_values(by=["#CHROM","POS"]).rename(columns={"#CHROM":"CHR","ID":"SNPID"})
    log.write("Finished clumping.",verbose=verbose)
    results_sumstats = insumstats.loc[insumstats["SNPID"].isin(results["SNPID"]),:].copy()
    finished(log=log, verbose=verbose, end_line=_end_line)

    return results_sumstats, results, plink_log



       