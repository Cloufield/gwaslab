import subprocess
import numpy as np
import os
import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.util.util_ex_process_ref import _process_plink_input_files
from gwaslab.extension import _checking_plink_version
from gwaslab.qc.qc_decorator import with_logging
@with_logging(
        start_to_msg="perfrom clumping",
        finished_msg="clumping",
        start_cols=["SNPID","CHR","POS"],
        start_function=".clump()"
)
def _clump(gls, vcf=None, scaled=False, out="clumping_plink2", 
           p="P",mlog10p="MLOG10P", overwrite=False, study=None, bfile=None, pfile=None,
           n_cores=1, memory=None, chrom=None, clump_p1=5e-8, clump_p2=5e-8, clump_r2=0.01, clump_kb=250,
           log=Log(),verbose=True,plink="plink",plink2="plink2"):
    """
    Perform LD clumping of GWAS summary statistics using PLINK2.

    Parameters
    ----------
    vcf : str or None, optional
        Path or prefix to reference VCF or genotype data compatible with PLINK2.
        Used when deriving `--pfile` inputs.
    bfile : str or None, optional
        Prefix to PLINK binary files (`.bed/.bim/.fam`). May include "@"
        as a chromosome placeholder.
    pfile : str or None, optional
        Prefix to PLINK2 files (`.pgen/.pvar/.psam`). May include "@"
        as a chromosome placeholder.
    scaled : bool, optional
        If True, clump on `mlog10p` using PLINK2 `--clump-log10`. If False,
        clump on `p`.
    out : str or None, optional
        Output prefix. If None, uses "./{study}_clumpping".
    p : str, optional
        Column name of p-values in `gls.data`.
    mlog10p : str, optional
        Column name of -log10(p) in `gls.data`.
    overwrite : bool, optional
        Whether to overwrite any intermediate reference files produced while
        preparing inputs.
    study : str or None, optional
        Study name used when `out` is None.
    n_cores : int, optional
        Number of threads to pass to PLINK2 via `--threads`.
    memory : int or None, optional
        Memory limit (MB) for PLINK2 via `--memory`.
    chrom : any, optional
        Unused parameter kept for API compatibility.
    clump_p1 : float, optional
        Primary p-value threshold (`--clump-p1` or `--clump-log10-p1`).
    clump_p2 : float, optional
        Secondary p-value threshold (`--clump-p2` or `--clump-log10-p2`).
    clump_r2 : float, optional
        LD threshold (`--clump-r2`).
    clump_kb : int, optional
        Window size in kilobases (`--clump-kb`).
    log : gwaslab.g_Log.Log, optional
        Logger instance used for progress reporting.
    verbose : bool, optional
        Whether to emit verbose log messages.
    plink : str, optional
        Path to PLINK (v1). Not used directly in clumping.
    plink2 : str, optional
        Path to PLINK2 binary.

    Returns
    -------
    results_sumstats : pandas.DataFrame
        Subset of input summary statistics for clumped lead variants.
    results : pandas.DataFrame
        Concatenated PLINK2 `.clumps` output across processed chromosomes.
    plink_log : str
        Combined PLINK2 log output captured during execution.

    Notes
    -----
    - Writes temporary files named "{out}_gwaslab_tmp.{sumstats_id}.{chr}.SNPIDP",
      which are removed after clumping.
    - Produces per-chromosome output files "{out}.{chr}.clumps".

    Examples
    --------
    >>> results_sumstats, results, logstr = _clump(
    ...     bfile="ref/chr@",
    ...     clump_p1=5e-8,
    ...     clump_p2=1e-5,
    ...     clump_r2=0.1,
    ...     clump_kb=250,
    ...     n_cores=4
    ... )
    """
    ##start function with col checking##########################################################
    
    if out is None:
        out = f"./{study}_clumpping".lstrip('/')
    else:
        out = out.lstrip('/')
    sumstats_id = gls.id
    sumstats = gls.data
    gls.offload()

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
        sumstats = sumstats.loc[sumstats[mlog10p]>min(clump_log10_p1,clump_log10_p2),:].copy()
    # extract lead variants
    else:
        log.write(" -Clumping will be performed using {}".format(p),verbose=verbose)
        sumstats = sumstats.loc[sumstats[p]<max(clump_p1,clump_p2),:].copy()

    if len(sumstats)==0:
        log.write(" -No significant variants after filtering.")
        gls.reload()
        return pd.DataFrame(), pd.DataFrame(), ""
    
    log.write(" -Significant variants on CHR: ",list(sumstats["CHR"].unique()),verbose=verbose)
    
    plink_log=""

    # process reference file
    bfile, plink_log, ref_bim,filetype = _process_plink_input_files(chrlist=sumstats["CHR"].unique(), 
                                                       bfile=bfile, 
                                                       pfile=pfile,
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
                sumstats.loc[is_avaialable_variant,["SNPID",mlog10p]].to_csv("{}_gwaslab_tmp.{}.{}.SNPIDP".format(out, sumstats_id, i),index=False,sep="\t")
            else:
                sumstats.loc[is_avaialable_variant,["SNPID",p]].to_csv("{}_gwaslab_tmp.{}.{}.SNPIDP".format(out, sumstats_id,i),index=False,sep="\t")
        except:
            log.write(" -Not available for: {}...".format(i),verbose=verbose)
        
    # create a empty dataframe for combining results from each CHR 
    results = pd.DataFrame()

    
    # clumping using plink
    for i in sumstats["CHR"].unique():
        chrom = i
        # temp file  
        clump = "{}_gwaslab_tmp.{}.{}.SNPIDP".format(out,sumstats_id,chrom)
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
        
        try:
            clumped = pd.read_csv("{}.clumps".format(out_single_chr),sep="\s+")
            results = pd.concat([results,clumped],ignore_index=True)
            
        except:
            log.write(f"Clumping failed for chr{i}")
        # remove temp SNPIDP file 
        os.remove(clump)
    
    results = results.sort_values(by=["#CHROM","POS"]).rename(columns={"#CHROM":"CHR","ID":"SNPID"})
    results_sumstats = sumstats.loc[sumstats["SNPID"].isin(results["SNPID"]),:].copy()
    gls.reload()
    
    return results_sumstats, results, plink_log



       
