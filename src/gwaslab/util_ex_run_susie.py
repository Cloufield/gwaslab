import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.g_version import _checking_r_version
from gwaslab.g_version import _check_susie_version
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished

def _run_susie_rss(filepath, r="Rscript", mode="bs",max_iter=100000,min_abs_corr=0.1,refine="TRUE",L=10, fillldna=True, n=None, delete=False, susie_args="", log=Log(),verbose=True):
    ##start function with col checking##########################################################
    _start_line = "run finemapping using SuSieR from command line"
    _end_line = "running finemapping using SuSieR from command line"
    _start_cols =[]
    _start_function = ".run_susie_rss()"
    _must_args ={}

    is_enough_info = start_to(sumstats=None,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: raise ValueError("Not enough columns for calculating LD matrix")
    ############################################################################################
    if filepath is None:
        log.write(" -File path is None.")
        log.write("Finished finemapping using SuSieR.")
        return pd.DataFrame()
        
    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()

    log = _checking_r_version(r, log)
    log = _check_susie_version(r,log)

    for index, row in filelist.iterrows(): 
        gc.collect()
        study = row["STUDY"]
        ld_r_matrix = row["LD_R_MATRIX"]
        sumstats = row["LOCUS_SUMSTATS"]
        output_prefix = sumstats.replace(".sumstats.gz","")
        log.write(" -Running for: {} - {}".format(row["SNPID"],row["STUDY"] ))
        log.write("  -Locus sumstats:{}".format(sumstats))
        log.write("  -LD r matrix:{}".format(ld_r_matrix))
        log.write("  -output_prefix:{}".format(output_prefix))
        
        rscript='''
        library(susieR)
        
        sumstats <- read.csv("{}")
        
        R <- as.matrix(read.csv("{}",sep="\t",header=FALSE))
        {}

        n <- floor(mean(sumstats$N))

        fitted_rss1 <- susie_rss({}, n = {}, R = R, max_iter = {}, min_abs_corr={}, refine = {}, L = {}{})

        susie_fitted_summary <- summary(fitted_rss1)

        output <- susie_fitted_summary$vars
        output$SNPID <- sumstats$SNPID[susie_fitted_summary$vars$variable]

        write.csv(output, "{}.pipcs", row.names = FALSE)
        '''.format(sumstats, 
                   ld_r_matrix,
                    "R[is.na(R)] <- 0" if fillldna==True else "",
                    "z= sumstats$Z," if mode=="z" else "bhat = sumstats$BETA,shat = sumstats$SE",
                    n if n is not None else "n", 
                    max_iter,
                    min_abs_corr, 
                    refine, 
                    L, 
                    susie_args, 
                    output_prefix)
        susier_line = "susie_rss({}, n = {}, R = R, max_iter = {}, min_abs_corr={}, refine = {}, L = {}{})".format("z= sumstats$Z," if mode=="z" else "bhat = sumstats$BETA,shat = sumstats$SE,",
                    n if n is not None else "n", 
                    max_iter,
                    min_abs_corr, 
                    refine, 
                    L, 
                    susie_args)
        log.write("  -SuSieR script: {}".format(susier_line))
        with open("_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"]),"w") as file:
                file.write(rscript)

        script_run_r = "{} _{}_{}_gwaslab_susie_temp.R".format(r, study,row["SNPID"])
        
        try:
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            #plink_process = subprocess.Popen("exec "+script_run_r, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
            #output1,output2 = plink_process.communicate()
            #output= output1 + output2+ "\n"
            #plink_process.kill()
            log.write("  -Running SuSieR from command line...")
            r_log+= output + "\n"
            pip_cs = pd.read_csv("{}.pipcs".format(output_prefix))
            pip_cs["LOCUS"] = row["SNPID"]
            pip_cs["STUDY"] = row["STUDY"]
            locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)
            	
            os.remove("_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"]))
            if delete == True:
                os.remove("{}.pipcs".format(output_prefix))
            else:
                log.write("  -SuSieR result summary to: {}".format("{}.pipcs".format(output_prefix)))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            os.remove("_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"]))
    
    locus_pip_cs = locus_pip_cs.rename(columns={"variable":"N_SNP","variable_prob":"PIP","cs":"CREDIBLE_SET_INDEX"})	
    finished(log=log, verbose=verbose, end_line=_end_line)
    return locus_pip_cs

