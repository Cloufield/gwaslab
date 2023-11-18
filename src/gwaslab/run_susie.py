import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.Log import Log
from gwaslab.getsig import getsig
from gwaslab.processreference import _process_vcf_and_bfile
from gwaslab.version import _checking_r_version
from gwaslab.version import _check_susie_version

def _run_susie_rss(filepath, r="Rscript", max_iter=100000,min_abs_corr=0.1,refine="TRUE",L=10, n=None, susie_args="", log=Log()):
    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()

    log = _checking_r_version(r, log)
    log = _check_susie_version(r,log)

    for index, row in filelist.iterrows(): 
        study = row["study"]
        ld_r_matrix = row["LD_r_matrix"]
        sumstats = row["Locus_sumstats"]
        output_prefix = sumstats.replace(".sumstats.gz","")
        log.write("Locus sumstats:{}".format(sumstats))
        log.write("LD r matrix:{}".format(ld_r_matrix))
        log.write("output_prefix:{}".format(output_prefix))
        
        rscript='''
        library(susieR)
        
        sumstats <- read.csv("{}")
        
        R <- as.matrix(read.csv("{}",sep="\t",header=FALSE))
        R[is.na(R)] <- 0

        n <- floor(mean(sumstats$N))

        fitted_rss1 <- susie_rss(bhat = sumstats$BETA, 
                                shat = sumstats$SE, 
                                n = {}, 
                                R = R, 
                                max_iter = {}, 
                                min_abs_corr={}, 
                                refine = {},
                                L = {}{})

        susie_fitted_summary <- summary(fitted_rss1)

        output <- susie_fitted_summary$vars
        output$SNPID <- sumstats$SNPID[susie_fitted_summary$vars$variable]

        write.csv(output, "{}.pipcs", row.names = FALSE)
        '''.format(sumstats, ld_r_matrix, n if n is not None else "n", max_iter,min_abs_corr, refine, L, susie_args, output_prefix)

        with open("_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"]),"w") as file:
                file.write(rscript)

        script_run_r = "{} _{}_{}_gwaslab_susie_temp.R".format(r, study,row["SNPID"])
        
        try:
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            r_log+= output + "\n"
            pip_cs = pd.read_csv("{}.pipcs".format(output_prefix))
            pip_cs["Locus"] = row["SNPID"]
            locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)
            os.remove("_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"]))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            os.remove("_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"]))

    return locus_pip_cs

