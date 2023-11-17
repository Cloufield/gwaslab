import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.Log import Log
from gwaslab.getsig import getsig
from gwaslab.processreference import _process_vcf_and_bfile

def _run_susie_rss(filepath, r="Rscript", max_iter=100000,min_abs_corr=0.1, refine="TRUE",L=10, susie_args="", log=Log()):
    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()
    for index, row in filelist.iterrows(): 
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
                                n = n, 
                                R = R, 
                                max_iter = {}, 
                                min_abs_corr={}, 
                                refine = {},
                                L = {}{})

        susie_fitted_summary <- summary(fitted_rss1)

        output <- susie_fitted_summary$vars
        output$SNPID <- sumstats$SNPID[susie_fitted_summary$vars$variable]

        write.csv(output, "{}.pipcs", row.names = FALSE)
        '''.format(sumstats, ld_r_matrix, max_iter,min_abs_corr, refine, L, susie_args, output_prefix)

        with open("_gwaslab_susie_temp.R","w") as file:
                file.write(rscript)

        script_run_r = "{} _gwaslab_susie_temp.R".format(r)
        
        try:
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            r_log+= output + "\n"
            pip_cs = pd.read_csv("{}.pipcs".format(output_prefix))
            pip_cs["Locus"] = row["SNPID"]
            locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)
            os.remove("_gwaslab_susie_temp.R")
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            os.remove("_gwaslab_susie_temp.R")

    return locus_pip_cs
