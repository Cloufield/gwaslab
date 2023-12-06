import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.g_version import _checking_r_version
from gwaslab.g_version import _check_susie_version

def _run_coloc_susie(filepath, r="Rscript",types=None, ns=None, fillldna=True, delete=False, coloc_args="", susie_args="", log=Log()):
    log.write(" Start to run coloc.susie from command line:")

    if types is None:
         types = ("cc","cc")

    if filepath is None:
        log.write(" -File path is None.")
        log.write("Finished finemapping using SuSieR.")
        return pd.DataFrame()
        
    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()

    log = _checking_r_version(r, log)
    #log = _check_susie_version(r,log)

    for index, row in filelist.iterrows(): 
        gc.collect()
        study = row["study"]
        ld_r_matrix = row["LD_r_matrix"]
        sumstats = row["Locus_sumstats"]
        output_prefix = sumstats.replace(".sumstats.gz","")
        log.write(" -Running for: {} - {}".format(row["SNPID"],row["study"] ))
        log.write("  -Locus sumstats:{}".format(sumstats))
        log.write("  -LD r matrix:{}".format(ld_r_matrix))
        log.write("  -output_prefix:{}".format(output_prefix))
        
        rscript='''
        library(coloc)
        
        df = read.csv("{}",header=TRUE)

        R <- as.matrix(read.csv("{}",sep="\t",header=FALSE))

        rownames(R)<-df[,"SNPID"]
        colnames(R)<-df[,"SNPID"]
        
        {}

        D1 <- list( "LD"=R, "beta"=df[,"BETA_1"],"varbeta"=df[,"SE_1"]**2,"snp"=df[,"SNPID"],"position"=df[,"POS"],"type"="{}","N"={})
        D2 <- list( "LD"=R, "beta"=df[,"BETA_2"],"varbeta"=df[,"SE_2"]**2,"snp"=df[,"SNPID"],"position"=df[,"POS"],"type"="{}","N"={})

        S1=runsusie(D1)
        S2=runsusie(D2)

        susie.res=coloc.susie(S1,S2{})

        write.csv(susie.res$summary, "{}.coloc.susie", row.names = FALSE)
        '''.format(sumstats, 
                   ld_r_matrix,
                    "R[is.na(R)] <- 0" if fillldna==True else "",
                    types[0], ns[0],
                    types[1], ns[1],
                    coloc_args,
                    output_prefix)
        
        log.write("  -coloc script: {}".format("coloc.susie(S1,S2)"))
        with open("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]),"w") as file:
                file.write(rscript)

        script_run_r = "{} _{}_{}_gwaslab_coloc_susie_temp.R".format(r, study,row["SNPID"])
        
        try:
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            #plink_process = subprocess.Popen("exec "+script_run_r, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
            #output1,output2 = plink_process.communicate()
            #output= output1 + output2+ "\n"
            #plink_process.kill()
            log.write(" Running coloc.SuSieR from command line...")
            r_log+= output + "\n"
            pip_cs = pd.read_csv("{}..coloc.susie".format(output_prefix))
            pip_cs["Locus"] = row["SNPID"]
            pip_cs["STUDY"] = row["study"]
            locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)
            os.remove("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]))
            if delete == True:
                os.remove("{}.pipcs".format(output_prefix))
            else:
                log.write("  -SuSieR result summary to: {}".format("{}.pipcs".format(output_prefix)))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            os.remove("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]))
    log.write("Finished finemapping using SuSieR.")
    return locus_pip_cs