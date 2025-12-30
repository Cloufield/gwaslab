from typing import TYPE_CHECKING, Optional, Union
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_susie_version
from gwaslab.qc.qc_decorator import with_logging

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

@with_logging(
        start_to_msg="run finemapping using SuSieR from command line",
        finished_msg="running finemapping using SuSieR from command line",
        start_cols=["SNPID","CHR","POS"],
        start_function=".run_susie_rss()"
)
def _run_susie_rss(gls: 'Sumstats',
                   filepath: Optional[str], 
                   r: str = "Rscript", 
                   mode: str = "bs",
                   out: Optional[str] = None,
                   max_iter: int = 100,
                   min_abs_corr: float = 0.5,
                   refine: str = "FALSE",
                   L: int = 10,
                   fillldna: bool = True, 
                   n: Optional[Union[int, str]] = None, 
                   delete: bool = False,  #if delete output file
                   susie_kwargs: str = "", 
                   log: Log = Log(),
                   verbose: bool = True) -> pd.DataFrame:
    if filepath is None:
        log.write(" -File path is None.")
        log.write("Finished finemapping using SuSieR.")
        return pd.DataFrame()
    
    gls.offload()

    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()

    log = _checking_r_version(r, log)
    log = _check_susie_version(r,log)
    
    for index, row in filelist.iterrows(): 
        gc.collect()
        study = row["STUDY"]
        ld_r_matrix = row["LD_R_MATRIX"] #ld matrix path
        sumstats = row["LOCUS_SUMSTATS"] #sumsttas path
        
        # out: directory for output files
        if out is None:
            output_prefix = sumstats.replace(".sumstats.gz","")
        else:
            output_prefix = os.path.join(out, os.path.basename(sumstats.replace(".sumstats.gz","")))
        
        log.write(" -Running for: {} - {}".format(row["SNPID"],row["STUDY"] ))
        log.write("  -Locus sumstats:{}".format(sumstats))
        log.write("  -LD r matrix:{}".format(ld_r_matrix))
        log.write("  -output_prefix:{}".format(output_prefix))
        
        rscript='''
library(susieR)

sumstats <- read.csv("{}",sep="\t")

R <- as.matrix(read.csv("{}",sep="\t",header=FALSE))
{}

n <- floor(mean(sumstats$N))

fitted_rss1 <- susie_rss({}, n = {}, R = R, max_iter = {}, min_abs_corr={}, refine = {}, L = {}{})

susie_fitted_summary <- summary(fitted_rss1)

output <- susie_fitted_summary$vars
output$SNPID <- sumstats$SNPID[susie_fitted_summary$vars$variable]
output$LOCUS <- "{}"
output$STUDY <- "{}"

write.csv(output, "{}.pipcs", row.names = FALSE)

png(filename="{}_diagnostic.png")
diagnostic <- kriging_rss({}, R, n=n)
diagnostic$plot
dev.off()
        '''.format(sumstats, 
                   ld_r_matrix,
                    "R[is.na(R)] <- 0" if fillldna==True else "",
                    "z= sumstats$Z," if mode=="z" else "bhat = sumstats$BETA,shat = sumstats$SE",
                    n if n is not None else "n", 
                    max_iter,
                    min_abs_corr, 
                    refine, 
                    L, 
                    susie_kwargs, 
                    row["SNPID"],
                    row["STUDY"],
                    output_prefix,
                    output_prefix,
                    "sumstats$Z" if mode=="z" else "sumstats$BETA/sumstats$SE")
        susier_line = "susie_rss({}, n = {}, R = R, max_iter = {}, min_abs_corr={}, refine = {}, L = {}{})".format("z= sumstats$Z," if mode=="z" else "bhat = sumstats$BETA,shat = sumstats$SE,",
                    n if n is not None else "n", 
                    max_iter,
                    min_abs_corr, 
                    refine, 
                    L, 
                    susie_kwargs)
        log.write("  -SuSieR script: {}".format(susier_line))
        
        # temporary R script path
        temp_r_path = "_{}_{}_{}_gwaslab_susie_temp.R".format(study,row["SNPID"],id(sumstats))
        if out is not None:
            temp_r_path = os.path.join(out, temp_r_path)
        
        
        log.write("  -Createing temp R script: {}".format(temp_r_path))
        with open(temp_r_path,"w") as file:
                file.write(rscript)

        script_run_r = "{} {}".format(r, temp_r_path)
        
        try:
            log.write("  -Running SuSieR from command line...")
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            #plink_process = subprocess.Popen("exec "+script_run_r, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
            #output1,output2 = plink_process.communicate()
            #output= output1 + output2+ "\n"
            #plink_process.kill()
            
            r_log+= output + "\n"
            pip_cs = pd.read_csv("{}.pipcs".format(output_prefix))
            pip_cs["LOCUS"] = row["SNPID"]
            pip_cs["STUDY"] = row["STUDY"]
            locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)
            	
            os.remove(temp_r_path)
            log.write("  -Removing temp R script: {}".format(temp_r_path))

            if delete == True:
                os.remove("{}.pipcs".format(output_prefix))
                log.write("  -Removing output file: {}".format(temp_r_path))
            else:
                log.write("  -SuSieR result summary to: {}".format("{}.pipcs".format(output_prefix)))
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            os.remove(temp_r_path)
            log.write("  -Removing temp R script: {}".format(temp_r_path))
    
    gls.reload()

    locus_pip_cs = locus_pip_cs.rename(columns={"variable":"N_SNP","variable_prob":"PIP","cs":"CREDIBLE_SET_INDEX"})	
    locus_pip_cs = pd.merge(locus_pip_cs, gls.data[["SNPID","CHR","POS"]], on="SNPID",how="left")
    
    return locus_pip_cs

def _get_cs_lead(pipcs: pd.DataFrame) -> pd.DataFrame:
    leads = pipcs.loc[pipcs["CREDIBLE_SET_INDEX"]>0,:]
    leads = leads.sort_values(by="PIP",ascending=False).drop_duplicates(subset=["STUDY","LOCUS","CREDIBLE_SET_INDEX"])
    return leads
