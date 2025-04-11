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
from gwaslab.viz_plot_stackedregional import _sort_args

def _run_mesusie(filepath, 
                 r="Rscript",
                 types=None, ns=None, 
                 fillldna=True, delete=False, 
                coloc_args="", 
                susie_args="", 
                ncols=None,
                d1_args="",
                d2_args="",
                log=Log(), 
                verbose=True):
    
    log.write(" Start to run mesusie from command line:", verbose=verbose)
    pass

    if ns is None:
        if ncols is not None:
            ns = ncols
    log.write(" -Ns: {} and {}".format(ns[0],ns[1]), verbose=verbose)

    if filepath is None:
        log.write(" -File path is None.", verbose=verbose)
        log.write("Finished finemapping using MESuSie.", verbose=verbose)
        return pd.DataFrame()
        
    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()

    log = _checking_r_version(r, log)
    #log = _check_susie_version(r,log)
    r_script_init='''
library(MESuSiE)
ld_list <- list()
summ_stat_list <- list()
    '''
    r_scripts_for_loading =[r_script_init]
    
    for index, row in filelist.iterrows(): 
        gc.collect()
        study = row["STUDY"]
        ld_r_matrix = row["LD_R_MATRIX"]
        sumstats = row["LOCUS_SUMSTATS"]
        locus=row["LOCUS"]

        log.write(" -Running for: {} - {}".format(row["SNPID"],row["STUDY"] ), verbose=verbose)
        log.write("  -Locus sumstats:{}".format(sumstats), verbose=verbose)
        log.write("  -LD r matrix:{}".format(ld_r_matrix), verbose=verbose)

        rscript='''
sum{index} <-  read.csv("{sumstats}",sep="\\t")
sum{index}$Z <- sum{index}$Beta/sum{index}$Se
sum{index}$N <- {n}
ld{index} <- read.csv("{ld_r_matrix}",sep="\\t",header=FALSE)
ld{index}[is.na(ld{index})]  <- 0
names(ld{index}) <- sum{index}$SNP
ld_list$study{index} <-  as.matrix(ld{index})
summ_stat_list$study{index} <- sum{index}

png(filename="./diagnostic_{study}_{locus}_{index}.png")
diagnostic <- kriging_rss(summ_stat_list$study{index}$Z, ld_list$study{index})
diagnostic$plot
dev.off()
        '''.format(
             index = index,
             study = study,
             locus = locus,
             n = ns[index],
             sumstats = sumstats,
             ld_r_matrix = ld_r_matrix
        )
        r_scripts_for_loading.append(rscript)
    
    rscript_loading = "".join(r_scripts_for_loading)
    

    rscript_computing='''
MESuSiE_res<-meSuSie_core(ld_list, summ_stat_list, L=10)'''

    rscript_output = '''
saveRDS(MESuSiE_res, file = "{study}_{locus}.rds")
pips <- cbind(summ_stat_list$study0$SNP, MESuSiE_res$pip_config)
colnames(pips)[1] <-"SNPID"
write.csv(pips, "{study}_{locus}.pipcs", row.names = FALSE)

write.csv(MESuSiE_res$cs$cs_index, "{study}_{locus}.cscs_index", row.names = FALSE)
write.csv(MESuSiE_res$cs$purity, "{study}_{locus}.cspurity", row.names = FALSE)
write.csv(MESuSiE_res$cs$cs_category, "{study}_{locus}.cscs_category", row.names = FALSE)
write.csv(MESuSiE_res$cs$cs, "{study}_{locus}.cscs", row.names = FALSE) 
    '''.format(study=study,locus=locus)
    
    rscript_plotting='''
png(filename="./{study}_{locus}_stacked_regions.png")
MESuSiE_Plot(MESuSiE_res, ld_list ,summ_stat_list)
dev.off()
'''.format(study=study,locus=locus)
    
    rscript = rscript_loading + rscript_computing + rscript_output + rscript_plotting

    log.write("  -MESuSie script: {}".format(rscript_computing), verbose=verbose)
    
    with open("_{}_{}_gwaslab_mesusie_temp.R".format(study,row["SNPID"]),"w") as file:
            file.write(rscript)

    script_run_r = "{} _{}_{}_gwaslab_mesusie_temp.R".format(r, study,row["SNPID"])
    
    try:
        output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
        #plink_process = subprocess.Popen("exec "+script_run_r, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
        #output1,output2 = plink_process.communicate()
        #output= output1 + output2+ "\n"
        #plink_process.kill()
        log.write(" Running MESuSie from command line...", verbose=verbose)
        r_log+= output + "\n"
        
        #os.remove("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]))
        
    except subprocess.CalledProcessError as e:
        log.write(e.output)
        #os.remove("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]))
    log.write("Finished cross ancestry finemapping using MESuSie.", verbose=verbose)
    return locus_pip_cs