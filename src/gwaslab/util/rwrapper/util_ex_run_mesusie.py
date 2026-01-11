from typing import Optional, Tuple, List
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_susie_version
from gwaslab.viz.viz_plot_stackedregional import _sort_kwargs

def _run_mesusie(
    filepath: Optional[str],
    r: str = "Rscript",
    types: Optional[Tuple[str, str]] = None,
    ns: Optional[List[int]] = None,
    fillldna: bool = True,
    delete: bool = False,
    coloc_kwargs: str = "",
    susie_kwargs: str = "",
    ncols: Optional[List[int]] = None,
    d1_kwargs: str = "",
    d2_kwargs: str = "",
    log: Log = Log(),
    verbose: bool = True
) -> str:
    
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
        if index==0:
            study0 = row["STUDY"]
        study = row["STUDY"]
        group = row["GROUP"]
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
ld_list${study} <-  as.matrix(ld{index})
summ_stat_list${study} <- sum{index}

png(filename="./diagnostic_{group}_{locus}_{index}.png")
diagnostic <- kriging_rss(summ_stat_list${study}$Z, ld_list${study})
diagnostic$plot
dev.off()
        '''.format(
             index = index,
             study = study,
             group=group,
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
saveRDS(MESuSiE_res, file = "{group}_{locus}.rds")
pips <- cbind(summ_stat_list${study0}$SNP, summ_stat_list${study0}$CHR, summ_stat_list${study0}$POS, MESuSiE_res$pip_config)
colnames(pips)[1] <-"SNPID"
colnames(pips)[2] <-"CHR"
colnames(pips)[3] <-"POS"
pips <- data.frame(pips)
pips[c("CREDIBLE_SET_INDEX")] <- 0 
pips[c("CS_CATEGORY")] <- NA
for (i in 1:length(MESuSiE_res$cs$cs)) {{
    pips[MESuSiE_res$cs$cs[[i]],c("CREDIBLE_SET_INDEX")]<-i
    pips[MESuSiE_res$cs$cs[[i]],c("CS_CATEGORY")] <- MESuSiE_res$cs$cs_category[[i]]
}}
write.csv(pips, "{group}_{locus}.pipcs", row.names = FALSE)

write.csv(MESuSiE_res$cs$cs_index, "{group}_{locus}.cscs_index", row.names = FALSE)
write.csv(MESuSiE_res$cs$purity, "{group}_{locus}.cspurity", row.names = FALSE)
write.csv(MESuSiE_res$cs$cs_category, "{group}_{locus}.cscs_category", row.names = FALSE)


for (p in MESuSiE_res$cs$cs) {{
  write(p,"{group}_{locus}.cscs_i", append=TRUE, sep="\t", ncolumns=10000000)
  write(summ_stat_list${study0}$SNP[p],"{group}_{locus}.cscs_snpid", append=TRUE, sep="\t", ncolumns=10000000)
}}

    '''.format(group=group,locus=locus,study0=study0)
    
    rscript_plotting='''
png(filename="./{group}_{locus}_stacked_regions.png")
MESuSiE_Plot(MESuSiE_res, ld_list ,summ_stat_list)
dev.off()
'''.format(group=group,locus=locus)
    
    rscript = rscript_loading + rscript_computing + rscript_output + rscript_plotting

    log.write("  -MESuSie script: {}".format(rscript_computing), verbose=verbose)
    
    with open("_{}_{}_gwaslab_mesusie_temp.R".format(group,row["SNPID"]),"w") as file:
            file.write(rscript)

    script_run_r = "{} _{}_{}_gwaslab_mesusie_temp.R".format(r, group,row["SNPID"])
    
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
    return "./{}_@.pipcs".format(group)
