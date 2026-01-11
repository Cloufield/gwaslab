from typing import TYPE_CHECKING, Optional, Tuple, List, Any
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_susie_version

if TYPE_CHECKING:
    from gwaslab.g_SumstatsPair import SumstatsPair

def _run_coloc_susie(
    glsp: 'SumstatsPair',
    filepath: Optional[str],
    r: str = "Rscript",
    types: Optional[Tuple[str, str]] = None,
    ns: Optional[Tuple[int, int]] = None,
    fillldna: bool = True,
    delete: bool = False,
    coloc_kwargs: str = "",
    susie_kwargs: str = "",
    ncols: Optional[Tuple[int, int]] = None,
    d1_kwargs: str = "",
    d2_kwargs: str = "",
    out: Optional[str] = None,
    log: Log = Log(),
    verbose: bool = True
) -> pd.DataFrame:
    
    log.write("Start to run coloc.susie from command line:", verbose=verbose)

    if filepath is None:
        log.write(" -File path is None.", verbose=verbose)
        log.write("Finished finemapping using SuSieR.", verbose=verbose)
        return pd.DataFrame()
    
    glsp.offload()

    if types is None:
        types = ("cc","cc")
    log.write(" -Phenotype types: {} and {}".format(types[0],types[1]), verbose=verbose)

    if ns is None:
        if ncols is not None:
            ns = ncols
    log.write(" -Ns: {} and {}".format(ns[0],ns[1]), verbose=verbose)

    filelist = pd.read_csv(filepath,sep="\t")
    r_log=""
    # write R script
    locus_pip_cs = pd.DataFrame()

    log = _checking_r_version(r, log)
    #log = _check_susie_version(r,log)

    for index, row in filelist.iterrows(): 
        gc.collect()
        study = row["STUDY"]
        ld_r_matrix = row["LD_R_MATRIX"]
        sumstats = row["LOCUS_SUMSTATS"]

        if out is None:
            output_prefix = sumstats.replace(".sumstats.gz","")
        else:
            output_prefix = os.path.join(out, os.path.basename(sumstats.replace(".sumstats.gz","")))
        
        log.write(" -Running for: {} - {}".format(row["SNPID"],row["STUDY"] ), verbose=verbose)
        log.write("  -Locus sumstats:{}".format(sumstats), verbose=verbose)
        log.write("  -LD r matrix:{}".format(ld_r_matrix), verbose=verbose)
        log.write("  -output_prefix:{}".format(output_prefix), verbose=verbose)
        
        rscript='''
        library(coloc)
        
        df = read.csv("{sumstats_path}",sep="\t",header=TRUE)

        R <- as.matrix(read.csv("{ld_r_matrix_path}",sep="\t",header=FALSE))

        rownames(R)<-df[,"SNPID"]
        colnames(R)<-df[,"SNPID"]
        
        {fillna_script}

        D1 <- list( "LD"=R, "beta"=df[,"BETA_1"],"varbeta"=df[,"SE_1"]**2,"snp"=df[,"SNPID"],"position"=df[,"POS"],"type"="{type1}","N"={n1}{d1_kwargs})
        D2 <- list( "LD"=R, "beta"=df[,"BETA_2"],"varbeta"=df[,"SE_2"]**2,"snp"=df[,"SNPID"],"position"=df[,"POS"],"type"="{type2}","N"={n2}{d2_kwargs})

        abf <- coloc.abf(dataset1=D1,dataset2=D2)
        write.csv(t(data.frame(abf$summary)) , "{output_prefix}.coloc.abf", row.names = FALSE)

        S1=runsusie(D1{susie_kwargs})
        S2=runsusie(D2{susie_kwargs})

        susie.res=coloc.susie(S1,S2{coloc_kwargs})

        write.csv(susie.res$summary, "{output_prefix}.coloc.susie", row.names = FALSE)

        '''.format(sumstats_path = sumstats, 
                   ld_r_matrix_path = ld_r_matrix,
                    fillna_script = "R[is.na(R)] <- 0" if fillldna==True else "",
                    type1 = types[0], 
                    n1 =ns[0],
                    d1_kwargs = d1_kwargs,
                    type2= types[1], 
                    d2_kwargs = d2_kwargs,
                    n2= ns[1],
                    susie_kwargs = susie_kwargs,
                    coloc_kwargs = coloc_kwargs,
                    output_prefix = output_prefix)
        
        log.write("  -coloc abf script: {}".format("coloc.abf(dataset1=D1,dataset2=D2)"), verbose=verbose)
        log.write("  -coloc susie script: {}".format("coloc.susie(S1,S2)"), verbose=verbose)
        
        with open("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]),"w") as file:
                file.write(rscript)

        script_run_r = "{} _{}_{}_gwaslab_coloc_susie_temp.R".format(r, study,row["SNPID"])
        
        try:
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            #plink_process = subprocess.Popen("exec "+script_run_r, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
            #output1,output2 = plink_process.communicate()
            #output= output1 + output2+ "\n"
            #plink_process.kill()
            log.write(" Running coloc.SuSieR from command line...", verbose=verbose)
            r_log+= output + "\n"
            
            pip_cs = pd.read_csv("{}.coloc.abf".format(output_prefix))
            if len(pip_cs)==0:
                 log.write("  -SuSieR result for {} is empty. Please check parameters.".format(output_prefix), verbose=verbose)
            else:
                pip_cs["LOCUS"] = row["SNPID"]
                pip_cs["STUDY"] = row["STUDY"]
                pip_cs["HIT1"] = row["SNPID"]
                pip_cs["METHOD"] = "abf"
                locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)

            pip_cs = pd.read_csv("{}.coloc.susie".format(output_prefix))
            if len(pip_cs)==0:
                 log.write("  -SuSieR result for {} is empty. Please check parameters.".format(output_prefix), verbose=verbose)
            else:
                pip_cs["LOCUS"] = row["SNPID"]
                pip_cs["STUDY"] = row["STUDY"]
                pip_cs["METHOD"] = "susie"
                locus_pip_cs = pd.concat([locus_pip_cs,pip_cs],ignore_index=True)
            
            os.remove("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]))
            
            if delete == True:
                os.remove("{}.coloc.susie".format(output_prefix))
                os.remove("{}.coloc.abf".format(output_prefix))
            else:
                log.write("  -coloc-abf result summary to: {}".format("{}.coloc.abf".format(output_prefix)), verbose=verbose)
                log.write("  -coloc-susie result summary to: {}".format("{}.coloc.susie".format(output_prefix)), verbose=verbose)
                
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            os.remove("_{}_{}_gwaslab_coloc_susie_temp.R".format(study,row["SNPID"]))
    
    log.write("Finished clocalization using coloc and SuSiE.", verbose=verbose)
    glsp.reload()
    return locus_pip_cs
