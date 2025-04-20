import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log

def _run_mtag(       sumstats_multi, 
                     python="Rscript",
                     mtag="",
                     study="Group1",
                     traits=None,
                     out_prefix=None,
                     types=None, 
                     n_min=0,
                     loci=None,
                     nstudy=2,
                     windowsizekb=1000,
                     build="99",
                     log=Log(), 
                     verbose=True):
    
    log.write("Start to run MTAG from command line:", verbose=verbose)

    if traits is None:
         traits_to_form_string = [ 'trait_{}'.format(i+1) for i in range(nstudy)]
    else:
         traits_to_form_string = ['{}'.format(i) for i in traits]
    
    res_combined = pd.DataFrame()
    # snpid    chr    bpos    a1    a2    freq    z    pval    n
    
    output_snp_info_cols =["rsID","CHR","POS","EA","NEA"]
    sumstats_paths = []
    for i in range(nstudy):
        output_stats_cols=[]
        for col in ["Z","P","EAF","N"]:
            output_stats_cols.append("{}_{}".format(col, i+1))

        rename_dict = {
             "rsID":"snpid",
             "CHR":"chr",
             "POS":"bpos",
             "EA":"a1",
             "NEA":"a2",
             "EAF_{}".format( i+1) :"freq",
             "Z_{}".format( i+1)   :"z",
             "P_{}".format( i+1)   :"pval",
             "N_{}".format( i+1)   :"n",

        }
    
        sumstats_multi[output_snp_info_cols+ output_stats_cols].rename(columns=rename_dict).to_csv("{}_{}.tsv.gz".format(study, traits_to_form_string[i]), index=None,sep="\t")
        sumstats_paths.append("{}_{}.tsv.gz".format(study, traits_to_form_string[i]))
    
    python_log=""
    if out_prefix is None:
        out_prefix = "./{study}_{nstudy}studies".format(study=study, nstudy=nstudy)
    
    script='''
{python} {mtag} \
--sumstats {sumstats_paths_string} \
--out {out_prefix} \
--n_min {n_min} \
--stream_stdout &
        '''.format(
            python=python,
            n_min=n_min,
            mtag=mtag,
            out_prefix=out_prefix,
            sumstats_paths_string = ",".join(sumstats_paths)
        )
    log.write(" MTAG script: {} ".format(script), verbose=verbose)

    
    with open("_{}_gwaslab_mtag_temp.sh".format(study),"w") as file:
            file.write(script)
    
    os.chmod("_{}_gwaslab_mtag_temp.sh".format(study), 0o700)   

    script_run = "./_{}_gwaslab_mtag_temp.sh".format(study)
    
    try:
        log.write(" Running MTAG from command line...", verbose=verbose)
        output = subprocess.check_output(script_run, stderr=subprocess.STDOUT, shell=True,text=True)
        log.write(output)
        python_log+= output + "\n"

    except subprocess.CalledProcessError as e:
        log.write(e.output)

    log.write("Finished MTAG.", verbose=verbose)
