from typing import TYPE_CHECKING, Optional, List
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_path_manager import _path

if TYPE_CHECKING:
    from gwaslab.g_SumstatsMulti import SumstatsMulti

def _run_mtag(sumstats_multi: 'SumstatsMulti', 
                     python: str = "python",
                     mtag: str = "",
                     study: str = "Group1",
                     special_flags: str = "",
                     ld_ref_panel: Optional[str] = None,
                     traits: Optional[List[str]] = None,
                     out_prefix: Optional[str] = None,
                     perfect_gencov: bool = False, 
                     equal_h2: bool = False, 
                     no_overlap: bool = False, 
                     fdr: bool = False,
                     n_min: int = 0,
                     nstudy: int = 2,
                     log: Log = Log(), 
                     verbose: bool = True) -> None:
    
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
        csv_path = _path(study = study,
                    trait = traits_to_form_string[i],
                    suffix="tsv.gz")

        sumstats_multi.data[output_snp_info_cols+ output_stats_cols].rename(columns=rename_dict).to_csv(csv_path, index=None,sep="\t")
        sumstats_paths.append(csv_path)
    
    sumstats_multi.offload()
    
    python_log=""
    if out_prefix is None:
        out_prefix = _path(study=study, 
                           nstudy = nstudy)
        
        #out_prefix = "./{study}_{nstudy}studies".format(study=study, nstudy=nstudy)
    if ld_ref_panel is not None:
        ld_ref_flag = "--ld_ref_panel {}".format(ld_ref_panel)
    else:
        ld_ref_flag=""

    if perfect_gencov == True:
        special_flags += "--perfect_gencov "
    if equal_h2 == True:
        special_flags += "--equal_h2 "
    if no_overlap == True:
        special_flags += "--no_overlap "
    if fdr == True:
        special_flags += "--fdr "

    script='''
{python} {mtag} {special_flags} {ld_ref_flag} \
--sumstats {sumstats_paths_string} \
--out {out_prefix} \
--n_min {n_min} \
--stream_stdout &
        '''.format(
            python=python,
            n_min=n_min,
            mtag=mtag,
            special_flags=special_flags,
            out_prefix=out_prefix,
            ld_ref_flag=ld_ref_flag,
            sumstats_paths_string = ",".join(sumstats_paths)
        )
    log.write("MTAG script: {} ".format(script), verbose=verbose)

    temp_script_path = _path(tmp=True,
                             study=study,
                             analysis="mtag", 
                             suffix="sh"
    )

    with open(temp_script_path,"w") as file:
            file.write(script)
    
    os.chmod(temp_script_path, 0o700)   
    
    try:
        log.write(" -Running MTAG from command line...", verbose=verbose)
        output = subprocess.check_output(os.path.join(temp_script_path)
                                         ,stderr=subprocess.STDOUT, shell=True,text=True)
        log.write(output)
        python_log+= output + "\n"

    except subprocess.CalledProcessError as e:
        log.write(e.output)

    sumstats_multi.reload()
    
    log.write("Finished MTAG.", verbose=verbose)
