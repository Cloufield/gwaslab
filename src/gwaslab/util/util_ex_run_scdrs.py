from typing import TYPE_CHECKING, Optional
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _run_scdrs(
    gls: 'Sumstats',
    scdrs: str = "scdrs",
    python: str = "python",
    study: str = "Study1",
    conda_env: Optional[str] = None,
    zscore_file: Optional[str] = None,
    out_file: Optional[str] = None,
    h5ad_file: Optional[str] = None,
    out_folder: Optional[str] = None,
    exclude_hla: bool = True,
    group_analysis: Optional[str] = None,
    gene_analysis: bool = False,
    gs_species: str = "human",
    h5ad_species: str = "human",
    flag_filter_data: bool = True,
    flag_raw_count: bool = True,
    munge_gs: bool = True,
    compute_score: bool = True,
    perform_downstream: bool = True,
    out: str = "./",
    delete: bool = True,
    ncol: str = "N",
    build: str = "19",
    log: Log = Log(),
    verbose: bool = True
) -> None:
    
    log.write(" Start to run scDRS from command line:", verbose=verbose)
    
    log.write(f" Output prefix: {out}", verbose=verbose)
    gls.offload()
    trait = study
    
    if out_file is None:
        out_file = f"./{trait}.gs"
        out_file = os.path.join(out, out_file)
    if out_folder is None:
        out_folder = out

    if conda_env is not None:
        conda_env_string = f"conda init bash\n conda activate {conda_env}\n"
    else:
        conda_env_string=""
    

    if group_analysis is not None:
        analysis_string = f"--group-analysis {group_analysis} "
    if gene_analysis == True:
        analysis_string += "--gene-analysis"
    
    bash_script=f'''#!/bin/bash
{conda_env_string} 
'''

    if munge_gs==True:
    
        bash_script+=f'''
{python} {scdrs} munge-gs \
    --out-file {out_file} \
    --zscore-file {zscore_file} \
    --weight zscore \
    --n-max 1000

'''
    
    if compute_score==True:
    
        bash_script+=f'''
{python} {scdrs} compute-score \
    --h5ad-file {h5ad_file} \
    --h5ad-species {h5ad_species} \
    --gs-file {out_file} \
    --gs-species {gs_species} \
    --out-folder {out_folder} \
    --flag-filter-data {flag_filter_data} \
    --flag-raw-count {flag_raw_count} \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

'''
    
    if perform_downstream==True:
        bash_script+=f'''
{python} {scdrs} perform-downstream {analysis_string} \
    --h5ad-file {h5ad_file} \
    --score-file ./{trait}.full_score.gz \
    --out-folder {out_folder} \
    --min_genes 250 \
    --min_cells 50 \
    --knn_n_neighbors 15 \
    --knn_n_pcs 20 \
    --flag-filter-data {flag_filter_data} \
    --flag-raw-count {flag_raw_count} 

'''
    log.write(f"Script: {bash_script}")

    try:
        log.write(" Running scDRS from command line...", verbose=verbose)
        output = subprocess.check_output(bash_script, stderr=subprocess.STDOUT, shell=True,text=True)
        output =  output + "\n"
        
    except subprocess.CalledProcessError as e:
        log.warning("ERROR!")
        log.write(e.output)
    gls.reload()
    log.write("Finished running scDRS.", verbose=verbose)