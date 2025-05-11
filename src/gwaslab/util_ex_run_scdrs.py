import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log

def _run_scdrs( scdrs="scdrs",
                python="python",
                study="Study1",
                conda_env=None,
                zscore_file=None,
                out_file=None,
                h5ad_file=None,
                out_folder=None,
                exclude_hla=True,
                group_analysis = None,
                gene_analysis = False,
                gs_species="human",
                h5ad_species="human",
                flag_filter_data=True,
                flag_raw_count=True,
                munge_gs=True,
                compute_score=True,
                perform_downstream=True,
                out="./",
                delete=True, 
                ncol="N",
                build="19",
                log=Log(), 
                verbose=True):
    
    log.write(" Start to run scDRS from command line:", verbose=verbose)
    
    trait = study
    if out_file is None:
        out_file = f"./{trait}.gs"
    if out_folder is None:
        out_folder = f"./"
    if conda_env is not None:
        conda_env_string = f"conda init bash\n conda activate {conda_env}\n"
    else:
        conda_env_string=""
    log.write(f" Output prefix: {out}", verbose=verbose)

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

    log.write("Finished running scDRS.", verbose=verbose)