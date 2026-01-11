from typing import TYPE_CHECKING, Optional, List, Dict, Any, Union
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_susie_version
from gwaslab.util.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.util.util_in_get_sig import _get_sig

if TYPE_CHECKING:
    from gwaslab.g_SumstatsPair import SumstatsPair

def _run_ccgwas(
    sumstats_pair: 'SumstatsPair',
    r: str = "Rscript",
    group: str = "Group1",
    studies: Optional[List[str]] = None,
    traits: Optional[List[str]] = None,
    meta: Optional[Dict[str, Any]] = None,
    ldsc: Optional[List[pd.DataFrame]] = None,
    ldsc_rg: Optional[pd.DataFrame] = None,
    nstudy: int = 2,
    K_A1A0: Optional[float] = None,
    K_A1A0_high: Optional[float] = None,
    K_A1A0_low: Optional[float] = None,
    K_B1B0: Optional[float] = None,
    K_B1B0_high: Optional[float] = None,
    K_B1B0_low: Optional[float] = None,
    h2l_A1A0: Optional[float] = None,
    h2l_B1B0: Optional[float] = None,
    rg_A1A0_B1B0: Optional[float] = None,
    intercept_A1A0_B1B0: Optional[float] = None,
    m: float = 1e4,
    N_A1: Optional[int] = None,
    N_B1: Optional[int] = None,
    N_A0: Optional[int] = None,
    N_B0: Optional[int] = None,
    N_overlap_A0B0: int = 0,
    log: Log = Log(),
    verbose: bool = True
) -> str:
    
    log.write(" Start to run CCGWAS from command line:", verbose=verbose)
    log.write(" -Methods: : {}...".format("Peyrot, W. J., & Price, A. L. (2021). Identifying loci with different allele frequencies among cases of eight psychiatric disorders using CC-GWAS. Nature genetics, 53(4), 445-454."),verbose=verbose)
    #"SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff"
    log.write(" -Running CCGWAS for: {}...".format(group),verbose=verbose)
    
    snp_info_cols=["SNPID","CHR","POS","EA","NEA"]
    stats_cols=["EAF","OR","SE","P","N_EFF"]

    if meta["gwaslab"]["objects"][0]["gwaslab"]["population_prevalence"] != "Unknown":
          K_A1A0 = float(meta["gwaslab"]["objects"][0]["gwaslab"]["population_prevalence"])
          K_A1A0_high = 1.1 * K_A1A0
          K_A1A0_low = K_A1A0/ 1.1
    
    if meta["gwaslab"]["objects"][1]["gwaslab"]["population_prevalence"] != "Unknown":
          K_B1B0 = float(meta["gwaslab"]["objects"][1]["gwaslab"]["population_prevalence"])
          K_B1B0_high = 1.1 * K_B1B0
          K_B1B0_low = K_B1B0/ 1.1

    if h2l_A1A0 is None: 
        h2l_A1A0 = ldsc[0].loc[0, "h2_liab"]
    if h2l_B1B0 is None: 
        h2l_B1B0 = ldsc[1].loc[0, "h2_liab"]
    if rg_A1A0_B1B0 is None:
        rg_A1A0_B1B0 = ldsc_rg.loc[(ldsc_rg["p1"]==studies[0])&(ldsc_rg["p2"]==studies[1]), :].iloc[0,ldsc_rg.columns.get_loc("rg")]
    if intercept_A1A0_B1B0 is None:
        intercept_A1A0_B1B0 = ldsc_rg.loc[(ldsc_rg["p1"]==studies[0])&(ldsc_rg["p2"]==studies[1]), :].iloc[0,ldsc_rg.columns.get_loc("gcov_int")]

    # prepare input files sumstats_multi
    for i in range(nstudy):
        output_cols = snp_info_cols + list(map(lambda x: x+"_{}".format(i+1), stats_cols))

        dic= {"SNPID":"SNP",
              "POS":"BP",
              "EAF_{}".format(i+1):"FRQ",
              "OR_{}".format(i+1):"OR",
              "SE_{}".format(i+1):"SE",
              "P_{}".format(i+1):"P",
              "N_EFF_{}".format(i+1):"Neff"
              }

        sumstats_pair[output_cols].rename(columns=dic).to_csv("./{}_{}.txt.gz".format(group, studies[i]),index=None,sep="\t")
    
    output_prefix = "{group}_ccgwas".format(group=group)
    r_log=""
    log = _checking_r_version(r, log)

    rscript='''
library(data.table)
library(R.utils)
library(CCGWAS)

CCGWAS( outcome_file = "{output_prefix}" , 
        A_name = "{study1}" , 
        B_name = "{study2}" , 
        sumstats_fileA1A0 = "./{group}_{study1}.txt.gz" ,
        sumstats_fileB1B0 = "./{group}_{study2}.txt.gz" ,
        K_A1A0 = {K_A1A0} , 
        K_A1A0_high = {K_A1A0_high} , 
        K_A1A0_low = {K_A1A0_low},  

        K_B1B0 ={K_B1B0} , 
        K_B1B0_high ={K_B1B0_high} , 
        K_B1B0_low = {K_B1B0_low} , 

        h2l_A1A0 ={h2l_A1A0}, 
        h2l_B1B0 = {h2l_B1B0} , 
        rg_A1A0_B1B0 = {rg_A1A0_B1B0} , 

        intercept_A1A0_B1B0 = {intercept_A1A0_B1B0} , 
        m = {m} ,  
        N_A1 =  {N_A1} , 
        N_B1 =  {N_B1} , 
        N_A0 =  {N_A0} , 
        N_B0 =  {N_B0} , 
        N_overlap_A0B0 =  {N_overlap_A0B0} )
        '''.format(
            output_prefix=output_prefix,
            study1=studies[0],
            study2=studies[1],
            group=group,
            K_A1A0 = K_A1A0 , 
            K_A1A0_high = K_A1A0_high , 
            K_A1A0_low = K_A1A0_low ,  
            K_B1B0 =K_B1B0 , 
            K_B1B0_high = K_B1B0_high , 
            K_B1B0_low =K_B1B0_low , 
            h2l_A1A0 = h2l_A1A0 , 
            h2l_B1B0 = h2l_B1B0 , 
            rg_A1A0_B1B0 = rg_A1A0_B1B0 , 
            intercept_A1A0_B1B0 = intercept_A1A0_B1B0 , 
            m = m,
            N_A1 = N_A1,
            N_B1 = N_B1,
            N_A0 = N_A0 , 
            N_B0 = N_B0 , 
            N_overlap_A0B0 = N_overlap_A0B0,
        )

    with open("_{}_gwaslab_ccgwas_temp.R".format( group),"w") as file:
            file.write(rscript)

    script_run_r = "{} _{}_gwaslab_ccgwas_temp.R".format(r, group )
    
    try:
        log.write(" -Running CCGWAS from command line...", verbose=verbose)
        output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
        r_log+= output + "\n"
        #os.remove("_{}_{}_gwaslab_hyprcoloc_temp.R".format(study,locus))
    except subprocess.CalledProcessError as e:
        log.write(e.output)
        #os.remove("_{}_{}_gwaslab_hyprcoloc_temp.R".format(study,locus))
    log.write(" -Finishing CCGWAS for {}...".format(group),verbose=verbose)
    log.write("Finished Case-case GWAS using CCGWAS.", verbose=verbose)
    return output_prefix
