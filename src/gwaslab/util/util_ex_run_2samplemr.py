import subprocess
import os
import gc
import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Optional, List, Tuple, Any, Union

if TYPE_CHECKING:
    from gwaslab.g_SumstatsPair import SumstatsPair

from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_susie_version
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2


def _run_two_sample_mr(
    sumstatspair_object: 'SumstatsPair',
    r: str,
    out: str = "./",
    clump: bool = False,
    f_check: float = 10,
    exposure1: str = "Trait1",
    outcome2: str = "Trait2",
    n1: Optional[int] = None,
    n2: Optional[int] = None,
    binary1: bool = False,
    cck1: Optional[Tuple[int, int, float]] = None,
    cck2: Optional[Tuple[int, int, float]] = None,
    ncase1: Optional[int] = None,
    ncontrol1: Optional[int] = None,
    prevalence1: Optional[float] = None,
    binary2: bool = False,
    ncase2: Optional[int] = None,
    ncontrol2: Optional[int] = None,
    prevalence2: Optional[float] = None,
    methods: Optional[List[str]] = None,
    log: Log = Log()
) -> None:

    log.write(" Start to run MR using twosampleMR from command line:")
    if methods is None:
        methods = ["mr_ivw","mr_simple_mode","mr_weighted_median","mr_egger_regression","mr_ivw_mre", "mr_weighted_mode"]
        methods_string = '"{}"'.format('","'.join(methods))
    
    if cck1 is not None:
        log.write(" - ncase1, ncontrol1, prevalence1:{}".format(cck1))
        binary1 = True
        ncase1 = cck1[0]
        ncontrol1 = cck1[1]
        prevalence1 =  cck1[2]
        n1 = ncase1 + ncontrol1
    if cck2 is not None:
        log.write(" - ncase2, ncontrol2, prevalence2:{}".format(cck2))
        binary2 = True
        ncase2 = cck2[0]
        ncontrol2 = cck2[1]
        prevalence2 =  cck2[2]
        n2 = ncase2 + ncontrol2

    if clump==True:
        sumstatspair = sumstatspair_object.clumps["clumps"]
    else: 
        sumstatspair = sumstatspair_object.data
    
    if n1 is not None:
        sumstatspair["N_1"] = n1
    if n2 is not None:
        sumstatspair["N_2"] = n2
    
    sumstatspair = _filter_by_f(sumstatspair, f_check, n1, binary1, ncase1, prevalence1)

    log = _checking_r_version(r, log)
    #log = _check_susie_version(r,log)       

    r_log=""

    cols_for_trait1, cols_for_trait2 = _sort_columns_to_load(sumstatspair)
    cols_for_trait1_script = _cols_list_to_r_script(cols_for_trait1)
    cols_for_trait2_script = _cols_list_to_r_script(cols_for_trait2)
    
    # Clumping       

    prefix = "{exposure}_{outcome}_{memory_id}".format(exposure = exposure1, outcome= outcome2, memory_id = id(sumstatspair))
    prefix = "{}{}".format(out.rstrip('/') + "/",prefix)
    temp_sumstats_path = "{out}twosample_mr_{exposure}_{outcome}_{memory_id}.csv.gz".format(out=out.rstrip('/') + "/",
                                                                                               exposure = exposure1, 
                                                                                               outcome= outcome2, 
                                                                                               memory_id = id(sumstatspair))
    if len(sumstatspair)>0:
        sumstatspair.to_csv(temp_sumstats_path ,index=None)
    else:
        return 0
    ###
    calculate_r_script = ""
    
    if binary1==True:
        calculate_r_script+= _make_script_for_calculating_r("exposure", ncase1, ncontrol1, prevalence1)
    else:
        calculate_r_script+= _make_script_for_calculating_r_quant("exposure")
    
    if binary2==True:
        calculate_r_script+= _make_script_for_calculating_r("outcome", ncase2, ncontrol2, prevalence2)
    else:
        calculate_r_script+= _make_script_for_calculating_r_quant("outcome")
    
    # create scripts
    directionality_test_script='''
    results_directionality <- directionality_test(harmonized_data)
    write.csv(results_directionality, "{prefix}.directionality", row.names = FALSE)
    '''.format(prefix = prefix)
    
    # Two Sample MR
    # Tests
    ## Pleiotropy
    ## Heterogeneity
    ## Reverse causality
    rscript='''
    library(TwoSampleMR)
    sumstats <- read.csv("{temp_sumstats_path}")

    {pheno1}
    {pheno2}

    exp_raw <-sumstats[,c({cols_for_trait1})]

    out_raw <-sumstats[,c({cols_for_trait2})] 

    exp_dat <- format_data( exp_raw,
                  type = "exposure",
                  snp_col = "SNPID",
                  beta_col = "BETA_1",
                  se_col = "SE_1",
                  effect_allele_col = "EA",
                  other_allele_col = "NEA",
                  eaf_col = "EAF_1",
                  pval_col = "P_1",
                  phenotype_col = "PHENO_1",
                  samplesize_col= "N_1"
                 )
    
    out_dat <- format_data( out_raw,
                  type = "outcome",
                  snp_col = "SNPID",
                  beta_col = "BETA_2",
                  se_col = "SE_2",
                  effect_allele_col = "EA",
                  other_allele_col = "NEA",
                  eaf_col = "EAF_2",
                  pval_col = "P_2",
                  phenotype_col = "PHENO_2",
                  samplesize_col= "N_2"
                  )
    
    harmonized_data <- harmonise_data(exp_dat,out_dat,action=1)
    
    results_mr <- mr(harmonized_data, method_list = c({methods_string}))
    write.csv(results_mr, "{prefix}.mr", row.names = FALSE)

    results_heterogeneity <- mr_heterogeneity(harmonized_data)
    write.csv(results_heterogeneity, "{prefix}.heterogeneity", row.names = FALSE)

    results_pleiotropy <- mr_pleiotropy_test(harmonized_data)
    write.csv(results_pleiotropy, "{prefix}.pleiotropy", row.names = FALSE)
    
    results_single <- mr_singlesnp(harmonized_data)
    write.csv(results_single, "{prefix}.singlesnp", row.names = FALSE)
    
    results_loo <- mr_leaveoneout(harmonized_data)
    write.csv(results_loo, "{prefix}.leaveoneout", row.names = FALSE)
    
    {calculate_r}
    {directionality_test}

    '''.format(
        temp_sumstats_path = temp_sumstats_path,
        pheno1 = 'sumstats$PHENO_1 <- "{}"'.format(exposure1) if exposure1 is not None else "", 
        pheno2 = 'sumstats$PHENO_2 <- "{}"'.format(outcome2) if outcome2 is not None else "", 
        cols_for_trait1 = cols_for_trait1_script,
        cols_for_trait2 = cols_for_trait2_script,
        methods_string=methods_string,
        prefix=prefix,
        calculate_r = calculate_r_script,
        directionality_test = directionality_test_script
    )
        
    temp_r_script_path = "{}_{}_{}_{}_gwaslab_2smr_temp.R".format(out.rstrip('/') + "/",
                                                                  exposure1,
                                                                  outcome2,
                                                                  id(sumstatspair))
    with open(temp_r_script_path,"w") as file:
            file.write(rscript)

    script_run_r = "{} {}".format(r, temp_r_script_path)
        
    try:
        log.write(" Running TwoSampleMR from command line...")
        output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
        #plink_process = subprocess.Popen("exec "+script_run_r, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
        #output1,output2 = plink_process.communicate()
        #output= output1 + output2+ "\n"
        #plink_process.kill()
        log.write(output)
        r_log+= output + "\n"
        os.remove(temp_r_script_path)
        try:
            for suffix in ["mr","pleiotropy","heterogeneity","singlesnp","leaveoneout","directionality"]:
                sumstatspair_object.mr[suffix] = pd.read_csv("{}.{}".format(prefix,suffix))
        except:
            pass     
        sumstatspair_object.mr["r_log"] = r_log
    except subprocess.CalledProcessError as e:
        log.write(" Error!")
        log.write(rscript)
        log.write(e.output)
        os.remove(temp_r_script_path)
    log.write(" Finished running MR using twosampleMR from command line.")



def _sort_columns_to_load(sumstatspair: pd.DataFrame) -> Tuple[List[str], List[str]]:
    cols_for_trait1=["SNPID","CHR","POS","EA","NEA","PHENO_1","N_1"]
    cols_for_trait2=["SNPID","CHR","POS","EA","NEA","PHENO_2","N_2"]
    for i in ["EAF","BETA","SE","P"]:
        if i+"_1" not in cols_for_trait1 and i+"_1" in sumstatspair.columns:
            cols_for_trait1.append(i+"_1")
        if i+"_2" not in cols_for_trait2 and i+"_2" in sumstatspair.columns:
            cols_for_trait2.append(i+"_2")
    return cols_for_trait1, cols_for_trait2




def _cols_list_to_r_script(cols_for_trait1: List[str]) -> str:
    script = '"{}"'.format('","'.join(cols_for_trait1))
    return script

def _make_script_for_calculating_r(
    exposure_or_outcome: str, 
    ncase: int, 
    ncontrol: int, 
    prevalence: float
) -> str:
        
        script = """
        harmonized_data$"r.{exposure_or_outcome}" <- get_r_from_lor(  harmonized_data$"beta.{exposure_or_outcome}",
                                                        harmonized_data$"eaf.{exposure_or_outcome}",
                                                        {ncase},
                                                        {ncontrol},
                                                        {prevalence},
                                                        model = "logit",
                                                        correction = FALSE
                                                        )
        """.format(
             exposure_or_outcome = exposure_or_outcome,
             ncase = ncase,
             ncontrol = ncontrol,
             prevalence = prevalence
        )
        return script


def _make_script_for_calculating_r_quant(exposure_or_outcome: str) -> str:
        script = """
        harmonized_data$"r.{exposure_or_outcome}" <- get_r_from_bsen(  harmonized_data$"beta.{exposure_or_outcome}",
                                                        harmonized_data$"se.{exposure_or_outcome}",
                                                        harmonized_data$"samplesize.{exposure_or_outcome}"
                                                        )
        """.format(
             exposure_or_outcome = exposure_or_outcome
        )
        return script


def _filter_by_f(
    sumstatspair: pd.DataFrame, 
    f_check: float, 
    n1: Optional[int], 
    binary1: Optional[bool] = None, 
    ncase1: Optional[int] = None, 
    ncontrol1: Optional[int] = None, 
    prevalence1: Optional[float] = None, 
    log: Log = Log()
) -> pd.DataFrame:

    if binary1==True:
        sumstatspair = _get_per_snp_r2(sumstatspair,
            beta="BETA_1",
            af="EAF_1",
            n = "N_1", 
            mode="b",
            se="SE_1",
            vary=1,
            ncase=ncase1, 
            ncontrol=ncontrol1, 
            prevalence=prevalence1, 
            k=1)
    else:
        sumstatspair = _get_per_snp_r2(sumstatspair,
            beta="BETA_1",
            af="EAF_1",
            n = "N_1", 
            mode="q",
            se="SE_1",
            vary=1,
            k=1)

    log.write("Filtered out {} variants with F < {}".format(len(sumstatspair) - sum(sumstatspair["F"]>f_check),f_check))
    sumstatspair = sumstatspair.loc[sumstatspair["F"]>f_check,:]

    return sumstatspair
