from typing import TYPE_CHECKING, Union, Optional
import pandas as pd
import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import norm
from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats
from gwaslab.qc.qc_decorator import with_logging
import gc

############################################################################################################################################################################

@with_logging(
        start_to_msg="perform meta-analysis",
        finished_msg="meta-analysis successfully"
)
def meta_analyze_multi(
    sumstats_multi: pd.DataFrame,
    random_effects: bool = False,
    nstudy: int = 1,
    match_allele: bool = True,
    log: Log = Log()
) -> Sumstats:
    ###########################################################################
    log.write(" -Initiating result DataFrame...")
    
    # Deduplicate BEFORE creating _INDEX to prevent same variant with different _INDEX
    # Always use CHR:POS:EA:NEA to determine duplicates (more reliable than SNPID)
    # This is much faster than groupby loops - uses vectorized pandas operations
    n_before = len(sumstats_multi)
    variant_cols = ["CHR", "POS", "EA", "NEA"]
    if all(col in sumstats_multi.columns for col in variant_cols):
        sumstats_multi = sumstats_multi.drop_duplicates(subset=variant_cols, keep="first")
    else:
        log.write(" -Warning: Missing required columns (CHR, POS, EA, NEA) for deduplication", verbose=True)
    
    n_after = len(sumstats_multi)
    if n_before > n_after:
        log.write(" -Removed {} duplicate variants before meta-analysis".format(n_before - n_after), verbose=True)
    
    # Now create _INDEX after deduplication - each variant gets unique _INDEX
    sumstats_multi["_INDEX"] = range(len(sumstats_multi))
    results_df = _init_result_df(sumstats_multi)
    ##########################################################################


    log.write(" -Iterating through {} datasets to compute statistics for fixed-effect model...".format(nstudy))
    for i in range(nstudy):
        n="N_{}".format(i+1)
        beta="BETA_{}".format(i+1)
        se="SE_{}".format(i+1)
        eaf="EAF_{}".format(i+1)
        single_study_cols=[n,beta,se,eaf,"SNPID","_INDEX"]
        # Filter for valid studies: BETA not null, SE not null and > 0, N and EAF not null
        valid_mask = (
            ~sumstats_multi["BETA_{}".format(i+1)].isna() &
            ~sumstats_multi[se].isna() &
            (sumstats_multi[se] > 0) &
            ~sumstats_multi[n].isna() &
            ~sumstats_multi[eaf].isna()
        )
        to_use_sumstats = sumstats_multi.loc[valid_mask, single_study_cols].drop_duplicates(subset="_INDEX").set_index("_INDEX")


        sumstats_index = to_use_sumstats.index

        results_df_not_in_sumstat_index = results_df.index[~results_df.index.isin(to_use_sumstats.index)]
        
        # N and DOF
        results_df.loc[sumstats_index, "N"]         += to_use_sumstats[n].fillna(0)
        results_df.loc[sumstats_index, "DOF"]       += 1        
        
        # BEAT and SE
        # Calculate weights (inverse variance) - SE is already validated to be > 0
        weights = 1 / (to_use_sumstats[se]**2)
        results_df.loc[sumstats_index,"_BETA2W_SUM"] += to_use_sumstats[beta]**2 * weights
        results_df.loc[sumstats_index,"_BETAW_SUM"]  += to_use_sumstats[beta] * weights
        results_df.loc[sumstats_index,"_W_SUM"]      += weights
        # _W2_SUM should be sum of squared individual weights, not square of cumulative sum
        # Correct: W2 = sum(w_i^2) where w_i = 1/SE_i^2
        results_df.loc[sumstats_index,"_W2_SUM"]     += weights ** 2  
        
        # EAF
        results_df.loc[sumstats_index,"_EA_N"] += to_use_sumstats[n]*to_use_sumstats[eaf]
        results_df.loc[sumstats_index,"_NEA_N"]  += to_use_sumstats[n]*(1 - to_use_sumstats[eaf])  

        # DIRECTION
        beta_index = to_use_sumstats[to_use_sumstats[beta]>0].index
        results_df.loc[beta_index, "DIRECTION"] += "+"
        beta_index = to_use_sumstats[to_use_sumstats[beta]==0].index
        results_df.loc[beta_index, "DIRECTION"] += "0"
        beta_index = to_use_sumstats[to_use_sumstats[beta]<0].index
        results_df.loc[beta_index, "DIRECTION"] += "-"
        results_df.loc[results_df_not_in_sumstat_index, "DIRECTION"] += "?"
        
        del to_use_sumstats
        gc.collect()

    ############################################################################## 
    # fixed - effect statistics
    results_df["BETA"] = results_df["_BETAW_SUM"] / results_df["_W_SUM"]
    results_df["EAF"] = results_df["_EA_N"] / (results_df["_EA_N"] + results_df["_NEA_N"])
    results_df["SE"] = np.sqrt(1/results_df["_W_SUM"])
    results_df["Z"] = results_df["BETA"] / results_df["SE"]
    results_df["P"] = norm.sf(abs(results_df["Z"]))*2
    results_df["Q"] = results_df["_BETA2W_SUM"] - (results_df["_BETAW_SUM"]**2 / results_df["_W_SUM"])
    
    for dof in results_df["DOF"].unique():
        results_df_dof_index = results_df["DOF"] == dof
        results_df.loc[results_df_dof_index,"P_HET"] = chi2.sf(results_df.loc[results_df_dof_index, "Q"].values,dof)
        gc.collect()
    
    # I2: I-squared statistic for heterogeneity
    # I2 = max(0, (Q - df) / Q) when Q > 0, else 0
    # Guard against division by zero when Q == 0
    results_df["I2"] = np.where(
        results_df["Q"] <= 0,
        0.0,
        (results_df["Q"] - results_df["DOF"]) / results_df["Q"]
    )
    results_df.loc[results_df["I2"]<0, "I2"] = 0
    
    results_df=results_df.drop(columns=["_EA_N","_NEA_N"])
    gc.collect()

    ###########################################################################
    if random_effects==True:        
        log.write(" -Iterating through {} datasets to compute statistics for random-effects model...".format(nstudy))
        # Calculate tau^2 (between-study variance) using DerSimonian-Laird estimator
        # tau^2 = max(0, (Q - df) / C) where C = W - W2/W
        # Guard against division by zero or negative denominator
        C = results_df["_W_SUM"] - (results_df["_W2_SUM"] / results_df["_W_SUM"])
        results_df["_R2"] = np.where(
            (C <= 0) | (results_df["Q"] <= results_df["DOF"]),
            0.0,
            (results_df["Q"] - results_df["DOF"]) / C
        )
        results_df.loc[results_df["_R2"]<0, "_R2"] = 0
        variant_index_random = results_df[results_df["_R2"]>0].index

        results_df["_BETAW_SUM_R"] = 0.0  
        results_df["_W_SUM_R"] = 0.0
        results_df["BETA_RANDOM"] = results_df["BETA"]
        results_df["SE_RANDOM"] = results_df["SE"]

        for i in range(nstudy):
            n="N_{}".format(i+1)
            beta="BETA_{}".format(i+1)
            se="SE_{}".format(i+1)
            eaf="EAF_{}".format(i+1)
            single_study_cols=[n,beta,se,eaf,"SNPID","_INDEX"]
            # Filter for valid studies: BETA not null, SE not null and > 0
            valid_mask = (
                ~sumstats_multi["BETA_{}".format(i+1)].isna() &
                ~sumstats_multi[se].isna() &
                (sumstats_multi[se] > 0)
            )
            to_use_sumstats = sumstats_multi.loc[valid_mask, single_study_cols].drop_duplicates(subset="_INDEX").set_index("_INDEX")
            sumstats_index = to_use_sumstats.index
            
            # BEAT and SE
            # Calculate random effects weights: 1 / (SE^2 + tau^2)
            # SE is already validated to be > 0
            weights_r = 1 / (to_use_sumstats[se]**2 + results_df.loc[sumstats_index,"_R2"])
            results_df.loc[sumstats_index,"_BETAW_SUM_R"]  += to_use_sumstats[beta] * weights_r
            results_df.loc[sumstats_index,"_W_SUM_R"]      += weights_r
            
            del to_use_sumstats
            del sumstats_index
            gc.collect()
            
        results_df.loc[variant_index_random,"BETA_RANDOM"] = results_df.loc[variant_index_random,"_BETAW_SUM_R"] / results_df.loc[variant_index_random,"_W_SUM_R"]
        results_df.loc[variant_index_random,"SE_RANDOM"] = np.sqrt(1/results_df.loc[variant_index_random,"_W_SUM_R"])
        results_df["Z_RANDOM"] = results_df["BETA_RANDOM"] / results_df["SE_RANDOM"]
        results_df["P_RANDOM"] = norm.sf(abs(results_df["Z_RANDOM"]))*2
        results_df = results_df.drop(columns=["_BETAW_SUM_R","_W_SUM_R"])

        gc.collect()
    ###########################################################################
    results_df = results_df.drop(columns=["_BETAW_SUM","_BETA2W_SUM","_W_SUM","_R2","_W2_SUM"]).sort_values(by=["CHR","POS"]).reset_index()
    gc.collect()

    if random_effects==True:
        other_cols = ["BETA_RANDOM","SE_RANDOM","Z_RANDOM","P_RANDOM"]
    else:
        other_cols = []
    
    results_df = results_df.drop(columns=["_INDEX"])    

    # Return as Sumstats object
    result_sumstats = Sumstats(results_df, fmt="gwaslab", other=other_cols)
    
    return result_sumstats

def _init_result_df(sumstats: pd.DataFrame) -> pd.DataFrame:

    results_df = sumstats[["_INDEX","SNPID","CHR","POS","EA","NEA"]]
    # First deduplicate by _INDEX to handle any duplicate _INDEX values
    results_df = results_df.drop_duplicates(subset="_INDEX")
    
    # Then deduplicate by variant identity (CHR:POS:EA:NEA) to prevent
    # the same variant with different _INDEX from being counted multiple times
    # Always use CHR:POS:EA:NEA (more reliable than SNPID)
    variant_cols = ["CHR", "POS", "EA", "NEA"]
    if all(col in results_df.columns for col in variant_cols):
        results_df = results_df.drop_duplicates(subset=variant_cols, keep="first")
    
    results_df = results_df.set_index("_INDEX")

    results_df["N"] = 0 
    results_df["_BETAW_SUM"] = 0.0  
    results_df["_BETA2W_SUM"] = 0.0  
    results_df["_W_SUM"] = 0.0 
    results_df["_W2_SUM"] = 0.0
    results_df["_EA_N"] = 0.0
    results_df["_NEA_N"] = 0.0 
    results_df["N"] = 0 
    results_df["DIRECTION"] = ""
    results_df["BETA"] = 0.0 
    results_df["SE"] = 0.0 
    results_df["DOF"] = -1
    results_df["_R2"] = 0

    dtype_dict ={
        "_BETAW_SUM":"float64",
        "_EA_N":"float64",
        "_NEA_N":"float64",
        "_BETA2W_SUM":"float64",
        "_W_SUM":"float64",
        "BETA":"float64",
        "SE":"float64",
        "N":"Int64",
        "DOF":"Int64"
    }
    results_df=results_df.astype(dtype_dict)
    return results_df