from typing import TYPE_CHECKING, List, Union, Optional, Tuple, Any, Dict
import pandas as pd
import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import norm
from gwaslab.info.g_Log import Log
from gwaslab.io.io_to_pickle import load_data_from_pickle
from gwaslab.g_Sumstats import Sumstats
from gwaslab.qc.qc_decorator import with_logging
import gc

@with_logging(
        start_to_msg="perform meta-analysis",
        finished_msg="meta-analysis successfully"
)
def meta_analyze(sumstats_list: List[Union[str, pd.DataFrame, Sumstats]],
                 random_effects: bool = False, 
                 match_allele: bool = True, 
                 log: Log = Log()) -> pd.DataFrame:
    
    ###########################################################################
    columns=["SNPID","CHR","POS","EA","NEA"]
    results_df = pd.DataFrame(columns=columns)

    log.write(" -Datasets:")

    for index,sumstats_path in enumerate(sumstats_list):
        if isinstance(sumstats_path, pd.DataFrame):
            log.write("  -Sumstats #{}: {} ".format(index, sumstats_path))
        elif isinstance(sumstats_path, Sumstats):
            log.write("  -Sumstats #{}: {} ".format(index, sumstats_path))
        else:
            log.write("  -Sumstats #{}: {} ".format(index, sumstats_path))
    
    
    # extract all variants information
    log.write(" -Iterating through {} datasets to determine variant list...".format(len(sumstats_list)))
    
    for index,sumstats_path in enumerate(sumstats_list):
        sumstats = get_sumstats(sumstats_path,usekeys=["SNPID","CHR","POS","EA","NEA"])
        new_rows = sumstats.loc[~sumstats["SNPID"].isin(results_df["SNPID"]),["SNPID","CHR","POS","EA","NEA"]]
        log.write("  -Sumstats #{}: {} new variants (out of {}) are being added to analysis...".format(index, len(new_rows),len(sumstats)))
        
        if len(new_rows)>0:
            if len(results_df) == 0:
                results_df = new_rows
            else:
                results_df = pd.concat([results_df, new_rows],ignore_index=True)
        del sumstats
        del new_rows
        gc.collect()
    
    ###########################################################################
    log.write(" -Initiating result DataFrame...")
    columns=["SNPID","CHR","POS","EA","NEA","_BETAW_SUM","_EA_N","_NEA_N","_BETA2W_SUM","_W_SUM","EAF","N","DIRECTION","BETA","SE","DOF"]
    results_df = results_df.set_index("SNPID")
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
    ###########################################################################
    
    log.write(" -Iterating through {} datasets to compute statistics for fixed-effect model...".format(len(sumstats_list)))
    for index,sumstats_path in enumerate(sumstats_list):
        to_use_sumstats = process_sumstats(sumstats_path, 
                                           results_df[["EA","NEA"]],
                                           index=index,
                                           match_allele=match_allele,)
        sumstats_index = to_use_sumstats.index
        results_df_not_in_sumstat_index = results_df.index[~results_df.index.isin(to_use_sumstats.index)]
        
        # N and DOF
        results_df.loc[sumstats_index, "N"]         += to_use_sumstats["N"]
        results_df.loc[sumstats_index, "DOF"]       += 1        
        
        # BEAT and SE
        results_df.loc[sumstats_index,"_BETA2W_SUM"] += to_use_sumstats["BETA"]**2 *(1/(to_use_sumstats["SE"]**2))
        results_df.loc[sumstats_index,"_BETAW_SUM"]  += to_use_sumstats["BETA"]*(1/(to_use_sumstats["SE"]**2))
        results_df.loc[sumstats_index,"_W_SUM"]      += 1/(to_use_sumstats["SE"]**2)  
        results_df.loc[sumstats_index,"_W2_SUM"]     += results_df.loc[sumstats_index,"_W_SUM"]**2  
        
        # EAF
        results_df.loc[sumstats_index,"_EA_N"] += to_use_sumstats["N"]*to_use_sumstats["EAF"]
        results_df.loc[sumstats_index,"_NEA_N"]  += to_use_sumstats["N"]*(1 - to_use_sumstats["EAF"])  

        # DIRECTION
        beta_index = to_use_sumstats[to_use_sumstats["BETA"]>0].index
        results_df.loc[beta_index, "DIRECTION"] += "+"
        beta_index = to_use_sumstats[to_use_sumstats["BETA"]==0].index
        results_df.loc[beta_index, "DIRECTION"] += "0"
        beta_index = to_use_sumstats[to_use_sumstats["BETA"]<0].index
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
    
    results_df["I2_HET"] = (results_df["Q"] - results_df["DOF"])/results_df["Q"]
    results_df.loc[results_df["I2_HET"]<0, "I2_HET"] = 0
    
    results_df=results_df.drop(columns=["_EA_N","_NEA_N"])
    gc.collect()

    ###########################################################################
    if random_effects==True:        
        log.write(" -Iterating through {} datasets to compute statistics for random-effects model...".format(len(sumstats_list)))
        results_df["_R2"] = (results_df["Q"] - results_df["DOF"])/(results_df["_W_SUM"] - (results_df["_W2_SUM"]/results_df["_W_SUM"]))
        results_df.loc[results_df["_R2"]<0, "_R2"] = 0
        variant_index_random = results_df[results_df["_R2"]>0].index

        results_df["_BETAW_SUM_R"] = 0.0  
        results_df["_W_SUM_R"] = 0.0
        results_df["BETA_RANDOM"] = results_df["BETA"]
        results_df["SE_RANDOM"] = results_df["SE"]

        for index,sumstats_path in enumerate(sumstats_list):
            to_use_sumstats = process_sumstats(sumstats_path, 
                                               results_df.loc[variant_index_random, ["EA","NEA"]], 
                                               index=index,
                                               match_allele=match_allele,
                                               extract_index=variant_index_random)
            
            sumstats_index = to_use_sumstats.index
            
            # BEAT and SE
            results_df.loc[sumstats_index,"_BETAW_SUM_R"]  += to_use_sumstats["BETA"]*(1/(to_use_sumstats["SE"]**2 + results_df.loc[sumstats_index,"_R2"]))
            results_df.loc[sumstats_index,"_W_SUM_R"]      += 1/(to_use_sumstats["SE"]**2 + results_df.loc[sumstats_index,"_R2"])
            
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
    results_df = results_df.drop(columns=["_BETAW_SUM","_BETA2W_SUM","_W_SUM","_R2","_W2_SUM"]).sort_values(by=["CHR","POS"])
    gc.collect()
    
    return results_df
     
def process_sumstats(
    sumstats_path: Union[str, pd.DataFrame, Sumstats, Tuple[Union[str, pd.DataFrame], Dict[str, Any]]],
    results_df: pd.DataFrame,
    index: int,
    extract_index: Optional[pd.Index] = None,
    match_allele: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    
    if extract_index is None:
        extract_index = results_df.index
        
    sumstats = get_sumstats(sumstats_path)
    
    to_use_sumstats = sumstats.loc[sumstats["SNPID"].isin(extract_index.values),["SNPID","EA","NEA","BETA","N","SE","EAF"]]
    
    if len(to_use_sumstats)>0:
        n_pre_dup = len(to_use_sumstats)
        log.write("  -Processing {} variants from sumstats #{}".format(len(to_use_sumstats), index))
        
        to_use_sumstats = to_use_sumstats.drop_duplicates(subset="SNPID").set_index("SNPID")
        n_post_dup = len(to_use_sumstats)
        
        if n_pre_dup - n_post_dup>0:
            log.write("  -Dropping {} duplicated variants from sumstats #{}".format(n_pre_dup - n_post_dup, index))
        
        if match_allele==True:
            sumstats_index = to_use_sumstats.index
            # drop not matched
            is_match = (to_use_sumstats.loc[sumstats_index,"EA"] == results_df.loc[sumstats_index, "EA"] )&(to_use_sumstats.loc[sumstats_index,"NEA"] == results_df.loc[sumstats_index, "NEA"])
            is_flip = (to_use_sumstats.loc[sumstats_index,"EA"] == results_df.loc[sumstats_index, "NEA"])&( to_use_sumstats.loc[sumstats_index,"NEA"] == results_df.loc[sumstats_index, "EA"])
            is_flip = is_flip | ((to_use_sumstats.loc[sumstats_index,"NEA"] == results_df.loc[sumstats_index, "EA"])&( to_use_sumstats.loc[sumstats_index,"EA"] == results_df.loc[sumstats_index, "NEA"])) 
            is_to_use = is_match|is_flip
            
            if sum(~is_to_use)>0:
                log.write("  -Dropping {} variants with unmatched alleles from sumstats #{}".format(sum(~is_to_use), index))
            
            to_use_sumstats.loc[is_flip[is_flip].index, "BETA"] =  -to_use_sumstats.loc[is_flip[is_flip].index, "BETA"] 
            to_use_sumstats.loc[is_flip[is_flip].index, "EAF"]  = 1-to_use_sumstats.loc[is_flip[is_flip].index, "EAF"] 
            to_use_sumstats = to_use_sumstats.loc[is_to_use[is_to_use].index,:]  
    
    gc.collect()

    return to_use_sumstats

def get_sumstats(
    input_path: Union[str, pd.DataFrame, Sumstats, Tuple[Union[str, pd.DataFrame], Dict[str, Any]]],
    usekeys: Optional[List[str]] = None
) -> pd.DataFrame:
    if isinstance(input_path, tuple):
        path = input_path[0]
        path_kwargs = input_path[1]
    else:
        path = input_path
        path_kwargs={}
        
    if isinstance(path, pd.DataFrame):
        sumstats = Sumstats(path,fmt="auto",verbose=False,usekeys=usekeys,**path_kwargs).data
    elif isinstance(path, Sumstats):
        sumstats = path.data
        if usekeys is not None:
            sumstats = sumstats[usekeys]
    elif path[-6:] == "pickle":
        sumstats = load_data_from_pickle(path)
        if usekeys is not None:
            sumstats = sumstats[usekeys]
    else:
        sumstats = Sumstats(path,fmt="auto",verbose=False,usekeys=usekeys,**path_kwargs).data
    return sumstats


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
        to_use_sumstats = sumstats_multi.loc[~sumstats_multi["BETA_{}".format(i+1)].isna(),single_study_cols].drop_duplicates(subset="_INDEX").set_index("_INDEX")


        sumstats_index = to_use_sumstats.index

        results_df_not_in_sumstat_index = results_df.index[~results_df.index.isin(to_use_sumstats.index)]
        
        # N and DOF
        results_df.loc[sumstats_index, "N"]         += to_use_sumstats[n].fillna(0)
        results_df.loc[sumstats_index, "DOF"]       += 1        
        
        # BEAT and SE
        results_df.loc[sumstats_index,"_BETA2W_SUM"] += to_use_sumstats[beta]**2 *(1/(to_use_sumstats[se]**2))
        results_df.loc[sumstats_index,"_BETAW_SUM"]  += to_use_sumstats[beta]*(1/(to_use_sumstats[se]**2))
        results_df.loc[sumstats_index,"_W_SUM"]      += 1/(to_use_sumstats[se]**2)  
        results_df.loc[sumstats_index,"_W2_SUM"]     += results_df.loc[sumstats_index,"_W_SUM"]**2  
        
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
    
    results_df["I2"] = (results_df["Q"] - results_df["DOF"])/results_df["Q"]
    results_df.loc[results_df["I2"]<0, "I2"] = 0
    
    results_df=results_df.drop(columns=["_EA_N","_NEA_N"])
    gc.collect()

    ###########################################################################
    if random_effects==True:        
        log.write(" -Iterating through {} datasets to compute statistics for random-effects model...".format(nstudy))
        results_df["_R2"] = (results_df["Q"] - results_df["DOF"])/(results_df["_W_SUM"] - (results_df["_W2_SUM"]/results_df["_W_SUM"]))
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
            to_use_sumstats = sumstats_multi.loc[~sumstats_multi["BETA_{}".format(i+1)].isna(),single_study_cols].drop_duplicates(subset="_INDEX").set_index("_INDEX")
            sumstats_index = to_use_sumstats.index
            
            # BEAT and SE
            results_df.loc[sumstats_index,"_BETAW_SUM_R"]  += to_use_sumstats[beta]*(1/(to_use_sumstats[se]**2 + results_df.loc[sumstats_index,"_R2"]))
            results_df.loc[sumstats_index,"_W_SUM_R"]      += 1/(to_use_sumstats[se]**2 + results_df.loc[sumstats_index,"_R2"])
            
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

    results_df = Sumstats(results_df, fmt="gwaslab", other = other_cols)
    
    return results_df

def _init_result_df(sumstats: pd.DataFrame) -> pd.DataFrame:

    results_df = sumstats[["_INDEX","SNPID","CHR","POS","EA","NEA"]]
    results_df = results_df.drop_duplicates(subset="_INDEX").set_index("_INDEX")

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