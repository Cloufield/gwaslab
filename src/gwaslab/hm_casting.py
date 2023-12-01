import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from pandas.api.types import CategoricalDtype


def _merge_mold_with_sumstats(mold, sumstats, log=Log()):
    
    cols_to_drop = []
    for i in sumstats.columns:
        if i in ["SNPID","rsID"]:
            cols_to_drop.append(i)
    
    if len(cols_to_drop)>0:
        log.write("Dropping old IDs:{}".format(cols_to_drop))
        sumstats = sumstats.drop(columns=cols_to_drop)
    
    mold_sumstats = pd.merge(mold, sumstats, on=["CHR","POS"], how="inner",suffixes=("_MOLD",""))
    log.write("After merging by CHR and POS:{}".format(len(mold_sumstats)))

    mold_sumstats = _keep_variants_with_same_allele_set(mold_sumstats)
    log.write("Matched variants:{}".format(len(mold_sumstats)))
    
    return mold_sumstats

def _keep_variants_with_same_allele_set(sumstats):

    all_alleles = set(list(sumstats["EA"].unique())+list(sumstats["NEA"].unique())+list(sumstats["EA_MOLD"].unique())+list(sumstats["NEA_MOLD"].unique()))
    allele_type = CategoricalDtype(categories=all_alleles, ordered=False)
    sumstats.loc[:, ["EA","EA_MOLD","NEA","NEA_MOLD"]] = sumstats.loc[:, ["EA","EA_MOLD","NEA","NEA_MOLD"]].astype(allele_type)
    
    is_perfect_match = (sumstats["EA"] == sumstats["EA_MOLD"]) & (sumstats["NEA"] == sumstats["NEA_MOLD"])
    is_flipped_match = (sumstats["EA"] == sumstats["NEA_MOLD"]) & (sumstats["NEA"] == sumstats["EA_MOLD"])
    is_allele_set_match = is_flipped_match | is_perfect_match
    
    return sumstats.loc[is_allele_set_match,:]

def _align_with_mold():
    pass