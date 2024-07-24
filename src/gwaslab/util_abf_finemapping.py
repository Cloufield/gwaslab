import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.util_in_filter_value import _get_flanking_by_chrpos
from gwaslab.util_in_filter_value import _get_flanking_by_id

# Calculate PIP based on approximate Bayesian factor (ABF)
# Wakefield, J. A bayesian measure of the probability of false discovery in genetic epidemiology studies. Am J Hum Genet 81, 208–227 (2007).


def calc_abf(insumstats,w=0.2,log=Log(),verbose=True,**kwargs):



    log.write("Start to calculate approximate Bayesian factor for {} variants".format(len(insumstats)),verbose=verbose)
    log.write(" - Reference: akefield, J. A bayesian measure of the probability of false discovery in genetic epidemiology studies. Am J Hum Genet 81, 208–227 (2007).",verbose=verbose)
    log.write(" - Priors for the standard deviation W of the effect size parameter β : {} ".format(w),verbose=verbose)
    # binary -> w=0.2
    # quant  -> w=0.15
    omega = w**2
    se = insumstats["SE"]
    v = se**2
    r = omega / (omega+v)
    beta = insumstats["BETA"]
    z = beta/se
    insumstats = insumstats.copy()

    # (6) ABF -> reciprocal
    insumstats.loc[:, "log_ABF"] = 1/2* (np.log(1-r) + (r * z**2))
    
    return insumstats

def calc_PIP(insumstats,log=Log(),verbose=True,**kwargs):
    # Calculate the logarithmic sum of each ABF to find the logarithm of total_abf
    log_total_abf = np.log(np.sum(np.exp(insumstats["log_ABF"] - np.max(insumstats["log_ABF"])))) + np.max(insumstats["log_ABF"])
    insumstats = insumstats.copy()
    log.write("Start to calculate PIP for {} variants".format(len(insumstats)),verbose=verbose)
    # Calculate PIP on a logarithmic scale by subtracting log_total_abf from each log_abf
    insumstats.loc[:, "log_PIP"] = insumstats['log_ABF'] - log_total_abf
    # Convert PIP on logarithmic scale to exponential and back to normal scale
    insumstats.loc[:, "PIP"] = np.exp(insumstats['log_PIP'])
    return insumstats

def abf_finemapping(insumstats,region=None,chrpos=None,snpid=None, log=Log(),**kwargs):

    if region is not None:
        region_data = insumstats[(insumstats["CHR"] == region[0]) & (insumstats["POS"] >= region[1]) & (insumstats["POS"] <= region[2])]
    elif chrpos is not None:
        region_data = _get_flanking_by_chrpos(insumstats, chrpos=chrpos,**kwargs)
    elif snpid is not None:
        region_data = _get_flanking_by_id(insumstats, snpid=snpid,**kwargs)

    region_data = calc_abf(region_data,log=log,**kwargs)
    region_data = calc_PIP(region_data,log=log,**kwargs)
    return region_data

def make_cs(insumstats,threshold=0.95,log=Log(),verbose=True):
    insumstats = insumstats.sort_values(by="PIP",ascending=False)
    pip_sum = 0
    cs = pd.DataFrame()
    for index, row in insumstats.iterrows():
        cs = pd.concat([cs,pd.DataFrame(row).T])
        pip_sum += row["PIP"]
        if pip_sum > threshold:
            break
    log.write("Finished constructing a {}% credible set with {} variant(s)".format(str(threshold * 100),str(len(cs))),verbose=verbose)
    return cs
