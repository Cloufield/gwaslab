import pandas as pd
import numpy as np
from gwaslab.g_Log import Log

# Calculate PIP based on approximate Bayesian factor (ABF)
# Wakefield, J. A bayesian measure of the probability of false discovery in genetic epidemiology studies. Am J Hum Genet 81, 208–227 (2007).


def calc_abf(insumstats,omega=0.04,log=Log(),verbose=True):
    log.write("Start to calculate approximate Bayesian factor for each {} variant...".format(str(len(insumstats))),verbose=verbose)
    se = insumstats["SE"]
    beta = insumstats["BETA"]
    insumstats = insumstats.copy()
    insumstats.loc[:, "log_ABF"] = 1/2*np.log(((se*se)/(se*se+omega))) + (omega*beta*beta)/(2*se*se*(se*se+omega))
    return insumstats

def calc_PIP(insumstats,log=Log(),verbose=True):
    # Calculate the logarithmic sum of each ABF to find the logarithm of total_abf
    log_total_abf = np.log(np.sum(np.exp(insumstats["log_ABF"] - np.max(insumstats["log_ABF"])))) + np.max(insumstats["log_ABF"])
    insumstats = insumstats.copy()
    log.write("Start to calculate PIP Bayesian factor for each {} variant...".format(str(len(insumstats))),verbose=verbose)
    # Calculate PIP on a logarithmic scale by subtracting log_total_abf from each log_abf
    insumstats.loc[:, "log_PIP"] = insumstats['log_ABF'] - log_total_abf
    # Convert PIP on logarithmic scale to exponential and back to normal scale
    insumstats.loc[:, "PIP"] = np.exp(insumstats['log_PIP'])
    return insumstats

def abf_finemapping(insumstats):
    insumstats = calc_abf(insumstats)
    insumstats = calc_PIP(insumstats)
    return insumstats

def make_cs(insumstats,threshold=0.95,log=Log(),verbose=True):
    if "PIP" not in list(insumstats.columns):
        insumstats = abf_finemapping(insumstats)
    insumstats = insumstats.sort_values(by="PIP",ascending=False)
    pip_sum = 0
    cs = pd.DataFrame()
    for index, row in insumstats.iterrows():
        cs = pd.concat([cs,pd.DataFrame(row).T])
        pip_sum += row["PIP"]
        if pip_sum > threshold:
            break
    log.write("Finish constructing {}% credible set with {} variant...".format(str(len(threshold)),str(len(cs))),verbose=verbose)
    return cs
