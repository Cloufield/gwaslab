import pandas as pd
import numpy as np
from gwaslab.g_Log import Log

# Calculate PIP based on approximate Bayesian factor (ABF)
# Wakefield, J. A bayesian measure of the probability of false discovery in genetic epidemiology studies. Am J Hum Genet 81, 208–227 (2007).


def calc_abf(df,omega=0.04):
    se = df["SE"]
    beta = df["BETA"]
    log_ABF = 1/2*np.log(((se*se)/(se*se+omega))) + (omega*beta*beta)/(2*se*se*(se*se+omega))
    df["log_ABF"] = 1/2*np.log(((se*se)/(se*se+omega))) + (omega*beta*beta)/(2*se*se*(se*se+omega))
    return df

def calc_PIP(df):
    # Calculate the logarithmic sum of each ABF to find the logarithm of total_abf
    log_total_abf = np.log(np.sum(np.exp(df["log_ABF"] - np.max(df["log_ABF"])))) + np.max(df["log_ABF"])
    # Calculate PIP on a logarithmic scale by subtracting log_total_abf from each log_abf
    df['log_PIP'] = df['log_ABF'] - log_total_abf
    # Convert PIP on logarithmic scale to exponential and back to normal scale
    df['PIP'] = np.exp(df['log_PIP'])
    return df

def abf_finemapping(df):
    df = calc_abf(df)
    df = calc_PIP(df)
    return df