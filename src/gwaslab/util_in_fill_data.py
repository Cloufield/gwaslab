import pandas as pd
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
from scipy import stats
from gwaslab.g_Log import Log
import gc
#from gwaslab.qc_fix_sumstats import sortcolumn
from gwaslab.g_version import _get_version
from gwaslab.qc_check_datatype import check_datatype


def filldata( 
    insumstats,
    to_fill=None,
    df=None,
    overwrite=False,
    verbose=True,
    only_sig=False,
    sig_level=5e-8,
    extreme=False,
    log = Log()
    ):
    
    # if a string is passed to to_fill, convert it to list
    if type(to_fill) is str:
        to_fill = [to_fill]
    sumstats = insumstats.copy()
    log.write("Start filling data using existing columns...{}".format(_get_version()), verbose=verbose)
    
    check_datatype(sumstats,verbose=verbose,log=log)
    
# check dupication ##############################################################################################
    skip_cols=[]
    log.write(" -Overwrite mode: ",overwrite, verbose=verbose)
    if overwrite is False:
        for i in to_fill:
            if i in sumstats.columns:
                skip_cols.append(i)
        for i in skip_cols:
            to_fill.remove(i)
        log.write("  -Skipping columns: ",skip_cols, verbose=verbose) 
    if len(set(to_fill) & set(["OR","OR_95L","OR_95U","BETA","SE","P","Z","CHISQ","MLOG10P","MAF","SIG"]))==0:
        log.write(" -No available columns to fill. Skipping.", verbose=verbose)
        log.write("Finished filling data using existing columns.", verbose=verbose)
        return sumstats
    log.write(" -Filling columns: ",to_fill, verbose=verbose)
    fill_iteratively(sumstats,to_fill,log,only_sig,df,extreme,verbose,sig_level)
       
# ###################################################################################
    #sumstats = sortcolumn(sumstats, verbose=verbose, log=log)
    gc.collect()
    log.write("Finished filling data using existing columns.", verbose=verbose)
    return sumstats
    
##########################################################################################################################    
    
def fill_p(sumstats,log,df=None,only_sig=False,sig_level=5e-8,overwrite=False,verbose=True,filled_count=0):
        # MLOG10P -> P
    if "MLOG10P" in sumstats.columns:    
        log.write("  - Filling P value using MLOG10P column...", verbose=verbose)
        sumstats["P"] = np.power(10,-sumstats["MLOG10P"])
        filled_count +=1

    # Z -> P
    elif "Z" in sumstats.columns:
        log.write("  - Filling P value using Z column...", verbose=verbose)
        stats.chisqprob = lambda chisq, degree_of_freedom: stats.chi2.sf(chisq, degree_of_freedom)
        sumstats["P"] = ss.chisqprob(sumstats["Z"]**2,1)
        filled_count +=1

    elif "CHISQ" in sumstats.columns:
    #CHISQ -> P
        log.write("  - Filling P value using CHISQ column...", verbose=verbose)
        stats.chisqprob = lambda chisq, degree_of_freedom: stats.chi2.sf(chisq, degree_of_freedom)
        if df is None:
            if only_sig is True and overwrite is True:
                sumstats.loc[sumstats["P"]<sig_level,"P"] = stats.chisqprob(sumstats.loc[sumstats["P"]<sig_level,"CHISQ"],1)
                filled_count +=1
            else:
                sumstats["P"] = stats.chisqprob(sumstats["CHISQ"],1)
                filled_count +=1
        else:
            if only_sig is True and overwrite is True:
                log.write("  - Filling P value using CHISQ column for variants:" , sum(sumstats["P"]<sig_level), verbose=verbose)
                sumstats.loc[sumstats["P"]<sig_level,"P"] = stats.chisqprob(sumstats.loc[sumstats["P"]<sig_level,"CHISQ"],sumstats.loc[sumstats["P"]<sig_level,df].astype("int"))
                filled_count +=1
            else:
                log.write("  - Filling P value using CHISQ column for all valid variants:", verbose=verbose)
                sumstats["P"] = stats.chisqprob(sumstats["CHISQ"],sumstats[df].astype("int"))
                filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_z(sumstats,log,verbose=True,filled_count=0):
    # BETA/SE -> Z
    if ("BETA" in sumstats.columns) and ("SE" in sumstats.columns):
        log.write("  - Filling Z using BETA/SE column...", verbose=verbose)
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_chisq(sumstats,log,verbose=True,filled_count=0):
    # Z -> CHISQ
    if "Z" in sumstats.columns:
        log.write("  - Filling CHISQ using Z column...", verbose=verbose)
        sumstats["CHISQ"] = (sumstats["Z"])**2
        filled_count +=1
    elif "P" in sumstats.columns:
    # P -> CHISQ
        log.write("  - Filling CHISQ using P column...", verbose=verbose)
        sumstats["CHISQ"] = ss.chi2.isf(sumstats["P"], 1)
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_or(sumstats,log,verbose=True,filled_count=0):
    # BETA -> OR
    if "BETA" in sumstats.columns:
        log.write("  - Filling OR using BETA column...", verbose=verbose)
        sumstats["OR"]   = np.exp(sumstats["BETA"])
        filled_count +=1
        # BETA/SE -> OR_95L / OR_95U
        # get confidence interval 95
        if ("BETA" in sumstats.columns) and ("SE" in sumstats.columns):
            log.write("  - Filling OR_95L/OR_95U using BETA/SE columns...", verbose=verbose)
            # beta - 1.96 x se , beta + 1.96 x se
            sumstats["OR_95L"] = np.exp(sumstats["BETA"]-ss.norm.ppf(0.975)*sumstats["SE"])
            sumstats["OR_95U"] = np.exp(sumstats["BETA"]+ss.norm.ppf(0.975)*sumstats["SE"])
            filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count
def fill_or95(sumstats,log,verbose=True,filled_count=0):
    # get confidence interval 95
    if ("BETA" in sumstats.columns) and ("SE" in sumstats.columns):
        log.write("  - Filling OR_95L/OR_95U using BETA/SE columns...", verbose=verbose)
        # beta - 1.96 x se , beta + 1.96 x se
        sumstats["OR_95L"] = np.exp(sumstats["BETA"]-ss.norm.ppf(0.975)*sumstats["SE"])
        sumstats["OR_95U"] = np.exp(sumstats["BETA"]+ss.norm.ppf(0.975)*sumstats["SE"])
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count
    
def fill_beta(sumstats,log,verbose=True,filled_count=0):
    # OR -> beta
    if "OR" in sumstats.columns:
        log.write("  - Filling BETA value using OR column...", verbose=verbose)
        sumstats["BETA"]  = np.log(sumstats["OR"])    
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_se(sumstats,log,verbose=True,filled_count=0):
    # OR / OR_95L /OR_95U -> SE
    if ("P" in sumstats.columns) and ("BETA" in sumstats.columns):
        log.write("  - Filling SE value using BETA and P column...", verbose=verbose)
        sumstats["SE"]= np.abs(sumstats["BETA"]/ ss.norm.ppf(1-sumstats["P"]/2))
        filled_count +=1
    elif ("OR" in sumstats.columns) and ("OR_95U" in sumstats.columns): 
        log.write("  - Filling SE value using OR/OR_95U column...", verbose=verbose)
        # 
        sumstats["SE"]=(np.log(sumstats["OR_95U"]) - np.log(sumstats["OR"]))/ss.norm.ppf(0.975)
        filled_count +=1
    elif ("OR" in sumstats.columns) and ("OR_95L" in sumstats.columns):
        log.write("  - Filling SE value using OR/OR_95L column...", verbose=verbose)
        sumstats["SE"]=(np.log(sumstats["OR"]) - np.log(sumstats["OR_95L"]))/ss.norm.ppf(0.975)
        filled_count +=1
    else:
        log.write("  - Not enough information to fill SE...", verbose=verbose)
        return 0,filled_count
    return 1,filled_count

def fill_mlog10p(sumstats,log,verbose=True,filled_count=0):
    if "P" in sumstats.columns:
        # P -> MLOG10P
        log.write("  - Filling MLOG10P using P column...", verbose=verbose)
        sumstats["MLOG10P"] = -np.log10(sumstats["P"])
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_extreme_mlog10p(sumstats,df,log,verbose=True,filled_count=0):
    # ref: https://stackoverflow.com/questions/46416027/how-to-compute-p-values-from-z-scores-in-r-when-the-z-score-is-large-pvalue-muc/46416222#46416222
    if "Z" in sumstats.columns:
        # P -> MLOG10P
        log.write("  - Filling MLOG10P using Z column...", verbose=verbose)
        sumstats = fill_extreme_mlog10(sumstats, "Z")
        filled_count +=1
    elif "BETA" in sumstats.columns and "SE" in sumstats.columns:
        log.write("  - Z column not available...", verbose=verbose)
        log.write("  - Filling Z using BETA/SE column...", verbose=verbose)
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
        log.write("  - Filling MLOG10P using Z column...", verbose=verbose)
        sumstats = fill_extreme_mlog10(sumstats, "Z")
        filled_count +=1
    elif "CHISQ" in sumstats.columns and "DOF" in sumstats.columns:
        log.write("  - Filling MLOG10P using CHISQ and DOF column...", verbose=verbose)
        sumstats = fill_extreme_mlog10_chisq(sumstats, "CHISQ", df)
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_maf(sumstats,log,verbose=True,filled_count=0):
    if "EAF" in sumstats.columns:
        # EAF -> MAF
        log.write("  - Filling MAF using EAF column...", verbose=verbose)
        sumstats["MAF"] =  sumstats["EAF"].apply(lambda x: min(x,1-x) if pd.notnull(x) else np.nan)
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

def fill_sig(sumstats,log,sig_level=5e-8, verbose=True,filled_count=0):
    if "P" in sumstats.columns or "MLOG10P" in sumstats.columns:
        log.write("  - Determining significant using P and MLOG10P with threshold:{}".format(sig_level), verbose=verbose)
        if "P" in sumstats.columns:
            is_sig = sumstats["P"]<sig_level 
        elif "MLOG10P" in sumstats.columns:
            is_sig = sumstats["MLOG10P"]>np.log10(sig_level) 
        sumstats["SIGNIFICANT"] = False
        sumstats.loc[is_sig, "SIGNIFICANT"] = True
        filled_count +=1
    else:
        return 0,filled_count
    return 1,filled_count

####################################################################################################################
def fill_extreme_mlog10(sumstats, z):
    log_pvalue = np.log(2) + ss.norm.logsf(np.abs(sumstats[z])) #two-sided
    log10_pvalue = log_pvalue/np.log(10)
    mantissa = 10**(log10_pvalue %1 )
    exponent = log10_pvalue // 1
    sumstats["MLOG10P"] = -log10_pvalue
    sumstats["P_MANTISSA"]= mantissa
    sumstats["P_EXPONENT"]= exponent
    return sumstats

def fill_extreme_mlog10_chisq(sumstats, chisq, df):
    #https://stackoverflow.com/a/46416222/199475
    log_pvalue = ss.chi2.logsf(sumstats[chisq], sumstats[df])

    log10_pvalue = log_pvalue/np.log(10)
    
    mantissa = 10**(log10_pvalue %1)
    exponent = log10_pvalue // 1
    sumstats["MLOG10P"] = -log10_pvalue
    sumstats["P_MANTISSA"]= mantissa
    sumstats["P_EXPONENT"]= exponent
    return sumstats

####################################################################################################################
def fill_iteratively(sumstats,raw_to_fill,log,only_sig,df,extreme,verbose,sig_level):
    to_fill = raw_to_fill.copy()
    log.write("  - Filling Columns iteratively...", verbose=verbose)

    filled_count=0
    for i in range(len(to_fill)+1):
    # beta to or ####################################################################################################     
        if "OR" in to_fill:
            status, filled_count = fill_or(sumstats,log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("OR")
    # or to beta #################################################################################################### 
        if "BETA" in to_fill:
            status,filled_count = fill_beta(sumstats,log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("BETA")
        if "SE" in to_fill:
            status,filled_count = fill_se(sumstats,log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("SE")
    # z/chi2 to p ##################################################################################################
        if "P" in to_fill:
            status,filled_count = fill_p(sumstats,log,only_sig=only_sig,df=df,sig_level=sig_level,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("P")
    # beta/se to z ##################################################################################################            
        if "Z" in to_fill:   
            status,filled_count = fill_z(sumstats,log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("Z")
    # z/p to chisq ##################################################################################################             
        if "CHISQ" in to_fill:
            status,filled_count = fill_chisq(sumstats,log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("CHISQ")
    # EAF to MAF ##################################################################################################   
        if "MAF" in to_fill:
            status,filled_count = fill_maf(sumstats,log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("MAF")
    # p to -log10(P)  ###############################################################################################
        if "MLOG10P" in to_fill:
            if extreme==True:
                status,filled_count = fill_extreme_mlog10p(sumstats,df, log,verbose=verbose,filled_count=filled_count)
                filled_count +=1
            elif "P" not in sumstats.columns:
                fill_p(sumstats,log,verbose=verbose)
                status,filled_count = fill_mlog10p(sumstats,log,verbose=verbose,filled_count=filled_count)
                sumstats = sumstats.drop(labels=["P"],axis=1)
            else:
                status,filled_count = fill_mlog10p(sumstats,log,verbose=verbose)
            if status == 1 : to_fill.remove("MLOG10P")

        if "SIG" in to_fill:
            status,filled_count = fill_sig(sumstats,sig_level=sig_level ,log=log,verbose=verbose,filled_count=filled_count)
            if status == 1 : to_fill.remove("SIG")
        if filled_count == 0:
            break
         
###Base functions########################################################################################

def _convert_betase_to_z(beta, se):
    return beta/se 

def _convert_betase_to_p(beta, se):
    z = _convert_betase_to_z(beta, se)
    p = _convert_z_to_p(z)
    return p

def _convert_betase_to_mlog10p(beta, se):
    z = _convert_betase_to_z(beta, se)
    mlog10p = _convert_z_to_mlog10p(z)
    return mlog10p

def _convert_p_to_chisq(p):
    return ss.chi2.isf(p, 1)

def _convert_z_to_chisq(z):
    return (z)**2

def _convert_z_to_p(z):
    return ss.chi2.sf(z**2,1) 

def _convert_z_to_mlog10p(z):
    log_pvalue = np.log(2) + ss.norm.logsf(np.abs(z)) #two-sided
    mlog10p = log_pvalue/np.log(10)
    return -mlog10p 

def _conver_chisq_to_p(chisq):
    return ss.chi2.sf(chisq,1)

def _convert_mlog10p_to_p(mlog10p):
    return np.power(10, -mlog10p) 

def _convert_or_to_beta(OR):
    return np.log(OR)  

def _convert_beta_to_or(beta):
    return np.exp(beta)   

def rank_based_int(series, c=3/8):
    #https://onlinelibrary.wiley.com/doi/10.1111/biom.13214
    n=sum(~series.isna())
    normalized_value = norm.ppf((series.rank()-c)/(n+1-2*c))
    return normalized_value