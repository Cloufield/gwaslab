import pandas as pd
import numpy as np
import scipy.stats as ss
from scipy import stats
from gwaslab.Log import Log

def filldata( 
    sumstats,
    to_fill=[],
    df=None,
    overwrite=False,
    verbose=True,
    only_sig=False,
    log = Log()
    ):
    
    if type(to_fill) is str:
        to_fill = [to_fill]

    if verbose: log.write("Start filling data using existing columns...")
    if verbose: log.write(" -Raw input columns: ",list(sumstats.columns))
# check dupication ##############################################################################################
    skip_cols=[]
    if verbose: log.write(" -Overwrite mode: ",overwrite)
    if overwrite is False:
        for i in to_fill:
            if i in sumstats.columns:
                skip_cols.append(i)
        for i in skip_cols:
            to_fill.remove(i)
        if verbose: log.write("  - Skipping columns: ",skip_cols) 
    if verbose: log.write("Filling columns: ",to_fill)
        
# beta to or ####################################################################################################     
    if "OR" in to_fill:
        # BETA -> OR
        if "BETA" in sumstats.columns:
            if verbose: log.write("  - Filling OR using BETA column...")
            sumstats["OR"]   = np.exp(sumstats["BETA"])
        
        # BETA/SE -> OR_95L / OR_95U
        # get confidence interval 95
        if ("BETA" in sumstats.columns) and ("SE" in sumstats.columns):
            if verbose: log.write("  - Filling OR_95L/OR_95U using BETA/SE columns...")
            sumstats["OR_95L"] = np.exp(sumstats["BETA"]-ss.norm.ppf(0.975)*sumstats["SE"])
            sumstats["OR_95U"] = np.exp(sumstats["BETA"]+ss.norm.ppf(0.975)*sumstats["SE"])
# or to beta #################################################################################################### 
    if "BETA" in to_fill:
        # OR -> beta
        if "OR" in sumstats.columns:
            if verbose: log.write("  - Filling BETA value using OR column...")
            sumstats["BETA"]  = np.log(sumstats["OR"])
    
    if "SE" in to_fill:
        # OR / OR_95L /OR_95U -> SE
        if ("OR" in sumstats.columns) and ("OR_95U" in sumstats.columns): 
            if verbose: log.write("  - Filling SE value using OR/OR_95U column...")
            sumstats["SE"]=(np.log(sumstats["OR_95U"]) - np.log(sumstats["OR"]))/ss.norm.ppf(0.975)
        elif ("OR" in sumstats.columns) and ("OR_95L" in sumstats.columns):
            if verbose: log.write("  - Filling SE value using OR/OR_95L column...")
            sumstats["SE"]=(np.log(sumstats["OR"]) - np.log(sumstats["OR_95L"]))/ss.norm.ppf(0.975)

# z/chi2 to p ##################################################################################################
    if "P" in to_fill:
        # Z -> P
        if "Z" in sumstats.columns:
            if verbose: log.write("  - Filling P value using Z column...")
            sumstats["P"] = ss.chisqprob(sumstats["Z"]**2,1)
        elif "CHISQ" in sumstats.columns:
        #CHISQ -> P
            if verbose: log.write("  - Filling P value using CHISQ column...")
            stats.chisqprob = lambda chisq, degree_of_freedom: stats.chi2.sf(chisq, degree_of_freedom)
            if df is None:
                if only_sig is True and overwrite is True:
                    sumstats.loc[sumstats["P"]<5e-8,"P"] = stats.chisqprob(sumstats.loc[sumstats["P"]<5e-8,"CHISQ"],1)
                else:
                    sumstats["P"] = stats.chisqprob(sumstats["CHISQ"],1)
            else:
                if only_sig is True and overwrite is True:
                    if verbose: log.write("  - Filling P value using CHISQ column for variants:" , sum(sumstats["P"]<5e-8))
                    sumstats.loc[sumstats["P"]<5e-8,"P"] = stats.chisqprob(sumstats.loc[sumstats["P"]<5e-8,"CHISQ"],sumstats.loc[sumstats["P"]<5e-8,df].astype("int"))
                else:
                    if verbose: log.write("  - Filling P value using CHISQ column for all valid variants:")
                    sumstats["P"] = stats.chisqprob(sumstats["CHISQ"],sumstats[df].astype("int"))
        # MLOG10P -> P
        elif "MLOG10P" in sumstats.columns:    
            if verbose: log.write("  - Filling P value using MLOG10P column...")
            sumstats["P"] = np.power(10,-sumstats["MLOG10P"])

# beta/se to z ##################################################################################################            
    if "Z" in to_fill:   
        # BETA/SE -> Z
        if ("BETA" in sumstats.columns) and ("SE" in sumstats.columns):
            if verbose: log.write("  - Filling Z using BETA/SE column...")
            sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

# z/p to chisq ##################################################################################################             
    if "CHISQ" in to_fill:
        # Z -> CHISQ
        if "Z" in sumstats.columns:
            if verbose: log.write("  - Filling CHISQ using Z column...")
            sumstats["CHISQ"] = (sumstats["Z"])**2
        elif "P" in sumstats.columns:
        # P -> CHISQ
            if verbose: log.write("  - Filling CHISQ using P column...")
            sumstats["CHISQ"] = ss.chi2.isf(sumstats["P"], 1)
            
# p to -log10(P)  ###############################################################################################
    if "MLOG10P" in to_fill:
        # P -> MLOG10P
        if verbose: log.write("  - Filling MLOG10P using P column...")
        sumstats["MLOG10P"] = -np.log10(sumstats["P"])
        
# ###################################################################################
    
    core_col = ["SNPID","rsID","CHR","POS","EA","NEA","EAF","N","BETA","SE","Z","CHISQ","P","MLOG10P","INFO", "OR","OR_SE","OR_95L","OR_95U"]
    other=[]
    for i in sumstats.columns:
        if i not in core_col:
                other.append(i)
    order = core_col + other
    
    output_columns=[]
    for i in order:
        if i in sumstats.columns:
            output_columns.append(i)
     
    sumstats = sumstats.loc[:,output_columns]
    if verbose: log.write("Finished filling data using existing columns.")
    return sumstats
    