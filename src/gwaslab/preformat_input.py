import pandas as pd
import numpy as np
import scipy.stats as ss

#20220425
def preformat(sumstats,
          snpid=None,
          rsid=None,
          chrom=None,
          pos=None,
          ea=None,
          nea=None,
          eaf=None,
          n=None,
          beta=None,
          se=None,
          chisq=None,
          z=None,
          p=None,
          mlog10p=None,
          info=None,
          OR=None,
          OR_se=None,
          OR_95L=None,
          OR_95U=None,
          direction=None,
          status=None,
          build=None,
          other=[],
          verbose=False,
          readargs={"sep": "\s+"},
          log=None):

    #renaming dictionary
    rename_dictionary = {
        snpid: "SNPID",
        rsid: "rsID",
        chrom: "CHR",
        pos: "POS",
        nea: "NEA",
        ea: "EA",
        eaf: "EAF",
        n: "N",
        beta: "BETA",
        se: "SE",
        z: "Z",
        chisq: "CHISQ",
        p: "P",
        mlog10p: "MLOG10P",
        info: "INFO",
        OR: "OR",
        OR_se: "OR_SE",
        OR_95L: "OR_95L",
        OR_95U: "OR_95U",
        direction:"DIRECTION",
        status:"STATUS"
    }
    
    usecols = []
    
    if snpid:
        usecols.append(snpid)
    if rsid:
        usecols.append(rsid)
    if chrom:
        usecols.append(chrom)
    if pos:
        usecols.append(pos)
    if ea:
        usecols.append(ea)
    if nea:
        usecols.append(nea)
    if eaf:
        usecols.append(eaf)
    if n and (type(n) is str):
        usecols.append(n)
    if beta:
        usecols.append(beta)
    if se:
        usecols.append(se)
    if chisq:
        usecols.append(chisq)
    if z:
        usecols.append(z)
    if p:
        usecols.append(p)
    if mlog10p:
        usecols.append(mlog10p)
    if info:
        usecols.append(info)
    if OR:
        usecols.append(OR)
    if OR_se:
        usecols.append(OR_se)
    if OR_95L:
        usecols.append(OR_95L)
    if OR_95U:
        usecols.append(OR_95U)
    if direction:
        usecols.append(direction)
    if status:
        usecols.append(status)
    if other:
        usecols = usecols + other
        for i in other:
            rename_dictionary[i] = i
            
 #loading data ##################################################################################
    
    try:
        if type(sumstats) is str:
            ## loading data from path
            inpath = sumstats
            if verbose: log.write("Initiating from file :" + inpath)
            sumstats = pd.read_table(inpath,
                             usecols=usecols,
                             dtype="string",
                             **readargs)
        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            if verbose: log.write("Initiating from pandas DataFrame ...")
            sumstats = sumstats.loc[:, usecols].astype("string")
    except ValueError:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure it contains the columns.")

    ## renaming columns
    converted_columns = list(map(lambda x: rename_dictionary[x], usecols))
    
    ## renaming log
    if verbose: log.write(" -Reading columns          :", ",".join(usecols))
    if verbose: log.write(" -Renaming columns to      :", ",".join(converted_columns))
    if verbose: log.write(" -Current dataframe shape  : Rows ", len(sumstats), " x ",len(sumstats.columns), " Columns")
    
    ## renaming
    sumstats = sumstats.rename(columns=rename_dictionary)

    ## if n was provided as int
    if type(n) is int:
        sumstats["N"] = n
        
    ## if markername was not provided
    #if (not snpid) and (not rsid):
    #    if ("CHR" in sumstats.columns) and ("POS" in sumstats.columns):
    #        sumstats["SNPID"] = sumstats["CHR"].astype(
    #            "string") + ":" + sumstats["POS"].astype("string")
    #    else:
    #        raise ValueError("Please input at least 1.rsID or 2.snpID or 3.CHR and POS")
    
    if status is None:
        if verbose: log.write(" -Initiating a status column ...")
        sumstats["STATUS"] = build +"99999"
        sumstats["STATUS"] = sumstats["STATUS"].astype("string")
    
    ##reodering 
    order = [
        "SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "Z",
        "CHISQ", "P", "MLOG10P", "OR", "OR_SE", "OR_95L", "OR_95U", "INFO", "N","DIRECTION","STATUS"
           ] + other
    
    output_columns = []
    for i in order:
        if i in sumstats.columns: output_columns.append(i)

    if verbose: log.write(" -Reordering columns to    :", ",".join(output_columns))
    sumstats = sumstats.loc[:, output_columns]
    if verbose: log.write("Loading data finished successfully!")
    

    return sumstats