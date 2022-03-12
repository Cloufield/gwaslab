import pandas as pd
import numpy as np
import scipy.stats as ss

#20220310
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
              other=[],
              verbose=False,
              readargs={"sep": "\s+"}):

    #renaming dictionary
    rename_dictionary = {
        snpid: "MARKERNAME",
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
        OR_95U: "OR_95U"
    }

    usecols = []
    datatype_dictionary = {}
    
    if snpid:
        usecols.append(snpid)
        datatype_dictionary[snpid] = "string"
    if rsid:
        usecols.append(rsid)
        datatype_dictionary[rsid] = "string"
    if chrom:
        usecols.append(chrom)
        datatype_dictionary[chrom] = "string"
    if pos:
        usecols.append(pos)
        datatype_dictionary[pos] = "int"
    if ea:
        usecols.append(ea)
        datatype_dictionary[ea] = "string"
    if nea:
        usecols.append(nea)
        datatype_dictionary[nea] = "string"
    if eaf:
        usecols.append(eaf)
        datatype_dictionary[eaf] = "float"
    if n and (type(n) is str):
        usecols.append(n)
        datatype_dictionary[n] = "int"
    if beta:
        usecols.append(beta)
        datatype_dictionary[beta] = "float"
    if se:
        usecols.append(se)
        datatype_dictionary[se] = "float"
    if chisq:
        usecols.append(chisq)
        datatype_dictionary[chisq] = "float"
    if z:
        usecols.append(z)
        datatype_dictionary[z] = "float"
    if p:
        usecols.append(p)
        datatype_dictionary[p] = "float"
    if mlog10p:
        usecols.append(mlog10p)
        datatype_dictionary[mlog10p] = "float"
    if info:
        usecols.append(info)
        datatype_dictionary[info] = "float"
    if OR:
        usecols.append(OR)
        datatype_dictionary[OR] = "float"
    if OR_se:
        usecols.append(OR_se)
        datatype_dictionary[OR_se] = "float"
    if OR_95L:
        usecols.append(OR_95L)
        datatype_dictionary[OR_95L] = "float"
    if OR_95U:
        usecols.append(OR_95U)
        datatype_dictionary[OR_95U] = "float"
    if other:
        usecols = usecols + other
        for i in other:
            rename_dictionary[i] = i
            datatype_dictionary[i] = "string"


#loading data ##################################################################################
    if type(sumstats) is str:
        ## loading data from path
        inpath = sumstats
        if verbose:
            print("Initiating from file :" + inpath)
        sumstats = pd.read_table(inpath,
                                 usecols=usecols,
                                 dtype=datatype_dictionary,
                                 **readargs)

    elif type(sumstats) is pd.DataFrame:
        ## loading data from dataframe
        if verbose:
            print("Initiating from pandas DataFrame ...")
        sumstats = sumstats.loc[:, usecols].astype(datatype_dictionary)

    ## renaming columns
    converted_columns = list(map(lambda x: rename_dictionary[x], usecols))
    datatype_columns = list(map(lambda x: datatype_dictionary[x], usecols))
    
    ## renaming log
    if verbose: print("  - Reading columns          :", ",".join(usecols))
    if verbose: print("  - Renaming columns to      :", ",".join(converted_columns))
    if verbose: print("  - Datatype of each columns :", ",".join(datatype_columns))
    if verbose:
        print("  - Current dataframe shape  : Rows ", len(sumstats), " x ",
              len(sumstats.columns), " Columns")
    ## renaming
    sumstats = sumstats.rename(columns=rename_dictionary)

    ## if n was provided as int
    if type(n) is int:
        sumstats["N"] = n
    ## if markername was not provided
    if (not snpid) and (not rsid):
        sumstats["MARKERNAME"] = sumstats["CHR"].astype(
            "string") + ":" + sumstats["POS"].astype("string")
    
    ##reodering 
    order = [
        "MARKERNAME","rsID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "Z",
        "CHISQ", "P", "MLOG10P", "OR", "OR_SE", "OR_95L", "OR_95U", "INFO", "N"
    ] + other
    output_columns = []
    for i in order:
        if i in sumstats.columns: output_columns.append(i)
            
    if verbose: print("  - Reordering columns to    :", ",".join(output_columns))
    sumstats = sumstats.loc[:, output_columns]
    if verbose: print("Loading data finished successfully!")
    return sumstats