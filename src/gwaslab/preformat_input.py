#import modin.pandas as pd
import pandas as pd
import numpy as np
import scipy.stats as ss
from gwaslab.CommonData import get_format_dict

#20220425
def preformat(sumstats,
          fmt=None,
          snpid=None,
          rsid=None,
          chrom=None,
          pos=None,
          ea=None,
          nea=None,
          eaf=None,
          neaf=None,
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
          readargs={"sep": "\t"},
          log=None):

    #renaming dictionary
    rename_dictionary = {}
    usecols = []
    dtype_dictionary ={}
    
    if fmt is not None:
        if verbose: log.write("Start to load format from formatbook....")
        meta_data,rename_dictionary = get_format_dict(fmt)
        if verbose: 
            log.write(" -"+fmt+" format meta info:")   
            for key,value in meta_data.items():
                log.write("  -",key," : ",value)  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            if fmt!="gwaslab":
                log.write(" - "+fmt+" format dictionary:")  
                log.write("  - "+fmt+" keys:",keys) 
                log.write("  - gwaslab values:",values)  
    try:
        if type(sumstats) is str:
            ## loading data from path
            inpath = sumstats
            readargs_header = readargs.copy()
            readargs_header["nrows"]=1
            readargs_header["dtype"]="string"
            raw_cols = pd.read_table(inpath,**readargs_header).columns
        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            raw_cols = sumstats.columns
        for key,value in rename_dictionary.items():
            if key in raw_cols:
                usecols.append(key)
                if value in ["EA","NEA"]:
                    dtype_dictionary[value]="category"
                if value in ["CHR","STATUS"]:
                    dtype_dictionary[value]="string"
                    
        
    except ValueError:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure it contains the columns.")

    if snpid:
        usecols.append(snpid)
        rename_dictionary[snpid]= "SNPID"
    if rsid:
        usecols.append(rsid)
        rename_dictionary[rsid]= "rsID"
    if chrom:
        usecols.append(chrom)
        rename_dictionary[chrom]= "CHR"
        dtype_dictionary[chrom]="string"
    if pos:
        usecols.append(pos)
        rename_dictionary[pos]= "POS"
    if ea:
        usecols.append(ea)
        rename_dictionary[ea]= "EA"
        dtype_dictionary[ea]="category"
    if nea:
        usecols.append(nea)
        rename_dictionary[nea]= "NEA"
        dtype_dictionary[nea]="category"
    if eaf:
        usecols.append(eaf)
        rename_dictionary[eaf]= "EAF"
    elif neaf:
        usecols.append(neaf)
        rename_dictionary[neaf]= "EAF"
    if n and (type(n) is str):
        usecols.append(n)
        rename_dictionary[n]= "N"
    if beta:
        usecols.append(beta)
        rename_dictionary[beta]= "BETA"
    if se:
        usecols.append(se)
        rename_dictionary[se]= "SE"
    if chisq:
        usecols.append(chisq)
        rename_dictionary[chisq]="CHISQ"
    if z:
        usecols.append(z)
        rename_dictionary[z]= "Z"
    if p:
        usecols.append(p)
        rename_dictionary[p]= "P"
    if mlog10p:
        usecols.append(mlog10p)
        rename_dictionary[mlog10p]= "MLOG10P"
    if info:
        usecols.append(info)
        rename_dictionary[info]= "INFO"
    if OR:
        usecols.append(OR)
        rename_dictionary[OR]= "OR"
    if OR_se:
        usecols.append(OR_se)
        rename_dictionary[OR_se]= "OR_SE"
    if OR_95L:
        usecols.append(OR_95L)
        rename_dictionary[OR_95L]= "OR_95L"
    if OR_95U:
        usecols.append(OR_95U)
        rename_dictionary[OR_95U]= "OR_95U"
    if direction:
        usecols.append(direction)
        rename_dictionary[direction]="DIRECTION"
    if status:
        usecols.append(status)
        rename_dictionary[status]="STATUS"
        dtype_dictionary[status]="string"
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
                             usecols=set(usecols),
                             dtype=dtype_dictionary,
                             **readargs)
        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            if verbose: log.write("Initiating from pandas DataFrame ...")
            sumstats = sumstats.loc[:, usecols]
    except ValueError:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure it contains the columns.")

    ## renaming columns
    converted_columns = list(map(lambda x: rename_dictionary[x], set(usecols)))
    
    ## renaming log
    if verbose: log.write(" -Reading columns          :", ",".join(set(usecols)))
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
        categories = {str(j+i) for j in [1900000,3800000,9700000,9800000,9900000] for i in range(0,100000)}
        sumstats["STATUS"] = pd.Categorical(sumstats["STATUS"],categories=categories)
    
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
    
    if neaf is not None :
        if verbose: log.write(" -NEAF is specified...") 
        pre_number=len(sumstats)
        if verbose: log.write(" -Checking if 0<= NEAF <=1 ...") 
        sumstats.loc[:,"EAF"] = pd.to_numeric(sumstats.loc[:,"EAF"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["EAF"]>=0) & (sumstats["EAF"]<=1),:]
        sumstats.loc[:,"EAF"] = 1- sumstats.loc[:,"EAF"].round(4)
        if verbose: log.write(" -Converted NEAF to EAF.") 
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad NEAF.")    
    
    if verbose: log.write("Loading data finished successfully!")
    

    return sumstats