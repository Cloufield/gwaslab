import pandas as pd
import numpy as np
import scipy.stats as ss
import gzip
import os
import gc
from gwaslab.CommonData import get_format_dict
from gwaslab.fixdata import sortcolumn

#20221030
def preformat(sumstats,
          fmt=None,
          snpid=None,
          rsid=None,
          chrom=None,
          pos=None,
          ea=None,
          nea=None,
          ref=None,
          alt=None,
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
          OR_95L=None,
          OR_95U=None,
          direction=None,
          status=None,
          study=None,
          build=None,
          other=[],
          verbose=False,
          readargs={"sep": "\t"},
          log=None):

    #renaming dictionary
    rename_dictionary = {}
    usecols = []
    dtype_dictionary ={}    
    
 #######################################################################################################################################################
    if fmt is not None:
        if verbose: log.write("Start to load format from formatbook....")
        
        # load format data
        meta_data,rename_dictionary = get_format_dict(fmt)
        
        ########## print format information################################################
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log)
                
#########################################################################################################################################################      
    
    # check chr-separated path / vcf / then print header.            
    try:
        if type(sumstats) is str:
            ## loading data from path #################################################
            inpath = sumstats
            ###load sumstats by each chromosome #################################################  
            if "@" in inpath:
                if verbose: log.write(" -Detected @ in path: load sumstats by each chromosome...")
                inpath_chr_list=[]
                inpath_chr_num_list=[]
                for chromosome in list(range(1,26))+["x","y","X","Y","MT","mt","m","M"]:
                    inpath_chr = inpath.replace("@",str(chromosome))  
                    if isfile_casesensitive(inpath_chr):
                        inpath_chr_num_list.append(str(chromosome))
                        inpath_chr_list.append(inpath_chr)        
                if verbose: log.write(" -Chromosomes detected:",",".join(inpath_chr_num_list))
                readargs_header = get_readargs_header(inpath = inpath_chr_list[0], readargs = readargs)
                row_one = pd.read_table(inpath_chr_list[0],**readargs_header)
                raw_cols = row_one.columns
            else:
            ##### loading data from tabular file#################################################
                readargs_header = get_readargs_header(inpath = inpath, readargs = readargs)
                row_one = pd.read_table(inpath,**readargs_header)
                raw_cols = row_one.columns
            
            if fmt=="vcf":
                # expanded
                format_cols = list(row_one["FORMAT"].str.split(":"))[0]
                # fixed + study1 + expanded
                raw_cols = meta_data["format_fixed"] + [raw_cols[9]] + format_cols

            ###################################################################################### 
        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            raw_cols = sumstats.columns
        
        ################################################
        for key,value in rename_dictionary.items():
            if key in raw_cols:
                usecols.append(key)
            if value in ["EA","NEA"]:
                dtype_dictionary[value]="category"
            if value in ["CHR","STATUS"]:
                dtype_dictionary[value]="string"     
    
    except ValueError:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure the columns you specified are in the file.")
        
    ###################################################################################################################################################
    ## check columns/datatype to use 
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
        dtype_dictionary[ea]="string"
    if nea:
        usecols.append(nea)
        rename_dictionary[nea]= "NEA"
        dtype_dictionary[nea]="string"
    if ref:
        usecols.append(ref)
        rename_dictionary[nea]= "REF"
        dtype_dictionary[nea]="string"
    if alt:
        usecols.append(alt)
        rename_dictionary[nea]= "ALT"
        dtype_dictionary[nea]="string"
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
    if fmt=="vcf":
        # store the final column list
        vcf_usecols = usecols.copy()
        # loading the fixed columns + study
        usecols = meta_data["format_fixed"]
        if study is not None:
            usecols =  usecols + [study]
        else:
            study = raw_cols[9]
            usecols =  usecols + [study]
 #loading data ##########################################################################################################
    
    try:
        if type(sumstats) is str:
            ## loading data from path
            inpath = sumstats
            if "@" in inpath:
                if verbose: log.write("Start to initiate from files with pattern :" + inpath)
                sumstats_chr_list=[]
                for i in inpath_chr_list:
                    if verbose: log.write(" -Loading:" + i)
                    skip_rows = get_skip_rows(i)
                    readargs["skiprows"] = skip_rows
                    sumstats_chr = pd.read_table(i,
                                        usecols=set(usecols),
                                        dtype=dtype_dictionary,
                                        **readargs)
                    sumstats_chr_list.append(sumstats_chr)
                if verbose: log.write(" -Merging sumstats for chromosomes:",",".join(inpath_chr_num_list))
                sumstats = pd.concat(sumstats_chr_list, axis=0, ignore_index=True) 
                del(sumstats_chr_list)
                gc.collect()
            else:
                skip_rows = get_skip_rows(inpath)
                readargs["skiprows"] = skip_rows
                if verbose: log.write("Start to initiate from file :" + inpath)
                sumstats = pd.read_table(inpath,
                                 usecols=set(usecols),
                                 dtype=dtype_dictionary,
                                 **readargs)
        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            if verbose: log.write("Start to initiate from pandas DataFrame ...")
            sumstats = sumstats.loc[:, usecols]
    except ValueError:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure it contain the columns.")

    ## renaming columns ###############################################################################################
    if fmt == "vcf":
        sumstats = parse_vcf_study(sumstats,format_cols,study,vcf_usecols,log=log,verbose=verbose)
        usecols = vcf_usecols
    
    converted_columns = list(map(lambda x: rename_dictionary[x], set(usecols)))
    
    ## renaming log
    if verbose: log.write(" -Reading columns          :", ",".join(set(usecols)))
    if verbose: log.write(" -Renaming columns to      :", ",".join(converted_columns))
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns)) 
    
    ## renaming  #####################################################################################
    sumstats = sumstats.rename(columns=rename_dictionary)

    ## if n was provided as int #####################################################################################
    if type(n) is int:
        sumstats["N"] = n
    
    ### status ######################################################################################################
    if status is None:
        sumstats = process_status(sumstats=sumstats,build=build,log=log,verbose=verbose)
    
    ## ea/nea, ref/alt ##############################################################################################
    sumstats = process_allele(sumstats=sumstats,log=log,verbose=verbose)
        
    ## NEAF to EAF ###########################################################################################################
    if neaf is not None :
        sumstats = process_neaf(sumstats=sumstats,log=log,verbose=verbose)

    ## reodering ###################################################################################################  
    sumstats = sortcolumn(sumstats=sumstats,log=log,verbose=verbose)    
    
    gc.collect()
    if verbose: log.write("Finished loading data successfully!")
    return sumstats


#### helper #######################################################################
def isfile_casesensitive(path):
    if not os.path.isfile(path): 
        return False   # exit early
    directory, filename = os.path.split(path)
    return filename in os.listdir(directory)

def get_readargs_header(inpath,readargs):
    if "vcf.gz" in inpath:
        with gzip.open(inpath,'r') as file:      
            skip=0
            for line in file:        
                if line.decode('utf-8').startswith('##'):
                    skip+=1
                else:
                    readargs["skiprows"]=skip
                    readargs["sep"]="\t"
                    break
    readargs_header = readargs.copy()
    readargs_header["nrows"]=1
    readargs_header["dtype"]="string"
    return readargs_header

def get_skip_rows(inpath):
    if "vcf.gz" in inpath:
        with gzip.open(inpath,'r') as file:      
            skip=0
            for line in file:        
                if line.decode('utf-8').startswith('##'):
                    skip+=1
                else:
                    return skip
    else:
        return 0

def parse_vcf_study(sumstats,format_cols,study,vcf_usecols,log,verbose=True):
    if verbose: log.write(" -Parsing based on FORMAT: ", format_cols)
    if verbose: log.write(" -Parsing vcf study : ", study)
    sumstats[format_cols] = sumstats[study].str.split(":",expand=True).values
    sumstats = sumstats.drop(["FORMAT",study],axis=1)
    sumstats = sumstats.loc[:, vcf_usecols]
    gc.collect()
    return sumstats

def print_format_info(fmt,meta_data, rename_dictionary, verbose, log):
    if verbose: log.write(" -"+fmt+" format meta info:")   
    for key,value in meta_data.items():
        if type(value) is str:
            if "\n" in value:
                value_first_line=value.split("\n")[0]
                log.write("  -",key," : "+value_first_line.strip()+"...")  
            else:
                log.write("  -",key," : "+value.strip())  
        else:
            log.write("  -",key," : ",value)  
    keys=[]
    values=[]
    for key,value in rename_dictionary.items():
        keys.append(key)
        values.append(value)
    if fmt!="gwaslab":
        log.write(" -"+fmt+" format dictionary:")  
        log.write("  - "+fmt+" keys:",",".join(keys)) 
        log.write("  - gwaslab values:",",".join(values))  

def process_neaf(sumstats,log,verbose):
    if verbose: log.write(" -NEAF is specified...") 
    pre_number=len(sumstats)
    if verbose: log.write(" -Checking if 0<= NEAF <=1 ...") 
    sumstats.loc[:,"EAF"] = pd.to_numeric(sumstats.loc[:,"EAF"], errors='coerce')
    sumstats = sumstats.loc[(sumstats["EAF"]>=0) & (sumstats["EAF"]<=1),:]
    sumstats.loc[:,"EAF"] = 1- sumstats.loc[:,"EAF"]
    if verbose: log.write(" -Converted NEAF to EAF.") 
    after_number=len(sumstats)
    if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad NEAF.") 
    return sumstats

def process_allele(sumstats,log,verbose):
    if "EA" in sumstats.columns:
        if "REF" in sumstats.columns and "ALT" in sumstats.columns:
            if "NEA" not in sumstats.columns:
                if verbose: log.write(" NEA not available: assigning REF to NEA...") 
                sumstats["NEA"]=sumstats["REF"]    
            if verbose: log.write(" -EA,REF and ALT columns are available: assigning NEA...") 
            ea_alt = sumstats["EA"]==sumstats["ALT"]
            if verbose: log.write(" -For variants with EA == ALT : assigning REF to NEA ...") 
            sumstats.loc[ea_alt,"NEA"] = sumstats.loc[ea_alt,"REF"]
            
            ea_not_alt = sumstats["EA"]!=sumstats["ALT"]
            if verbose: log.write(" -For variants with EA != ALT : assigning ALT to NEA ...") 
            sumstats.loc[ea_not_alt,"NEA"] = sumstats.loc[ea_not_alt,"ALT"]

            #sumstats = sumstats.drop(labels=["REF","ALT"],axis=1)
            sumstats["REF"]=sumstats["REF"].astype("category") 
            sumstats["ALT"]=sumstats["ALT"].astype("category") 
        sumstats["EA"]=sumstats["EA"].astype("category")     
    if "NEA" in sumstats.columns:
        sumstats["NEA"]=sumstats["NEA"].astype("category")  
    return sumstats

def process_status(sumstats,build,log,verbose):
    if verbose: log.write(" -Initiating a status column: STATUS ...")
    #sumstats["STATUS"] = int(build)*(10**5) +99999
    sumstats["STATUS"] = build +"99999"
    categories = {str(j+i) for j in [1900000,3800000,9700000,9800000,9900000] for i in range(0,100000)}
    sumstats["STATUS"] = pd.Categorical(sumstats["STATUS"],categories=categories)
    return sumstats
