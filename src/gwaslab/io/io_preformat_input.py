import pandas as pd
import numpy as np
import scipy.stats as ss
import gzip
import os
import re
import gc
from gwaslab.bd.bd_common_data import get_format_dict
from gwaslab.qc.qc_fix_sumstats import sortcolumn
from gwaslab.qc.qc_fix_sumstats import _process_build
from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import quick_convert_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_memory_usage
from gwaslab.g_headers import _check_overlap_with_reserved_keys
from gwaslab.g_vchange_status import STATUS_CATEGORIES
from gwaslab.g_Log import Log


#20221030
def preformat(sumstats,
          fmt=None,
          tab_fmt="tsv",
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
          maf=None,
          n=None,
          beta=None,
          se=None,
          chisq=None,
          z=None,
          f=None,
          t=None,
          p=None,
          q=None,
          mlog10p=None,
          test=None,
          info=None,
          OR=None,
          OR_95L=None,
          OR_95U=None,
          beta_95L=None,
          beta_95U=None,
          HR=None,
          HR_95L=None,
          HR_95U=None,
          i2=None,
          snpr2=None,
          phet=None,
          dof=None,
          ncase=None,
          ncontrol=None,
          neff=None,
          direction=None,
          status=None,
          study=None,
          trait=None,
          build=None,
          other=None,
          exclude=None,
          include=None,
          chrom_pat=None,
          snpid_pat=None,
          verbose=False,
          log=None,
          readargs=None,
          **kwreadargs):
    """
    Load and preformat summary statistics data into standardized GWASLab format.

    Parameters
    ----------
    sumstats : str or pandas.DataFrame
        Input summary statistics, provided either as a file path or a DataFrame.
        When summary statistics are split by chromosome, a single path pattern
        may be supplied using the `@` symbol as a placeholder for the chromosome
        number. For example: "gwas/chr@.sumstats.gz" will load
        "gwas/chr1.sumstats.gz", "gwas/chr2.sumstats.gz", ... automatically.
    fmt : str, optional
        Format name to get predefined mapping if provided. (e.g., 'gwaslab', 'vcf').
    tab_fmt : str, default: 'tsv'
        Table format ('tsv', 'parquet').
    snpid : str, optional
        Column name for SNP identifiers in the input data. Expected format is CHR:POS:NEA:EA (e.g., 1:123:A:G), although the delimiter may vary depending on the source.
    rsid : str, optional
        Column name for rsID in the input data. Values should follow the standard "rs" prefix (e.g., rs12345).
    chrom : str, optional
        Column name for chromosome in input data.
    pos : str, optional
        Column name for position in input data.
    ea : str, optional
        Column name for effect allele in input data. (assuming alternative allele)
    nea : str, optional
        Column name for non-effect allele in input data. (assuming reference allele)
    eaf : str, optional
        Column name for effect allele frequency in input data.
    neaf : str, optional
        Column name for non-effect allele frequency in input data.
    maf : str, optional
        Column name for minor allele frequency in input data.
    n : str or int, optional
        Column name or constant value for sample size.
    beta : str, optional
        Column name for beta in input data.
    se : str, optional
        Column name for standard error in input data.
    chisq : str, optional
        Column name for chi-square in input data.
    z : str, optional
        Column name for z-score in input data.
    f : str, optional
        Column name for F-statistic in input data.
    t : str, optional
        Column name for T-statistic in input data.
    p : str, optional
        Column name for p-value in input data.
    q : str, optional
        Column name for Q-statistic in input data.
    mlog10p : str, optional
        Column name for -log10(p) in input data.
    test : str, optional
        Column name for test type in input data.
    info : str, optional
        Column name for imputation info in input data.
    OR : str, optional
        Column name for odds ratio in input data.
    OR_95L : str, optional
        Column name for lower 95% CI of OR in input data.
    OR_95U : str, optional
        Column name for upper 95% CI of OR in input data.
    beta_95L : str, optional
        Column name for lower 95% CI of beta in input data.
    beta_95U : str, optional
        Column name for upper 95% CI of beta in input data.
    HR : str, optional
        Column name for hazard ratio in input data.
    HR_95L : str, optional
        Column name for lower 95% CI of HR in input data.
    HR_95U : str, optional
        Column name for upper 95% CI of HR in input data.
    i2 : str, optional
        Column name for I2 statistic in input data.
    snpr2 : str, optional
        Column name for SNP R2 in input data.
    phet : str, optional
        Column name for p-heterogeneity in input data.
    dof : str, optional
        Column name for degrees of freedom in input data.
    ncase : str or int, optional
        Column name or constant value for case count.
    ncontrol : str or int, optional
        Column name or constant value for control count.
    neff : str or int, optional
        Column name or constant value for effective sample size.
    direction : str, optional
        Column name for meta-analysis effect-direction strings, where each character
        ("+", "-", "0") represents the direction of effect for one cohort (e.g., "++-0+").
    status : str, optional
        Column name for status in input data.
    study : str, optional, default="Study_1"
        Column name for study ID in input data.
    trait : str, optional, default="Trait_1"
        Column name for trait in input data.
    build : str, optional
        Genome build version (e.g., '19' and '38').
    species : str, default="homo sapiens"
        species
    other : list, optional
        Additional columns in the raw file to load.
    chrom_pat : str, optional
        Regex pattern to filter chromosomes like'chrX'
    snpid_pat : str, optional
        Regex pattern to filter SNPs based on snpid like'chrX:'.
    verbose : bool, default: False
        Enable verbose output.
    readargs : dict, optional
        Additional arguments for reading files using pd.read_csv() like `nrows`, `comment`. 
        Example, {"nrows": 1000} means to load first 1000 rows.

    Returns
    -------
    pandas.DataFrame
        Formatted summary statistics with standardized column names.
    dict
        Updated readargs dictionary.

    Raises
    ------
    ValueError
        If input is not a path or DataFrame, or if columns are missing.

    Less used parameters
    -------------------------
    exclude : list, optional
        Columns to exclude explicitly. Columns should be passed as GWASLab built-in HEADER keywords in uppercase like BETA, DIRECTION. Not original headers.
    include : list, optional
        Columns to include explicitly. Columns should be passed as GWASLab built-in HEADER keywords in uppercase like BETA, DIRECTION. Not original headers.
    """
    
    if readargs is None:
        readargs = dict()

    readargs = readargs | kwreadargs
    
    if log is None:
        log = Log()
    rename_dictionary = {}
    usecols = list()
    if other is None:
        other = list()
    if exclude is None:
        exclude = list()
    if include is None:
        include = list()
    
    dtype_dictionary = {}    
    if readargs is None:
        readargs={}
 #######################################################################################################################################################
    # workflow: 
    # 1. formatbook
    # 2. user specified header
    # 3. include & exclude
    if tab_fmt=="parquet":
        if type(sumstats) is str:
            log.write("Start to load data from parquet file....",verbose=verbose)
            log.write(" -path: {}".format(sumstats),verbose=verbose)
            sumstats = pd.read_parquet(sumstats,**readargs) 
            log.write("Finished loading parquet file into pd.DataFrame....",verbose=verbose)
        else:
            raise ValueError("Please input a path for parquet file.")
    
    if fmt is not None:
        # loading format parameters
        log.write("Start to load format from formatbook....",verbose=verbose)
        
        # load format data
        meta_data,rename_dictionary = get_format_dict(fmt)
        
        ########## print format information################################################
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log)
        
        if "format_separator" in meta_data.keys():
            if "sep" not in readargs.keys():
                readargs["sep"] = meta_data["format_separator"]
            else:
                if readargs["sep"] != meta_data["format_separator"]:
                    log.write('  - format_separator will be changed to: "{}"'.format(readargs["sep"]),verbose=verbose)

        if "format_na" in meta_data.keys():
            readargs["na_values"] = meta_data["format_na"]
        
        if "format_comment" in meta_data.keys():
            readargs["comment"] = meta_data["format_comment"]
        
        if "format_other_cols" in meta_data.keys():
            other += meta_data["format_other_cols"]
        
        if "sep" not in readargs.keys():
            readargs["sep"] = "\t"
    else:
        meta_data = None
        
#########################################################################################################################################################      
    
    # check chr-separated path / vcf / then print header.     
    inpath, inpath_chr_list, inpath_chr_num_list, format_cols, raw_cols, usecols, dtype_dictionary  = check_path_and_header(sumstats, 
                                                                                                                    fmt, 
                                                                                                                    meta_data, 
                                                                                                                    readargs, 
                                                                                                                    usecols, 
                                                                                                                    dtype_dictionary, 
                                                                                                                    rename_dictionary, 
                                                                                                                    log, 
                                                                                                                    verbose)     

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
        rename_dictionary[ref]= "REF"
        dtype_dictionary[ref]="string"
    if alt:
        usecols.append(alt)
        rename_dictionary[alt]= "ALT"
        dtype_dictionary[alt]="string"
    if eaf:
        usecols.append(eaf)
        rename_dictionary[eaf]= "EAF"
    elif neaf:
        # neaf will be converted to eaf
        usecols.append(neaf)
        rename_dictionary[neaf]= "EAF"
    if maf:
        usecols.append(maf)
        rename_dictionary[maf]= "MAF"
    if n and (type(n) is str):
        usecols.append(n)
        rename_dictionary[n]= "N"
    if ncase and (type(ncase) is str):
        usecols.append(ncase)
        rename_dictionary[ncase]= "N_CASE"    
    if ncontrol and (type(ncontrol) is str):
        usecols.append(ncontrol)
        rename_dictionary[ncontrol]= "N_CONTROL"    
    if neff and (type(neff) is str):
        usecols.append(neff)
        rename_dictionary[neff]= "N_EFF" 
    if beta:
        usecols.append(beta)
        rename_dictionary[beta]= "BETA"
    if beta_95L:
        usecols.append(beta_95L)
        rename_dictionary[beta_95L]= "BETA_95L"
    if beta_95U:
        usecols.append(beta_95U)
        rename_dictionary[beta_95U]= "BETA_95U"  
    if se:
        usecols.append(se)
        rename_dictionary[se]= "SE"
    if chisq:
        usecols.append(chisq)
        rename_dictionary[chisq]="CHISQ"
    if z:
        usecols.append(z)
        rename_dictionary[z]= "Z"
    if q:
        usecols.append(q)
        rename_dictionary[q]= "Q"
    if p:
        usecols.append(p)
        rename_dictionary[p]= "P"
    if t:
        usecols.append(t)
        rename_dictionary[t]= "T"
    if f:
        usecols.append(f)
        rename_dictionary[f]= "F"
    if mlog10p:
        usecols.append(mlog10p)
        rename_dictionary[mlog10p]= "MLOG10P"
    if test:
        usecols.append(test)
        rename_dictionary[test]= "TEST"        
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
    if HR:
        usecols.append(HR)
        rename_dictionary[HR]= "HR"
    if HR_95L:
        usecols.append(HR_95L)
        rename_dictionary[HR_95L]= "HR_95L"
    if HR_95U:
        usecols.append(HR_95U)
        rename_dictionary[HR_95U]= "HR_95U"        
    if phet:
        usecols.append(phet)
        rename_dictionary[phet]= "P_HET"      
    if i2:
        usecols.append(i2)
        rename_dictionary[i2]= "I2"    
    if snpr2:
        usecols.append(snpr2)
        rename_dictionary[snpr2]= "SNPR2"    
    if dof:
        usecols.append(dof)
        rename_dictionary[dof]= "DOF"    
    if direction:
        usecols.append(direction)
        rename_dictionary[direction]="DIRECTION"
    if status:
        usecols.append(status)
        rename_dictionary[status]="STATUS"
        dtype_dictionary[status]="string"
    if other:
        overlapped = _check_overlap_with_reserved_keys(other)
        log.warning("Columns with headers overlapping with GWASLab reserved keywords:{}".format(overlapped),verbose=verbose)
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

    if len(include)>0:
    # extract only specified keys
        usecols_new =[]
        for i in include:
            # rename_dictionary: sumstats to gwaslab
            for k, v in rename_dictionary.items():
                if i == v:
                    # get list of sumstats header
                    usecols_new.append(k)
        usecols_valid =[]
        for i in usecols_new:
            if i in usecols:
                usecols_valid.append(i)
        log.write(f' -Include columns :{",".join(usecols_valid)}' ,verbose=verbose)
        usecols = usecols_valid

    if len(exclude)>0:
    # exclude specified keys
        exclude_cols =[]
        for i in exclude:
            # rename_dictionary: sumstats to gwaslab
            for k, v in rename_dictionary.items():
                if i == v:
                    # get list of sumstats header
                    exclude_cols.append(k)
        log.write(f' -Exclude columns :{",".join(exclude_cols)}' ,verbose=verbose)
        for i in exclude_cols:
            if i in usecols:
                usecols.remove(i)

 #loading data ##########################################################################################################
    
    try:
        if type(sumstats) is str:
            ## loading data from path
            #inpath = sumstats
            if "@" in inpath:
                log.write("Start to initialize gl.Sumstats from files with pattern :" + inpath,verbose=verbose)
                sumstats_chr_list=[]
                for i in inpath_chr_list:
                    log.write(" -Loading:" + i)
                    skip_rows = get_skip_rows(i)
                    readargs["skiprows"] = skip_rows
                    explicit = {"usecols","dtype_dictionary"}
                    readargs = {k: v for k, v in readargs.items() if k not in explicit}
                    sumstats_chr = pd.read_table(i,
                                        usecols=set(usecols),
                                        dtype=dtype_dictionary,
                                        **readargs)
                    sumstats_chr_list.append(sumstats_chr)
                log.write(" -Merging sumstats for chromosomes:",",".join(inpath_chr_num_list),verbose=verbose)
                sumstats = pd.concat(sumstats_chr_list, axis=0, ignore_index=True) 
                del(sumstats_chr_list)
                gc.collect()
            else:
                skip_rows = get_skip_rows(inpath)
                readargs["skiprows"] = skip_rows
                log.write("Start to initialize gl.Sumstats from file :" + inpath,verbose=verbose)
                if chrom_pat is not None:
                    sumstats = _load_single_chr(inpath,
                                                usecols,
                                                dtype_dictionary,
                                                readargs=readargs,
                                                rename_dictionary=rename_dictionary,
                                                chrom_pat=chrom_pat,
                                                log=log,
                                                verbose=verbose)
                elif snpid_pat is not None: 
                                
                    sumstats = _load_variants_with_pattern(inpath,
                                                usecols,
                                                dtype_dictionary,
                                                readargs=readargs,
                                                rename_dictionary=rename_dictionary,
                                                snpid_pat=snpid_pat,
                                                log=log,
                                                verbose=verbose)
                else:
                    explicit = {"usecols","dtype_dictionary"}
                    readargs = {k: v for k, v in readargs.items() if k not in explicit}
                    sumstats = pd.read_table(inpath,
                                    usecols=set(usecols),
                                    dtype=dtype_dictionary,
                                    **readargs)

        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            log.write("Start to initialize gl.Sumstats from pandas DataFrame ...",verbose=verbose)
            sumstats = sumstats[usecols].copy()
            for key,value in dtype_dictionary.items():
                if key in usecols:
                    astype = value
                    if rename_dictionary[key]=="CHR":
                        astype ="Int64"  
                    try:
                        sumstats[key] = sumstats[key].astype(astype)
                    except:
                        sumstats[key] = sumstats[key].astype("string")
    except ValueError:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure it contain the columns.")

    ## renaming columns ###############################################################################################
    if fmt == "vcf":
        sumstats = parse_vcf_study(sumstats,format_cols,study,vcf_usecols,log=log,verbose=verbose)
        usecols = vcf_usecols
    
    converted_columns = list(map(lambda x: rename_dictionary[x], set(usecols)))
    
    ## renaming log
    log.write(" -Reading columns          :", ",".join(set(usecols)),verbose=verbose)
    log.write(" -Renaming columns to      :", ",".join(converted_columns),verbose=verbose)
    log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns),verbose=verbose) 
    
    ## renaming  #####################################################################################
    sumstats = sumstats.rename(columns=rename_dictionary)

    ## if n was provided as int #####################################################################################
    if type(n) is int:
        sumstats["N"] = n
    if type(ncase) is int:
        sumstats["N_CASE"] = ncase
    if type(ncontrol) is int:
        sumstats["N_CONTROL"] = ncontrol
    
    ### status ######################################################################################################
    
    sumstats = process_status(sumstats=sumstats,build=build,status=status,log=log,verbose=verbose)
    
    ## ea/nea, ref/alt ##############################################################################################
    sumstats = process_allele(sumstats=sumstats,log=log,rename_dictionary=rename_dictionary,verbose=verbose)
        
    ## NEAF to EAF ###########################################################################################################
    if neaf is not None or ("NEAF" in sumstats.columns and "EAF" not in sumstats.columns):
        sumstats = process_neaf(sumstats=sumstats,log=log,verbose=verbose)

    ## reodering ###################################################################################################  
    sumstats = sortcolumn(sumstats=sumstats,log=log,verbose=verbose)    
    sumstats = quick_convert_datatype(sumstats,log=log,verbose=verbose)
    
    check_datatype(sumstats,log=log,verbose=verbose)
    gc.collect()
    check_dataframe_memory_usage(sumstats,log=log,verbose=verbose)
    
    log.write("Finished loading data successfully!",verbose=verbose)
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
    log.write(" -Parsing based on FORMAT: ", format_cols,verbose=verbose)
    log.write(" -Parsing vcf study : ", study,verbose=verbose)
    sumstats[format_cols] = sumstats[study].str.split(":",expand=True).values
    sumstats = sumstats.drop(["FORMAT",study],axis=1)
    sumstats = sumstats[ vcf_usecols]
    gc.collect()
    return sumstats

def print_format_info(fmt,meta_data, rename_dictionary, verbose, log,output=False, skip_meta_records=None):
    log.write(" -"+fmt+" format meta info:",verbose=verbose)   
    if skip_meta_records is None:
        skip_meta_records =[]
    for key,value in meta_data.items():
        if key in skip_meta_records:
            continue
        if value is None:
            continue
        if type(value) is str:
            if "\n" in value:
                value_first_line=value.split("\n")[0]
                log.write("  -",key," : "+value_first_line.strip()+"...",verbose=verbose)
            elif value==" ":
                log.write('  -',key,' : \\s ',verbose=verbose)     
            elif value=="\t":
                log.write('  -',key,' : \\t',verbose=verbose)    
            else:
                log.write("  -",key," : "+value.strip(),verbose=verbose)  
        elif type(value) is list:
            log.write("  -",key," : "+','.join(value),verbose=verbose)  
        else:
            log.write("  -",key," : ",value,verbose=verbose)  
    keys=[]
    values=[]
    for key,value in rename_dictionary.items():
        keys.append(key)
        values.append(value)
    if fmt!="gwaslab":
        if output == False:
            if fmt!="auto":
                log.write(" -"+fmt+" to gwaslab format dictionary:",verbose=verbose)  
                log.write("  - "+fmt+" keys:",",".join(keys),verbose=verbose) 
                log.write("  - gwaslab values:",",".join(values),verbose=verbose) 
            else:
                log.write("  - Auto-detection mode. Note: auto-detection assumes A1=EA; Alt=EA and Frq=EAF...",verbose=verbose)
                log.write("  - Header conversion source: https://github.com/Cloufield/formatbook/blob/main/formats/auto.json",verbose=verbose)  
        else:
            log.write(" -gwaslab to "+fmt+" format dictionary:",verbose=verbose)  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            log.write("  - gwaslab keys:",  ','.join(keys),verbose=verbose) 
            log.write("  - "+fmt+" values:"  , ','.join(values),verbose=verbose)       

def process_neaf(sumstats,log,verbose):
    log.write(" -NEAF is specified...",verbose=verbose) 
    pre_number=len(sumstats)
    log.write(" -Checking if 0<= NEAF <=1 ...",verbose=verbose) 
    if "NEAF" in sumstats.columns:
        sumstats["NEAF"] = pd.to_numeric(sumstats["NEAF"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["NEAF"]>=0) & (sumstats["NEAF"]<=1),:]
        sumstats["EAF"] = 1- sumstats["NEAF"]
        sumstats.drop(columns=["NEAF"], inplace=True)
    else:
        sumstats["EAF"] = pd.to_numeric(sumstats["EAF"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["EAF"]>=0) & (sumstats["EAF"]<=1),:]
        sumstats["EAF"] = 1- sumstats["EAF"]
    log.write(" -Converted NEAF to EAF.",verbose=verbose) 
    after_number=len(sumstats)
    log.write(" -Removed "+str(pre_number - after_number)+" variants with bad NEAF.",verbose=verbose) 
    return sumstats

def process_allele(sumstats,
                   log,
                   rename_dictionary,
                   verbose):
    
    if "EA" in sumstats.columns:

        if "REF" in sumstats.columns and "ALT" in sumstats.columns:

            if "NEA" not in sumstats.columns:
                log.write(" NEA not available: assigning REF to NEA...",verbose=verbose) 
                sumstats["NEA"]=sumstats["REF"]    
            
            log.write(" -EA,REF and ALT columns are available: assigning NEA...",verbose=verbose) 
            ea_alt = sumstats["EA"]==sumstats["ALT"]
            
            log.write(" -For variants with EA == ALT : assigning REF to NEA ...",verbose=verbose) 
            sumstats.loc[ea_alt,"NEA"] = sumstats.loc[ea_alt,"REF"]
            
            ea_not_alt = sumstats["EA"]!=sumstats["ALT"]
            log.write(" -For variants with EA != ALT : assigning ALT to NEA ...",verbose=verbose) 
            sumstats.loc[ea_not_alt,"NEA"] = sumstats.loc[ea_not_alt,"ALT"]

            #sumstats = sumstats.drop(labels=["REF","ALT"],axis=1)
            sumstats["REF"]=sumstats["REF"].astype("category") 
            sumstats["ALT"]=sumstats["ALT"].astype("category") 
            
        sumstats["EA"]=sumstats["EA"].astype("category")     
    if "NEA" in sumstats.columns:
        sumstats["NEA"]=sumstats["NEA"].astype("category")  

    if "NEA" not in sumstats.columns and "EA" not in sumstats.columns:
        if "REF" in sumstats.columns and "ALT" in sumstats.columns:
            log.write(" -Converting REF and ALT to NEA and EA ...")
            sumstats["REF"]=sumstats["REF"].astype("category") 
            sumstats["ALT"]=sumstats["ALT"].astype("category") 
            sumstats["NEA"]=sumstats["REF"].astype("category") 
            sumstats["EA"]=sumstats["ALT"].astype("category") 
            sumstats = sumstats.drop(columns=["ALT","REF"])
    return sumstats

def process_status(sumstats,build,status, log,verbose):
    if status is None:
        log.write(" -Initiating a status column: STATUS ...",verbose=verbose)
        #sumstats["STATUS"] = int(build)*(10**5) +99999
        build = _process_build(build,log,verbose)
        sumstats["STATUS"] = build +"99999"

    sumstats["STATUS"] = pd.Categorical(sumstats["STATUS"],categories=STATUS_CATEGORIES)
    return sumstats


def _load_single_chr(inpath,usecols,dtype_dictionary,readargs,rename_dictionary,chrom_pat,log,verbose):
    explicit = {"usecols","dtype_dictionary","chunksize","iterator"}
    readargs = {k: v for k, v in readargs.items() if k not in explicit}
    sumstats_iter = pd.read_table(inpath,
                usecols=set(usecols),
                dtype=dtype_dictionary, 
                iterator=True, 
                chunksize=500000,
                **readargs)
    # get chr 
    for k,v in rename_dictionary.items():
        if v=="CHR":
            if k in usecols:
                log.write(" -Columns used to filter variants: {}".format(k),verbose=verbose)
                chunk_chrom = k
                break

    log.write(" -Loading only variants on chromosome with pattern : {} ...".format(chrom_pat),verbose=verbose)
    sumstats_filtered = pd.concat([chunk[chunk[chunk_chrom].str.match(chrom_pat, case=False,na=False) ] for chunk in sumstats_iter])
    log.write(" -Loaded {} variants on chromosome with pattern :{} ...".format(len(sumstats_filtered), chrom_pat),verbose=verbose)
    return sumstats_filtered

def _load_variants_with_pattern(inpath,usecols,dtype_dictionary,readargs,rename_dictionary,snpid_pat,log,verbose):
    explicit = {"usecols","dtype_dictionary","chunksize","iterator"}
    readargs = {k: v for k, v in readargs.items() if k not in explicit}
    sumstats_iter = pd.read_table(inpath,
                usecols=set(usecols),
                dtype=dtype_dictionary, 
                iterator=True, 
                chunksize=500000,
                **readargs)
    # get chr 
    for k,v in rename_dictionary.items():
        if v=="SNPID":
            if k in usecols:
                log.write(" -Columns used to filter variants: {}".format(k),verbose=verbose)
                chunk_snpid = k
                break

    log.write(" -Loading only variants with pattern :  {} ...".format(snpid_pat),verbose=verbose)
    sumstats_filtered = pd.concat([chunk[chunk[chunk_snpid].str.match(snpid_pat, case=False,na=False) ] for chunk in sumstats_iter])
    log.write(" -Loaded {} variants with pattern : {} ...".format(len(sumstats_filtered), snpid_pat),verbose=verbose)
    return sumstats_filtered


def check_path_and_header(sumstats=None, 
                          fmt=None, 
                          meta_data=None, 
                          readargs=None, 
                          usecols=None, 
                          dtype_dictionary=None, 
                          rename_dictionary=None, 
                          log=None, 
                          verbose=None):
    

    if type(sumstats) is str:
        ## loading data from path #################################################
        inpath = sumstats
        
        try:
            format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list = process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose)
        
        except (FileNotFoundError, IndexError):
            log.warning("Loading {} failed...Tesing if compressed/uncompressed...".format(inpath),verbose=verbose)
            try:
                if inpath[-3:]==".gz":
                    inpath = inpath[:-3]
                    log.write(" -Trying to load {}...".format(inpath),verbose=verbose)
                    format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list =process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose)
                else:
                    inpath = inpath+".gz"
                    log.write(" -Trying to load {}...".format(inpath),verbose=verbose)
                    format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list = process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose)
            except:
                raise ValueError("Please input a valid path, and make sure the separator is correct and the columns you specified are in the file.") 

        ###################################################################################### 
    elif type(sumstats) is pd.DataFrame:
        inpath = None
        format_cols = None
        inpath_chr_list = None
        inpath_chr_num_list = None
        ## loading data from dataframe
        raw_cols = sumstats.columns

    ################################################
    for key,value in rename_dictionary.items():
        # check avaiable keys  key->raw header
        # usecols : a list of raw headers to load from file/DataFrame 
        if key in raw_cols:
            usecols.append(key)
        if value in ["EA","NEA"]:
            dtype_dictionary[key]="category"
        if value in ["STATUS"]:
            dtype_dictionary[key]="string"     
        if value in ["CHR"]:
            dtype_dictionary[key]="string"  

    return inpath, inpath_chr_list, inpath_chr_num_list, format_cols, raw_cols, usecols, dtype_dictionary

def process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose):
    
    format_cols = None
    inpath_chr_list = None
    inpath_chr_num_list = None

    if "@" in inpath:
        log.write(" -Detected @ in path: load sumstats by each chromosome...",verbose=verbose)
        inpath_chr_list=[]
        inpath_chr_num_list=[]
        
        # create a regex pattern for matching
        pat = os.path.basename(inpath).replace("@","(\w+)")
        
        # get dir
        dirname = os.path.dirname(inpath)
        
        # all files in the directory
        files = os.listdir(dirname)
        
        files.sort()

        for file in files:
            # match
            result = re.match(pat, file)
            if result:
                # get chr
                chr_matched = str(result.group(1))
                inpath_chr_num_list.append(chr_matched)
                inpath_chr_list.append(inpath.replace("@",str(chr_matched))  )
        
        log.write(" -Chromosomes detected:",",".join(inpath_chr_num_list),verbose=verbose)

        #if inpath_chr_list is empty-> IndexError
        readargs_header = get_readargs_header(inpath = inpath_chr_list[0], readargs = readargs)
        row_one = pd.read_table(inpath_chr_list[0],**readargs_header)
        # columns in the sumstats
        raw_cols = row_one.columns
    else:
    ##### loading data from tabular file#################################################
    #if file not found, FileNotFoundError
        readargs_header = get_readargs_header(inpath = inpath, readargs = readargs)
        row_one = pd.read_table(inpath,**readargs_header)
        raw_cols = row_one.columns

    if fmt=="vcf":
        # expanded
        format_cols = list(row_one["FORMAT"].str.split(":"))[0]
        # fixed + study1 + expanded
        raw_cols = meta_data["format_fixed"] + [raw_cols[9]] + format_cols

    return format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list
