import pandas as pd
import polars as pl
import numpy as np
import scipy.stats as ss
import gzip
import os
import re
import gc
from gwaslab.bd.bd_common_data import get_format_dict
from gwaslab.qc.qc_fix_sumstats import _process_build
from gwaslab.qc.qc_check_datatype_polars import check_datatype_polars
from gwaslab.qc.qc_check_datatype_polars import quick_convert_datatype
from gwaslab.qc.qc_check_datatype_polars import check_dataframe_memory_usage
from gwaslab.qc.qc_reserved_headers import _check_overlap_with_reserved_keys
from gwaslab.info.g_vchange_status import STATUS_CATEGORIES
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_reserved_headers import DEFAULT_COLUMN_ORDER
#20221030
def preformatp(sumstats,
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
 #######################################################################################################################################################
    # workflow: 
    # 1. formatbook
    # 2. user specified header
    # 3. include & exclude
    if tab_fmt=="parquet":
        if type(sumstats) is str:
            log.write("Start to load data from parquet file....",verbose=verbose)
            log.write(" -path: {}".format(sumstats),verbose=verbose)
            sumstats = pl.read_parquet(sumstats,**readargs)
            log.write("Finished loading parquet file into DataFrame....",verbose=verbose)
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
            if "separator" not in readargs.keys():
                readargs["separator"] = meta_data["format_separator"]
            else:
                if readargs["separator"] != meta_data["format_separator"]:
                    log.write('  - format_separator will be changed to: "{}"'.format(readargs["separator"]),verbose=verbose)

        if "format_na" in meta_data.keys():
            readargs["null_values"] = meta_data["format_na"]
        
        if "format_comment" in meta_data.keys():
            readargs["comment_prefix"] = meta_data["format_comment"]
        
        if "format_other_cols" in meta_data.keys():
            other += meta_data["format_other_cols"]
        
        if "separator" not in readargs.keys():
            readargs["separator"] = "\t"
    else:
        meta_data = None
        
        if "separator" not in readargs.keys():
            readargs["separator"] = "\t"
        
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
    if snpid and snpid not in rename_dictionary:
        usecols.append(snpid)
        rename_dictionary[snpid]= "SNPID"
    if rsid and rsid not in rename_dictionary:
        usecols.append(rsid)
        rename_dictionary[rsid]= "rsID"
    if chrom and chrom not in rename_dictionary:
        usecols.append(chrom)
        rename_dictionary[chrom]= "CHR"
        dtype_dictionary[chrom]=pl.String()
    if pos and pos not in rename_dictionary:
        usecols.append(pos)
        rename_dictionary[pos]= "POS"
    if ea and ea not in rename_dictionary:
        usecols.append(ea)
        rename_dictionary[ea]= "EA"
        dtype_dictionary[ea]=pl.String()
    if nea and nea not in rename_dictionary:
        usecols.append(nea)
        rename_dictionary[nea]= "NEA"
        dtype_dictionary[nea]=pl.String()
    if ref and ref not in rename_dictionary:
        usecols.append(ref)
        rename_dictionary[ref]= "REF"
        dtype_dictionary[ref]=pl.String()
    if alt and alt not in rename_dictionary:
        usecols.append(alt)
        rename_dictionary[alt]= "ALT"
        dtype_dictionary[alt]=pl.String()
    if eaf and eaf not in rename_dictionary:
        usecols.append(eaf)
        rename_dictionary[eaf]= "EAF"
    elif neaf and neaf not in rename_dictionary:
        # neaf will be converted to eaf
        usecols.append(neaf)
        rename_dictionary[neaf]= "EAF"
    if maf and maf not in rename_dictionary:
        usecols.append(maf)
        rename_dictionary[maf]= "MAF"
    if n and (type(n) is str) and n not in rename_dictionary:
        usecols.append(n)
        rename_dictionary[n]= "N"
    if ncase and (type(ncase) is str) and ncase not in rename_dictionary:
        usecols.append(ncase)
        rename_dictionary[ncase]= "N_CASE"    
    if ncontrol and (type(ncontrol) is str) and ncontrol not in rename_dictionary:
        usecols.append(ncontrol)
        rename_dictionary[ncontrol]= "N_CONTROL"    
    if neff and (type(neff) is str) and neff not in rename_dictionary:
        usecols.append(neff)
        rename_dictionary[neff]= "N_EFF" 
    if beta and beta not in rename_dictionary:
        usecols.append(beta)
        rename_dictionary[beta]= "BETA"
    if beta_95L and beta_95L not in rename_dictionary:
        usecols.append(beta_95L)
        rename_dictionary[beta_95L]= "BETA_95L"
    if beta_95U and beta_95U not in rename_dictionary:
        usecols.append(beta_95U)
        rename_dictionary[beta_95U]= "BETA_95U"  
    if se and se not in rename_dictionary:
        usecols.append(se)
        rename_dictionary[se]= "SE"
    if chisq and chisq not in rename_dictionary:
        usecols.append(chisq)
        rename_dictionary[chisq]="CHISQ"
    if z and z not in rename_dictionary:
        usecols.append(z)
        rename_dictionary[z]= "Z"
    if q and q not in rename_dictionary:
        usecols.append(q)
        rename_dictionary[q]= "Q"
    if p and p not in rename_dictionary:
        usecols.append(p)
        rename_dictionary[p]= "P"
    if t and t not in rename_dictionary:
        usecols.append(t)
        rename_dictionary[t]= "T"
    if f and f not in rename_dictionary:
        usecols.append(f)
        rename_dictionary[f]= "F"
    if mlog10p and mlog10p not in rename_dictionary:
        usecols.append(mlog10p)
        rename_dictionary[mlog10p]= "MLOG10P"
    if test and test not in rename_dictionary:
        usecols.append(test)
        rename_dictionary[test]= "TEST"        
    if info and info not in rename_dictionary:
        usecols.append(info)
        rename_dictionary[info]= "INFO"
    if OR and OR not in rename_dictionary:
        usecols.append(OR)
        rename_dictionary[OR]= "OR"
    if OR_95L and OR_95L not in rename_dictionary:
        usecols.append(OR_95L)
        rename_dictionary[OR_95L]= "OR_95L"
    if OR_95U and OR_95U not in rename_dictionary:
        usecols.append(OR_95U)
        rename_dictionary[OR_95U]= "OR_95U"
    if HR and HR not in rename_dictionary:
        usecols.append(HR)
        rename_dictionary[HR]= "HR"
    if HR_95L and HR_95L not in rename_dictionary:
        usecols.append(HR_95L)
        rename_dictionary[HR_95L]= "HR_95L"
    if HR_95U and HR_95U not in rename_dictionary:
        usecols.append(HR_95U)
        rename_dictionary[HR_95U]= "HR_95U"        
    if phet and phet not in rename_dictionary:
        usecols.append(phet)
        rename_dictionary[phet]= "P_HET"      
    if i2 and i2 not in rename_dictionary:
        usecols.append(i2)
        rename_dictionary[i2]= "I2"    
    if snpr2 and snpr2 not in rename_dictionary:
        usecols.append(snpr2)
        rename_dictionary[snpr2]= "SNPR2"    
    if dof and dof not in rename_dictionary:
        usecols.append(dof)
        rename_dictionary[dof]= "DOF"    
    if direction and direction not in rename_dictionary:
        usecols.append(direction)
        rename_dictionary[direction]="DIRECTION"
    if status and status not in rename_dictionary:
        usecols.append(status)
        rename_dictionary[status]="STATUS"
        dtype_dictionary[status]=pl.Int64()
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

    # Remove duplicates from usecols while preserving order
    seen = set()
    usecols = [x for x in usecols if not (x in seen or seen.add(x))]

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
                    readargs["skip_rows"] = skip_rows
                    explicit = {"columns","schema_overrides"}
                    readargs_clean = {k: v for k, v in readargs.items() if k not in explicit}
                    sumstats_chr = pl.read_csv(i,
                                        columns=usecols,
                                        schema_overrides=dtype_dictionary,
                                        **readargs_clean)
                    sumstats_chr_list.append(sumstats_chr)
                log.write(" -Merging sumstats for chromosomes:",",".join(inpath_chr_num_list),verbose=verbose)
                sumstats = pl.concat(sumstats_chr_list, rechunk=True) 
                del(sumstats_chr_list)
                gc.collect()
            else:
                skip_rows = get_skip_rows(inpath)
                readargs["skip_rows"] = skip_rows
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
                    explicit = {"columns","schema_overrides"}
                    readargs_clean = {k: v for k, v in readargs.items() if k not in explicit}
                    sumstats = pl.read_csv(inpath,
                                    columns=usecols,
                                    schema_overrides=dtype_dictionary,
                                    **readargs_clean)

        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            log.write("Start to initialize gl.Sumstats from pandas DataFrame ...",verbose=verbose)
            sumstats = sumstats.copy()
            for key,value in dtype_dictionary.items():
                if key in sumstats.columns:
                    astype = value
                    if key in rename_dictionary and rename_dictionary[key]=="CHR":
                        astype ="Int64"  
                    try:
                        # Convert pandas types - for polars String, use pandas string
                        if isinstance(value, pl.String):
                            sumstats[key] = sumstats[key].astype("string")
                        else:
                            # For other types, try to convert
                            sumstats[key] = sumstats[key].astype(astype)
                    except:
                        sumstats[key] = sumstats[key].astype("string")
            sumstats = pl.from_pandas(sumstats)
            # Get unique columns that exist in the dataframe
            available_cols = [c for c in usecols if c in sumstats.columns]
            # Remove duplicates
            seen = set()
            available_cols = [x for x in available_cols if not (x in seen or seen.add(x))]
            sumstats = sumstats.select([pl.col(c) for c in available_cols])
            if len(dtype_dictionary)>0:
                sumstats = sumstats.with_columns([
                    pl.col(k).cast(dtype_dictionary[k], strict=False) for k in dtype_dictionary.keys() if k in sumstats.columns
                ])
        elif isinstance(sumstats, pl.DataFrame):
            ## loading data from polars dataframe
            log.write("Start to initialize gl.Sumstats from polars DataFrame ...",verbose=verbose)
            # Get unique columns that exist in the dataframe
            available_cols = [c for c in usecols if c in sumstats.columns]
            # Remove duplicates
            seen = set()
            available_cols = [x for x in available_cols if not (x in seen or seen.add(x))]
            sumstats = sumstats.select([pl.col(c) for c in available_cols])
            if len(dtype_dictionary)>0:
                sumstats = sumstats.with_columns([
                    pl.col(k).cast(dtype_dictionary[k], strict=False) for k in dtype_dictionary.keys() if k in sumstats.columns
                ])
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
    sumstats = sumstats.rename(rename_dictionary)

    ## if n was provided as int #####################################################################################
    if type(n) is int:
        sumstats = sumstats.with_columns(pl.lit(n).alias("N"))
    if type(ncase) is int:
        sumstats = sumstats.with_columns(pl.lit(ncase).alias("N_CASE"))
    if type(ncontrol) is int:
        sumstats = sumstats.with_columns(pl.lit(ncontrol).alias("N_CONTROL"))
    
    ### status ######################################################################################################
    
    sumstats = process_status(sumstats=sumstats,build=build,status=status,log=log,verbose=verbose)
    
    ## ea/nea, ref/alt ##############################################################################################
    sumstats = process_allele(sumstats=sumstats,log=log,rename_dictionary=rename_dictionary,verbose=verbose)
        
    ## NEAF to EAF ###########################################################################################################
    if neaf is not None or ("NEAF" in sumstats.columns and "EAF" not in sumstats.columns):
        sumstats = process_neaf(sumstats=sumstats,log=log,verbose=verbose)

    ## reodering ###################################################################################################  
    sumstats = _sort_column(sumstats=sumstats,log=log,verbose=verbose)    
    sumstats = quick_convert_datatype(sumstats,log=log,verbose=verbose)
    
    # Force create IDs if both rsID and SNPID are absent
    if ("rsID" not in sumstats.columns) and ("SNPID" not in sumstats.columns):
        if ("CHR" in sumstats.columns) and ("POS" in sumstats.columns):
            if ("EA" in sumstats.columns) and ("NEA" in sumstats.columns):
                sumstats = sumstats.with_columns(
                    (pl.col("CHR").cast(pl.String) + ":" +
                     pl.col("POS").cast(pl.String) + ":" +
                     pl.col("NEA").cast(pl.String) + ":" +
                     pl.col("EA").cast(pl.String)).alias("SNPID")
                )
            else:
                sumstats = sumstats.with_columns(
                    (pl.col("CHR").cast(pl.String) + ":" +
                     pl.col("POS").cast(pl.String)).alias("SNPID")
                )
            log.write(" -No rsID/SNPID found; created SNPID from CHR:POS[:NEA:EA]", verbose=verbose)
        else:
            sumstats = sumstats.with_columns(
                pl.lit(None).cast(pl.String).alias("SNPID")
            )
            log.warning(" -No rsID/SNPID and missing CHR/POS; created empty SNPID", verbose=verbose)

    check_datatype_polars(sumstats,log=log,verbose=verbose)
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
                    readargs["skip_rows"]=skip
                    readargs["separator"]="\t"
                    break
    readargs_header = readargs.copy()
    readargs_header["n_rows"] = 1
    readargs_header["infer_schema_length"] = 0
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
    # Split the study column by ":" and expand into columns
    format_cols_list = format_cols if isinstance(format_cols, list) else format_cols.split(":")
    # Split the study column and create new columns
    split_cols = []
    for i, col_name in enumerate(format_cols_list):
        split_cols.append(
            pl.col(study).str.split(":").list.get(i).alias(col_name)
        )
    sumstats = sumstats.with_columns(split_cols)
    sumstats = sumstats.drop(["FORMAT", study])
    sumstats = sumstats.select([pl.col(c) for c in vcf_usecols if c in sumstats.columns])
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
    pre_number=sumstats.height
    log.write(" -Checking if 0<= NEAF <=1 ...",verbose=verbose) 
    if "NEAF" in sumstats.columns:
        sumstats = sumstats.with_columns(
            pl.col("NEAF").cast(pl.Float64, strict=False)
        )
        sumstats = sumstats.filter((pl.col("NEAF")>=0) & (pl.col("NEAF")<=1))
        sumstats = sumstats.with_columns(
            (1 - pl.col("NEAF")).alias("EAF")
        )
        sumstats = sumstats.drop("NEAF")
    else:
        sumstats = sumstats.with_columns(
            pl.col("EAF").cast(pl.Float64, strict=False)
        )
        sumstats = sumstats.filter((pl.col("EAF")>=0) & (pl.col("EAF")<=1))
        sumstats = sumstats.with_columns(
            (1 - pl.col("EAF")).alias("EAF")
        )
    log.write(" -Converted NEAF to EAF.",verbose=verbose) 
    after_number=sumstats.height
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
                sumstats = sumstats.with_columns(pl.col("REF").alias("NEA"))
            
            log.write(" -EA,REF and ALT columns are available: assigning NEA...",verbose=verbose) 
            ea_alt = pl.col("EA") == pl.col("ALT")
            
            log.write(" -For variants with EA == ALT : assigning REF to NEA ...",verbose=verbose) 
            log.write(" -For variants with EA != ALT : assigning ALT to NEA ...",verbose=verbose) 
            # Single operation: if EA==ALT then REF, else ALT
            sumstats = sumstats.with_columns(
                pl.when(ea_alt)
                .then(pl.col("REF"))
                .otherwise(pl.col("ALT"))
                .alias("NEA")
            )

    if "NEA" not in sumstats.columns and "EA" not in sumstats.columns:
        if "REF" in sumstats.columns and "ALT" in sumstats.columns:
            log.write(" -Converting REF and ALT to NEA and EA ...")
            sumstats = sumstats.with_columns([
                pl.col("REF").alias("NEA"),
                pl.col("ALT").alias("EA")
            ]).drop(["ALT","REF"])
    return sumstats

def process_status(sumstats,build,status, log,verbose):
    if status is None:
        log.write(" -Initiating a status column: STATUS ...",verbose=verbose)
        #sumstats["STATUS"] = int(build)*(10**5) +99999
        build = _process_build(build,log,verbose)
        # Create integer status code: build (2 digits) + 99999 (5 digits) = 7 digits
        sumstats = sumstats.with_columns(
            pl.lit(int(build + "99999")).alias("STATUS")
        )

    # Convert STATUS to Int64 (integer type)
    if "STATUS" in sumstats.columns:
        sumstats = sumstats.with_columns(
            pl.col("STATUS").cast(pl.Int64)
        )
    return sumstats


def _load_single_chr(inpath,usecols,dtype_dictionary,readargs,rename_dictionary,chrom_pat,log,verbose):
    explicit = {"columns","schema_overrides","chunk_size"}
    readargs_clean = {k: v for k, v in readargs.items() if k not in explicit}
    # get chr 
    for k,v in rename_dictionary.items():
        if v=="CHR":
            if k in usecols:
                log.write(" -Columns used to filter variants: {}".format(k),verbose=verbose)
                chunk_chrom = k
                break

    log.write(" -Loading only variants on chromosome with pattern : {} ...".format(chrom_pat),verbose=verbose)
    
    # Read in chunks and filter
    sumstats_iter = pl.scan_csv(inpath,
                columns=usecols,
                schema_overrides=dtype_dictionary, 
                **readargs_clean)
    
    # Filter by chromosome pattern
    sumstats_filtered = sumstats_iter.filter(
        pl.col(chunk_chrom).str.contains(chrom_pat, literal=False)
    ).collect()
    
    log.write(" -Loaded {} variants on chromosome with pattern :{} ...".format(sumstats_filtered.height, chrom_pat),verbose=verbose)
    return sumstats_filtered

def _load_variants_with_pattern(inpath,usecols,dtype_dictionary,readargs,rename_dictionary,snpid_pat,log,verbose):
    explicit = {"columns","schema_overrides","chunk_size"}
    readargs_clean = {k: v for k, v in readargs.items() if k not in explicit}
    # get chr 
    for k,v in rename_dictionary.items():
        if v=="SNPID":
            if k in usecols:
                log.write(" -Columns used to filter variants: {}".format(k),verbose=verbose)
                chunk_snpid = k
                break

    log.write(" -Loading only variants with pattern :  {} ...".format(snpid_pat),verbose=verbose)
    
    # Read in chunks and filter
    sumstats_iter = pl.scan_csv(inpath,
                columns=usecols,
                schema_overrides=dtype_dictionary, 
                **readargs_clean)
    
    # Filter by snpid pattern
    sumstats_filtered = sumstats_iter.filter(
        pl.col(chunk_snpid).str.contains(snpid_pat, literal=False)
    ).collect()
    
    log.write(" -Loaded {} variants with pattern : {} ...".format(sumstats_filtered.height, snpid_pat),verbose=verbose)
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
    elif isinstance(sumstats, pl.DataFrame):
        inpath = None
        format_cols = None
        inpath_chr_list = None
        inpath_chr_num_list = None
        ## loading data from polars dataframe
        raw_cols = sumstats.columns

    ################################################
    for key,value in rename_dictionary.items():
        # check available keys  key->raw header
        # usecols : a list of raw headers to load from file/DataFrame 
        if key in raw_cols:
            usecols.append(key)
        if value in ["EA","NEA"]:
            dtype_dictionary[key]=pl.String()
        if value in ["STATUS"]:
            dtype_dictionary[key]=pl.String()     
        if value in ["CHR"]:
            dtype_dictionary[key]=pl.String()  

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
        row_one = pl.read_csv(inpath_chr_list[0],**readargs_header)
        # columns in the sumstats
        raw_cols = row_one.columns
    else:
    ##### loading data from tabular file#################################################
    #if file not found, FileNotFoundError
        readargs_header = get_readargs_header(inpath = inpath, readargs = readargs)
        row_one = pl.read_csv(inpath,**readargs_header)
        raw_cols = row_one.columns

    if fmt=="vcf":
        # expanded
        format_cols = list(row_one["FORMAT"].str.split(":"))[0]
        # fixed + study1 + expanded
        raw_cols = meta_data["format_fixed"] + [raw_cols[9]] + format_cols

    return format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list


def _sort_column(sumstats,verbose=True,log=Log(),order = None):
    '''
    Reorder columns according to a specified order.
    
    Reorders the dataframe columns to match a predefined standard order, placing standard
    GWAS columns first (SNPID, rsID, CHR, POS, EA, NEA, statistics, etc.) followed by
    any additional columns not in the standard list.

    Parameters
    ----------
    verbose : bool, optional
        Whether to print progress. Default is True.

    Returns
    -------
    pl.DataFrame
        Modified sumstats with reordered columns.
    '''

    if order is None:
        order = [
        "SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "MAF", "BETA", "SE","BETA_95L","BETA_95U", "Z","T","F",
        "CHISQ", "P", "MLOG10P", "OR", "OR_95L", "OR_95U","HR", "HR_95L", "HR_95U","INFO", "N","N_CASE","N_CONTROL","DIRECTION","I2","P_HET","DOF","SNPR2","STATUS"]
    output_columns = []
    for i in order:
        if i in sumstats.columns: output_columns.append(i)
    for i in sumstats.columns:
        if i not in order: output_columns.append(i)
    log.write(" -Reordering columns to    :", ",".join(output_columns), verbose=verbose)
    sumstats = sumstats.select(output_columns)

    return sumstats
