import pandas as pd
import yaml
import hashlib
import copy
from pysam import tabix_compress 
from pysam import tabix_index
from datetime import datetime
from datetime import date
from gwaslab.io_preformat_input import print_format_info
from gwaslab.bd_common_data import get_formats_list
from gwaslab.g_Log import Log
from gwaslab.bd_common_data import get_format_dict
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.g_version import gwaslab_info
from gwaslab.bd_get_hapmap3 import gethapmap3
from gwaslab.util_in_filter_value import _exclude_hla
from gwaslab.util_in_filter_value import _exclude
from gwaslab.util_in_filter_value import _extract
# to vcf
# to fmt
    ## vcf
    ## vep
    ## bed
    ## annovar
    ## general : ldsc, plink, plink2, saige, regenie 
###################################################################################################################################################
def _to_format(sumstats,
              path="./sumstats",
              fmt="gwaslab",   
              extract=None,
              exclude=None,
              cols=None,
              id_use="rsID",
              hapmap3=False,
              exclude_hla=False,
              hla_range=(25,34),  
              build=None, 
              n=None,
              no_status=False,
              output_log=True,
              to_csvargs=None,
              float_formats=None,
              xymt_number=False,
              xymt=None,
              chr_prefix="",
              meta=None,
              ssfmeta=False,
              md5sum=False,
              bgzip=False,
              tabix=False,
              tabix_indexargs={},
              log=Log(),
              verbose=True):
    
    if  to_csvargs is None:
        to_csvargs = {}
    if  float_formats is None:
        float_formats={}
    if cols is None:
        cols=[]
    if xymt is None:
        xymt = ["X","Y","MT"]
    onetime_log = copy.deepcopy(log)

    #######################################################################################################
    
    formatlist= get_formats_list() + ["vep","bed","annovar","vcf"]
    if fmt in formatlist:
        onetime_log.write("Start to convert the output sumstats in: ",fmt, " format",verbose=verbose)
    else:
        raise ValueError("Please select a format to output")
    suffix=fmt
    
    #######################################################################################################
    # filter
    output = sumstats.copy()
    
    #######################################################################################################

    if extract is not None:
        output = _extract(output, extract, id_use=id_use, log=onetime_log, verbose=verbose)
    
    if exclude is not None:
        output = _exclude(output, exclude, id_use=id_use, log=onetime_log, verbose=verbose)
    
    #hla and hapmap3 #######################################################################################
    
    #exclude hla
    if exclude_hla==True:
        output = _exclude_hla(output, lower=hla_range[0]*1000000 ,upper=hla_range[1]*1000000 ,log=onetime_log, verbose=verbose)
        suffix = "noMHC."+suffix
    
    #extract hapmap3 SNPs
    if hapmap3==True:
        output = gethapmap3(output,build=build,verbose=verbose)
        after = len(output)
        onetime_log.write(" -Extract {} variants in Hapmap3 datasets for build {}.".format(after, build ),verbose=verbose)
        suffix = "hapmap3."+suffix
    
    # add a n column
    if n is not None:
        output["N"] = n
            
    #######################################################################################################
    #formatting float statistics
    onetime_log.write(" -Formatting statistics ...",verbose=verbose)
    
    formats = {
            'EAF': '{:.4g}', 
            'MAF': '{:.4g}', 
            'BETA': '{:.4f}', 
            'SE': '{:.4f}',
            'BETA_95U': '{:.4f}',
            'BETA_95L': '{:.4f}',
            'Z': '{:.4f}',
            'CHISQ': '{:.4f}',
            'F': '{:.4f}',
            'OR': '{:.4f}',
            'OR_95U': '{:.4f}',
            'OR_95L': '{:.4f}',
            'HR': '{:.4f}',
            'HR_95U': '{:.4f}',
            'HR_95L': '{:.4f}',
            'INFO': '{:.4f}',
            'P': '{:.4e}',
            'MLOG10P': '{:.4f}',
            'DAF': '{:.4f}'}
    
    for col, f in float_formats.items():
        if col in output.columns: 
            formats[col]=f
    
    for col, f in formats.items():
        if col in output.columns: 
            if str(output[col].dtype) in ["Float32","Float64","float64","float32","float16","float"]:
                output[col] = output[col].map(f.format)
    
    onetime_log.write(" -Float statistics formats:",verbose=verbose)  
    keys=[]
    values=[]
    for key,value in formats.items():
        if key in output.columns: 
            keys.append(key)
            values.append(value)
    
    onetime_log.write("  - Columns       :",keys,verbose=verbose) 
    onetime_log.write("  - Output formats:",values,verbose=verbose) 
        
    ##########################################################################################################          
    # output, mapping column names
    
    if fmt in get_formats_list() + ["vep","bed","annovar","vcf"]:
        tofmt(output,
              path=path,
              fmt=fmt,
              cols=cols,
              suffix=suffix,
              build=build,
              verbose=verbose,
                no_status=no_status,
                log=onetime_log,
                to_csvargs=to_csvargs,
                chr_prefix=chr_prefix,
                meta=meta,
                ssfmeta=ssfmeta,
                bgzip=bgzip,
                tabix=tabix,
                tabix_indexargs=tabix_indexargs,
                md5sum=md5sum,
                xymt_number=xymt_number,
                xymt=xymt)
    
    if output_log is True:
        log_path = path + "."+ suffix + ".log"
        onetime_log.write(" -Saving log file to: {}".format(log_path),verbose=verbose)
        onetime_log.write("Finished outputting successfully!",verbose=verbose)
        try:
            onetime_log.save(log_path, verbose=False)
        except:
            pass

###################################################################################################################################################
def tofmt(sumstats,
          meta,
          path=None,
          suffix=None,
          fmt=None,
          cols=[],
          xymt_number=False,
          xymt=["X","Y","MT"],
          build="19",
          chr_prefix="",
          ssfmeta=False,
          md5sum=False,
          bgzip=False,
          tabix=False,
          tabix_indexargs={},
          verbose=True,
          no_status=False,
          log=Log(),
          to_csvargs=None):
    
    if to_csvargs is None:
        to_csvargs=dict()
    
    if fmt in ["ssf"]: 
        xymt_number=True
        if "SNPID" in sumstats.columns:
            log.write(' -Replacing SNPID separator from ":" to "_"...')
            sumstats["SNPID"] = sumstats["SNPID"].str.replace(":","_")
    log.write(" -Start outputting sumstats in "+fmt+" format...")
    
    if "CHR" in sumstats.columns:
        if xymt_number is False and pd.api.types.is_integer_dtype(sumstats["CHR"]):
            sumstats["CHR"]= sumstats["CHR"].map(get_number_to_chr(xymt=xymt,prefix=chr_prefix))
        elif chr_prefix is not None:
            sumstats["CHR"]= chr_prefix + sumstats["CHR"].astype("string")

    ####################################################################################################################
    if fmt=="bed":
        # bed-like format, 0-based, 
        # first 3 columns : chromosome, start, end
        # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
        is_snp,is_insert,is_delete = _check_indel(sumstats,log,verbose)
        log.write(" -formatting to 0-based bed-like file...")
        log.write(" -format description: {}".format("https://genome.ucsc.edu/FAQ/FAQformat.html#format1"))
        
        sumstats = _adjust_position(sumstats, fmt, is_snp, is_insert, is_delete, log, verbose )

        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"] + cols

        _output_bed_like(sumstats,  path, "bed", suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose)
    ####################################################################################################################   
    elif fmt=="vep":
        # bed-like format, 1-based
        # first 6 columns : chromosome, start, end, allele, strand, identifier
        # https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html
        
        is_snp,is_insert,is_delete = _check_indel(sumstats,log,verbose)
            
        log.write(" -formatting to 1-based bed-like file (for vep)...")
        log.write(" -format description: {}".format("http://asia.ensembl.org/info/docs/tools/vep/vep_formats.html"))
        sumstats = _adjust_position(sumstats, fmt, is_snp, is_insert, is_delete , log, verbose)
        
        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"]+ cols

        _output_bed_like(sumstats,  path,"vep", suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose)

    ####################################################################################################################
    elif fmt=="annovar":
        # bed-like format, 1-based, 
        # first 3 columns : Chromosome ("chr" prefix is optional), Start, End, Reference Allelel, Alternative Allele
        # https://annovar.openbioinformatics.org/en/latest/user-guide/input/
        is_snp,is_insert,is_delete = _check_indel(sumstats,log,verbose)

        log.write(" -formatting to 1-based bed-like file...")
        log.write(" -format description: {}".format("https://annovar.openbioinformatics.org/en/latest/user-guide/input/"))
        
        sumstats = _adjust_position(sumstats, fmt, is_snp, is_insert, is_delete, log, verbose )

        ouput_cols=["CHR","START","END","NEA_out","EA_out","SNPID"]+ cols
        
        _output_bed_like(sumstats, path, fmt, suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose)
    
    ####################################################################################################################       
    elif fmt=="vcf":
        # GWAS-VCF
        log.write(" -"+fmt+" format will be loaded...",verbose=verbose)
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True, skip_meta_records=["format_fixed_header","format_contig_19","format_contig_38"])
        
        # determine which ID to use
        if "rsID" in sumstats.columns:
            rename_dictionary["rsID"]="ID"
        else:
            rename_dictionary["SNPID"]="ID" 
        
        # get the columns to output
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        ouput_cols = ouput_cols +["STATUS"]+ cols
        sumstats = sumstats[ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        
        # replace : with _
        sumstats["ID"] = sumstats["ID"].str.replace(":","_")
        
        # process Allele frequency data
        if "AF" in sumstats.columns:
            sumstats["INFO"] = "AF="+sumstats["AF"].astype("string")
        else:
            sumstats["INFO"] = "."

        # sumstats columns placed in vcf-FORMAT column     
        output_format=[]
        for i in sumstats.columns:
            if i in meta_data["format_format"]:
                output_format.append(i)  

        # determine path
        path = path + "."+suffix
        
        
        vcf_header =  _process_vcf_header(sumstats, meta, meta_data, build, log, verbose)

        log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)
        # output header
        with open(path,"w") as file:
            file.write(vcf_header)
        
        with open(path,"a") as file:
            log.write(" -Output columns:"," ".join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]]))
            file.write("\t".join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]])+"\n")
            log.write(" -Outputing data...")
            QUAL="."
            FILTER="PASS"
            for index,row in sumstats.iterrows():
                CHROM=str(row["#CHROM"])
                POS=str(row["POS"])
                ID=str(row["ID"])
                REF=str(row["REF"])
                ALT=str(row["ALT"])
                INFO=str(row["INFO"])
                FORMAT=":".join(output_format)
                DATA=":".join(row[output_format].astype("string"))
                file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, DATA))
        _bgzip_tabix_md5sum(path, fmt, bgzip, md5sum, tabix, tabix_indexargs, log, verbose)
    
    ####################################################################################################################
    elif fmt in get_formats_list():      
        # tabular 
        log.write(" -"+fmt+" format will be loaded...",verbose=verbose)
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True)
        
        ymal_path = path + "."+suffix+".tsv-meta.ymal"
        path = path + "."+suffix+".tsv.gz"
        log.write(" -Output path:",path, verbose=verbose) 
        
        sumstats,to_csvargs = _configure_output_cols_and_args(sumstats, rename_dictionary, cols, no_status, path, meta_data, to_csvargs, log, verbose)
        
        log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)
        sumstats.to_csv(path, index=None,**to_csvargs)

        if md5sum == True: 
            md5_value = md5sum_file(path,log,verbose)
        else:
            md5_value = calculate_md5sum_file(path)
        
        ## update ssf-style meta data and export to yaml file
        _configure_ssf_meta(sumstats, fmt, ssfmeta, meta, meta_data, path, md5_value, ymal_path, log, verbose)
        
        return sumstats  
####################################################################################################################
def _configure_output_cols_and_args(sumstats, rename_dictionary, cols, no_status, path, meta_data, to_csvargs, log, verbose):
    # grab format cols that exist in sumstats
    ouput_cols=[]
    for i in sumstats.columns:
        if i in rename_dictionary.keys():
            ouput_cols.append(i)  
    
    # + additional cols and remove duplicated
    ouput_cols = list(set(ouput_cols + cols)) 
    
    # remove STATUS
    try:
        if no_status == True:
            ouput_cols.remove("STATUS")
    except:
        pass
    
    #filter and rename to target fromat headers
    sumstats = sumstats[ouput_cols]
    sumstats = sumstats.rename(columns=rename_dictionary) 

    # configure target format args and reorder columns
    if "format_separator" in meta_data.keys():
        to_csvargs["sep"] = meta_data["format_separator"]
    else:
        to_csvargs["sep"]="\t"
    if "format_na" in meta_data.keys():
        to_csvargs["na_rep"] = meta_data["format_na"]
    if "format_col_order" in meta_data.keys():
        fixed_col =[]
        other_col=[]
        for i in meta_data["format_col_order"]:
            if i in sumstats.columns:
                fixed_col.append(i)
        for i in sumstats.columns:
            if i not in meta_data["format_col_order"]:
                other_col.append(i)
        sumstats = sumstats[fixed_col + other_col]
    log.write(" -Output columns: {}".format(",".join(sumstats.columns)),verbose=verbose)
    return sumstats, to_csvargs


def _configure_ssf_meta(sumstats, fmt, ssfmeta, meta, meta_data, path, md5_value, ymal_path, log, verbose):
    ### calculate meta data
    if "EAF" in sumstats.columns:
        min_maf = sumstats["EAF"].min()
    else:
        min_maf = "Unknown"
    
    if "N" in sumstats.columns:
        n_median =  sumstats["N"].median()
        n_max = sumstats["N"].max()
        n_min = sumstats["N"].min()
    else:
        n_median = "Unknown"
        n_max = "Unknown"
        n_min = "Unknown"

    if ssfmeta==True:
        sumstats_meta_copy = meta.copy()
        if "format_cite_name" in meta_data.keys():
            sumstats_meta_copy["file_type"] = meta_data["format_cite_name"]
        else:
            sumstats_meta_copy["file_type"] = fmt
        sumstats_meta_copy["minor_allele_freq_lower_limit"] = min_maf
        sumstats_meta_copy["data_file_name"] = path
        sumstats_meta_copy["data_file_md5sum"] = md5_value
        sumstats_meta_copy["date_last_modified"] = get_format_date_and_time()
        sumstats_meta_copy["samples"]["sample_size"] = n_max
        sumstats_meta_copy["gwaslab"]["samples"]["sample_size_min"] = n_min
        sumstats_meta_copy["gwaslab"]["samples"]["sample_size_median"] = n_median
        sumstats_meta_copy["gwaslab"]["variants"]["variant_number"] = len(sumstats)
        log.write(" -Exporting SSF-style meta data to {}".format(ymal_path),verbose=verbose) 
        with open(ymal_path, 'w') as outfile:
            yaml.dump(sumstats_meta_copy, outfile)



def _output_bed_like(sumstats, path, fmt, suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose):
    sumstats = sumstats[ouput_cols]
    path = path + "."+suffix
    log.write(" -Output columns: {}".format(",".join(sumstats.columns)),verbose=verbose)
    log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)
    sumstats.to_csv(path,sep="\t",index=None,header=None,**to_csvargs)
    _bgzip_tabix_md5sum(path, fmt, bgzip, md5sum, tabix, tabix_indexargs, log, verbose)


def _bgzip_tabix_md5sum(path, fmt, bgzip, md5sum, tabix, tabix_indexargs, log, verbose):
    if bgzip == True:
        log.write(" -bgzip compressing : {}...".format(path+".gz"),verbose=verbose) 
        tabix_compress(path, path+".gz",force=True)
    if md5sum == True: 
        if bgzip == True:
            md5sum_file(path+".gz",log,verbose)
        else:
            md5sum_file(path,log,verbose)
    if tabix == True and bgzip == True:
        log.write(" -tabix indexing : : {}...".format(path+".gz.tbi"),verbose=verbose)
        if "preset" not in  tabix_indexargs:
            tabix_indexargs["preset"] = fmt
        if "force" not in tabix_indexargs:
            tabix_indexargs["force"] = True
        tabix_index(path+".gz", **tabix_indexargs)    


def _check_indel(sumstats,log,verbose):
    is_snp = (sumstats["EA"].str.len() == sumstats["NEA"].str.len())
    is_insert = (sumstats["EA"].str.len()>1) &(sumstats["NEA"].str.len()==1)
    is_delete = (sumstats["EA"].str.len()==1) &(sumstats["NEA"].str.len()>1)
    
    log.write(" -Number of SNPs :",sum(is_snp))
    log.write(" -Number of Insertions :",sum(is_insert)) 
    log.write(" -Number of Deletions :",sum(is_delete)) 
    return is_snp,is_insert,is_delete


def md5sum_file(filename,log,verbose):
    log.write(" -md5sum hashing for the file:",filename,verbose=verbose) 
    md5_hash = hashlib.md5()
    with open(filename,"rb") as f:
        # Read and update hash in chunks
        for byte_block in iter(lambda: f.read(4096*1000),b""):
            md5_hash.update(byte_block)
    with open(filename+".md5sum","w") as f:
        out = str(md5_hash.hexdigest())
        f.write(out+"\n")
        log.write(" -md5sum path:",filename+".md5sum",verbose=verbose)    
        log.write(" -md5sum: {}".format(out),verbose=verbose) 
        return out

def calculate_md5sum_file(filename):
    md5_hash = hashlib.md5()
    with open(filename,"rb") as f:
        # Read and update hash in chunks
        for byte_block in iter(lambda: f.read(4096*1000),b""):
            md5_hash.update(byte_block)
        out = str(md5_hash.hexdigest())
        return out    

def get_format_date_and_time():
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d-%H:%M:%S")
    return dt_string


def _adjust_position(sumstats, fmt,is_snp, is_insert, is_delete, log, verbose):
    log.write(" -Adjusting positions in format-specific manner..",verbose=verbose) 
    if fmt=="bed":
        sumstats.loc[is_snp,"START"]  = sumstats.loc[is_snp,"POS"]-1 
        sumstats.loc[is_snp,"END"]    = sumstats.loc[is_snp,"POS"]-1 + sumstats.loc[is_snp,"NEA"].str.len()
        sumstats.loc[is_snp,"NEA/EA"] = sumstats.loc[is_snp,"NEA"].astype("string")+"/"+sumstats.loc[is_snp,"EA"].astype("string")
        
        # for insertion
        # start = pos : end = pos
        # A/ATC -> -/TC
        sumstats.loc[is_insert,"START"]  = sumstats.loc[is_insert,"POS"]
        sumstats.loc[is_insert,"END"]    = sumstats.loc[is_insert,"POS"]
        sumstats.loc[is_insert,"NEA/EA"] = "-/"+sumstats.loc[is_insert,"EA"].str.slice(start=1)
        
        # for deletion 
        # start = pos - 1 +1; end = pos -1 +1+ len(Ref)
        # ATC/A -> TC/-
        sumstats.loc[is_delete,"START"]  = sumstats.loc[is_delete,"POS"]
        sumstats.loc[is_delete,"END"]    = sumstats.loc[is_delete,"POS"] + sumstats.loc[is_delete,"NEA"].str.len() - 1
        sumstats.loc[is_delete,"NEA/EA"] = sumstats.loc[is_delete,"NEA"].str.slice(start=1)+"/-"
        sumstats["STRAND"]="+"
    elif fmt=="vep":
        sumstats.loc[is_snp,"START"]  = sumstats.loc[is_snp,"POS"] + (sumstats.loc[is_snp,"NEA"].str.len() - 1 )
        sumstats.loc[is_snp,"END"]    = sumstats.loc[is_snp,"POS"] + (sumstats.loc[is_snp,"NEA"].str.len() - 1 )
        sumstats.loc[is_snp,"NEA/EA"] = sumstats.loc[is_snp,"NEA"].astype("string")+"/"+sumstats.loc[is_snp,"EA"].astype("string")
        
        # for insertion
        # start = pos+1 ; end = pos
        # A/ATC -> -/TC
        sumstats.loc[is_insert,"START"]  = sumstats.loc[is_insert,"POS"] + 1
        sumstats.loc[is_insert,"END"]    = sumstats.loc[is_insert,"POS"]
        sumstats.loc[is_insert,"NEA/EA"] = "-/" + sumstats.loc[is_insert,"EA"].str.slice(start=1)
        
        # for deletion 
        # start = pos ; end = pos + len(Ref) -1
        # ATC/A -> TC/-
        sumstats.loc[is_delete,"START"]  = sumstats.loc[is_delete,"POS"] + 1
        sumstats.loc[is_delete,"END"]    = sumstats.loc[is_delete,"POS"] + (sumstats.loc[is_delete,"NEA"].str.len() -1)
        sumstats.loc[is_delete,"NEA/EA"] = sumstats.loc[is_delete,"NEA"].str.slice(start=1)+"/-"
        sumstats["STRAND"]="+"
    elif fmt=="annovar":
        # for snp  
        # start = pos ; end = pos
        # A/G
        # AT/CG
        sumstats.loc[is_snp,"START"]  = sumstats.loc[is_snp,"POS"]
        sumstats.loc[is_snp,"END"]    = sumstats.loc[is_snp,"POS"]-1 + sumstats.loc[is_snp,"NEA"].str.len()
        sumstats.loc[is_snp,"NEA_out"] = sumstats.loc[is_snp,"NEA"].astype("string")
        sumstats.loc[is_snp,"EA_out"] = sumstats.loc[is_snp,"EA"].astype("string")
        # for insertion
        # start = pos : end = pos
        # A/ATC -> -/TC
        sumstats.loc[is_insert,"START"]  = sumstats.loc[is_insert,"POS"]+1
        sumstats.loc[is_insert,"END"]   = sumstats.loc[is_insert,"POS"]+1
        sumstats.loc[is_insert,"NEA_out"] = "-"
        sumstats.loc[is_insert,"EA_out"] = sumstats.loc[is_insert,"EA"].str.slice(start=1)
        
        # for deletion 
        # start = pos - 1 +1; end = pos -1 +1+ len(Ref)
        # ATC/A -> TC/-
        sumstats.loc[is_delete,"START"] = sumstats.loc[is_delete,"POS"]
        sumstats.loc[is_delete,"END"]  = sumstats.loc[is_delete,"POS"]- 1 + sumstats.loc[is_delete,"NEA"].str.len() 
        sumstats.loc[is_delete,"NEA_out"] = sumstats.loc[is_delete,"NEA"].str.slice(start=1)
        sumstats.loc[is_delete,"EA_out"] = "-"

    
    sumstats["START"] = sumstats["START"].astype("Int64")
    sumstats["END"] = sumstats["END"].astype("Int64")
    return sumstats

def _process_vcf_header(sumstats, meta, meta_data, build, log, verbose):
    
    log.write(" -Creating VCF file header...",verbose=verbose)
    log.write("  -VCF header contig build:"+str(build),verbose=verbose) 
    
    # calculate meta data
    harmonised = sum(sumstats["STATUS"].str.match( r"\w\w\w[0][0123][012][01234]", case=False, flags=0, na=False ) )
    switchedalleles = sum(sumstats["STATUS"].str.match( r"\w\w\w[0][0123][12][24]", case=False, flags=0, na=False ) )

    # Create vcf header
    vcf_header = meta_data["format_fixed_header"] +"\n"+ meta_data["format_contig_"+str(build)]+"\n"

    # Create sample header
    vcf_header+="##SAMPLE=<ID={},TotalVariants={},VariantsNotRead=0,HarmonisedVariants={},VariantsNotHarmonised={},SwitchedAlleles={},StudyType={}>\n".format(
                    meta["gwaslab"]["study_name"], len(sumstats), harmonised, len(sumstats)-harmonised, switchedalleles, meta["gwaslab"]["study_type"])
    vcf_header+="##gwaslab_version="+gwaslab_info()["version"]+"\n"

    log.write("  -ID:{}".format( meta["gwaslab"]["study_name"]),verbose=verbose) 
    log.write("  -StudyType:{}".format(meta["gwaslab"]["study_type"]),verbose=verbose) 
    log.write("  -TotalVariants:{}".format(len(sumstats)),verbose=verbose) 
    log.write("  -HarmonisedVariants:{}".format(harmonised),verbose=verbose) 
    log.write("  -VariantsNotHarmonised:{}".format(len(sumstats)-harmonised),verbose=verbose) 
    log.write("  -SwitchedAlleles:{}".format(switchedalleles),verbose=verbose) 

    return vcf_header
        