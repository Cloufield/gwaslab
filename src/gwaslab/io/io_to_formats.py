import pandas as pd
import yaml
import hashlib
import copy
import os
"""
Convert summary statistics to specified output format with various processing options.
"""

from pysam import tabix_compress 
from pysam import tabix_index
from datetime import datetime
from datetime import date
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import gwaslab_info

from gwaslab.io.io_preformat_input import _print_format_info

from gwaslab.bd.bd_common_data import get_format_dict
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_formats_list
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap3
from gwaslab.bd.bd_sex_chromosomes import Chromosomes

from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.util.util_in_filter_value import _exclude
from gwaslab.util.util_in_filter_value import _extract
# to vcf
# to fmt
    ## vcf
    ## vep
    ## bed
    ## annovar
    ## general : ldsc, plink, plink2, saige, regenie 
###################################################################################################################################################
def _to_format(sumstats_or_dataframe,
              path="./sumstats",
              fmt="gwaslab",   
              tab_fmt="tsv",
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
              float_formats=None,
              xymt_number=False,
              xymt=None,
              chr_prefix="",
              meta=None,
              ssfmeta=False,
              md5sum=False,
              gzip=True,
              bgzip=False,
              tabix=False,
              tabix_indexargs={},
              to_csvargs=None,
              to_tabular_kwargs=None,
              validate=False,
              gwas_ssf_path=None,
              log=Log(),
              verbose=True):
    """
    Output the sumstats. Convert summary statistics to specified output format with various processing options. VCF, VEP, BED, tsv-like formats.
    
    Parameters
    ----------
    path : str, optional
        Output file path prefix. Default is "./sumstats"
    fmt : str, optional, default="gwaslab"
        Output format. Supported formats include 'gwaslab', 'vcf', 'bed', 'annovar', etc.  The list of available formats can be checked using `get_formats_list()`.  
        Default is 'gwaslab'.
    tab_fmt : str, optional
        Tabular output format type used when `fmt` is not one of 'vcf', 'bed', or 'annovar'.  
        Supported options are 'tsv', 'csv', and 'parquet'.  
        Default is 'tsv'.
    extract : list or str, optional
        SNPs to extract. Default is None
    exclude : list or str, optional
        SNPs to exclude. Default is None
    cols : list, optional
        Additional columns to include in output. Default is empty list
    hapmap3 : bool, optional
        Whether to extract Hapmap3 SNPs. Default is False
    exclude_hla : bool, optional
        Whether to exclude HLA region. Default is False
    hla_range : tuple, optional
        HLA exclusion range in Mbp. Default is (25,34)
    build : str, optional
        Genome build version. Default is None
    n : float, optional
        Sample size to add as 'N' column. Default is None
    no_status : bool, optional
        Whether to exclude 'STATUS' column. Default is False
    output_log : bool, optional
        Whether to save a log file. Default is True
    float_formats : dict, optional
        Dictionary of float formatting strings. Default is None
    xymt_number : bool, optional
        Whether to use numeric codes for X/Y/MT chromosomes. Default is False
    xymt : list, optional
        Chromosome names for X, Y, MT. Default is ["X","Y","MT"]
    chr_prefix : str, optional
        Prefix for chromosome numbers. Default is empty string
    ssfmeta : bool, optional
        Whether to create SSF-style metadata. Default is False
    md5sum : bool, optional
        Whether to generate MD5 checksum. Default is False
    gzip : bool, optional
        Whether to gzip compress output. Default is True
    bgzip : bool, optional
        Whether to bgzip compress output. Default is False
    tabix : bool, optional
        Whether to create Tabix index. Default is False
    tabix_indexargs : dict, optional
        Dictionary of Tabix indexing arguments. Default is empty dict
    to_csvargs : dict, optional
        Dictionary of CSV writing arguments. Default is None
    to_tabular_kwargs : dict, optional
        Dictionary of tabular format arguments. Default is None
    verbose : bool, optional
        Whether to print progress messages. Default is True
    validate : bool, optional
         Whether to use the gwas-ssf CLI tool for validation (only for SSF format). Default is False
    gwas_ssf_path : str, optional
         Path to gwas-ssf CLI executable. If None, will try default paths. Default is None

    Returns
    -------
    None
        Output is written to file(s) specified by path parameter

    Less used parameters
    ---------------------------------------------------------------
    log : gwaslab.Log, optional
        Log object for tracking progress. Default is new Log instance
    meta : dict, optional
        Metadata dictionary. Default is None
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
        chromosomes_obj = None
    else:
        sumstats = sumstats_or_dataframe.data
        # Get Chromosomes instance from Sumstats object if available
        chromosomes_obj = getattr(sumstats_or_dataframe, 'chromosomes', None)
             
    if to_csvargs is None:
        to_csvargs=dict()
    if tabix_indexargs is None:
        tabix_indexargs=dict()
    if to_tabular_kwargs is None:
        to_tabular_kwargs=dict()
    if  float_formats is None:
        float_formats=dict()
    if cols is None:
        cols=[]
    if xymt is None:
        # Use chromosomes from Sumstats object if available, otherwise default to human
        if chromosomes_obj is not None:
            # Build xymt list from chromosomes object
            xymt_list = []
            if len(chromosomes_obj.sex_chromosomes) >= 1:
                xymt_list.append(chromosomes_obj.sex_chromosomes[0])
            if len(chromosomes_obj.sex_chromosomes) >= 2:
                xymt_list.append(chromosomes_obj.sex_chromosomes[1])
            if chromosomes_obj.mitochondrial:
                xymt_list.append(chromosomes_obj.mitochondrial)
            xymt = xymt_list if xymt_list else ["X", "Y", "MT"]  # Fallback if no sex chromosomes
        else:
            xymt = ["X","Y","MT"]  # Default for DataFrame case
    non_gzip_tab_fmt = ["parquet"]
    non_md5sum_tab_fmt = ["parquet"]

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
    if id_use =="rsID":
        if id_use not in sumstats.columns:
            id_use="SNPID"

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
        output = _get_hapmap3(output,build=build,verbose=verbose)
        after = len(output)
        onetime_log.write(" -Extract {} variants in Hapmap3 datasets for build {}.".format(after, build ),verbose=verbose)
        suffix = "hapmap3."+suffix
    
    # add a n column
    if n is not None:
        output["N"] = n
            
    #######################################################################################################
    #formatting float statistics
    
    if tab_fmt!="parquet":
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
        output_file = tofmt(output,
              path=path,
              fmt=fmt,
              tab_fmt=tab_fmt,
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
                gzip=gzip,
                bgzip=bgzip,
                non_gzip_tab_fmt=non_gzip_tab_fmt,
                non_md5sum_tab_fmt=non_md5sum_tab_fmt,
                tabix=tabix,
                tabix_indexargs=tabix_indexargs,
                to_tabular_kwargs=to_tabular_kwargs,
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

    # Validate SSF format output if requested
    if fmt == "ssf":
        try:
            import os
            import subprocess
            from pathlib import Path
            
            output_path = Path(output_file)
            if not output_path.exists():
                log.write(f"Output file not found for validation: {output_file}", verbose=verbose)
                return
            
            log.write("Validating SSF format output...", verbose=verbose)
            
            # Try to use gwas-ssf CLI first
            if gwas_ssf_path:
                # Use user-specified path
                gwas_ssf_paths = [gwas_ssf_path]
            else:
                # Default: try system PATH
                gwas_ssf_paths = ["gwas-ssf"]
            
            cli_available = False
            for cli_path in gwas_ssf_paths:
                if os.path.exists(cli_path) or cli_path == "gwas-ssf":
                    try:
                        result = subprocess.run(
                            [cli_path, "validate", str(output_path)],
                            capture_output=True,
                            text=True,
                            timeout=120
                        )
                        cli_available = True
                        
                        if result.returncode == 0:
                            log.write("✓ SSF validation successful (via gwas-ssf CLI)", verbose=verbose)
                            if result.stdout:
                                log.write(result.stdout.strip(), verbose=verbose)
                        else:
                            log.warning("✗ SSF validation failed (via gwas-ssf CLI)", verbose=verbose)
                            if result.stderr:
                                log.warning(result.stderr.strip()[:500], verbose=verbose)
                            if result.stdout:
                                log.write(result.stdout.strip()[:500], verbose=verbose)
                        break
                    except (FileNotFoundError, subprocess.TimeoutExpired, Exception) as e:
                        continue
            
            # Fallback to new validator if CLI not available
            if not cli_available:
                try:
                    from gwaslab.extension.gwas_sumstats_tools.validate_ssf import validate_ssf_file
                    
                    log.write("Using built-in SSF validator (gwas-ssf CLI not available)...", verbose=verbose)
                    is_valid, message, errors = validate_ssf_file(
                        filename=output_path,
                        log=log,
                        verbose=verbose,
                        minimum_rows=100_000,
                        pval_zero=False
                    )
                    
                    if is_valid:
                        log.write(f"✓ SSF validation successful: {message}", verbose=verbose)
                    else:
                        log.warning(f"✗ SSF validation failed: {message}", verbose=verbose)
                        if errors:
                            for error in errors[:10]:  # Show first 10 errors
                                log.warning(f"  - {error}", verbose=verbose)
                            if len(errors) > 10:
                                log.warning(f"  ... and {len(errors) - 10} more errors", verbose=verbose)
                except ImportError as e:
                    log.write(f"Warning: SSF validation not available (CLI and built-in validator both unavailable): {e}", verbose=verbose)
                except Exception as e:
                    log.warning(f"Error during SSF validation: {e}", verbose=verbose)
                    
        except Exception as e:
            log.warning(f"Error during SSF validation: {e}", verbose=verbose)


###################################################################################################################################################
def tofmt(sumstats,
          meta,
          path=None,
          suffix=None,
          fmt=None,
          tab_fmt="csv",
          cols=[],
          xymt_number=False,
          xymt=["X","Y","MT"],
          build="19",
          chr_prefix="",
          ssfmeta=False,
          md5sum=False,
          bgzip=False,
          gzip=True,
          non_gzip_tab_fmt=None,
          non_md5sum_tab_fmt=None,
          tabix=False,
          tabix_indexargs=None,
          verbose=True,
          no_status=False,
          log=Log(),
          to_csvargs=None,
          to_tabular_kwargs=None):
    
    if fmt in ["ssf"]: 
        xymt_number=True
        if "SNPID" in sumstats.columns:
            log.write(' -Replacing SNPID separator from ":" to "_"...')
            sumstats["SNPID"] = sumstats["SNPID"].str.replace(":","_")
    log.write(" -Start outputting sumstats in "+fmt+" format...")
    
    if "CHR" in sumstats.columns:
        # output X,Y,MT instead of 23,24,25
        if xymt_number is False and pd.api.types.is_integer_dtype(sumstats["CHR"]):
            # Get species from metadata if available, otherwise use default
            species = None
            if meta is not None and isinstance(meta, dict):
                species = meta.get("gwaslab", {}).get("species", None)
            # Use species-aware get_number_to_chr (xymt is already set correctly from _to_format)
            sumstats["CHR"]= sumstats["CHR"].map(get_number_to_chr(xymt=xymt, prefix=chr_prefix, species=species))
        # add prefix to CHR
        elif len(chr_prefix)>0:
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
        _print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True, skip_meta_records=["format_fixed_header","format_contig_19","format_contig_38"])
        
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
        try:
            fast_to_vcf(sumstats, path, vcf_header, output_format, meta_data, meta)
        except:
            log.write(f"Error in using fast_to_vcf. Falling back to original implementation.",verbose=verbose)
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
    elif fmt in get_formats_list() :      
        # tabular 
        log.write(" -"+fmt+" format will be loaded...",verbose=verbose)
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        _print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True)
        
        # determine if gzip or not / create path for output
        if gzip ==True and tab_fmt not in non_gzip_tab_fmt:
            path = path + "."+suffix+".{}.gz".format(tab_fmt)
        else:
            path = path + "."+suffix+".{}".format(tab_fmt)

        yaml_path = path + "."+suffix+".{}-meta.yaml".format(tab_fmt)
        log.write(" -Output path:",path, verbose=verbose) 

        sumstats,to_csvargs = _configure_output_cols_and_kwargs(sumstats, rename_dictionary, cols, no_status, path, meta_data, to_csvargs, log, verbose)
        
        log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)
        
        _write_tabular(sumstats,rename_dictionary, path, tab_fmt, to_csvargs, to_tabular_kwargs, log, verbose, gzip)
        
        if tab_fmt not in non_md5sum_tab_fmt and "@" not in path:
            if md5sum == True: 
                # write a md5sum file
                md5_value = md5sum_file(path,log,verbose)
            else:
                # calculate md5sum without saveing a file
                md5_value = calculate_md5sum_file(path)
        else:
            md5_value = "NA"

        ## update ssf-style meta data and export to yaml file
        _configure_ssf_meta(sumstats, fmt, ssfmeta, meta, meta_data, path, md5_value, yaml_path, log, verbose)
        
        return path  
    
####################################################################################################################
def _write_tabular(sumstats,rename_dictionary, path, tab_fmt, to_csvargs, to_tabular_kwargs, log, verbose, gzip):
    if tab_fmt=="tsv" or tab_fmt=="csv":
        # pandas automatically detects compression when path ends with .gz
        if "@" in path:
            chr_header = rename_dictionary["CHR"]
            log.write(f"  -@ detected: writing each chromosome to a single file...",verbose=verbose)
            log.write("  -Chromosomes:{}...".format(list(sumstats["CHR"].unique())),verbose=verbose)
            for single_chr in list(sumstats["CHR"].unique()):
                single_path = path.replace("@","{}".format(single_chr))
                sumstats.loc[sumstats[chr_header]==single_chr,:].to_csv(single_path, index=False, **to_csvargs)
        else:
            sumstats.to_csv(path, index=False, **to_csvargs)

    elif tab_fmt=="parquet":
        # When partition_cols is provided, PyArrow expects a directory path, not a file path
        if to_tabular_kwargs and "partition_cols" in to_tabular_kwargs and to_tabular_kwargs["partition_cols"]:
            # Convert file path to directory path by removing .parquet extension
            if path.endswith(".parquet"):
                path = path[:-8]  # Remove ".parquet" extension
            # Create directory if it doesn't exist
            os.makedirs(path, exist_ok=True)
            log.write("  -Partitioned parquet will be written to directory: {}".format(path), verbose=verbose)
        sumstats.to_parquet(path, index=None, **to_tabular_kwargs)


def fast_to_vcf(dataframe, path, vcf_header, output_format, meta_data, meta):
    # Get the columns in the right order and convert to numpy
    df_numpy = dataframe[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO'] + output_format].to_numpy()

    sep = '\t'
    QUAL = "."
    FILTER = "PASS"
    FORMAT = ":".join(output_format)
    format_format = ':'.join(['%s']*len(output_format))

    single_row_format = f'%s %s %s %s %s {QUAL} {FILTER} %s {FORMAT} {format_format}'

    out_string = vcf_header
    out_string += sep.join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]]) + "\n"
    fmt = sep.join(single_row_format.split(' ')) # build formatting for one single row
    fmt = '\n'.join([fmt]*dataframe.shape[0]) # add newline and replicate the formatting for all rows
    out_string += fmt % tuple(df_numpy.ravel()) # flatten the array and then apply formatting
    out_string += '\n'

    with open(path, 'w') as f:
        f.write(out_string)

####################################################################################################################
def _configure_output_cols_and_kwargs(sumstats, rename_dictionary, cols, no_status, path, meta_data, to_csvargs, log, verbose):
    # grab format cols that exist in sumstats
    ouput_cols=[]
    for i in sumstats.columns:
        if i in rename_dictionary.keys():
            ouput_cols.append(i)  
    
    # + additional cols and remove duplicated
    ouput_cols_final = []
    for i in ouput_cols + cols:
        if i not in ouput_cols_final:
            ouput_cols_final.append(i)
    
    # remove STATUS
    try:
        if no_status == True:
            ouput_cols_final.remove("STATUS")
    except:
        pass
    
    #filter and rename to target fromat headers
    sumstats = sumstats[ouput_cols_final]
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


def _configure_ssf_meta(sumstats, fmt, ssfmeta, meta, meta_data, path, md5_value, yaml_path, log, verbose):
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
        log.write(" -Exporting SSF-style meta data to {}".format(yaml_path),verbose=verbose) 
        with open(yaml_path, 'w') as outfile:
            yaml.dump(sumstats_meta_copy, outfile)



def _output_bed_like(sumstats, path, fmt, suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose):
    sumstats = sumstats[ouput_cols]
    path = path + "."+suffix
    log.write(" -Output columns: {}".format(",".join(sumstats.columns)),verbose=verbose)
    log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)
    if to_csvargs is None:
        to_csvargs = {}
    sumstats.to_csv(path,sep="\t",index=False,header=None,**to_csvargs)
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
    
    # Pre-compute common values to avoid repeated operations
    pos = sumstats["POS"]
    nea_len = sumstats["NEA"].str.len()
    nea_str = sumstats["NEA"].astype("string")
    ea_str = sumstats["EA"].astype("string")
    ea_slice = sumstats["EA"].str.slice(start=1)
    nea_slice = sumstats["NEA"].str.slice(start=1)
    
    # Initialize START and END columns
    sumstats["START"] = pd.NA
    sumstats["END"] = pd.NA
    
    if fmt=="bed":
        # SNPs: 0-based coordinates
        sumstats.loc[is_snp,"START"] = pos[is_snp] - 1
        sumstats.loc[is_snp,"END"] = pos[is_snp] - 1 + nea_len[is_snp]
        sumstats.loc[is_snp,"NEA/EA"] = nea_str[is_snp] + "/" + ea_str[is_snp]
        
        # Insertions: START = END = POS (0-based)
        sumstats.loc[is_insert,"START"] = pos[is_insert]
        sumstats.loc[is_insert,"END"] = pos[is_insert]
        sumstats.loc[is_insert,"NEA/EA"] = "-/" + ea_slice[is_insert]
        
        # Deletions: START = POS, END = POS + len(NEA) - 1 (0-based)
        sumstats.loc[is_delete,"START"] = pos[is_delete]
        sumstats.loc[is_delete,"END"] = pos[is_delete] + nea_len[is_delete] - 1
        sumstats.loc[is_delete,"NEA/EA"] = nea_slice[is_delete] + "/-"
        
        sumstats["STRAND"] = "+"
        
    elif fmt=="vep":
        # SNPs: 1-based coordinates
        nea_len_snp = nea_len[is_snp] - 1
        sumstats.loc[is_snp,"START"] = pos[is_snp] + nea_len_snp
        sumstats.loc[is_snp,"END"] = pos[is_snp] + nea_len_snp
        sumstats.loc[is_snp,"NEA/EA"] = nea_str[is_snp] + "/" + ea_str[is_snp]
        
        # Insertions: START = POS + 1, END = POS (VEP convention: START > END)
        sumstats.loc[is_insert,"START"] = pos[is_insert] + 1
        sumstats.loc[is_insert,"END"] = pos[is_insert]
        sumstats.loc[is_insert,"NEA/EA"] = "-/" + ea_slice[is_insert]
        
        # Deletions: START = POS + 1, END = POS + len(NEA) - 1 (1-based)
        nea_len_del = nea_len[is_delete] - 1
        sumstats.loc[is_delete,"START"] = pos[is_delete] + 1
        sumstats.loc[is_delete,"END"] = pos[is_delete] + nea_len_del
        sumstats.loc[is_delete,"NEA/EA"] = nea_slice[is_delete] + "/-"
        
        sumstats["STRAND"] = "+"
        
    elif fmt=="annovar":
        # SNPs: 1-based coordinates
        sumstats.loc[is_snp,"START"] = pos[is_snp]
        sumstats.loc[is_snp,"END"] = pos[is_snp] - 1 + nea_len[is_snp]
        sumstats.loc[is_snp,"NEA_out"] = nea_str[is_snp]
        sumstats.loc[is_snp,"EA_out"] = ea_str[is_snp]
        
        # Insertions: START = END = POS (ANNOVAR convention)
        sumstats.loc[is_insert,"START"] = pos[is_insert]
        sumstats.loc[is_insert,"END"] = pos[is_insert]
        sumstats.loc[is_insert,"NEA_out"] = "-"
        sumstats.loc[is_insert,"EA_out"] = ea_slice[is_insert]
        
        # Deletions: START = POS, END = POS - 1 + len(NEA) (1-based)
        sumstats.loc[is_delete,"START"] = pos[is_delete]
        sumstats.loc[is_delete,"END"] = pos[is_delete] - 1 + nea_len[is_delete]
        sumstats.loc[is_delete,"NEA_out"] = nea_slice[is_delete]
        sumstats.loc[is_delete,"EA_out"] = "-"
    
    # Convert to Int64 at the end (more efficient than per-assignment)
    sumstats["START"] = sumstats["START"].astype("Int64")
    sumstats["END"] = sumstats["END"].astype("Int64")
    return sumstats

def _process_vcf_header(sumstats, meta, meta_data, build, log, verbose):
    
    log.write(" -Creating VCF file header...",verbose=verbose)
    log.write("  -VCF header contig build:"+str(build),verbose=verbose) 
    
    # calculate meta data
    # Match patterns using integer arithmetic
    from gwaslab.info.g_vchange_status import status_match
    # Pattern: digit 4=0, digit 5=0-3, digit 6=0-2, digit 7=0-4
    harmonised = sum(status_match(sumstats["STATUS"], 4, [0]) & 
                     status_match(sumstats["STATUS"], 5, [0,1,2,3]) &
                     status_match(sumstats["STATUS"], 6, [0,1,2]) &
                     status_match(sumstats["STATUS"], 7, [0,1,2,3,4]))
    # Pattern: digit 4=0, digit 5=0-3, digit 6=1-2, digit 7=2 or 4
    switchedalleles = sum(status_match(sumstats["STATUS"], 4, [0]) &
                          status_match(sumstats["STATUS"], 5, [0,1,2,3]) &
                          status_match(sumstats["STATUS"], 6, [1,2]) &
                          status_match(sumstats["STATUS"], 7, [2,4]))

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
