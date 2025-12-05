import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
# pandas.api.types.is_int64_dtype
# pandas.api.types.is_categorical_dtype

dtype_dict ={
    "SNPID":["string","object"],
    "rsID": ["string","object"],
    "CHR":  ["Int64","int64","int32","Int32","int"],
    "POS":  ["int64","Int64","int"],
    "EA":   ["category"],  
    "NEA":["category"],  
    "REF":["category"],  
    "ALT":["category"],  
    "BETA":["float64"],
    "BETA_95L":["float64"],
    "BETA_95U":["float64"],
    "SE":["float64"],
    "N":["Int64","int64","int32","Int32","int"],
    "N_CASE":["Int64","int64","int32","Int32","int"],
    "N_CONTROL":["Int64","int64","int32","Int32","int"],
    "OR":["float64"],
    "OR_95L":["float64"],
    "OR_95U":["float64"],
    "HR":["float64"],
    "HR_95L":["float64"],
    "HR_95U":["float64"],
    "P":["float64"],
    "MLOG10P":["float64"],
    "Z":["float64"],
    "F":["float64"],
    "T":["float64"],
    "TEST":["string","object","category"],
    "CHISQ":["float64"],
    "I2":["float64"],
    "P_HET":["float64"],
    "SNPR2":["float64"],
    "EAF":["float64","float","float32"],
    "NEAF":["float64","float","float32"],
    "MAF":["float64","float","float32"],
    "INFO":["float64","float","float32"],
    "DOF":["Int64","int64","int32","Int32","int"],  
    "STATUS":["category"],  
    "DIRECTION":["string","object"],
    'PIP'               :["float64","float","float32"],
    'CREDIBLE_SET_INDEX':["Int64","int64","int32","Int32","int"],
    'N_SNP'             :["Int64","int64","int32","Int32","int"],
    'LOCUS'             :["string","object","category"],
    'STUDY'             :["string","object","category"],
    'BETA_RANDOM' :["float64"],
    'SE_RANDOM' :["float64"],
    'Z_RANDOM' :["float64"],
    'P_RANDOM' :["float64"]
    }

def check_datatype(sumstats, verbose=True, log=Log()):
    """
    Inspect and report data types for all columns in a summary statistics
    DataFrame, comparing each against expected dtypes.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table to check.
    verbose : bool, default True
        Whether to print log messages to stdout.
    log : gwaslab.g_Log.Log
        Logger used for structured messages and warnings.

    Behavior
    --------
    - Logs aligned lists of column names, their dtypes, and verification flags
      ("T" for match, "F" for mismatch, "NA" for unknown header).
    - Emits targeted fix suggestions for CHR, POS, allele, and ID columns when
      their dtypes are incompatible.
    - For statistics columns present (e.g., BETA, SE, Z, P, OR, HR, CHISQ,
      INFO), warns to use `Sumstats.check_sanity()` only when those columns have
      incompatible dtypes.
    """
    
    try:
        headers = []
        dtypes = []
        verified = []
        raw_verified =[]

        for header,dtype in sumstats.dtypes.items():
            width = max(len(header),len(str(dtype)))
            
            header_fix_length = header + " "*(width- len(header) )
            
            dtype_fix_length  = str(dtype) + " "*(width- len(str(dtype)))
            
            verified_str = verify_datatype(header, dtype)
            verified_fix_length  = verified_str + " " *(width- len(verified_str))
            
            headers.append(format(header_fix_length))
            dtypes.append((str(dtype_fix_length)))
            verified.append(verified_fix_length)
            if verified_str == "F":
                raw_verified.append(header)

        log.write(" -Column  :", " ".join(headers), verbose=verbose)
        log.write(" -DType   :", " ".join(dtypes), verbose=verbose)
        log.write(" -Verified:", " ".join(verified), verbose=verbose)

        if len(raw_verified)>0:
            log.warning("Columns with possibly incompatible dtypes: {}".format(",".join(raw_verified)), verbose=verbose)
            try:
                if "CHR" in raw_verified:
                    log.warning("Consider using Sumstats.fix_chr() to fix CHR dtype", verbose=verbose)
                if "POS" in raw_verified:
                    log.warning("Consider using Sumstats.fix_pos() to fix POS dtype", verbose=verbose)
                if ("EA" in raw_verified) or ("NEA" in raw_verified) or ("REF" in raw_verified) or ("ALT" in raw_verified):
                    log.warning("Consider using Sumstats.fix_allele() to fix allele dtypes", verbose=verbose)
                if ("SNPID" in raw_verified) or ("rsID" in raw_verified):
                    log.warning("Consider using Sumstats.fix_id() to fix ID dtypes", verbose=verbose)
                stats_candidates = [
                    "BETA","SE","Z","P","MLOG10P","OR","OR_95L","OR_95U",
                    "HR","HR_95L","HR_95U","CHISQ","F","T","SNPR2","INFO"
                ]
                present_stats = [c for c in stats_candidates if c in sumstats.columns]
                bad_stats = [c for c in present_stats if verify_datatype(c, sumstats.dtypes[c]) == "F"]
                if len(bad_stats) > 0:
                    for c in bad_stats:
                        log.warning("Consider using Sumstats.check_sanity() to check {} statistics".format(c), verbose=verbose)
            except:
                pass
    except:
        pass

def verify_datatype(header, dtype):
    """
    Verify a single column's dtype against expected dtypes.

    Parameters
    ----------
    header : str
        Column name to check.
    dtype : object
        The pandas dtype object or its string representation.

    Returns
    -------
    str
        "T" if dtype matches one of the expected entries, "F" if it does not,
        and "NA" if the column is not recognized.
    """

    if header in dtype_dict.keys():
        if str(dtype) in dtype_dict[header]:
            return "T"
        else:
            return "F"
    else:
        return "NA"

def quick_convert_datatype(sumstats, log, verbose):
    """
    Attempt to convert columns with incompatible dtypes to the first expected
    dtype for that column, logging successes and failures.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table to convert.
    log : gwaslab.g_Log.Log
        Logger used for messages.
    verbose : bool
        Whether to print log messages to stdout.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns converted where possible.
    """
    for col in sumstats.columns:
        if col in dtype_dict.keys():
            if str(sumstats[col].dtypes) not in dtype_dict[col]:
                datatype=dtype_dict[col][0]
                log.write(" -Trying to convert datatype for {}: {} -> {}...".format(col, str(sumstats[col].dtypes), datatype), end="" ,verbose=verbose)
                try:
                    sumstats[col] = sumstats[col].astype(datatype)
                    log.write("{}".format(datatype),show_time=False, verbose=verbose)
                except:
                    log.write("Failed...",show_time=False,verbose=verbose)
                    pass
    return sumstats

def check_dataframe_shape(sumstats, log, verbose):
    """
    Log DataFrame shape and estimated memory usage in megabytes.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table.
    log : gwaslab.g_Log.Log
        Logger used for messages.
    verbose : bool
        Whether to print log messages to stdout.
    """
    memory_in_mb = sumstats.memory_usage().sum()/1024/1024
    try:
        log.write(" -Current Dataframe shape : {} x {} ; Memory usage: {:.2f} MB".format(len(sumstats),len(sumstats.columns),memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Dataframe shape...")
    
def check_dataframe_memory_usage(sumstats, log, verbose):
    """
    Log DataFrame memory usage in megabytes.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table.
    log : gwaslab.g_Log.Log
        Logger used for messages.
    verbose : bool
        Whether to print log messages to stdout.
    """
    memory_in_mb = sumstats.memory_usage().sum()/1024/1024
    try:
        log.write(" -Current Dataframe memory usage: {:.2f} MB".format(memory_in_mb), verbose=verbose)
    except:
        log.warning("Error: cannot get Memory usage...")

def check_datatype_for_cols(sumstats, cols=None, verbose=True, log=Log(), fix=False, **fix_kwargs):
    """
    Verify dtypes for a specified subset of columns and emit fix suggestions.

    Parameters
    ----------
    sumstats : pandas.DataFrame
        Summary statistics table to check.
    cols : list[str] or None, default None
        Column names to verify. If None, no columns are checked.
    verbose : bool, default True
        Whether to print log messages to stdout.
    log : gwaslab.g_Log.Log
        Logger used for messages and warnings.

    Behavior
    --------
    - Collects failing columns and logs targeted fix suggestions for CHR, POS,
      allele, and ID columns.
    - For present statistics columns among the specified `cols`, warns to use
      `Sumstats.check_sanity()` only when those columns have incompatible
      dtypes.
    - Raises `ValueError` listing failing columns to encourage corrective
      actions via Sumstats methods.
    """
    if cols is None:
        cols = []
    try:
        failed = []

        for header,dtype in sumstats.dtypes.items():
            if header in cols:
                verified_str = verify_datatype(header, dtype)
                if verified_str=="F":
                    failed.append(header)
    except:
        pass

    if len(failed) >0:
        log.warning("Data types were not fixed for : {} ".format(",".join(failed)))
        try:
            if fix is True:
                try:
                    from gwaslab.qc.qc_fix_sumstats import fixchr, fixpos, fixallele, fixID
                except Exception:
                    fixchr = None; fixpos = None; fixallele = None; fixID = None

                if "CHR" in failed and fixchr is not None:
                    sumstats = fixchr(sumstats, log=log, verbose=verbose, **{k:v for k,v in fix_kwargs.items() if k in ["remove","add_prefix","x","y","mt","chrom_list","minchr"]})
                if "POS" in failed and fixpos is not None:
                    sumstats = fixpos(sumstats, log=log, verbose=verbose, **{k:v for k,v in fix_kwargs.items() if k in ["remove","lower_limit","upper_limit","limit"]})
                if (set(["EA","NEA","REF","ALT"]) & set(failed)) and fixallele is not None:
                    sumstats = fixallele(sumstats, log=log, verbose=verbose, **{k:v for k,v in fix_kwargs.items() if k in ["remove"]})
                if ("SNPID" in failed or "rsID" in failed) and fixID is not None:
                    sumstats = fixID(sumstats, log=log, verbose=verbose, **{k:v for k,v in fix_kwargs.items() if k in ["fixprefix","fixchrpos","fixid","fixeanea","fixeanea_flip","fixsep","reversea","overwrite","forcefixid"]})

                # Re-verify targeted columns after attempted fixes
                re_failed = []
                for header,dtype in sumstats.dtypes.items():
                    if header in cols:
                        verified_str = verify_datatype(header, dtype)
                        if verified_str=="F":
                            re_failed.append(header)

                if len(re_failed) == 0:
                    log.write(" -Dtype fixes applied successfully for requested columns.", verbose=verbose)
                    return sumstats
                else:
                    log.warning("Some columns still have incompatible dtypes after fixes: {}".format(",".join(re_failed)), verbose=verbose)

            # Suggest manual next steps if fix was False or fixes failed
            if "CHR" in failed:
                log.warning("Consider using Sumstats.fix_chr() to fix CHR dtype", verbose=verbose)
            if "POS" in failed:
                log.warning("Consider using Sumstats.fix_pos() to fix POS dtype", verbose=verbose)
            if ("EA" in failed) or ("NEA" in failed) or ("REF" in failed) or ("ALT" in failed):
                log.warning("Consider using Sumstats.fix_allele() to fix allele dtypes", verbose=verbose)
            if ("SNPID" in failed) or ("rsID" in failed):
                log.warning("Consider using Sumstats.fix_id() to fix ID dtypes", verbose=verbose)

            stats_candidates = [
                "BETA","SE","Z","P","MLOG10P","OR","OR_95L","OR_95U",
                "HR","HR_95L","HR_95U","CHISQ","F","T","SNPR2","INFO"
            ]
            present_stats = [c for c in stats_candidates if c in sumstats.columns]
            bad_stats = [c for c in present_stats if verify_datatype(c, sumstats.dtypes[c]) == "F"]
            if len(bad_stats) > 0:
                for c in bad_stats:
                    log.warning("Consider using Sumstats.check_sanity() to check {} statistics".format(c), verbose=verbose)
        except:
            pass

        if fix is True:
            raise ValueError("Failed to fix dtypes for requested columns: {}".format(",".join(failed)))
        else:
            raise ValueError("Please fix dtypes using Sumstats methods (fix_chr, fix_pos, fix_allele, fix_id) for: {}. If statistics columns are present, consider Sumstats.check_sanity().".format(",".join(failed)))
