from typing import TYPE_CHECKING, Union, Tuple, List, Any, Optional
import re
import gc
import pandas as pd
import numpy as np
from multiprocessing import  Pool
from functools import partial

from gwaslab.info.g_vchange_status import vchange_status, set_status_digit, status_match
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version
from gwaslab.info.g_vchange_status import STATUS_CATEGORIES
from gwaslab.info.g_vchange_status import match_status

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chain
from gwaslab.bd.bd_common_data import NA_STRINGS
from gwaslab.bd.bd_sex_chromosomes import Chromosomes

from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_pattern import RSID_PATTERN, SNPID_PATTERN_EXTRACT, CHR_PATTERN_EXTRACT, FLAGS, SNPID_SEP_PATTERN, CHR_PREFIX_PATTERN, SNPID_PATTERN_STRIP
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.qc.qc_reserved_headers import DEFAULT_COLUMN_ORDER
from gwaslab.util.util_in_fill_data import fill_extreme_mlog10p

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

# Main QC and data fixing functions in this module:
# - _process_build: Process and validate genome build information
# - _set_build: Set genome build for sumstats
# - _fix_ID: Fix SNPID/rsID format, extract CHR/POS/EA/NEA from IDs
# - _rsid_to_chrpos: Convert rsID to CHR:POS coordinates (if implemented)
# - _remove_dup: Remove duplicate or multiallelic variants
# - _fix_chr: Standardize chromosome notation (handle X, Y, MT, prefixes)
# - _fix_pos: Validate and standardize genomic positions
# - _fix_allele: Validate and standardize allele representations (EA/NEA)
# - _flip_allele_stats: Adjust statistics when alleles are flipped/reverse complemented
# - _sort_coordinate: Sort variants by genomic coordinates (CHR, then POS)
# - _sort_column: Reorder DataFrame columns according to standard order

@with_logging(
        start_to_msg= "check SNPID/rsID",
        finished_msg= "checking SNPID/rsID",
        start_function= ".fix_id()",
        start_cols=[],
        check_dtype=False,
        fix=False
)
def _fix_ID(sumstats_obj: Union['Sumstats', pd.DataFrame],
       snpid: str = "SNPID", rsid: str = "rsID", chrom: str = "CHR", pos: str = "POS", nea: str = "NEA", ea: str = "EA", status: str = "STATUS", fixprefix: bool = False,
       fixchrpos: bool = False, fixid: bool = False, fixeanea: bool = False, fixeanea_flip: bool = False, fixsep: bool = False, reversea: bool = False,
       overwrite: bool = False, verbose: bool = True, forcefixid: bool = False, log: Log = Log()) -> pd.DataFrame:  
    
    '''
    Fix various aspects of genomic data including SNPID, rsID, chromosome positions, and allele information.
    
    This function performs multiple data quality checks and fixes: (1) validates and fixes SNPID format
    (CHR:POS:NEA:EA pattern), (2) extracts and fixes chromosome and position from SNPID or rsID, (3) validates
    rsID format (rsxxxxxx pattern), (4) extracts and fixes EA and NEA from SNPID, (5) standardizes separators
    and removes prefixes in SNPID, and (6) generates new SNPID from available CHR, POS, EA, NEA data.
    
    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to fix.
    fixprefix : bool
        Whether to remove 'chr' prefix in SNPID.
    fixchrpos : bool
        Whether to fix chromosome and position from SNPID.
    fixid : bool
        Whether to generate new SNPID from available data.
    fixeanea : bool
        Whether to fix EA and NEA from SNPID.
    fixeanea_flip : bool
        Whether to flip EA and NEA during fixing.
    fixsep : bool
        Whether to standardize separators in SNPID.
    reversea : bool
        Whether to reverse alleles in SNPID.
    overwrite : bool
        Whether to overwrite existing values.
    verbose : bool, optional
        Whether to print progress.
    forcefixid : bool
        Whether to force fix even without status check.

    Returns
    -------
    pd.DataFrame
        Modified sumstats.data with fixed data.
    '''
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called during initialization - no Sumstats object yet
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False

    ############################  checking datatype ###################################################  
    if rsid in sumstats.columns:
        log.log_operation("Checking rsID data type...", verbose=verbose)
        if sumstats[rsid].dtype != "string":
            log.log_datatype_change("rsID", str(sumstats[rsid].dtype), "string", success=True, verbose=verbose)
            sumstats[rsid] = sumstats[rsid].astype("string")

    if snpid in sumstats.columns:
        log.log_operation("Checking SNPID data type...", verbose=verbose)
        if sumstats[snpid].dtype != "string":
            log.log_datatype_change("SNPID", str(sumstats[snpid].dtype), "string", success=True, verbose=verbose)
            sumstats[snpid] = sumstats[snpid].astype("string")
            
    ############################  checking string NA ###################################################
    log.write(" -Checking NA strings :{}".format(",".join(NA_STRINGS)),verbose=verbose)  
    if snpid in sumstats.columns:  
        log.write(" -Checking if SNPID contains NA strings...",verbose=verbose)
        is_snpid_string_na = sumstats[snpid].isin(NA_STRINGS)
        if is_snpid_string_na.sum() > 0:
            log.log_operation("Converting {} NA strings in SNPID to pd.NA".format(is_snpid_string_na.sum()), indent=1, verbose=verbose)
            sumstats.loc[is_snpid_string_na ,snpid] = pd.NA

    if rsid in sumstats.columns: 
        log.write(" -Checking if rsID contains NA strings...",verbose=verbose)
        is_rsid_string_na = sumstats[rsid].isin(NA_STRINGS)
        if is_rsid_string_na.sum() > 0:
            log.log_operation("Converting {} NA strings in rsID to pd.NA".format(is_rsid_string_na.sum()), indent=1, verbose=verbose)
            sumstats.loc[is_rsid_string_na ,rsid] = pd.NA
    ############################  checking ###################################################  
    if snpid in sumstats.columns:  
        log.write(" -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)",verbose=verbose)
        # check if SNPID is CHR:POS:EA:NEA
        is_chrposrefalt = sumstats[snpid].str.match(SNPID_PATTERN_EXTRACT, flags=FLAGS, na=False)
        # check if SNPID is NA
        is_snpid_na = sumstats[snpid].isna()

        # Ensure STATUS is integer (not Categorical) before assignment
        if sumstats[status].dtype.name == 'category':
            sumstats[status] = sumstats[status].astype(str).astype(int).astype('Int64')
        
        # change STATUS code
        sumstats.loc[ is_chrposrefalt,status] = vchange_status(sumstats.loc[ is_chrposrefalt,status],3 ,"975" ,"630")
        sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status] = vchange_status(sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status],3 ,"975" ,"842")
        
    if rsid in sumstats.columns: 
        log.write(" -Checking if rsID is rsxxxxxx...", verbose=verbose)
        is_rsid = sumstats[rsid].str.match(RSID_PATTERN, flags=FLAGS, na=False)
        
        sumstats.loc[ is_rsid,status] = vchange_status(sumstats.loc[ is_rsid,status], 3, "986","520")
        sumstats.loc[~is_rsid,status] = vchange_status(sumstats.loc[~is_rsid,status], 3, "986","743")
        
        log.write(" -Checking if CHR:POS:NEA:EA is mixed in rsID column ...", verbose=verbose)
        is_rs_chrpos = sumstats[rsid].str.match(SNPID_PATTERN_EXTRACT, flags=FLAGS, na=False)
        
        log.write(" -Number of CHR:POS:NEA:EA mixed in rsID column :", is_rs_chrpos.sum(), verbose=verbose)
        log.write(" -Number of Unrecognized rsID :", len(sumstats) - is_rs_chrpos.sum() - is_rsid.sum(), verbose=verbose)
        log.write(" -A look at the unrecognized rsID :",set(sumstats.loc[(~is_rsid)&(~is_rs_chrpos),rsid].head()),"...", verbose=verbose) 
      
    ############################  fixing chr pos###################################################  
    if reversea == True:
        if snpid in sumstats.columns: 
            log.write(" -Reversing Alleles in SNPID...", verbose=verbose)
            to_fix = is_chrposrefalt
            to_fix_num = to_fix.sum()
            if to_fix_num>0: log.log_variants_count(to_fix_num, "variants could be reversed", verbose=verbose)
            extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
            sumstats.loc[to_fix, snpid] = extracted["CHR"] + ":" + extracted["POS"] + ":" + extracted["EA"] + ":" + extracted["NEA"]
            
    if fixchrpos == True:
    # from snpid or rsid, extract CHR:POS to fix CHR and POS    
        if snpid in sumstats.columns: 
            log.write(" -Fixing CHR and POS...", verbose=verbose)
            if overwrite is True: 
                log.write(" -Overwrite is applied...", verbose=verbose)
                # fix all
                to_fix = is_chrposrefalt
                
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                #fix variants with chr and pos being NA
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                to_fix_num = to_fix.sum()
                if to_fix_num: 
                    log.log_variants_count(to_fix_num, "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)
            
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                log.log_column_added("CHR", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                to_fix_num = to_fix.sum()
                if to_fix_num>0: 
                    log.log_variants_count(to_fix_num, "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)
            
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                log.log_column_added("CHR and POS", verbose=verbose)
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna() 
                to_fix_num = to_fix.sum()
                if to_fix_num>0: 
                    log.log_variants_count(to_fix_num, "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)     
                
            else:
                log.log_column_added("CHR and POS", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt
                to_fix_num = to_fix.sum()
                if to_fix_num>0: 
                    log.log_variants_count(to_fix_num, "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)   
                    
            if to_fix.sum()>0:
                log.log_formula("CHR and POS", "from SNPID", source_columns=["SNPID"], verbose=verbose)
                # format and qc filled chr and pos
                extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
                # Convert to regular pandas Series to avoid PyArrow type issues
                sumstats.loc[to_fix, chrom] = extracted["CHR"].astype("string")
                sumstats.loc[to_fix, pos] = pd.to_numeric(extracted["POS"], errors='coerce').astype('Int64')
                
                #sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,snpid].str.split(':|_|-').str[0].str.strip("chrCHR").astype("string")
                #sumstats.loc[to_fix,pos] =np.floor(pd.to_numeric(sumstats.loc[to_fix,snpid].str.split(':|_|-').str[1], errors='coerce')).astype('Int64')
                # change status
                #sumstats.loc[to_fix,status] = vchange_status(sumstats.loc[to_fix,status], 4, "98765432","00000000")    
        
        if rsid in sumstats.columns:
            log.write(" -Fixing CHR and POS using chr:pos:ref:alt format variants in rsID column...", verbose=verbose)
            if overwrite is True: 
                log.write(" -Overwrite is applied...", verbose=verbose)
                to_fix = is_rs_chrpos
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if to_fix.sum()>0: 
                    log.log_variants_count(to_fix.sum(), "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                log.log_column_added("CHR", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if to_fix.sum()>0: 
                    log.log_variants_count(to_fix.sum(), "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                log.log_column_added("CHR and POS", verbose=verbose)
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna() 
                if to_fix.sum()>0: 
                    log.log_variants_count(to_fix.sum(), "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)
            else:
                log.log_column_added("CHR and POS", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos
                if to_fix.sum()>0: 
                    log.log_variants_count(to_fix.sum(), "variants could be fixed", verbose=verbose)
                else:
                    log.log_operation("No fixable variants", verbose=verbose)   
            
            if to_fix.sum()>0:

                log.log_formula("CHR and POS", "from rsID", source_columns=["rsID"], verbose=verbose)
                # Convert to regular pandas Series to avoid PyArrow type issues
                sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,rsid].str.split(':|_|-',n=2).str[0].astype("string")
                sumstats.loc[to_fix,pos] = pd.to_numeric(sumstats.loc[to_fix,rsid].str.split(':|_|-',n=2).str[1], errors='coerce').astype('Int64')
                #sumstats.loc[to_fix,pos] = np.floor(pd.to_numeric(sumstats.loc[to_fix,rsid].str.split(':|_|-',x).get(1), errors='coerce')).astype('Int64')
                #sumstats.loc[to_fix,status] = vchange_status(sumstats.loc[to_fix,status], 4, "98765432","00000000").astype("string")  
                
    ############################  fixing chr pos###################################################   
    if fixeanea == True:
        log.warning("gwaslab assumes SNPID is in the format of CHR:POS:NEA:EA / CHR:POS:REF:ALT", verbose=verbose)
        if overwrite is True:
            log.write(" -Overwrite mode is applied...", verbose=verbose)
            to_fix = is_chrposrefalt
        elif (nea in sumstats.columns) and (ea in sumstats.columns):
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
        elif (nea in sumstats.columns) and (ea not in sumstats.columns):
            log.log_column_added("EA", verbose=verbose)
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
        elif (nea not in sumstats.columns) and (ea in sumstats.columns):
            log.log_column_added("NEA", verbose=verbose)
            sumstats[nea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
        else:
            log.log_column_added("EA and NEA", verbose=verbose)
            sumstats[nea]=pd.Series(dtype="string")
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt
            if to_fix.sum()>0:
                log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...", verbose=verbose)
    #                
        if to_fix.sum()>0:
            log.log_formula("EA and NEA", "from SNPID", source_columns=["SNPID"], verbose=verbose)
    #        
            if fixeanea_flip == True:
                log.write(" -Flipped : CHR:POS:NEA:EA -> CHR:POS:EA:NEA ", verbose=verbose)
                extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
                # Convert to regular pandas Series to avoid PyArrow type issues
                sumstats.loc[to_fix, ea] = extracted["NEA"].astype("string")
                sumstats.loc[to_fix, nea] = extracted["EA"].astype("string")
            else:
                log.write(" -Chr:pos:a1:a2...a1->EA , a2->NEA ", verbose=verbose)
                extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
                # Convert to regular pandas Series to avoid PyArrow type issues
                sumstats.loc[to_fix, ea] = extracted["EA"].astype("string")
                sumstats.loc[to_fix, nea] = extracted["NEA"].astype("string")
    
    #        #to_change_status = sumstats[status].str.match(r"\w\w\w[45]\w\w\w")
    #        #sumstats.loc[to_fix&to_change_status,status] = vchange_status(sumstats.loc[to_fix&to_change_status,status],4,"2")  
    #        #sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]).astype("string")
    #        #sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[1]).astype("Int64")
    
    ############################  fixing id ################################################### 
    if fixsep == True:
        if snpid in sumstats.columns: 
            log.write(' -Replacing separators in SNPID with ":" ...', verbose=verbose)
            sumstats[snpid] = sumstats[snpid].str.replace(SNPID_SEP_PATTERN, ":", regex=True, flags=FLAGS)
    
    if fixprefix == True:
        if snpid in sumstats.columns: 
            log.write(' -Removing /^chr/ in SNPID ...', verbose=verbose)
            sumstats[snpid] = sumstats[snpid].str.replace(CHR_PREFIX_PATTERN, '', regex=True, flags=FLAGS)

    if fixid == True:
        if snpid not in sumstats.columns: 
        # initiate a SNPID column
            sumstats[snpid]=pd.Series(dtype="string")
        
        if (rsid in sumstats.columns) and (sum(is_rs_chrpos)>0) :
            sumstats[snpid]= sumstats.loc[is_rs_chrpos,rsid]
            
        if (chrom in sumstats.columns) and (pos in sumstats.columns):
            #only fix when CHR and POS is available
            pre_number=sum(sumstats[snpid].isna())                                
            
            if overwrite is False:
                #fix empty - digit 4 = 0
                to_fix = sumstats[snpid].isna() & status_match(sumstats[status], 4, [0])
            else:
                #fix all - digit 4 = 0
                to_fix = status_match(sumstats[status], 4, [0])
            
            if (ea in sumstats.columns) and (nea in sumstats.columns):
            # when ea and nea is available  -> check status -> fix to chr:pos:nea:ea 
                
                # digit 4 = 0
                matched_index = status_match(sumstats[status], 4, [0])
                to_part_fix = matched_index & to_fix 
                
                # digit 5 in [0,1,2,3] AND digit 6 in [0,1,2,6,7] AND digit 7 in [0,1,2,3,4]
                digit_5_match = status_match(sumstats[status], 5, [0, 1, 2, 3])
                digit_6_match = status_match(sumstats[status], 6, [0, 1, 2, 6, 7])
                digit_7_match = status_match(sumstats[status], 7, [0, 1, 2, 3, 4])
                matched_index = digit_5_match & digit_6_match & digit_7_match
                if forcefixid is True:
                    matched_index = to_fix
                to_full_fix = matched_index & to_fix 
                
                
                if to_part_fix.sum() > 0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                if to_full_fix.sum() > 0:
                    sumstats.loc[to_full_fix,snpid] = sumstats.loc[to_full_fix,chrom].astype("string") + ":"+sumstats.loc[to_full_fix,pos].astype("string") +":"+ sumstats.loc[to_full_fix,nea].astype("string") +":"+ sumstats.loc[to_full_fix,ea].astype("string")
                log.log_formula("SNPID", "from CHR:POS", source_columns=["CHR", "POS"], verbose=verbose)
                log.log_formula("SNPID", "from CHR:POS:NEA:EA", source_columns=["CHR", "POS", "EA", "NEA"], verbose=verbose)
                sumstats.loc[(to_full_fix),status] = vchange_status(sumstats.loc[(to_full_fix),status],3,"975","630") 
                sumstats.loc[(to_part_fix),status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")  
                
            else:
            #when these is no ea or ena, just fix to chr:pos
                to_part_fix = to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull()
                log.log_formula("SNPID", "from CHR POS", source_columns=["CHR", "POS"], verbose=verbose)
                if to_part_fix.sum() > 0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                    sumstats.loc[to_part_fix,status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")
                    
            after_number = sumstats[snpid].isna().sum()
            log.write(" -Fixed "+ str(pre_number - after_number) +" variants ID...", verbose=verbose)
        else:
            log.write(" -ID unfixable: no CHR and POS columns or no SNPID. ", verbose=verbose)

    # Update QC status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_qc_step
            id_kwargs = {
                'snpid': snpid, 'rsid': rsid, 'chrom': chrom, 'pos': pos, 'nea': nea, 'ea': ea, 'status': status,
                'fixprefix': fixprefix, 'fixchrpos': fixchrpos, 'fixid': fixid, 'fixeanea': fixeanea,
                'fixeanea_flip': fixeanea_flip, 'fixsep': fixsep, 'reversea': reversea, 'overwrite': overwrite,
                'forcefixid': forcefixid
            }
            _update_qc_step(sumstats_obj, "id", id_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

@with_logging(
        start_to_msg= "strip SNPID",
        finished_msg= "stripping SNPID",
        start_function= ".strip_snpid()",
        start_cols=["SNPID"]
)
def _strip_SNPID(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], snpid: str = "SNPID", overwrite: bool = False, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:  
    '''
    Strip non-standard characters from SNPID values to standardize format.
    
    Removes leading and trailing non-standard characters from SNPID values that match
    the pattern (xxx:)CHR:POS:ATCG_Allele:ATCG_Allele(:xxx), keeping only the core
    CHR:POS:NEA:EA format.

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    overwrite : bool
        Whether to overwrite existing values.
    verbose : bool, optional
        Whether to print progress.
    Returns
    -------
    pd.DataFrame
        Modified sumstats with stripped SNPIDs.
    '''
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    log.write(" -Checking if SNPID is (xxx:)CHR:POS:ATCG_Allele:ATCG_Allele(:xxx)...(separator: - ,: , _)",verbose=verbose)
    is_chrposrefalt = sumstats[snpid].str.match(SNPID_PATTERN_STRIP, flags=FLAGS, na=False)
    # check if SNPID is NA
    log.write(" -Stripping {} non-NA fixable SNPIDs...".format(sum(is_chrposrefalt)),verbose=verbose)

    # flip 
    extracted = sumstats.loc[is_chrposrefalt, snpid].str.extract(SNPID_PATTERN_STRIP, flags=FLAGS)
    sumstats.loc[is_chrposrefalt, snpid] = (
        extracted["CHR"].astype("string") + ":" +
        extracted["POS"].astype("string") + ":" +
        extracted["NEA"].astype("string") + ":" +
        extracted["EA"].astype("string")
    )

    return sumstats

@with_logging(
        start_to_msg= "flip SNPID from CHR:POS:A1:A2 to CHR:POS:A2:A1",
        finished_msg= "flipping SNPID",
        start_function= ".flip_snpid()",
        start_cols=["SNPID"]
)
def _flip_SNPID(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], snpid: str = "SNPID", overwrite: bool = False, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:  
    '''
    Flip alleles in SNPID values without changing status codes or statistics.
    
    Converts SNPID from CHR:POS:EA:NEA format to CHR:POS:NEA:EA format by swapping
    the last two allele components. This function only modifies the SNPID column and
    does not affect EA, NEA, STATUS, or any statistics.

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    snpid : str
        Column name for SNPID.
    overwrite : bool
        Whether to overwrite existing values.
    verbose : bool, optional
        Whether to print progress.

    Returns
    -------
    pd.DataFrame
        Modified sumstats with flipped alleles.
    '''
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    ##start function with col checking##########################################################
    log.warning("This function only flips alleles in SNPID without changing EA, NEA, STATUS or any statistics.")
    log.write(" -Checking if SNPID is CHR:POS:ATCG_Allele:ATCG_Allele...(separator: - ,: , _)",verbose=verbose)
    is_chrposrefalt = sumstats[snpid].str.match(SNPID_PATTERN_EXTRACT, flags=FLAGS, na=False)
    
    log.write(" -Flipping {} non-NA fixable SNPIDs...".format(sum(is_chrposrefalt)),verbose=verbose)

    # flip 
    parts = sumstats.loc[is_chrposrefalt, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
    parts = parts.astype({"CHR":"string", "POS":"string", "NEA":"string", "EA":"string"})
    sumstats.loc[is_chrposrefalt, snpid] = parts["CHR"].str.cat([parts["POS"], parts["EA"], parts["NEA"]], sep=":")

    return sumstats

###############################################################################################################
# 20230128 
@with_logging(
        start_to_msg= "remove duplicated/multiallelic variants",
        finished_msg= "removing duplicated/multiallelic variants",
        start_function= ".remove_dup()",
        start_cols=None
)
def _remove_dup(sumstats_obj: Union['Sumstats', pd.DataFrame], mode: str = "dm", chrom: str = "CHR", pos: str = "POS", snpid: str = "SNPID", ea: str = "EA", nea: str = "NEA", rsid: str = "rsID", keep: str = 'first', keep_col: str = "P", remove_na: bool = False, keep_ascend: bool = True, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Remove duplicate or multiallelic variants based on user-selected criteria.

    Supports multiple duplicate-identification strategies depending on variant identifiers
    (SNPID, rsID) or allele and coordinate combinations. Can also collapse multi-allelic
    sites by retaining a single representative variant. Variants are sorted by a specified
    column (e.g., P-value) before removal to ensure the best variant is kept.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to process.
    mode : str
        String encoding the deduplication rules; may include one or more of:
        - 'ds' : Identify duplicates using SNPID.
        - 'dr' : Identify duplicates using rsID.
        - 'dc' : Identify duplicates using chromosome, position, effect allele, and non-effect allele.
        - 'm' : Identify multi-allelic variants (same chromosome + position).
    keep : {'first', 'last', False}, default 'first'
        Which record to retain when duplicates are detected.
    keep_col : str, optional
        Column to sort by prior to duplicate removal; used only when `keep` is
        not False.
    remove_na : bool, default False
        If True, remove rows containing missing values in deduplication-relevant columns.
    keep_ascend : bool, default True
        If True, sort in ascending order when determining which duplicate to keep.
    verbose : bool, default False
        If True, print progress and summary information.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with duplicates and multi-allelic variants removed
        according to the specified mode.

    """
    sumstats = sumstats_obj.data

    log.log_operation("Removing mode:{}".format(mode), verbose=verbose)
    if keep_col == "P" and ("P" in sumstats.columns):
        zero_count = (sumstats["P"] == 0).sum()
        if zero_count > 1:
            log.log_operation("Detected {} variants with P=0; converting to MLOG10P and sorting descending".format(zero_count), verbose=verbose)
            status, _ = fill_extreme_mlog10p(sumstats, None, log, verbose=verbose, filled_count=0)
            if status == 1 and ("MLOG10P" in sumstats.columns):
                keep_col = "MLOG10P"
                keep_ascend = False
    # sort the variants using the specified column before removing
    if keep_col is not None : 
        if keep_col in sumstats.columns:
            log.log_operation_start("sort the sumstats using {}".format(keep_col), verbose=verbose)
            sumstats = sumstats.sort_values(by=keep_col,ascending=keep_ascend)
        else:
            log.write("Column" + keep_col +" was not detected... skipping... ", verbose=verbose)
    total_number = len(sumstats)  
    
    # remove by duplicated SNPID
    if (snpid in sumstats.columns) and ("d" in mode or "s" in mode):
        log.log_operation_start("remove duplicated variants based on snpid", version=_get_version(), verbose=verbose)
        check_dataframe_shape(sumstats, log, verbose)
        log.write(" -Which variant to keep: ",  keep , verbose=verbose)   
        pre_number =len(sumstats)   
        if snpid in sumstats.columns:
            # keep na and remove duplicated
            sumstats = sumstats.loc[sumstats[snpid].isna() | (~sumstats.duplicated(subset=[snpid], keep=keep)),:]
            after_number=len(sumstats)   
            log.log_variants_removed(pre_number - after_number, reason="based on SNPID", verbose=verbose)
    
    # remove by duplicated rsID
    if (rsid in sumstats.columns) and ("d" in mode or "r" in mode):
        # keep na and remove duplicated
        pre_number =len(sumstats)
        log.log_operation_start("remove duplicated variants based on rsID", verbose=verbose)
        check_dataframe_shape(sumstats, log, verbose)
        sumstats = sumstats.loc[sumstats[rsid].isna() | (~sumstats.duplicated(subset=rsid, keep=keep)),:]
        after_number=len(sumstats)   
        log.log_variants_removed(pre_number - after_number, reason="based on rsID", verbose=verbose)
    
    # remove by duplicated variants by CHR:POS:NEA:EA
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and (nea in sumstats.columns) and (ea in sumstats.columns) and ("d" in mode or "c" in mode):
        log.log_operation_start("remove duplicated variants based on CHR,POS,EA and NEA", verbose=verbose)
        check_dataframe_shape(sumstats, log, verbose)
        log.write(" -Which variant to keep: ",  keep , verbose=verbose)   
        pre_number =len(sumstats)   
        if snpid in sumstats.columns:
            # keep na and remove duplicated
            sumstats = sumstats.loc[(~sumstats[[chrom,pos,ea,nea]].all(axis=1)) | (~sumstats.duplicated(subset=[chrom,pos,ea,nea], keep=keep)),:]
            after_number=len(sumstats)   
            log.log_variants_removed(pre_number - after_number, reason="based on CHR,POS,EA and NEA", verbose=verbose) 
    
    # remove by multiallelic variants by CHR:POS
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and "m" in mode:
        # keep na and remove duplicated
        pre_number =len(sumstats) 
        log.log_operation_start("remove multiallelic variants based on chr:pos", verbose=verbose)    
        check_dataframe_shape(sumstats, log, verbose)
        log.write(" -Which variant to keep: ",  keep , verbose=verbose) 
        sumstats = sumstats.loc[(~sumstats[[chrom,pos]].all(axis=1)) | (~sumstats.duplicated(subset=[chrom,pos], keep=keep)),:]
        after_number=len(sumstats)  
        log.log_variants_removed(pre_number - after_number, reason="multiallelic variants", verbose=verbose)   
    after_number=len(sumstats)   
    
    # resort the coordinates
    log.log_variants_removed(total_number - after_number, reason="in total", verbose=verbose)
    if keep_col is not None : 
        log.write(" -Sort the coordinates based on CHR and POS...", verbose=verbose)
        sumstats = _sort_coordinate(sumstats,verbose=False)
    
    if "n" in mode or remove_na==True:
        # if remove==True, remove NAs
        log.write(" -Removing NAs...", verbose=verbose)
        pre_number =len(sumstats) 
        specified_columns = []
        if "d" in mode:
            if rsid in sumstats.columns: specified_columns.append(rsid)
            if snpid in sumstats.columns: specified_columns.append(snpid)
            if chrom in sumstats.columns: specified_columns.append(chrom)
            if pos in sumstats.columns: specified_columns.append(pos)
            if ea in sumstats.columns: specified_columns.append(ea)
            if nea in sumstats.columns: specified_columns.append(nea)
        if "r" in mode:
            if rsid in sumstats.columns:specified_columns.append(rsid)
        if "s" in mode:
            if snpid in sumstats.columns:specified_columns.append(snpid)
        if "m" in mode:
            if chrom in sumstats.columns:specified_columns.append(chrom)
            if pos in sumstats.columns:specified_columns.append(pos)
        if "c" in mode:
            if chrom in sumstats.columns:specified_columns.append(chrom)
            if pos in sumstats.columns:specified_columns.append(pos)
            if ea in sumstats.columns:specified_columns.append(ea)
            if nea in sumstats.columns:specified_columns.append(nea)
        specified_columns = list(set(specified_columns))
        sumstats = sumstats.loc[~sumstats[specified_columns].isna().any(axis=1),:]
        after_number=len(sumstats) 
        log.log_variants_removed(pre_number - after_number, reason="with NA values in {}".format(specified_columns), verbose=verbose)  

    # Update the sumstats_obj.data with the filtered dataframe
    sumstats_obj.data = sumstats
    
    # Update QC status
    try:
        from gwaslab.info.g_meta import _update_qc_step
        remove_dup_kwargs = {
            'mode': mode, 'chrom': chrom, 'pos': pos, 'snpid': snpid, 'ea': ea, 'nea': nea, 'rsid': rsid,
            'keep': keep, 'keep_col': keep_col, 'remove_na': remove_na, 'keep_ascend': keep_ascend
        }
        _update_qc_step(sumstats_obj, "remove_dup", remove_dup_kwargs, True)
    except:
        pass

    return sumstats_obj.data

def _build_chrom_list_with_numeric(chrom_list, chromosomes_obj):
    """
    Build a chromosome list that includes both string and numeric representations.
    This is needed after sex chromosomes are converted to numeric values.
    
    Parameters:
    -----------
    chrom_list : list
        Original chromosome list with string identifiers
    chromosomes_obj : Chromosomes or None
        Chromosomes instance to get numeric mappings
        
    Returns:
    --------
    list
        Extended chromosome list with numeric sex chromosome values
    """
    chrom_list_with_numeric = chrom_list.copy()
    
    if chromosomes_obj is not None:
        all_sex_chr_numeric = chromosomes_obj.get_all_sex_chromosomes_numeric()
        # Add numeric values as strings (since sumstats[chrom] is still string at this point)
        for num_val in all_sex_chr_numeric:
            chrom_list_with_numeric.append(str(num_val))
    else:
        # Default numeric values for X, Y, MT
        chrom_list_with_numeric.extend(["23", "24", "25"])
    
    return chrom_list_with_numeric


def _convert_sex_chromosomes_to_numeric(sumstats, fixable_indices, 
                                       x, y, mt, chrom, log, verbose):
    """
    Convert sex chromosome labels (X, Y, MT) to numeric values.
    
    Optimized version that:
    - Uses set operations for faster lookups
    - Pre-computes case variations
    - Directly maps only the subset that needs conversion
    - Checks sumstats directly after values have been assigned
    
    Parameters:
    -----------
    sumstats : pd.DataFrame
        Summary statistics dataframe (already contains extracted chromosome values)
    fixable_indices : pd.Index
        Indices of variants that were fixable (where extracted values were assigned)
    x, y, mt : tuple
        Chromosome mappings as (label, numeric_value)
    chrom : str
        Column name for chromosome
    log : Log
        Logging object
    verbose : bool
        Verbosity flag
    """
    import pandas as pd
    
    # Early return if no fixable indices
    if len(fixable_indices) == 0:
        return
    
    # Pre-compute case variations for faster lookups (use set for O(1) membership testing)
    x_label, x_num = x[0], str(x[1])
    y_label, y_num = y[0], str(y[1])
    mt_label, mt_num = mt[0], str(mt[1])
    
    x_lower, x_upper = x_label.lower(), x_label.upper()
    y_lower, y_upper = y_label.lower(), y_label.upper()
    mt_lower, mt_upper = mt_label.lower(), mt_label.upper()
    
    # Build set of all case variations for fast membership testing
    xymt_set = {x_lower, x_upper, y_lower, y_upper, mt_lower, mt_upper}
    
    # Check which fixable values in sumstats are sex chromosomes (using set for O(1) lookup)
    # Note: sumstats already contains the extracted values at these indices
    fixable_chr_values = sumstats.loc[fixable_indices, chrom].astype(str)
    sex_chr_mask = fixable_chr_values.isin(xymt_set)
    
    # Early return if no sex chromosomes found
    sex_chr_count = sex_chr_mask.sum()
    if sex_chr_count == 0:
        return
    
    log.write(" -Identifying non-autosomal chromosomes : {}, {}, and {} ...".format(x_label, y_label, mt_label), verbose=verbose)
    log.write(" -Identified {} variants on sex chromosomes...".format(sex_chr_count), verbose=verbose)
    
    # Get unique sex chromosome values present in data (as lowercase for comparison)
    sex_chr_values = fixable_chr_values[sex_chr_mask]
    sex_chr_values_lower = set(sex_chr_values.str.lower())
    
    # Build conversion dictionary only for chromosomes actually present in data
    # Since we convert to lowercase before mapping, we only need lowercase keys
    convert_num_to_xymt = {}
    chromosome_mappings = [
        (x_lower, x_num, x_label, "first sex chromosome"),
        (y_lower, y_num, y_label, "second sex chromosome"),
        (mt_lower, mt_num, mt_label, "mitochondrial")
    ]
    
    for label_lower, numeric_str, label_orig, description in chromosome_mappings:
        if label_lower in sex_chr_values_lower:
            convert_num_to_xymt[label_lower] = numeric_str
            log.write(" -Standardizing {} chromosome notations: {} to {}...".format(description, label_orig, numeric_str), verbose=verbose)
    
    # Directly map only the subset that needs conversion
    sex_chr_indices = fixable_indices[sex_chr_mask]
    if len(sex_chr_indices) > 0:
        # Use vectorized string operations for case-insensitive mapping
        # Convert to lowercase, map to numeric values, then assign back
        # All values should map successfully since they're all identified sex chromosomes
        sex_chr_values_lower = sumstats.loc[sex_chr_indices, chrom].astype(str).str.lower()
        sumstats.loc[sex_chr_indices, chrom] = sex_chr_values_lower.map(convert_num_to_xymt)


###############################################################################################################
# 20230128
@with_logging(
        start_to_msg= "fix chromosome notation (CHR)",
        finished_msg= "fixing chromosome notation (CHR)",
        start_function= ".fix_chr()",
        start_cols=["CHR","STATUS"],
        check_dtype=False,
        fix=False
)
def _fix_chr(sumstats_obj: Union['Sumstats', pd.DataFrame], chrom: str = "CHR", status: str = "STATUS", add_prefix: str = "", remove: bool = False, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Standardize chromosome notation and handle special chromosome cases (X, Y, MT).
    
    All chromosome notations are converted to string type first. After fix, all chromosome
    notations will be int. This function normalizes chromosome labels to a consistent format,
    extracts chromosome numbers from various formats (e.g., "chr1", "1", "chrX"), maps special
    chromosomes (X, Y, mitochondrial) to standardized numeric identifiers, and optionally
    removes invalid chromosome values.
    
    Chromosome mappings (x, y, mt), chrom_list, and minchr are automatically derived from
    the Sumstats object's chromosomes attribute (Chromosomes instance).

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to fix.
    add_prefix : str, optional, default=""
        Prefix to prepend to chromosome labels (e.g., "chr").
    remove : bool, default False
        If True, remove records with invalid or unrecognized chromosome labels.
    verbose : bool, default False
        If True, print progress or diagnostic messages.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table with standardized chromosome identifiers.
    """
    import pandas as pd
    
    # ============================================================================
    # Step 1: Initialize and get chromosome parameters
    # ============================================================================
    if isinstance(sumstats_obj, pd.DataFrame):
        sumstats = sumstats_obj
        is_dataframe = True
        chromosomes_obj = None
    else:
        sumstats = sumstats_obj.data
        is_dataframe = False
        chromosomes_obj = getattr(sumstats_obj, 'chromosomes', None)
    
    # Get chromosome mappings from Chromosomes object or use defaults
    if chromosomes_obj is not None:
        x, y, mt = chromosomes_obj.get_chromosome_mappings()
        chrom_list = chromosomes_obj.chromosomes.copy()
        minchr = chromosomes_obj.get_min_chromosome()
    else:
        x, y, mt = ("X", 23), ("Y", 24), ("MT", 25)
        chrom_list = get_chr_list()
        minchr = 1
    
    # ============================================================================
    # Step 2: Convert CHR column to string type
    # ============================================================================
    try:
        log.write(" -Checking CHR data type...", verbose=verbose)
        if sumstats[chrom].dtype != "string":
            sumstats[chrom] = sumstats[chrom].astype("string")
    except:
        log.log_datatype_change("CHR", str(sumstats[chrom].dtype), "string", status="attempt", verbose=verbose)
        sumstats[chrom] = sumstats[chrom].astype("string")
    
    # ============================================================================
    # Step 3: Identify which chromosomes need fixing
    # ============================================================================
    is_chr_fixed = sumstats[chrom].str.isnumeric().fillna(False)
    log.log_variants_with_condition("standardized chromosome notation", is_chr_fixed.sum(), verbose=verbose)
    
    # If all chromosomes are already numeric, skip extraction
    if is_chr_fixed.sum() == len(sumstats):
        log.write(" -All CHR are already fixed...", verbose=verbose)
        sumstats.loc[is_chr_fixed, status] = vchange_status(sumstats.loc[is_chr_fixed, status], 4, "986", "520")
    else:
        # ========================================================================
        # Step 4: Extract and fix chromosome notations
        # ========================================================================
        chr_extracted = sumstats.loc[~is_chr_fixed, chrom].str.extract(CHR_PATTERN_EXTRACT, flags=FLAGS)[1]
        is_chr_fixable = ~chr_extracted.isna()
        log.log_variants_with_condition("fixable chromosome notations", is_chr_fixable.sum(), verbose=verbose)
        
        # Check for NA and invalid chromosomes
        is_chr_na = sumstats.loc[~is_chr_fixed, chrom].isna()
        if is_chr_na.sum() > 0:
            log.log_variants_with_condition("NA chromosome notations", is_chr_na.sum(), verbose=verbose)
        
        is_chr_invalid = (~is_chr_fixable) & (~is_chr_na)
        if is_chr_invalid.sum() > 0:
            log.log_variants_with_condition("invalid chromosome notations", is_chr_invalid.sum(), verbose=verbose)
            try:
                invalid_examples = set(sumstats.loc[~is_chr_fixed, chrom][is_chr_invalid].head())
                log.write(" -A look at invalid chromosome notations:", invalid_examples, verbose=verbose)
            except:
                pass
        else:
            log.write(" -No unrecognized chromosome notations...", verbose=verbose)
        
        # ========================================================================
        # Step 5: Assign extracted values back to sumstats
        # ========================================================================
        # The regex extraction already provides normalized chromosome values (no chr prefix)
        # Just convert to string type and assign directly
        chr_extracted_series = chr_extracted[is_chr_fixable.index].astype("string")
        sumstats.loc[is_chr_fixable.index, chrom] = chr_extracted_series
        
        # ========================================================================
        # Step 6: Convert sex chromosomes to numeric values
        # ========================================================================
        _convert_sex_chromosomes_to_numeric(sumstats, is_chr_fixable.index, x, y, mt, chrom, log, verbose)
        
        # ========================================================================
        # Step 7: Update status codes
        # ========================================================================
        sumstats.loc[is_chr_fixed, status] = vchange_status(sumstats.loc[is_chr_fixed, status], 4, "986", "520")
        if len(is_chr_fixable.index) > 0:
            sumstats.loc[is_chr_fixable.index, status] = vchange_status(
                sumstats.loc[is_chr_fixable.index, status], 4, "986", "520"
            )
            sumstats.loc[is_chr_invalid.index, status] = vchange_status(
                sumstats.loc[is_chr_invalid.index, status], 4, "986", "743"
            )
        
        # ========================================================================
        # Step 8: Remove invalid chromosomes if requested
        # ========================================================================
        if remove:
            chrom_list_with_numeric = _build_chrom_list_with_numeric(chrom_list, chromosomes_obj)
            good_chr = sumstats[chrom].isin(chrom_list_with_numeric)
            unrecognized_num = (~good_chr).sum()
            
            if unrecognized_num > 0:
                try:
                    numeric_chrs = [int(x) for x in chrom_list_with_numeric if x.isnumeric()]
                    if numeric_chrs:
                        log.write(" -Valid CHR list: {} - {}".format(min(numeric_chrs), max(numeric_chrs)), verbose=verbose)
                except:
                    pass
                log.log_variants_removed(unrecognized_num, reason="with chromosome notations not in CHR list", verbose=verbose)
                try:
                    invalid_examples = set(sumstats.loc[~good_chr, chrom].head())
                    log.write(" -A look at chromosome notations not in CHR list:", invalid_examples, verbose=verbose)
                except:
                    pass
                sumstats = sumstats.loc[good_chr, :].copy()

    # ============================================================================
    # Step 9: Convert chromosome column to integer type
    # ============================================================================
    try:
        sumstats[chrom] = sumstats[chrom].astype('Int64')
    except:
        sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    
    # ============================================================================
    # Step 10: Filter out variants with CHR < minchr
    # ============================================================================
    out_of_range_chr = (sumstats[chrom] < minchr).fillna(False)
    out_of_range_num = out_of_range_chr.sum()
    if out_of_range_num > 0:
        log.write(" -Sanity check for CHR...", verbose=verbose)
        log.log_variants_removed(out_of_range_num, reason="with CHR < {}".format(minchr), verbose=verbose)
        sumstats = sumstats.loc[~out_of_range_chr, :]
    
    # ============================================================================
    # Step 11: Update Sumstats object and return
    # ============================================================================
    if not is_dataframe:
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_qc_step
            chr_kwargs = {
                'chrom': chrom, 'status': status, 'add_prefix': add_prefix,
                'remove': remove
            }
            _update_qc_step(sumstats_obj, "chr", chr_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

###############################################################################################################    
# 20230128
@with_logging(
        start_to_msg= "fix basepair positions (POS)",
        finished_msg= "fixing basepair positions (POS)",
        start_function= ".fix_pos()",
        start_cols=["POS","STATUS"],
        check_dtype=False,
        fix=False
)
def _fix_pos(sumstats_obj: Union['Sumstats', pd.DataFrame], pos: str = "POS", status: str = "STATUS", remove: bool = False, verbose: bool = True, lower_limit: int = 0, upper_limit: Optional[int] = None, limit: int = 250000000, log: Log = Log()) -> pd.DataFrame:
    '''
    Standardize and validate genomic base-pair positions.
    
    This function checks that reported genomic positions fall within valid chromosomal bounds
    and optionally removes invalid entries. It handles string-formatted positions with thousands
    separators, converts positions to Int64 type, and filters out positions outside the specified
    range. If explicit limits are not provided, a default maximum bound is applied.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to fix.
    remove : bool, default False
        If True, remove records with invalid or out-of-range positions.
    verbose : bool, default False
        If True, print progress or diagnostic messages.
    lower_limit : int, optional
        Minimum acceptable genomic position. Default is 0.
    upper_limit : int, optional
        Maximum acceptable genomic position.
    limit : int, default 250000000
        Default upper limit applied when `upper_limit` is not provided.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with standardized and validated base-pair positions.
    '''
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called during initialization - no Sumstats object yet
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False
    # Set default upper limit if not provided
    if upper_limit is None:
        upper_limit = limit
    
    # Track initial number of variants for reporting
    all_var_num = len(sumstats)
    
    # Check for missing positions before processing
    is_pos_na = sumstats[pos].isna()
    
    # Handle string/object types: remove thousands separators
    try:
        dtype_str = str(sumstats[pos].dtype)
        if dtype_str == "string" or dtype_str == "object":
            sumstats[pos] = sumstats[pos].astype('string')
            log.write(' -Removing thousands separator "," or underbar "_" ...', verbose=verbose)
            sumstats.loc[~is_pos_na, pos] = sumstats.loc[~is_pos_na, pos].str.replace(r'[,_]', '', regex=True)
    except Exception:
        pass

    # Convert POS to integer type
    try:
        log.log_datatype_change("POS", str(sumstats[pos].dtype), "Int64", status="attempt", verbose=verbose)
        sumstats[pos] = sumstats[pos].astype('Int64')
    except Exception:
        log.log_datatype_change("POS", str(sumstats[pos].dtype), "Int64", status="attempt", verbose=verbose)
        sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    
    # Identify fixed and invalid positions
    # Performance optimization: Compute is_pos_na once after conversion and reuse
    is_pos_na_after = sumstats[pos].isna()
    is_pos_fixed = ~is_pos_na_after
    is_pos_invalid = (~is_pos_na) & (~is_pos_fixed)
    
    # Update status codes for fixed and invalid positions
    sumstats.loc[is_pos_fixed, status] = vchange_status(sumstats.loc[is_pos_fixed, status], 4, "975", "630")
    sumstats.loc[is_pos_invalid, status] = vchange_status(sumstats.loc[is_pos_invalid, status], 4, "975", "842")
    
    # Remove outliers outside the specified bounds
    log.write(" -Position bound:({} , {:,})".format(lower_limit, upper_limit), verbose=verbose)
    # Performance optimization: Reuse is_pos_na_after instead of recomputing, and compute outlier count once
    is_outlier = ((sumstats[pos] <= lower_limit) | (sumstats[pos] >= upper_limit)) & (~is_pos_na_after)
    outlier_num = is_outlier.sum()
    log.log_variants_removed(outlier_num, reason="outliers", verbose=verbose)
    sumstats = sumstats.loc[~is_outlier, :]
    
    # Optionally remove remaining NA positions
    if remove is True:
        # Performance optimization: Recompute is_pos_na only after outlier filtering
        is_pos_na_final = sumstats[pos].isna()
        sumstats = sumstats.loc[~is_pos_na_final, :]
        remain_var_num = len(sumstats)
        log.log_variants_removed(all_var_num - remain_var_num, reason="with bad positions", verbose=verbose)
 
    # Update QC status only if called with Sumstats object
    if not is_dataframe:
        # Assign filtered dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_qc_step
            pos_kwargs = {
                'pos': pos, 'status': status, 'remove': remove, 'lower_limit': lower_limit,
                'upper_limit': upper_limit, 'limit': limit
            }
            _update_qc_step(sumstats_obj, "pos", pos_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

###############################################################################################################    
# 20220514
@with_logging(
        start_to_msg= "fix alleles (EA and NEA)",
        finished_msg= "fixing alleles (EA and NEA)",
        start_function= ".fix_allele()",
        start_cols=["EA","NEA","STATUS"],
        check_dtype=False,
        fix=False
)
def _fix_allele(sumstats_obj: Union['Sumstats', pd.DataFrame], ea: str = "EA", nea: str = "NEA", status: str = "STATUS", remove: bool = False, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    """
    Validate and standardize allele representations.
    
    This function checks allele fields for valid nucleotide characters (A, T, C, G), converts
    all alleles to uppercase, standardizes their format using categorical data types, and
    classifies variants as SNPs, indels, normalized, or not normalized. Optionally, rows with
    invalid allele values can be removed.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to fix.
    remove : bool, default False
        If True, remove variants with invalid allele representations.
    verbose : bool, default False
        If True, print progress or warning messages.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table with validated and standardized allele values.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called during initialization - no Sumstats object yet
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False

    ############################################################################################
    #try:
    #    ea_missing = sum(sumstats[ea].isna())
    #    nea_missing = sum(sumstats[nea].isna())
    #    if sum(ea_missing)>0:
    #        log.write(" -Converting {} missing EA to letter N.".format(ea_missing))
    #        sumstats[ea] = sumstats[ea].add_categories("N").fillna("N")
    #    if sum(sumstats[nea].isna())>0:
    #        log.write(" -Converting {} missing NEA to letter N.".format(nea_missing))
    #        sumstats[nea] = sumstats[nea].add_categories("N").fillna("N")
    #except:
    #    pass

    log.log_operation("Converted all bases to string datatype and UPPERCASE", verbose=verbose)
    # Performance optimization: Compute uppercase once and reuse for categories
    ea_upper = sumstats[ea].str.upper()
    nea_upper = sumstats[nea].str.upper()
    categories = set(ea_upper) | set(nea_upper) | set("N")
    categories = {x for x in categories if pd.notna(x)}
    sumstats[ea] = pd.Categorical(ea_upper, categories=categories) 
    sumstats[nea] = pd.Categorical(nea_upper, categories=categories) 
    all_var_num = len(sumstats)
    
    ## check ATCG
    bad_ea = sumstats[ea].str.contains("[^ATCG]", na=True, regex=True)
    bad_nea = sumstats[nea].str.contains("[^ATCG]", na=True, regex=True)
    good_ea = ~bad_ea
    good_nea = ~bad_nea
    
    # Performance optimization: Compute sums once and reuse
    bad_ea_num = bad_ea.sum()
    bad_nea_num = bad_nea.sum()
    log.log_variants_with_condition("bad EA", bad_ea_num, verbose=verbose)
    log.log_variants_with_condition("bad NEA", bad_nea_num, verbose=verbose)
    
    ## check NA
    is_eanea_na = sumstats[ea].isna() |  sumstats[nea].isna()
    log.log_variants_with_condition("NA for EA or NEA", is_eanea_na.sum(), verbose=verbose)
    
    ## check same alleles
    not_variant = sumstats[nea] == sumstats[ea]
    log.log_variants_with_condition("same EA and NEA", not_variant.sum(), verbose=verbose)

    ## sum up invalid variants
    is_invalid = bad_ea | bad_nea | not_variant
    
    exclude = bad_nea | bad_ea
    # Performance optimization: Compute exclude sum once
    exclude_num = exclude.sum()
    
    if len(set(sumstats.loc[bad_ea,ea].head())) >0:
        log.write(" -A look at the non-ATCG EA:",set(sumstats.loc[bad_ea,ea].head()),"...", verbose=verbose) 
    if len(set(sumstats.loc[bad_nea,nea].head())) >0:
        log.write(" -A look at the non-ATCG NEA:",set(sumstats.loc[bad_nea,nea].head()),"...", verbose=verbose) 
    
    # Compute is_eanea_fixed BEFORE filtering (consistent with previous code)
    is_eanea_fixed = good_ea | good_nea
    
    if remove == True:
        sumstats = sumstats.loc[(good_ea & good_nea),:].copy()
        good_eanea_num = len(sumstats)
        log.log_variants_removed(all_var_num - good_eanea_num, reason="with NA alleles or alleles that contain bases other than A/C/T/G", verbose=verbose)  
        # Filter masks to align with filtered dataframe
        is_eanea_fixed = is_eanea_fixed.loc[sumstats.index]
        is_invalid = is_invalid.loc[sumstats.index]
        is_eanea_na = is_eanea_na.loc[sumstats.index]
        # Recalculate not_variant after filtering to ensure index alignment
        not_variant = sumstats[nea] == sumstats[ea]
        sumstats = sumstats.loc[~not_variant,:].copy()
        good_eanea_notsame_num = len(sumstats)
        log.log_variants_removed(good_eanea_num - good_eanea_notsame_num, reason="with same allele for EA and NEA", verbose=verbose) 
        # Filter masks again after second filtering
        is_eanea_fixed = is_eanea_fixed.loc[sumstats.index]
        is_invalid = is_invalid.loc[sumstats.index]
        is_eanea_na = is_eanea_na.loc[sumstats.index]
    else:
        sumstats[[ea,nea]] = sumstats[[ea,nea]].fillna("N")
        log.write(" -Detected "+str(exclude_num)+" variants with alleles that contain bases other than A/C/T/G .", verbose=verbose) 
    # Performance optimization: Compute uppercase once and reuse for categories
    ea_upper_after = sumstats[ea].str.upper()
    nea_upper_after = sumstats[nea].str.upper()
    categories = set(ea_upper_after) | set(nea_upper_after) | set("N")
    sumstats[ea] = pd.Categorical(ea_upper_after, categories=categories) 
    sumstats[nea] = pd.Categorical(nea_upper_after, categories=categories) 
    
    # Performance optimization: Compute lengths once and reuse
    ea_len = sumstats[ea].str.len()
    nea_len = sumstats[nea].str.len()
    is_snp = (ea_len == 1) & (nea_len == 1)
    is_indel = (ea_len != nea_len)
    is_not_normalized = (ea_len > 1) & (nea_len > 1)
    # Match previous code logic: normalized = indel AND (one length is 1 and other is >1)
    is_normalized = is_indel & ((ea_len == 1) & (nea_len > 1) | (ea_len > 1) & (nea_len == 1))
        
    is_eanea_fixed_not_norm = is_eanea_fixed & is_not_normalized
    is_eanea_fixed_snp = is_eanea_fixed & is_snp
    is_eanea_fixed_indel = is_eanea_fixed & is_indel
    is_eanea_fixed_norm = is_eanea_fixed & is_normalized
    
    # Status updates should always happen (consistent with previous code)
    if sum(is_invalid) > 0:
        sumstats.loc[is_invalid, status] = vchange_status(sumstats.loc[is_invalid, status], 5, "9", "6") 
    if sum(is_eanea_na) > 0:
        sumstats.loc[is_eanea_na, status] = vchange_status(sumstats.loc[is_eanea_na, status], 5, "9", "7")
    if sum(is_eanea_fixed_not_norm) > 0:
        sumstats.loc[is_eanea_fixed_not_norm, status] = vchange_status(sumstats.loc[is_eanea_fixed_not_norm, status], 5, "9", "5")
    if sum(is_eanea_fixed_snp) > 0:
        sumstats.loc[is_eanea_fixed_snp, status] = vchange_status(sumstats.loc[is_eanea_fixed_snp, status], 5, "9", "0")
    if sum(is_eanea_fixed_indel) > 0:
        sumstats.loc[is_eanea_fixed_indel, status] = vchange_status(sumstats.loc[is_eanea_fixed_indel, status], 5, "9", "4")
    if sum(is_eanea_fixed_norm) > 0:
        sumstats.loc[is_eanea_fixed_norm, status] = vchange_status(sumstats.loc[is_eanea_fixed_norm, status], 5, "4", "3")

    # Update QC status only if called with Sumstats object
    if not is_dataframe:
        # Assign filtered dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_qc_step
            allele_kwargs = {
                'ea': ea, 'nea': nea, 'status': status, 'remove': remove
            }
            _update_qc_step(sumstats_obj, "allele", allele_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

###############################################################################################################   
# 20220721

@with_logging(
        start_to_msg= "normalize indels",
        finished_msg= "normalizing indels",
        start_function= ".normalize()",
        start_cols=["EA","NEA","STATUS"]
)
def _parallelize_normalize_allele(sumstats_obj: Union['Sumstats', pd.DataFrame], snpid: str = "SNPID", rsid: str = "rsID", pos: str = "POS", nea: str = "NEA", ea: str = "EA", status: str = "STATUS", chunk: int = 3000000, threads: int = 1, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    '''
    Normalize indels in parallel using left-alignment and parsimony principles.
    
    This function standardizes allele representations for insertion/deletion variants by
    left-aligning and trimming shared sequence context. It removes common suffixes and
    prefixes from both alleles and adjusts positions accordingly, following the VCF
    normalization standard.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to normalize.
    chunk : int, default 3000000
        Size of chunks for parallel processing.
    threads : int, default 1
        Number of threads used for parallel processing.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with normalized indel allele representations.
    '''
    sumstats = sumstats_obj.data
    ############################################################################################

    # digit 5 in [4,5]
    variants_to_check = status_match(sumstats[status], 5, [4, 5])
    # Performance optimization: Compute sum once and reuse
    variants_to_check_num = variants_to_check.sum()
    if variants_to_check_num == 0:
        log.write(" -No available variants to normalize..", verbose=verbose)
        return sumstats_obj.data
    ###############################################################################################################
    # Performance optimization: Extract subset once and reuse
    variants_subset = sumstats.loc[variants_to_check, [pos, nea, ea, status]]
    
    # Use vectorized fast normalization (previously mode="v")
    if variants_to_check_num < 100000:
        threads=1  
    if threads==1:
        normalized_pd, changed_index = fastnormalizeallele(variants_subset, pos=pos, nea=nea, ea=ea, status=status, chunk=chunk, log=log, verbose=verbose)
    else:
        map_func = partial(fastnormalizeallele, pos=pos, nea=nea, ea=ea, status=status)
        df_split = _df_split(variants_subset, threads)
        with Pool(threads) as pool:
            results = pool.map(map_func,df_split)
        normalized_pd = pd.concat([i[0] for i in results])
        changed_index = np.concatenate([i[1] for i in results])
        del results
        gc.collect()
    ###############################################################################################################
    try:
        # Performance optimization: changed_index may not align with sumstats index, use it directly
        changed_num = len(changed_index)
        if changed_num>0:
            # Get example variants from the changed_index (which are indices from variants_subset)
            # Map back to original sumstats indices
            if hasattr(changed_index, '__len__') and len(changed_index) > 0:
                variant_indices = variants_subset.index[changed_index]
                if len(variant_indices) > 0:
                    example_sumstats = sumstats.loc[variant_indices[:5], :]  # Get first 5 examples
                    if snpid in example_sumstats.columns:
                        before_normalize_id = example_sumstats[snpid]
                    elif rsid in example_sumstats.columns:
                        before_normalize_id = example_sumstats[rsid]
                    else:
                        before_normalize_id = example_sumstats.index
                    
                    log.write(" -Not normalized allele IDs:",end="", verbose=verbose)
                    for i in before_normalize_id.values:
                        log.write(i,end=" ",show_time=False)
                    log.write("... \n",end="",show_time=False, verbose=verbose) 
                
                    log.write(" -Not normalized allele:",end="", verbose=verbose)
                    for i in example_sumstats[[ea,nea]].values:
                        log.write(i,end="",show_time=False, verbose=verbose)
                    log.write("... \n",end="",show_time=False, verbose=verbose)     
            log.write(" -Modified "+str(changed_num) +" variants according to parsimony and left alignment principal.", verbose=verbose)
        else:
            log.write(" -All variants are already normalized..", verbose=verbose)
    except:
        pass
    
    # Performance optimization: Compute categories more efficiently by combining sets in one operation
    categories = set(sumstats[ea]) | set(sumstats[nea]) | set(normalized_pd[ea]) | set(normalized_pd[nea])
    sumstats[ea] = pd.Categorical(sumstats[ea], categories=categories) 
    sumstats[nea] = pd.Categorical(sumstats[nea], categories=categories) 
    sumstats.loc[variants_to_check,[pos,nea,ea,status]] = normalized_pd.values
    
    # Optimize POS conversion: check dtype first, then convert only if needed
    if sumstats[pos].dtype.name != 'Int64':
        try:
            sumstats[pos] = sumstats[pos].astype('Int64')
        except (ValueError, TypeError):
            sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
  
    # Assign modified dataframe back to the Sumstats object
    sumstats_obj.data = sumstats
  
    # Update QC status
    try:
        from gwaslab.info.g_meta import _update_qc_step
        # Extract relevant kwargs for QC status tracking (exclude internal variables)
        normalize_kwargs = {
            'snpid': snpid,
            'rsid': rsid,
            'pos': pos,
            'nea': nea,
            'ea': ea,
            'status': status,
            'chunk': chunk,
            'threads': threads
        }
        _update_qc_step(sumstats_obj, "normalize", normalize_kwargs, True)
    except:
        pass
  
    return sumstats_obj.data

def fastnormalizeallele(insumstats: pd.DataFrame, pos: str = "POS", nea: str = "NEA", ea: str = "EA", status: str = "STATUS", chunk: int = 3000000, log: Log = Log(), verbose: bool = False) -> Tuple[pd.DataFrame, np.ndarray]:
    log.write(" -Number of variants to check:{}".format(len(insumstats)), verbose=verbose)
    log.write(" -Chunk size:{}".format(chunk), verbose=verbose)
    log.write(" -Processing in chunks:",end="", verbose=verbose)
    # Performance optimization: Convert to string once before loop instead of every iteration
    insumstats[nea] = insumstats[nea].astype("string")
    insumstats[ea] = insumstats[ea].astype("string")
    # Performance optimization: Use list to collect indices, then concatenate once (O(n) instead of O(n²))
    changed_index_list = []
    num_chunks = len(insumstats) // chunk + 1
    for part_n in range(num_chunks):
        log.write(part_n, end=" ",show_time=False, verbose=verbose)
        chunk_slice = slice(part_n * chunk, (part_n + 1) * chunk)
        insumstats.iloc[chunk_slice, :], changed_index_single = normalizae_chunk(
            insumstats.iloc[chunk_slice, :].copy(), pos=pos, nea=nea, ea=ea, status=status
        )
        changed_index_list.append(changed_index_single)
    # Performance optimization: Concatenate once at the end
    changed_index = np.concatenate(changed_index_list) if changed_index_list else np.array([], dtype=int)
    log.write("\n",end="",show_time=False, verbose=verbose)   
    return insumstats, changed_index

def normalizae_chunk(sumstats: pd.DataFrame, pos: str = "POS", nea: str = "NEA", ea: str = "EA", status: str = "STATUS") -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Normalize indels in a chunk by removing common suffixes and prefixes.
    Optimized version using vectorized pandas operations.
    """
    # Early exit if empty
    if len(sumstats) == 0:
        return sumstats, np.array([], dtype=int)
    
    # Pre-compute boolean masks and lengths
    is_same = sumstats[nea] == sumstats[ea]
    nea_len = sumstats[nea].str.len()
    ea_len = sumstats[ea].str.len()
    is_normalized = ((nea_len == 1) | (ea_len == 1)) & (~is_same)
    
    # Early exit if all are already normalized
    if is_normalized.all():
        sumstats.loc[is_normalized, status] = vchange_status(sumstats.loc[is_normalized, status], 5, "4", "0")
        sumstats.loc[is_same, status] = vchange_status(sumstats.loc[is_same, status], 5, "4", "3")
        return sumstats, np.array([], dtype=int)
    
    # Track changed variants
    changed = pd.Series(False, index=sumstats.index, dtype=bool)
    
    # Right side - remove common suffixes
    max_length = max(nea_len.max(), ea_len.max()) if len(sumstats) > 0 else 0
    
    for i in range(1, max_length):
        # Use vectorized string operations - more efficient than individual string access
        is_pop = (sumstats[nea].str[-1] == sumstats[ea].str[-1]) & (~is_normalized)
        is_pop_num = is_pop.sum()
        if is_pop_num == 0:
            break
        
        if i == 1:
            changed |= is_pop
        
        # Update lengths incrementally
        nea_len[is_pop] = nea_len[is_pop] - 1
        ea_len[is_pop] = ea_len[is_pop] - 1
        
        # Vectorized string slicing - pandas handles this efficiently
        sumstats.loc[is_pop, nea] = sumstats.loc[is_pop, nea].str[:-1]
        sumstats.loc[is_pop, ea] = sumstats.loc[is_pop, ea].str[:-1]
        
        # Update is_normalized mask using cached lengths
        is_normalized = ((nea_len == 1) | (ea_len == 1)) & (~is_same)
    
    # Left side - remove common prefixes and adjust position
    max_length = max(nea_len.max(), ea_len.max()) if len(sumstats) > 0 else 0
    
    for i in range(1, max_length):
        # Use vectorized string operations
        is_pop = (sumstats[nea].str[0] == sumstats[ea].str[0]) & (~is_normalized)
        is_pop_num = is_pop.sum()
        if is_pop_num == 0:
            break
        
        if i == 1:
            changed |= is_pop
        
        # Vectorized string slicing and position update
        sumstats.loc[is_pop, nea] = sumstats.loc[is_pop, nea].str[1:]
        sumstats.loc[is_pop, ea] = sumstats.loc[is_pop, ea].str[1:]
        sumstats.loc[is_pop, pos] = sumstats.loc[is_pop, pos] + 1
        
        # Update lengths incrementally
        nea_len[is_pop] = nea_len[is_pop] - 1
        ea_len[is_pop] = ea_len[is_pop] - 1
        
        # Update is_normalized mask using cached lengths
        is_normalized = ((nea_len == 1) | (ea_len == 1)) & (~is_same)
    
    # Update status
    # Update variants that were previously classified as indel (digit 5 = 4) to SNP (0) or normalized indel (3)
    sumstats.loc[is_normalized, status] = vchange_status(sumstats.loc[is_normalized, status], 5, "4", "0")
    sumstats.loc[is_same, status] = vchange_status(sumstats.loc[is_same, status], 5, "4", "3")
    # Also update variants that were previously classified as not normalized (digit 5 = 5) to normalized indel (3)
    # This handles cases where variants like AT/GT get normalized to A/G
    # Note: These variants get status 3 (normalized indel) even if they become SNPs after normalization,
    # because they were originally not normalized and got normalized
    sumstats.loc[is_normalized, status] = vchange_status(sumstats.loc[is_normalized, status], 5, "5", "3")
    
    # Return changed indices
    changed_index = sumstats.index[changed]
    return sumstats, changed_index.values

###############################################################################################################
# 20220426
def get_reverse_complementary_allele(a: str) -> str:
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a[::-1].translate(dic)

def flip_direction(series: pd.Series) -> pd.Series:
    """
    Flip direction string by swapping '+' and '-' characters (vectorized).
    
    Converts '+' to '-' and '-' to '+', while preserving '?' and other characters.
    This is used to flip the direction of effect when alleles are swapped.
    
    Parameters
    ----------
    series : pd.Series
        Series of direction strings to flip (e.g., "++-?", "+-+", etc.)
        
    Returns
    -------
    pd.Series
        Series with flipped direction strings (e.g., "--+?", "-+-", etc.)
    """
    # Use vectorized string operations with translate for maximum performance
    # Create translation table once and apply to all values
    translation_table = str.maketrans("+-", "-+")
    return series.astype(str).apply(lambda x: x.translate(translation_table) if pd.notna(x) else x)

def flip_by_swap(sumstats: pd.DataFrame, matched_index: pd.Series, log: Log, verbose: bool) -> pd.DataFrame:
    """
    Swap NEA and EA columns for matched variants.
    """
    # Early exit if no matches
    if not matched_index.any():
        return sumstats
    
    if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns):
        log.write(" -Swapping column: NEA <=> EA...", verbose=verbose) 
        # Efficient swap using values to avoid index alignment issues
        sumstats.loc[matched_index, ['NEA', 'EA']] = sumstats.loc[matched_index, ['EA', 'NEA']].values
    return sumstats

def flip_by_inverse(sumstats: pd.DataFrame, matched_index: pd.Series, log: Log, verbose: bool, cols: Optional[List[str]] = None, factor: float = 1) -> pd.DataFrame:
    """
    Flip ratio statistics (OR, HR) by taking inverse (1/x) for matched variants.
    """
    # Early exit if no matches
    if not matched_index.any():
        return sumstats
    
    # Cache column checks
    has_or = "OR" in sumstats.columns
    has_or_95l = "OR_95L" in sumstats.columns
    has_or_95u = "OR_95U" in sumstats.columns
    has_hr = "HR" in sumstats.columns
    has_hr_95l = "HR_95L" in sumstats.columns
    has_hr_95u = "HR_95U" in sumstats.columns
    
    if has_or:
        log.write(" -Flipping column: OR = 1 / OR...", verbose=verbose) 
        sumstats.loc[matched_index, "OR"] = factor / sumstats.loc[matched_index, "OR"].values
    
    if has_or_95l:
        log.write(" -Flipping column: OR_95U = 1 / OR_95L...", verbose=verbose) 
        sumstats.loc[matched_index, "OR_95U"] = factor / sumstats.loc[matched_index, "OR_95L"].values
    
    if has_or_95u:
        log.write(" -Flipping column: OR_95L = 1 / OR_95U...", verbose=verbose) 
        sumstats.loc[matched_index, "OR_95L"] = factor / sumstats.loc[matched_index, "OR_95U"].values
    
    if has_hr:
        log.write(" -Flipping column: HR = 1 / HR...", verbose=verbose) 
        sumstats.loc[matched_index, "HR"] = factor / sumstats.loc[matched_index, "HR"].values
    
    if has_hr_95l:
        log.write(" -Flipping column: HR_95U = 1 / HR_95L...", verbose=verbose) 
        sumstats.loc[matched_index, "HR_95U"] = factor / sumstats.loc[matched_index, "HR_95L"].values
    
    if has_hr_95u:
        log.write(" -Flipping column: HR_95L = 1 / HR_95U...", verbose=verbose) 
        sumstats.loc[matched_index, "HR_95L"] = factor / sumstats.loc[matched_index, "HR_95U"].values
    
    return sumstats

def flip_by_subtract(sumstats: pd.DataFrame, matched_index: pd.Series, log: Log, verbose: bool, cols: Optional[List[str]] = None, factor: float = 1) -> pd.DataFrame:
    """
    Flip frequency statistics (EAF) by subtracting from factor (1 - EAF) for matched variants.
    """
    # Early exit if no matches
    if not matched_index.any():
        return sumstats
    
    if "EAF" in sumstats.columns:
        log.write(" -Flipping column: EAF = 1 - EAF...", verbose=verbose) 
        sumstats.loc[matched_index, "EAF"] = factor - sumstats.loc[matched_index, "EAF"].values
    return sumstats

def flip_by_sign(sumstats: pd.DataFrame, matched_index: pd.Series, log: Log, verbose: bool, cols: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Flip sign of effect size statistics (BETA, Z, T) and swap confidence intervals for matched variants.
    """
    # Early exit if no matches
    if not matched_index.any():
        return sumstats
    
    # Cache column checks
    has_beta = "BETA" in sumstats.columns
    has_beta_95l = "BETA_95L" in sumstats.columns
    has_beta_95u = "BETA_95U" in sumstats.columns
    has_z = "Z" in sumstats.columns
    has_t = "T" in sumstats.columns
    has_direction = "DIRECTION" in sumstats.columns
    
    if has_beta:
        log.write(" -Flipping column: BETA = - BETA...", verbose=verbose) 
        sumstats.loc[matched_index, "BETA"] = -sumstats.loc[matched_index, "BETA"].values
    
    if has_beta_95l:
        log.write(" -Flipping column: BETA_95U = - BETA_95L...", verbose=verbose) 
        sumstats.loc[matched_index, "BETA_95U"] = -sumstats.loc[matched_index, "BETA_95L"].values
    
    if has_beta_95u:
        log.write(" -Flipping column: BETA_95L = - BETA_95U...", verbose=verbose) 
        sumstats.loc[matched_index, "BETA_95L"] = -sumstats.loc[matched_index, "BETA_95U"].values
    
    if has_z:
        log.write(" -Flipping column: Z = - Z...", verbose=verbose) 
        sumstats.loc[matched_index, "Z"] = -sumstats.loc[matched_index, "Z"].values
    
    if has_t:
        log.write(" -Flipping column: T = - T...", verbose=verbose) 
        # Bug fix: was assigning to "Z" instead of "T"
        sumstats.loc[matched_index, "T"] = -sumstats.loc[matched_index, "T"].values
    
    if has_direction:
        log.write(" -Flipping column: DIRECTION +-?0 <=> -+?0 ...", verbose=verbose) 
        # Use vectorized flip_direction function
        sumstats.loc[matched_index, "DIRECTION"] = flip_direction(sumstats.loc[matched_index, "DIRECTION"])
    
    return sumstats

@with_logging(
        start_to_msg= "adjust statistics based on STATUS code",
        finished_msg= "adjusting statistics based on STATUS code",
        start_function= ".flip_allele_stats()",
        start_cols=None
)
def _flip_allele_stats(sumstats_obj: Union['Sumstats', pd.DataFrame], status: str = "STATUS", verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    '''
    Adjust statistics when allele direction has changed based on STATUS codes.
    
    This function adjusts effect sizes and allele-specific statistics when variants have been
    flipped or converted to reverse complement. It handles multiple scenarios: reverse
    complement conversion for SNPs, allele swapping for REF/ALT mismatches, flipping for
    standardized indels, and strand flipping for palindromic variants. Run after checking
    with reference sequence.

    Parameters
    ----------
    sumstats_obj : Sumstats or pandas.DataFrame
        Sumstats object or DataFrame containing the data to process.
    verbose : bool, default False
        If True, print progress messages during processing.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with effect sizes and alleles flipped where required.
    '''
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Cache column checks for performance
    has_nea_ea = ("NEA" in sumstats.columns) and ("EA" in sumstats.columns)
    
    if_stats_flipped = False
    
    ###################get reverse complementary####################
    # Use status_match for integer status codes (digit 6, values 4 or 5)
    matched_index = status_match(sumstats[status], 6, [4, 5])
    matched_count = matched_index.sum()
    if matched_count > 0:
        log.log_operation_start("convert alleles to reverse complement for SNPs with status xxxxx[45]x", version=_get_version(), verbose=verbose) 
        log.write(f" -Flipping {matched_count:,} variants...", verbose=verbose) 
        if has_nea_ea:
            log.log_operation("Converting to reverse complement: EA and NEA", verbose=verbose) 
            # Vectorized reverse complement conversion
            matched_subset = sumstats.loc[matched_index, ['NEA', 'EA']]
            reverse_complement_nea = matched_subset['NEA'].apply(get_reverse_complementary_allele)
            reverse_complement_ea = matched_subset['EA'].apply(get_reverse_complementary_allele)
            
            # Convert to strings once to avoid categorical dtype conflicts when assigning
            reverse_complement_nea = reverse_complement_nea.astype(str)
            reverse_complement_ea = reverse_complement_ea.astype(str)
            
            # Update categories to include reverse complement values
            # Convert full columns to strings once (reused for both unique extraction and categorical conversion)
            ea_str = sumstats['EA'].astype(str)
            nea_str = sumstats['NEA'].astype(str)
            
            # Get unique values efficiently using .unique() (faster than set conversion from Series)
            # Combine all unique values and filter invalid ones
            categories = set(ea_str.unique()) | set(nea_str.unique()) | set(reverse_complement_nea.unique()) | set(reverse_complement_ea.unique())
            categories = {c for c in categories if pd.notna(c) and c not in ('nan', '<NA>', 'None')}
            
            # Convert to categorical with updated categories (reusing already-converted strings)
            sumstats['EA'] = pd.Categorical(ea_str, categories=categories)
            sumstats['NEA'] = pd.Categorical(nea_str, categories=categories)
            
            # Assign the reverse complement values (already strings, pandas handles categorical conversion)
            sumstats.loc[matched_index, 'NEA'] = reverse_complement_nea
            sumstats.loc[matched_index, 'EA'] = reverse_complement_ea
            sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 6, "4", "2")
            log.write(" -Changed the status for flipped variants : xxxxx4x -> xxxxx2x", verbose=verbose)
        if_stats_flipped = True
    
    ###################flip ref####################
    # digit 6 in [3,5]
    matched_index = status_match(sumstats[status], 6, [3, 5])
    matched_count = matched_index.sum()
    if matched_count > 0:
        log.write(f"Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...{_get_version()}", verbose=verbose) 
        log.write(f" -Flipping {matched_count:,} variants...", verbose=verbose) 
        
        flip_by_swap(sumstats, matched_index, log, verbose)
        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants : xxxxx[35]x -> xxxxx[12]x", verbose=verbose) 
        sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 6, "35", "12")
        if_stats_flipped = True
        
    ###################flip ref for undistingushable indels####################
    # Pattern: xxxx[123][67]6 means digit 5 in [1,2,3] AND digit 6 in [6,7] AND digit 7 = 6
    # Use status_match to ensure all conditions are met (more reliable than regex pattern)
    digit_5_match = status_match(sumstats[status], 5, [1, 2, 3])
    digit_6_match = status_match(sumstats[status], 6, [6, 7])
    digit_7_match = status_match(sumstats[status], 7, [6])
    matched_index = digit_5_match & digit_6_match & digit_7_match
    matched_count = matched_index.sum()
    if matched_count > 0:
        log.write(f"Start to flip allele-specific stats for standardized indels with status xxxx[123][67][6]: ALT->EA , REF->NEA...{_get_version()}", verbose=verbose) 
        log.write(f" -Flipping {matched_count:,} variants...", verbose=verbose) 
        
        flip_by_swap(sumstats, matched_index, log, verbose)
        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4", verbose=verbose) 
        sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 7, "6", "4")
        if_stats_flipped = True
    
    ###################flip statistics for reverse strand panlindromic variants####################
    # Pattern: xxxxx[012]5 means digit 6 in [0,1,2] AND digit 7 = 5
    # Use status_match to ensure both conditions are met (more reliable than regex pattern with literal digits)
    digit_6_match = status_match(sumstats[status], 6, [0, 1, 2])
    digit_7_match = status_match(sumstats[status], 7, [5])
    matched_index = digit_6_match & digit_7_match
    matched_count = matched_index.sum()
    if matched_count > 0:
        log.write(f"Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...{_get_version()}", verbose=verbose) 
        log.write(f" -Flipping {matched_count:,} variants...", verbose=verbose) 

        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2", verbose=verbose) 
        sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 7, "5", "2")
        if_stats_flipped = True

    if not if_stats_flipped:
        log.write(" -No statistics have been changed.", verbose=verbose)
    
    # Update harmonization status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            flip_kwargs = {'status': status}
            _update_harmonize_step(sumstats_obj, "flip_allele_stats", flip_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats


###############################################################################################################
# 20220426
@with_logging(
        start_to_msg= "sort the genome coordinates",
        finished_msg= "sorting coordinates",
        start_function= ".sort_coordinate()",
        start_cols=["CHR","POS"],
        show_shape=False
)
def _sort_coordinate(sumstats_obj: Union['Sumstats', pd.DataFrame], chrom: str = "CHR", pos: str = "POS", reindex: bool = True, verbose: bool = True, log: Log = Log()) -> pd.DataFrame:
    '''
    Sort variants by genomic coordinates (chromosome, then position).
    
    Sorts the dataframe first by chromosome number, then by position in ascending order.
    The index is reset to sequential integers after sorting.

    Parameters
    ----------
    sumstats_obj : Sumstats or pandas.DataFrame
        Sumstats object or DataFrame containing the data to sort.
    verbose : bool, default False
        If True, print progress messages.

    Returns
    -------
    pandas.DataFrame
        DataFrame with sorted genomic coordinates.
    '''
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called with DataFrame
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    # Performance optimization: Check if POS needs conversion only if not already Int64
    if sumstats[pos].dtype != "Int64":
        try:
            log.log_datatype_change("POS", str(sumstats[pos].dtype), "Int64", status="attempt", verbose=verbose)
            sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
        except Exception:
            pass
    
    # Performance optimization: Use reindex parameter to control index reset
    # If reindex=False, preserve original index; if True, reset to sequential integers
    sumstats = sumstats.sort_values(by=[chrom, pos], ascending=True, ignore_index=reindex)
    
    # Update QC status and metadata only if called with Sumstats object
    if not is_dataframe:
        # Assign sorted dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_qc_step
            sort_coord_kwargs = {'chrom': chrom, 'pos': pos, 'reindex': reindex}
            _update_qc_step(sumstats_obj, "sort_coord", sort_coord_kwargs, True)
            # Set metadata
            sumstats_obj.meta["is_sorted"] = True
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats
###############################################################################################################
# 20230430 added HR HR_95 BETA_95 N_CASE N_CONTROL
@with_logging(
        start_to_msg= "reorder the columns",
        finished_msg= "reordering the columns",
        start_function= ".sort_column()",
        start_cols=None,
        show_shape=False
)
def _sort_column(sumstats_obj: Union['Sumstats', pd.DataFrame], verbose: bool = True, log: Log = Log(), order: Optional[List[str]] = None) -> pd.DataFrame:
    '''
    Reorder columns according to a specified order.
    
    Reorders the dataframe columns to match a predefined standard order, placing standard
    GWAS columns first (SNPID, rsID, CHR, POS, EA, NEA, statistics, etc.) followed by
    any additional columns not in the standard list.

    Parameters
    ----------
    sumstats_obj : Sumstats or pd.DataFrame
        Sumstats object or DataFrame containing the data to reorder.
    verbose : bool, optional
        Whether to print progress. Default is True.

    Returns
    -------
    pd.DataFrame
        Modified sumstats with reordered columns.
    '''
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called during initialization - no Sumstats object yet
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False

    if order is None:
        order = DEFAULT_COLUMN_ORDER.copy()
    # Performance optimization: Create set from order for O(1) lookup
    order_set = set(order)
    
    # Performance optimization: Convert to sets for O(1) lookup instead of O(n) list membership
    sumstats_columns_set = set(sumstats.columns)
    
    # Performance optimization: Use list comprehension with set lookup for better performance
    # First, get columns in order that exist in sumstats
    output_columns = [col for col in order if col in sumstats_columns_set]
    # Then, append remaining columns not in order
    output_columns.extend([col for col in sumstats.columns if col not in order_set])
    
    log.write(" -Reordering columns to    :", ",".join(output_columns), verbose=verbose)
    sumstats = sumstats[output_columns]

    # Update QC status and metadata only if called with Sumstats object
    if not is_dataframe:
        # Assign reordered dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_qc_step
            sort_column_kwargs = {'order': order}
            _update_qc_step(sumstats_obj, "sort_column", sort_column_kwargs, True)
            # Set metadata
            sumstats_obj.meta["is_sorted"] = True
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats


###############################################################################################################
def _df_split(dataframe: pd.DataFrame, n: int) -> List[pd.DataFrame]:
    k, m = divmod(len(dataframe), n)
    return [dataframe.iloc[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]
