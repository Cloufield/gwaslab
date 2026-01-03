import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from functools import partial
from gwaslab.info.g_vchange_status_polars import vchange_statusp
from gwaslab.info.g_vchange_status import status_match
from gwaslab.info.g_vchange_status import change_status
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chain

from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_build import _set_build
from gwaslab.qc.qc_pattern import CHR_PATTERN_EXTRACT, FLAGS

from gwaslab.util.util_in_fill_data import _convert_betase_to_mlog10p
from gwaslab.util.util_in_fill_data import _convert_betase_to_p
from gwaslab.util.util_in_fill_data import _convert_mlog10p_to_p

import polars as pl
from typing import TYPE_CHECKING, Union, Optional

if TYPE_CHECKING:
    from gwaslab.g_Sumstats_polars import Sumstatsp
###############################################################################################################
# 20220426
def get_reverse_complementary_allele(a):
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a[::-1].translate(dic)

def flip_direction(string):
    flipped_string=""
    for char in string:
        if char=="?":
            flipped_string+="?"
        elif char=="+":
            flipped_string+="-"
        elif char=="-":
            flipped_string+="+"
        else: #sometime it is 0
            flipped_string+=char
    return flipped_string

def flip_by_swap(sumstats, matched_index, log, verbose):
    if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
        log.write(" -Swapping column: NEA <=> EA...", verbose=verbose)

        sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  pl.col("EA")  )  
                .otherwise( pl.col("NEA") )
                .alias("NEA"),

                pl.when( matched_index )  
                .then(   pl.col("NEA")  )  
                .otherwise( pl.col("EA") )
                .alias("EA"),
                )

    return sumstats

def flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1):
    for header in ["OR","OR_95L","OR_95U","HR","HR_95L","HR_95U"]:
        if header in sumstats.columns:
                log.write(" -Flipping column: {header} = 1 / {header}...".format(header = header), verbose=verbose) 
                sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  1/ pl.col(header) )  
                .otherwise( pl.col(header) )
                .alias(header)
                )
    return sumstats

def flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1):
    header="EAF"
    if header in sumstats.columns:
        log.write(" -Flipping column: EAF = 1 - EAF...", verbose=verbose) 
        sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  1 - pl.col(header) )  
                .otherwise( pl.col(header) )
                .alias(header)
                )
    return sumstats

def flip_by_sign(sumstats, matched_index, log, verbose, cols=None):
    for header in ["BETA","BETA_95L","BETA_95U","T","Z"]:
        if header in sumstats.columns:
                log.write(" -Flipping column: {header} = - {header}...".format(header = header), verbose=verbose) 
                sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  - pl.col(header) )  
                .otherwise( pl.col(header) )
                .alias(header)
                )
    
    if "DIRECTION" in sumstats.columns:
        sumstats = sumstats.with_columns(
                pl.when( matched_index )  
                .then(  pl.col("DIRECTION").map_batches(lambda x: pl.Series(flip_direction(x))) )  
                .otherwise( pl.col("DIRECTION") )
                .alias("DIRECTION")
                )
    return sumstats

def flipallelestatsp(sumstats,status="STATUS",verbose=True,log=Log()):
    ##start function with col checking#########################################################

    if_stats_flipped = False
    ###################get reverse complementary####################
    pattern = r"\w\w\w\w\w[45]\w"  
    #matched_index = status_match(sumstats[status],6,[4,5]) #
    #matched_index = sumstats[status].str[5].str.match(r"4|5")

    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w\w[45]\w")

    if len(sumstats.filter(matched_index))>0:
        log.write("Start to convert alleles to reverse complement for SNPs with status xxxxx[45]x...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            log.write(" -Converting to reverse complement : EA and NEA...", verbose=verbose) 
            
            sumstats = sumstats.filter(matched_index).with_columns(
                NEA = pl.col("NEA").map_batches(lambda x: pl.Series(get_reverse_complementary_allele(x))),
                EA = pl.col("EA").map_batches(lambda x: pl.Series(get_reverse_complementary_allele(x)))
                                                                   )
        
            sumstats  = vchange_statusp(sumstats, matched_index, status,6, ["4"], ["2"])
            log.write(" -Changed the status for flipped variants : xxxxx4x -> xxxxx2x", verbose=verbose)
        if_stats_flipped = True
    
    ###################flip ref####################
    pattern = r"\w\w\w\w\w[35]\w"  
    #matched_index = status_match(sumstats[status],6,[3,5]) #sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w\w[35]\w")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        
        sumstats = flip_by_swap(sumstats, matched_index, log, verbose)
        sumstats = flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        sumstats = flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        sumstats = flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants : xxxxx[35]x -> xxxxx[12]x", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status,6, ["35"], ["12"])
        if_stats_flipped = True
        
    ###################flip ref for undistingushable indels####################
    pattern = r"\w\w\w\w[123][67]6"  
    #matched_index = status_match(sumstats[status],6,[1,2,3])|status_match(sumstats[status],6,[6,7])|status_match(sumstats[status],7,6) #sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w[123][67]6")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for standardized indels with status xxxx[123][67][6]: ALT->EA , REF->NEA...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 
        
        sumstats = flip_by_swap(sumstats, matched_index, log, verbose)
        sumstats = flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        sumstats = flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        sumstats = flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status, 7, ["6"], ["4"])
        if_stats_flipped = True
         # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    pattern = r"\w\w\w\w\w[012]5"  
    #matched_index = status_match(sumstats[status],6,[0,1,2]) | status_match(sumstats[status],7,[5])#sumstats[status].str.match(pattern)
    matched_index = pl.col(status).cast(pl.String).str.contains("^\w\w\w\w\w[012]5")
    if len(sumstats.filter(matched_index))>0:
        log.write("Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(len(sumstats.filter(matched_index))) +" variants...", verbose=verbose) 

        sumstats = flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        sumstats = flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        sumstats = flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2", verbose=verbose) 
        sumstats  = vchange_statusp(sumstats, matched_index,status,7, ["5"], ["2"])
        if_stats_flipped = True

    if if_stats_flipped != True:
        log.write(" -No statistics have been changed.")
    return sumstats

###############################################################################################################
# 20230128
def _fix_chrp(sumstats_obj: Union['Sumstatsp', pl.DataFrame], chrom: str = "CHR", status: str = "STATUS", add_prefix: str = "", remove: bool = False, verbose: bool = True, log: Log = Log()) -> pl.DataFrame:
    """
    Standardize chromosome notation and handle special chromosome cases (X, Y, MT) using polars.
    
    All chromosome notations are converted to string type first. After fix, all chromosome
    notations will be int. This function normalizes chromosome labels to a consistent format,
    extracts chromosome numbers from various formats (e.g., "chr1", "1", "chrX"), maps special
    chromosomes (X, Y, mitochondrial) to standardized numeric identifiers, and optionally
    removes invalid chromosome values.
    
    Chromosome mappings (x, y, mt), chrom_list, and minchr are automatically derived from
    the Sumstats object's chromosomes attribute (Chromosomes instance).

    Parameters
    ----------
    sumstats_obj : Sumstatsp or pl.DataFrame
        Sumstats object or polars DataFrame containing the data to fix.
    chrom : str, default "CHR"
        Column name for chromosome.
    status : str, default "STATUS"
        Column name for status.
    add_prefix : str, optional, default=""
        Prefix to prepend to chromosome labels (e.g., "chr").
    remove : bool, default False
        If True, remove records with invalid or unrecognized chromosome labels.
    verbose : bool, default False
        If True, print progress or diagnostic messages.
    log : Log, optional
        Logging object.

    Returns
    -------
    pl.DataFrame
        Polars DataFrame with standardized chromosome identifiers.
    """
    # ============================================================================
    # Step 1: Initialize and get chromosome parameters
    # ============================================================================
    if isinstance(sumstats_obj, pl.DataFrame):
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
        if sumstats[chrom].dtype != pl.String:
            sumstats = sumstats.with_columns(pl.col(chrom).cast(pl.String).alias(chrom))
    except:
        log.log_datatype_change("CHR", str(sumstats[chrom].dtype), "string", status="attempt", verbose=verbose)
        sumstats = sumstats.with_columns(pl.col(chrom).cast(pl.String).alias(chrom))
    
    # ============================================================================
    # Step 3: Identify which chromosomes need fixing
    # ============================================================================
    is_chr_fixed = pl.col(chrom).cast(pl.String).str.contains(r'^\d+$', strict=False).fill_null(False)
    fixed_count = sumstats.filter(is_chr_fixed).height
    log.log_variants_with_condition("standardized chromosome notation", fixed_count, verbose=verbose)
    
    # If all chromosomes are already numeric, skip extraction
    if fixed_count == sumstats.height:
        log.write(" -All CHR are already fixed...", verbose=verbose)
        sumstats = vchange_statusp(sumstats, is_chr_fixed, status, 4, ["986"], ["520"])
    else:
        # ========================================================================
        # Step 4: Extract and fix chromosome notations
        # ========================================================================
        # Extract chromosome using regex pattern
        # CHR_PATTERN_EXTRACT = r'^(chr)?(\d{1,3}|[XYM]|MT)$' - group 0 is full match, group 1 is chr prefix, group 2 is chromosome value
        chr_extracted = pl.col(chrom).cast(pl.String).str.extract(CHR_PATTERN_EXTRACT, 2)
        is_chr_fixable = chr_extracted.is_not_null()
        fixable_count = sumstats.filter(is_chr_fixable).height
        log.log_variants_with_condition("fixable chromosome notations", fixable_count, verbose=verbose)
        
        # Check for NA and invalid chromosomes
        is_chr_na = pl.col(chrom).is_null()
        na_count = sumstats.filter(is_chr_na).height
        if na_count > 0:
            log.log_variants_with_condition("NA chromosome notations", na_count, verbose=verbose)
        
        is_chr_invalid = (~is_chr_fixable) & (~is_chr_na)
        invalid_count = sumstats.filter(is_chr_invalid).height
        if invalid_count > 0:
            log.log_variants_with_condition("invalid chromosome notations", invalid_count, verbose=verbose)
            try:
                invalid_examples = sumstats.filter(is_chr_invalid).select(chrom).unique().head(5)
                invalid_set = set(invalid_examples[chrom].to_list())
                log.write(" -A look at invalid chromosome notations:", invalid_set, verbose=verbose)
            except:
                pass
        else:
            log.write(" -No unrecognized chromosome notations...", verbose=verbose)
        
        # ========================================================================
        # Step 5: Assign extracted values back to sumstats
        # ========================================================================
        # Update chromosome values where fixable
        sumstats = sumstats.with_columns(
            pl.when(is_chr_fixable)
            .then(chr_extracted.cast(pl.String))
            .otherwise(pl.col(chrom))
            .alias(chrom)
        )
        
        # ========================================================================
        # Step 6: Convert sex chromosomes to numeric values
        # ========================================================================
        x_label, x_num = x[0], str(x[1])
        y_label, y_num = y[0], str(y[1])
        mt_label, mt_num = mt[0], str(mt[1])
        
        # Build mapping for sex chromosomes
        sex_chr_map = {
            x_label.lower(): x_num,
            x_label.upper(): x_num,
            y_label.lower(): y_num,
            y_label.upper(): y_num,
            mt_label.lower(): mt_num,
            mt_label.upper(): mt_num,
        }
        
        # Check if any sex chromosomes exist
        chr_lower = pl.col(chrom).cast(pl.String).str.to_lowercase()
        is_sex_chr = chr_lower.is_in(list(sex_chr_map.keys()))
        sex_chr_count = sumstats.filter(is_sex_chr).height
        
        if sex_chr_count > 0:
            log.write(" -Identifying non-autosomal chromosomes : {}, {}, and {} ...".format(x_label, y_label, mt_label), verbose=verbose)
            log.write(" -Identified {} variants on sex chromosomes...".format(sex_chr_count), verbose=verbose)
            
            # Convert sex chromosomes to numeric using when/then chain
            chr_lower_col = pl.col(chrom).cast(pl.String).str.to_lowercase()
            sumstats = sumstats.with_columns(
                pl.when(chr_lower_col == x_label.lower())
                .then(pl.lit(x_num))
                .when(chr_lower_col == y_label.lower())
                .then(pl.lit(y_num))
                .when(chr_lower_col == mt_label.lower())
                .then(pl.lit(mt_num))
                .otherwise(pl.col(chrom))
                .alias(chrom)
            )
        
        # ========================================================================
        # Step 7: Update status codes
        # ========================================================================
        sumstats = vchange_statusp(sumstats, is_chr_fixed, status, 4, ["986"], ["520"])
        if fixable_count > 0:
            sumstats = vchange_statusp(sumstats, is_chr_fixable, status, 4, ["986"], ["520"])
        if invalid_count > 0:
            sumstats = vchange_statusp(sumstats, is_chr_invalid, status, 4, ["986"], ["743"])
        
        # ========================================================================
        # Step 8: Remove invalid chromosomes if requested
        # ========================================================================
        if remove:
            # Build chromosome list with numeric sex chromosomes
            chrom_list_with_numeric = chrom_list.copy()
            if chromosomes_obj is not None:
                all_sex_chr_numeric = chromosomes_obj.get_all_sex_chromosomes_numeric()
                for num_val in all_sex_chr_numeric:
                    chrom_list_with_numeric.append(str(num_val))
            else:
                chrom_list_with_numeric.extend(["23", "24", "25"])
            
            good_chr = pl.col(chrom).is_in(chrom_list_with_numeric)
            unrecognized_num = sumstats.filter(~good_chr).height
            
            if unrecognized_num > 0:
                try:
                    numeric_chrs = [int(x) for x in chrom_list_with_numeric if x.isnumeric()]
                    if numeric_chrs:
                        log.write(" -Valid CHR list: {} - {}".format(min(numeric_chrs), max(numeric_chrs)), verbose=verbose)
                except:
                    pass
                log.log_variants_removed(unrecognized_num, reason="with chromosome notations not in CHR list", verbose=verbose)
                try:
                    invalid_examples = sumstats.filter(~good_chr).select(chrom).unique().head(5)
                    invalid_set = set(invalid_examples[chrom].to_list())
                    log.write(" -A look at chromosome notations not in CHR list:", invalid_set, verbose=verbose)
                except:
                    pass
                sumstats = sumstats.filter(good_chr)

    # ============================================================================
    # Step 9: Convert chromosome column to integer type
    # ============================================================================
    try:
        sumstats = sumstats.with_columns(pl.col(chrom).cast(pl.Int64, strict=False).alias(chrom))
    except:
        # Try converting via string first
        sumstats = sumstats.with_columns(
            pl.col(chrom).cast(pl.String).str.strip_chars().cast(pl.Int64, strict=False).alias(chrom)
        )
    
    # ============================================================================
    # Step 10: Filter out variants with CHR < minchr
    # ============================================================================
    out_of_range_chr = (pl.col(chrom) < minchr) & (pl.col(chrom).is_not_null())
    out_of_range_num = sumstats.filter(out_of_range_chr).height
    if out_of_range_num > 0:
        log.write(" -Sanity check for CHR...", verbose=verbose)
        log.log_variants_removed(out_of_range_num, reason="with CHR < {}".format(minchr), verbose=verbose)
        sumstats = sumstats.filter(~out_of_range_chr)
    
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
def _fix_posp(sumstats_obj: Union['Sumstatsp', pl.DataFrame], pos: str = "POS", status: str = "STATUS", remove: bool = False, verbose: bool = True, lower_limit: int = 0, upper_limit: Optional[int] = None, limit: int = 250000000, log: Log = Log()) -> pl.DataFrame:
    '''
    Standardize and validate genomic base-pair positions using polars.
    
    This function checks that reported genomic positions fall within valid chromosomal bounds
    and optionally removes invalid entries. It handles string-formatted positions with thousands
    separators, converts positions to Int64 type, and filters out positions outside the specified
    range. If explicit limits are not provided, a default maximum bound is applied.

    Parameters
    ----------
    sumstats_obj : Sumstatsp or pl.DataFrame
        Sumstats object or polars DataFrame containing the data to fix.
    pos : str, default "POS"
        Column name for position.
    status : str, default "STATUS"
        Column name for status.
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
    log : Log, optional
        Logging object.

    Returns
    -------
    pl.DataFrame
        Polars DataFrame with standardized and validated base-pair positions.
    '''
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pl.DataFrame):
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    # Set default upper limit if not provided
    if upper_limit is None:
        upper_limit = limit
    
    # Track initial number of variants for reporting
    all_var_num = sumstats.height
    
    # Check for missing positions before processing
    is_pos_na = pl.col(pos).is_null()
    
    # Handle string types: remove thousands separators
    try:
        if sumstats[pos].dtype == pl.String:
            log.write(' -Removing thousands separator "," or underbar "_" ...', verbose=verbose)
            sumstats = sumstats.with_columns(
                pl.when(~is_pos_na)
                .then(pl.col(pos).str.replace_all(r'[,_]', ''))
                .otherwise(pl.col(pos))
                .alias(pos)
            )
    except Exception:
        pass

    # Convert POS to integer type
    try:
        log.log_datatype_change("POS", str(sumstats[pos].dtype), "Int64", status="attempt", verbose=verbose)
        sumstats = sumstats.with_columns(pl.col(pos).cast(pl.Int64, strict=False).alias(pos))
    except Exception:
        log.log_datatype_change("POS", str(sumstats[pos].dtype), "Int64", status="attempt", verbose=verbose)
        # Try converting via string first
        sumstats = sumstats.with_columns(
            pl.col(pos).cast(pl.String).str.strip_chars().cast(pl.Int64, strict=False).alias(pos)
        )
    
    # Identify fixed and invalid positions
    is_pos_na_after = pl.col(pos).is_null()
    is_pos_fixed = ~is_pos_na_after
    is_pos_invalid = (~is_pos_na) & (~is_pos_fixed)
    
    # Update status codes for fixed and invalid positions
    sumstats = vchange_statusp(sumstats, is_pos_fixed, status, 4, ["975"], ["630"])
    sumstats = vchange_statusp(sumstats, is_pos_invalid, status, 4, ["975"], ["842"])
    
    # Remove outliers outside the specified bounds
    log.write(" -Position bound:({} , {:,})".format(lower_limit, upper_limit), verbose=verbose)
    is_outlier = ((pl.col(pos) <= lower_limit) | (pl.col(pos) >= upper_limit)) & (~is_pos_na_after)
    outlier_num = sumstats.filter(is_outlier).height
    if outlier_num == 0:
        log.write(" -No outlier variants were removed.", verbose=verbose)
    else:
        log.log_variants_removed(outlier_num, reason="outliers", verbose=verbose)
    sumstats = sumstats.filter(~is_outlier)
    
    # Optionally remove remaining NA positions
    if remove is True:
        is_pos_na_final = pl.col(pos).is_null()
        sumstats = sumstats.filter(~is_pos_na_final)
        remain_var_num = sumstats.height
        log.log_variants_removed(all_var_num - remain_var_num, reason="with bad positions", verbose=verbose)
 
    # Update QC status only if called with Sumstats object
    if not is_dataframe:
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