from typing import TYPE_CHECKING, Optional, List, Tuple, Union, Dict, Any
import pandas as pd
import numpy as np
from pysam import VariantFile
import gzip
from gwaslab.io.io_fasta import load_fasta_auto, build_fasta_records, load_fasta_filtered, load_and_build_fasta_records
from itertools import repeat
from multiprocessing import Pool
from functools import partial
import re
import os
import gc
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
    from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.qc.qc_fix_sumstats import _fix_chr
from gwaslab.qc.qc_fix_sumstats import _fix_pos
from gwaslab.qc.qc_fix_sumstats import _sort_column
from gwaslab.qc.qc_fix_sumstats import _df_split
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.qc.qc_fix_sumstats import _sort_coordinate
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_NC
from gwaslab.bd.bd_common_data import _maketrans
from gwaslab.info.g_vchange_status import vchange_status, set_status_digit, status_match, ensure_status_int
from gwaslab.info.g_version import _get_version
from gwaslab.cache_manager import CacheManager, PALINDROMIC_INDEL, NON_PALINDROMIC
from gwaslab.info.g_vchange_status import STATUS_CATEGORIES
from gwaslab.io.io_vcf import check_vcf_chr_prefix, check_vcf_chr_NC
#rsidtochrpos
#checkref
#parallelizeassignrsid
#inferstrand
#parallelecheckaf

### CONSTANTS AND MAPPINGS ###
# Epsilon for floating point precision in frequency comparisons
# Used to handle floating point precision issues when comparing MAF values at threshold
FREQ_COMPARISON_EPSILON = 1e-6

### HELPER FUNCTIONS FOR STATUS OPTIMIZATION ###
def _extract_status_digit(status_series: pd.Series, digit: int) -> pd.Series:
    """
    Extract a specific digit from status codes (optimized version).
    
    Parameters
    ----------
    status_series : pd.Series
        Series of status codes (7-digit integers)
    digit : int
        Digit position (1-indexed from left: 1=leftmost, 7=rightmost)
    
    Returns
    -------
    pd.Series
        Series of extracted digits (0-9)
    """
    status_int = status_series.astype('int64')
    # For 7-digit number: digit 7 (rightmost) = status_int % 10
    #                    digit 6 = (status_int // 10) % 10
    #                    digit 5 = (status_int // 100) % 10
    #                    digit 4 = (status_int // 1000) % 10
    power = 10 ** (7 - digit)
    return (status_int // power) % 10

PADDING_VALUE = 100

# chr(0) should not be used in the mapping dict because it's a reserved value.
# Instead of starting from chr(1), we start from chr(2) because this could be useful in the future
# to compute the complementary allele with a simple XOR operation (e.g. 2 ^ 1 = 3, 3 ^ 1 = 2, 4 ^ 1 = 5, 5 ^ 1 = 4, ...)
MAPPING = {
    "A": chr(2),
    "T": chr(3),
    "C": chr(4),
    "G": chr(5),
    "N": chr(6),
}
assert all(value != chr(0) for value in MAPPING.values()), "Mapping in the dictionary should not be equal to chr(0). This is a reserved value"

_COMPLEMENTARY_MAPPING = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
}
COMPLEMENTARY_MAPPING = {k: MAPPING[v] for k,v in _COMPLEMENTARY_MAPPING.items()}

TRANSLATE_TABLE = _maketrans(MAPPING)
TRANSLATE_TABLE_COMPL = _maketrans(COMPLEMENTARY_MAPPING)

#20220808
#################################################################################################################

####################################################################################################################

def _process_chrpos_task(task: Tuple[Any, ...]) -> pd.DataFrame:
    """
    Helper function for parallel processing of chromosome position tasks.
    This must be a module-level function to be picklable for multiprocessing.
    
    Parameters
    ----------
    task : tuple
        (dataframe, h5_file_path, chr_num, group, build, status, chrom, pos, log, task_num, verbose)
    
    Returns
    -------
    DataFrame
        Updated subset of summary statistics with merged data.
    """
    df, h5_file, chr_num, group, build, status, chrom, pos, log, task_num, verbose = task
    return merge_chrpos(df, h5_file, chr_num, group, build, status, chrom, pos, log, task_num, verbose)
    
def merge_chrpos(
    sumstats_part: pd.DataFrame,
    h5_file_path: str,
    chr_num: Optional[int],
    group: int,
    build: str,
    status: str,
    chrom: str = "CHR",
    pos: str = "POS",
    log: Optional[Log] = None,
    task_num: Optional[int] = None,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Merge chromosome position data from HDF5 group into sumstats_part using index-based matching.
    
    Simple and efficient version using pandas index operations:
    - Loads reference data from HDF5 with rsn as index
    - Uses pandas index intersection for fast matching
    - Vectorized assignment
    - Supports modulo 10 grouping (group_0 through group_9)
    - CHR extracted from filename (not stored in HDF5 since each file is chromosome-specific)
    - Uses rsn as index for efficient matching

    Parameters
    ----------
    sumstats_part : DataFrame
        Subset of summary statistics to update. Must have "rsn" as index.
    h5_file_path : str
        Path to the chromosome-specific HDF5 file.
    chr_num : int or None
        Chromosome number (extracted from filename). None for old single-file format.
    group : int
        Group number (0-9 for mod10 grouping).
    build : str
        Genome build version.
    status : str
        Column name for status codes.
    chrom : str, default="CHR"
        Column name for chromosome values.
    pos : str, default="POS"
        Column name for position values.
    log : Log, optional
        Logging object for progress output.
    task_num : int, optional
        Task number for progress tracking.
    verbose : bool, default=True
        Verbose output.

    Returns
    -------
    DataFrame
        Updated subset of summary statistics with merged data (rsn as index).
    """
    # Print task number if provided (no time, space at end)
    if log is not None and task_num is not None:
        log.write("{}".format(task_num), verbose=verbose, show_time=False, end=" ")
    
    # Validate group is within expected range (0-9 for mod10 grouping)
    if group is None or not (0 <= group <= 9):
        return sumstats_part
    
    # Ensure rsn is the index
    if sumstats_part.index.name != "rsn":
        if "rsn" in sumstats_part.columns:
            sumstats_part = sumstats_part.set_index("rsn")
        else:
            # Cannot proceed without rsn
            return sumstats_part
    
    try:
        # Use HDFStore context manager for efficient reading
        with pd.HDFStore(h5_file_path, mode='r') as store:
            key = "group_{}".format(group)
            if key in store:
                # Read reference data as DataFrame
                ref_df = store[key]
                
                # Ensure rsn is the index in reference DataFrame
                if ref_df.index.name != "rsn":
                    if "rsn" in ref_df.columns:
                        ref_df = ref_df.set_index("rsn")
                    else:
                        # Cannot find rsn, skip this group
                        return sumstats_part
                
                # Find intersection of indices (matched variants)
                matched_rsn = sumstats_part.index.intersection(ref_df.index)
                
                if len(matched_rsn) > 0:
                    # Ensure STATUS is integer before assignment
                    sumstats_part = ensure_status_int(sumstats_part, status)
                    
                    # Update status codes for fixable variants
                    sumstats_part.loc[matched_rsn, status] = vchange_status(
                        sumstats_part.loc[matched_rsn, status], 1, "139", 3*int(build[0]))
                    sumstats_part.loc[matched_rsn, status] = vchange_status(
                        sumstats_part.loc[matched_rsn, status], 2, "987", 3*int(build[1]))
                    
                    # Assign CHR from parameter (only if chr_num is provided and CHR is missing)
                    if chr_num is not None:
                        if chrom in sumstats_part.columns:
                            # Only assign CHR if it's missing (preserve existing values)
                            # Get CHR values for matched variants
                            matched_chr = sumstats_part.loc[matched_rsn, chrom]
                            # Find which ones are missing
                            missing_chr_indices = matched_chr[matched_chr.isna()].index
                            if len(missing_chr_indices) > 0:
                                sumstats_part.loc[missing_chr_indices, chrom] = chr_num
                        else:
                            sumstats_part[chrom] = pd.Series(dtype="Int64", index=sumstats_part.index)
                            sumstats_part.loc[matched_rsn, chrom] = chr_num
                    
                    # Assign POS values from reference DataFrame (only if POS is missing)
                    if pos in sumstats_part.columns:
                        # Only assign POS if it's missing (preserve existing values)
                        # Get POS values for matched variants
                        matched_pos = sumstats_part.loc[matched_rsn, pos]
                        # Find which ones are missing
                        missing_pos_indices = matched_pos[matched_pos.isna()].index
                        if len(missing_pos_indices) > 0:
                            sumstats_part.loc[missing_pos_indices, pos] = ref_df.loc[missing_pos_indices, "POS"].astype("Int64")
                    else:
                        sumstats_part[pos] = pd.Series(dtype="Int64", index=sumstats_part.index)
                        sumstats_part.loc[matched_rsn, pos] = ref_df.loc[matched_rsn, "POS"].astype("Int64")
    except Exception as e:
        # Silently pass on errors (group might not exist or other issues)
        # In production, you might want to log this for debugging
        pass
    return sumstats_part

@with_logging(
    start_to_msg="assign CHR and POS using rsIDs",
    finished_msg="assigning CHR and POS using rsIDs",
    start_cols=["rsID"],
    start_function=".rsid_to_chrpos2()"
)
def _parallelize_rsid_to_chrpos(
    sumstats: Union['Sumstats', pd.DataFrame],
    rsid: str = "rsID",
    chrom: str = "CHR",
    pos: str = "POS",
    path: Optional[str] = None,
    ref_rsid_to_chrpos_vcf: Optional[str] = None,
    ref_rsid_to_chrpos_hdf5: Optional[str] = None,
    build: str = "99",
    status: str = "STATUS",
    threads: int = 4,
    block_size: Optional[int] = None,
    verbose: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    """
    Assign CHR and POS using rsIDs (uses fast HDF5-based parallel processing).

    This function uses optimized HDF5-based parallel processing which is much faster than
    the old TSV-based approach. It requires an HDF5 reference file (generated from VCF using
    `process_vcf_to_hfd5()`) or will auto-generate the path from a VCF reference file.

    This function assigns CHR and POS values to summary statistics by matching rsIDs against
    reference HDF5 files (one per chromosome) containing rsID to POS mappings. The HDF5 files
    are typically generated from a VCF file using `process_vcf_to_hfd5()` and contain
    precomputed POS values grouped by modulo 10 (rsID % 10). CHR is extracted from the filename.
    The function processes data in parallel using multiple CPU cores for improved performance
    on large datasets.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Input summary statistics DataFrame containing variant data.
    rsid : str, default="rsID"
        Column name containing rsID values in `sumstats`.
    chrom : str, default="CHR"
        Column name for chromosome values to be updated.
    pos : str, default="POS"
        Column name for position values to be updated.
    path : str, optional
        Path to the HDF5 reference file. If not provided, must specify either `ref_rsid_to_chrpos_vcf` 
        or `ref_rsid_to_chrpos_hdf5`.
    ref_rsid_to_chrpos_vcf : str, optional
        Path to VCF file containing rsID to CHR:POS mappings. If provided, the corresponding HDF5 file 
        path will be automatically generated using the mod10 naming convention.
    ref_rsid_to_chrpos_hdf5 : str, optional
        Path to pre-generated HDF5 file containing rsID to CHR:POS mappings. Takes precedence over 
        `ref_rsid_to_chrpos_vcf`.
    build : str, default="99"
        Reference genome build identifier. "99" indicates unknown or unspecified build.
    status : str, default="STATUS"
        Column name for status codes in `sumstats`.
    threads : int, default=4
        Number of threads to use for parallel processing.
    block_size : int, optional
        **Deprecated**: This parameter is ignored. The function now uses modulo 10 grouping 
        (rsID % 10) to match the HDF5 file structure created by `process_vcf_to_hfd5()`.
    verbose : bool, default=True
        If True, print progress messages.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object for recording process information.

    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with CHR and POS values assigned based on rsID matches. 
        Variants without valid rsIDs or matches will retain original values (or NaN if new columns were 
        created).

    Notes
    -----
    - The HDF5 reference files are organized as one file per chromosome with groups named
      "group_0" through "group_9" (based on rsID % 10), each containing a DataFrame with
      columns "rsn" (int64) and "POS" (int32). CHR is extracted from the filename.
    - This function uses modulo 10 grouping to match the HDF5 structure created by
      `process_vcf_to_hfd5()`.
    - The HDF5 files follow the naming convention:
      `{vcf_file_name}.chr{chr_num}.rsID_CHR_POS_mod10.h5`
    - Non-valid rsIDs (containing non-numeric characters after "rs") are processed separately.
    - Optimized with vectorized operations and efficient HDF5 access using pandas index operations.
    """
    import pandas as pd
    import os
    import glob
    import re
    from multiprocessing import Pool
    
    # ============================================================================
    # Step 1: Handle input format (DataFrame or Sumstats object)
    # ============================================================================
    if isinstance(sumstats, pd.DataFrame):
        is_dataframe = True
    else:
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # ============================================================================
    # Step 2: Determine HDF5 file/directory path
    # ============================================================================
    if ref_rsid_to_chrpos_hdf5 is not None:
        path = ref_rsid_to_chrpos_hdf5
    elif ref_rsid_to_chrpos_vcf is not None:
        vcf_file_name = os.path.basename(ref_rsid_to_chrpos_vcf)
        vcf_dir_path = os.path.dirname(ref_rsid_to_chrpos_vcf)
        if vcf_dir_path == "":
            vcf_dir_path = "."
        # process_vcf_to_hfd5 returns a directory with separate files per chromosome
        path = vcf_dir_path
    
    if path is None:
        raise ValueError(
            "Please provide path to HDF5 file or VCF file "
            "(via ref_rsid_to_chrpos_hdf5 or ref_rsid_to_chrpos_vcf)."
        )
    
    # ============================================================================
    # Step 3: Preprocess rsIDs - extract numeric part
    # ============================================================================
    # Strip first 2 characters and convert to numeric (simple and fast)
    rsid_processed = sumstats[rsid].astype(str).str[2:]
    sumstats["rsn"] = pd.to_numeric(rsid_processed, errors="coerce").astype("Int64")
    
    # Log processing parameters
    log.write(" -Source hdf5 file: ", path, verbose=verbose)
    log.write(" -Threads to use: ", threads, verbose=verbose)
    log.write(" -Grouping method: modulo 10 (rsID % 10)", verbose=verbose)
    
    # ============================================================================
    # Step 4: Separate valid and invalid rsIDs
    # ============================================================================
    input_columns = sumstats.columns
    rsn_isna = sumstats["rsn"].isna()
    
    sumstats_nonrs = sumstats.loc[rsn_isna, :].copy()
    sumstats_rs = sumstats.loc[~rsn_isna, :].copy()
    
    log.write(" -Non-Valid rsIDs: ", sumstats_nonrs.shape[0], verbose=verbose)
    log.write(" -Valid rsIDs: ", len(sumstats_rs), verbose=verbose)
    
    # Early return if no valid rsIDs
    if len(sumstats_rs) == 0:
        log.write(" -No valid rsIDs to process", verbose=verbose)
        sumstats_nonrs = sumstats_nonrs.drop(columns=["rsn"], errors='ignore')
        if not is_dataframe:
            sumstats_obj.data = sumstats_nonrs
            return sumstats_obj.data
        else:
            return sumstats_nonrs
    
    # ============================================================================
    # Step 5: Find HDF5 chromosome files and create mapping
    # ============================================================================
    # Build mapping of chromosome number -> HDF5 file path
    chr_to_h5_file = {}
    if os.path.isdir(path):
        # Look for chromosome-specific files (new design)
        h5_files = glob.glob(os.path.join(path, "*.chr*.rsID_CHR_POS_mod10.h5"))
        if h5_files:
            # New format: extract chromosome number from each filename
            for h5_file in h5_files:
                match = re.search(r'\.chr(\d+)\.', os.path.basename(h5_file))
                if match:
                    chr_to_h5_file[int(match.group(1))] = h5_file
        else:
            # Fallback: try old single-file format
            h5_files = glob.glob(os.path.join(path, "*.rsID_CHR_POS_mod10.h5"))
            if h5_files:
                chr_to_h5_file[None] = h5_files[0]
        if not chr_to_h5_file:
            raise ValueError(f"No HDF5 files found in directory: {path}")
    else:
        # Path is a single file (backward compatibility)
        chr_to_h5_file[None] = path
    
    log.write(" -Found HDF5 files for {} chromosome(s)".format(len(chr_to_h5_file)), verbose=verbose)
    
    # ============================================================================
    # Step 6: Initialize columns and prepare for processing
    # ============================================================================
    # Initialize CHR and POS columns if they don't exist
    if chrom not in input_columns:
        log.write(" -Initiating CHR ... ", verbose=verbose)
        sumstats_rs[chrom] = pd.Series(dtype="Int64", index=sumstats_rs.index)
    if pos not in input_columns:
        log.write(" -Initiating POS ... ", verbose=verbose)
        sumstats_rs[pos] = pd.Series(dtype="Int64", index=sumstats_rs.index)
    
    # Assign groups and set rsn as index for efficient matching
    sumstats_rs["group"] = sumstats_rs["rsn"] % 10
    sumstats_rs = sumstats_rs.set_index("rsn")
    
    # ============================================================================
    # Step 7: Create processing tasks
    # ============================================================================
    # Check if CHR column has values
    has_chr = chrom in sumstats_rs.columns and sumstats_rs[chrom].notna().any()
    
    tasks = []
    if has_chr and None not in chr_to_h5_file:
        # Process by chromosome: group data by chromosome, then by group
        log.write(" -Processing by chromosome (CHR column available)", verbose=verbose)
        for chr_num, chr_df in sumstats_rs.groupby(chrom):
            if chr_num in chr_to_h5_file:
                h5_file = chr_to_h5_file[chr_num]
                for group_id, group_df in chr_df.groupby("group"):
                    tasks.append((group_df, h5_file, int(chr_num), int(group_id), build, status, chrom, pos, log, len(tasks) + 1, verbose))
    else:
        # Search across all chromosomes: group by group, try all chromosome files
        log.write(" -Searching across all chromosomes (CHR column not available or single-file format)", verbose=verbose)
        for group_id, group_df in sumstats_rs.groupby("group"):
            for chr_num, h5_file in chr_to_h5_file.items():
                tasks.append((group_df, h5_file, chr_num, int(group_id), build, status, chrom, pos, log, len(tasks) + 1, verbose))
    
    log.write(" -Created {} processing tasks: ".format(len(tasks)), verbose=verbose)
    
    # ============================================================================
    # Step 8: Parallel lookup and assignment
    # ============================================================================
    if tasks:
        with Pool(threads) as pool:
            sumstats_rs_list = pool.map(_process_chrpos_task, tasks)
        
        # Print newline after all tasks complete
        log.write("", verbose=verbose, show_time=False, end="\n")
        
        # Concatenate preserving rsn index
        sumstats_rs = pd.concat(sumstats_rs_list)
        
        # Deduplicate if we searched across all chromosomes (keep matches with CHR/POS)
        if not has_chr or None in chr_to_h5_file:
            # Sort by match status (matched variants first), then deduplicate
            if chrom in sumstats_rs.columns and pos in sumstats_rs.columns:
                sumstats_rs['_match_priority'] = (sumstats_rs[chrom].notna() & sumstats_rs[pos].notna()).astype(int)
                sumstats_rs = sumstats_rs.sort_values(by='_match_priority', ascending=False, kind='mergesort')
                sumstats_rs = sumstats_rs.drop(columns=['_match_priority'])
            sumstats_rs = sumstats_rs[~sumstats_rs.index.duplicated(keep='first')]
        
        # Log match statistics
        if chrom in sumstats_rs.columns and pos in sumstats_rs.columns:
            matched = (sumstats_rs[chrom].notna() & sumstats_rs[pos].notna()).sum()
            log.write(" -Variants matched: {} / {}".format(matched, len(sumstats_rs)), verbose=verbose)
    else:
        log.write(" -No tasks to process", verbose=verbose)
    
    # ============================================================================
    # Step 9: Merge results and cleanup
    # ============================================================================
    log.write(" -Merging group data... ", verbose=verbose)
    
    # Drop temporary columns and reset index
    sumstats_rs = sumstats_rs.drop(columns=["group"], errors='ignore').reset_index()
    sumstats_nonrs = sumstats_nonrs.drop(columns=["rsn"], errors='ignore')
    
    # Merge back valid and invalid rsIDs
    log.write(" -Append data... ", verbose=verbose)
    sumstats = pd.concat([sumstats_rs, sumstats_nonrs], ignore_index=True)
    
    # ============================================================================
    # Step 10: Post-processing and validation
    # ============================================================================
    sumstats = _fix_chr(sumstats, verbose=verbose)
    sumstats = _fix_pos(sumstats, verbose=verbose)
    sumstats = _sort_column(sumstats_obj=sumstats, verbose=verbose)

    # ============================================================================
    # Step 11: Update metadata and return results
    # ============================================================================
    if not is_dataframe:
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            rsid_to_chrpos_kwargs = {
                'rsid': rsid,
                'chrom': chrom,
                'pos': pos,
                'path': path,
                'ref_rsid_to_chrpos_vcf': ref_rsid_to_chrpos_vcf,
                'ref_rsid_to_chrpos_hdf5': ref_rsid_to_chrpos_hdf5,
                'build': build,
                'status': status,
                'threads': threads,
                'block_size': block_size
            }
            _update_harmonize_step(sumstats_obj, "rsid_to_chrpos", rsid_to_chrpos_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats
####################################################################################################################
# old version
#20240320 check if non-effect allele is aligned with reference genome         
def _fast_check_status(
    x: pd.DataFrame,
    record: np.ndarray,
    starting_positions: np.ndarray,
    records_len: np.ndarray,
    powers_of_10: np.ndarray,
    chrom_order: np.ndarray
) -> np.ndarray:
    """
    Check if non-effect allele (NEA) is aligned with reference genome using vectorized operations.

    This function efficiently checks alignment status for multiple variants against a reference genome
    using numpy vectorization. It updates the status codes based on whether alleles match the reference,
    need flipping, or are not present in the reference.

    Parameters
    ----------
    x : pd.DataFrame
        DataFrame containing variant data with columns in order: ['CHR', 'POS', 'EA', 'NEA', 'STATUS'].
        Note: Column names are not explicitly checked - the function expects specific column ordering.
    record : np.array
        Concatenated reference genome sequence as a 1D numpy array of integers (after translation).
        Note: This should already be padded in the caller to avoid out-of-bounds errors.
    starting_positions : np.array
        1D array containing the starting index in `record` for each chromosome in `x`.
        Must correspond to chromosomes in `chrom_order` and be ordered accordingly (sorted).
    records_len : np.array
        1D array containing the length of each chromosome's sequence in `record`.
        Must correspond to chromosomes in `chrom_order` and be ordered accordingly (sorted).
    powers_of_10 : np.array
        Pre-computed array of powers of 10 for status digit extraction: [1000000, 100000, ..., 1].
        This avoids recreating it for each call.
    chrom_order : np.array
        Sorted array of unique chromosome values in `x`. This ensures alignment between
        chromosome values and array indices in `starting_positions` and `records_len`.
        This avoids recomputing np.unique() and fixes a correctness bug.

    Returns
    -------
    np.array
        Array of updated status strings after checking against the reference genome.

    Notes
    -----
    Status codes are stored as the 6th character (0-based index 5) of the status string:
    0: Alleles match reference (NEA == ref)
    1: Flipped Fixed
    2: Reverse_complementary Fixed
    3: Flipped (EA == ref but NEA != ref)
    4: Reverse_complementary (rev_NEA == ref)
    5: Reverse_complementary + Flipped
    6: Both alleles on genome but indistinguishable (indel)
    7: Reverse_complementary + both alleles on genome + indistinguishable
    8: Not on reference genome
    9: Unchecked
    """
    if x.empty:
        return np.array([])
    
    # x is expected to be a DataFrame with these columns in that order: ['CHR', 'POS', 'EA', 'NEA', 'STATUS']
    # In this way, we don't need to specify the columns names
    _chrom = x.iloc[:, 0]
    _pos = x.iloc[:, 1]
    _ea = x.iloc[:, 2]
    _nea = x.iloc[:, 3]
    _status = x.iloc[:, 4]

    # position of the status (i.e. x['STATUS']) that will be modified
    status_flip_idx = 5 

    pos = _pos.values.astype(np.int64) # convert to int64 because they could be of type 'object'

    # Rebase the chromosome numbers to 0-based indexing
    # e.g. ['1', '2', '4', '2'] -> [0, 1, 2, 1]
    # This is needed because record is a single 1D array containing all the records for all the selected chromosomes,
    # so for instance if record contains the records for chr1, chr2, chr4 ([...chr1...chr2...chr4...]), we need to 
    # rebase the chromosome numbers to 0-based indexing to index the correct record portion when we do starting_positions[chrom]
    # Note that in x there are only the rows for the same chromosomes for which we have the records in record
    # (i.e. we don't have rows for chr3 if we don't have the record for chr3). This filtering is done in the caller function
    # Use precomputed chrom_order from caller to avoid recomputing np.unique() and ensure correct alignment
    _chrom = _chrom.values
    chrom = np.searchsorted(chrom_order, _chrom) # Replace each value in '_chrom' with its corresponding index in chrom_order
 
    max_len_nea = _nea.str.len().max()
    max_len_ea = _ea.str.len().max()

    ########################################## mask for variants with out of range POS
    mask_outlier = pos > records_len[chrom]

    #########################################

    # Let's apply the same magic used for the fasta records (check build_fasta_records() for details) to convert the NEA and EA to
    # a numpy array of integers in a very fast way.
    # In that case we start from a pd.Series to we can apply some built-in methods.
    # Also, when doing nea.view('<u4'), each row will be automatically right-padded with zeros to reach the max_len_nea.
    # For this reason, we then replace the zeros with out padding value
    # (and that's why the mapping dict can't have chr(0) as a value, otherwise we would have zeros for both padding and a character)
    # Reshaping is needed because .view('<u4') will create a flattened array    
    nea = _nea.str.translate(TRANSLATE_TABLE).to_numpy().astype(f'<U{max_len_nea}')
    nea = nea.view('<u4').reshape(-1, max_len_nea).astype(np.uint8)
    nea[nea == 0] = PADDING_VALUE # padding value
    ###########################################
    
    ###########################################
    # Create a mask holding True at the position of non-padding values
    mask_nea = nea != PADDING_VALUE

    # Let's do everything again for EA
    ea = _ea.str.translate(TRANSLATE_TABLE).to_numpy().astype(f'<U{max_len_ea}')
    ea = ea.view('<u4').reshape(-1, max_len_ea).astype(np.uint8)
    ea[ea == 0] = PADDING_VALUE # padding value
    ###########################################
    
    ###########################################
    mask_ea = ea != PADDING_VALUE


    # Optimized status conversion: work directly with integers instead of string conversion
    # Status codes are always 7-digit integers
    # Extract digits directly using integer operations (much faster than string conversion)
    assert (_status >= 1000000).all() and (_status <= 9999999).all(), "All status codes should be 7-digit integers"
    status_len = len(powers_of_10)  # Use length from pre-computed powers_of_10
    
    # Convert status to numpy array and extract digits using vectorized operations
    # powers_of_10 is pre-computed in check_status to avoid duplication
    status_vals = _status.values.astype(np.int64)
    # Extract all digits at once using broadcasting: (status // power) % 10
    status = ((status_vals[:, np.newaxis] // powers_of_10[np.newaxis, :]) % 10).astype(np.uint8)


    # Expand the position to a 2D array and subtract 1 to convert to 0-based indexing
    # e.g. [2, 21, 46] -> [[1], [20], [45]]
    pos = np.expand_dims(pos, axis=-1) - 1

    # Create a modified indices array specifying the starting position of each chromosome in the concatenated record array
    modified_indices = starting_positions[chrom]
    modified_indices = modified_indices[:, np.newaxis] # Add a new axis to modified_indices to align with the dimensions of pos

    # Create the range of indices: [0, ..., max_len_nea-1]
    indices_range = np.arange(max_len_nea)

    # Add the range of indices to the starting indices
    # e.g. pos = [[1], [20], [45]], indices_range = [0, 1, 2], indices = [[1, 2, 3], [20, 21, 22], [45, 46, 47]]
    indices = pos + indices_range

    # Modify indices to select the correct absolute position in the concatenated record array
    indices = indices + modified_indices

    # Note: record is already padded in check_status() to avoid out-of-bounds errors
    # This avoids duplicating the padding operation for each call
    
    # Index the record array using the computed indices.
    # Since we use np.take, indices must all have the same length, and this is why we added the padding to NEA
    # and we create the indices using max_len_nea (long story short, we can't obtain a scattered/ragged array)
    output_nea = np.take(record, indices, mode="clip")
    ##################################################################
    output_nea[mask_outlier] = PADDING_VALUE
    ##################################################################
    
    # Check if the NEA is equal to the reference sequence at the given position
    # In a non-matrix way, this is equivalent (for one single element) to:
    # nea == record[pos-1: pos+len(nea)-1]
    # where for example:
    #  a) nea = "AC", record = "ACTG", pos = 1 -> True
    #  b) nea = "T", record = "ACTG", pos = 3 -> True
    #  c) nea = "AG", record = "ACTG", pos = 1 -> False
    # Since we want to do everything in a vectorized way, we will compare the padded NEA with the output 
    # and then we use the mask to focus only on the non-padded elements
    # Pseudo example (X represents the padding value):
    #  nea = ['AC', 'T'], record = 'ACTGAAG', pos = [1, 3]
    #  -> nea = ['AC', 'TX'], indices = [[1, 2], [3, 4]], mask = [[True, True], [True, False]], output_nea = [['A', 'C'], ['T', 'G']]
    #  -> nea == output_nea: [[True, True], [True, False]], mask: [[True, True], [True, False]]
    #  -> nea == output_nea + ~mask: [[True, True], [True, True]]
    #  -> np.all(nea == output_nea + ~mask, 1): [True, True]

    nea_eq_ref = np.all((nea == output_nea) + ~mask_nea, 1)

    # Let's do everything again for EA
    indices_range = np.arange(max_len_ea)
    indices = pos + indices_range
    indices = indices + modified_indices
    output_ea = np.take(record, indices, mode="clip")
    ##################################################################
    output_ea[mask_outlier] = PADDING_VALUE
    ##################################################################

    ea_eq_ref = np.all((ea == output_ea) + ~mask_ea, 1)

    # Only compute reverse complements for rows where both nea and ea don't match reference
    # This optimization avoids unnecessary computation for most straightforward SNP matches
    need_rev = (~nea_eq_ref) & (~ea_eq_ref)
    
    # Initialize reverse complement arrays with False (will be updated only for need_rev rows)
    rev_nea_eq_ref = np.zeros(len(nea_eq_ref), dtype=bool)
    rev_ea_eq_ref = np.zeros(len(ea_eq_ref), dtype=bool)
    
    if need_rev.any():
        # Build reverse complements only for rows that need it
        _nea_need_rev = _nea[need_rev]
        _ea_need_rev = _ea[need_rev]
        output_nea_need_rev = output_nea[need_rev]
        output_ea_need_rev = output_ea[need_rev]
        mask_nea_need_rev = mask_nea[need_rev]
        mask_ea_need_rev = mask_ea[need_rev]
        
        # Create the reverse complement of NEA for rows that need it
        # In this case, we manually left-pad the translated string with the padding value, since the padding done by view('<u4') would be right-padded
        # and that will make hard the reverse operation (because we would have e.g. [2, 2, 4, 100, ..., 100] which will be hard to convert into [4, 2, 2, 100, ..., 100])
        rev_nea_subset = _nea_need_rev.str.translate(TRANSLATE_TABLE_COMPL).str.pad(max_len_nea, 'left', chr(PADDING_VALUE)).to_numpy().astype(f'<U{max_len_nea}')
        rev_nea_subset = rev_nea_subset.view('<u4').reshape(-1, max_len_nea).astype(np.uint8)
        rev_nea_subset = rev_nea_subset[:, ::-1]
        
        # Create the reverse complement of EA for rows that need it
        rev_ea_subset = _ea_need_rev.str.translate(TRANSLATE_TABLE_COMPL).str.pad(max_len_ea, 'left', chr(PADDING_VALUE)).to_numpy().astype(f'<U{max_len_ea}')
        rev_ea_subset = rev_ea_subset.view('<u4').reshape(-1, max_len_ea).astype(np.uint8)
        rev_ea_subset = rev_ea_subset[:, ::-1]
        
        # Check reverse complements only for the subset
        rev_nea_eq_ref_subset = np.all((rev_nea_subset == output_nea_need_rev) + ~mask_nea_need_rev, 1)
        rev_ea_eq_ref_subset = np.all((rev_ea_subset == output_ea_need_rev) + ~mask_ea_need_rev, 1)
        
        # Update the full arrays with results for the subset
        rev_nea_eq_ref[need_rev] = rev_nea_eq_ref_subset
        rev_ea_eq_ref[need_rev] = rev_ea_eq_ref_subset

    masks_max_len = max(mask_nea.shape[1], mask_ea.shape[1])

    len_nea_eq_len_ea = np.all(
        np.pad(mask_nea, ((0,0),(0, masks_max_len-mask_nea.shape[1])), constant_values=False) == 
        np.pad(mask_ea, ((0,0),(0, masks_max_len-mask_ea.shape[1])), constant_values=False)
        , axis=1) # pad masks with False to reach same shape
    len_rev_nea_eq_rev_len_ea = len_nea_eq_len_ea

    # The following conditions replicates the if-else statements of the original check_status function:
    # https://github.com/Cloufield/gwaslab/blob/f6b4c4e58a26e5d67d6587141cde27acf9ce2a11/src/gwaslab/hm_harmonize_sumstats.py#L238

    # nea == ref && ea == ref && len(nea) != len(ea)
    status[nea_eq_ref * ea_eq_ref * ~len_nea_eq_len_ea, status_flip_idx] = 6

    # nea == ref && ea != ref
    status[nea_eq_ref * ~ea_eq_ref, status_flip_idx] = 0

    # nea != ref && ea == ref
    status[~nea_eq_ref * ea_eq_ref, status_flip_idx] = 3

    # nea != ref && ea != ref && rev_nea == ref && rev_ea == ref && len(rev_nea) != len(rev_ea)
    status[~nea_eq_ref * ~ea_eq_ref * rev_nea_eq_ref * rev_ea_eq_ref * ~len_rev_nea_eq_rev_len_ea, status_flip_idx] = 8

    # nea != ref && ea != ref && rev_nea == ref && rev_ea != ref
    status[~nea_eq_ref * ~ea_eq_ref * rev_nea_eq_ref * ~rev_ea_eq_ref, status_flip_idx] = 4

    # nea != ref && ea != ref && rev_nea != ref && rev_ea == ref
    status[~nea_eq_ref * ~ea_eq_ref * ~rev_nea_eq_ref * rev_ea_eq_ref, status_flip_idx] = 5

    # nea != ref && ea != ref && rev_nea != ref && rev_ea != ref
    status[~nea_eq_ref * ~ea_eq_ref * ~rev_nea_eq_ref * ~rev_ea_eq_ref, status_flip_idx] = 8

    # Convert back the (now modified) 2D status digit array to integers efficiently
    # Since 'status' is a 2D array of integers ranging from 0 to 9, we can build the integer representation
    # of each row using the efficient operation below (e.g. [1, 2, 3, 4, 5, 6, 7] -> 1234567).
    # This avoids the string conversion overhead we had before
    # Then convert to string array for compatibility with existing code
    status_flat = np.sum(status * powers_of_10[np.newaxis, :], axis=1)
    status_arr = status_flat.astype(f'<U{status_len}')

    return status_arr


def check_status(
    sumstats: pd.DataFrame,
    fasta_records_dict: Optional[Dict[int, Any]] = None,
    record: Optional[np.ndarray] = None,
    starting_positions_dict: Optional[Dict[int, int]] = None,
    records_len_dict: Optional[Dict[int, int]] = None,
    log: Log = Log(),
    verbose: bool = True
) -> np.ndarray:
    """
    Check status of variants against reference genome.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame
    fasta_records_dict : dict, optional
        Dictionary of FastaRecord objects (for backward compatibility)
    record : np.ndarray, optional
        Pre-built concatenated numpy array of sequences (if provided, skips building)
    starting_positions_dict : dict, optional
        Pre-built dictionary of starting positions (if provided, skips building)
    records_len_dict : dict, optional
        Pre-built dictionary of record lengths (if provided, skips building)
    log : Log, optional
        Logging object
    verbose : bool, optional
        Verbose output
        
    Returns
    -------
    np.ndarray
        Status values
    """
    chrom,pos,ea,nea,status = sumstats.columns

    # Use pre-built records if provided, otherwise build from dict
    if record is None or starting_positions_dict is None or records_len_dict is None:
        # Backward compatibility: build from dict
        record, starting_positions_dict, records_len_dict = build_fasta_records(fasta_records_dict, pos_as_dict=True, log=log, verbose=verbose)

    # Pre-compute maximum padding needed across both groups to avoid padding twice
    # This is more efficient than padding the record separately for each call
    max_len_nea_all = sumstats[nea].str.len().max()
    max_len_ea_all = sumstats[ea].str.len().max()
    max_padding_needed = max(max_len_nea_all, max_len_ea_all)
    
    # Pad the record once based on maximum needed padding
    # This avoids duplicating the padding operation in _fast_check_status
    record = np.pad(record, (0, max_padding_needed), constant_values=PADDING_VALUE)
    
    # Pre-compute powers_of_10 for status conversion (used in both calls)
    # Status codes are always 7-digit integers
    status_len = 7
    powers_of_10 = 10 ** np.arange(status_len - 1, -1, -1, dtype=np.int64)

    # In _fast_check_status(), several 2D numpy arrays are created and they are padded to have shape[1] == max_len_nea or max_len_ea
    # Since most of the NEA and EA strings are short, we perform the check first on the records having short NEA and EA strings,
    # and then we perform the check on the records having long NEA and EA strings. In this way we can speed up the process (since the 
    # arrays are smaller) and save memory.
    max_len = 4 # this is a chosen value, we could compute it using some stats about the length and count of NEA and EA strings
    condition = (sumstats[nea].str.len() <= max_len) & (sumstats[ea].str.len() <= max_len)

    log.write(f"   -Checking records for ( len(NEA) <= {max_len} and len(EA) <= {max_len} )", verbose=verbose)
    sumstats_cond = sumstats[condition]
    # Sort chromosomes to ensure alignment with np.unique() inside _fast_check_status
    # This fixes a correctness bug where unsorted unique_chrom_cond could mismatch with sorted np.unique(_chrom)
    chrom_order_cond = np.sort(sumstats_cond[chrom].unique())
    starting_pos_cond = np.array([starting_positions_dict[k] for k in chrom_order_cond], dtype=np.int64)
    records_len_cond = np.array([records_len_dict[k] for k in chrom_order_cond], dtype=np.int64)

    sumstats.loc[condition, status] = _fast_check_status(sumstats_cond, record=record, starting_positions=starting_pos_cond, records_len=records_len_cond, powers_of_10=powers_of_10, chrom_order=chrom_order_cond)

    log.write(f"   -Checking records for ( len(NEA) > {max_len} or len(EA) > {max_len} )", verbose=verbose)
    sumstats_not_cond = sumstats[~condition]
    # Sort chromosomes to ensure alignment with np.unique() inside _fast_check_status
    chrom_order_not_cond = np.sort(sumstats_not_cond[chrom].unique())
    starting_not_pos_cond = np.array([starting_positions_dict[k] for k in chrom_order_not_cond], dtype=np.int64)
    records_len_not_cond = np.array([records_len_dict[k] for k in chrom_order_not_cond], dtype=np.int64)
    sumstats.loc[~condition, status] = _fast_check_status(sumstats_not_cond, record=record, starting_positions=starting_not_pos_cond, records_len=records_len_not_cond, powers_of_10=powers_of_10, chrom_order=chrom_order_not_cond)

    return sumstats[status].values

# load_fasta_auto is now imported from gwaslab.io.io_fasta       

@with_logging(
    start_to_msg="check if NEA is aligned with reference sequence",
    finished_msg="checking if NEA is aligned with reference sequence",
    start_cols=["CHR","POS","EA","NEA","STATUS"],
    start_function=".check_ref()"
)
def _check_ref(
    sumstats: Union['Sumstats', pd.DataFrame],
    ref_seq: str,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    status: str = "STATUS",
    mapper: Optional['ChromosomeMapper'] = None,
    remove: bool = False,
    verbose: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    """
    Check if non-effect allele (NEA) is aligned with reference genome.

    This function checks whether the non-effect allele (NEA) in the summary statistics
    matches the reference genome sequence. It updates the status codes in the summary
    statistics DataFrame to reflect the alignment status.

    Parameters
    ----------
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object containing variant data.
    ref_seq : str
        Path to the reference genome FASTA file.
    chrom : str, default="CHR"
        Column name for chromosome information in sumstats.
    pos : str, default="POS"
        Column name for position information in sumstats.
    ea : str, default="EA"
        Column name for effect allele in sumstats.
    nea : str, default="NEA"
        Column name for non-effect allele in sumstats.
    status : str, default="STATUS"
        Column name for status codes in sumstats.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper.
    remove : bool, default=False
        If True, remove variants not on the reference genome.
    verbose : bool, default=True
        If True, print progress messages.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object for recording process information.

    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with status codes reflecting
        alignment with reference genome.

    Notes
    -----
    The function uses the following status codes (6th character of status string):
    0: Alleles match reference (NEA == ref)
    3: Flipped (EA == ref but NEA != ref)
    4: Reverse_complementary (rev_NEA == ref)
    5: Reverse_complementary + Flipped
    6: Both alleles on genome but indistinguishable (indel)
    8: Not on reference genome
    9: Unchecked
    """
    # Handle both DataFrame and Sumstats object
    import pandas as pd
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
        sumstats_obj = None
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    # Get mapper from Sumstats object if available, otherwise create one
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            # Use the Sumstats object's mapper
            mapper = sumstats_obj.mapper
        else:
            # Create default mapper (will auto-detect format when needed)
            from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(
                species=species,
                build=build,
                log=log,
                verbose=verbose
            )
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chrom in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chrom])
    
    # Ensure STATUS is integer type before any operations
    sumstats = ensure_status_int(sumstats, status)
    
    log.write(" -Reference genome FASTA file: "+ ref_seq,verbose=verbose)  
    
    sumstats = _sort_coordinate(sumstats,verbose=False)
    
    # Convert chromlist to set for O(1) lookup instead of O(n)
    chromlist = get_chr_list(add_number=True)
    chromlist_set = set(chromlist)
    
    # Convert chroms_in_sumstats to set for O(1) lookup (compute set directly from unique)
    chroms_in_sumstats_set = set(sumstats[chrom].unique())  # load records from Fasta file only for the chromosomes present in the sumstats
    
    # Load, filter, and build FASTA records in a single pass for maximum performance
    # This avoids creating intermediate FastaRecord objects
    record, starting_positions_dict, records_len_dict = load_and_build_fasta_records(
        ref_seq,
        chromlist_set,
        chroms_in_sumstats_set,
        mapper=mapper,
        pos_as_dict=True,
        log=log,
        verbose=verbose
    )

    # Early filtering: combine all conditions in one pass
    # Cache these masks for reuse later in status counting and available_to_check computation
    valid_pos_mask = ~sumstats[pos].isna() & (sumstats[pos] > 0)
    valid_alleles_mask = ~sumstats[nea].isna() & ~sumstats[ea].isna()
    
    # Track if we actually checked any records
    checked_any_records = False
    
    if len(records_len_dict) > 0:
        log.write(" -Checking records", verbose=verbose)
        # Convert dict keys to set for faster isin() operation
        all_records_keys_set = set(records_len_dict.keys())
        valid_chrom_mask = sumstats[chrom].isin(all_records_keys_set)
        
        to_check_ref = valid_chrom_mask & valid_pos_mask & valid_alleles_mask
        
        if to_check_ref.any():
            sumstats_to_check = sumstats.loc[to_check_ref,[chrom,pos,ea,nea,status]]
            # Pass pre-built records directly to check_status for better performance
            sumstats.loc[to_check_ref,status] = check_status(
                sumstats_to_check, 
                record=record,
                starting_positions_dict=starting_positions_dict,
                records_len_dict=records_len_dict,
                log=log, 
                verbose=verbose
            )
            checked_any_records = True
        log.write(" -Finished checking records", verbose=verbose) 
    
    # Convert STATUS to integer only if we checked records or need to compute statistics
    # Optimize: check dtype once and convert directly, reuse status_int for digit extraction
    if checked_any_records or not remove:
        status_dtype = sumstats[status].dtype
        if status_dtype.name == 'category':
            status_int = sumstats[status].astype(str).astype('int64')
            sumstats[status] = status_int.astype('Int64')
        elif status_dtype not in ['int64', 'Int64', 'int32', 'Int32']:
            status_int = sumstats[status].astype('int64')
            sumstats[status] = status_int.astype('Int64')
        else:
            # Already integer type, just ensure Int64
            status_int = sumstats[status].astype('int64')
            if status_dtype != 'Int64':
                sumstats[status] = status_int.astype('Int64')
    else:
        # No records checked, but we still need status_int for remove operation
        status_int = sumstats[status].astype('int64')

    # Compute available_to_check: variants with valid pos, nea, and ea
    # Reuse the masks we already computed
    available_to_check = (valid_pos_mask & valid_alleles_mask).sum()
    
    # Optimize: Extract 6th digit once instead of calling status_match 7 times
    # Digit 6 is the 2nd from right (1-indexed: 7=rightmost, 6=2nd from right)
    # For 7-digit number: extract using (status // 10) % 10
    digit_6 = (status_int // 10) % 10  # Extract 6th digit (2nd from right)
    
    # Use vectorized comparisons instead of multiple status_match calls
    status_0 = (digit_6 == 0).sum()
    status_3 = (digit_6 == 3).sum()
    status_4 = (digit_6 == 4).sum()
    status_5 = (digit_6 == 5).sum()
    status_6 = (digit_6 == 6).sum()
    #status_7 = (digit_6 == 7).sum()
    status_8 = (digit_6 == 8).sum()
    
    log.write(" -Variants allele on given reference sequence : ",status_0,verbose=verbose)
    log.write(" -Variants flipped : ",status_3,verbose=verbose)
    
    # Avoid division by zero
    if available_to_check > 0:
        raw_matching_rate = (status_3+status_0)/available_to_check
        flip_rate = status_3/available_to_check
        log.write("  -Raw Matching rate : ","{:.2f}%".format(raw_matching_rate*100),verbose=verbose)
        if raw_matching_rate <0.8:
            log.warning("Matching rate is low, please check if the right reference genome is used.")
        if flip_rate > 0.85 :
            log.write("  -Flipping variants rate > 0.85, it is likely that the EA is aligned with REF in the original dataset.",verbose=verbose)
    else:
        log.write("  -No variants available to check against reference.",verbose=verbose)
    
    log.write(" -Variants inferred reverse_complement : ",status_4,verbose=verbose)
    log.write(" -Variants inferred reverse_complement_flipped : ",status_5,verbose=verbose)
    log.write(" -Both allele on genome + unable to distinguish : ",status_6,verbose=verbose)
    #log.write(" -Reverse_complementary + both allele on genome + unable to distinguish: ",status_7)
    log.write(" -Variants not on given reference sequence : ",status_8,verbose=verbose)
    
    if remove is True:
        # Use the pre-computed digit_6 instead of calling status_match again
        sumstats = sumstats.loc[digit_6 != 8, :]
        log.write(" -Variants not on given reference sequence were removed.",verbose=verbose)

    # Update harmonization status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            check_ref_kwargs = {
                'ref_seq': ref_seq, 'chrom': chrom, 'pos': pos, 'ea': ea, 'nea': nea,
                'status': status, 'mapper': mapper, 'remove': remove
            }
            _update_harmonize_step(sumstats_obj, "check_ref", check_ref_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats


#######################################################################################################################################

#20220721
def chrposref_rsid(
    chr: Union[str, int],
    end: int,
    ref: str,
    alt: str,
    vcf_reader: VariantFile,
    mapper: Optional['ChromosomeMapper'] = None
) -> Union[str, pd.NA]:
    ## single record assignment
    start=end-1
    if mapper is not None:
        # Convert sumstats chr to reference format for VCF lookup
        # as_string=True for VCF compatibility (pysam requires string)
        chr = mapper.sumstats_to_reference(chr, reference_file=None, as_string=True)
    
    try:
        chr_seq = vcf_reader.fetch(chr,start,end)
    except:
        return pd.NA
    
    for record in chr_seq:
        if record.pos==end: 
            if record.alts is None:
                return pd.NA
            if record.ref==ref and (alt in record.alts):
                return record.id
            elif (ref in record.alts) and record.ref==alt:
                return record.id
    return pd.NA

def assign_rsid_single(
    sumstats: pd.DataFrame,
    path: str,
    rsid: str = "rsID",
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    mapper: Optional['ChromosomeMapper'] = None
) -> pd.DataFrame:
    ## single df assignment
    vcf_reader = VariantFile(path)
    # Auto-detect reference format from VCF file if mapper provided
    if mapper is not None:
        mapper.detect_reference_format(path)
    def rsid_helper(x,vcf_reader,mapper):
         return chrposref_rsid(x.iloc[0],x.iloc[1],x.iloc[2],x.iloc[3],vcf_reader,mapper)
    map_func=partial(rsid_helper,vcf_reader=vcf_reader,mapper=mapper)
    rsID = sumstats.apply(map_func,axis=1)
    return rsID

@with_logging(
    start_to_msg="assign rsID using reference file",
    finished_msg="assign rsID using reference file",
    start_cols=["CHR","POS","EA","NEA","STATUS"],
    start_function=".assign_rsid()"
)
def _parallelize_assign_rsid(
    sumstats: Union['Sumstats', pd.DataFrame],
    path: str,
    ref_mode: str = "vcf",
    snpid: str = "SNPID",
    rsid: str = "rsID",
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    status: str = "STATUS",
    threads: int = 1,
    chunksize: int = 5000000,
    ref_snpid: str = "SNPID",
    ref_rsid: str = "rsID",
    overwrite: str = "empty",
    verbose: bool = True,
    log: Log = Log(),
    mapper: Optional['ChromosomeMapper'] = None
) -> pd.DataFrame:
    """
    Assign rsID to variants by matching with reference file.

    This function assigns rsID to variants in the summary statistics by matching
    chromosome position and alleles with a reference file (VCF or TSV format).
    It supports different overwrite modes and can process data in parallel.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame containing variant data.
    path : str
        Path to the reference file (VCF or TSV).
    ref_mode : str, default="vcf"
        Reference file format ("vcf" for VCF files, "tsv" for TSV files).
    snpid : str, default="SNPID"
        Column name for SNP IDs in the summary statistics.
    rsid : str, default="rsID"
        Column name for rsIDs in the summary statistics.
    chr : str, default="CHR"
        Column name for chromosome information.
    pos : str, default="POS"
        Column name for position information.
    ref : str, default="NEA"
        Column name for reference allele (non-effect allele).
    alt : str, default="EA"
        Column name for alternative allele (effect allele).
    status : str, default="STATUS"
        Column name for status codes.
    threads : int, default=1
        Number of CPU cores to use for parallel processing.
    chunksize : int, default=5000000
        Size of chunks for processing large reference files.
    ref_snpid : str, default="SNPID"
        Column name for SNP IDs in the reference TSV file.
    ref_rsid : str, default="rsID"
        Column name for rsIDs in the reference TSV file.
    overwrite : str, default="empty"
        Overwrite mode for rsID assignment:
        - "all": overwrite rsID for all available rsID
        - "invalid": only assign rsID for variants with invalid rsID
        - "empty": only assign rsID for variants with NA rsID
    verbose : bool, default=True
        If True, print progress messages.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object for recording process information.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper.

    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with rsID assignments.

    Notes
    -----
    The function first checks if the required columns are present in the summary
    statistics. For VCF reference files, it matches variants based on
    CHR:POS:REF:ALT/ALT:REF. For TSV reference files, it matches based on SNPID.
    The function handles parallel processing for large datasets and provides
    detailed logging of the assignment process.
    """
    # Handle both DataFrame and Sumstats object
    import pandas as pd
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
        sumstats_obj = None
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    # Get mapper from Sumstats object if available
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chr in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chr])

    if ref_mode=="vcf":
        # Auto-detect reference format from VCF file
        mapper.detect_reference_format(path)

        log.write(" -Assigning rsID based on CHR:POS and REF:ALT/ALT:REF...",verbose=verbose)
        ##############################################
        if rsid not in sumstats.columns:
            sumstats[rsid]=pd.Series(dtype="string")

        ###############################################
        total_number= len(sumstats)
        pre_number = (~sumstats[rsid].isna()).sum()

        ##################################################################################################################
        # Match: digit 4 is 0 and digit 5 is 0-4 - Optimized: extract digits once
        digit_4 = _extract_status_digit(sumstats[status], 4)
        digit_5 = _extract_status_digit(sumstats[status], 5)
        standardized_normalized = (digit_4 == 0) & (digit_5.isin([0,1,2,3,4]))
        if overwrite=="all":
            to_assign = standardized_normalized
        if overwrite=="invalid":
            to_assign = (~sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)) & standardized_normalized
        if overwrite=="empty":
            to_assign = sumstats[rsid].isna()& standardized_normalized
        ##################################################################################################################
        # multicore arrangement

        if to_assign.sum()>0:
            if to_assign.sum()<10000: threads=1
            #df_split = np.array_split(sumstats.loc[to_assign, [chr,pos,ref,alt]], threads)
            df_split = _df_split(sumstats.loc[to_assign, [chr,pos,ref,alt]], threads)
            with Pool(threads) as pool:
                map_func = partial(assign_rsid_single,path=path,chr=chr,pos=pos,ref=ref,alt=alt,mapper=mapper) 
                assigned_rsid = pd.concat(pool.map(map_func,df_split))
                sumstats.loc[to_assign,rsid] = assigned_rsid.values
        gc.collect()
        ##################################################################################################################

        after_number = (~sumstats[rsid].isna()).sum()
        log.write(" -rsID annotation for "+str(total_number - after_number) +" variants need to be fixed!",verbose=verbose)
        log.write(" -Annotated "+str(after_number - pre_number) +" rsID successfully!",verbose=verbose)
    
    ##################################################################################################################
    elif ref_mode=="tsv":
        
        #standardized_normalized = sumstats[status].str.match("\w\w\w[0][01234]\w\w", case=False, flags=0, na=False)
        standardized_normalized = sumstats[status] == sumstats[status]

        if rsid not in sumstats.columns:
            sumstats[rsid]=pd.Series(dtype="string")
            
        if overwrite == "empty":
            to_assign = sumstats[rsid].isna() & standardized_normalized
        if overwrite=="all":
            to_assign = standardized_normalized
        if overwrite=="invalid":
            to_assign = (~sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)) & standardized_normalized
            
        total_number= len(sumstats)
        pre_number = (~sumstats[rsid].isna()).sum()
        log.write(" -"+str(to_assign.sum()) +" rsID could be possibly fixed...", verbose=verbose)
        if to_assign.sum()>0:
            sumstats = sumstats.set_index(snpid)  
            dic_chuncks = pd.read_csv(path,sep="\t",usecols=[ref_snpid,ref_rsid],
                              chunksize=chunksize,index_col=ref_snpid,
                              dtype={ref_snpid:"string",ref_rsid:"string"})

            log.write(" -Setting block size: ",chunksize,verbose=verbose)
            log.write(" -Loading block: ",end="",verbose=verbose)     
            for i,dic in enumerate(dic_chuncks):
                gc.collect()
                log.write(i," ",end=" ",show_time=False)  
                dic = dic.rename(index={ref_snpid:snpid})
                dic = dic.rename(columns={ref_rsid:rsid})  
                dic = dic.loc[~dic.index.duplicated(keep=False),:]
                sumstats.update(dic,overwrite=True)

            log.write("\n",end="",show_time=False,verbose=verbose) 
            sumstats = sumstats.reset_index()
            sumstats = sumstats.rename(columns = {'index':snpid})

            after_number = (~sumstats[rsid].isna()).sum()
            log.write(" -rsID annotation for "+str(total_number - after_number) +" variants needed to be fixed!",verbose=verbose)
            log.write(" -Annotated "+str(after_number - pre_number) +" rsID successfully!",verbose=verbose)
        else:
            log.write(" -No rsID can be fixed...skipping...",verbose=verbose)
        ################################################################################################################
    
    # Update harmonization status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _append_meta_record, _update_harmonize_step
            if ref_mode == "tsv":
                sumstats_obj.meta["gwaslab"]["references"]["ref_rsid_tsv"] = path
            elif ref_mode == "vcf":
                sumstats_obj.meta["gwaslab"]["references"]["ref_rsid_vcf"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_rsid_vcf"], path)
            assign_rsid_kwargs = {
                'path': path, 'ref_mode': ref_mode, 'snpid': snpid, 'rsid': rsid, 'chr': chr, 'pos': pos,
                'ref': ref, 'alt': alt, 'status': status, 'threads': threads, 'chunksize': chunksize,
                'overwrite': overwrite
            }
            _update_harmonize_step(sumstats_obj, "assign_rsid", assign_rsid_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats
#################################################################################################################################################
#single record assignment

def check_strand_status(
    chr: Union[str, int],
    start: int,
    end: int,
    ref: str,
    alt: str,
    eaf: float,
    vcf_reader: VariantFile,
    alt_freq: str,
    ref_maf_threshold: float,
    status: Union[int, str],
    mapper: Optional['ChromosomeMapper'] = None
) -> Union[int, str]:
    ### 0 : not palindromic
    ### 1 : palindromic +strand 
    ### 2 : palindromic -strand -> need to flip -> flipped
    ### 5 : palindromic -strand -> need to flip
    ### 8 : no ref data
    if mapper is not None:
        # Convert sumstats chr to reference format for VCF lookup
        # as_string=True for VCF compatibility (pysam requires string)
        chr = mapper.sumstats_to_reference(chr, reference_file=None, as_string=True)
    # Convert status to integer
    status_val = int(status) if isinstance(status, (int, np.integer)) else int(str(status))
    
    try:
        chr_seq = vcf_reader.fetch(chr,start,end)
    except:
        return set_status_digit(status_val, 7, 8)  # Set digit 7 to 8
        
    
    for record in chr_seq:
        if record.pos==end and record.ref==ref and (alt in record.alts):
            if min(record.info[alt_freq][0], 1- record.info[alt_freq][0]) > ref_maf_threshold + FREQ_COMPARISON_EPSILON:
                return set_status_digit(status_val, 7, 8)  # Set digit 7 to 8
            # After check_ref() and flip_allele_stats(), EA should match ALT
            # So we compare EAF with ALT AF directly
            if  (record.info[alt_freq][0]<0.5) and (eaf<0.5):
                return set_status_digit(status_val, 7, 1)  # Set digit 7 to 1
            elif (record.info[alt_freq][0]>0.5) and (eaf>0.5):
                return set_status_digit(status_val, 7, 1)  # Set digit 7 to 1
            else:
                return set_status_digit(status_val, 7, 5)  # Set digit 7 to 5
    return set_status_digit(status_val, 7, 8)  # Set digit 7 to 8

def check_strand_status_cache(
    data: Union[pd.DataFrame, np.ndarray],
    cache: Any,
    ref_infer: Optional[str] = None,
    ref_alt_freq: Optional[str] = None,
    ref_maf_threshold: float = 0.4,
    mapper: Optional['ChromosomeMapper'] = None,
    trust_cache: bool = True,
    log: Log = Log(),
    verbose: bool = True
) -> List[Union[int, str]]:
    if not trust_cache:
        assert ref_infer is not None, "If trust_cache is False, ref_infer must be provided"
        log.warning("You are not trusting the cache, this will slow down the process. Please consider building a complete cache.")

    if ref_infer is not None and not trust_cache:
        vcf_reader = VariantFile(ref_infer)

    if isinstance(data, pd.DataFrame):
        data = data.values
    
    in_cache = 0
    new_statuses = []
    
    for i in range(data.shape[0]):
        _chrom, pos, ref, alt, eaf, status = data[i]
        chrom = _chrom
        start = pos - 1
        end = pos
        
        if mapper is not None:
            # as_string=True for VCF compatibility (pysam requires string)
            chrom = mapper.sumstats_to_reference(chrom, reference_file=ref_infer if ref_infer else None, as_string=True)
        
        # Convert status to integer
        status_val = int(status) if isinstance(status, (int, np.integer)) else int(str(status))
        
        new_status = set_status_digit(status_val, 7, 8)  # default value: set digit 7 to 8
        
        cache_key = f"{chrom}:{pos}:{ref}:{alt}"
        if cache_key in cache:
            in_cache += 1
            record = cache[cache_key]


            if record is None:
                new_status = set_status_digit(status_val, 7, 8)
            else:
                if min(record, 1- record) > ref_maf_threshold + FREQ_COMPARISON_EPSILON:
                    new_status = set_status_digit(status_val, 7, 8)
                elif (record<0.5) and (eaf<0.5):
                    new_status = set_status_digit(status_val, 7, 1)
                elif (record>0.5) and (eaf>0.5):
                    new_status = set_status_digit(status_val, 7, 1)
                else:
                    new_status = set_status_digit(status_val, 7, 5)
        else:
            if not trust_cache:
                # If we don't trust the cache as a not complete cache, we should perform the check reading from the VCF file
                new_status = check_strand_status(_chrom, start, end, ref, alt, eaf, vcf_reader, ref_alt_freq, ref_maf_threshold, status, mapper)
        
        new_statuses.append(new_status)
        
    log.write(f"  -Elements in cache: {in_cache}", verbose=verbose)
    return new_statuses


def check_unkonwn_indel(
    chr: Union[str, int],
    start: int,
    end: int,
    ref: str,
    alt: str,
    eaf: float,
    vcf_reader: VariantFile,
    alt_freq: str,
    ref_maf_threshold: float,
    status: Union[int, str],
    mapper: Optional['ChromosomeMapper'] = None,
    daf_tolerance: float = 0.2
) -> Union[int, str]:
    ### input : unknown indel, both on genome (xx1[45]x)
    ### 3 no flip
    ### 4 unknown indel,fixed   (6->5)
    ### 6 flip
    
    if mapper is not None:
        # Convert sumstats chr to reference format for VCF lookup
        chr = mapper.sumstats_to_reference(chr, reference_file=None)
    # Convert status to integer
    status_val = int(status) if isinstance(status, (int, np.integer)) else int(str(status))
    
    try:
        chr_seq = vcf_reader.fetch(chr,start,end)
    except:
        return set_status_digit(status_val, 7, 8)

    for record in chr_seq:
        if min(record.info[alt_freq][0], 1- record.info[alt_freq][0]) > ref_maf_threshold + FREQ_COMPARISON_EPSILON:
                return set_status_digit(status_val, 7, 8)
        if record.pos==end and record.ref==ref and (alt in record.alts):
            if  abs(record.info[alt_freq][0] - eaf)<daf_tolerance:
                return set_status_digit(status_val, 7, 3)
   
        elif record.pos==end and record.ref==alt and (ref in record.alts):
            if  abs(record.info[alt_freq][0] - (1 - eaf))<daf_tolerance:
                return set_status_digit(status_val, 7, 6)

    return set_status_digit(status_val, 7, 8)


def check_unkonwn_indel_cache(
    data: Union[pd.DataFrame, np.ndarray],
    cache: Any,
    ref_infer: Optional[str] = None,
    ref_alt_freq: Optional[str] = None,
    ref_maf_threshold: Optional[float] = None,
    mapper: Optional['ChromosomeMapper'] = None,
    daf_tolerance: float = 0.2,
    trust_cache: bool = True,
    log: Log = Log(),
    verbose: bool = True
) -> List[Union[int, str]]:
    if not trust_cache:
        assert ref_infer is not None, "If trust_cache is False, ref_infer must be provided"
        log.warning("You are not trusting the cache, this will slow down the process. Please consider building a complete cache.")

    if ref_infer is not None:
        vcf_reader = VariantFile(ref_infer)

    if isinstance(data, pd.DataFrame):
        data = data.values
    
    in_cache = 0
    new_statuses = []
    
    for i in range(data.shape[0]):
        _chrom, pos, ref, alt, eaf, status = data[i]
        chrom = _chrom
        
        if mapper is not None:
            # as_string=True for VCF compatibility (pysam requires string)
            chrom = mapper.sumstats_to_reference(chrom, reference_file=ref_infer if ref_infer else None, as_string=True)
        start = pos - 1
        end = pos
        
        # Convert status to integer
        status_val = int(status) if isinstance(status, (int, np.integer)) else int(str(status))
        
        new_status = set_status_digit(status_val, 7, 8)  # default value: set digit 7 to 8
        
        cache_key_ref_alt = f"{chrom}:{pos}:{ref}:{alt}"
        cache_key_alt_ref = f"{chrom}:{pos}:{alt}:{ref}"

        if cache_key_ref_alt in cache:
            in_cache += 1
            record = cache[cache_key_ref_alt]

            if record is None:
                new_status = set_status_digit(status_val, 7, 8)
            else:
                if min(record, 1- record) > ref_maf_threshold + FREQ_COMPARISON_EPSILON:
                    return set_status_digit(status_val, 7, 8)
                if  abs(record - eaf)<daf_tolerance:
                    new_status = set_status_digit(status_val, 7, 3)

        elif cache_key_alt_ref in cache:
            in_cache += 1
            record = cache[cache_key_alt_ref]
            if record is None:
                new_status = set_status_digit(status_val, 7, 8)
            else:
                if min(record, 1- record) > ref_maf_threshold + FREQ_COMPARISON_EPSILON:
                    return set_status_digit(status_val, 7, 8)
                if  abs(record - (1 - eaf))<daf_tolerance:
                    new_status = set_status_digit(status_val, 7, 6)

        else:
            if not trust_cache:
                # If we don't trust the cache as a not complete cache, we should perform the check reading from the VCF file
                new_status = check_unkonwn_indel(_chrom, start, end, ref, alt, eaf, vcf_reader, ref_alt_freq, ref_maf_threshold,status, mapper, daf_tolerance)
                
        new_statuses.append(new_status)
        
    log.write(f"  -Elements in cache: {in_cache}", verbose=verbose)
    return new_statuses

                                               
def get_reverse_complementary_allele(a: str) -> str:
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a[::-1].translate(dic)
                                                 
def is_palindromic(sumstats: pd.DataFrame, a1: str = "EA", a2: str = "NEA") -> pd.Series:
    gc= (sumstats[a1]=="G") & (sumstats[a2]=="C")
    cg= (sumstats[a1]=="C") & (sumstats[a2]=="G")
    at= (sumstats[a1]=="A") & (sumstats[a2]=="T")
    ta= (sumstats[a1]=="T") & (sumstats[a2]=="A")
    palindromic = gc | cg | at | ta 
    return palindromic
##################################################################################################################################################
#single df assignment

def check_strand(
    sumstats: pd.DataFrame,
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    ref_maf_threshold: Optional[float] = None,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    mapper: Optional['ChromosomeMapper'] = None,
    status: str = "STATUS"
) -> pd.Series:
    vcf_reader = VariantFile(ref_infer)
    status_part = sumstats.apply(lambda x:check_strand_status(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],x.iloc[4],vcf_reader,ref_alt_freq,ref_maf_threshold,x.iloc[5],mapper),axis=1) 
    return status_part

def check_strand_cache(
    sumstats: pd.DataFrame,
    cache: Any,
    ref_infer: Optional[str] = None,
    ref_alt_freq: Optional[float] = None,
    ref_maf_threshold: Optional[float] = None,
    mapper: Optional['ChromosomeMapper'] = None,
    trust_cache: bool = True,
    log: Log = Log(),
    verbose: bool = True
) -> pd.Series:
    assert cache is not None, "Cache must be provided"
    status_part = check_strand_status_cache(sumstats,cache,ref_infer,ref_alt_freq,ref_maf_threshold,mapper,trust_cache,log,verbose)
    return status_part

def check_indel(
    sumstats: pd.DataFrame,
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    ref_maf_threshold: Optional[float] = None,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    mapper: Optional['ChromosomeMapper'] = None,
    status: str = "STATUS",
    daf_tolerance: float = 0.2
) -> pd.Series:
    vcf_reader = VariantFile(ref_infer)
    status_part = sumstats.apply(lambda x:check_unkonwn_indel(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],x.iloc[4],vcf_reader,ref_alt_freq,ref_maf_threshold,x.iloc[5],mapper,daf_tolerance),axis=1)
    return status_part

def check_indel_cache(
    sumstats: pd.DataFrame,
    cache: Any,
    ref_infer: Optional[str] = None,
    ref_alt_freq: Optional[str] = None,
    ref_maf_threshold: Optional[float] = None,
    mapper: Optional['ChromosomeMapper'] = None,
    daf_tolerance: float = 0.2,
    trust_cache: bool = True,
    log: Log = Log(),
    verbose: bool = True
) -> pd.Series:
    assert cache is not None, "Cache must be provided"
    status_part = check_unkonwn_indel_cache(sumstats,cache,ref_infer,ref_alt_freq,ref_maf_threshold,mapper,daf_tolerance,trust_cache,log,verbose)
    return status_part

##################################################################################################################################################
@with_logging(
    start_to_msg="infer strand for palindromic SNPs/align indistinguishable indels",
    finished_msg="inferring strand for palindromic SNPs/align indistinguishable indels",
    start_cols=["CHR","POS","EA","NEA","STATUS"],
    start_function=".infer_strand()",
    must_kwargs=["ref_alt_freq"]
)
def _parallelize_infer_strand(
    sumstats: Union['Sumstats', pd.DataFrame],
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    maf_threshold: float = 0.40,
    ref_maf_threshold: float = 0.4,
    daf_tolerance: float = 0.20,
    remove_snp: str = "",
    mode: str = "pi",
    threads: int = 1,
    remove_indel: str = "",
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    status: str = "STATUS",
    mapper: Optional['ChromosomeMapper'] = None,
    cache_options: Dict[str, Any] = {},
    verbose: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    """
    Args:
    cache_options : A dictionary with the following keys:
        - cache_manager: CacheManager object or None. If any between cache_loader and cache_process is not None, or use_cache is True, a CacheManager object will be created automatically.
        - trust_cache: bool (optional, default: True). Whether to completely trust the cache or not. Trusting the cache means that any key not found inside the cache will be considered as a missing value even in the VCF file.
        - cache_loader: Object with a get_cache() method or None.
        - cache_process: Object with an apply_fn() method or None.
        - use_cache: bool (optional, default: False). If any of the cache_manager, cache_loader or cache_process is not None, this will be set to True automatically.
                     If set to True and all between cache_manager, cache_loader and cache_process are None, the cache will be loaded (or built) on the spot.

        The usefulness of a cache_loader or cache_process object is to pass a custom object which already has the cache loaded. This can be useful if the cache is loaded in background in another thread/process while other operations are performed.
        The cache_manager is a CacheManager object is used to expose the API to interact with the cache.
    """
    # Handle both DataFrame and Sumstats object
    import pandas as pd
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Get mapper from Sumstats object if available
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chr in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chr])
    
    # Auto-detect reference format from VCF file
    mapper.detect_reference_format(ref_infer)
    
    # Setup cache variables
    cache_manager = cache_options.get("cache_manager", None)
    if cache_manager is not None:
        assert isinstance(cache_manager, CacheManager), "cache_manager must be a CacheManager object"
    trust_cache = cache_options.get("trust_cache", True)
    cache_loader = cache_options.get("cache_loader", None)
    cache_process = cache_options.get("cache_process", None)
    use_cache = any(c is not None for c in [cache_manager, cache_loader, cache_process]) or cache_options.get('use_cache', False)
    _threads = threads # backup threads
    
    log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)  

    if "p" in mode:
        ## checking \w\w\w\w[0]\w\w -> standardized and normalized snp
        # Match: digit 4 is 0 and digit 5 is 0 - Optimized: extract digits once
        digit_4 = _extract_status_digit(sumstats[status], 4)
        digit_5 = _extract_status_digit(sumstats[status], 5)
        good_chrpos = (digit_4 == 0) & (digit_5 == 0) 
        palindromic = good_chrpos & is_palindromic(sumstats[[ref,alt]],a1=ref,a2=alt)   
        not_palindromic_snp = good_chrpos & (~palindromic)

        ##not palindromic : change status
        sumstats.loc[not_palindromic_snp,status] = vchange_status(sumstats.loc[not_palindromic_snp,status], 7 ,"9","0")  
        log.write(" -Identified ", palindromic.sum(), " palindromic SNPs...", verbose=verbose)
        
        #palindromic but can not infer
        
        maf_can_infer   = (sumstats[eaf] <= maf_threshold + FREQ_COMPARISON_EPSILON) | (sumstats[eaf] >= 1 - maf_threshold - FREQ_COMPARISON_EPSILON)
        
        sumstats.loc[palindromic&(~maf_can_infer),status] = vchange_status(sumstats.loc[palindromic&(~maf_can_infer),status],7,"9","7")
        
        #palindromic WITH UNKNWON OR UNCHECKED STATUS
        # Match: digit 6 is 0,1,2 and digit 7 is 8,9 - Optimized: extract digits once
        digit_6 = _extract_status_digit(sumstats[status], 6)
        digit_7 = _extract_status_digit(sumstats[status], 7)
        unknow_palindromic = (digit_6.isin([0,1,2])) & (digit_7.isin([8,9])) 

        unknow_palindromic_to_check = palindromic & maf_can_infer & unknow_palindromic
        
        log.write(" -After filtering by MAF<= {} , {} palindromic SNPs with unknown strand will be inferred...".format(maf_threshold, unknow_palindromic_to_check.sum()), verbose=verbose)

        ######################################################################################### 
        if unknow_palindromic_to_check.sum()>0:
            if unknow_palindromic_to_check.sum()<10000:
                threads=1

            if use_cache and cache_manager is None:
                cache_manager = CacheManager(base_path=ref_infer, cache_loader=cache_loader, cache_process=cache_process,
                                             ref_alt_freq=ref_alt_freq, category=PALINDROMIC_INDEL,
                                             threads=_threads, log=log, verbose=verbose)

            log.write(" -Starting strand inference for palindromic SNPs...",verbose=verbose)
            df_to_check = sumstats.loc[unknow_palindromic_to_check,[chr,pos,ref,alt,eaf,status]]
           
            if use_cache and cache_manager.cache_len > 0:
                log.write("  -Using cache for strand inference",verbose=verbose)
                status_inferred = cache_manager.apply_fn(check_strand_cache, sumstats=df_to_check, ref_infer=ref_infer, ref_alt_freq=ref_alt_freq, ref_maf_threshold=ref_maf_threshold, mapper=mapper, trust_cache=trust_cache, log=log, verbose=verbose)
                sumstats.loc[unknow_palindromic_to_check,status] = status_inferred
            else:
                #df_split = np.array_split(df_to_check, threads)
                df_split = _df_split(df_to_check, threads)
                with Pool(threads) as pool:
                    map_func = partial(check_strand,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,status=status,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,ref_maf_threshold=ref_maf_threshold,mapper=mapper) 
                    status_inferred = pd.concat(pool.map(map_func,df_split))
                    sumstats.loc[unknow_palindromic_to_check,status] = status_inferred.values
            log.write(" -Finished strand inference.",verbose=verbose)
        else:
            log.warning("No palindromic variants available for checking.")
        #########################################################################################
        #0 Not palindromic SNPs
        #1 Palindromic +strand  -> no need to flip
        #2 palindromic -strand  -> need to flip -> fixed
        #3 Indel no need flip
        #4 Unknown Indel -> fixed
        #5 Palindromic -strand -> need to flip
        #6 Indel need flip
        #7 indistinguishable
        #8 Not matching or No information
        #9 Unchecked

        # Match digit 7 - Optimized: extract digit once, then use vectorized comparisons
        digit_7 = _extract_status_digit(sumstats[status], 7)
        status0 = (digit_7 == 0)
        status1 = (digit_7 == 1)
        status5 = (digit_7 == 5)
        status7 = (digit_7 == 7)
        status8 = (digit_7 == 8)  

        log.write("  -Non-palindromic : ",sum(status0),verbose=verbose)
        log.write("  -Palindromic SNPs on + strand: ",sum(status1),verbose=verbose)
        log.write("  -Palindromic SNPs on - strand and needed to be flipped:",sum(status5),verbose=verbose)   
        log.write("  -Palindromic SNPs with MAF not available to infer : ",sum(status7),verbose=verbose)  
        log.write("  -Palindromic SNPs lacking reference matches, or with reference ALT unavailable or above MAF threshold for inference : ",sum(status8),verbose=verbose)  

        if ("7" in remove_snp) and ("8" in remove_snp) :
            log.write("  -Palindromic SNPs with MAF not available to infer and with no matches or no information will be removed",verbose=verbose) 
            sumstats = sumstats.loc[~(status7 | status8),:].copy()
        elif "8" in remove_snp:
            log.write("  -Palindromic SNPs with no matches or no information will be removed",verbose=verbose)
            sumstats = sumstats.loc[~status8,:].copy()
        elif "7" in remove_snp:
            log.write("  -Palindromic SNPs with MAF not available to infer will be removed",verbose=verbose) 
            sumstats = sumstats.loc[~status7,:].copy()

    ### unknow_indel
    if "i" in mode:
        # Match: digit 6 is 6 and digit 7 is 8,9 - Optimized: extract digits once
        digit_6 = _extract_status_digit(sumstats[status], 6)
        digit_7 = _extract_status_digit(sumstats[status], 7)
        unknow_indel = (digit_6 == 6) & (digit_7.isin([8,9]))   
        log.write(" -Identified ", unknow_indel.sum(), " indistinguishable Indels...", verbose=verbose)
        if unknow_indel.sum()>0:
            log.write(" -Indistinguishable indels will be inferred from reference vcf REF and ALT...",verbose=verbose)
            #########################################################################################  
            #with maf can not infer
            #maf_can_infer   = (sumstats[eaf] <= maf_threshold + FREQ_COMPARISON_EPSILON) | (sumstats[eaf] >= 1 - maf_threshold - FREQ_COMPARISON_EPSILON) 
            #sumstats.loc[unknow_indel&(~maf_can_infer),status] = vchange_status(sumstats.loc[unknow_indel&(~maf_can_infer),status],7,"9","8") 
            log.write(" -Difference in allele frequency (DAF) tolerance: {}".format(daf_tolerance),verbose=verbose)
                         
            if unknow_indel.sum()>0:
                if unknow_indel.sum()<10000:
                    threads=1

                if use_cache and cache_manager is None:
                    cache_manager = CacheManager(base_path=ref_infer, cache_loader=cache_loader, cache_process=cache_process,
                                                ref_alt_freq=ref_alt_freq, category=PALINDROMIC_INDEL,
                                                threads=_threads, log=log, verbose=verbose)

                log.write(" -Starting indistinguishable indel inference...",verbose=verbose)
                df_to_check = sumstats.loc[unknow_indel,[chr,pos,ref,alt,eaf,status]]
            
                if use_cache and cache_manager.cache_len > 0:
                    log.write("  -Using cache for indel inference",verbose=verbose)
                    status_inferred = cache_manager.apply_fn(check_indel_cache, sumstats=df_to_check, ref_infer=ref_infer, ref_alt_freq=ref_alt_freq,ref_maf_threshold=ref_maf_threshold, mapper=mapper, daf_tolerance=daf_tolerance, trust_cache=trust_cache, log=log, verbose=verbose)
                    sumstats.loc[unknow_indel,status] = status_inferred
                else:
                    #df_split = np.array_split(sumstats.loc[unknow_indel, [chr,pos,ref,alt,eaf,status]], threads)
                    df_split = _df_split(sumstats.loc[unknow_indel, [chr,pos,ref,alt,eaf,status]], threads)
                    with Pool(threads) as pool:
                        map_func = partial(check_indel,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,status=status,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,ref_maf_threshold=ref_maf_threshold,mapper=mapper,daf_tolerance=daf_tolerance) 
                        status_inferred = pd.concat(pool.map(map_func,df_split))
                        sumstats.loc[unknow_indel,status] = status_inferred.values
                log.write(" -Finished indistinguishable indel inference.",verbose=verbose)

            #########################################################################################

            # Optimized: extract digits once, then use vectorized comparisons
            digit_6 = _extract_status_digit(sumstats[status], 6)
            digit_7 = _extract_status_digit(sumstats[status], 7)
            status3 = (digit_7 == 3)
            status6 = (digit_7 == 6)
            status8 = (digit_6 == 6) & (digit_7 == 8)  

            log.write("  -Indels ea/nea match reference : ", status3.sum(), verbose=verbose)
            log.write("  -Indels ea/nea need to be flipped : ", status6.sum(), verbose=verbose)
            log.write("  -Indels with no matches or no information : ", status8.sum(), verbose=verbose)
            if "8" in remove_indel:
                log.write("  -Indels with no matches or no information will be removed",verbose=verbose)
                sumstats = sumstats.loc[~status8,:].copy()    
        else:
            log.warning("No indistinguishable indels available for checking.") 
    
    # Update harmonization status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            infer_strand_kwargs = {
                'ref_infer': ref_infer, 'ref_alt_freq': ref_alt_freq, 'maf_threshold': maf_threshold,
                'ref_maf_threshold': ref_maf_threshold, 'daf_tolerance': daf_tolerance,
                'remove_snp': remove_snp, 'mode': mode, 'threads': threads, 'remove_indel': remove_indel
            }
            _update_harmonize_step(sumstats_obj, "infer_strand", infer_strand_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats




















################################################################################################################

@with_logging(
    start_to_msg="check the difference between EAF (sumstats) and ALT frequency (reference VCF)",
    finished_msg="checking the difference between EAF (sumstats) and ALT frequency (reference VCF)",
    start_cols=["CHR","POS","EA","NEA","EAF","STATUS"],
    start_function=".check_daf()",
    must_kwargs=["ref_alt_freq"]
)
def _parallelize_check_af(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    maf_threshold: float = 0.4,
    column_name: str = "DAF",
    suffix: str = "",
    threads: int = 1,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    status: str = "STATUS",
    mapper: Optional['ChromosomeMapper'] = None,
    force: bool = False,
    verbose: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    """
    Check the difference between effect allele frequency (EAF) in summary statistics and alternative allele frequency in reference VCF.

    This function calculates the difference between the effect allele frequency (EAF) in the summary 
    statistics and the alternative allele frequency (ALT_AF) in the reference VCF file. The difference 
    (DAF) is stored in a new column in the summary statistics DataFrame.

    **Purpose:**
    - Quality control: Identify variants with large differences in allele frequency between sumstats 
      and reference, which may indicate:
      - Population differences
      - Allele mismatches or strand flips
      - Data quality issues
    - Validation: Verify that EAF values in sumstats are consistent with reference population frequencies.

    **DAF Calculation:**
    - DAF = EAF (sumstats) - ALT_AF (reference VCF)
    - Positive DAF: EAF in sumstats is higher than reference
    - Negative DAF: EAF in sumstats is lower than reference
    - Large |DAF| values (> 0.2) may indicate issues requiring investigation

    **When to use:**
    - After `infer_af()` to validate inferred EAF values
    - Before harmonization to identify potential allele mismatches
    - For quality control to flag variants with unusual frequency differences

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame containing variant data. Must have columns: CHR, POS, EA, NEA, EAF.
    ref_infer : str
        Path to the reference VCF file. Must be indexed (tabix) for efficient querying.
    ref_alt_freq : str, optional
        Field name for alternative allele frequency in VCF INFO section (e.g., "AF", "AF_popmax", 
        "gnomAD_AF"). If None, the function will attempt to auto-detect common AF field names.
    maf_threshold : float, default=0.4
        Minor allele frequency threshold for filtering variants. Only variants with MAF ≤ threshold 
        in both sumstats and reference are included in DAF statistics. This helps focus on common 
        variants where frequency differences are more meaningful.
    column_name : str, default="DAF"
        Name of the column to store the difference values. The final column name will be 
        `column_name + suffix`.
    suffix : str, default=""
        Suffix to append to the column name (e.g., "_pop1", "_gnomad"). Useful when comparing 
        multiple reference populations.
    threads : int, default=1
        Number of CPU cores to use for parallel processing. Set to 1 if processing < 10,000 variants.
    chr : str, default="CHR"
        Column name for chromosome in sumstats.
    pos : str, default="POS"
        Column name for position in sumstats.
    ref : str, default="NEA"
        Column name for reference/non-effect allele in sumstats.
    alt : str, default="EA"
        Column name for alternative/effect allele in sumstats.
    eaf : str, default="EAF"
        Column name for effect allele frequency in sumstats.
    status : str, default="STATUS"
        Column name for status codes. By default, only processes variants with STATUS digit 4 = 0 
        (standardized and normalized), unless `force=True`.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper with automatic format detection.
    force : bool, default=False
        If True, check all variants regardless of STATUS codes. If False, only processes variants 
        with valid harmonization status (STATUS digit 4 = 0).
    verbose : bool, default=True
        If True, print progress messages and DAF statistics (max, min, mean, std, abs statistics).
    log : gwaslab.g_Log.Log, default=Log()
        Logging object for recording process information.

    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with a new column (`column_name + suffix`) containing 
        the difference between EAF (sumstats) and ALT_AF (reference VCF). Variants that could not 
        be matched in the reference VCF will have DAF = NaN.

    Notes
    -----
    - The difference in allele frequency (DAF) is calculated as: DAF = EAF (sumstats) - ALT_AF (reference VCF)
    - **Important**: This DAF is NOT the derived allele frequency. It is simply the difference in 
      allele frequencies between two datasets.
    - The function requires the reference VCF to be indexed (tabix) for efficient chromosome-based 
      querying. Use `tabix -p vcf reference.vcf.gz` to create an index.
    - By default, only processes variants with valid harmonization status (STATUS digit 4 = 0) to 
      ensure alleles are standardized and normalized.
    - The function provides comprehensive statistics about DAF values including:
      - Maximum and minimum DAF
      - Mean and standard deviation of DAF
      - Statistics for absolute DAF values
      - Count of variants with |DAF| > 0.2 (potential issues)
    - Only variants with valid chromosome position and allele information are checked by default.
    - For small datasets (< 10,000 variants), the function automatically sets `threads=1` to avoid overhead.
    - Large |DAF| values (> 0.2) may indicate:
      - Population differences (expected for population-specific variants)
      - Allele mismatches (check EA/NEA alignment)
      - Strand flips (use `infer_strand()` to resolve)
      - Data quality issues (verify EAF calculation in sumstats)

    See Also
    --------
    infer_af : Infer EAF from reference VCF before checking differences.
    plot_daf : Visualize DAF distribution to identify outliers.
    infer_strand : Resolve strand orientation issues that may cause large DAF values.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
        is_dataframe = True
        sumstats_obj = None
    else:
        sumstats_obj = sumstats_or_dataframe
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Get mapper from Sumstats object if available
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chr in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chr])
    
    # Auto-detect reference format from VCF file
    mapper.detect_reference_format(ref_infer)

    column_name = column_name + suffix
    


    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)  
        if not force:
            # Match: digit 4 is 0 - Optimized: extract digit once
            digit_4 = _extract_status_digit(sumstats[status], 4)
            good_chrpos = (digit_4 == 0)  
        else:
            good_chrpos = pd.Series([True] * len(sumstats), index=sumstats.index)
        log.write(" -Checking variants:", good_chrpos.sum(), verbose=verbose)
        sumstats[column_name]=np.nan
    
    ########################  
        if (~sumstats[eaf].isna()).sum()<10000:
            threads=1       
        #df_split = np.array_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]], threads)
        df_split = _df_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]], threads)
        with Pool(threads) as pool:
            if (~sumstats[eaf].isna()).sum()>0:
                map_func = partial(checkaf,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,column_name=column_name,mapper=mapper) 
                result = pd.concat(pool.map(map_func,df_split))
                # Cast to match original column dtype to avoid FutureWarning
                if column_name in sumstats.columns:
                    original_dtype = sumstats[column_name].dtype
                    result[column_name] = result[column_name].astype(original_dtype)
                sumstats.loc[good_chrpos,[column_name]] = result[column_name]
    ###########################
        #status_inferred = sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]].apply(lambda x:check_daf(x[0],x[1]-1,x[1],x[2],x[3],x[4],vcf_reader,ref_alt_freq,chr_dict),axis=1)
        log.write(" -Difference in allele frequency (DAF) = EAF (sumstats) - ALT_AF (reference VCF)", verbose=verbose)  
        log.write(" -Note: this DAF is not the derived allele frequency.", verbose=verbose)  
        #sumstats.loc[good_chrpos,"DAF"] = status_inferred.values
        #sumstats["DAF"]=sumstats["DAF"].astype("float")     
        log.write(" - {} max:".format(column_name), np.nanmax(sumstats[column_name]),verbose=verbose)
        log.write(" - {} min:".format(column_name), np.nanmin(sumstats[column_name]),verbose=verbose)
        log.write(" - {} sd:".format(column_name), np.nanstd(sumstats[column_name]),verbose=verbose)
        log.write(" - abs({}) min:".format(column_name), np.nanmin(np.abs(sumstats[column_name])),verbose=verbose) 
        log.write(" - abs({}) max:".format(column_name), np.nanmax(np.abs(sumstats[column_name])),verbose=verbose)
        log.write(" - abs({}) sd:".format(column_name), np.nanstd(np.abs(sumstats[column_name])),verbose=verbose) 
        log.write("Finished allele frequency checking!") 
    
    # Set metadata if Sumstats object is available
    try:
        sumstats_obj = getattr(log, '_sumstats_obj', None)
        if sumstats_obj is not None:
            from gwaslab.info.g_meta import _append_meta_record
            sumstats_obj.meta["gwaslab"]["references"]["ref_infer_daf"] = _append_meta_record(
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_daf"], ref_infer)
    except:
        pass
    
    return sumstats

def checkaf(
    sumstats: pd.DataFrame,
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    column_name: str = "DAF",
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    mapper: Optional['ChromosomeMapper'] = None
) -> pd.DataFrame:
    #vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    vcf_reader = VariantFile(ref_infer)
    # Auto-detect reference format from VCF file if mapper provided
    if mapper is not None:
        mapper.detect_reference_format(ref_infer)
    def afapply(x,vcf,alt_freq,mapper):
            return check_daf(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],x.iloc[4],vcf_reader,ref_alt_freq,mapper)
    map_func = partial(afapply,vcf=vcf_reader,alt_freq=ref_alt_freq,mapper=mapper)
    status_inferred = sumstats.apply(map_func,axis=1)
    sumstats[column_name] = status_inferred.values
    sumstats[column_name]=sumstats[column_name].astype("float") 
    return sumstats

def check_daf(
    chr: Union[str, int],
    start: int,
    end: int,
    ref: str,
    alt: str,
    eaf: float,
    vcf_reader: VariantFile,
    alt_freq: str,
    mapper: Optional['ChromosomeMapper'] = None
) -> float:
    if mapper is not None:
        # Convert sumstats chr to reference format for VCF lookup
        chr = mapper.sumstats_to_reference(chr, reference_file=None)
    chr_seq = vcf_reader.fetch(chr,start,end)
    
    for record in chr_seq:
        if record.pos==end:
            if record.ref==ref and (alt in record.alts):
                return eaf - record.info[alt_freq][0]
    return np.nan
################################################################################################################

@with_logging(
    start_to_msg="infer sumstats EAF using reference VCF ALT frequency",
    finished_msg="inferring sumstats EAF using reference VCF ALT frequency",
    start_cols=["CHR","POS","EA","NEA","STATUS"],
    start_function=".check_daf()",
    must_kwargs=["ref_alt_freq"]
)
def _parallelize_infer_af(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    threads: int = 1,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    status: str = "STATUS",
    mapper: Optional['ChromosomeMapper'] = None,
    force: bool = False,
    verbose: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    """
    Infer effect allele frequency (EAF) in summary statistics using reference VCF ALT frequency.

    This function infers the effect allele frequency (EAF) in the summary statistics by matching
    variants with a reference VCF file. It extracts the alternative allele frequency (ALT_AF) from 
    the reference VCF INFO field and updates the EAF values in the summary statistics DataFrame.

    **Workflow:**
    1. Matches variants in sumstats with reference VCF by CHR:POS:EA:NEA.
    2. Extracts ALT frequency from VCF INFO field (specified by `ref_alt_freq`).
    3. Handles allele matching: If EA matches ALT in VCF, uses ALT_AF directly. If EA matches REF 
       in VCF, uses 1 - ALT_AF.
    4. Updates EAF column in sumstats with inferred values.

    **When to use:**
    - When sumstats are missing EAF values but have valid CHR, POS, EA, NEA.
    - When you want to fill in EAF from a reference population (e.g., 1000 Genomes, gnomAD).
    - As a preprocessing step before strand inference or allele frequency comparison.

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame containing variant data. Must have columns: CHR, POS, EA, NEA.
    ref_infer : str
        Path to the reference VCF file. Must be indexed (tabix) for efficient querying.
    ref_alt_freq : str, optional
        Field name for alternative allele frequency in VCF INFO section (e.g., "AF", "AF_popmax", 
        "gnomAD_AF"). If None, the function will attempt to auto-detect common AF field names.
    threads : int, default=1
        Number of CPU cores to use for parallel processing. Set to 1 if processing < 10,000 variants.
    chr : str, default="CHR"
        Column name for chromosome in sumstats.
    pos : str, default="POS"
        Column name for position in sumstats.
    ref : str, default="NEA"
        Column name for reference/non-effect allele in sumstats.
    alt : str, default="EA"
        Column name for alternative/effect allele in sumstats.
    eaf : str, default="EAF"
        Column name for effect allele frequency. This column will be created if it doesn't exist, 
        and existing values will be updated where inference is successful.
    status : str, default="STATUS"
        Column name for status codes. By default, only processes variants with STATUS digit 4 = 0 
        (standardized and normalized), unless `force=True`.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper with automatic format detection.
    force : bool, default=False
        If True, infer EAF for all variants regardless of STATUS codes. If False, only processes 
        variants with valid harmonization status (STATUS digit 4 = 0).
    verbose : bool, default=True
        If True, print progress messages and statistics about inference success rate.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object for recording process information.

    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with inferred EAF values. Variants that could not be 
        matched in the reference VCF will have EAF = NaN.

    Notes
    -----
    - The function requires the reference VCF to be indexed (tabix) for efficient chromosome-based 
      querying. Use `tabix -p vcf reference.vcf.gz` to create an index.
    - By default, only processes variants with valid harmonization status (STATUS digit 4 = 0) to 
      ensure alleles are standardized and normalized.
    - The function uses parallel processing to improve performance on large datasets. For small 
      datasets (< 10,000 variants), it automatically sets `threads=1` to avoid overhead.
    - After inference, the function reports statistics about:
      - Number of variants with EAF successfully inferred
      - Number of variants still missing EAF (not found in reference or missing ref_alt_freq field)
    - The inferred EAF values are stored in the specified EAF column, overwriting existing values 
      where inference is successful.
    - The function automatically handles chromosome name mapping using ChromosomeMapper. 
      for proper matching.

    See Also
    --------
    check_af : Calculate difference between EAF (sumstats) and ALT_AF (reference) after inference.
    infer_eaf_from_maf : Infer EAF from MAF using reference VCF (alternative method).
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
        is_dataframe = True
        sumstats_obj = None
    else:
        sumstats_obj = sumstats_or_dataframe
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Get mapper from Sumstats object if available
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chr in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chr])
    
    # Auto-detect reference format from VCF file
    mapper.detect_reference_format(ref_infer)
    
    if eaf not in sumstats.columns:
        sumstats[eaf]=np.nan
    
    prenumber = sumstats[eaf].isna().sum()

    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)   
        if not force:
            # Match: digit 4 is 0 - Optimized: extract digit once
            digit_4 = _extract_status_digit(sumstats[status], 4)
            good_chrpos = (digit_4 == 0)  
        else:
            good_chrpos = pd.Series([True] * len(sumstats), index=sumstats.index)
        log.write(" -Checking variants:", good_chrpos.sum(), verbose=verbose) 
    
    ########################  
        if sumstats[eaf].isna().sum()<10000:
            threads=1       
        #df_split = np.array_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], threads)
        df_split = _df_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], threads)
        with Pool(threads) as pool:
            map_func = partial(inferaf,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,mapper=mapper) 
            result = pd.concat(pool.map(map_func,df_split))
            # Cast to match original column dtype to avoid FutureWarning
            original_dtype = sumstats[eaf].dtype
            result[eaf] = result[eaf].astype(original_dtype)
            sumstats.loc[good_chrpos,[eaf]] = result[eaf]
    ###########################
        
        afternumber = sumstats[eaf].isna().sum()
        log.write(" -Inferred EAF for {} variants.".format(prenumber - afternumber),verbose=verbose) 
        log.write(" -EAF is still missing for {} variants.".format(afternumber),verbose=verbose) 
    
    # Set metadata if Sumstats object is available
    try:
        sumstats_obj = getattr(log, '_sumstats_obj', None)
        if sumstats_obj is not None:
            from gwaslab.info.g_meta import _append_meta_record
            sumstats_obj.meta["gwaslab"]["references"]["ref_infer_af"] = _append_meta_record(
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_af"], ref_infer)
    except:
        pass
    
    return sumstats

def inferaf(
    sumstats: pd.DataFrame,
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    mapper: Optional['ChromosomeMapper'] = None
) -> pd.DataFrame:
    #vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    vcf_reader = VariantFile(ref_infer)
    # Auto-detect reference format from VCF file if mapper provided
    if mapper is not None:
        mapper.detect_reference_format(ref_infer)
    def afapply(x,vcf,alt_freq,mapper):
            return infer_af(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],vcf_reader,ref_alt_freq,mapper)
    map_func = partial(afapply,vcf=vcf_reader,alt_freq=ref_alt_freq,mapper=mapper)
    status_inferred = sumstats.apply(map_func,axis=1)
    sumstats[eaf] = status_inferred.values
    sumstats[eaf]=sumstats[eaf].astype("float") 
    return sumstats

def infer_af(
    chr: Union[str, int],
    start: int,
    end: int,
    ref: str,
    alt: str,
    vcf_reader: VariantFile,
    alt_freq: str,
    mapper: Optional['ChromosomeMapper'] = None
) -> float:
    if mapper is not None:
        # Convert sumstats chr to reference format for VCF lookup
        chr = mapper.sumstats_to_reference(chr, reference_file=None)
    chr_seq = vcf_reader.fetch(chr,start,end)
    
    for record in chr_seq:
        if record.pos==end:
            if record.ref==ref and (alt in record.alts):
                return record.info[alt_freq][0]
            elif record.ref==alt and (ref in record.alts):
                return 1 - record.info[alt_freq][0]
    return np.nan
##############################################################################################################################################################################################

################################################################################################################
@with_logging(
    start_to_msg="infer sumstats EAF from sumstats MAF using reference VCF ALT frequency",
    finished_msg="inferring sumstats EAF from sumstats MAF using reference VCF ALT frequency",
    start_cols=["CHR","POS","EA","NEA","MAF","STATUS"],
    start_function=".infer_af()",
    must_kwargs=["ref_alt_freq"]
)
def _parallele_infer_af_with_maf(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    threads: int = 1,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    maf: str = "MAF",
    ref_eaf: str = "_REF_EAF",
    status: str = "STATUS",
    mapper: Optional['ChromosomeMapper'] = None,
    force: bool = False,
    verbose: bool = True,
    log: Log = Log()
) -> pd.DataFrame:
    """
    Infer effect allele frequency (EAF) in summary statistics from MAF using reference VCF ALT frequency.

    This function infers the effect allele frequency (EAF) in the summary statistics by first
    extracting the reference allele frequency from a VCF file, then using this information along
    with the summary statistics MAF to calculate the correct EAF. It handles cases where the
    effect allele might need to be flipped based on the reference data.

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame containing variant data.
    ref_infer : str
        Path to the reference VCF file.
    ref_alt_freq : str, optional
        Field name for alternative allele frequency in VCF INFO section.
    threads : int, default=1
        Number of CPU cores to use for parallel processing.
    chr : str, default="CHR"
        Column name for chromosome information.
    pos : str, default="POS"
        Column name for position information.
    ref : str, default="NEA"
        Column name for reference allele (non-effect allele).
    alt : str, default="EA"
        Column name for alternative allele (effect allele).
    eaf : str, default="EAF"
        Column name for effect allele frequency.
    maf : str, default="MAF"
        Column name for minor allele frequency.
    ref_eaf : str, default="_REF_EAF"
        Temporary column name for storing reference allele frequency.
    status : str, default="STATUS"
        Column name for status codes.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper.
    force : bool, default=False
        If True, infer EAF for all variants regardless of status.
    verbose : bool, default=True
        If True, print progress messages.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object for recording process information.

    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with inferred EAF values.

    Notes
    -----
    - The function first extracts reference allele frequencies from the VCF file.
    - It then compares these reference frequencies with the MAF in the summary statistics
      to determine if flipping is needed.
    - If the reference allele frequency and MAF suggest different major/minor alleles,
      the EAF is calculated as 1 - MAF (flipping the allele).
    - The temporary reference EAF column (_REF_EAF) is dropped before returning the DataFrame.
    - The function provides statistics about the number of variants for which EAF was
      successfully inferred and those still missing EAF values.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
        is_dataframe = True
        sumstats_obj = None
    else:
        sumstats_obj = sumstats_or_dataframe
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Get mapper from Sumstats object if available
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chr in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chr])
    
    # Auto-detect reference format from VCF file
    mapper.detect_reference_format(ref_infer)
    
    if eaf not in sumstats.columns:
        sumstats[eaf]=np.nan
    if ref_eaf not in sumstats.columns:
        sumstats[ref_eaf]=np.nan

    prenumber = sumstats[eaf].isna().sum()

    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)   
        if not force:
            # Match: digit 4 is 0 - Optimized: extract digit once
            digit_4 = _extract_status_digit(sumstats[status], 4)
            good_chrpos = (digit_4 == 0)
        else:
            good_chrpos = pd.Series([True] * len(sumstats), index=sumstats.index)
        
        # Only process variants with valid MAF (required for EAF inference)
        good_chrpos = good_chrpos & sumstats[maf].notna()
        log.write(" -Checking variants:", good_chrpos.sum(), verbose=verbose)
    
        ########################  
        #extract ref af
        if good_chrpos.sum() < 10000:
            threads = 1       
        #df_split = np.array_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], threads)
        df_split = _df_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], threads)
        with Pool(threads) as pool:
            map_func = partial(inferaf,chr=chr,pos=pos,ref=ref,alt=alt,eaf=ref_eaf,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,mapper=mapper) 
            result = pd.concat(pool.map(map_func,df_split))
            # Cast to match original column dtype to avoid FutureWarning
            original_dtype = sumstats[ref_eaf].dtype
            result[ref_eaf] = result[ref_eaf].astype(original_dtype)
            sumstats.loc[good_chrpos,[ref_eaf]] = result[ref_eaf]
        
        ###########################
        # infer sumstats EAF 
        # based on sumstats MAF and reference EAF
        # Optimized: Use XOR logic - flip when ref_eaf and maf indicate different major/minor alleles
        # (ref_eaf >= 0.5) means ref allele is major, (maf <= 0.5) means maf is minor
        # If ref is major and maf is minor, or ref is minor and maf is major, we need to flip
        is_flipped = (sumstats[ref_eaf] >= 0.5) != (sumstats[maf] > 0.5)
        # Vectorized assignment: use np.where for efficient conditional assignment
        sumstats.loc[good_chrpos, eaf] = np.where(
            is_flipped[good_chrpos],
            1 - sumstats.loc[good_chrpos, maf],
            sumstats.loc[good_chrpos, maf]
        )
        log.write(" -Flipping MAF to obtain EAF for {} variants".format(is_flipped[good_chrpos].sum()),verbose=verbose)

        ###########################    
        afternumber = sumstats[eaf].isna().sum()
        log.write(" -Inferred EAF for {} variants.".format(prenumber - afternumber),verbose=verbose) 
        log.write(" -EAF is still missing for {} variants.".format(afternumber),verbose=verbose) 
        sumstats = sumstats.drop(columns=[ref_eaf])
    
    # Set metadata if Sumstats object is available
    try:
        sumstats_obj = getattr(log, '_sumstats_obj', None)
        if sumstats_obj is not None:
            from gwaslab.info.g_meta import _append_meta_record
            sumstats_obj.meta["gwaslab"]["references"]["ref_infer_maf"] = _append_meta_record(
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_maf"], ref_infer)
    except:
        pass
    
    return sumstats

def inferaf(
    sumstats: pd.DataFrame,
    ref_infer: str,
    ref_alt_freq: Optional[str] = None,
    chr: str = "CHR",
    pos: str = "POS",
    ref: str = "NEA",
    alt: str = "EA",
    eaf: str = "EAF",
    mapper: Optional['ChromosomeMapper'] = None
) -> pd.DataFrame:
    #vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    vcf_reader = VariantFile(ref_infer)
    # Auto-detect reference format from VCF file if mapper provided
    if mapper is not None:
        mapper.detect_reference_format(ref_infer)
    def afapply(x,vcf,alt_freq,mapper):
            return infer_af(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],vcf_reader,ref_alt_freq,mapper)
    map_func = partial(afapply,vcf=vcf_reader,alt_freq=ref_alt_freq,mapper=mapper)
    status_inferred = sumstats.apply(map_func,axis=1)
    sumstats[eaf] = status_inferred.values
    sumstats[eaf]=sumstats[eaf].astype("float") 
    return sumstats

def infer_af(
    chr: Union[str, int],
    start: int,
    end: int,
    ref: str,
    alt: str,
    vcf_reader: VariantFile,
    alt_freq: str,
    mapper: Optional['ChromosomeMapper'] = None
) -> float:
    if mapper is not None:
        # Convert sumstats chr to reference format for VCF lookup
        chr = mapper.sumstats_to_reference(chr, reference_file=None)
    chr_seq = vcf_reader.fetch(chr,start,end)
    
    for record in chr_seq:
        if record.pos==end:
            if record.ref==ref and (alt in record.alts):
                return record.info[alt_freq][0]
            elif record.ref==alt and (ref in record.alts):
                return 1 - record.info[alt_freq][0]
    return np.nan
