from __future__ import annotations

import gc
from typing import Tuple, List, Any
from multiprocessing import Pool
from functools import partial

import pandas as pd
import numpy as np

from liftover import get_lifter
from liftover import ChainFile

from gwaslab.info.g_Log import Log
from gwaslab.info.g_vchange_status import match_status, ensure_status_int
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chain
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_fix_sumstats import _fix_chr, _fix_pos

try:
    from sumstats_liftover import liftover_df
except ImportError:
    raise ImportError(
        "sumstats_liftover package is required. Install it with: pip install sumstats-liftover"
    )


def _is_standard_chromosome(chrom: str) -> bool:
    """
    Check if a chromosome name represents a standard chromosome.
    
    Standard chromosomes are: 1-22, X, Y, M/MT (and their 'chr' prefixed versions).
    Filters out alternate contigs, unplaced sequences, etc.
    
    Parameters
    ----------
    chrom : str
        Chromosome name to check
        
    Returns
    -------
    bool
        True if it's a standard chromosome, False otherwise
    """
    chrom_str = str(chrom).strip()
    # Remove 'chr' prefix (case-insensitive)
    if chrom_str.lower().startswith('chr'):
        chrom_str = chrom_str[3:]
    
    # Check for standard chromosomes
    # Standard autosomes: 1-22
    # Standard sex chromosomes: X, Y
    # Standard mitochondrial: M, MT
    # Also accept numeric special chromosomes: 23, 24, 25
    if chrom_str in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                     "21", "22", "23", "24", "25",
                     "X", "Y", "M", "MT", "x", "y", "m", "mt"]:
        return True
    
    # Check if it's a numeric chromosome 1-22
    try:
        chrom_num = int(chrom_str)
        if 1 <= chrom_num <= 22:
            return True
    except (ValueError, TypeError):
        pass
    
    return False


def _normalize_chrom_name(chrom: str) -> str:
    """
    Normalize chromosome name by removing 'chr' prefix and handling special cases.
    Converts numeric special chromosomes to string format for matching.
    
    Examples:
        'chr1' -> '1'
        '1' -> '1'
        'chrX' -> 'X'
        'X' -> 'X'
        '23' -> 'X' (for matching with chain files)
        'chrM' -> 'M'
        'M' -> 'M'
        '25' -> 'M' (for matching with chain files)
    """
    chrom_str = str(chrom).strip()
    # Remove 'chr' prefix (case-insensitive)
    if chrom_str.lower().startswith('chr'):
        chrom_str = chrom_str[3:]
    # Convert numeric special chromosomes to string format for matching
    # (chain files use "X", "Y", "M", not "23", "24", "25")
    if chrom_str == "23":
        chrom_str = "X"
    elif chrom_str == "24":
        chrom_str = "Y"
    elif chrom_str == "25":
        chrom_str = "M"
    return chrom_str


# ----------------------------
# Helper functions
# ----------------------------

def _df_split(dataframe: pd.DataFrame, n: int) -> List[pd.DataFrame]:
    """Split a DataFrame into n approximately equal parts for parallel processing."""
    k, m = divmod(len(dataframe), n)
    return [dataframe.iloc[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]


# ----------------------------
# Original liftover functions (using liftover package)
# ----------------------------

def liftover_snv(row, chrom, converter, to_build) -> Tuple[Any, Any, int]:
    """
    Perform liftover for a single variant (SNV).
    
    Converts position from one genome build to another and updates status code.
    Returns (chrom, pos, status) tuple or (NA, NA, status) if unmapped.
    
    Parameters
    ----------
    row : pd.Series
        Row containing POS (iloc[0]) and STATUS (iloc[1])
    chrom : str
        Chromosome name for conversion
    converter : ChainFile or Lifter
        Converter object for coordinate transformation
    to_build : str
        Target genome build number
        
    Returns
    -------
    Tuple[Any, Any, int]
        (chrom, pos, status) tuple with lifted coordinates
    """
    # Convert status to integer and extract digits efficiently
    status_val = row.iloc[1]
    if not isinstance(status_val, (int, np.integer)):
        status_val = int(str(status_val))
    
    # Extract digits: digit 3, set digit 4 to 9, digit 5, set digits 6-7 to 99
    # Optimize: use bit operations where possible, but division is clearer here
    digit3 = (status_val // 10000) % 10
    digit5 = (status_val // 100) % 10
    # Build status code: build (2 digits) + digit3 + 9 + digit5 + 99
    status_end_int = digit3 * 10000 + 9000 + digit5 * 100 + 99
    
    # Get input position (1-based from GWAS data)
    pos_input = int(row.iloc[0])
    
    # Convert position to 0-based for converter (chain files use 0-based coordinates)
    # Converter is created with one_based=False, so it expects 0-based input and returns 0-based output
    # 
    # WHY THE BUG ONLY AFFECTED ~0.01% OF VARIANTS:
    # The old implementation used one_based=True, which caused a 2 bp offset only for variants where:
    # 1. The exception handler was triggered (e.g., certain chain file versions, specific chromosomes)
    # 2. The fallback converter was created without one_based parameter (default behavior may differ)
    # 3. OR in specific chain file segments (inverted regions, boundaries) where coordinate handling
    #    behaved differently with one_based=True vs False
    # 
    # The coordinate mismatch occurred because:
    # - Input: 1-based position (e.g., 1000)
    # - Converted to 0-based: 999
    # - Passed to converter with one_based=True: converter expects 1-based, treats 999 as position 1000
    # - Converter returns 1-based result (e.g., 2000)
    # - We add 1 again: 2000 + 1 = 2001 (should be 2000) → 2 bp offset
    # 
    # This only manifested when the exception path was used or in edge cases where the converter's
    # internal coordinate handling differed, explaining why only ~1000 out of 10M variants were affected.
    pos_0_based = pos_input - 1
    
    # Get conversion result (cache to avoid double lookup)
    try:
        results = converter[chrom][pos_0_based]
    except (KeyError, IndexError):
        # Chromosome or position not in converter
        return pd.NA, pd.NA, 9700000 + status_end_int
    
    # Check if conversion succeeded
    if not results:
        return pd.NA, pd.NA, 9700000 + status_end_int
    
    # Extract and process results
    result_chrom = results[0][0]
    result_pos = results[0][1]  # This is 0-based (since converter uses one_based=False)
    
    # Strip "chr" prefix and compare (cache stripped value)
    chrom_stripped = result_chrom.strip("chr")
    if chrom_stripped != chrom:
        # Chromosome mismatch - unmapped
        return pd.NA, pd.NA, 9700000 + status_end_int
    
    # Successful conversion: return lifted coordinates
    # Convert 0-based result back to 1-based for GWAS output
    # FIXED: Previously had double conversion issue - converter with one_based=True was returning
    # 1-based, but we were adding 1 again, causing 2 bp offset. Now using one_based=False
    # consistently: input 1-based -> convert to 0-based -> converter returns 0-based -> convert to 1-based
    return chrom_stripped, result_pos + 1, int(to_build) * 100000 + status_end_int


def liftover_variant(sumstats, 
             chrom="CHR", 
             pos="POS", 
             status="STATUS",
             from_build="19", 
             to_build="38",
             chain=None) -> pd.DataFrame:
    """
    Perform liftover conversion for variants in a dataframe chunk.
    
    This function is called in parallel by _parallelize_liftover_variant.
    It converts genomic coordinates from one build to another using UCSC chain files.
    """
    # Early exit if empty
    if len(sumstats) == 0:
        return sumstats
    
    # Create converter (this is per-chunk, but converter creation is relatively fast)
    # Use one_based=False to work with 0-based coordinates internally (chain files are 0-based)
    # This avoids coordinate system confusion
    # 
    # NOTE: The old implementation used one_based=True, which caused a 2 bp offset for a small
    # subset of variants (~0.01%). This happened because:
    # 1. When one_based=True, the converter expects 1-based input and returns 1-based output
    # 2. But we were converting 1-based input to 0-based before passing it
    # 3. The converter treated the 0-based position as 1-based (off by 1)
    # 4. Then we added 1 again to convert back to 1-based (total off by 2)
    # 5. The issue only affected variants where the exception handler triggered and used
    #    a fallback converter, or in specific chain file segments where coordinate handling
    #    differed (e.g., inverted regions, boundary conditions, or certain chromosome regions)
    try:
        if chain is None:
            converter = get_lifter("hg{}".format(from_build), "hg{}".format(to_build), one_based=False)
        else:
            converter = ChainFile(chain, target="", query="", one_based=False)
    except Exception:
        # Fallback: try without one_based parameter (defaults to False/0-based)
        # This ensures consistent 0-based coordinate handling even if one_based parameter fails
        if chain is None:
            converter = get_lifter("hg{}".format(from_build), "hg{}".format(to_build))
        else:
            converter = ChainFile(chain, target="", query="")

    # Get chromosome dictionaries once
    dic = get_number_to_chr(in_chr=False, xymt=["X", "Y", "M"])
    dic2 = get_chr_to_number(out_chr=False)
    
    # Process each chromosome separately
    for chrom_val in sumstats[chrom].unique():
        chrom_mask = sumstats[chrom] == chrom_val
        if not chrom_mask.any():
            continue
        
        chrom_to_convert = dic[chrom_val]
        chrom_subset = sumstats.loc[chrom_mask, [pos, status]]
        
        # Apply liftover to all variants on this chromosome
        # liftover_snv returns (chrom, pos, status) tuple
        lifted = chrom_subset.apply(
            lambda row: liftover_snv(row, chrom_to_convert, converter, to_build),
            axis=1
        )
        
        # Extract tuple elements efficiently using list comprehension (faster than multiple apply calls)
        # Handle both tuple and NA cases
        lifted_chrom = [x[0] if isinstance(x, tuple) and len(x) > 0 else pd.NA for x in lifted]
        lifted_pos = [x[1] if isinstance(x, tuple) and len(x) > 1 else pd.NA for x in lifted]
        lifted_status = [x[2] if isinstance(x, tuple) and len(x) > 2 else pd.NA for x in lifted]
        
        # Update dataframe with lifted values (preserve index)
        sumstats.loc[chrom_mask, pos] = lifted_pos
        sumstats.loc[chrom_mask, status] = pd.Series(lifted_status, index=sumstats.loc[chrom_mask].index, dtype='Int64')
        sumstats.loc[chrom_mask, chrom] = pd.Series(lifted_chrom, index=sumstats.loc[chrom_mask].index).map(dic2).astype("Int64")
    
    return sumstats


@with_logging(
        start_to_msg= "perform liftover",
        finished_msg= "liftover",
        start_function= ".liftover()",
        start_cols=["CHR","POS","STATUS"]
)
def _parallelize_liftover_variant(sumstats_obj,threads=1,chrom="CHR", pos="POS", from_build="19", to_build="38",status="STATUS",remove=True,chain=None, verbose=True,log=Log()) -> pd.DataFrame:
    '''
    Perform parallelized liftover of variants to a new genome build.
    
    Converts genomic coordinates from one genome build (e.g., hg19/GRCh37) to another
    (e.g., hg38/GRCh38) using UCSC chain files. Only processes variants with valid
    coordinates (status code xxx0xxx). After liftover, validates and fixes chromosome
    and position columns.

    Parameters
    ----------
    sumstats_obj : Sumstats
        Sumstats object containing the data to liftover.
    threads : int
        Number of CPU cores to use for parallel processing.
    from_build : str
        Name of the original genome build (e.g., "19").
    to_build : str
        Name of the target genome build (e.g., "38").
    remove : bool, default False
        If True, remove variants that fail to map.
    verbose : bool, default False
        If True, print progress messages during processing.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table with updated genomic coordinates.
    '''
    # Helper function to update metadata for Sumstats object (defined early for use in early exit)
    def _update_liftover_metadata(sumstats_obj_ref, to_build, from_build, remove, chain, threads, log):
        """Update metadata after liftover operation."""
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            from gwaslab.qc.qc_build import _process_build
            sumstats_obj_ref.meta["is_sorted"] = False
            sumstats_obj_ref.meta["is_harmonised"] = False
            sumstats_obj_ref.meta["gwaslab"]["genome_build"] = _process_build(to_build, log=log, verbose=False)
            sumstats_obj_ref.build = to_build
            liftover_kwargs = {
                'from_build': from_build, 'to_build': to_build, 'remove': remove, 'chain': chain, 'threads': threads
            }
            _update_harmonize_step(sumstats_obj_ref, "liftover", liftover_kwargs, True)
        except Exception:
            # Silently fail if metadata update is not available
            pass
    
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        # Called during initialization - no Sumstats object yet
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    lifter_from_build = _process_build(from_build, log=log, verbose=False)
    lifter_to_build = _process_build(to_build, log=log, verbose=False)

    if chain is not None:
        log.log_operation("Creating converter using ChainFile: {}".format(chain), verbose=verbose)
    else:
        try:
            chain = get_chain(from_build=from_build, to_build=to_build)
            if chain is None or chain==False:
                raise ValueError("No available chain file for {} -> {}".format(lifter_from_build, lifter_to_build))
            log.log_operation("Creating converter using provided ChainFile: {}".format(chain), verbose=verbose)
        except:
            chain = None
            log.log_operation("Try creating converter using liftover package", verbose=verbose)

    log.log_operation("Creating converter: {} -> {}".format(lifter_from_build, lifter_to_build), verbose=verbose)
    
    # Find variants with valid coordinates (status code xxx0xxx)
    pattern = r"\w\w\w0\w\w\w"  
    to_lift = match_status(sumstats[status], pattern)
    to_lift_count = to_lift.sum()
    
    # Early exit if no variants to lift
    if to_lift_count == 0:
        log.log_operation("No variants with status code xxx0xxx to liftover", verbose=verbose)
        if not is_dataframe:
            _update_liftover_metadata(sumstats_obj, to_build, from_build, remove, chain, threads, log)
        return sumstats_obj.data if not is_dataframe else sumstats
    
    # Filter to variants that need lifting
    sumstats = sumstats.loc[to_lift, :].copy()
    log.log_operation("Converting variants with status code xxx0xxx: {:,}".format(to_lift_count), verbose=verbose)
    
    ###########################################################################
    # Auto-adjust threads for small datasets
    if to_lift_count < 10000:
        threads = 1
    
    # Parallel liftover processing
    df_split = _df_split(sumstats[[chrom, pos, status]], threads)
    pool = Pool(threads)
    try:
        func = liftover_variant
        map_func = partial(
            func, chrom=chrom, pos=pos, from_build=lifter_from_build,
            to_build=lifter_to_build, status=status, chain=chain
        )
        lifted_results = pool.map(map_func, df_split)
        # Efficiently concatenate and assign results
        sumstats[[chrom, pos, status]] = pd.concat(lifted_results, ignore_index=False)
    finally:
        pool.close()
        pool.join()
        del df_split, lifted_results
        gc.collect()
    ############################################################################
    
    # Count unmapped variants efficiently
    is_unmapped = sumstats[pos].isna()
    unmap_num = is_unmapped.sum()
    
    if remove and unmap_num > 0:
        log.log_variants_removed(unmap_num, reason="unmapped variants", verbose=verbose)
        sumstats = sumstats.loc[~is_unmapped, :]
    elif unmap_num > 0:
        log.write(" -Unmapped variants: {:,} (not removed)".format(unmap_num), verbose=verbose)
    
    # After liftover, validate and fix chr and pos
    if is_dataframe:
        sumstats = _fix_chr(sumstats, chrom=chrom, add_prefix="", remove=remove, verbose=verbose)
        sumstats = _fix_pos(sumstats, pos=pos, remove=remove, verbose=verbose)
    else:
        # Assign filtered dataframe back before calling fix functions
        sumstats_obj.data = sumstats
        sumstats_obj.data = _fix_chr(sumstats_obj, chrom=chrom, add_prefix="", remove=remove, verbose=verbose)
        sumstats_obj.data = _fix_pos(sumstats_obj, pos=pos, remove=remove, verbose=verbose)

    # Set metadata and update harmonization status if called with Sumstats object
    # This handles both normal execution and mocked/test scenarios
    if not is_dataframe:
        # Update metadata on the Sumstats object
        _update_liftover_metadata(sumstats_obj, to_build, from_build, remove, chain, threads, log)
        # Always return DataFrame (sumstats_obj.data contains the processed DataFrame)
        return sumstats_obj.data
    else:
        # Always return DataFrame
        return sumstats


# ----------------------------
# Fast liftover2 functions (using sumstats_liftover package)
# ----------------------------


@with_logging(
        start_to_msg= "perform liftover2",
        finished_msg= "liftover2",
        start_function= ".liftover2()",
        start_cols=["CHR","POS","STATUS"]
)
def _liftover2_variant(sumstats_obj,
                      chain_path: str,
                      chrom: str = "CHR",
                      pos: str = "POS",
                      status: str = "STATUS",
                      to_build: str = None,
                      out_chrom_col: str = "CHR_LIFT",
                      out_pos_col: str = "POS_LIFT",
                      out_strand_col: str = "STRAND_LIFT",
                      one_based_input: bool = True,
                      one_based_output: bool = True,
                      remove: bool = True,
                      verbose: bool = True,
                      log: Log = Log()) -> pd.DataFrame:
    """
    Perform liftover2 conversion for variants using a UCSC chain file.
    
    This is a fast chain-based liftover implementation that directly parses
    UCSC chain files and performs vectorized coordinate conversion. It converts
    genomic coordinates from one genome build (e.g., hg19/GRCh37) to another
    (e.g., hg38/GRCh38) using UCSC chain files.
    
    **Step-by-Step Process:**
    
    1. **Chain File Parsing** (`parse_chain_to_segments`)
       - Reads the UCSC chain file (supports both .chain and .chain.gz formats)
       - Parses chain headers containing alignment information:
         * Target chromosome (tName), query chromosome (qName)
         * Alignment scores, strand information, chromosome sizes
       - Extracts ungapped alignment blocks from each chain
       - Normalizes chromosome names (removes 'chr' prefix for consistent matching)
       - Builds a best-disjoint-cover index for fast point lookups
       
       Example chain header:
       ```
       chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
       ```
       This indicates: score=20851231461, target=chr1, query=chr1, both forward strand
       
       Example alignment block:
       ```
       167376  50041   80290
       ```
       This means: size=167376 bases aligned, then gap of 50041 on target, 80290 on query
    
    2. **Chromosome Name Normalization**
       - Input chromosomes from sumstats are normalized for matching:
         * Autosomal: "1" or "chr1" → "1"
         * Special chromosomes: "X"/"chrX"/"23" → "X", "Y"/"chrY"/"24" → "Y", "M"/"MT"/"chrM"/"chrMT"/"25" → "M"
       - Chain file chromosomes are normalized when stored (e.g., "chr1" → "1", "chrX" → "X")
       - This ensures matching works regardless of naming convention
       - Special chromosomes are converted to numeric values in the output:
         * X → 23
         * Y → 24
         * M or MT → 25
       
       Examples (matching phase):
       - Sumstats has CHR="1" → normalized to "1" → matches chain "chr1" (normalized to "1")
       - Sumstats has CHR="chr1" → normalized to "1" → matches chain "chr1" (normalized to "1")
       - Sumstats has CHR="X" or "chrX" or "23" → normalized to "X" for matching → matches chain "chrX" (normalized to "X")
       - Sumstats has CHR="Y" or "chrY" or "24" → normalized to "Y" for matching → matches chain "chrY" (normalized to "Y")
       - Sumstats has CHR="M", "MT", "chrM", "chrMT", or "25" → normalized to "M" for matching → matches chain "chrM" (normalized to "M")
       
       Examples (output phase):
       - All special chromosomes in output are converted to numeric: "X" → "23", "Y" → "24", "M"/"MT" → "25"
    
    3. **Position Conversion** (`liftover_sumstats_df`)
       - For each variant, converts input position to 0-based if needed:
         * 1-based input: pos_0based = pos - 1
         * 0-based input: pos_0based = pos
       - Looks up the position in the best-disjoint-cover index using binary search
       - Finds the matching alignment segment that contains the position
       
       Example:
       - Input: CHR=1, POS=725932 (1-based)
       - Converted: CHR="1" (normalized), POS=725931 (0-based)
       - Lookup: Finds segment where 725931 falls within [t0, t1) interval
    
    4. **Coordinate Mapping**
       - Calculates offset within the alignment segment:
         offset = position_0based - segment_start
       - Maps to query coordinates:
         * Forward strand: query_pos = segment_query_start + offset
         * Reverse strand: query_pos = query_size - 1 - (segment_query_start + offset)
       - Converts back to requested coordinate system (0-based or 1-based)
       
       Example (forward strand):
       - Segment: target [725000, 726000) → query [730000, 731000)
       - Variant at target position 725932 (0-based)
       - Offset: 725932 - 725000 = 932
       - Query position: 730000 + 932 = 730932 (0-based)
       - If 1-based output: 730932 + 1 = 730933
       
       Example (reverse strand):
       - Segment: target [725000, 726000) → query reverse [730000, 731000)
       - Variant at target position 725932 (0-based)
       - Offset: 725932 - 725000 = 932
       - Query reverse position: 730000 + 932 = 730932
       - Query forward position: query_size - 1 - 730932
       - If query_size = 248956422: 248956422 - 1 - 730932 = 248225489 (0-based)
       - If 1-based output: 248225490
    
    5. **Strand Assignment**
       - Records whether the mapping is on forward (+) or reverse (-) strand
       - Forward strand: query coordinates increase in same direction as target
       - Reverse strand: query coordinates are on reverse complement
       
       Example:
       - If chain has qStrand='+': output strand = "+"
       - If chain has qStrand='-': output strand = "-"
    
    6. **Unmapped Variant Handling**
       - Variants that don't fall within any alignment segment are marked as unmapped
       - Unmapped variants get: CHR_LIFT=None, POS_LIFT=-1 or NaN, STRAND_LIFT=None
       - For unmapped variants: CHR and POS are set to NA (matching original liftover behavior)
       - Status codes are updated: 9700000 + (digit3 * 10000 + 9000 + digit5 * 100 + 99)
       - If remove=True, unmapped variants are filtered out
       
       Example:
       - Variant at CHR=1, POS=999999999 (outside all alignment blocks)
       - Result: CHR=NA, POS=NA, STATUS=9700000 + status_end_int
       - If remove=True: variant is removed from output
    
    7. **Status Code Updates** (if STATUS column exists)
       - Extracts digit3 and digit5 from original status code
       - For mapped variants: status = to_build * 100000 + (digit3 * 10000 + 9000 + digit5 * 100 + 99)
       - For unmapped variants: status = 9700000 + (digit3 * 10000 + 9000 + digit5 * 100 + 99)
       - This matches the behavior of the original liftover function
       
       Example:
       - Original status: 1900000 (build 19, valid coordinates)
       - After liftover to build 38 (mapped): status = 3800000 + status_end_int
       - After liftover (unmapped): status = 9700000 + status_end_int
    
    8. **Output Column Creation and Overwrite**
       - Creates temporary columns with lifted coordinates:
         * CHR_LIFT: Lifted chromosome name (with 'chr' prefix removed, special chromosomes as 23/24/25)
         * POS_LIFT: Lifted position (in requested coordinate system)
         * STRAND_LIFT: Strand of the lifted position ("+" or "-")
       - Overwrites original CHR and POS columns with lifted values (only for mapped variants)
       - For unmapped variants: CHR and POS are set to NA (matching original liftover behavior)
       - Drops temporary columns (CHR_LIFT, POS_LIFT, STRAND_LIFT) from final output
    
    **Performance Features:**
    - Vectorized operations for batch processing
    - Binary search on sorted intervals for O(log n) lookup
    - Best-disjoint-cover index eliminates overlap resolution during lookup
    - Per-chromosome processing for memory efficiency
    
    **Parameters**
    ----------
    sumstats_obj : Sumstats or pd.DataFrame
        Sumstats object or DataFrame containing the data to liftover.
        Must contain columns specified by `chrom` and `pos` parameters.
        
        Example input:
        ```
        CHR  POS      EA  NEA
        1    725932   G   A
        1    725933   A   G
        2    100000   C   T
        ```
    
    chain_path : str
        Path to UCSC chain file (can be .chain or .chain.gz).
        Chain files can be downloaded from UCSC Genome Browser.
        
        Examples:
        - "/path/to/hg19ToHg38.over.chain.gz"
        - "/path/to/hg38ToHg19.over.chain"
        - Common chain files:
          * hg19ToHg38.over.chain.gz (hg19 → hg38)
          * hg38ToHg19.over.chain.gz (hg38 → hg19)
    
    chrom : str, default "CHR"
        Column name for chromosome in the input dataframe.
        Can be numeric (1, 2, 3..., 23, 24, 25) or string ("1", "chr1", "X", "chrX", etc.).
        Special chromosomes can be:
        - X chromosome: "X", "chrX", or 23
        - Y chromosome: "Y", "chrY", or 24
        - Mitochondrial: "M", "MT", "chrM", "chrMT", or 25
        Will be normalized internally for matching (numeric special chromosomes converted to string).
        
        Example: If your dataframe has column "Chromosome", set chrom="Chromosome"
    
    pos : str, default "POS"
        Column name for position in the input dataframe.
        Must be integer positions.
        
        Example: If your dataframe has column "BP", set pos="BP"
    
    status : str, default "STATUS"
        Column name for status code in the input dataframe.
        Status codes will be updated to reflect the new build and mapping status.
        - Mapped variants: status = to_build * 100000 + (digit3 * 10000 + 9000 + digit5 * 100 + 99)
        - Unmapped variants: status = 9700000 + (digit3 * 10000 + 9000 + digit5 * 100 + 99)
        If STATUS column doesn't exist, status codes will not be updated.
        
        Example: If your dataframe has column "STAT", set status="STAT"
    
    to_build : str, optional
        Target genome build number (e.g., "19", "38").
        If not provided, will be inferred from chain file name.
        Examples: "19" (hg19/GRCh37), "38" (hg38/GRCh38)
        
        Example: to_build="38"
    
    out_chrom_col : str, default "CHR_LIFT"
        Temporary column name for lifted chromosome (will be dropped after overwriting CHR).
        Chromosome names will have 'chr' prefix removed and special chromosomes
        converted to numeric values:
        - Autosomal: "chr1" → "1", "chr2" → "2", etc.
        - X chromosome: "chrX" or "X" → "23"
        - Y chromosome: "chrY" or "Y" → "24"
        - Mitochondrial: "chrM", "chrMT", "M", or "MT" → "25"
        
        Note: This column is used internally and then dropped. The lifted values
        overwrite the original CHR column.
    
    out_pos_col : str, default "POS_LIFT"
        Temporary column name for lifted position (will be dropped after overwriting POS).
        Position will be in the coordinate system specified by `one_based_output`.
        Unmapped variants will have -1 or NaN.
        
        Note: This column is used internally and then dropped. The lifted values
        overwrite the original POS column.
    
    out_strand_col : str, default "STRAND_LIFT"
        Temporary column name for lifted strand (will be dropped from final output).
        Values will be "+" (forward) or "-" (reverse).
        Unmapped variants will have None.
        
        Note: This column is created during liftover but dropped from the final output.
        Strand information is not preserved in the final dataframe.
    
    one_based_input : bool, default True
        Whether input positions are 1-based (GWAS standard) or 0-based (BED format).
        - True: Input positions are 1-based (first base = 1)
        - False: Input positions are 0-based (first base = 0)
        
        Examples:
        - GWAS sumstats typically use 1-based: one_based_input=True
        - BED files use 0-based: one_based_input=False
    
    one_based_output : bool, default True
        Whether output positions should be 1-based or 0-based.
        - True: Output positions will be 1-based (GWAS standard)
        - False: Output positions will be 0-based (BED format)
        
        Examples:
        - For GWAS output: one_based_output=True
        - For BED file output: one_based_output=False
    
    remove : bool, default True
        If True, remove variants that fail to map (unmapped variants).
        If False, keep unmapped variants with CHR_LIFT=None, POS_LIFT=-1.
        
        Examples:
        - remove=True: Only mapped variants in output
        - remove=False: All variants kept, unmapped marked with None/-1
    
    verbose : bool, default True
        If True, print progress messages including:
        - Chain file path
        - Coordinate system information
        - Mapping statistics (mapped/unmapped counts)
        - Removal statistics
    
    log : Log, default Log()
        Logging object for output messages.
        Used internally for consistent logging format.
        
    **Returns**
    -------
    pd.DataFrame
        Summary statistics table with lifted coordinates added as new columns.
        Original columns are preserved. New columns are added:
        - `out_chrom_col`: Lifted chromosome
        - `out_pos_col`: Lifted position  
        - `out_strand_col`: Lifted strand
        
        Example output (CHR and POS are overwritten, temporary columns dropped):
        ```
        CHR  POS      EA  NEA
        1    730933   G   A
        1    730934   A   G
        2    150000   C   T
        23   550000   A   G
        24   120000   C   T
        25   1500     A   G
        ```
        
        Note: Original CHR and POS values are replaced with lifted coordinates.
        Temporary columns (CHR_LIFT, POS_LIFT, STRAND_LIFT) are not in the final output.
        
        If remove=True and variants are unmapped, they are excluded from output.
        If remove=False, unmapped variants have:
        - CHR = NA (set to NA, not original value)
        - POS = NA (set to NA, not original value)
        - STATUS = 9700000 + status_end_int (if STATUS column exists)
    
    **Examples**
    --------
    Basic usage:
    >>> mysumstats.liftover2(chain_path="/path/to/hg19ToHg38.over.chain.gz")
    
    Custom column names:
    >>> mysumstats.liftover2(
    ...     chain_path="/path/to/hg19ToHg38.over.chain.gz",
    ...     chrom="Chromosome",
    ...     pos="BP",
    ...     out_chrom_col="CHR_hg38",
    ...     out_pos_col="POS_hg38"
    ... )
    
    Keep unmapped variants:
    >>> mysumstats.liftover2(
    ...     chain_path="/path/to/hg19ToHg38.over.chain.gz",
    ...     remove=False
    ... )
    
    **Notes**
    -----
    - Chain files use 0-based, half-open intervals [start, end)
    - Chromosome name normalization handles both "1" and "chr1" formats
    - Special chromosomes are converted to numeric values in output:
      X → 23, Y → 24, M/MT → 25 (consistent with GWASLab standard)
    - Reverse strand mappings convert coordinates to forward strand
    - Unmapped variants occur when positions fall outside alignment blocks
    - This implementation is faster than the original liftover for large datasets
      due to vectorized operations and optimized indexing
    """
    import pandas as pd
    
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    if len(sumstats) == 0:
        log.write(" -Empty dataframe, skipping liftover2.", verbose=verbose)
        return sumstats if is_dataframe else sumstats_obj.data
    
    # Infer to_build from chain_path if not provided
    if to_build is None:
        import re
        # Try to extract build number from chain file name (e.g., hg19ToHg38 -> 38)
        chain_name = chain_path.split('/')[-1]  # Get filename
        match = re.search(r'ToHg(\d+)|To(\d+)', chain_name, re.IGNORECASE)
        if match:
            to_build = match.group(1) or match.group(2)
        else:
            # Try to infer from common patterns
            if 'hg38' in chain_name.lower() or 'grch38' in chain_name.lower():
                to_build = "38"
            elif 'hg19' in chain_name.lower() or 'grch37' in chain_name.lower():
                to_build = "19"
            else:
                to_build = "99"  # Unknown build
        log.write(f" -Inferred target build: {to_build} from chain file name", verbose=verbose)
    
    # Process build number
    from gwaslab.qc.qc_build import _process_build
    to_build_processed = _process_build(to_build, log=log, verbose=False)
    
    log.write(f" -Using chain file: {chain_path}", verbose=verbose)
    log.write(f" -Target build: {to_build_processed}", verbose=verbose)
    log.write(f" -Input positions are {'1-based' if one_based_input else '0-based'}", verbose=verbose)
    log.write(f" -Output positions will be {'1-based' if one_based_output else '0-based'}", verbose=verbose)
    
    # Store original chromosome values for comparison (normalized for matching)
    # This is needed to detect chromosome mismatches (variants mapped to different chromosomes)
    original_chrom_normalized = sumstats[chrom].astype(str).apply(_normalize_chrom_name)
    
    # Perform liftover using sumstats_liftover package
    sumstats = liftover_df(
        df=sumstats,
        chain_path=chain_path,
        chrom_col=chrom,
        pos_col=pos,
        out_chrom_col=out_chrom_col,
        out_pos_col=out_pos_col,
        out_strand_col=out_strand_col,
        one_based_input=one_based_input,
        one_based_output=one_based_output,
        remove_unmapped=False,  # We handle removal separately below
        convert_special_chromosomes=False,  # We handle conversion separately below
    )
    
    # Count unmapped variants (positions that failed to map)
    unmapped = sumstats[out_pos_col].isna() | (sumstats[out_pos_col] == -1)
    
    # Check for chromosome mismatches (variants mapped to different chromosomes)
    # This mimics the old implementation behavior where chromosome mismatches are treated as unmapped
    if out_chrom_col in sumstats.columns:
        # Normalize lifted chromosomes for comparison
        lifted_chrom_normalized = sumstats[out_chrom_col].astype(str).apply(
            lambda x: _normalize_chrom_name(x) if pd.notna(x) and str(x) != '-1' else None
        )
        
        # Mark variants with chromosome mismatches as unmapped
        # Compare normalized chromosome names (handles "1" vs "chr1", "X" vs "23", etc.)
        chrom_mismatch = lifted_chrom_normalized.notna() & (lifted_chrom_normalized != original_chrom_normalized)
        chrom_mismatch_count = chrom_mismatch.sum()
        
        if chrom_mismatch_count > 0:
            log.write(f" -Chromosome mismatches detected: {chrom_mismatch_count} variants (treated as unmapped)", verbose=verbose)
            # Mark chromosome mismatches as unmapped
            unmapped = unmapped | chrom_mismatch
    
    unmapped_count = unmapped.sum()
    mapped_count = len(sumstats) - unmapped_count
    
    log.write(f" -Mapped: {mapped_count} variants", verbose=verbose)
    log.write(f" -Unmapped: {unmapped_count} variants", verbose=verbose)
    
    # Update status codes if STATUS column exists (mimic original liftover behavior)
    if status in sumstats.columns:
        # Ensure STATUS column is integer type before processing
        sumstats = ensure_status_int(sumstats, status)
        # Extract status digits for all variants
        status_vals = sumstats[status].astype(str)
        
        def extract_status_digits(status_str):
            """Extract digit3 and digit5 from status code."""
            try:
                status_int = int(str(status_str))
                digit3 = (status_int // 10000) % 10
                digit5 = (status_int // 100) % 10
                status_end_int = digit3 * 10000 + 9000 + digit5 * 100 + 99
                return status_end_int
            except (ValueError, TypeError):
                return None
        
        status_end_ints = status_vals.apply(extract_status_digits)
        
        # Update status for mapped variants: to_build * 100000 + status_end_int
        if mapped_count > 0:
            mapped_mask = ~unmapped
            mapped_status_end = status_end_ints[mapped_mask]
            valid_mapped = mapped_status_end.notna()
            if valid_mapped.any():
                new_status_mapped = (int(to_build_processed) * 100000 + mapped_status_end[valid_mapped]).astype('Int64')
                sumstats.loc[mapped_mask & valid_mapped, status] = new_status_mapped
        
        # Update status for unmapped variants: 9700000 + status_end_int (build 97 = unmapped)
        if unmapped_count > 0:
            unmapped_status_end = status_end_ints[unmapped]
            valid_unmapped = unmapped_status_end.notna()
            if valid_unmapped.any():
                new_status_unmapped = (9700000 + unmapped_status_end[valid_unmapped]).astype('Int64')
                sumstats.loc[unmapped & valid_unmapped, status] = new_status_unmapped
    
    # Update chromosome format if needed (strip 'chr' prefix and convert special chromosomes)
    if out_chrom_col in sumstats.columns:
        # Convert chromosome names to match expected format
        # First, strip 'chr' prefix
        sumstats[out_chrom_col] = sumstats[out_chrom_col].astype(str).str.replace("^chr", "", regex=True)
        # Only convert special chromosomes for standard chromosomes (filter out alternate contigs)
        # Check which values are standard chromosomes before converting
        is_standard = sumstats[out_chrom_col].apply(_is_standard_chromosome)
        if is_standard.any():
            # Convert special chromosomes to numeric: X→23, Y→24, M/MT→25
            # Only apply to standard chromosomes
            standard_mask = is_standard
            sumstats.loc[standard_mask, out_chrom_col] = sumstats.loc[standard_mask, out_chrom_col].replace({
                "X": "23", "x": "23", "chrX": "23", "chrx": "23",
                "Y": "24", "y": "24", "chrY": "24", "chry": "24",
                "M": "25", "m": "25", "MT": "25", "mt": "25", "chrM": "25", "chrm": "25", "chrMT": "25", "chrmt": "25"
            })
    
    # Overwrite CHR and POS with lifted values, drop temporary columns
    # For unmapped variants, set CHR and POS to NA (mimic original liftover behavior)
    # Note: unmapped mask now includes both position failures and chromosome mismatches
    if out_chrom_col in sumstats.columns and chrom in sumstats.columns:
        # Recalculate mapped_mask after chromosome mismatch detection
        mapped_mask = ~unmapped
        
        # Update mapped variants with lifted chromosome
        if mapped_mask.any():
            sumstats.loc[mapped_mask, chrom] = sumstats.loc[mapped_mask, out_chrom_col]
        
        # Set unmapped variants to NA (before removal, if not removing)
        if unmapped_count > 0 and not remove:
            sumstats.loc[unmapped, chrom] = pd.NA
        
        # Convert to numeric if possible (for special chromosomes 23, 24, 25)
        # Use errors='coerce' to handle non-numeric chromosome names (alternate contigs, etc.)
        # Only convert standard chromosomes that can be safely converted
        try:
            # First, identify which values are standard chromosomes and can be converted
            chrom_series = sumstats[chrom].astype(str)
            # Check which are standard chromosomes (excluding NA)
            is_standard = chrom_series.apply(lambda x: _is_standard_chromosome(x) if pd.notna(x) else False)
            
            # Only try to convert standard chromosomes
            if is_standard.any():
                # For standard chromosomes, try to convert to numeric
                # This handles: "1"-"22" -> 1-22, "23"->23, "24"->24, "25"->25
                # Create a copy for conversion
                chrom_to_convert = chrom_series.copy()
                # Only convert standard chromosomes
                numeric_chrom = pd.to_numeric(chrom_to_convert, errors='coerce')
                # Only update values that are standard AND successfully converted (not NaN)
                valid_mask = is_standard & (~numeric_chrom.isna())
                if valid_mask.any():
                    sumstats.loc[valid_mask, chrom] = numeric_chrom[valid_mask].astype('Int64')
            # Non-standard chromosomes (alternate contigs, etc.) remain as strings
        except Exception:
            # If conversion fails entirely, keep original values as strings
            pass
    
    if out_pos_col in sumstats.columns and pos in sumstats.columns:
        # Recalculate mapped_mask after chromosome mismatch detection
        mapped_mask = ~unmapped
        
        # Update mapped variants with lifted position
        if mapped_mask.any():
            sumstats.loc[mapped_mask, pos] = sumstats.loc[mapped_mask, out_pos_col]
        
        # Set unmapped variants to NA (before removal, if not removing)
        if unmapped_count > 0 and not remove:
            sumstats.loc[unmapped, pos] = pd.NA
        
        # Convert to integer type
        try:
            sumstats[pos] = pd.to_numeric(sumstats[pos], errors='coerce').astype('Int64')
        except:
            pass
    
    # Remove unmapped variants if requested (after setting NA values)
    if remove and unmapped_count > 0:
        sumstats = sumstats[~unmapped].copy()
        log.write(f" -Removed {unmapped_count} unmapped variants", verbose=verbose)
    
    # Drop temporary output columns
    columns_to_drop = [out_chrom_col, out_pos_col, out_strand_col]
    for col in columns_to_drop:
        if col in sumstats.columns:
            sumstats = sumstats.drop(columns=[col])
    
    # After liftover, validate and fix chr and pos (mimic original liftover behavior)
    # This ensures chromosome and position values are properly formatted and validated
    if is_dataframe:
        sumstats = _fix_chr(sumstats, chrom=chrom, add_prefix="", remove=remove, verbose=verbose, log=log)
        sumstats = _fix_pos(sumstats, pos=pos, remove=remove, verbose=verbose, log=log)
    else:
        # Assign filtered dataframe back before calling fix functions
        sumstats_obj.data = sumstats
        sumstats_obj.data = _fix_chr(sumstats_obj, chrom=chrom, add_prefix="", remove=remove, verbose=verbose, log=log)
        sumstats_obj.data = _fix_pos(sumstats_obj, pos=pos, remove=remove, verbose=verbose, log=log)
        sumstats = sumstats_obj.data
    
    # Update QC status and metadata only if called with Sumstats object
    if not is_dataframe:
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            liftover2_kwargs = {
                'chain_path': chain_path,
                'chrom': chrom,
                'pos': pos,
                'status': status,
                'to_build': to_build,
                'out_chrom_col': out_chrom_col,
                'out_pos_col': out_pos_col,
                'out_strand_col': out_strand_col,
                'one_based_input': one_based_input,
                'one_based_output': one_based_output,
                'remove': remove
            }
            _update_harmonize_step(sumstats_obj, "liftover2", liftover2_kwargs, True)
            sumstats_obj.meta["is_sorted"] = False
            sumstats_obj.meta["is_harmonised"] = False
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

