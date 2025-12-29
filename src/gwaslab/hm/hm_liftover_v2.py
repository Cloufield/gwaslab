from __future__ import annotations

from typing import List

import pandas as pd

from gwaslab.info.g_Log import Log
from gwaslab.info.g_vchange_status import status_match, ensure_status_int, vchange_status
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_common_data import get_chain
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from typing import Optional
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_fix_sumstats import _fix_chr, _fix_pos

try:
    from sumstats_liftover import liftover_df
except ImportError:
    raise ImportError(
        "sumstats_liftover package is required. Install it with: pip install sumstats-liftover"
    )


def _normalize_chrom_name_vectorized(chrom_series: pd.Series, mapper: Optional[ChromosomeMapper] = None) -> pd.Series:
    """
    Vectorized version of _normalize_chrom_name for processing pandas Series.
    Uses ChromosomeMapper for normalization.
    
    Parameters
    ----------
    chrom_series : pd.Series
        Series of chromosome names (can be strings, ints, or mixed)
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use. If None, creates a new one.
    
    Returns
    -------
    pd.Series
        Series of normalized chromosome names (string format for X/Y/MT, numeric for autosomes)
    """
    # Use provided mapper or create a new one
    if mapper is None:
        mapper = ChromosomeMapper()
    else:
        # Use the provided mapper (it already has species/build info)
        mapper = ChromosomeMapper(
            species=mapper.species,
            build=mapper.build,
            log=mapper.log,
            verbose=mapper.verbose
        )
    # Detect format and convert to string
    mapper.detect_sumstats_format(chrom_series)
    return mapper.to_string(chrom_series)


# ----------------------------
# Helper functions
# ----------------------------

def _df_split(dataframe: pd.DataFrame, n: int) -> List[pd.DataFrame]:
    """Split a DataFrame into n approximately equal parts for parallel processing."""
    k, m = divmod(len(dataframe), n)
    return [dataframe.iloc[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]


# ----------------------------
# Fast liftover functions (using sumstats_liftover package)
# ----------------------------


@with_logging(
        start_to_msg= "perform liftover",
        finished_msg= "liftover",
        start_function= ".liftover()",
        start_cols=["CHR","POS","STATUS"]
)
def _liftover_variant(sumstats_obj,
                      chain_path: str = None,
                      from_build: str = None,
                      to_build: str = None,
                      chrom: str = "CHR",
                      pos: str = "POS",
                      status: str = "STATUS",
                      out_chrom_col: str = "CHR_LIFT",
                      out_pos_col: str = "POS_LIFT",
                      out_strand_col: str = "STRAND_LIFT",
                      one_based_input: bool = True,
                      one_based_output: bool = True,
                      remove: bool = True,
                      filter_by_status: bool = True,
                      verbose: bool = True,
                      log: Log = Log()) -> pd.DataFrame:
    """
    Perform liftover of variants to a new genome build.
    
    Can be called with either:
    - from_build and to_build (chain file will be automatically found)
    - chain_path (direct path to chain file)
    
    This is a fast chain-based liftover implementation that directly parses
    UCSC chain files and performs vectorized coordinate conversion. It converts
    genomic coordinates from one genome build (e.g., hg19/GRCh37) to another
    (e.g., hg38/GRCh38) using UCSC chain files.
    
    Chain files are automatically obtained from built-in data (for hg19<->hg38)
    or downloaded from UCSC if not available. Built-in chain files are preferred
    for better performance and reliability.
    
    Parameters
    ----------
    chain_path : str, optional
        Path to UCSC chain file. If provided, from_build and to_build are optional.
    from_build : str, optional
        Source genome build (e.g., "19"). If None, uses sumstats_obj.build if available.
    to_build : str, optional
        Target genome build (e.g., "38"). Required if chain_path is not provided.
    remove : bool, default=True
        Whether to remove unmapped variants
    verbose : bool, default=True
        Whether to print progress messages
    
    Returns
    -------
    pd.DataFrame
        DataFrame with lifted coordinates (or updated sumstats_obj.data if Sumstats object).
        The returned DataFrame contains:
        - Updated CHR and POS columns with lifted coordinates (from the target build)
        - Updated STATUS column with new status codes reflecting the liftover results
        - Unmapped variants are either removed (if remove=True) or kept with NA values
          for CHR and POS (if remove=False)
    
    Notes
    -----
    If called with a Sumstats object, the function will update the object's data
    and metadata in place. If called with a DataFrame, it returns a new DataFrame.
    """
    import pandas as pd
    from pathlib import Path
    
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_obj, pd.DataFrame):
        sumstats = sumstats_obj
        is_dataframe = True
    else:
        sumstats = sumstats_obj.data
        is_dataframe = False
    
    if len(sumstats) == 0:
        log.write(" -Empty dataframe, skipping liftover.", verbose=verbose)
        return sumstats if is_dataframe else sumstats_obj.data
    
    # Get chain file path if from_build/to_build provided instead
    if chain_path is None:
        if from_build is None or to_build is None:
            raise ValueError("Either chain_path or both from_build and to_build must be provided")
        try:
            chain_path = get_chain(from_build=from_build, to_build=to_build)
            if chain_path is None or chain_path == False:
                lifter_from_build = _process_build(from_build, log=log, verbose=False)
                lifter_to_build = _process_build(to_build, log=log, verbose=False)
                raise ValueError("No available chain file for {} -> {}".format(lifter_from_build, lifter_to_build))
            
            # Check if using built-in chain file
            builtin_chains_dir = Path(__file__).parents[1] / "data" / "chains"
            is_builtin = Path(chain_path).parent == builtin_chains_dir
            
            if is_builtin:
                log.log_operation("Using built-in chain file: {}".format(chain_path), verbose=verbose)
            else:
                log.log_operation("Using chain file: {}".format(chain_path), verbose=verbose)
        except Exception as e:
            lifter_from_build = _process_build(from_build, log=log, verbose=False)
            lifter_to_build = _process_build(to_build, log=log, verbose=False)
            raise ValueError("No available chain file for {} -> {}. Error: {}".format(lifter_from_build, lifter_to_build, str(e)))
    else:
        log.log_operation("Using provided chain file: {}".format(chain_path), verbose=verbose)
    
    # Filter by status if requested
    # Only copy if we might need to merge later (avoid unnecessary copy overhead)
    original_sumstats = None
    to_lift = None
    needs_merging = False
    if filter_by_status and status in sumstats.columns:
        # digit 4 = 0
        to_lift = status_match(sumstats[status], 4, [0])
        to_lift_count = to_lift.sum()
        
        if to_lift_count == 0:
            log.log_operation("No variants with status code xxx0xxx to liftover", verbose=verbose)
            # Update metadata if Sumstats object
            if not is_dataframe:
                try:
                    from gwaslab.info.g_meta import _update_harmonize_step
                    sumstats_obj.meta["is_sorted"] = False
                    sumstats_obj.meta["is_harmonised"] = False
                    if to_build is not None:
                        # Use setter to ensure consistency between self.build and meta
                        sumstats_obj.build = to_build
                    liftover_kwargs = {
                        'from_build': from_build, 'to_build': to_build, 'remove': remove, 'chain_path': chain_path
                    }
                    _update_harmonize_step(sumstats_obj, "liftover", liftover_kwargs, True)
                except Exception:
                    pass
            return sumstats_obj.data if not is_dataframe else sumstats
        
        # Only copy when we know we need to merge (reduce memory overhead)
        original_sumstats = sumstats.copy()
        # Use view instead of copy for filtering (more efficient)
        sumstats_to_lift = sumstats.loc[to_lift, :]
        log.log_operation("Converting variants with status code xxx0xxx: {:,}".format(to_lift_count), verbose=verbose)
        # Create a copy only when we modify it (lazy copy)
        sumstats = sumstats_to_lift.copy()
        needs_merging = True
    
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
    
    log.write(f" -Target build: {to_build_processed}", verbose=verbose)
    log.write(f" -Input positions are {'1-based' if one_based_input else '0-based'}", verbose=verbose)
    log.write(f" -Output positions will be {'1-based' if one_based_output else '0-based'}", verbose=verbose)
    
    # Get mapper from Sumstats object if available
    mapper = None
    if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
        mapper = sumstats_obj.mapper
    
    # Store original chromosome values for comparison (normalized for matching)
    # This is needed to detect chromosome mismatches (variants mapped to different chromosomes)
    # Use vectorized normalization for better performance
    original_chrom_normalized = _normalize_chrom_name_vectorized(sumstats[chrom], mapper=mapper)
    
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
        # Normalize lifted chromosomes for comparison (only for non-NA values)
        # Use vectorized normalization for better performance
        # First filter out NA and -1 values, then normalize
        valid_mask = sumstats[out_chrom_col].notna() & (sumstats[out_chrom_col].astype(str) != '-1')
        lifted_chrom_normalized = pd.Series(index=sumstats.index, dtype=object)
        if valid_mask.any():
            lifted_chrom_normalized.loc[valid_mask] = _normalize_chrom_name_vectorized(sumstats.loc[valid_mask, out_chrom_col], mapper=mapper)
        
        # Mark variants with chromosome mismatches as unmapped
        # Compare normalized chromosome names (handles "1" vs "chr1", "X" vs "23", etc.)
        chrom_mismatch = lifted_chrom_normalized.notna() & (lifted_chrom_normalized != original_chrom_normalized)
        chrom_mismatch_count = chrom_mismatch.sum()
        
        if chrom_mismatch_count > 0:
            log.write(f" -Chromosome mismatches detected: {chrom_mismatch_count} variants (treated as unmapped)", verbose=verbose)
            
            # Show examples of chromosome mismatches
            if verbose:
                chrom_mismatch_variants = sumstats[chrom_mismatch].head(5)  # Show up to 5 examples
                example_cols = []
                # Add SNPID if available
                if 'SNPID' in sumstats.columns:
                    example_cols.append('SNPID')
                # Add original CHR and POS
                example_cols.extend([chrom, pos])
                # Add lifted CHR and POS
                if out_chrom_col in sumstats.columns:
                    example_cols.append(out_chrom_col)
                if out_pos_col in sumstats.columns:
                    example_cols.append(out_pos_col)
                # Add REF/ALT if available
                if 'REF' in sumstats.columns:
                    example_cols.append('REF')
                if 'ALT' in sumstats.columns:
                    example_cols.append('ALT')
                # Add STATUS if available
                if status in sumstats.columns:
                    example_cols.append(status)
                
                # Only show columns that exist
                example_cols = [col for col in example_cols if col in chrom_mismatch_variants.columns]
                
                if len(example_cols) > 0:
                    examples = chrom_mismatch_variants[example_cols]
                    log.write(" -Examples of chromosome mismatches:", verbose=verbose)
                    for idx, row in examples.iterrows():
                        example_str = "   " + " | ".join([f"{col}={row[col]}" for col in example_cols])
                        log.write(example_str, verbose=verbose)
            
            unmapped = unmapped | chrom_mismatch
    
    unmapped_count = unmapped.sum()
    mapped_count = len(sumstats) - unmapped_count
    mapped_mask = ~unmapped  # Compute mapped_mask early for use in status updates
    
    log.write(f" -Mapped: {mapped_count} variants", verbose=verbose)
    log.write(f" -Unmapped: {unmapped_count} variants", verbose=verbose)
    
    # Show examples of unmapped variants
    if unmapped_count > 0 and verbose:
        unmapped_variants = sumstats[unmapped].head(5)  # Show up to 5 examples
        example_cols = [chrom, pos]
        # Add SNPID if available
        if 'SNPID' in sumstats.columns:
            example_cols.insert(0, 'SNPID')
        # Add REF/ALT if available
        if 'REF' in sumstats.columns:
            example_cols.append('REF')
        if 'ALT' in sumstats.columns:
            example_cols.append('ALT')
        # Add STATUS if available
        if status in sumstats.columns:
            example_cols.append(status)
        
        # Only show columns that exist
        example_cols = [col for col in example_cols if col in unmapped_variants.columns]
        
        if len(example_cols) > 0:
            examples = unmapped_variants[example_cols]
            log.write(" -Examples of unmapped variants:", verbose=verbose)
            for idx, row in examples.iterrows():
                example_str = "   " + " | ".join([f"{col}={row[col]}" for col in example_cols])
                log.write(example_str, verbose=verbose)
    
    # Update status codes after liftover conversion
    # Set digit 4 to "9" (unchecked) for mapped variants since coordinates are from liftover
    # fix_chr and fix_pos will later validate and update digit 4 appropriately
    if status in sumstats.columns:
        # Ensure STATUS column is integer type before processing
        sumstats = ensure_status_int(sumstats, status)
        
        # Extract status digits for all variants (vectorized operation)
        status_vals = sumstats[status].values
        # Extract digit3 and digit5: digit3 = (status // 10000) % 10, digit5 = (status // 100) % 10
        digit3 = (status_vals // 10000) % 10
        digit5 = (status_vals // 100) % 10
        # Build status_end_ints with digit 4 = 9 (unchecked) since coordinates are from liftover
        # fix_chr and fix_pos will validate and update digit 4 later
        status_end_ints = digit3 * 10000 + 9000 + digit5 * 100 + 99
        
        # Update status for mapped variants: to_build * 100000 + status_end_int (digit 4 = 9)
        # mapped_mask already computed above
        if mapped_count > 0:
            new_status_mapped = (int(to_build_processed) * 100000 + status_end_ints[mapped_mask]).astype('Int64')
            sumstats.loc[mapped_mask, status] = new_status_mapped
        
        # Update status for unmapped variants: 9700000 + status_end_int (build 97 = unmapped)
        if unmapped_count > 0:
            new_status_unmapped = (9700000 + status_end_ints[unmapped]).astype('Int64')
            sumstats.loc[unmapped, status] = new_status_unmapped
    # mapped_mask already computed above for use in CHR/POS updates
    
    # Overwrite CHR and POS with lifted values, drop temporary columns
    # For unmapped variants, set CHR and POS to NA (mimic original liftover behavior)
    # Combine operations to reduce overhead
    if out_chrom_col in sumstats.columns and chrom in sumstats.columns:
        # For chromosomes, don't use .values to allow pandas to handle type conversion
        # (chromosomes can be strings like 'X', 'Y', 'M' or integers)
        # Convert column to object dtype first if it's categorical or numeric and source has object dtype
        if mapped_mask.any():
            # If target is categorical or numeric but source is object (contains strings), convert to object
            if (isinstance(sumstats[chrom].dtype, pd.CategoricalDtype) or 
                (sumstats[chrom].dtype != 'object' and sumstats[out_chrom_col].dtype == 'object')):
                sumstats[chrom] = sumstats[chrom].astype('object')
            sumstats.loc[mapped_mask, chrom] = sumstats.loc[mapped_mask, out_chrom_col]
        if unmapped_count > 0 and not remove:
            sumstats.loc[unmapped, chrom] = pd.NA
        sumstats = sumstats.drop(columns=[out_chrom_col])
    
    if out_pos_col in sumstats.columns and pos in sumstats.columns:
        # For positions, use .values for better performance (always numeric)
        if mapped_mask.any():
            sumstats.loc[mapped_mask, pos] = sumstats.loc[mapped_mask, out_pos_col].values
        if unmapped_count > 0 and not remove:
            sumstats.loc[unmapped, pos] = pd.NA
        sumstats = sumstats.drop(columns=[out_pos_col])
    
    if out_strand_col in sumstats.columns:
        sumstats = sumstats.drop(columns=[out_strand_col])
    
    # Remove unmapped variants if requested (after setting NA values)
    if remove and unmapped_count > 0:
        sumstats = sumstats[~unmapped].copy()
        log.write(f" -Removed {unmapped_count} unmapped variants", verbose=verbose)
    
    # Merge results back if we filtered by status
    if needs_merging:
        if remove:
            # If removing unmapped, concatenate non-lifted variants with lifted (unmapped removed)
            non_lifted = original_sumstats.loc[~to_lift, :]
            if len(non_lifted) > 0:
                sumstats = pd.concat([non_lifted, sumstats], ignore_index=False)
            # else: sumstats already contains only lifted variants
        else:
            # If keeping unmapped, update the original dataframe with lifted coordinates
            original_sumstats.loc[to_lift, chrom] = sumstats[chrom]
            original_sumstats.loc[to_lift, pos] = sumstats[pos]
            if status in sumstats.columns:
                original_sumstats.loc[to_lift, status] = sumstats[status]
            sumstats = original_sumstats
    
    # After liftover and merging, validate and fix chr and pos
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
        try:
            from gwaslab.info.g_meta import _update_harmonize_step
            if to_build is not None:
                # Use setter to ensure consistency between self.build and meta
                sumstats_obj.build = to_build
            sumstats_obj.meta["is_sorted"] = False
            sumstats_obj.meta["is_harmonised"] = False
            liftover_kwargs = {
                'from_build': from_build, 'to_build': to_build, 'remove': remove, 'chain_path': chain_path
            }
            _update_harmonize_step(sumstats_obj, "liftover", liftover_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

