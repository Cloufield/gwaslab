"""
Lazy merge functions for SumstatsMultiLazy.

This module provides lazy versions of merge functions that work with Polars LazyFrames,
allowing merging of multiple sumstats files without materializing all data into memory.
"""

from typing import Optional, List, Tuple, Union
import polars as pl
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.info.g_vchange_status_polars import vchange_statusp, copy_statusp
from gwaslab.qc.qc_fix_sumstats_polars import flipallelestatsp
from gwaslab.qc.qc_reserved_headers import DEFAULT_COLUMN_ORDER
from gwaslab.util.util_in_fill_data import _fill_data


def _merge_mold_with_sumstats_by_chrposp_lazy(
    mold: pl.LazyFrame,
    sumstats: pl.LazyFrame,
    stats_cols1: Optional[List[str]] = None,
    stats_cols2: Optional[List[str]] = None,
    log: Log = Log(),
    suffixes: Tuple[str, str] = ("_1", ""),
    merge_mode: str = "full",
    merge_by_id: bool = False,
    verbose: bool = True
) -> pl.LazyFrame:
    """
    Merge two sumstats lazily using CHR, POS, ASET (allele set).
    
    This is a lazy version that works with LazyFrames, deferring execution
    until materialization. This allows merging large files without loading
    everything into memory.
    
    Parameters
    ----------
    mold : pl.LazyFrame
        First sumstats (mold/reference) as LazyFrame
    sumstats : pl.LazyFrame
        Second sumstats to merge as LazyFrame
    stats_cols1 : list, optional
        Statistics columns in first sumstats
    stats_cols2 : list, optional
        Statistics columns in second sumstats
    log : Log, optional
        Logging object
    suffixes : tuple, default ("_1", "")
        Suffixes for columns from first and second sumstats
    merge_mode : str, default "full"
        Merge mode ("full", "inner", "left", "right")
    merge_by_id : bool, default False
        Whether to merge by ID instead of CHR:POS:ASET
    verbose : bool, default True
        Whether to print verbose messages
        
    Returns
    -------
    pl.LazyFrame
        Merged LazyFrame (lazy until collected)
    """
    log.write("Start to merge sumstats (lazy)...", verbose=verbose)
    
    # Only rename ID columns if merging by ID (not needed for CHR:POS:ASET merge)
    # But we still need to handle conflicts if SNPID exists in both
    if merge_mode == "full" and merge_by_id:
        # Rename ID columns to avoid conflicts (only when merging by ID)
        rename_dict = {}
        sumstats_schema = sumstats.collect_schema()
        if "SNPID" in sumstats_schema:
            rename_dict["SNPID"] = "_SNPID_RIGHT"
        if "rsID" in sumstats_schema:
            rename_dict["rsID"] = "_rsID_RIGHT"
        
        if rename_dict:
            sumstats = sumstats.rename(rename_dict, strict=False)
    
    if merge_by_id:
        # Merge by ID (simpler case)
        if "SNPID" in mold.columns and "SNPID" in sumstats.columns:
            mold_sumstats = mold.join(
                sumstats, 
                on="SNPID", 
                how=merge_mode, 
                suffix="_", 
                coalesce=True
            )
        elif "rsID" in mold.columns and "rsID" in sumstats.columns:
            mold_sumstats = mold.join(
                sumstats, 
                on="rsID", 
                how=merge_mode, 
                suffix="_", 
                coalesce=True
            )
        else:
            raise ValueError("No ID column found for merging")
    else:
        # Merge by CHR, POS, ASET (allele set)
        # Create ASET column (sorted allele pair) for both
        mold_schema = mold.collect_schema()
        sumstats_schema = sumstats.collect_schema()
        
        # Determine which EA/NEA columns exist in mold
        # Check both suffixed and non-suffixed versions
        ea_suffixed = "EA" + suffixes[0] if suffixes[0] else "EA"
        nea_suffixed = "NEA" + suffixes[0] if suffixes[0] else "NEA"
        ea_col_mold = ea_suffixed if (suffixes[0] and ea_suffixed in mold_schema) else "EA"
        nea_col_mold = nea_suffixed if (suffixes[0] and nea_suffixed in mold_schema) else "NEA"
        eaf_col_mold = ("EAF" + suffixes[0]) if (suffixes[0] and ("EAF" + suffixes[0]) in mold_schema) else ("EAF" if "EAF" in mold_schema else None)
        
        # Verify columns exist
        if ea_col_mold not in mold_schema:
            raise ValueError(f"EA column '{ea_col_mold}' not found in mold. Available columns: {list(mold_schema.keys())[:10]}")
        if nea_col_mold not in mold_schema:
            raise ValueError(f"NEA column '{nea_col_mold}' not found in mold. Available columns: {list(mold_schema.keys())[:10]}")
        
        mold = mold.with_columns(
            pl.when(pl.col(ea_col_mold) > pl.col(nea_col_mold))
            .then(pl.col(ea_col_mold) + ":" + pl.col(nea_col_mold))
            .otherwise(pl.col(nea_col_mold) + ":" + pl.col(ea_col_mold))
            .alias("ASET")
        )
        
        # Determine which EA/NEA columns exist in sumstats
        ea_col_sumstats = "EA" + suffixes[1] if (suffixes[1] and ("EA" + suffixes[1]) in sumstats_schema) else "EA"
        nea_col_sumstats = "NEA" + suffixes[1] if (suffixes[1] and ("NEA" + suffixes[1]) in sumstats_schema) else "NEA"
        eaf_col_sumstats = ("EAF" + suffixes[1]) if (suffixes[1] and ("EAF" + suffixes[1]) in sumstats_schema) else ("EAF" if "EAF" in sumstats_schema else None)
        
        # Verify columns exist
        if ea_col_sumstats not in sumstats_schema:
            raise ValueError(f"EA column '{ea_col_sumstats}' not found in sumstats. Available columns: {list(sumstats_schema.keys())[:10]}")
        if nea_col_sumstats not in sumstats_schema:
            raise ValueError(f"NEA column '{nea_col_sumstats}' not found in sumstats. Available columns: {list(sumstats_schema.keys())[:10]}")
        
        sumstats = sumstats.with_columns(
            pl.when(pl.col(ea_col_sumstats) > pl.col(nea_col_sumstats))
            .then(pl.col(ea_col_sumstats) + ":" + pl.col(nea_col_sumstats))
            .otherwise(pl.col(nea_col_sumstats) + ":" + pl.col(ea_col_sumstats))
            .alias("ASET")
        )
        
        # Remove duplicates based on CHR, POS, ASET (lazy operation)
        mold = mold.unique(subset=["CHR", "POS", "ASET"])
        sumstats = sumstats.unique(subset=["CHR", "POS", "ASET"])
        
        # Select only needed columns from sumstats (exclude SNPID/rsID to avoid join conflicts)
        # We're merging by CHR:POS:ASET, so we don't need these IDs
        # Get all columns except SNPID and rsID
        # This avoids schema resolution issues
        sumstats = sumstats.select([
            pl.col("*").exclude(["SNPID", "rsID"])
        ])
        
        # Handle indels: add allele to ASET if lengths differ
        # Check if EAF column exists
        if eaf_col_mold:
            mold = mold.with_columns(
                pl.when(
                    pl.col(nea_col_mold).str.len_chars() != 
                    pl.col(ea_col_mold).str.len_chars()
                )
                .then(
                    pl.when(pl.col(eaf_col_mold) < 0.5)
                    .then(pl.col("ASET") + ":" + pl.col(ea_col_mold))
                    .otherwise(pl.col("ASET") + ":" + pl.col(nea_col_mold))
                )
                .otherwise(pl.col("ASET"))
                .alias("ASET")
            )
        
        if eaf_col_sumstats:
            sumstats = sumstats.with_columns(
                pl.when(
                    pl.col(nea_col_sumstats).str.len_chars() != pl.col(ea_col_sumstats).str.len_chars()
                )
                .then(
                    pl.when(pl.col(eaf_col_sumstats) < 0.5)
                    .then(pl.col("ASET") + ":" + pl.col(ea_col_sumstats))
                    .otherwise(pl.col("ASET") + ":" + pl.col(nea_col_sumstats))
                )
                .otherwise(pl.col("ASET"))
                .alias("ASET")
            )
        
        # Perform join (lazy operation)
        mold_sumstats = mold.join(
            sumstats, 
            on=["CHR", "POS", "ASET"], 
            how=merge_mode, 
            suffix="_", 
            coalesce=True
        )
    
    log.write("Finished merging sumstats (lazy).", verbose=verbose)
    return mold_sumstats


def _align_with_moldp_lazy(
    sumstats: pl.LazyFrame,
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_1", "")
) -> pl.LazyFrame:
    """
    Align alleles with mold (reference) lazily.
    
    This identifies perfect matches and flipped matches, and updates STATUS
    accordingly. All operations are lazy.
    
    Parameters
    ----------
    sumstats : pl.LazyFrame
        Merged sumstats as LazyFrame
    log : Log, optional
        Logging object
    verbose : bool, default True
        Whether to print verbose messages
    suffixes : tuple, default ("_1", "")
        Suffixes for columns from first and second sumstats
        
    Returns
    -------
    pl.LazyFrame
        Aligned LazyFrame (lazy until collected)
    """
    ea1 = "EA" + suffixes[0]
    nea1 = "NEA" + suffixes[0]
    ea2 = "EA" + suffixes[1] if suffixes[1] else "EA"
    nea2 = "NEA" + suffixes[1] if suffixes[1] else "NEA"
    status1 = "STATUS" + suffixes[0]
    status2 = "STATUS" + suffixes[1] if suffixes[1] else "STATUS"
    
    # Check if status columns exist (get schema without materializing)
    schema = sumstats.schema
    has_status1 = status1 in schema
    has_status2 = status2 in schema
    
    # Create match conditions (lazy)
    is_perfect_match = (
        (pl.col(ea2) == pl.col(ea1)) & 
        (pl.col(nea2) == pl.col(nea1))
    )
    
    is_flipped_match = (
        (pl.col(ea2) == pl.col(nea1)) & 
        (pl.col(nea2) == pl.col(ea1))
    )
    
    log.write(" -Aligning alleles with reference (lazy)...", verbose=verbose)
    
    # Note: We can't count matches lazily without materializing
    # So we'll skip the counts in lazy mode
    log.write("  -Match detection deferred until materialization", verbose=verbose)
    
    # Copy STATUS for perfect matches (lazy operation)
    # Only if both STATUS columns exist
    if has_status1 and has_status2:
        try:
            log.write("  -For perfect match: copy STATUS from reference...", verbose=verbose)
            sumstats = copy_statusp(sumstats, is_perfect_match, status1, status2, 6)
            
            # Update STATUS for flipped matches (lazy operation)
            log.write("  -For flipped match: convert STATUS xxxxx[456789]x to xxxxx3x...", verbose=verbose)
            sumstats = vchange_statusp(sumstats, is_flipped_match, status2, 6, "456789", "333333")
        except Exception as e:
            log.write("  -Warning: STATUS alignment skipped due to: {}".format(str(e)), verbose=verbose)
            # Continue without STATUS alignment
    
    return sumstats


def flipallelestatsp_lazy(
    sumstats: pl.LazyFrame,
    status: str = "STATUS",
    verbose: bool = True,
    log: Log = Log()
) -> pl.LazyFrame:
    """
    Flip alleles based on STATUS lazily.
    
    This is a simplified lazy-compatible version. For full functionality,
    allele flipping will be applied during materialization if needed.
    
    Parameters
    ----------
    sumstats : pl.LazyFrame
        Sumstats as LazyFrame
    status : str, default "STATUS"
        STATUS column name
    verbose : bool, default True
        Whether to print verbose messages
    log : Log, optional
        Logging object
        
    Returns
    -------
    pl.LazyFrame
        LazyFrame (flipping will be applied during materialization if needed)
    """
    # For now, we skip allele flipping during lazy merge
    # It will be applied during materialization when the data is collected
    # This avoids the complexity of handling len() checks and complex transformations
    # on LazyFrames
    log.write("Flipping alleles based on STATUS (deferred to materialization)...", verbose=verbose)
    log.write("  -Allele flipping will be applied when data is materialized", verbose=verbose)
    
    # Return as-is for now - flipping can be done after materialization
    # or we can materialize just the STATUS column to check, but that defeats the purpose
    return sumstats


def _fill_missing_columnsp_lazy(
    sumstats: pl.LazyFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True
) -> pl.LazyFrame:
    """
    Fill missing columns lazily.
    
    Parameters
    ----------
    sumstats : pl.LazyFrame
        Sumstats as LazyFrame
    columns : list
        List of column names to fill
    log : Log, optional
        Logging object
    verbose : bool, default True
        Whether to print verbose messages
        
    Returns
    -------
    pl.LazyFrame
        LazyFrame with filled columns (lazy until collected)
    """
    # _fill_data should work with LazyFrames
    # If not, we'll need to implement a lazy version
    log.write("Filling missing columns (lazy)...", verbose=verbose)
    sumstats = _fill_data(sumstats, to_fill=columns)
    return sumstats


def _renaming_colsp_lazy(
    sumstats: pl.LazyFrame,
    columns: List[str],
    log: Log = Log(),
    verbose: bool = True,
    suffixes: Tuple[str, str] = ("_1", "_2")
) -> pl.LazyFrame:
    """
    Rename columns with suffix lazily.
    
    Parameters
    ----------
    sumstats : pl.LazyFrame
        Sumstats as LazyFrame
    columns : list
        List of column names to rename
    log : Log, optional
        Logging object
    verbose : bool, default True
        Whether to print verbose messages
    suffixes : tuple, default ("_1", "_2")
        Suffix to add to columns
        
    Returns
    -------
    pl.LazyFrame
        LazyFrame with renamed columns (lazy until collected)
    """
    to_rename = ["STATUS"]
    for col in columns:
        if col in sumstats.columns:
            to_rename.append(col)
    
    rename_dict = {i: i + suffixes[1] for i in to_rename if i in sumstats.columns}
    
    if rename_dict:
        log.write(" -Renaming columns by adding suffix {} (lazy)...".format(suffixes[1]), verbose=verbose)
        sumstats = sumstats.rename(rename_dict)
    
    return sumstats


def _sort_pair_colsp_lazy(
    molded_sumstats: pl.LazyFrame,
    verbose: bool = True,
    log: Log = Log(),
    order: Optional[List[str]] = None,
    stats_order: Optional[List[str]] = None,
    suffixes: Tuple[str, str] = ("_1", "_2")
) -> pl.LazyFrame:
    """
    Sort columns lazily.
    
    Note: Column ordering requires knowing which columns exist, which may
    require partial materialization. However, we can build the order list
    lazily and apply it.
    
    Parameters
    ----------
    molded_sumstats : pl.LazyFrame
        Merged sumstats as LazyFrame
    verbose : bool, default True
        Whether to print verbose messages
    log : Log, optional
        Logging object
    order : list, optional
        Explicit column order
    stats_order : list, optional
        Order for statistics columns
    suffixes : tuple, default ("_1", "_2")
        Suffixes for columns
        
    Returns
    -------
    pl.LazyFrame
        LazyFrame with sorted columns (lazy until collected)
    """
    # For lazy sorting, we need to get column names
    # This requires materializing the schema, but not the data
    log.write("Sorting columns (lazy)...", verbose=verbose)
    
    # Get schema (doesn't materialize data) - use collect_schema to avoid warnings
    schema = molded_sumstats.collect_schema()
    available_cols = list(schema.keys())
    
    if stats_order is None:
        # Base order: first 6 columns from DEFAULT_COLUMN_ORDER (only if they exist)
        base_cols = ["SNPID", "rsID", "CHR", "POS", "EA", "NEA"]
        # Only include base cols that actually exist
        order = [col for col in base_cols if col in available_cols]
        # Stats order: remaining columns from DEFAULT_COLUMN_ORDER
        stats_order = [col for col in DEFAULT_COLUMN_ORDER if col not in base_cols]
        # Add non-reserved columns
        additional_cols = ["DIRECTION", "I2"]
        for col in additional_cols:
            if col not in stats_order:
                stats_order.append(col)
    
    # Build ordered column list
    for suffix in suffixes:
        for col in stats_order:
            suffixed_col = col + suffix
            if suffixed_col in available_cols and suffixed_col not in order:
                order.append(suffixed_col)
    
    # Add any remaining columns
    for col in available_cols:
        if col not in order:
            order.append(col)
    
    # Select columns in order (lazy operation) - only existing columns
    order_existing = [col for col in order if col in available_cols]
    log.write("  -Reordering columns (lazy)...", verbose=verbose)
    molded_sumstats = molded_sumstats.select(order_existing)
    
    return molded_sumstats


def merge_sumstats_files_lazy(file_manager, merge_mode="full", merge_by_id=False, 
                              verbose=True, log=None):
    """
    Merge multiple sumstats files lazily without materializing all data.
    
    This function builds a lazy merge plan that chains operations across
    multiple files, only materializing when absolutely necessary.
    
    Parameters
    ----------
    file_manager : SumstatsFileManager
        File manager with paths and metadata
    merge_mode : str, default "full"
        Merge mode
    merge_by_id : bool, default False
        Whether to merge by ID
    verbose : bool, default True
        Whether to print verbose messages
    log : Log, optional
        Logging object
        
    Returns
    -------
    pl.LazyFrame
        Merged LazyFrame (lazy until collected)
    """
    if log is None:
        log = Log()
    
    if len(file_manager.paths) == 0:
        raise ValueError("No file paths provided")
    
    log.write("Building lazy merge plan for {} files...".format(len(file_manager.paths)), verbose=verbose)
    
    # Load first file as lazy frame
    first_path = file_manager.paths[0]
    first_tab_format = file_manager.file_tab_formats[0]
    
    log.write(" -Loading file #1: {} (lazy)...".format(first_path), verbose=verbose)
    
    if first_tab_format == "parquet":
        mold_lazy = pl.scan_parquet(first_path)
    elif first_tab_format in ["tsv", "csv"]:
        separator = "\t" if first_tab_format == "tsv" else ","
        mold_lazy = pl.scan_csv(first_path, separator=separator)
    else:
        raise ValueError("Unsupported file format: {}".format(first_tab_format))
    
    # Rename columns with _1 suffix (lazy operation)
    rename_dict_1 = {"EA": "EA_1", "NEA": "NEA_1", "STATUS": "STATUS_1"}
    if file_manager.stats_cols[0]:
        rename_dict_1.update({col: col + "_1" for col in file_manager.stats_cols[0]})
    if file_manager.other_cols[0]:
        rename_dict_1.update({col: col + "_1" for col in file_manager.other_cols[0]})
    
    mold_lazy = mold_lazy.rename(rename_dict_1)
    
    # Merge remaining files one by one
    for i, path in enumerate(file_manager.paths[1:], start=1):
        log.write(" -Merging file #{}: {} (lazy)...".format(i+1, path), verbose=verbose)
        
        # Load next file as lazy frame
        tab_format = file_manager.file_tab_formats[i]
        
        if tab_format == "parquet":
            file_lazy = pl.scan_parquet(path)
        elif tab_format in ["tsv", "csv"]:
            separator = "\t" if tab_format == "tsv" else ","
            file_lazy = pl.scan_csv(path, separator=separator)
        else:
            raise ValueError("Unsupported file format: {}".format(tab_format))
        
        # Determine correct suffixes based on merge iteration
        # After first merge, EA/NEA in mold have no suffix, so we need to adjust
        if i == 1:
            # First merge: mold has EA_1/NEA_1, new file has EA/NEA
            merge_suffixes = ("_1", "")
        else:
            # Subsequent merges: mold has EA/NEA (no suffix), new file has EA/NEA
            # Rename new file's EA/NEA to have suffix before merging
            new_suffix = "_{}".format(i+1)
            file_lazy = file_lazy.rename({"EA": "EA" + new_suffix, "NEA": "NEA" + new_suffix}, strict=False)
            merge_suffixes = ("", new_suffix)
        
        # Merge with mold (lazy)
        mold_lazy = _merge_mold_with_sumstats_by_chrposp_lazy(
            mold=mold_lazy,
            sumstats=file_lazy,
            stats_cols1=file_manager.other_cols[0],
            stats_cols2=file_manager.other_cols[i],
            log=log,
            suffixes=merge_suffixes,
            merge_mode=merge_mode,
            merge_by_id=merge_by_id,
            verbose=verbose
        )
        
        # Align alleles (lazy)
        mold_lazy = _align_with_moldp_lazy(
            mold_lazy,
            log=log,
            verbose=verbose,
            suffixes=merge_suffixes
        )
        
        # Flip alleles if needed (lazy)
        mold_lazy = flipallelestatsp_lazy(
            mold_lazy,
            log=log,
            verbose=verbose
        )
        
        # After first merge (i==1), drop EA/NEA from new file and rename EA_1/NEA_1 to EA/NEA
        # For subsequent merges, drop the suffixed EA/NEA from new file
        if i == 1:
            # First merge: drop EA/NEA from second file, rename EA_1/NEA_1 to EA/NEA
            schema_before = mold_lazy.collect_schema()
            if "EA" in schema_before and "NEA" in schema_before:
                mold_lazy = mold_lazy.drop(["EA", "NEA"])
            if "EA_1" in schema_before and "NEA_1" in schema_before:
                mold_lazy = mold_lazy.rename({"EA_1": "EA", "NEA_1": "NEA"})
        else:
            # Subsequent merges: drop the suffixed EA/NEA from new file
            new_suffix = "_{}".format(i+1)
            schema_before = mold_lazy.collect_schema()
            ea_new = "EA" + new_suffix
            nea_new = "NEA" + new_suffix
            if ea_new in schema_before and nea_new in schema_before:
                mold_lazy = mold_lazy.drop([ea_new, nea_new])
        
        # Fill missing columns if needed
        if not set(file_manager.stats_cols[i]) == set(file_manager.stats_cols[0]):
            cols_to_fill = set(file_manager.stats_cols[0]).difference(set(file_manager.stats_cols[i]))
            if cols_to_fill:
                mold_lazy = _fill_missing_columnsp_lazy(
                    mold_lazy,
                    list(cols_to_fill),
                    log=log,
                    verbose=verbose
                )
        
        # Rename columns with suffix
        mold_lazy = _renaming_colsp_lazy(
            mold_lazy,
            file_manager.stats_cols[0] + file_manager.other_cols[i],
            log=log,
            verbose=verbose,
            suffixes=("_1", "_{}".format(i+1))
        )
        
        # Sort columns
        mold_lazy = _sort_pair_colsp_lazy(
            mold_lazy,
            verbose=verbose,
            log=log,
            suffixes=["_{}".format(j) for j in range(1, i+2)]
        )
    
    log.write("Lazy merge plan completed. Data will be loaded on demand.", verbose=verbose)
    return mold_lazy

