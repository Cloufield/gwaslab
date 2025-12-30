"""
Module for assigning rsIDs and annotating GWAS summary statistics from reference files.

This module provides functionality to:
- Assign rsIDs to variants in GWAS summary statistics by matching against reference VCF/BCF/TSV files
- Annotate summary statistics with additional fields (e.g., rsID, AF) from reference lookup tables
- Extract lookup tables from VCF/BCF files for efficient annotation
- Convert VCF files to BCF format for improved performance
- Handle allele-aware matching and STATUS-based filtering for variant assignment

Key Functions:
- _assign_rsid(): Main function for assigning rsIDs with allele matching and overwrite control
- _annotate_sumstats(): General annotation function for assigning multiple fields from lookup tables
- _extract_lookup_table_from_vcf_bcf(): Extract variant information from VCF/BCF to create lookup tables
- _assign_from_lookup(): Core function that performs the actual assignment using lookup tables
- _convert_vcf_to_bcf(): Convert VCF files to BCF format for faster processing

The module supports multiple reference formats (VCF, BCF, TSV), handles chromosome name mapping,
and provides flexible overwrite modes for existing rsID values. It uses parallel processing
where applicable and includes comprehensive logging and error handling.
"""

import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
import subprocess
import shutil
import os
import tempfile
from gwaslab.io.io_vcf import is_vcf_file
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

@with_logging(
    start_to_msg="assign rsID from reference",
    finished_msg="assigning rsID from reference",
    start_cols=["CHR","POS","EA","NEA","STATUS"],
    start_function=".assign_rsid2()"
)
def _assign_rsid(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    lookup_path: str | None = None,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    threads: int = 6,
    rsid: str = "rsID",
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    overwrite: str = "empty",
    mapper: ChromosomeMapper | None = None,
    log: "Log" = Log(),
    verbose: bool = True,
):
    """
    Assign rsIDs to GWAS summary statistics using reference data with allele matching and STATUS filtering.

    This function assigns rsIDs to a GWAS summary statistics DataFrame by matching variants against a reference
    VCF or TSV file. It performs allele-aware matching and applies STATUS-based filtering to determine which
    variants should be assigned rsIDs. The function handles various reference formats and allows control over
    overwrite behavior for existing rsID values.

    Parameters
    ----------
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object.
    path : str or None, optional
        Path to reference file (VCF/BCF or TSV). Overrides `tsv_path`.
    vcf_path : str or None, optional
        Path to VCF/BCF file. Overrides `path` and `tsv_path`.
    tsv_path : str or None, optional
        Path to precomputed lookup TSV file. If not provided, generated from VCF.
    reuse_lookup : bool, optional
        If True, reuse existing lookup TSV if available.
    convert_to_bcf : bool, optional
        If True, convert VCF to BCF before processing.
    strip_info : bool, optional
        If True, strip INFO fields when converting VCF to BCF.
    threads : int, optional
        Number of threads for bcftools operations.
    overwrite : str, optional
        Overwrite mode: "all", "invalid", or "empty". Determines which existing rsID values to overwrite.
        Default is "empty".
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper with automatic format detection.
    log : gwaslab.g_Log.Log, optional
        Log object for recording progress. Default is a new Log instance.
    verbose : bool, optional
        If True, log detailed progress messages. Default is True.

    Returns
    -------
    pd.DataFrame
        The input `sumstats` DataFrame with rsID column updated.

    Raises
    ------
    ValueError
        If required columns are missing in `sumstats` or invalid `overwrite` value is provided.
    FileNotFoundError
        If specified reference file is not found.

    Notes
    -----
    - The function first checks for required columns in `sumstats`.
    - STATUS filtering uses a regex pattern to identify variants eligible for rsID assignment.
    - Overwrite modes:
      * "all": overwrite all rsIDs for eligible variants
      * "invalid": overwrite only non-rsID formatted values (e.g., not matching "rs[0-9]+")
      * "empty": only fill missing rsID values
    - If `ref_mode` is "auto", the function determines whether to use VCF or TSV based on file extension.
    """
    import os
    import re

    # Handle both DataFrame and Sumstats object
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
            mapper = sumstats_obj.mapper
        else:
            # Fix: safely get species and build
            if not is_dataframe and hasattr(sumstats_obj, 'meta'):
                species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens")
            else:
                species = "homo sapiens"
            if not is_dataframe and hasattr(sumstats_obj, 'build'):
                build = sumstats_obj.build
            else:
                build = None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chrom in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chrom])

    # Ensure rsID exists
    if rsid not in sumstats.columns:
        sumstats[rsid] = pd.Series(pd.NA, dtype="string")

    total_before = len(sumstats)
    had_rsid_before = (~sumstats[rsid].isna()).sum()

    # ---------------------------
    # STATUS filter - optimized: compute once and reuse
    # ---------------------------
    # Match: digit 4 is 0 and digit 5 is 0-4
    from gwaslab.info.g_vchange_status import status_match
    # Optimize: extract STATUS once and compute both conditions
    status_series = sumstats["STATUS"]
    digit4_match = status_match(status_series, 4, [0])
    digit5_match = status_match(status_series, 5, [0, 1, 2, 3, 4])
    standardized_normalized = digit4_match & digit5_match
    to_assign_mask = standardized_normalized

    # ---------------------------
    # Overwrite modes - optimized: early exit if no variants to assign
    # ---------------------------
    if overwrite == "all":
        sumstats.loc[to_assign_mask, rsid] = pd.NA
        # After setting to NA, count how many eligible variants need assignment
        variants_to_assign = (to_assign_mask & sumstats[rsid].isna()).sum()
    elif overwrite == "invalid":
        # Optimize: only convert to string for non-NA values that need checking
        rsid_series = sumstats[rsid]
        # For invalid check, only need to check non-NA values
        non_na_mask = rsid_series.notna()
        if non_na_mask.any():
            # Only convert non-NA values to string for regex matching
            rsid_str = rsid_series[non_na_mask].astype("string")
            invalid_mask = pd.Series(False, index=sumstats.index)
            invalid_mask[non_na_mask] = ~rsid_str.str.match(r"^rs[0-9]+$", na=False)
        else:
            invalid_mask = pd.Series(False, index=sumstats.index)
        sumstats.loc[to_assign_mask & invalid_mask, rsid] = pd.NA
        # After setting invalid ones to NA, count how many eligible variants need assignment
        variants_to_assign = (to_assign_mask & sumstats[rsid].isna()).sum()
    elif overwrite == "empty":
        # only fill missing → check which variants need assignment
        variants_to_assign = (to_assign_mask & sumstats[rsid].isna()).sum()
    else:
        raise ValueError("overwrite must be 'all', 'invalid', or 'empty'")

    # Early exit if no variants need assignment
    if variants_to_assign == 0:
        log.write(" -No variants need rsID assignment, skipping lookup...", verbose=verbose)
        had_rsid_after = (~sumstats[rsid].isna()).sum()
        log.write(" -rsID count: {} → {} / {}...".format(had_rsid_before, had_rsid_after, total_before), verbose=verbose)
        log.write(" -Finished assigning rsID from reference.", verbose=verbose)
        
        # Drop ALLELE_FLIPPED as it's an internal temporary column
        if "ALLELE_FLIPPED" in sumstats.columns:
            sumstats = sumstats.drop(columns=["ALLELE_FLIPPED"])
        
        if not is_dataframe:
            sumstats_obj.data = sumstats
            return sumstats_obj.data
        else:
            return sumstats

    before_missing = sumstats[rsid].isna().sum()

    # ---------------------------
    # Determine lookup TSV - only if we have variants to assign
    # ---------------------------
    if vcf_path is not None:
        if is_vcf_file(vcf_path):
            ref_mode = "vcf/bcf"
            path_to_use = vcf_path
        else:
            ref_mode = "tsv"
            path_to_use = vcf_path
    elif path is not None:
        if is_vcf_file(path):
            ref_mode = "vcf/bcf" 
            path_to_use = path
        else:
            ref_mode = "tsv"
            path_to_use = path
    else:
        ref_mode = "tsv"
        path_to_use = tsv_path

    log.write(" -Determining reference mode: {}...".format(ref_mode), verbose=verbose)

    if ref_mode == "tsv":
        # path_to_use is tsv
        if path_to_use is None or not os.path.exists(path_to_use):
            raise FileNotFoundError(f"Lookup TSV not found: {path_to_use}")
        log.write(" -Using TSV directly for lookup: {}...".format(path_to_use), verbose=verbose)
        lookup_tsv = path_to_use
        rm_tmp_lookup = False
    else:  
        # path_to_use is vcf/bcf
        if convert_to_bcf:
            if len(path_to_use) < 4 or path_to_use[-4:] != ".bcf":
                log.write(" -Converting VCF to BCF (strip_info={})...".format(strip_info), verbose=verbose)
                path_to_use = _convert_vcf_to_bcf(path_to_use, threads=threads, strip=strip_info, log=log, verbose=verbose)
            else:
                log.write(" -Already bcf", verbose=verbose)

        log.write(" -Extracting new lookup TSV from: {}...".format(path_to_use), verbose=verbose)

        # Optimize: only extract lookup for variants that need assignment
        # Create a subset of sumstats with only variants that need rsID assignment
        variants_needing_rsid = sumstats.loc[to_assign_mask & sumstats[rsid].isna(), [chrom, pos]].rename(
            columns={chrom: "CHR", pos: "POS"}
        )
        
        # Auto-detect reference format from VCF file
        mapper.detect_reference_format(path_to_use)
        
        lookup_tsv, rm_tmp_lookup = _extract_lookup_table_from_vcf_bcf(
            vcf_path=path_to_use,
            sumstats=variants_needing_rsid,
            mapper=mapper,
            assign_cols=["rsID"],  # ensures the lookup carries ID/rsID usable by assigner
            out_lookup=lookup_path,
            threads=threads,
            verbose=verbose,
            log=log,
        )

    sumstats = _assign_from_lookup(
        sumstats=sumstats,
        lookup_table=lookup_tsv,
        assign_cols=("rsID",),   # function will map ID → rsID automatically
        chrom=chrom,
        pos=pos,
        ea=ea,
        nea=nea,
        verbose=verbose,
        log=log,
        rm_tmp_lookup=rm_tmp_lookup
    )

    after_missing = sumstats[rsid].isna().sum()
    filled = before_missing - after_missing

    had_rsid_after = (~sumstats[rsid].isna()).sum()

    log.write(" -Filled {} rsIDs...".format(filled), verbose=verbose)
    log.write(" -rsID count: {} → {} / {}...".format(had_rsid_before, had_rsid_after, total_before), verbose=verbose)
    log.write(" -Finished assigning rsID from reference.", verbose=verbose)

    # Update harmonization status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _append_meta_record, _update_harmonize_step
            # Determine which path was used (path_to_use is already determined above)
            if path_to_use is not None:
                if ref_mode == "vcf/bcf":
                    sumstats_obj.meta["gwaslab"]["references"]["ref_rsid_vcf"] = _append_meta_record(
                        sumstats_obj.meta["gwaslab"]["references"]["ref_rsid_vcf"], path_to_use)
                else:
                    sumstats_obj.meta["gwaslab"]["references"]["ref_rsid_tsv"] = path_to_use
                assign_rsid_kwargs = {
                    'path': path, 'vcf_path': vcf_path, 'tsv_path': tsv_path, 'lookup_path': lookup_path,
                    'threads': threads, 'rsid': rsid, 'chrom': chrom, 'pos': pos, 'ea': ea, 'nea': nea,
                    'overwrite': overwrite, 'convert_to_bcf': convert_to_bcf, 'strip_info': strip_info
                }
                _update_harmonize_step(sumstats_obj, "assign_rsid", assign_rsid_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

@with_logging(
    start_to_msg="annotate GWAS summary statistics",
    finished_msg="annotating GWAS summary statistics",
    start_cols=["CHR","POS","EA","NEA"],
    start_function=".annotate_sumstats()"
)
def _annotate_sumstats(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    lookup_path: str | None = None,
    assign_cols=("rsID",),
    mapper: ChromosomeMapper | None = None,
    threads=6,
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    reuse_lookup=True,
    convert_to_bcf=False,
    strip_info=True,
    verbose=True,
    log=Log()
):
    """
    Annotate GWAS summary statistics by assigning fields (e.g., rsID, AF)
    from a lookup table extracted from a VCF/BCF.

    Two modes:
      (1) If tsv_path exists and reuse_lookup=True → skip extraction.
      (2) Otherwise extract from VCF → create tsv_path → annotate.

    Parameters
    ----------
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object. Must contain CHR, POS, EA, NEA.
    vcf_path : str or None
        Required if lookup table needs to be generated.
    tsv_path : str
        Lookup table file (tsv or tsv.gz).
    assign_cols : tuple[str]
        Columns to assign (e.g., ("rsID","AF")).
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper with automatic format detection.
    threads : int
        bcftools threads.
    chrom, pos, ea, nea : str
        Column names in sumstats.
    reuse_lookup : bool
        If True, reuse lookup if exists.

    Returns
    -------
    sumstats : pd.DataFrame
    tsv_path : str
    """
    import os

    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # ----------------------------------------------
    # Initialize ALLELE_FLIPPED column early (before early exit check)
    # ----------------------------------------------
    # Initialize ALLELE_FLIPPED column to track which variants had swapped alleles
    # This must be done before the early exit check to ensure it's always present
    if "ALLELE_FLIPPED" not in sumstats.columns:
        sumstats["ALLELE_FLIPPED"] = False
        log.write(f" -Initialized ALLELE_FLIPPED column (early initialization)", verbose=verbose)
    
    # ----------------------------------------------
    # Early check: determine which columns need annotation
    # ----------------------------------------------
    # Initialize missing columns and check which need annotation
    cols_to_assign = []
    for col in assign_cols:
        if col not in sumstats.columns:
            sumstats[col] = pd.NA
            cols_to_assign.append(col)
        elif sumstats[col].isna().any():
            # Column exists but has missing values
            cols_to_assign.append(col)
    
    # Early exit if no columns need annotation
    # However, we still need to check for ALLELE_FLIPPED if it's all False
    # (meaning it was just initialized and hasn't been set based on matching yet)
    if not cols_to_assign:
        # Check if ALLELE_FLIPPED needs to be set based on matching
        # If ALLELE_FLIPPED exists but is all False, we should still do matching
        # to determine which variants are flipped
        need_flip_check = (
            "ALLELE_FLIPPED" in sumstats.columns and 
            sumstats["ALLELE_FLIPPED"].fillna(False).sum() == 0 and
            sumstats["ALLELE_FLIPPED"].notna().all()
        )
        if not need_flip_check:
            log.write(" -All annotation columns already present, skipping annotation...", verbose=verbose)
            if not is_dataframe:
                sumstats_obj.data = sumstats
                return sumstats_obj.data
            else:
                return sumstats
        else:
            log.write(" -All annotation columns already present, but checking for allele flips...", verbose=verbose)
            # Continue to matching logic to set ALLELE_FLIPPED
            # We'll skip annotation updates but still do matching
            # Set a flag to indicate we only need to update ALLELE_FLIPPED, not annotation columns
            only_update_flip = True

    # ----------------------------------------------
    # Step 1 — reuse or extract lookup table
    # ----------------------------------------------
    if vcf_path is not None:
        if is_vcf_file(vcf_path):
            ref_mode = "vcf/bcf"
            path_to_use = vcf_path
        else:
            ref_mode = "tsv"
            path_to_use = vcf_path
    elif path is not None:
        if is_vcf_file(path):
            ref_mode = "vcf/bcf" 
            path_to_use = path
        else:
            ref_mode = "tsv"
            path_to_use = path
    else:
        ref_mode = "tsv"
        path_to_use = tsv_path

    log.write(" -Annotating {} from {}".format(",".join(assign_cols), path_to_use if path_to_use else "lookup"))
    log.write(" -Determining reference mode: {}...".format(ref_mode), verbose=verbose)

    if ref_mode == "tsv":
        # path_to_use is tsv
        if path_to_use is None or not os.path.exists(path_to_use):
            raise FileNotFoundError(f"Lookup TSV not found: {path_to_use}")
        log.write(" -Using TSV directly for lookup: {}...".format(path_to_use), verbose=verbose)
        lookup_tsv = path_to_use
        rm_tmp_lookup = False
    else:  
        # path_to_use is vcf/bcf
        if convert_to_bcf:
            # Fix: use proper string check instead of indexing
            if len(path_to_use) < 4 or path_to_use[-4:] != ".bcf":
                log.write(" -Converting VCF to BCF (strip_info={})...".format(strip_info), verbose=verbose)
                path_to_use = _convert_vcf_to_bcf(path_to_use, threads=threads, strip=strip_info, log=log, verbose=verbose)
            else:
                log.write(" -Already bcf", verbose=verbose)

        log.write(" -Creating lookup table from VCF: {}...".format(path_to_use), verbose=verbose)
        
        # Optimize: only extract lookup for variants that need annotation
        # Determine which variants need annotation (those with missing values in any assign_col)
        # Since we already initialized missing columns above, all assign_cols should exist
        assign_cols_list = list(assign_cols)
        missing_mask = sumstats[assign_cols_list].isna().any(axis=1)
        if missing_mask.any():
            # Only extract lookup for variants that need annotation
            variants_needing_annotation = sumstats.loc[missing_mask, [chrom, pos]].rename(
                columns={chrom: "CHR", pos: "POS"}
            )
        else:
            # All variants already have all annotations - this shouldn't happen after our early check
            # but handle gracefully by using all variants
            variants_needing_annotation = sumstats[[chrom, pos]].rename(
                columns={chrom: "CHR", pos: "POS"}
            )
        
        # Auto-detect reference format from VCF file
        mapper.detect_reference_format(path_to_use)
        
        lookup_tsv, rm_tmp_lookup = _extract_lookup_table_from_vcf_bcf(
            vcf_path   = path_to_use,
            sumstats   = variants_needing_annotation,
            mapper     = mapper,
            assign_cols= assign_cols,
            out_lookup = lookup_path,
            threads    = threads,
            log        = log,
            verbose    = verbose
        )

    # ----------------------------------------------
    # Step 2 — assign annotation fields
    # ----------------------------------------------
    sumstats = _assign_from_lookup(
        sumstats     = sumstats,
        lookup_table = lookup_tsv,
        assign_cols  = assign_cols,
        chrom        = chrom,
        pos          = pos,
        ea           = ea,
        nea          = nea,
        log          = log,
        verbose      = verbose,
        rm_tmp_lookup=rm_tmp_lookup
    )

    # Update only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        return sumstats_obj.data
    else:
        return sumstats

@with_logging(
    start_to_msg="extract lookup table from vcf/bcf",
    finished_msg="extracting lookup table from vcf/bcf",
    start_cols=["CHR","POS"],
    start_function="._extract_lookup_table_from_vcf_bcf()"
)
def _extract_lookup_table_from_vcf_bcf_old(
    vcf_path,
    sumstats,
    chr_dict = None,
    assign_cols   = None,
    out_lookup=None,
    threads=6,
    verbose=True,
    rm_out_lookup=False,
    log=Log()
    ):
    # Early exit if sumstats is empty
    if sumstats.empty or len(sumstats) == 0:
        log.write(" -No variants to extract, creating empty lookup table...", verbose=verbose)
        if out_lookup is None:
            out_lookup_tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".gwaslab.lookup.txt")
            out_lookup = out_lookup_tmp.name
            out_lookup_tmp.close()
        
        # Create empty lookup with proper headers
        info_tags = [col for col in (assign_cols or []) if col not in ("ID","rsID")]
        id_col = "rsID" if (assign_cols and "rsID" in assign_cols) else "ID"
        header_cols = ["CHR", "POS", "REF", "ALT", id_col] + info_tags
        header_line = "\t".join(header_cols) + "\n"
        
        if out_lookup.endswith(".gz"):
            import gzip
            with gzip.open(out_lookup, "wt") as f:
                f.write(header_line)
        else:
            with open(out_lookup, "w") as f:
                f.write(header_line)
        
        log.write(" -Empty lookup table created: {}...".format(out_lookup), verbose=verbose)
        return out_lookup, rm_out_lookup

    # Note: This is a deprecated function. Use _extract_lookup_table_from_vcf_bcf instead.
    # Keeping chr_dict for backward compatibility but it's deprecated.
    if chr_dict is None:
        # Create mapper and detect format
        mapper = ChromosomeMapper(log=log, verbose=verbose)
        if not sumstats.empty and "CHR" in sumstats.columns:
            mapper.detect_sumstats_format(sumstats["CHR"])
        mapper.detect_reference_format(vcf_path)
        # Create chr_dict from mapper for backward compatibility
        # This is a temporary bridge - should use mapper directly
        from gwaslab.bd.bd_common_data import get_number_to_chr
        chr_dict = get_number_to_chr(prefix="chr" if mapper._reference_prefix == "chr" else "")
        log.write(" -Auto-determined chr_dict: {}...".format(chr_dict), verbose=verbose)

    if assign_cols is None:
        assign_cols = []
    if not shutil.which("bcftools"):
        raise RuntimeError("bcftools not found in PATH")

    if out_lookup is None:
        out_lookup_tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".gwaslab.lookup.txt")
        out_lookup = out_lookup_tmp.name
        out_lookup_tmp.close()

    # Optimize: Create targets more efficiently
    # Pre-compute inverse dict if needed
    inv_chr_dict = {v: k for k, v in chr_dict.items()} if chr_dict else None
    
    # Optimize: Avoid unnecessary copy - work with selected columns directly
    targets_df = sumstats[["CHR", "POS"]].copy()
    
    # Apply chr_dict mapping if needed
    if chr_dict is not None:
        log.write(" -Converting chromosome notation using chr_dict to reference notation...", verbose=verbose)
        targets_df["CHR"] = targets_df["CHR"].map(chr_dict)
    
    # Create targets file - optimize: drop NA and duplicates before sorting
    tmp_targets = tempfile.NamedTemporaryFile(delete=False, suffix=".gwaslab.targets.sorted.txt.gz")
    _targets_path = tmp_targets.name
    tmp_targets.close()
    
    # Optimize: Chain operations more efficiently
    targets_df = targets_df.dropna(subset=["CHR"]).drop_duplicates()
    if targets_df.empty:
        # No valid targets after filtering - create empty lookup table
        log.write(" -No valid targets after filtering, creating empty lookup table...", verbose=verbose)
        info_tags = [col for col in assign_cols if col not in ("ID","rsID")]
        id_col = "rsID" if ("rsID" in assign_cols) else "ID"
        header_cols = ["CHR", "POS", "REF", "ALT", id_col] + info_tags
        header_line = "\t".join(header_cols) + "\n"
        
        if out_lookup.endswith(".gz"):
            import gzip
            with gzip.open(out_lookup, "wt") as f:
                f.write(header_line)
        else:
            with open(out_lookup, "w") as f:
                f.write(header_line)
        
        log.write(" -Empty lookup table created: {}...".format(out_lookup), verbose=verbose)
        return out_lookup, rm_out_lookup
    
    # Sort and write targets
    targets_df = targets_df.sort_values(["CHR", "POS"])
    targets_df.to_csv(_targets_path, sep="\t", header=False, index=False, compression="gzip")

    log.write(" -Created target list: {}...".format(_targets_path), verbose=verbose)

    # ---- extract from VCF/BCF ----
    _tmp_filtered_bcf = tempfile.NamedTemporaryFile(delete=False, suffix=".gwaslab.filtered.bcf").name
    cmd_filter = [
        "bcftools", "view",
        "-T", _targets_path,
        "-Ob", "-o", _tmp_filtered_bcf,
        "--threads", str(threads),
        vcf_path
    ]
    log.write(" -Extracting target sites...", verbose=verbose)
    log.write(" -Calling: {}".format(" ".join(cmd_filter)), verbose=verbose)
    
    subprocess.check_call(cmd_filter)
    subprocess.check_call(["bcftools", "index", "-f", _tmp_filtered_bcf])

    # Optimize: Build fmt string more efficiently
    info_tags = [col for col in assign_cols if col not in ("ID","rsID")]
    # Build fmt string in one go instead of loop
    info_fmt = "".join(f"\t%INFO/{tag}" for tag in info_tags)
    fmt = f"%CHROM\t%POS\t%REF\t%ALT\t%ID{info_fmt}\n"

    id_col = "rsID" if ("rsID" in assign_cols) else "ID"
    header_cols = ["CHR", "POS", "REF", "ALT", id_col] + info_tags
    header_line = "\t".join(header_cols) + "\n"
    cmd_query = ["bcftools", "query", "-f", fmt, _tmp_filtered_bcf]

    log.write(" -Writing lookup to {}...".format(out_lookup), verbose=verbose)
    
    # Optimize: Process data in streaming fashion to avoid double I/O
    # First, write header and query output to a temporary location for processing
    _tmp_query_output = tempfile.NamedTemporaryFile(delete=False, suffix=".gwaslab.query.txt")
    _tmp_query_path = _tmp_query_output.name
    _tmp_query_output.close()
    
    # Write query output to temp file
    with open(_tmp_query_path, "w") as tmp_f:
        log.write(" -Calling: {}".format(" ".join(cmd_query)), verbose=verbose)
        p = subprocess.Popen(cmd_query, stdout=tmp_f, text=True)
        p.wait()
        if p.returncode != 0:
            raise RuntimeError(f"bcftools query failed with return code {p.returncode}")

    # Read and process in one go
    # Handle empty results gracefully
    try:
        df = pd.read_csv(_tmp_query_path, sep="\t", header=None, 
                         names=["CHR", "POS", "REF", "ALT", "ID"] + info_tags)
        if df.empty:
            # Create empty dataframe with correct columns
            df = pd.DataFrame(columns=["CHR", "POS", "REF", "ALT", "ID"] + info_tags)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        # Handle completely empty file
        df = pd.DataFrame(columns=["CHR", "POS", "REF", "ALT", "ID"] + info_tags)

    # Convert ID → rsID if needed
    if id_col == "rsID" and "ID" in df.columns:
        df = df.rename(columns={"ID": "rsID"})

    # Convert CHR back using inverse dictionary
    if inv_chr_dict is not None:
        log.write(" -Converting CHR back to original sumstats notation...", verbose=verbose)
        df["CHR"] = df["CHR"].astype(str).map(inv_chr_dict).astype("category")

    # Write final output
    df.to_csv(out_lookup, sep="\t", index=False, compression="infer")
    
    # Cleanup temp query file
    try:
        os.remove(_tmp_query_path)
    except:
        pass
    
    # Optimize cleanup: collect all files to clean up
    to_clean_up = [
        _targets_path,
        _tmp_filtered_bcf,
        _tmp_filtered_bcf + ".csi", 
        _tmp_filtered_bcf + ".tbi"
    ]
    
    # Cleanup with better error handling
    for f in to_clean_up:
        if isinstance(f, str) and os.path.exists(f):
            try:
                os.remove(f)
                log.write(" -Cleaning up: {}...".format(f), verbose=verbose)
            except OSError as e:
                log.write(" -Warning: Could not remove {}: {}...".format(f, e), verbose=verbose)

    log.write(" -Lookup table created: {}...".format(out_lookup), verbose=verbose)
    return out_lookup, rm_out_lookup
from io import StringIO
import pandas as pd
import subprocess, tempfile, os

def _worker_bcf_lookup(args):
    (chr_val, df_chr, vcf_path, fmt, split_by_chr,
     mapper, info_tags, id_col, log, verbose) = args

    import tempfile, subprocess, os
    import pandas as pd
    from io import StringIO
    # Convert reference chromosome back to sumstats format for logging
    orig_chr = mapper.reference_to_sumstats(chr_val, reference_file=vcf_path)

    # Create per-chr target list
    tmp_targets = tempfile.NamedTemporaryFile(
        delete=False, suffix=f".{orig_chr}.targets.txt"
    )
    targets_path = tmp_targets.name
    tmp_targets.close()

    # Don't sort targets - preserve original order to match VCF file order
    # Sorting targets causes bcftools to output in sorted order, not VCF file order
    # This is critical for matching the old method which processes variants in VCF order
    # Sort targets by CHR:POS for efficient bcftools processing
    # Note: This may cause bcftools to output in position-sorted order rather than
    # original VCF file order, which can lead to different rsID selection when
    # multiple rsIDs exist at the same position compared to the old method
    # (which uses vcf_reader.fetch() and preserves VCF file order)
    df_chr[["CHR","POS"]].sort_values(["CHR","POS"]) \
        .to_csv(targets_path, sep="\t", header=False, index=False)

    # FAST bcftools pipeline: uncompressed streaming (-Ou)
    if split_by_chr:
        cmd = (
            f"bcftools view -r {chr_val} -T {targets_path} -Ou {vcf_path} | "
            f'bcftools query -f "{fmt}"'
        )
    else:
        cmd = (
            f"bcftools view -T {targets_path} -Ou {vcf_path} | "
            f"bcftools query -f '{fmt}'"
        )

    
    # Run pipeline and capture text output
    data = subprocess.check_output(cmd, shell=True, text=True)
    log.write(f" {orig_chr}", verbose=verbose, end="", show_time=False)

    # Convert to DataFrame
    if data.strip() == "":
        df_out = pd.DataFrame(columns=["CHR","POS","REF","ALT",id_col] + info_tags)
    else:
        df_out = pd.read_csv(StringIO(data), sep="\t", header=None)
        df_out.columns = ["CHR","POS","REF","ALT",id_col] + info_tags

    # cleanup only target file
    try:
        os.remove(targets_path)
    except:
        pass
    
    return df_out
from multiprocessing import Pool

@with_logging(
    start_to_msg="extract lookup table from vcf/bcf",
    finished_msg="extracting lookup table from vcf/bcf",
    start_cols=["CHR","POS"],
    start_function="._extract_lookup_table_from_vcf_bcf()"
)
def _extract_lookup_table_from_vcf_bcf(
    vcf_path,
    sumstats,
    mapper: ChromosomeMapper | None = None,
    assign_cols=None,
    out_lookup=None,
    verbose=True,
    rm_out_lookup=False,
    split_by_chr=True,
    threads=6,     # <-- threads for Pool
    log=Log()
):
    import os, shutil, tempfile, pandas as pd

    def is_indexed(p):
        return os.path.exists(p + ".tbi") or os.path.exists(p + ".csi")

    if split_by_chr is None:
        split_by_chr = not is_indexed(vcf_path)

    # Get or create mapper
    if mapper is None:
        # Create default mapper
        mapper = ChromosomeMapper(log=log, verbose=verbose)
        # Auto-detect sumstats format
        if not sumstats.empty and "CHR" in sumstats.columns:
            mapper.detect_sumstats_format(sumstats["CHR"])
    
    # Auto-detect reference format from VCF file
    mapper.detect_reference_format(vcf_path)

    if assign_cols is None:
        assign_cols = []

    if out_lookup is None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".lookup.txt.gz")
        out_lookup = tmp.name
        tmp.close()

    # Convert sumstats chromosomes to reference format
    targets = sumstats[["CHR","POS"]].copy()
    log.write(" -Converting chromosome notation to reference notation...", verbose=verbose)
    targets["CHR"] = mapper.sumstats_to_reference_series(targets["CHR"], reference_file=vcf_path)
    targets = targets.dropna().drop_duplicates()
    # Group by chromosome and sort by chromosome to ensure consistent processing order
    # Note: Sorting chromosomes may affect the order in which bcftools processes records,
    # potentially leading to different rsID selection when multiple rsIDs exist at the same
    # position compared to the old method (which processes variants in sumstats order)
    chr_groups = {chr_val: df_chr for chr_val, df_chr in sorted(targets.groupby("CHR"))}

    info_tags = [c for c in assign_cols if c not in ("ID","rsID")]
    id_col = "rsID" if ("rsID" in assign_cols) else "ID"
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID" + "".join([f"\t%INFO/{t}" for t in info_tags]) + "\n"

    # Prepare args for pool
    tasks = []
    for chr_val, df_chr in chr_groups.items():
        tasks.append((chr_val, df_chr, vcf_path, fmt, split_by_chr,
                      mapper, info_tags, id_col, log, verbose))

    log.write(f" -Running multiprocessing: {threads} workers, {len(tasks)} chromosomes",
              verbose=verbose)
    
    if split_by_chr:
        log.write(" -Calling: "
                  "bcftools view -r <CHR> -T <TARGETS> -Ou <VCF>| "
                  "bcftools query -f '<FMT>'", verbose=verbose)
    else:
        # TEMPLATE log
        log.write(" -Calling: bcftools view -T <TARGETS> -Ou <VCF> | "
                  "bcftools query -f '<FMT>'", verbose=verbose)

    with Pool(threads) as pool:
        log.write(f" -Finished:",end="", verbose=verbose)
        dfs = pool.map(_worker_bcf_lookup, tasks)
        log.write(f"", verbose=verbose, show_time=False)

    # Merge results - filter out empty DataFrames to avoid FutureWarning
    dfs_filtered = [d for d in dfs if not d.empty]
    if dfs_filtered:
        df = pd.concat(dfs_filtered, ignore_index=True)
    else:
        # If all DataFrames are empty, create an empty DataFrame with the expected columns
        df = pd.DataFrame(columns=["CHR","POS","REF","ALT",id_col] + info_tags)
    
    if not df.empty:
        # Convert reference chromosomes back to sumstats format
        log.write(" -Converting chromosome notation back to sumstats notation...", verbose=verbose)
        df["CHR"] = mapper.reference_to_sumstats_series(df["CHR"], reference_file=vcf_path)

    if id_col == "rsID" and "ID" in df.columns:
        df.rename(columns={"ID":"rsID"}, inplace=True)

    df.sort_values(["CHR","POS"], inplace=True)
    df.to_csv(out_lookup, sep="\t", index=False, compression="infer")

    log.write(f" -Lookup table created: {out_lookup}", verbose=verbose)
    return out_lookup, rm_out_lookup

import pandas as pd
import numpy as np

@with_logging(
    start_to_msg="assign from lookup table",
    finished_msg="assigning from lookup table",
    start_cols=["CHR","POS"],
    start_function="._assign_from_lookup()",
    must_kwargs=["lookup_table"]
)
def _assign_from_lookup(
    sumstats,
    lookup_table,
    assign_cols=("rsID",),
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    verbose=True,
    log=Log(),
    rm_tmp_lookup=False
):
    """
    Assign annotation values from a lookup table to sumstats by matching CHR:POS:EA:NEA.
    
    This function performs allele-aware matching, handling both forward (EA:NEA) and 
    reverse (NEA:EA) allele orders. It processes the lookup table in chunks for memory
    efficiency and tracks which variants were updated and which had allele flips.
    """
    import pandas as pd
    import numpy as np
    chunksize = 5_000_000  # Process lookup table in chunks of 5M rows to manage memory

    # ============================================================================
    # Step 1: Validate lookup table file exists
    # ============================================================================
    import os
    if isinstance(lookup_table, str) and not os.path.exists(lookup_table):
        log.warning(f"Lookup table not found: {lookup_table}, skipping assignment...", verbose=verbose)
        return sumstats

    # ============================================================================
    # Step 2: Read lookup table header to understand structure
    # ============================================================================
    # Read only the header (0 rows) to get column names without loading data
    lookup_header = pd.read_csv(lookup_table, sep="\t", nrows=0).columns.tolist()

    # ============================================================================
    # Step 3: Initialize ALLELE_FLIPPED column if not present
    # ============================================================================
    # ALLELE_FLIPPED is used to track which variants had swapped alleles during matching
    # It should be initialized before we start processing
    if "ALLELE_FLIPPED" not in sumstats.columns:
        sumstats["ALLELE_FLIPPED"] = False
        log.write(" -Initialized ALLELE_FLIPPED column", verbose=verbose)
    
    # ============================================================================
    # Step 4: Check if lookup table has any data
    # ============================================================================
    # Try to read first row to verify file contains data (not just headers)
    try:
        first_row = pd.read_csv(lookup_table, sep="\t", nrows=1)
        if first_row.empty:
            log.write(" -Lookup table is empty, skipping assignment...", verbose=verbose)
            return sumstats
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        # Handle case where file exists but is empty or malformed
        log.write(" -Lookup table is empty, skipping assignment...", verbose=verbose)
        return sumstats

    # ============================================================================
    # Step 4: Detect allele column naming convention in lookup table
    # ============================================================================
    # Lookup tables can use either REF/ALT (VCF convention) or EA/NEA (sumstats convention)
    # We need to identify which columns represent effect/non-effect alleles
    if ("REF" in lookup_header) and ("ALT" in lookup_header):
        # VCF convention: REF = non-effect allele, ALT = effect allele
        lookup_nea_col = "REF"
        lookup_ea_col  = "ALT"
        mode = "REF_ALT"
    elif ("NEA" in lookup_header) and ("EA" in lookup_header):
        # Sumstats convention: NEA = non-effect allele, EA = effect allele
        lookup_nea_col = "NEA"
        lookup_ea_col  = "EA"
        mode = "EA_NEA"
    else:
        raise ValueError(f"Lookup must contain REF/ALT or NEA/EA. Found: {lookup_header}")

    log.write(f" -Detected allele mode: {mode} using {lookup_ea_col}(EA/ALT) / {lookup_nea_col}(NEA/REF)...",
              verbose=verbose)

    # ============================================================================
    # Step 5: Normalize ID column names (rsID vs ID)
    # ============================================================================
    # Save original requested columns for reporting dropped fields
    original_assign_cols = assign_cols

    # Determine which ID column exists in lookup (rsID preferred over ID)
    # Case 1: lookup has both → prefer rsID
    if "rsID" in lookup_header and "ID" in lookup_header:
        id_col = "rsID"
    # Case 2: lookup has only one
    elif "rsID" in lookup_header:
        id_col = "rsID"
    elif "ID" in lookup_header:
        id_col = "ID"
    else:
        id_col = None  # no rsID or ID in lookup

    # ============================================================================
    # Step 6: Map requested columns to available columns in lookup
    # ============================================================================
    # Convert requested assign_cols to columns that actually exist in lookup
    # - If user requests "rsID" or "ID", use whichever exists (id_col)
    # - For other columns, only include if they exist in lookup
    lookup_header_set = set(lookup_header)  # Use set for O(1) lookup performance
    normalized_assign_cols = []
    for col in assign_cols:
        if col in ("rsID", "ID"):
            # User requested ID/rsID - use whichever column exists in lookup
            if id_col is not None:
                normalized_assign_cols.append(id_col)
        elif col in lookup_header_set:
            # Column exists in lookup, include it
            normalized_assign_cols.append(col)

    # Remove duplicates but preserve order (dict.fromkeys maintains insertion order)
    assign_cols = tuple(dict.fromkeys(normalized_assign_cols))

    # ============================================================================
    # Step 7: Report any requested columns that are not available
    # ============================================================================
    dropped = set(original_assign_cols) - set(assign_cols)
    if dropped:
        log.warning("Annotation columns not available in lookup, skipped: {}...".format(dropped), verbose=verbose)

    # ============================================================================
    # Step 8: Early exit if no columns to assign
    # ============================================================================
    if not assign_cols:
        log.write(" -No columns to assign after normalization, skipping...", verbose=verbose)
        return sumstats

    # ============================================================================
    # Step 9: Initialize missing annotation columns in sumstats
    # ============================================================================
    # Create columns in sumstats if they don't exist, fill with NA
    for col in assign_cols:
        if col not in sumstats.columns:
            sumstats[col] = pd.NA
    # ALLELE_FLIPPED column is already initialized earlier (before early exit check)

    # ============================================================================
    # Step 10: Prepare column selection and data types for chunked reading
    # ============================================================================
    # Define which columns to read from lookup table
    usecols = [chrom, pos, lookup_ea_col, lookup_nea_col] + list(assign_cols)
    # Define data types for efficient memory usage
    dtype = {chrom:"int64", pos:"int64",
             lookup_ea_col:"category", lookup_nea_col:"category"}

    # ============================================================================
    # Step 11: Track initial state of missing annotations
    # ============================================================================
    # Record which rows have missing annotations BEFORE assignment
    # This allows us to track which rows were newly filled (not just overwritten)
    assign_cols_list = list(assign_cols)
    initial_missing = sumstats[assign_cols_list].isna()
    
    # ============================================================================
    # Step 12: Early exit if all annotations already filled
    # ============================================================================
    if not initial_missing.any().any():
        log.write(" -All annotation columns already filled, skipping lookup...", verbose=verbose)
        return sumstats
    
    # ============================================================================
    # Step 13: Initialize tracking sets for statistics
    # ============================================================================
    # Track which rows were updated and which had allele flips for final reporting
    updated_rows = set()  # Set of row indices that were newly annotated
    flipped_rows = set()  # Set of row indices that had allele flips

    # ============================================================================
    # Step 14: Chunk-wise streaming processing
    # ============================================================================
    # Process lookup table in chunks to manage memory for large files
    # Pre-compute unique chromosomes in sumstats for faster filtering
    sumstats_chrs = set(sumstats[chrom].unique())
    
    # Iterate over chunks of the lookup table
    for chunk in pd.read_csv(
        lookup_table, sep="\t", usecols=usecols,
        dtype=dtype, chunksize=chunksize
    ):
        # ========================================================================
        # Step 14a: Skip empty chunks
        # ========================================================================
        if chunk.empty:
            continue
            
        # ========================================================================
        # Step 14b: Expand multi-allelic variants
        # ========================================================================
        # Some variants have multiple alternate alleles (e.g., REF=A, ALT=T,G)
        # Expand these into separate rows for each ALT allele
        try:
            chunk = _expand_multiallelic_fast(chunk, lookup_ea_col, lookup_nea_col)
        except Exception:
            log.warning(" -Multiallelic expansion failed: using original rows...",
                        verbose=verbose)

        lookup_rows = len(chunk)
        if lookup_rows == 0:
            continue
            
        log.write(f" -Loaded {lookup_rows:,} lookup rows...", verbose=verbose)

        # ========================================================================
        # Step 14c: Filter by chromosome for efficiency
        # ========================================================================
        # Only process chromosomes that exist in both sumstats and this chunk
        # This avoids unnecessary processing of non-matching chromosomes
        chunk_chrs = set(chunk[chrom].unique())
        matching_chrs = chunk_chrs & sumstats_chrs  # Set intersection
        if not matching_chrs:
            log.write(" -No matching CHR in this chunk, skipping...", verbose=verbose)
            continue
            
        # Filter sumstats to only rows with matching chromosomes
        ss_sub = sumstats[sumstats[chrom].isin(matching_chrs)]
        if ss_sub.empty:
            log.write(" -No matching CHR in this chunk, skipping...", verbose=verbose)
            continue
        
        # Create a copy because we'll modify allele columns (convert to categorical)
        ss_sub = ss_sub.copy()

        # ========================================================================
        # Step 14d: Create unified allele symbol space
        # ========================================================================
        # Collect all unique allele symbols from both sumstats and lookup chunk
        # This ensures categorical dtype includes all possible alleles for matching
        allele_space = set()
        allele_space.update(ss_sub[ea].dropna().unique())
        allele_space.update(ss_sub[nea].dropna().unique())
        allele_space.update(chunk[lookup_ea_col].dropna().unique())
        allele_space.update(chunk[lookup_nea_col].dropna().unique())
        
        if not allele_space:
            log.write(" -No valid alleles in chunk, skipping...", verbose=verbose)
            continue
            
        # Convert to categorical dtype for efficient matching
        # Categorical dtype allows fast equality comparisons
        alleles = pd.CategoricalDtype(categories=list(allele_space), ordered=False)

        # Convert all allele columns to the same categorical dtype
        # This ensures allele symbols match exactly (e.g., "A" == "A" not "A" == "A ")
        ss_sub[ea]  = ss_sub[ea].astype(alleles)
        ss_sub[nea] = ss_sub[nea].astype(alleles)
        chunk[lookup_ea_col]  = chunk[lookup_ea_col].astype(alleles)
        chunk[lookup_nea_col] = chunk[lookup_nea_col].astype(alleles)

        # ========================================================================
        # Step 14e: Create lookup index for fast matching
        # ========================================================================
        # Build a MultiIndex lookup table keyed by CHR:POS:NEA:EA
        # Remove duplicates first (keep first occurrence if same variant appears multiple times)
        # Note: The order of deduplication may differ from the old method (which uses
        # vcf_reader.fetch() and preserves VCF file order), potentially leading to minor
        # differences (0.12% in test cases) when multiple rsIDs exist at the same position.
        lookup = (
            chunk.drop_duplicates([chrom, pos, lookup_nea_col, lookup_ea_col])
                 .set_index([chrom, pos, lookup_nea_col, lookup_ea_col])[list(assign_cols)]
        )
        
        if lookup.empty:
            log.write(" -Lookup index is empty after deduplication, skipping...", verbose=verbose)
            continue

        # ========================================================================
        # Step 14f: Create MultiIndex keys for sumstats variants
        # ========================================================================
        # Create two sets of keys to match against lookup:
        # - key_fwd: CHR:POS:NEA:EA (forward orientation, matches if alleles are in same order)
        # - key_rev: CHR:POS:EA:NEA (reverse orientation, matches if alleles are swapped)
        # This allows us to handle cases where sumstats and lookup have opposite allele orders
        key_fwd = pd.MultiIndex.from_arrays(
            [ss_sub[chrom], ss_sub[pos], ss_sub[nea], ss_sub[ea]],
            names=[chrom, pos, nea, ea]
        )
        key_rev = pd.MultiIndex.from_arrays(
            [ss_sub[chrom], ss_sub[pos], ss_sub[ea], ss_sub[nea]],
            names=[chrom, pos, ea, nea]
        )

        # ========================================================================
        # Step 14g: Lookup annotation values for both orientations
        # ========================================================================
        # Try to find matches in lookup table for both forward and reverse allele orders
        # vals_fwd: annotation values when alleles match in forward order
        # vals_rev: annotation values when alleles match in reverse order (flipped)
        vals_fwd = lookup.reindex(key_fwd).to_numpy()
        vals_rev = lookup.reindex(key_rev).to_numpy()

        # ========================================================================
        # Step 14h: Determine which variants had allele flips
        # ========================================================================
        # A variant is "flipped" if reverse match succeeded.
        # This matches the old method's logic which checks:
        #   1. if record.ref == NEA and EA in record.alts → return ALT_AF (forward)
        #   2. elif record.ref == EA and NEA in record.alts → return 1-ALT_AF (reverse/flipped)
        # The old method returns the FIRST match, so if reverse matches, it uses 1-ALT_AF.
        # To match this behavior, we mark as flipped if reverse match succeeded,
        # regardless of whether forward also succeeded (we'll prefer forward if both succeed,
        # but the ALLELE_FLIPPED flag indicates that reverse was a valid match).
        # However, to match the old method exactly, we should mark as flipped ONLY if
        # forward failed and reverse succeeded (old method's elif condition).
        # But actually, if both succeed, the old method would use forward (first condition),
        # so we should NOT mark as flipped. Only mark as flipped if forward failed.
        vals_fwd_na = pd.isna(vals_fwd)
        vals_rev_notna = pd.notna(vals_rev)
        # Mark as flipped if forward failed AND reverse succeeded
        # This matches the old method's elif condition (only checked if first condition failed)
        flipped = vals_fwd_na.all(axis=1) & vals_rev_notna.any(axis=1)
        
        # ========================================================================
        # Step 14i: Combine forward and reverse matches
        # ========================================================================
        # Use forward values where available (element-wise), otherwise use reverse values
        # This handles cases where some columns matched forward and others matched reverse
        assigned = np.where(pd.notna(vals_fwd), vals_fwd, vals_rev)

        # ========================================================================
        # Step 14j: Identify rows that need annotation updates
        # ========================================================================
        # Only update rows where annotations are currently missing
        # This prevents overwriting existing annotations
        missing_mask = sumstats.loc[ss_sub.index, assign_cols].isna().any(axis=1)
        rows_now = ss_sub.index[missing_mask]
        
        # Also identify rows that have matches in lookup (for ALLELE_FLIPPED update)
        # This includes rows with matches even if annotations are already present
        has_match = pd.notna(vals_fwd).any(axis=1) | pd.notna(vals_rev).any(axis=1)
        rows_with_match = ss_sub.index[has_match]
        
        # Debug: Count flipped variants in this chunk
        flipped_count = flipped.sum()
        if flipped_count > 0:
            log.write(f" -Found {flipped_count} flipped variants in this chunk", verbose=verbose)
        
        if len(rows_now) == 0:
            log.write(" -No missing annotations in this chunk, skipping annotation updates...", verbose=verbose)
            # Even if no annotations to update, we should still update ALLELE_FLIPPED
            # for rows that have matches in the lookup table
            if len(rows_with_match) > 0:
                sumstats.loc[rows_with_match, "ALLELE_FLIPPED"] |= flipped[has_match]
                log.write(f" -Updated ALLELE_FLIPPED for {len(rows_with_match)} variants (no missing annotations)", verbose=verbose)
            continue

        # ========================================================================
        # Step 14k: Update sumstats with annotation values
        # ========================================================================
        # Assign annotation values to rows that were missing them
        sumstats.loc[rows_now, assign_cols] = assigned[missing_mask]
        # Update ALLELE_FLIPPED flag for variants that had swapped alleles
        # Update for all rows with matches, not just rows with missing annotations
        if len(rows_with_match) > 0:
            sumstats.loc[rows_with_match, "ALLELE_FLIPPED"] |= flipped[has_match]
            flipped_updated = flipped[has_match].sum()
            if flipped_updated > 0:
                log.write(f" -Updated ALLELE_FLIPPED=True for {flipped_updated} variants in this chunk", verbose=verbose)

        # ========================================================================
        # Step 14l: Track which rows were newly filled (not just overwritten)
        # ========================================================================
        # Determine which rows were updated for the FIRST time (were missing before)
        # This is different from rows_now because rows_now includes all missing rows,
        # but we only want to count rows that were actually filled (not NA after assignment)
        newly_filled = (
            initial_missing.loc[rows_now, :] &  # Was missing before
            sumstats.loc[rows_now, assign_cols].notna()  # Is filled now
        )
        new_rows = newly_filled.index[newly_filled.any(axis=1)]
        updated_rows.update(new_rows.tolist())

        # ========================================================================
        # Step 14m: Track which newly filled rows had allele flips
        # ========================================================================
        # Map flipped flags to the corresponding sumstats indices for this chunk
        # Only count flips for rows that were newly filled (not pre-existing annotations)
        flipped_now = pd.Series(flipped[missing_mask], index=rows_now)
        flipped_new = flipped_now.loc[new_rows]
        flipped_rows.update(flipped_new.index[flipped_new].tolist())

        # ========================================================================
        # Step 14n: Update baseline missing map for next iteration
        # ========================================================================
        # Update the initial_missing tracking to reflect current state
        # This ensures we don't double-count rows in subsequent chunks
        initial_missing.loc[rows_now, :] = sumstats.loc[rows_now, assign_cols].isna()

        # ========================================================================
        # Step 14o: Log progress for this chunk
        # ========================================================================
        log.write(
            f" -Newly annotated sumstats rows: {len(new_rows):,} "
            f"(chunk lookup rows: {lookup_rows:,}) "
            f"| New flips: {len(flipped_new):,}",
            verbose=verbose
        )

    # ============================================================================
    # Step 15: Final statistics and reporting
    # ============================================================================
    # Count unique rows that were annotated and had allele flips
    processed_variants = len(updated_rows)
    flipped_count = len(flipped_rows)

    log.write(f" -Total unique sumstats rows annotated: {processed_variants:,}",
              verbose=verbose)
    log.write(f" -Total unique rows with allele flips: {flipped_count:,}",
              verbose=verbose)

    # ============================================================================
    # Step 16: Cleanup temporary lookup file if requested
    # ============================================================================
    # If lookup table was a temporary file created during processing,
    # remove it to free up disk space
    if rm_tmp_lookup:
        if isinstance(lookup_table, str) and os.path.exists(lookup_table):
            try:
                os.remove(lookup_table)
                log.write(f" -Cleaning lookup: {lookup_table}...", verbose=verbose)
            except:
                pass

    # ============================================================================
    # Step 17: Return updated sumstats
    # ============================================================================
    return sumstats

import subprocess
from pathlib import Path

def _convert_vcf_to_bcf(reference, 
                        threads=6, 
                        strip=True, 
                        ref_fa=None, 
                        log=Log(), 
                        verbose=True):
    """
    Normalize a reference VCF (multi-allelic splitting) with optional INFO/FORMAT stripping
    and optional left-normalization using a reference FASTA.

    Parameters
    ----------
    reference : str or Path
        Path to the reference VCF (bgzipped). If missing '.gz', it will be appended.
    threads : int
        Number of threads for bcftools processing.
    strip : bool, default=True
        If True, remove all INFO and FORMAT fields and name output as <reference>.strip.bcf.
        If False, keep INFO/FORMAT and name output <reference>.bcf.
    ref_fa : str or Path, optional
        Reference FASTA for left-normalization (`bcftools norm -f ref.fa`). If None, skip.
        FASTA must be indexed (.fai).
    log : Log
    verbose : bool

    Returns
    -------
    str
        Path to the generated .bcf file.
    """

    reference = str(reference)
    if not reference.endswith(".gz"):
        reference_vcf = reference + ".gz"
    else:
        reference_vcf = reference

    # Output name
    out_bcf = reference.replace(".vcf.gz","") + (".strip.bcf" if strip else ".bcf")

    # ---- Build normalization command ----
    # Always split multi-allelics
    cmd_norm = f"bcftools norm -m - {reference_vcf} -Ou"

    # Add FASTA for left normalization if provided
    if ref_fa is not None:
        cmd_norm += f" -f {ref_fa}"

    # ---- Strip or keep INFO/FORMAT ----
    if strip:
        cmd_annot = f" | bcftools annotate -x INFO,FORMAT -Ob --threads {threads} -o {out_bcf}"
    else:
        cmd_annot = f" | bcftools view -Ob --threads {threads} -o {out_bcf}"

    cmd1 = ["bash", "-c", cmd_norm + cmd_annot]
    if verbose:
        log.write(" -Running: {}...".format(cmd1[2]), verbose=verbose)
    subprocess.check_call(cmd1)

    # ---- Index ----
    cmd2 = ["bcftools", "index", "-f", out_bcf]
    if verbose:
        log.write(" -Running: {}...".format(" ".join(cmd2)), verbose=verbose)
    subprocess.check_call(cmd2)

    if verbose:
        log.write(" -Done. Output: {} and index...".format(out_bcf), verbose=verbose)
    return out_bcf


import os
import shutil
import subprocess
import tempfile
import pandas as pd
from multiprocessing import Pool


def _run_bcftools_extract(args):
    chr_name, chr_target_file, vcf_path, threads = args
    out_bcf = tempfile.NamedTemporaryFile(delete=False,
                                          suffix=f".gwaslab.chr{chr_name}.bcf").name

    cmd = [
        "bcftools", "view",
        "-T", chr_target_file,
        "-Ob", "-o", out_bcf,
        "--threads", str(threads),
        vcf_path
    ]
    subprocess.check_call(cmd)
    subprocess.check_call(["bcftools", "index", "-f", out_bcf])

    return out_bcf





def _expand_multiallelic_fast(df, ea_col, nea_col):
    """
    Expand multi-allelic lookup rows in a fast and memory-efficient way.
    A,T,G → 3 biallelic rows.
    """

    # Detect if any multiallelic rows exist
    if not df[ea_col].str.contains(",").any():
        return df  # nothing to do

    # ------------------------------
    # Split EA (ALT) column
    # ------------------------------
    ea_split = df[ea_col].str.split(",", expand=True)

    # Pre-split NEA/REF only if necessary
    if df[nea_col].str.contains(",").any():
        nea_split = df[nea_col].str.split(",", expand=True)
    else:
        # Avoid unnecessary expansion → more memory efficient
        nea_split = None

    out_frames = []

    # ------------------------------
    # Expand each ALT allele column
    # ------------------------------
    for i in range(ea_split.shape[1]):
        ea_i = ea_split[i]

        if ea_i.isna().all():
            continue  # no more alternate alleles

        df_i = df.copy()

        # Assign EA/ALT
        df_i[ea_col] = ea_i

        # Assign NEA/REF
        if nea_split is None:
            # Avoid allocating new arrays: repeat scalar NEA values
            pass
        else:
            df_i[nea_col] = nea_split[i]

        # Keep only rows where this allele exists
        df_i = df_i[ea_i.notna()]

        out_frames.append(df_i)

    # Concatenate efficiently
    return pd.concat(out_frames, ignore_index=True)
