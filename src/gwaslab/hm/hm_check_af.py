import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from .hm_assign_rsid import _annotate_sumstats
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

@with_logging(
    start_to_msg="check the difference between EAF (sumstats) and ALT frequency (reference VCF) using sweep mode",
    finished_msg="checking the difference between EAF (sumstats) and ALT frequency (reference VCF) using sweep mode",
    start_cols=["CHR","POS","EA","NEA","EAF","STATUS"],
    start_function=".check_af2()",
    must_kwargs=["ref_alt_freq"]
)
def _check_af_with_annotation(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    assign_cols: tuple = ("AF",),
    mapper: ChromosomeMapper | None = None,
    threads: int = 6,
    reuse_lookup: bool = True,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    eaf: str = "EAF",
    ref_alt_freq: str = "AF",
    column_name: str = "DAF",
    suffix: str = "",
    status: str = "STATUS",
    force: bool = False,
    log: "Log" = Log(),
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Check the difference between effect allele frequency (EAF) in summary statistics and 
    alternative allele frequency in reference VCF using sweep mode.
    
    This function performs a two-step process:
    1. **Annotation**: Annotates summary statistics with reference allele frequency (AF) from a 
       reference VCF/BCF or TSV file using `_annotate_sumstats`. The AF is stored in a column 
       (default: same as `ref_alt_freq`) for use in DAF calculation.
    2. **DAF Calculation**: Calculates the difference between EAF (sumstats) and ALT_AF (reference) 
       with proper allele matching. The difference (DAF) is stored in a new column.
    
    **DAF Calculation:**
    - DAF = EAF (sumstats) - ALT_AF (reference VCF)
    - Positive DAF: EAF in sumstats is higher than reference
    - Negative DAF: EAF in sumstats is lower than reference
    - Large |DAF| values (> 0.2) may indicate issues requiring investigation
    
    **When to use:**
    - After `infer_af2()` to validate inferred EAF values
    - Before harmonization to identify potential allele mismatches
    - For quality control to flag variants with unusual frequency differences
    
    Parameters
    ----------
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object. Must contain columns: CHR, POS, EA, NEA, EAF.
    path : str or None, optional
        Path to reference file (VCF/BCF or TSV). If provided, overrides `tsv_path`. 
        The function will automatically detect the file format.
    vcf_path : str or None, optional
        Path to VCF/BCF file. If provided, overrides both `path` and `tsv_path`.
    tsv_path : str or None, optional
        Path to precomputed lookup TSV file. If not provided and `vcf_path` is given, 
        a lookup table will be generated from the VCF/BCF file.
    assign_cols : tuple or str, default=("AF",)
        Column names to extract from reference file during annotation. The first column 
        will be used as the reference AF for DAF calculation.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper.
    threads : int, default=6
        Number of threads for parallel processing during annotation and lookup table generation.
    reuse_lookup : bool, default=True
        If True, reuse existing lookup table TSV file if found. If False, regenerate from VCF/BCF.
    convert_to_bcf : bool, default=False
        If True, convert VCF to BCF format for faster processing. Note: `strip_info` will be 
        set to False if converting to BCF (INFO fields are needed for AF extraction).
    strip_info : bool, default=True
        If True, strip INFO fields from VCF during lookup table generation to reduce file size.
        Set to False if INFO fields are needed for other purposes.
    chrom : str, default="CHR"
        Column name for chromosome in sumstats.
    pos : str, default="POS"
        Column name for position in sumstats.
    ea : str, default="EA"
        Column name for effect allele in sumstats.
    nea : str, default="NEA"
        Column name for non-effect allele in sumstats.
    eaf : str, default="EAF"
        Column name for effect allele frequency in sumstats.
    ref_alt_freq : str, default="AF"
        Field name for alternative allele frequency in VCF INFO section (e.g., "AF", "AF_popmax", 
        "gnomAD_AF"). This should match the field name in the reference VCF.
    column_name : str, default="DAF"
        Name of the column to store the difference values. The final column name will be 
        `column_name + suffix`.
    suffix : str, default=""
        Suffix to append to the column name (e.g., "_pop1", "_gnomad"). Useful when comparing 
        multiple reference populations.
    status : str, default="STATUS"
        Column name for status codes. By default, only processes variants with STATUS digit 4 = 0 
        (standardized and normalized), unless `force=True`.
    force : bool, default=False
        If True, check all variants regardless of STATUS codes. If False, only processes variants 
        with valid harmonization status (STATUS digit 4 = 0).
    verbose : bool, default=True
        If True, print progress messages and DAF statistics.
    log : Log, default=Log()
        Logging object for recording process information.
    
    Returns
    -------
    pd.DataFrame or Sumstats
        If input is a DataFrame, returns updated DataFrame with DAF column.
        If input is a Sumstats object, returns the Sumstats object with updated data.
    
    Notes
    -----
    - This function uses optimized bulk lookup methods for faster processing compared to 
      per-variant VCF queries.
    - The function automatically handles chromosome name mapping using ChromosomeMapper.
    - Lookup tables are cached as TSV files for faster subsequent runs when `reuse_lookup=True`.
    - The difference in allele frequency (DAF) is calculated as: DAF = EAF (sumstats) - ALT_AF (reference VCF)
    - **Important**: This DAF is NOT the derived allele frequency. It is simply the difference in 
      allele frequencies between two datasets.
    - Large |DAF| values (> 0.2) may indicate:
      - Population differences (expected for population-specific variants)
      - Allele mismatches (check EA/NEA alignment)
      - Strand flips (use `infer_strand()` to resolve)
      - Data quality issues (verify EAF calculation in sumstats)
    
    See Also
    --------
    _infer_af_with_annotation : Function that infers EAF from reference VCF before checking differences.
    _infer_strand_with_annotation : Function that infers strand orientation which may affect DAF values.
    """
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Normalize assign_cols to tuple/list
    if isinstance(assign_cols, str):
        assign_cols = tuple([assign_cols])
    elif isinstance(assign_cols, list):
        assign_cols = tuple(assign_cols)
    else:
        assign_cols = tuple(assign_cols)
    
    # Ensure assign_cols contains ref_alt_freq
    if ref_alt_freq not in assign_cols:
        assign_cols = tuple(list(assign_cols) + [ref_alt_freq])
    
    # Optimize: If converting to BCF, don't strip info (needed for AF extraction)
    if convert_to_bcf:
        strip_info = False

    # Get mapper from Sumstats object if available, otherwise create one
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chrom in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chrom])

    # Auto-detect reference format from VCF file
    if vcf_path is not None:
        mapper.detect_reference_format(vcf_path)
    elif path is not None:
        from gwaslab.hm.hm_assign_rsid import is_vcf_file
        if is_vcf_file(path):
            mapper.detect_reference_format(path)

    # First annotate with AF from reference
    annotated_sumstats = _annotate_sumstats(
        sumstats=sumstats,
        path=path,
        vcf_path=vcf_path,
        tsv_path=tsv_path,
        assign_cols=assign_cols,
        mapper=mapper,
        threads=threads,
        chrom=chrom,
        pos=pos,
        ea=ea,
        nea=nea,
        reuse_lookup=reuse_lookup,
        convert_to_bcf=convert_to_bcf,
        strip_info=strip_info,
        verbose=verbose,
        log=log,
    )
    
    # Calculate DAF = EAF - ALT_AF (reference)
    # The annotated AF column should contain ALT_AF from reference
    column_name = column_name + suffix
    
    # Initialize DAF column
    if column_name not in annotated_sumstats.columns:
        annotated_sumstats[column_name] = np.nan
    
    # Filter variants to process based on STATUS if not force
    if not force:
        from gwaslab.hm.hm_harmonize_sumstats import _extract_status_digit
        digit_4 = _extract_status_digit(annotated_sumstats[status], 4)
        good_chrpos = (digit_4 == 0)
        log.write(" -Checking variants:", good_chrpos.sum(), verbose=verbose)
    else:
        good_chrpos = pd.Series(True, index=annotated_sumstats.index)
        log.write(" -Checking all variants (force=True):", good_chrpos.sum(), verbose=verbose)
    
    # Only process variants with EAF available
    has_eaf = annotated_sumstats[eaf].notna()
    has_ref_af = annotated_sumstats[ref_alt_freq].notna()
    to_process = good_chrpos & has_eaf & has_ref_af
    
    if to_process.sum() > 0:
        # Calculate DAF = EAF - ALT_AF (reference)
        # The ref_alt_freq column contains ALT_AF from reference VCF
        # We need to match alleles: if EA matches ALT, use ALT_AF directly
        # If EA matches REF (ALLELE_FLIPPED=True), use 1 - ALT_AF
        # Check for ALLELE_FLIPPED flag to determine allele orientation
        flipped_col = "ALLELE_FLIPPED"
        if flipped_col in annotated_sumstats.columns:
            # If alleles were flipped, EA in sumstats matches REF in VCF
            # So we need to use 1 - ALT_AF to get the frequency of EA
            ref_af_matched = annotated_sumstats.loc[to_process, ref_alt_freq].copy()
            flipped_mask = annotated_sumstats.loc[to_process, flipped_col].fillna(False)
            ref_af_matched[flipped_mask] = 1 - ref_af_matched[flipped_mask]
        else:
            # No flip flag, assume EA matches ALT (use AF directly)
            ref_af_matched = annotated_sumstats.loc[to_process, ref_alt_freq]
        
        # DAF = EAF (sumstats) - ALT_AF (reference, matched to EA)
        eaf_values = annotated_sumstats.loc[to_process, eaf].values
        ref_af_values = ref_af_matched.values if hasattr(ref_af_matched, 'values') else ref_af_matched
        daf_values = eaf_values - ref_af_values
        annotated_sumstats.loc[to_process, column_name] = daf_values.astype(float)
        
        log.write(" -Difference in allele frequency (DAF) = EAF (sumstats) - ALT_AF (reference VCF)", verbose=verbose)
        log.write(" -Note: this DAF is not the derived allele frequency.", verbose=verbose)
        log.write(" - {} max:".format(column_name), np.nanmax(annotated_sumstats[column_name]), verbose=verbose)
        log.write(" - {} min:".format(column_name), np.nanmin(annotated_sumstats[column_name]), verbose=verbose)
        log.write(" - {} mean:".format(column_name), np.nanmean(annotated_sumstats[column_name]), verbose=verbose)
        log.write(" - {} sd:".format(column_name), np.nanstd(annotated_sumstats[column_name]), verbose=verbose)
        log.write(" - abs({}) min:".format(column_name), np.nanmin(np.abs(annotated_sumstats[column_name])), verbose=verbose)
        log.write(" - abs({}) max:".format(column_name), np.nanmax(np.abs(annotated_sumstats[column_name])), verbose=verbose)
        log.write(" - abs({}) sd:".format(column_name), np.nanstd(np.abs(annotated_sumstats[column_name])), verbose=verbose)
        log.write("Finished allele frequency checking!")
    
    # Drop ALLELE_FLIPPED as it's an internal temporary column
    if "ALLELE_FLIPPED" in annotated_sumstats.columns:
        annotated_sumstats = annotated_sumstats.drop(columns=["ALLELE_FLIPPED"])
    
    # Set metadata and update harmonization status if Sumstats object is available
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = annotated_sumstats
        try:
            from gwaslab.info.g_meta import _append_meta_record
            if path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_daf"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer_daf"], path)
            elif vcf_path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_daf"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer_daf"], vcf_path)
        except:
            pass
        return sumstats_obj.data
    else:
        return annotated_sumstats

@with_logging(
    start_to_msg="infer sumstats EAF using reference VCF ALT frequency in sweep mode",
    finished_msg="inferring sumstats EAF using reference VCF ALT frequency in sweep mode",
    start_cols=["CHR","POS","EA","NEA","STATUS"],
    start_function=".infer_af2()",
    must_kwargs=["ref_alt_freq"]
)
def _infer_af_with_annotation(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    assign_cols: tuple = ("AF",),
    mapper: ChromosomeMapper | None = None,
    threads: int = 6,
    reuse_lookup: bool = True,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    eaf: str = "EAF",
    ref_alt_freq: str = "AF",
    status: str = "STATUS",
    force: bool = False,
    log: "Log" = Log(),
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Infer effect allele frequency (EAF) in summary statistics using reference VCF ALT frequency 
    in sweep mode.
    
    This function performs a two-step process:
    1. **Annotation**: Annotates summary statistics with reference allele frequency (AF) from a 
       reference VCF/BCF or TSV file using `_annotate_sumstats`. The AF is stored in a column 
       (default: same as `ref_alt_freq`) for use in EAF inference.
    2. **EAF Inference**: Infers EAF values by matching variants and handling allele orientation. 
       If EA matches ALT in VCF, uses ALT_AF directly. If EA matches REF in VCF, uses 1 - ALT_AF.
    
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
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object. Must contain columns: CHR, POS, EA, NEA.
    path : str or None, optional
        Path to reference file (VCF/BCF or TSV). If provided, overrides `tsv_path`. 
        The function will automatically detect the file format.
    vcf_path : str or None, optional
        Path to VCF/BCF file. If provided, overrides both `path` and `tsv_path`.
    tsv_path : str or None, optional
        Path to precomputed lookup TSV file. If not provided and `vcf_path` is given, 
        a lookup table will be generated from the VCF/BCF file.
    assign_cols : tuple or str, default=("AF",)
        Column names to extract from reference file during annotation. The first column 
        will be used as the reference AF for EAF inference.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper.
    threads : int, default=6
        Number of threads for parallel processing during annotation and lookup table generation.
    reuse_lookup : bool, default=True
        If True, reuse existing lookup table TSV file if found. If False, regenerate from VCF/BCF.
    convert_to_bcf : bool, default=False
        If True, convert VCF to BCF format for faster processing. Note: `strip_info` will be 
        set to False if converting to BCF (INFO fields are needed for AF extraction).
    strip_info : bool, default=True
        If True, strip INFO fields from VCF during lookup table generation to reduce file size.
        Set to False if INFO fields are needed for other purposes.
    chrom : str, default="CHR"
        Column name for chromosome in sumstats.
    pos : str, default="POS"
        Column name for position in sumstats.
    ea : str, default="EA"
        Column name for effect allele in sumstats.
    nea : str, default="NEA"
        Column name for non-effect allele in sumstats.
    eaf : str, default="EAF"
        Column name for effect allele frequency. This column will be created if it doesn't exist, 
        and existing values will be updated where inference is successful.
    ref_alt_freq : str, default="AF"
        Field name for alternative allele frequency in VCF INFO section (e.g., "AF", "AF_popmax", 
        "gnomAD_AF"). This should match the field name in the reference VCF.
    status : str, default="STATUS"
        Column name for status codes. By default, only processes variants with STATUS digit 4 = 0 
        (standardized and normalized), unless `force=True`.
    force : bool, default=False
        If True, infer EAF for all variants regardless of STATUS codes. If False, only processes 
        variants with valid harmonization status (STATUS digit 4 = 0).
    verbose : bool, default=True
        If True, print progress messages and statistics about inference success rate.
    log : Log, default=Log()
        Logging object for recording process information.
    
    Returns
    -------
    pd.DataFrame or Sumstats
        If input is a DataFrame, returns updated DataFrame with inferred EAF values.
        If input is a Sumstats object, returns the Sumstats object with updated data.
    
    Notes
    -----
    - This function uses optimized bulk lookup methods for faster processing compared to 
      per-variant VCF queries.
    - The function automatically handles chromosome name mapping using ChromosomeMapper.
    - Lookup tables are cached as TSV files for faster subsequent runs when `reuse_lookup=True`.
    - By default, only processes variants with valid harmonization status (STATUS digit 4 = 0) to 
      ensure alleles are standardized and normalized.
    - After inference, the function reports statistics about:
      - Number of variants with EAF successfully inferred
      - Number of variants still missing EAF (not found in reference or missing ref_alt_freq field)
    - The inferred EAF values are stored in the specified EAF column, overwriting existing values 
      where inference is successful.
    
    See Also
    --------
    _check_af_with_annotation : Calculate difference between EAF (sumstats) and ALT_AF (reference) after inference.
    _infer_strand_with_annotation : Function that infers strand orientation which may affect EAF values.
    """
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Normalize assign_cols to tuple/list
    if isinstance(assign_cols, str):
        assign_cols = tuple([assign_cols])
    elif isinstance(assign_cols, list):
        assign_cols = tuple(assign_cols)
    else:
        assign_cols = tuple(assign_cols)
    
    # Ensure assign_cols contains ref_alt_freq
    if ref_alt_freq not in assign_cols:
        assign_cols = tuple(list(assign_cols) + [ref_alt_freq])
    
    # Optimize: If converting to BCF, don't strip info (needed for AF extraction)
    if convert_to_bcf:
        strip_info = False

    # Get mapper from Sumstats object if available, otherwise create one
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chrom in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chrom])

    # Auto-detect reference format from VCF file
    if vcf_path is not None:
        mapper.detect_reference_format(vcf_path)
    elif path is not None:
        from gwaslab.hm.hm_assign_rsid import is_vcf_file
        if is_vcf_file(path):
            mapper.detect_reference_format(path)

    # Initialize EAF column if it doesn't exist
    if eaf not in sumstats.columns:
        sumstats[eaf] = np.nan
    
    prenumber = sumstats[eaf].isna().sum()
    
    # Filter variants to process based on STATUS if not force
    if not force:
        from gwaslab.hm.hm_harmonize_sumstats import _extract_status_digit
        digit_4 = _extract_status_digit(sumstats[status], 4)
        good_chrpos = (digit_4 == 0)
        log.write(" -Checking variants:", good_chrpos.sum(), verbose=verbose)
    else:
        good_chrpos = pd.Series(True, index=sumstats.index)
        log.write(" -Checking all variants (force=True):", good_chrpos.sum(), verbose=verbose)
    
    # Only process variants that need EAF inference
    needs_eaf = sumstats[eaf].isna()
    to_process = good_chrpos & needs_eaf
    
    log.write(f" -Variants needing EAF inference: {to_process.sum()} / {len(sumstats)}", verbose=verbose)
    
    if to_process.sum() > 0:
        # Annotate with AF from reference for variants that need EAF
        # We'll annotate all variants first, then update only those that need EAF
        annotated_sumstats = _annotate_sumstats(
            sumstats=sumstats,
            path=path,
            vcf_path=vcf_path,
            tsv_path=tsv_path,
            assign_cols=assign_cols,
            mapper=mapper,
            threads=threads,
            chrom=chrom,
            pos=pos,
            ea=ea,
            nea=nea,
            reuse_lookup=reuse_lookup,
            convert_to_bcf=convert_to_bcf,
            strip_info=strip_info,
            verbose=verbose,
            log=log,
        )
        
        # Debug: Check if ALLELE_FLIPPED exists after annotation
        if "ALLELE_FLIPPED" not in annotated_sumstats.columns:
            log.write(" -WARNING: ALLELE_FLIPPED column not found after _annotate_sumstats!", verbose=verbose)
        else:
            flipped_count = annotated_sumstats["ALLELE_FLIPPED"].fillna(False).sum()
            log.write(f" -ALLELE_FLIPPED column found, {flipped_count} variants marked as flipped", verbose=verbose)
        
        # The annotated AF column contains ALT_AF from reference VCF
        # _annotate_sumstats handles allele matching, so the annotated AF should be
        # the frequency of the allele that matches EA in sumstats
        # For inference: if EA matches ALT, use ALT_AF directly; if EA matches REF, use 1 - ALT_AF
        # Since _annotate_sumstats already handles matching, we can use the annotated AF directly
        # as the inferred EAF for variants where EAF was missing
        
        # Update EAF for variants that were missing EAF and now have annotated AF
        has_ref_af = annotated_sumstats[ref_alt_freq].notna()
        to_update = to_process & has_ref_af
        
        if to_update.sum() > 0:
            # Use the annotated AF as inferred EAF, but adjust for allele orientation
            # If EA matches ALT, use ALT_AF directly
            # If EA matches REF (ALLELE_FLIPPED=True), use 1 - ALT_AF
            flipped_col = "ALLELE_FLIPPED"
            if flipped_col in annotated_sumstats.columns:
                # If alleles were flipped, EA in sumstats matches REF in VCF
                # So we need to use 1 - ALT_AF to get the frequency of EA
                inferred_eaf = annotated_sumstats.loc[to_update, ref_alt_freq].copy()
                flipped_mask = annotated_sumstats.loc[to_update, flipped_col].fillna(False)
                inferred_eaf[flipped_mask] = 1 - inferred_eaf[flipped_mask]
                annotated_sumstats.loc[to_update, eaf] = inferred_eaf
            else:
                # No flip flag, assume EA matches ALT (use AF directly)
                annotated_sumstats.loc[to_update, eaf] = annotated_sumstats.loc[to_update, ref_alt_freq]
        
        afternumber = annotated_sumstats[eaf].isna().sum()
        log.write(" -Inferred EAF for {} variants.".format(prenumber - afternumber), verbose=verbose)
        log.write(" -EAF is still missing for {} variants.".format(afternumber), verbose=verbose)
        
        sumstats = annotated_sumstats
    else:
        log.write(" -No variants need EAF inference (all have EAF or don't meet criteria).", verbose=verbose)
    
    # Drop ALLELE_FLIPPED as it's an internal temporary column
    if "ALLELE_FLIPPED" in sumstats.columns:
        sumstats = sumstats.drop(columns=["ALLELE_FLIPPED"])
    
    # Set metadata and update harmonization status if Sumstats object is available
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _append_meta_record
            if path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_af"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer_af"], path)
            elif vcf_path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_af"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer_af"], vcf_path)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

@with_logging(
    start_to_msg="infer sumstats EAF from sumstats MAF using reference VCF ALT frequency in sweep mode",
    finished_msg="inferring sumstats EAF from sumstats MAF using reference VCF ALT frequency in sweep mode",
    start_cols=["CHR","POS","EA","NEA","MAF","STATUS"],
    start_function=".infer_eaf_from_maf2()",
    must_kwargs=["ref_alt_freq"]
)
def _infer_af_with_maf_annotation(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    assign_cols: tuple = ("AF",),
    mapper: ChromosomeMapper | None = None,
    threads: int = 6,
    reuse_lookup: bool = True,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    eaf: str = "EAF",
    maf: str = "MAF",
    ref_alt_freq: str = "AF",
    ref_eaf: str = "_REF_EAF",
    status: str = "STATUS",
    force: bool = False,
    log: "Log" = Log(),
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Infer effect allele frequency (EAF) in summary statistics from MAF using reference VCF ALT frequency 
    in sweep mode.
    
    This function performs a two-step process:
    1. **Annotation**: Annotates summary statistics with reference allele frequency (AF) from a 
       reference VCF/BCF or TSV file using `_annotate_sumstats`. The AF is stored in a temporary 
       column (default: "_REF_EAF") for use in EAF inference.
    2. **EAF Inference**: Infers EAF values by comparing reference AF with MAF from sumstats. 
       If the reference allele frequency and MAF suggest different major/minor alleles, the EAF is 
       calculated as 1 - MAF (flipping the allele), otherwise MAF is used directly.
    
    **Workflow:**
    1. Matches variants in sumstats with reference VCF by CHR:POS:EA:NEA.
    2. Extracts ALT frequency from VCF INFO field (specified by `ref_alt_freq`).
    3. Compares reference AF with MAF to determine if flipping is needed:
       - If ref_AF >= 0.5 (ref allele is major) != MAF > 0.5 (MAF is minor), then flip (use 1 - MAF)
       - Otherwise use MAF directly
    4. Updates EAF column in sumstats with inferred values.
    
    **When to use:**
    - When sumstats have MAF values but are missing EAF values.
    - When you want to infer EAF from MAF using a reference population (e.g., 1000 Genomes, gnomAD).
    - As an alternative to `infer_af2()` when you have MAF but not direct EAF information.
    
    Parameters
    ----------
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object. Must contain columns: CHR, POS, EA, NEA, MAF.
    path : str or None, optional
        Path to reference file (VCF/BCF or TSV). If provided, overrides `tsv_path`. 
        The function will automatically detect the file format.
    vcf_path : str or None, optional
        Path to VCF/BCF file. If provided, overrides both `path` and `tsv_path`.
    tsv_path : str or None, optional
        Path to precomputed lookup TSV file. If not provided and `vcf_path` is given, 
        a lookup table will be generated from the VCF/BCF file.
    assign_cols : tuple or str, default=("AF",)
        Column names to extract from reference file during annotation. The first column 
        will be used as the reference AF for EAF inference.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper.
    threads : int, default=6
        Number of threads for parallel processing during annotation and lookup table generation.
    reuse_lookup : bool, default=True
        If True, reuse existing lookup table TSV file if found. If False, regenerate from VCF/BCF.
    convert_to_bcf : bool, default=False
        If True, convert VCF to BCF format for faster processing. Note: `strip_info` will be 
        set to False if converting to BCF (INFO fields are needed for AF extraction).
    strip_info : bool, default=True
        If True, strip INFO fields from VCF during lookup table generation to reduce file size.
        Set to False if INFO fields are needed for other purposes.
    chrom : str, default="CHR"
        Column name for chromosome in sumstats.
    pos : str, default="POS"
        Column name for position in sumstats.
    ea : str, default="EA"
        Column name for effect allele in sumstats.
    nea : str, default="NEA"
        Column name for non-effect allele in sumstats.
    eaf : str, default="EAF"
        Column name for effect allele frequency. This column will be created if it doesn't exist, 
        and existing values will be updated where inference is successful.
    maf : str, default="MAF"
        Column name for minor allele frequency in sumstats. This is required for EAF inference.
    ref_alt_freq : str, default="AF"
        Field name for alternative allele frequency in VCF INFO section (e.g., "AF", "AF_popmax", 
        "gnomAD_AF"). This should match the field name in the reference VCF.
    ref_eaf : str, default="_REF_EAF"
        Temporary column name for storing reference allele frequency. This column will be created 
        during annotation and dropped before returning.
    status : str, default="STATUS"
        Column name for status codes. By default, only processes variants with STATUS digit 4 = 0 
        (standardized and normalized), unless `force=True`.
    force : bool, default=False
        If True, infer EAF for all variants regardless of STATUS codes. If False, only processes 
        variants with valid harmonization status (STATUS digit 4 = 0).
    verbose : bool, default=True
        If True, print progress messages and statistics about inference success rate.
    log : Log, default=Log()
        Logging object for recording process information.
    
    Returns
    -------
    pd.DataFrame or Sumstats
        If input is a DataFrame, returns updated DataFrame with inferred EAF values.
        If input is a Sumstats object, returns the Sumstats object with updated data.
    
    Notes
    -----
    - This function uses optimized bulk lookup methods for faster processing compared to 
      per-variant VCF queries.
    - The function automatically handles chromosome name mapping using ChromosomeMapper.
    - Lookup tables are cached as TSV files for faster subsequent runs when `reuse_lookup=True`.
    - By default, only processes variants with valid harmonization status (STATUS digit 4 = 0) to 
      ensure alleles are standardized and normalized.
    - The function requires MAF values in sumstats to infer EAF.
    - After inference, the function reports statistics about:
      - Number of variants with EAF successfully inferred
      - Number of variants still missing EAF (not found in reference or missing ref_alt_freq field)
    - The inferred EAF values are stored in the specified EAF column, overwriting existing values 
      where inference is successful.
    - The temporary reference EAF column (`ref_eaf`) is dropped before returning.
    
    See Also
    --------
    _infer_af_with_annotation : Function that infers EAF directly from reference ALT_AF.
    _check_af_with_annotation : Calculate difference between EAF (sumstats) and ALT_AF (reference) after inference.
    """
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats, pd.DataFrame):
        # Called with DataFrame
        is_dataframe = True
    else:
        # Called with Sumstats object
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Normalize assign_cols to tuple/list
    if isinstance(assign_cols, str):
        assign_cols = tuple([assign_cols])
    elif isinstance(assign_cols, list):
        assign_cols = tuple(assign_cols)
    else:
        assign_cols = tuple(assign_cols)
    
    # Ensure assign_cols contains ref_alt_freq
    if ref_alt_freq not in assign_cols:
        assign_cols = tuple(list(assign_cols) + [ref_alt_freq])
    
    # Optimize: If converting to BCF, don't strip info (needed for AF extraction)
    if convert_to_bcf:
        strip_info = False

    # Get mapper from Sumstats object if available, otherwise create one
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect sumstats format if data is available
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chrom in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chrom])

    # Auto-detect reference format from VCF file
    if vcf_path is not None:
        mapper.detect_reference_format(vcf_path)
    elif path is not None:
        from gwaslab.hm.hm_assign_rsid import is_vcf_file
        if is_vcf_file(path):
            mapper.detect_reference_format(path)

    # Initialize EAF and ref_eaf columns if they don't exist
    if eaf not in sumstats.columns:
        sumstats[eaf] = np.nan
    if ref_eaf not in sumstats.columns:
        sumstats[ref_eaf] = np.nan
    
    prenumber = sumstats[eaf].isna().sum()
    
    # Filter variants to process based on STATUS if not force
    if not force:
        from gwaslab.hm.hm_harmonize_sumstats import _extract_status_digit
        digit_4 = _extract_status_digit(sumstats[status], 4)
        good_chrpos = (digit_4 == 0)
    else:
        good_chrpos = pd.Series(True, index=sumstats.index)
    
    # Only process variants with valid MAF (required for EAF inference)
    good_chrpos = good_chrpos & sumstats[maf].notna()
    log.write(" -Checking variants:", good_chrpos.sum(), verbose=verbose)
    
    if good_chrpos.sum() > 0:
        # Annotate with AF from reference for variants that need EAF inference
        annotated_sumstats = _annotate_sumstats(
            sumstats=sumstats,
            path=path,
            vcf_path=vcf_path,
            tsv_path=tsv_path,
            assign_cols=assign_cols,
            mapper=mapper,
            threads=threads,
            chrom=chrom,
            pos=pos,
            ea=ea,
            nea=nea,
            reuse_lookup=reuse_lookup,
            convert_to_bcf=convert_to_bcf,
            strip_info=strip_info,
            verbose=verbose,
            log=log,
        )
        
        # Store the annotated AF in the ref_eaf column
        # Handle allele matching: if EA matches ALT, use ALT_AF directly
        # If EA matches REF (ALLELE_FLIPPED=True), use 1 - ALT_AF
        flipped_col = "ALLELE_FLIPPED"
        if flipped_col in annotated_sumstats.columns:
            # If alleles were flipped, EA in sumstats matches REF in VCF
            # So we need to use 1 - ALT_AF to get the frequency of EA
            ref_af_matched = annotated_sumstats[ref_alt_freq].copy()
            flipped_mask = annotated_sumstats[flipped_col].fillna(False)
            ref_af_matched[flipped_mask] = 1 - ref_af_matched[flipped_mask]
            annotated_sumstats[ref_eaf] = ref_af_matched
        else:
            # No flip flag, assume EA matches ALT (use AF directly)
            annotated_sumstats[ref_eaf] = annotated_sumstats[ref_alt_freq]
        
        # Infer sumstats EAF based on sumstats MAF and reference EAF
        # Use XOR logic - flip when ref_eaf and maf indicate different major/minor alleles
        # (ref_eaf >= 0.5) means ref allele is major, (maf <= 0.5) means maf is minor
        # If ref is major and maf is minor, or ref is minor and maf is major, we need to flip
        has_ref_eaf = annotated_sumstats[ref_eaf].notna()
        to_process = good_chrpos & has_ref_eaf
        
        if to_process.sum() > 0:
            is_flipped = (annotated_sumstats.loc[to_process, ref_eaf] >= 0.5) != (annotated_sumstats.loc[to_process, maf] > 0.5)
            # Vectorized assignment: use np.where for efficient conditional assignment
            maf_values = annotated_sumstats.loc[to_process, maf].values
            inferred_eaf = np.where(
                is_flipped,
                1 - maf_values,
                maf_values
            ).astype(float)
            # Cast to match original column dtype to avoid FutureWarning
            original_dtype = annotated_sumstats[eaf].dtype
            inferred_eaf = inferred_eaf.astype(original_dtype)
            annotated_sumstats.loc[to_process, eaf] = inferred_eaf
            log.write(" -Flipping MAF to obtain EAF for {} variants".format(is_flipped.sum()), verbose=verbose)
        
        afternumber = annotated_sumstats[eaf].isna().sum()
        log.write(" -Inferred EAF for {} variants.".format(prenumber - afternumber), verbose=verbose)
        log.write(" -EAF is still missing for {} variants.".format(afternumber), verbose=verbose)
        
        # Drop the temporary ref_eaf column
        if ref_eaf in annotated_sumstats.columns:
            annotated_sumstats = annotated_sumstats.drop(columns=[ref_eaf])
        
        sumstats = annotated_sumstats
    else:
        log.write(" -No variants need EAF inference (all have EAF or don't meet criteria).", verbose=verbose)
    
    # Drop ALLELE_FLIPPED as it's an internal temporary column
    if "ALLELE_FLIPPED" in sumstats.columns:
        sumstats = sumstats.drop(columns=["ALLELE_FLIPPED"])
    
    # Set metadata and update harmonization status if Sumstats object is available
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = sumstats
        try:
            from gwaslab.info.g_meta import _append_meta_record
            if path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_maf"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer_maf"], path)
            elif vcf_path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer_maf"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer_maf"], vcf_path)
        except:
            pass
        return sumstats_obj.data
    else:
        return sumstats

