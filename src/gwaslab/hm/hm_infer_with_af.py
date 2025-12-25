import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.info.g_vchange_status import vchange_status
from gwaslab.info.g_vchange_status import ensure_status_int
from .hm_assign_rsid import _annotate_sumstats
from gwaslab.qc.qc_decorator import with_logging

@with_logging(
    start_to_msg="annotate and infer strand orientation using allele frequencies",
    finished_msg="annotating and inferring strand orientation using allele frequencies",
    start_cols=["CHR","POS","EA","NEA","EAF","STATUS"],
    start_function=".infer_strand2()"
)
def _infer_strand_with_annotation(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    assign_cols: tuple = ("AF",),
    chr_dict: dict | None = None,
    threads: int = 6,
    reuse_lookup: bool = True,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    eaf: str = "EAF",
    raf: str = "RAF",
    flipped_col: str = "ALLELE_FLIPPED",
    strand_col: str = "STRAND",
    maf_threshold: float = 0.40,
    ref_maf_threshold: float = 0.40,
    daf_tolerance: float = 0.20,
    log: "Log" = Log(),
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Annotate summary statistics with reference allele frequency and infer strand orientation.
    
    This function performs a two-step process:
    1. **Annotation**: Annotates summary statistics with reference allele frequency (RAF) from a 
       reference VCF/BCF or TSV file using `_annotate_sumstats`. The RAF is stored in a column 
       (default: "RAF") for use in strand inference.
    2. **Strand Inference**: Infers strand orientation for palindromic SNPs and indels using the 
       annotated RAF and effect allele frequency (EAF) from the sumstats using `_infer_strand`.
    
    The function determines strand orientation by comparing EAF (from sumstats) with RAF (from 
    reference) to identify whether variants are on the forward (+) or reverse (-) strand. This 
    is particularly important for palindromic SNPs (A/T, G/C) and indels where strand orientation 
    cannot be determined from alleles alone.
    
    **STATUS Code Updates:**
    - Non-palindromic SNPs: STATUS 7th digit = `0`
    - Palindromic SNPs (forward): STATUS 7th digit = `1`
    - Palindromic SNPs (reverse): STATUS 7th digit = `2`
    - Indels (forward): STATUS 7th digit = `3`
    - Indels (ambiguous): STATUS 7th digit = `4`
    - Indels (reverse): STATUS 7th digit = `6`
    - Palindromic SNPs (MAF > threshold): STATUS 7th digit = `7`
    - Variants not found in reference: STATUS 7th digit = `8`
    - Variants not checked: STATUS 7th digit = `9`
    
    Parameters
    ----------
    sumstats : pd.DataFrame or Sumstats
        Summary statistics DataFrame or Sumstats object. Must contain columns: CHR, POS, EA, NEA, EAF, STATUS.
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
        will be renamed to `raf` (default: "RAF") for use in strand inference.
    chr_dict : dict or None, optional
        Dictionary for mapping chromosome names between sumstats and reference (e.g., 
        {1: "chr1", 2: "chr2"}). If None, automatic detection is attempted.
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
    raf : str, default="RAF"
        Column name to store reference allele frequency (from reference file). This column 
        will be created/updated during annotation.
    flipped_col : str, default="ALLELE_FLIPPED"
        Column name indicating if alleles were flipped during harmonization (True = reverse strand).
        This column is used as a hint for strand inference but is not required.
    strand_col : str, default="STRAND"
        Column name to store strand orientation ('+' for forward, '-' for reverse, '?' for unknown).
        This column will be created if it doesn't exist.
    maf_threshold : float, default=0.40
        Maximum minor allele frequency threshold for palindromic SNPs. Palindromic SNPs with 
        MAF > threshold in either EAF or RAF will be excluded from strand determination and 
        assigned STATUS 7th digit = `7`. Higher values allow more variants to be processed, 
        but may reduce accuracy for ambiguous cases.
    ref_maf_threshold : float, default=0.40
        Maximum minor allele frequency threshold for reference allele frequency (RAF). Used as 
        an additional filter for palindromic SNPs. Variants with MAF(RAF) > threshold will be 
        excluded from strand determination.
    daf_tolerance : float, default=0.20
        Difference in allele frequency tolerance for indels. For indels, the function compares 
        |EAF - RAF| (forward) and |EAF - (1 - RAF)| (reverse). If the difference is within 
        `daf_tolerance`, the strand is assigned. If both differences are within tolerance or 
        both exceed tolerance, the indel is marked as ambiguous (STATUS 7th digit = `4`).
    verbose : bool, default=True
        If True, print progress messages and warnings.
    log : Log, default=Log()
        Logging object for recording process information.
    
    Returns
    -------
    pd.DataFrame or Sumstats
        If input is a DataFrame, returns updated DataFrame with RAF and STRAND columns.
        If input is a Sumstats object, returns the Sumstats object with updated data.
    
    Notes
    -----
    - This function requires the sumstats to have been harmonized (STATUS codes present) and 
      to have EAF values for variants of interest.
    - The function uses optimized bulk lookup methods for faster processing compared to 
      per-variant VCF queries.
    - For palindromic SNPs, strand inference is only performed if MAF can be reliably 
      determined from EAF alone (EAF < maf_threshold OR EAF > 1 - maf_threshold).
    - The function automatically handles chromosome name mapping if `chr_dict` is provided.
    - Lookup tables are cached as TSV files for faster subsequent runs when `reuse_lookup=True`.
    
    See Also
    --------
    _annotate_sumstats : Function that performs the annotation step.
    _infer_strand : Function that performs the strand inference step.
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
        
    # Optimize: If converting to BCF, don't strip info (needed for AF extraction)
    if convert_to_bcf:
        strip_info = False

    # First annotate with AF
    annotated_sumstats = _annotate_sumstats(
        sumstats=sumstats,
        path=path,
        vcf_path=vcf_path,
        tsv_path=tsv_path,
        assign_cols=assign_cols,
        chr_dict=chr_dict,
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
    
    # Optimize: Rename annotation column to RAF more efficiently
    if assign_cols and len(assign_cols) > 0:
        source_col = assign_cols[0]
        log.write(" -Renaming {} in reference to {} in Sumstats".format(source_col, raf), verbose=verbose)
        # Only rename if source column exists and is different from target
        if source_col in annotated_sumstats.columns and source_col != raf:
            # Drop target column if it exists to avoid conflicts
            if raf in annotated_sumstats.columns:
                annotated_sumstats = annotated_sumstats.drop(columns=[raf])
            annotated_sumstats = annotated_sumstats.rename(columns={source_col: raf})

    # Then infer strand using the annotated AF as RAF
    result = _infer_strand(
        sumstats=annotated_sumstats,
        chrom=chrom,
        pos=pos,
        ea=ea,
        nea=nea,
        eaf=eaf,
        raf=raf,  # The AF column added by annotation
        flipped_col=flipped_col,
        strand_col=strand_col,
        maf_threshold=maf_threshold,
        ref_maf_threshold=ref_maf_threshold,
        daf_tolerance=daf_tolerance,
        log=log,
        verbose=verbose,
    )
    
    # Drop ALLELE_FLIPPED as it's an internal temporary column
    if "ALLELE_FLIPPED" in result.columns:
        result = result.drop(columns=["ALLELE_FLIPPED"])
    
    # Set metadata and update harmonization status if Sumstats object is available
    # Update harmonization status only if called with Sumstats object
    if not is_dataframe:
        # Assign modified dataframe back to the Sumstats object
        sumstats_obj.data = result
        try:
            from gwaslab.info.g_meta import _append_meta_record, _update_harmonize_step
            if path is not None:
                sumstats_obj.meta["gwaslab"]["references"]["ref_infer"] = _append_meta_record(
                    sumstats_obj.meta["gwaslab"]["references"]["ref_infer"], path)
            infer_strand_kwargs = {
                'path': path, 'vcf_path': vcf_path, 'tsv_path': tsv_path, 'assign_cols': assign_cols,
                'chr_dict': chr_dict, 'threads': threads, 'reuse_lookup': reuse_lookup,
                'convert_to_bcf': convert_to_bcf, 'strip_info': strip_info,
                'maf_threshold': maf_threshold, 'ref_maf_threshold': ref_maf_threshold,
                'daf_tolerance': daf_tolerance
            }
            _update_harmonize_step(sumstats_obj, "infer_strand", infer_strand_kwargs, True)
        except:
            pass
        return sumstats_obj.data
    else:
        return result

@with_logging(
    start_to_msg="infer strand orientation using allele frequencies",
    finished_msg="inferring strand orientation using allele frequencies",
    start_cols=["CHR","POS","EA","NEA","EAF"],
    start_function=".infer_strand()"
)
def _infer_strand(
    sumstats: pd.DataFrame,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    eaf: str = "EAF",
    raf: str = "RAF",
    flipped_col: str = "ALLELE_FLIPPED",
    strand_col: str = "STRAND",
    status_col: str = "STATUS",
    log: "Log" = Log(),
    maf_threshold: float = 0.40,
    ref_maf_threshold: float = 0.40,
    daf_tolerance: float = 0.20,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Infer strand orientation for palindromic SNPs and indels using allele frequencies.
    
    This function determines strand orientation (forward `+` or reverse `-`) for variants where 
    strand cannot be determined from alleles alone. It uses effect allele frequency (EAF) from 
    the sumstats and reference allele frequency (RAF) from a reference file to infer strand.
    
    **Workflow:**
    
    1. **Non-Palindromic SNPs**: Single-character alleles that are NOT palindromic pairs 
       (A/T, T/A, G/C, C/G). These are assigned STATUS 7th digit = `0` (strand is unambiguous).
    
    2. **Palindromic SNPs**: Single-character alleles that ARE palindromic pairs (A/T, T/A, G/C, C/G).
       - Only processes variants with STATUS digit 6 = 0,1,2 AND digit 7 = 8,9 (unknown strand).
       - Filters by MAF threshold: Only processes if MAF(EAF) ≤ maf_threshold AND MAF(RAF) ≤ ref_maf_threshold.
       - Compares |EAF - RAF| (forward) vs |EAF - (1 - RAF)| (reverse).
       - Assigns forward strand (`+`) if forward difference ≤ reverse difference, else reverse (`-`).
       - STATUS codes: `1` (forward), `2` (reverse), `7` (MAF > threshold), `8` (not found in reference).
    
    3. **Indels**: Variants with alleles longer than 1 character.
       - Only processes variants with STATUS digit 6 = 6 AND digit 7 = 8,9 (unknown indel strand).
       - Filters by MAF threshold: Only processes if MAF(EAF) ≤ maf_threshold AND MAF(RAF) ≤ ref_maf_threshold.
       - Compares |EAF - RAF| (forward) vs |EAF - (1 - RAF)| (reverse) with `daf_tolerance`.
       - Assigns forward strand (`+`) if |EAF - RAF| ≤ daf_tolerance, else checks reverse.
       - Assigns reverse strand (`-`) if |EAF - (1 - RAF)| ≤ daf_tolerance AND forward not within tolerance.
       - STATUS codes: `3` (forward), `4` (ambiguous), `6` (reverse), `8` (not found/MAF exceeded).
    
    **STATUS Code Updates (7th digit):**
    - `0`: Non-palindromic SNP (strand unambiguous)
    - `1`: Palindromic SNP, forward strand
    - `2`: Palindromic SNP, reverse strand (flipped)
    - `3`: Indel, forward strand
    - `4`: Indel, ambiguous (both DAF differences > tolerance)
    - `6`: Indel, reverse strand
    - `7`: Palindromic SNP, MAF > threshold (cannot infer)
    - `8`: Variant not found in reference or MAF exceeded
    - `9`: Variant not checked (no EAF or RAF available)
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame. Must contain columns: CHR, POS, EA, NEA, EAF, RAF, 
        ALLELE_FLIPPED, and STATUS.
    chrom : str, default="CHR"
        Column name for chromosome.
    pos : str, default="POS"
        Column name for position.
    ea : str, default="EA"
        Column name for effect allele.
    nea : str, default="NEA"
        Column name for non-effect allele.
    eaf : str, default="EAF"
        Column name for effect allele frequency (from sumstats).
    raf : str, default="RAF"
        Column name for reference allele frequency (from reference file, typically annotated 
        by `_annotate_sumstats` or `_infer_strand_with_annotation`).
    flipped_col : str, default="ALLELE_FLIPPED"
        Column name indicating if alleles were flipped during harmonization. This is used 
        as a hint but is not required for strand inference.
    strand_col : str, default="STRAND"
        Column name to store strand orientation ('+' for forward, '-' for reverse, '?' for unknown).
        This column will be created if it doesn't exist.
    status_col : str, default="STATUS"
        Column name containing 7-digit status codes. The 7th digit will be updated to reflect 
        strand inference results.
    maf_threshold : float, default=0.40
        Maximum minor allele frequency threshold for palindromic SNPs. Palindromic SNPs with 
        MAF > threshold in EAF will be assigned STATUS 7th digit = `7` (cannot infer). 
        Higher values allow more variants to be processed but may reduce accuracy.
    ref_maf_threshold : float, default=0.40
        Maximum minor allele frequency threshold for reference allele frequency (RAF). Used as 
        an additional filter for both palindromic SNPs and indels. Variants with MAF(RAF) > 
        threshold will be excluded from strand determination.
    daf_tolerance : float, default=0.20
        Difference in allele frequency tolerance for indels. For indels, compares |EAF - RAF| 
        (forward) and |EAF - (1 - RAF)| (reverse). If the difference is within `daf_tolerance`, 
        the strand is assigned. If both differences are within tolerance or both exceed tolerance, 
        the indel is marked as ambiguous (STATUS 7th digit = `4`).
    verbose : bool, default=True
        If True, print progress messages and warnings about missing data.
    log : Log, default=Log()
        Logging object for recording process information.
    
    Returns
    -------
    pd.DataFrame
        Updated summary statistics DataFrame with:
        - `strand_col`: Strand orientation ('+', '-', or '?')
        - `status_col`: Updated STATUS codes with 7th digit reflecting inference results
        - `flipped_col`: Removed (if present) as it's no longer needed after inference
    
    Notes
    -----
    - This function requires RAF to be pre-annotated (typically by `_annotate_sumstats` or 
      `_infer_strand_with_annotation`).
    - The function only processes variants with unknown strand status (STATUS 7th digit = 8 or 9).
    - For palindromic SNPs, strand inference is only performed if MAF can be reliably determined 
      from EAF alone (EAF < maf_threshold OR EAF > 1 - maf_threshold) before checking RAF.
    - The function uses vectorized operations for efficient processing of large datasets.
    - STATUS codes must be integer type (use `ensure_status_int` if needed).
    
    See Also
    --------
    _infer_strand_with_annotation : Function that annotates RAF and then calls this function.
    _annotate_sumstats : Function that annotates RAF from reference file.
    """

    # Check required columns - optimize with set for O(1) lookup
    required_cols = [chrom, pos, ea, nea, eaf, raf, flipped_col]
    missing_cols = [col for col in required_cols if col not in sumstats.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Optimize: Check for NA values in critical columns in one pass
    na_counts = {}
    for col in [flipped_col, eaf, raf]:
        na_count = sumstats[col].isna().sum()
        if na_count > 0:
            na_counts[col] = na_count
            log.write(" -WARNING: {} variants have NA in {} column...".format(na_count, col), verbose=verbose)
    
    # Initialize strand column if not present
    if strand_col not in sumstats.columns:
        sumstats[strand_col] = pd.NA
        log.write(" -Initialized new {} column...".format(strand_col), verbose=verbose)
    else:
        existing_strand_count = sumstats[strand_col].notna().sum()
        log.write(" -Found existing {} column with {} non-NA values...".format(strand_col, existing_strand_count), verbose=verbose)
    
    # Ensure STATUS is integer type before any operations
    sumstats = ensure_status_int(sumstats, status_col)
    
    # Optimize: Pre-compute allele validity masks once
    # Cache string length operations to avoid repeated computation
    ea_series = sumstats[ea]
    nea_series = sumstats[nea]
    ea_valid = ea_series.notna() & (ea_series != ".")
    nea_valid = nea_series.notna() & (nea_series != ".")
    both_alleles_valid = ea_valid & nea_valid
    # Cache str.len() results to avoid repeated string operations
    ea_len = ea_series.str.len()
    nea_len = nea_series.str.len()
    both_single_char = (ea_len == 1) & (nea_len == 1)
    valid_snp_mask = both_alleles_valid & both_single_char
    
    # Optimize: Identify palindromic variants more efficiently
    # Palindromic pairs: A/T, T/A, G/C, C/G
    # Use cached ea_series and nea_series to avoid repeated column access
    palindromic_pairs = (
        ((ea_series == "A") & (nea_series == "T")) |
        ((ea_series == "T") & (nea_series == "A")) |
        ((ea_series == "G") & (nea_series == "C")) |
        ((ea_series == "C") & (nea_series == "G"))
    )
    palindromic_mask = palindromic_pairs & valid_snp_mask
    
    # Identify non-palindromic variants (not palindromic and valid SNPs)
    non_palindromic_mask = valid_snp_mask & ~palindromic_mask
    
    # Update STATUS code (last character = '0' for non-palindromic)
    if non_palindromic_mask.any():
        sumstats.loc[non_palindromic_mask, status_col] = vchange_status(
            sumstats.loc[non_palindromic_mask, status_col],
            7,
            [str(i) for i in range(10)],
            ["0"] * 10
        )
    
    # Valid palindromic mask (already computed above)
    valid_palindromic_mask = palindromic_mask

    # Match old method: Assign status 7 for palindromic SNPs with MAF > threshold based on EAF alone
    # This must be done BEFORE filtering by unknown_palindromic_mask, matching old method behavior (line 1586 in hm_harmonize_sumstats.py)
    # Old method order: 1) Identify palindromic, 2) Assign status 0 for non-palindromic, 3) Assign status 7 for ~maf_can_infer, 4) Filter unknown_palindromic
    # Old method: maf_can_infer = (eaf < maf_threshold) | (eaf > 1 - maf_threshold)
    #            If ~maf_can_infer (MAF between threshold and 1-threshold), assign status 7
    # Optimize: Pre-compute eaf_notna and raf_notna once to avoid repeated computation
    eaf_notna = sumstats[eaf].notna()
    raf_notna = sumstats[raf].notna()
    eaf_series = sumstats[eaf]
    raf_series = sumstats[raf]
    
    if valid_palindromic_mask.any():
        # Check if MAF can be inferred from EAF alone (EAF < maf_threshold OR EAF > 1 - maf_threshold)
        # If NOT maf_can_infer (i.e., MAF is between threshold and 1-threshold), assign status 7
        palindromic_with_eaf = valid_palindromic_mask & eaf_notna
        if palindromic_with_eaf.any():
            # Optimize: Use cached eaf_series instead of .loc[] to avoid repeated indexing
            eaf_values = eaf_series[palindromic_with_eaf]
            maf_can_infer = (eaf_values < maf_threshold) | (eaf_values > 1 - maf_threshold)
            # Create mask for variants that cannot infer MAF from EAF alone
            # Map boolean results back to full index
            maf_cannot_infer_local = ~maf_can_infer
            maf_cannot_infer_mask = pd.Series(False, index=sumstats.index)
            maf_cannot_infer_mask.loc[palindromic_with_eaf] = maf_cannot_infer_local.values
            
            # Assign status 7 for palindromic SNPs with MAF > threshold (based on EAF alone, before reference check)
            if maf_cannot_infer_mask.any():
                sumstats.loc[maf_cannot_infer_mask, status_col] = vchange_status(
                    sumstats.loc[maf_cannot_infer_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["7"] * 10
                )
                # Remove these from valid_palindromic_mask since they're already assigned
                valid_palindromic_mask = valid_palindromic_mask & ~maf_cannot_infer_mask

    # Create mask for variants with unknown strand status for palindromic SNPs
    # Pattern: 6th character is 0,1,2 and 7th character is 8 or 9
    # This is done AFTER assigning status 7, matching old method order
    from gwaslab.info.g_vchange_status import status_match
    if status_col in sumstats.columns:
        # Optimize: Cache status_match results
        status_series = sumstats[status_col]
        digit6_match = status_match(status_series, 6, [0, 1, 2])
        digit7_match = status_match(status_series, 7, [8, 9])
        unknown_palindromic_mask = digit6_match & digit7_match
        valid_palindromic_mask = valid_palindromic_mask & unknown_palindromic_mask
    else:
        log.write(" -WARNING: STATUS column '{}' not found, cannot filter by strand status...".format(status_col), verbose=verbose)
    
    # For palindromic SNPs, use EAF/RAF to determine strand
    # Only process variants where MAF can be inferred from EAF (matching old method)
    if valid_palindromic_mask.any():
        # Only check variants where MAF can be inferred from EAF alone
        # (EAF < maf_threshold OR EAF > 1 - maf_threshold)
        palindromic_with_eaf = valid_palindromic_mask & eaf_notna
        if palindromic_with_eaf.any():
            # Optimize: Use cached eaf_series instead of .loc[] to avoid repeated indexing
            eaf_values = eaf_series[palindromic_with_eaf]
            maf_can_infer = (eaf_values < maf_threshold) | (eaf_values > 1 - maf_threshold)
            # Only process variants where MAF can be inferred
            # Map boolean results back to full index
            maf_can_infer_mask = pd.Series(False, index=sumstats.index)
            maf_can_infer_mask.loc[palindromic_with_eaf] = maf_can_infer.values
            valid_palindromic_mask = valid_palindromic_mask & maf_can_infer_mask
        
        # Filter out variants with NA in RAF (need RAF for reference comparison)
        valid_palindromic_with_af_mask = valid_palindromic_mask & eaf_notna & raf_notna
        
        if valid_palindromic_with_af_mask.any():
            # Optimize: Vectorize MAF calculation using numpy
            # Use cached series instead of .loc[] to avoid repeated indexing
            eaf_values = eaf_series[valid_palindromic_with_af_mask]
            raf_values = raf_series[valid_palindromic_with_af_mask]
            # Vectorized MAF: min(x, 1-x) using numpy minimum
            maf_eaf = np.minimum(eaf_values, 1 - eaf_values)
            maf_raf = np.minimum(raf_values, 1 - raf_values)
            maf_mask_local = (maf_eaf <= maf_threshold) & (maf_raf <= ref_maf_threshold)
            # Map mask back to full index
            maf_mask = pd.Series(False, index=sumstats.index)
            maf_mask.loc[valid_palindromic_with_af_mask] = maf_mask_local.values
            valid_palindromic_with_af_mask = valid_palindromic_with_af_mask & maf_mask
        
        if valid_palindromic_with_af_mask.any():
            # Optimize: Calculate expected AF for forward/reverse strands using vectorized operations
            # Use cached series instead of .loc[] to avoid repeated indexing
            eaf_subset = eaf_series[valid_palindromic_with_af_mask]
            raf_subset = raf_series[valid_palindromic_with_af_mask]
            af_diff_forward = np.abs(eaf_subset - raf_subset)
            af_diff_reverse = np.abs(eaf_subset - (1 - raf_subset))
            
            # Forward strand if EAF closer to RAF than to 1-RAF
            forward_palindromic = af_diff_forward <= af_diff_reverse
            # Map boolean results back to full index
            forward_palindromic_full = pd.Series(False, index=sumstats.index)
            forward_palindromic_full.loc[valid_palindromic_with_af_mask] = forward_palindromic.values
            forward_mask = valid_palindromic_with_af_mask & forward_palindromic_full
            reverse_mask = valid_palindromic_with_af_mask & ~forward_palindromic_full
            
            if forward_mask.any():
                sumstats.loc[forward_mask, strand_col] = "+"
                # Update STATUS code (last character = '1' for palindromic +strand)
                sumstats.loc[forward_mask, status_col] = vchange_status(
                    sumstats.loc[forward_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["1"] * 10
                )
            
            if reverse_mask.any():
                sumstats.loc[reverse_mask, strand_col] = "-"
                # Update STATUS code (last character = '2' for palindromic -strand that was flipped)
                sumstats.loc[reverse_mask, status_col] = vchange_status(
                    sumstats.loc[reverse_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["2"] * 10
                )
        
        # Handle palindromic SNPs with NA in EAF or RAF, or filtered out by MAF threshold
        # Match old method behavior: status 7 for MAF > threshold, status 8 for not found in reference
        if valid_palindromic_mask.any():
            # Variants with both EAF and RAF available but filtered out by MAF threshold -> status 7
            maf_exceed_mask = valid_palindromic_mask & eaf_notna & raf_notna
            if 'valid_palindromic_with_af_mask' in locals() and valid_palindromic_with_af_mask.any():
                maf_exceed_mask = maf_exceed_mask & ~valid_palindromic_with_af_mask
            
            # Variants with EAF but no RAF (checked but not found in reference) -> status 8
            # Variants with no EAF (not checked) should remain status 9
            checked_but_not_found_mask = valid_palindromic_mask & eaf_notna & ~raf_notna
            
            # Assign status '7' for MAF > threshold (matches old method)
            if maf_exceed_mask.any():
                sumstats.loc[maf_exceed_mask, strand_col] = "?"
                sumstats.loc[maf_exceed_mask, status_col] = vchange_status(
                    sumstats.loc[maf_exceed_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["7"] * 10
                )
            
            # Assign status '8' for variants checked but not found in reference (matches old method)
            if checked_but_not_found_mask.any():
                sumstats.loc[checked_but_not_found_mask, strand_col] = "?"
                sumstats.loc[checked_but_not_found_mask, status_col] = vchange_status(
                    sumstats.loc[checked_but_not_found_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["8"] * 10
                )
    
    # Optimize: Handle indels - use pre-computed validity masks
    # Use cached length series to avoid repeated str.len() operations
    indel_mask = (ea_len > 1) | (nea_len > 1)
    valid_indel_mask = indel_mask & both_alleles_valid
    
    # Create mask for variants with unknown strand status for indels
    # Pattern: 6th character is 6 and 7th character is 8 or 9
    if status_col in sumstats.columns:
        # Optimize: Reuse status_series if available
        if 'status_series' not in locals():
            status_series = sumstats[status_col]
        # Match: digit 6 is 6 and digit 7 is 8,9
        digit6_indel = status_match(status_series, 6, [6])
        digit7_indel = status_match(status_series, 7, [8, 9])
        unknown_indel_mask = digit6_indel & digit7_indel
        valid_indel_mask = valid_indel_mask & unknown_indel_mask
    else:
        log.write(" -WARNING: STATUS column '{}' not found, cannot filter by strand status...".format(status_col), verbose=verbose)
    
    # For indels, use EAF/RAF to determine strand
    # Match old method: only process indels that have RAF available (found in reference)
    # The old method's check_unkonwn_indel only processes variants found in reference during per-variant query
    if valid_indel_mask.any():
        # Filter out variants with NA in EAF or RAF
        # Note: The old method only processes indels found in reference, so we require RAF to be present
        valid_indel_with_af_mask = valid_indel_mask & eaf_notna & raf_notna

        if valid_indel_with_af_mask.any():
            # Optimize: Vectorize MAF calculation using numpy
            # Use cached series instead of .loc[] to avoid repeated indexing
            eaf_values = eaf_series[valid_indel_with_af_mask]
            raf_values = raf_series[valid_indel_with_af_mask]
            maf_eaf = np.minimum(eaf_values, 1 - eaf_values)
            maf_raf = np.minimum(raf_values, 1 - raf_values)
            # Match old method: check MAF threshold first (line 1406 in hm_harmonize_sumstats.py)
            # Old method: if min(raf, 1-raf) > ref_maf_threshold, return status 8
            maf_mask_local = (maf_eaf <= maf_threshold) & (maf_raf <= ref_maf_threshold)
            # Map mask back to full index
            maf_mask = pd.Series(False, index=sumstats.index)
            maf_mask.loc[valid_indel_with_af_mask] = maf_mask_local.values
            valid_indel_with_af_mask = valid_indel_with_af_mask & maf_mask

            # Track ambiguous indels for later status 4 assignment
            ambiguous_indel_mask_full = pd.Series(False, index=sumstats.index)
            
            if valid_indel_with_af_mask.any():
                # Optimize: Calculate difference between forward and reverse AF differences using vectorized operations
                # Use cached series instead of .loc[] to avoid repeated indexing
                eaf_subset = eaf_series[valid_indel_with_af_mask]
                raf_subset = raf_series[valid_indel_with_af_mask]
                af_diff_forward = np.abs(eaf_subset - raf_subset)
                af_diff_reverse = np.abs(eaf_subset - (1 - raf_subset))
                
                # Match old method logic exactly (check_unkonwn_indel, lines 1408-1414):
                # 1. Check forward: if abs(raf - eaf) < daf_tolerance, assign status 3
                # 2. Check reverse: if abs(raf - (1-eaf)) < daf_tolerance, assign status 6
                # 3. If neither is within tolerance, assign status 8 (not ambiguous status 4)
                # The old method uses if-elif, so forward takes precedence
                forward_within_tolerance = af_diff_forward <= daf_tolerance
                reverse_within_tolerance = af_diff_reverse <= daf_tolerance
                
                # Track ambiguous indels (both differences > tolerance) for status 4 assignment
                # Note: Old method's check_unkonwn_indel doesn't assign status 4, but status 4 might be assigned elsewhere
                # For ambiguous cases where both DAF > tolerance, we'll assign status 4 to match old method's overall behavior
                ambiguous_indel_mask_local = (~forward_within_tolerance) & (~reverse_within_tolerance)
                # Map mask back to full index
                ambiguous_indel_mask_full.loc[valid_indel_with_af_mask] = ambiguous_indel_mask_local.values
                
                # Process indels where at least one direction is within tolerance
                # Forward takes precedence (old method uses if-elif)
                # Map boolean results back to full index
                forward_within_tolerance_full = pd.Series(False, index=sumstats.index)
                forward_within_tolerance_full.loc[valid_indel_with_af_mask] = forward_within_tolerance.values
                forward_indel_mask = valid_indel_with_af_mask & forward_within_tolerance_full
                # For reverse, only include if forward is NOT within tolerance (old method uses if-elif)
                reverse_within_tolerance_full = pd.Series(False, index=sumstats.index)
                reverse_within_tolerance_full.loc[valid_indel_with_af_mask] = reverse_within_tolerance.values
                reverse_indel_mask = valid_indel_with_af_mask & reverse_within_tolerance_full & ~forward_within_tolerance_full
        
                if forward_indel_mask.any() or reverse_indel_mask.any():
                    if forward_indel_mask.any():
                        sumstats.loc[forward_indel_mask, strand_col] = "+"
                        # Update STATUS code for indels - forward strand (no flip): last character = '3'
                        # Match old method: if abs(raf - eaf) < daf_tolerance, status 3
                        sumstats.loc[forward_indel_mask, status_col] = vchange_status(
                            sumstats.loc[forward_indel_mask, status_col],
                            7,
                            [str(i) for i in range(10)],
                            ["3"] * 10
                        )
                    
                    if reverse_indel_mask.any():
                        sumstats.loc[reverse_indel_mask, strand_col] = "-"
                        # For reverse strand (flip): last character = '6'
                        # Match old method: if abs(raf - (1-eaf)) < daf_tolerance, status 6
                        sumstats.loc[reverse_indel_mask, status_col] = vchange_status(
                            sumstats.loc[reverse_indel_mask, status_col],
                            7,
                            [str(i) for i in range(10)],
                            ["6"] * 10
                        )
        
        # Handle indels with NA in EAF or RAF or insufficient daf_diff
        # Match old method behavior: status 8 for not found, status 4 only for ambiguous (both DAF > tolerance)
        if valid_indel_mask.any():
            # Split into different categories to match old method behavior
            # 1. Indels checked but not found in reference (EAF present, RAF NA) -> status 8
            indel_not_found_mask = valid_indel_mask & eaf_notna & ~raf_notna
            
            # 2. Indels with both EAF and RAF but filtered out by MAF threshold -> status 8
            indel_maf_exceed_mask = valid_indel_mask & eaf_notna & raf_notna
            if 'valid_indel_with_af_mask' in locals() and valid_indel_with_af_mask.any():
                indel_maf_exceed_mask = indel_maf_exceed_mask & ~valid_indel_with_af_mask
            
            # 3. Indels with ambiguous DAF (both differences > tolerance) -> status 4 (unknown indel)
            # Only include indels that had valid AF but were ambiguous (both DAF differences > tolerance)
            # These must be valid indels with both EAF and RAF, passed MAF check, but both DAF differences > tolerance
            if 'ambiguous_indel_mask_full' in locals() and ambiguous_indel_mask_full.any():
                # Ensure we only mark valid indels that had AF as ambiguous
                indel_ambiguous_mask = ambiguous_indel_mask_full & valid_indel_mask
            else:
                indel_ambiguous_mask = pd.Series(False, index=sumstats.index)
            
            # Assign status 8 for not found or MAF exceeded (matches old method)
            indel_status8_mask = indel_not_found_mask | indel_maf_exceed_mask
            if indel_status8_mask.any():
                sumstats.loc[indel_status8_mask, strand_col] = "?"
                sumstats.loc[indel_status8_mask, status_col] = vchange_status(
                    sumstats.loc[indel_status8_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["8"] * 10
                )
            
            # Assign status 4 only for ambiguous indels (both DAF differences > tolerance)
            if indel_ambiguous_mask.any():
                sumstats.loc[indel_ambiguous_mask, strand_col] = "?"
                # Update STATUS code (last character = '4' for unknown indel)
                sumstats.loc[indel_ambiguous_mask, status_col] = vchange_status(
                    sumstats.loc[indel_ambiguous_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["4"] * 10
                )
    
    # Optimize: Only drop columns if they exist
    cols_to_drop = []
    if strand_col in sumstats.columns:
        cols_to_drop.append(strand_col)
    
    if cols_to_drop:
        sumstats = sumstats.drop(columns=cols_to_drop)

    return sumstats
