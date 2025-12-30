import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.info.g_vchange_status import vchange_status
from gwaslab.info.g_vchange_status import ensure_status_int
from .hm_assign_rsid import _annotate_sumstats
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

# Epsilon for floating point precision in frequency comparisons
# Used to handle floating point precision issues when comparing MAF values at threshold
FREQ_COMPARISON_EPSILON = 1e-6

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
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided and sumstats is a Sumstats object, uses sumstats.mapper.
        If not provided, creates a default mapper with automatic format detection.
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
    - The function automatically handles chromosome name mapping using ChromosomeMapper.
    - Lookup tables are cached as TSV files for faster subsequent runs when `reuse_lookup=True`.
    
    See Also
    --------
    _annotate_sumstats : Function that performs the annotation step.
    _infer_strand : Function that performs the strand inference step.
    """
    # Handle both DataFrame and Sumstats object inputs
    if isinstance(sumstats, pd.DataFrame):
        is_dataframe = True
    else:
        sumstats_obj = sumstats
        sumstats = sumstats_obj.data
        is_dataframe = False

    # Initialize chromosome mapper for coordinate conversion
    if mapper is None:
        if not is_dataframe and hasattr(sumstats_obj, 'mapper'):
            mapper = sumstats_obj.mapper
        else:
            species = sumstats_obj.meta.get("gwaslab", {}).get("species", "homo sapiens") if not is_dataframe else "homo sapiens"
            build = sumstats_obj.build if not is_dataframe and hasattr(sumstats_obj, 'build') else None
            mapper = ChromosomeMapper(species=species, build=build, log=log, verbose=verbose)
            # Auto-detect chromosome format from sumstats data
            if not is_dataframe and hasattr(sumstats_obj, 'data') and not sumstats_obj.data.empty and chrom in sumstats_obj.data.columns:
                mapper.detect_sumstats_format(sumstats_obj.data[chrom])

    # Normalize assign_cols parameter to tuple format
    if isinstance(assign_cols, str):
        assign_cols = tuple([assign_cols])
    elif isinstance(assign_cols, list):
        assign_cols = tuple(assign_cols)
    else:
        assign_cols = tuple(assign_cols)
        
    # If converting to BCF, preserve INFO fields (needed for AF extraction)
    if convert_to_bcf:
        strip_info = False

    # Auto-detect reference file format for chromosome mapping
    if vcf_path is not None:
        mapper.detect_reference_format(vcf_path)
    elif path is not None:
        from gwaslab.hm.hm_assign_rsid import is_vcf_file
        if is_vcf_file(path):
            mapper.detect_reference_format(path)

    # Step 1: Annotate sumstats with reference allele frequency (RAF)
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
    
    # Rename annotated column to RAF (if different from source column name)
    if assign_cols and len(assign_cols) > 0:
        source_col = assign_cols[0]
        if source_col in annotated_sumstats.columns and source_col != raf:
            log.write(" -Renaming {} in reference to {} in Sumstats".format(source_col, raf), verbose=verbose)
            # Drop target column if it exists to avoid conflicts
            if raf in annotated_sumstats.columns:
                annotated_sumstats = annotated_sumstats.drop(columns=[raf])
            annotated_sumstats = annotated_sumstats.rename(columns={source_col: raf})

    # Step 2: Infer strand orientation using EAF and RAF
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
    
    # Remove temporary ALLELE_FLIPPED column (used internally for strand inference)
    if "ALLELE_FLIPPED" in result.columns:
        result = result.drop(columns=["ALLELE_FLIPPED"])
    
    # Update metadata and harmonization status if input was a Sumstats object
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
                'threads': threads, 'reuse_lookup': reuse_lookup,
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

    # Validate required columns are present
    required_cols = [chrom, pos, ea, nea, eaf, raf, flipped_col]
    missing_cols = [col for col in required_cols if col not in sumstats.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Check for missing values in critical columns and log warnings
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
    
    # Pre-compute allele validity and type masks for efficient processing
    # Extract allele series once to avoid repeated column access
    ea_series = sumstats[ea]
    nea_series = sumstats[nea]
    ea_valid = ea_series.notna() & (ea_series != ".")
    nea_valid = nea_series.notna() & (nea_series != ".")
    both_alleles_valid = ea_valid & nea_valid
    
    # Identify SNPs (single-character alleles) vs indels (multi-character)
    ea_len = ea_series.str.len()
    nea_len = nea_series.str.len()
    both_single_char = (ea_len == 1) & (nea_len == 1)
    valid_snp_mask = both_alleles_valid & both_single_char
    
    # Identify palindromic SNPs: A/T, T/A, G/C, C/G
    # These require strand inference since they look the same on both strands
    palindromic_pairs = (
        ((ea_series == "A") & (nea_series == "T")) |
        ((ea_series == "T") & (nea_series == "A")) |
        ((ea_series == "G") & (nea_series == "C")) |
        ((ea_series == "C") & (nea_series == "G"))
    )
    palindromic_mask = palindromic_pairs & valid_snp_mask
    
    # Non-palindromic SNPs: strand is unambiguous, no inference needed
    non_palindromic_mask = valid_snp_mask & ~palindromic_mask
    
    # Set STATUS 7th digit = '0' for non-palindromic SNPs (strand unambiguous)
    if non_palindromic_mask.any():
        sumstats.loc[non_palindromic_mask, status_col] = vchange_status(
            sumstats.loc[non_palindromic_mask, status_col],
            7,
            [str(i) for i in range(10)],
            ["0"] * 10
        )
    
    # Start with all palindromic SNPs for further processing
    valid_palindromic_mask = palindromic_mask

    # Pre-compute frequency series and validity masks for efficient processing
    eaf_notna = sumstats[eaf].notna()
    raf_notna = sumstats[raf].notna()
    eaf_series = sumstats[eaf]
    raf_series = sumstats[raf]
    
    # Step 1: Filter palindromic SNPs by MAF threshold (based on EAF alone)
    # For palindromic SNPs, we can only infer strand if MAF is clearly defined
    # MAF can be inferred if: EAF <= maf_threshold OR EAF >= (1 - maf_threshold)
    # If MAF is between threshold and (1-threshold), strand cannot be reliably inferred
    if valid_palindromic_mask.any():
        palindromic_with_eaf = valid_palindromic_mask & eaf_notna
        if palindromic_with_eaf.any():
            eaf_values = eaf_series[palindromic_with_eaf]
            # Check if MAF can be determined from EAF: either minor or major allele is clear
            maf_can_infer = (eaf_values <= maf_threshold + FREQ_COMPARISON_EPSILON) | (eaf_values >= 1 - maf_threshold - FREQ_COMPARISON_EPSILON)
            
            # Variants where MAF cannot be inferred (ambiguous frequency) -> STATUS 7
            maf_cannot_infer_local = ~maf_can_infer
            maf_cannot_infer_mask = pd.Series(False, index=sumstats.index)
            maf_cannot_infer_mask.loc[palindromic_with_eaf] = maf_cannot_infer_local.values
            
            if maf_cannot_infer_mask.any():
                sumstats.loc[maf_cannot_infer_mask, strand_col] = "?"
                sumstats.loc[maf_cannot_infer_mask, status_col] = vchange_status(
                    sumstats.loc[maf_cannot_infer_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["7"] * 10
                )
                # Exclude these from further processing
                valid_palindromic_mask = valid_palindromic_mask & ~maf_cannot_infer_mask

    # Step 2: Filter to palindromic SNPs with unknown strand status
    # Only process variants where STATUS indicates unknown strand (7th digit = 8 or 9)
    # and variant type is SNP (6th digit = 0, 1, or 2)
    from gwaslab.info.g_vchange_status import status_match
    if status_col in sumstats.columns:
        status_series = sumstats[status_col]
        digit6_match = status_match(status_series, 6, [0, 1, 2])  # SNP types
        digit7_match = status_match(status_series, 7, [8, 9])      # Unknown strand
        unknown_palindromic_mask = digit6_match & digit7_match
        valid_palindromic_mask = valid_palindromic_mask & unknown_palindromic_mask
    else:
        log.write(" -WARNING: STATUS column '{}' not found, cannot filter by strand status...".format(status_col), verbose=verbose)
    
    # Step 3: Infer strand for palindromic SNPs using EAF and RAF comparison
    if valid_palindromic_mask.any():
        # Require both EAF and RAF to be available for comparison
        valid_palindromic_with_af_mask = valid_palindromic_mask & eaf_notna & raf_notna
        
        if valid_palindromic_with_af_mask.any():
            # Filter by MAF threshold: both EAF and RAF must have MAF <= threshold
            # This ensures frequencies are unambiguous enough for reliable inference
            eaf_values = eaf_series[valid_palindromic_with_af_mask]
            raf_values = raf_series[valid_palindromic_with_af_mask]
            maf_eaf = np.minimum(eaf_values, 1 - eaf_values)
            maf_raf = np.minimum(raf_values, 1 - raf_values)
            maf_mask_local = (maf_eaf <= maf_threshold + FREQ_COMPARISON_EPSILON) & (maf_raf <= ref_maf_threshold + FREQ_COMPARISON_EPSILON)
            
            # Map mask back to full index
            maf_mask = pd.Series(False, index=sumstats.index)
            maf_mask.loc[valid_palindromic_with_af_mask] = maf_mask_local.values
            valid_palindromic_with_af_mask = valid_palindromic_with_af_mask & maf_mask
        
        if valid_palindromic_with_af_mask.any():
            # Compare EAF with RAF to determine strand orientation
            # Forward strand: EAF should match RAF (same strand)
            # Reverse strand: EAF should match (1 - RAF) (opposite strand, alleles flipped)
            eaf_subset = eaf_series[valid_palindromic_with_af_mask]
            raf_subset = raf_series[valid_palindromic_with_af_mask]
            af_diff_forward = np.abs(eaf_subset - raf_subset)
            af_diff_reverse = np.abs(eaf_subset - (1 - raf_subset))
            
            # Assign forward strand if EAF is closer to RAF than to (1-RAF)
            forward_palindromic = af_diff_forward <= af_diff_reverse
            
            # Map results back to full index
            forward_palindromic_full = pd.Series(False, index=sumstats.index)
            forward_palindromic_full.loc[valid_palindromic_with_af_mask] = forward_palindromic.values
            forward_mask = valid_palindromic_with_af_mask & forward_palindromic_full
            reverse_mask = valid_palindromic_with_af_mask & ~forward_palindromic_full
            
            # Assign forward strand: STATUS 7th digit = '1'
            if forward_mask.any():
                sumstats.loc[forward_mask, strand_col] = "+"
                sumstats.loc[forward_mask, status_col] = vchange_status(
                    sumstats.loc[forward_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["1"] * 10
                )
            
            # Assign reverse strand: STATUS 7th digit = '2'
            if reverse_mask.any():
                sumstats.loc[reverse_mask, strand_col] = "-"
                sumstats.loc[reverse_mask, status_col] = vchange_status(
                    sumstats.loc[reverse_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["2"] * 10
                )
        
        # Step 4: Handle remaining palindromic SNPs that couldn't be processed
        if valid_palindromic_mask.any():
            # Category 1: Variants with both EAF and RAF but filtered out by MAF threshold
            maf_exceed_mask = valid_palindromic_mask & eaf_notna & raf_notna
            if 'valid_palindromic_with_af_mask' in locals() and valid_palindromic_with_af_mask.any():
                maf_exceed_mask = maf_exceed_mask & ~valid_palindromic_with_af_mask
                
                # Sub-category 1a: MAF(RAF) > ref_maf_threshold -> STATUS 8 (reference MAF too high)
                if maf_exceed_mask.any():
                    raf_exceed_values = raf_series[maf_exceed_mask]
                    maf_raf_exceed = np.minimum(raf_exceed_values, 1 - raf_exceed_values)
                    maf_raf_exceed_bool = maf_raf_exceed > ref_maf_threshold + FREQ_COMPARISON_EPSILON
                    maf_raf_exceed_mask = pd.Series(False, index=sumstats.index)
                    maf_raf_exceed_mask.loc[maf_exceed_mask] = maf_raf_exceed_bool.values
                    
                    if maf_raf_exceed_mask.any():
                        sumstats.loc[maf_raf_exceed_mask, strand_col] = "?"
                        sumstats.loc[maf_raf_exceed_mask, status_col] = vchange_status(
                            sumstats.loc[maf_raf_exceed_mask, status_col],
                            7,
                            [str(i) for i in range(10)],
                            ["8"] * 10
                        )
                        maf_exceed_mask = maf_exceed_mask & ~maf_raf_exceed_mask
            
            # Category 2: Variants checked but not found in reference (EAF present, RAF missing) -> STATUS 8
            checked_but_not_found_mask = valid_palindromic_mask & eaf_notna & ~raf_notna
            
            # Category 3: Variants with MAF(EAF) > threshold but MAF(RAF) <= threshold -> STATUS 7
            # (These passed RAF check but failed EAF check)
            if maf_exceed_mask.any():
                sumstats.loc[maf_exceed_mask, strand_col] = "?"
                sumstats.loc[maf_exceed_mask, status_col] = vchange_status(
                    sumstats.loc[maf_exceed_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["7"] * 10
                )
            
            # Assign STATUS 8 for variants not found in reference
            if checked_but_not_found_mask.any():
                sumstats.loc[checked_but_not_found_mask, strand_col] = "?"
                sumstats.loc[checked_but_not_found_mask, status_col] = vchange_status(
                    sumstats.loc[checked_but_not_found_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["8"] * 10
                )
    
    # Process indels: variants with multi-character alleles
    # Use pre-computed length masks to identify indels
    indel_mask = (ea_len > 1) | (nea_len > 1)
    valid_indel_mask = indel_mask & both_alleles_valid
    
    # Filter to indels with unknown strand status
    # Only process variants where STATUS indicates unknown indel strand (6th digit = 6, 7th digit = 8 or 9)
    if status_col in sumstats.columns:
        if 'status_series' not in locals():
            status_series = sumstats[status_col]
        digit6_indel = status_match(status_series, 6, [6])      # Indel type
        digit7_indel = status_match(status_series, 7, [8, 9])   # Unknown strand
        unknown_indel_mask = digit6_indel & digit7_indel
        valid_indel_mask = valid_indel_mask & unknown_indel_mask
    else:
        log.write(" -WARNING: STATUS column '{}' not found, cannot filter by strand status...".format(status_col), verbose=verbose)
    
    # Infer strand for indels using EAF and RAF comparison with tolerance
    # Only process indels found in reference (RAF available)
    if valid_indel_mask.any():
        # Require both EAF and RAF to be available
        valid_indel_with_af_mask = valid_indel_mask & eaf_notna & raf_notna

        if valid_indel_with_af_mask.any():
            # Filter by MAF threshold: both EAF and RAF must have MAF <= threshold
            eaf_values = eaf_series[valid_indel_with_af_mask]
            raf_values = raf_series[valid_indel_with_af_mask]
            maf_eaf = np.minimum(eaf_values, 1 - eaf_values)
            maf_raf = np.minimum(raf_values, 1 - raf_values)
            maf_mask_local = (maf_eaf <= maf_threshold + FREQ_COMPARISON_EPSILON) & (maf_raf <= ref_maf_threshold + FREQ_COMPARISON_EPSILON)
            
            # Map mask back to full index
            maf_mask = pd.Series(False, index=sumstats.index)
            maf_mask.loc[valid_indel_with_af_mask] = maf_mask_local.values
            valid_indel_with_af_mask = valid_indel_with_af_mask & maf_mask

            # Track ambiguous indels (both forward and reverse differences exceed tolerance)
            ambiguous_indel_mask_full = pd.Series(False, index=sumstats.index)
            
            if valid_indel_with_af_mask.any():
                # Compare EAF with RAF using tolerance-based approach
                # For indels, we use daf_tolerance instead of exact matching
                eaf_subset = eaf_series[valid_indel_with_af_mask]
                raf_subset = raf_series[valid_indel_with_af_mask]
                af_diff_forward = np.abs(eaf_subset - raf_subset)
                af_diff_reverse = np.abs(eaf_subset - (1 - raf_subset))
                
                # Determine strand based on which difference is within tolerance
                # Forward strand: |EAF - RAF| < daf_tolerance (matches old method: abs(raf - eaf) < daf_tolerance)
                # Reverse strand: |EAF - (1-RAF)| < daf_tolerance (matches old method: abs(raf - (1-eaf)) < daf_tolerance)
                # Ambiguous: both differences >= tolerance
                forward_within_tolerance = af_diff_forward < daf_tolerance
                reverse_within_tolerance = af_diff_reverse < daf_tolerance
                
                # Track ambiguous cases (both differences exceed tolerance)
                ambiguous_indel_mask_local = (~forward_within_tolerance) & (~reverse_within_tolerance)
                ambiguous_indel_mask_full.loc[valid_indel_with_af_mask] = ambiguous_indel_mask_local.values
                
                # Map results back to full index
                forward_within_tolerance_full = pd.Series(False, index=sumstats.index)
                forward_within_tolerance_full.loc[valid_indel_with_af_mask] = forward_within_tolerance.values
                forward_indel_mask = valid_indel_with_af_mask & forward_within_tolerance_full
                
                # Reverse strand only if forward is not within tolerance (forward takes precedence)
                reverse_within_tolerance_full = pd.Series(False, index=sumstats.index)
                reverse_within_tolerance_full.loc[valid_indel_with_af_mask] = reverse_within_tolerance.values
                reverse_indel_mask = valid_indel_with_af_mask & reverse_within_tolerance_full & ~forward_within_tolerance_full
        
                # Assign forward strand: STATUS 7th digit = '3'
                if forward_indel_mask.any():
                    sumstats.loc[forward_indel_mask, strand_col] = "+"
                    sumstats.loc[forward_indel_mask, status_col] = vchange_status(
                        sumstats.loc[forward_indel_mask, status_col],
                        7,
                        [str(i) for i in range(10)],
                        ["3"] * 10
                    )
                
                # Assign reverse strand: STATUS 7th digit = '6'
                if reverse_indel_mask.any():
                    sumstats.loc[reverse_indel_mask, strand_col] = "-"
                    sumstats.loc[reverse_indel_mask, status_col] = vchange_status(
                        sumstats.loc[reverse_indel_mask, status_col],
                        7,
                        [str(i) for i in range(10)],
                        ["6"] * 10
                    )
        
        # Handle remaining indels that couldn't be processed
        if valid_indel_mask.any():
            # Category 1: Indels checked but not found in reference (EAF present, RAF missing) -> STATUS 8
            indel_not_found_mask = valid_indel_mask & eaf_notna & ~raf_notna
            
            # Category 2: Indels with both EAF and RAF but filtered out by MAF threshold -> STATUS 8
            indel_maf_exceed_mask = valid_indel_mask & eaf_notna & raf_notna
            if 'valid_indel_with_af_mask' in locals() and valid_indel_with_af_mask.any():
                indel_maf_exceed_mask = indel_maf_exceed_mask & ~valid_indel_with_af_mask
            
            # Category 3: Ambiguous indels (both forward and reverse differences > tolerance) -> STATUS 4
            # These are indels with valid AF that passed MAF check but couldn't be assigned to either strand
            if 'ambiguous_indel_mask_full' in locals() and ambiguous_indel_mask_full.any():
                indel_ambiguous_mask = ambiguous_indel_mask_full & valid_indel_mask
            else:
                indel_ambiguous_mask = pd.Series(False, index=sumstats.index)
            
            # Assign STATUS 8 for not found or MAF exceeded
            indel_status8_mask = indel_not_found_mask | indel_maf_exceed_mask
            if indel_status8_mask.any():
                sumstats.loc[indel_status8_mask, strand_col] = "?"
                sumstats.loc[indel_status8_mask, status_col] = vchange_status(
                    sumstats.loc[indel_status8_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["8"] * 10
                )
            
            # Assign STATUS 4 for ambiguous indels (both DAF differences exceed tolerance)
            if indel_ambiguous_mask.any():
                sumstats.loc[indel_ambiguous_mask, strand_col] = "?"
                sumstats.loc[indel_ambiguous_mask, status_col] = vchange_status(
                    sumstats.loc[indel_ambiguous_mask, status_col],
                    7,
                    [str(i) for i in range(10)],
                    ["4"] * 10
                )
    
    # Note: strand_col is kept in the output for user reference
    # It is not dropped here (unlike ALLELE_FLIPPED which is dropped in _infer_strand_with_annotation)
    
    return sumstats
