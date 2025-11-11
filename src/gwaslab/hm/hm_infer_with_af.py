import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.g_vchange_status import vchange_status
from .hm_assign_rsid import _annotate_sumstats

def _infer_strand_with_annotation(
    sumstats: pd.DataFrame,
    vcf_path: str | None = None,
    tsv_path: str = "gwaslab_sumstats_lookup_table.txt.gz",
    assign_cols: tuple = ("AF",),
    chr_dict: dict | None = None,
    threads: int = 6,
    reuse_lookup: bool = True,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    overwrite: bool = False,
    is_id_rsid: bool = True,
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    eaf: str = "EAF",
    raf: str = "AF",
    flipped_col: str = "ALLELE_FLIPPED",
    strand_col: str = "STRAND",
    log: "Log" = Log(),
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Annotate summary statistics with reference allele frequency and infer strand orientation.
    
    This function first annotates the summary statistics with reference allele frequency
    using _annotate_sumstats, then infers strand orientation using infer_strand.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        GWAS summary statistics DataFrame.
    vcf_path : str or None
        Path to VCF file for annotation.
    tsv_path : str
        Path to TSV lookup table.
    assign_cols : tuple
        Columns to assign during annotation (default: ("AF",)).
    chr_dict : dict or None
        Chromosome dictionary for reference conversion.
    threads : int
        Number of threads for processing.
    reuse_lookup : bool
        Whether to reuse existing lookup table.
    convert_to_bcf : bool
        Whether to convert VCF to BCF.
    strip_info : bool
        Whether to strip INFO fields.
    overwrite : bool
        Whether to overwrite existing annotations.
    is_id_rsid : bool
        Whether ID is rsID.
    chrom, pos, ea, nea : str
        Column names for genomic coordinates and alleles.
    eaf : str
        Column name for effect allele frequency in sumstats.
    flipped_col : str
        Column name indicating if alleles were flipped.
    strand_col : str
        Column name to store strand information.
    log : Log
        Logging object.
    verbose : bool
        Verbosity flag.
    
    Returns
    -------
    pd.DataFrame
        Annotated and strand-inferred summary statistics.
    """
    # First annotate with AF
    annotated_sumstats = _annotate_sumstats(
        sumstats=sumstats,
        vcf_path=vcf_path,
        tsv_path=tsv_path,
        assign_cols=assign_cols,
        chr_dict=chr_dict,
        threads=threads,
        chrom=chrom,
        pos=pos,
        ea=ea,
        nea=nea,
        overwrite=overwrite,
        is_id_rsid=is_id_rsid,
        reuse_lookup=reuse_lookup,
        convert_to_bcf=convert_to_bcf,
        strip_info=strip_info,
        verbose=verbose,
        log=log,
    )
    
    # Then infer strand using the annotated AF as RAF
    return _infer_strand(
        sumstats=annotated_sumstats,
        chrom=chrom,
        pos=pos,
        ea=ea,
        nea=nea,
        eaf=eaf,
        raf=raf,  # The AF column added by annotation
        flipped_col=flipped_col,
        strand_col=strand_col,
        log=log,
        verbose=verbose,
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
    daf_tolerance: float = 0.20,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Infer strand orientation using FLIPPED column and allele frequencies.
    
    This function determines strand orientation using the FLIPPED column which indicates
    whether alleles were flipped during harmonization, and uses EAF and RAF for additional validation.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        GWAS summary statistics DataFrame containing at least CHR, POS, EA, NEA, EAF, RAF, and FLIPPED columns.
    chrom, pos, ea, nea : str
        Column names for chromosome, position, effect allele, and non-effect allele.
    eaf : str
        Column name for effect allele frequency (from sumstats).
    raf : str
        Column name for reference allele frequency (from reference file).
    flipped_col : str
        Column name indicating if alleles were flipped (True = reverse strand).
    strand_col : str
        Name of the column to store strand information ('+' for forward, '-' for reverse).
    status_col : str
        Name of the column containing status codes (default: "STATUS").
    log : Log
        Logging object for progress messages.
    maf_threshold : float
        Maximum minor allele frequency threshold for palindromic SNPs (default: 0.40).
        Palindromic SNPs with MAF > threshold will be excluded from strand determination.
    daf_tolerance : float
        Difference in allele frequency tolerance for indels (default: 0.20).
        Indels with |daf_forward - daf_reverse| < tolerance will be marked as ambiguous.
    verbose : bool
        Whether to log progress messages.
    
    Returns
    -------
    pd.DataFrame
        Updated summary statistics with strand information and updated STATUS codes.
    """

    # Check required columns
    required_cols = [chrom, pos, ea, nea, eaf, raf, flipped_col]
    for col in required_cols:
        if col not in sumstats.columns:
            raise ValueError(f"Missing required column: {col}")
    
    # Check for NA values in critical columns
    flipped_na_count = sumstats[flipped_col].isna().sum()
    if flipped_na_count > 0:
        log.write(f"[infer_strand] WARNING: {flipped_na_count} variants have NA in {flipped_col} column", verbose=verbose)
    
    eaf_na_count = sumstats[eaf].isna().sum()
    if eaf_na_count > 0:
        log.write(f"[infer_strand] WARNING: {eaf_na_count} variants have NA in {eaf} column", verbose=verbose)
    
    raf_na_count = sumstats[raf].isna().sum()
    if raf_na_count > 0:
        log.write(f"[infer_strand] WARNING: {raf_na_count} variants have NA in {raf} column", verbose=verbose)
    
    # Initialize strand column if not present
    if strand_col not in sumstats.columns:
        sumstats[strand_col] = pd.NA
        log.write(f"[infer_strand] Initialized new {strand_col} column", verbose=verbose)
    else:
        existing_strand_count = sumstats[strand_col].notna().sum()
        log.write(f"[infer_strand] Found existing {strand_col} column with {existing_strand_count} non-NA values", verbose=verbose)
    
    # Identify non-palindromic variants (not A/T, T/A, G/C, C/G) and not indels
    non_palindromic_mask = ~(
        ((sumstats[ea] == "A") & (sumstats[nea] == "T")) |
        ((sumstats[ea] == "T") & (sumstats[nea] == "A")) |
        ((sumstats[ea] == "G") & (sumstats[nea] == "C")) |
        ((sumstats[ea] == "C") & (sumstats[nea] == "G"))
    ) & sumstats[ea].notna() & sumstats[nea].notna() & (sumstats[ea].str.len() == 1) & (sumstats[nea].str.len() == 1)
    
    # Update STATUS code (last character = '0' for non-palindromic)
    sumstats.loc[non_palindromic_mask, status_col] = vchange_status(
        sumstats.loc[non_palindromic_mask, status_col],
        7,
        [str(i) for i in range(10)],
        ["0"] * 10
    )
    
    # Handle palindromic SNPs
    palindromic_mask = (
        ((sumstats[ea] == "A") & (sumstats[nea] == "T")) |
        ((sumstats[ea] == "T") & (sumstats[nea] == "A")) |
        ((sumstats[ea] == "G") & (sumstats[nea] == "C")) |
        ((sumstats[ea] == "C") & (sumstats[nea] == "G"))
    )
    
    # Handle NA values in allele columns for palindromic detection
    palindromic_na_mask = (
        sumstats[ea].isna() | 
        sumstats[nea].isna() |
        (sumstats[ea] == ".") |
        (sumstats[nea] == ".")
    )
    valid_palindromic_mask = palindromic_mask & ~palindromic_na_mask
    
    # Create mask for variants with unknown strand status for palindromic SNPs
    # Pattern: 6th character is 0,1,2 and 7th character is 8 or 9
    if status_col in sumstats.columns:
        unknown_palindromic_mask = sumstats[status_col].str.match(r'\w\w\w\w\w[012][89]', case=False, na=False)
        valid_palindromic_mask = valid_palindromic_mask & unknown_palindromic_mask
    else:
        log.write(f"[infer_strand] WARNING: STATUS column '{status_col}' not found, cannot filter by strand status", verbose=verbose)
    
    # For palindromic SNPs, use EAF/RAF to determine strand
    if valid_palindromic_mask.any():
        # Filter out variants with NA in EAF or RAF
        valid_palindromic_with_af_mask = valid_palindromic_mask & sumstats[eaf].notna() & sumstats[raf].notna()
        
        # Apply MAF threshold filter for palindromic SNPs
        eaf_values = sumstats.loc[valid_palindromic_with_af_mask, eaf]
        raf_values = sumstats.loc[valid_palindromic_with_af_mask, raf]
        maf_eaf = eaf_values.apply(lambda x: min(x, 1-x))
        maf_raf = raf_values.apply(lambda x: min(x, 1-x))
        maf_mask = (maf_eaf <= maf_threshold) & (maf_raf <= maf_threshold)
        valid_palindromic_with_af_mask = valid_palindromic_with_af_mask & maf_mask
        
        # Calculate expected AF for forward/reverse strands
        af_diff_forward = abs(sumstats.loc[valid_palindromic_with_af_mask, eaf] - sumstats.loc[valid_palindromic_with_af_mask, raf])
        af_diff_reverse = abs(sumstats.loc[valid_palindromic_with_af_mask, eaf] - (1 - sumstats.loc[valid_palindromic_with_af_mask, raf]))
        
        # Forward strand if EAF closer to RAF than to 1-RAF
        forward_palindromic = af_diff_forward <= af_diff_reverse
        sumstats.loc[valid_palindromic_with_af_mask & forward_palindromic, strand_col] = "+"
        
        # Update STATUS code (last character = '1' for palindromic +strand)
        sumstats.loc[valid_palindromic_with_af_mask & forward_palindromic, status_col] = vchange_status(
            sumstats.loc[valid_palindromic_with_af_mask & forward_palindromic, status_col],
            7,
            [str(i) for i in range(10)],
            ["1"] * 10
        )
        
        # Reverse strand otherwise
        reverse_palindromic = ~forward_palindromic
        sumstats.loc[valid_palindromic_with_af_mask & reverse_palindromic, strand_col] = "-"
        # Update STATUS code (last character = '2' for palindromic -strand that was flipped)
        sumstats.loc[valid_palindromic_with_af_mask & reverse_palindromic, status_col] = vchange_status(
            sumstats.loc[valid_palindromic_with_af_mask & reverse_palindromic, status_col],
            7,
            [str(i) for i in range(10)],
            ["2"] * 10
        )
        
        # Handle palindromic SNPs with NA in EAF or RAF, or filtered out by MAF threshold
        # Split into MAF > threshold and NA in EAF/RAF
        maf_exceed_mask = valid_palindromic_mask & sumstats[eaf].notna() & sumstats[raf].notna() & ~maf_mask
        na_raf_mask = valid_palindromic_mask & (sumstats[eaf].isna() | sumstats[raf].isna())
        
        # Assign status '7' for MAF > threshold
        if maf_exceed_mask.any():
            sumstats.loc[maf_exceed_mask, strand_col] = "?"
            sumstats.loc[maf_exceed_mask, status_col] = vchange_status(
                sumstats.loc[maf_exceed_mask, status_col],
                7,
                [str(i) for i in range(10)],
                ["7"] * 10
            )
        
        # Assign status '8' for NA in EAF/RAF
        if na_raf_mask.any():
            sumstats.loc[na_raf_mask, strand_col] = "?"
            sumstats.loc[na_raf_mask, status_col] = vchange_status(
                sumstats.loc[na_raf_mask, status_col],
                7,
                [str(i) for i in range(10)],
                ["8"] * 10
            )
    
    # Handle indels
    indel_mask = (
        (sumstats[ea].str.len() > 1) | 
        (sumstats[nea].str.len() > 1)
    )
    
    # Handle NA values in allele columns for indel detection
    indel_na_mask = (
        sumstats[ea].isna() | 
        sumstats[nea].isna() |
        (sumstats[ea] == ".") |
        (sumstats[nea] == ".")
    )
    valid_indel_mask = indel_mask & ~indel_na_mask
    
    # Create mask for variants with unknown strand status for indels
    # Pattern: 6th character is 6 and 7th character is 8 or 9
    if status_col in sumstats.columns:
        unknown_indel_mask = sumstats[status_col].str.match(r'\w\w\w\w\w[6][89]', case=False, na=False)
        valid_indel_mask = valid_indel_mask & unknown_indel_mask
    else:
        log.write(f"[infer_strand] WARNING: STATUS column '{status_col}' not found, cannot filter by strand status", verbose=verbose)
    
    # For indels, use EAF/RAF to determine strand
    if valid_indel_mask.any():
        # Filter out variants with NA in EAF or RAF
        valid_indel_with_af_mask = valid_indel_mask & sumstats[eaf].notna() & sumstats[raf].notna()
            
        # Calculate difference between forward and reverse AF differences
        af_diff_forward = abs(sumstats.loc[valid_indel_with_af_mask, eaf] - sumstats.loc[valid_indel_with_af_mask, raf])
        af_diff_reverse = abs(sumstats.loc[valid_indel_with_af_mask, eaf] - (1 - sumstats.loc[valid_indel_with_af_mask, raf]))
        
        # Mark as ambiguous if both differences exceed tolerance
        ambiguous_indel_mask = (af_diff_forward > daf_tolerance) & (af_diff_reverse > daf_tolerance)
        valid_indel_with_af_mask = valid_indel_with_af_mask & ~ambiguous_indel_mask
        
        # Only process indels with sufficient daf_diff
        if valid_indel_with_af_mask.any():
            forward_indel = af_diff_forward <= af_diff_reverse
            sumstats.loc[valid_indel_with_af_mask & forward_indel, strand_col] = "+"
            sumstats.loc[valid_indel_with_af_mask & ~forward_indel, strand_col] = "-"
            
            # Update STATUS code for indels
            # For forward strand (no flip): last character = '3'
            sumstats.loc[valid_indel_with_af_mask & forward_indel, status_col] = vchange_status(
                sumstats.loc[valid_indel_with_af_mask & forward_indel, status_col],
                7,
                [str(i) for i in range(10)],
                ["3"] * 10
            )
            
            # For reverse strand (flip): last character = '6'
            sumstats.loc[valid_indel_with_af_mask & ~forward_indel, status_col] = vchange_status(
                sumstats.loc[valid_indel_with_af_mask & ~forward_indel, status_col],
                7,
                [str(i) for i in range(10)],
                ["6"] * 10
            )
        
        # Handle indels with NA in EAF or RAF or insufficient daf_diff
        indel_ambiguous_mask = valid_indel_mask & ~valid_indel_with_af_mask
        if indel_ambiguous_mask.any():
            sumstats.loc[indel_ambiguous_mask, strand_col] = "?"
            # Update STATUS code (last character = '4' for unknown indel, fixed)
            sumstats.loc[indel_ambiguous_mask, status_col] = vchange_status(
                sumstats.loc[indel_ambiguous_mask, status_col],
                7,
                [str(i) for i in range(10)],
                ["4"] * 10
            )
    
    # Handle ambiguous cases (NA in FLIPPED column)
    ambiguous_mask = sumstats[raf].isna() & valid_indel_with_af_mask & valid_palindromic_with_af_mask
    if ambiguous_mask.any():
        sumstats.loc[ambiguous_mask, strand_col] = "?"
        # Update STATUS code (last character = '8' for no ref data)
        sumstats.loc[ambiguous_mask, status_col] = vchange_status(
            sumstats.loc[ambiguous_mask, status_col],
            7,
            [str(i) for i in range(10)],
            ["8"] * 10
        )
    
    return sumstats
