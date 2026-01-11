from __future__ import annotations
import numpy as np
import pandas as pd
from math import sqrt
from typing import Optional, Tuple, List, Union
from allel import read_vcf, GenotypeArray
from scipy.special import erfc
import scipy.stats as ss
from gwaslab.g_Sumstats import Sumstats
from gwaslab.io.io_vcf import auto_check_vcf_chr_dict
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log
from shutil import which

# Polars is required for performance
import polars as pl


# =============================================================================
# Small math helpers (NO SciPy dependency)
# =============================================================================

def _norm_sf(x: np.ndarray) -> np.ndarray:
    """
    Survival function of standard normal: sf(x) = P(Z > x), Z~N(0,1)

    sf(x) = 0.5 * erfc(x/sqrt(2))
    Uses scipy.special.erfc for vectorized computation (much faster than np.vectorize).
    
    Handles edge cases:
    - Very large |x| values to prevent underflow
    - NaN/Inf inputs
    """
    # Clip extreme values to prevent numerical issues
    # For |x| > 37, erfc is effectively 0 or 1
    x_clipped = np.clip(x, -37.0, 37.0)
    result = 0.5 * erfc(x_clipped / sqrt(2.0))
    # Ensure result is finite
    result = np.where(np.isfinite(result), result, 0.0)
    return result


def _p_from_z(z: np.ndarray) -> np.ndarray:
    """
    Two-sided p-value from Z:
        P = 2 * sf(|Z|)
    """
    return 2.0 * _norm_sf(np.abs(z))


def _z_to_mlog10p(z: np.ndarray) -> np.ndarray:
    """
    Convert Z-score to -log10(P-value) using log-space for numerical precision.
    Uses the same approach as util_in_fill_data.py for consistency.
    """
    z_arr = np.asarray(z, dtype=np.float64)
    # Two-sided test: log_pvalue = log(2) + norm.logsf(|z|)
    log_pvalue = np.log(2) + ss.norm.logsf(np.abs(z_arr))
    mlog10p = log_pvalue / np.log(10)
    return -mlog10p


def _make_blocks_by_bp(pos: np.ndarray, window_bp: int) -> List[Tuple[int, int]]:
    """
    Partition sorted variant positions into blocks of ~window_bp basepairs.
    Returns half-open index intervals [start, end).
    """
    blocks: List[Tuple[int, int]] = []
    n = len(pos)
    i = 0
    while i < n:
        p0 = pos[i]
        j = i + 1
        while j < n and (pos[j] - p0) < window_bp:
            j += 1
        blocks.append((i, j))
        i = j
    return blocks


# =============================================================================
# VCF reading helpers
# =============================================================================

def _read_vcf_variants(
    vcf_path: str,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    vcf_chr_dict: Optional[dict] = None,
    tabix: Optional[bool] = None,
    mapper: Optional[ChromosomeMapper] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> pl.DataFrame:
    """
    Read variants from VCF file and return Polars DataFrame with CHR, POS, EA, NEA, EAF, INFO.
    
    Parameters
    ----------
    vcf_path : str
        Path to VCF file
    region : tuple, optional
        (chromosome, start, end) region to extract
    vcf_chr_dict : dict, optional
        Chromosome mapping dictionary
    tabix : bool, optional
        Whether to use tabix indexing
    mapper : ChromosomeMapper, optional
        Chromosome mapper instance
    log : Log, optional
        Logging object
    verbose : bool
        Verbose flag
        
    Returns
    -------
    pl.DataFrame
        Polars DataFrame with columns: CHR, POS, EA, NEA, EAF, INFO (if available)
    """
    if log is None:
        log = Log()
    
    if tabix is None:
        tabix = which("tabix")
    
    # Get or create mapper
    if mapper is None:
        mapper = ChromosomeMapper(log=log, verbose=verbose)
        mapper.detect_reference_format(vcf_path)
    
    # Auto-detect chromosome dictionary if not provided
    if vcf_chr_dict is None:
        vcf_chr_dict = auto_check_vcf_chr_dict(vcf_path, None, verbose, log)
    
    # Read VCF data
    # Only load variant metadata fields (no genotypes) for fast loading
    # Note: scikit-allel automatically extracts standard INFO fields like AF as variants/AF.
    # We request both variants/INFO (for structured access) and variants/AF (for direct access).
    fields_variant_only = [
        "variants/CHROM",
        "variants/POS",
        "variants/REF",
        "variants/ALT",
        "variants/ID",
        "variants/INFO",  # Load entire INFO field for structured array access
        "variants/AF",    # Also request AF directly (scikit-allel auto-extracts standard INFO fields)
    ]
    
    if region is not None:
        chr_, start, end = region
        # Convert chromosome to reference format
        region_chr_ref = mapper.sumstats_to_reference(chr_, reference_file=vcf_path, as_string=True)
        region_str = f"{region_chr_ref}:{start}-{end}"
        log.write(f" -Loading VCF region: {region_str}", verbose=verbose)
        vcf_data = read_vcf(vcf_path, region=region_str, tabix=tabix, fields=fields_variant_only)
    else:
        log.write(" -Loading all variants from VCF", verbose=verbose)
        vcf_data = read_vcf(vcf_path, tabix=tabix, fields=fields_variant_only)
    
    if vcf_data is None or len(vcf_data["variants/POS"]) == 0:
        raise ValueError("No variants found in VCF for the requested region.")
    
    # Extract variant information
    n_variants = len(vcf_data["variants/POS"])
    log.write(f" -Found {n_variants} variants", verbose=verbose)
    
    # Get chromosome (convert from reference format if needed)
    # Vectorized string conversion
    chr_data = vcf_data["variants/CHROM"]
    if isinstance(chr_data[0], bytes):
        # Vectorized bytes decoding
        chr_data = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in chr_data], dtype=object)
    else:
        chr_data = np.array([str(c) for c in chr_data], dtype=object)
    
    # Convert chromosome names back to sumstats format (vectorized)
    if vcf_chr_dict:
        reverse_dict = {v: k for k, v in vcf_chr_dict.items()}
        # Vectorized mapping using numpy
        chr_array = np.array(chr_data, dtype=object)
        chr_data = np.array([reverse_dict.get(c, c) for c in chr_array], dtype=object)
    
    # Extract positions
    pos_data = vcf_data["variants/POS"].astype(np.int64)
    
    # Extract ID (SNP ID) from VCF if available (vectorized bytes decoding)
    snpid_data = None
    if "variants/ID" in vcf_data:
        snpid_data = vcf_data["variants/ID"]
        if isinstance(snpid_data[0], (bytes, np.bytes_)):
            snpid_data = np.char.decode(snpid_data.astype("S"), "utf-8")
        else:
            snpid_data = snpid_data.astype(str)
        # Replace missing IDs ('.' or None) with None
        snpid_data = np.where((snpid_data == '.') | (snpid_data == 'None') | (snpid_data == ''), None, snpid_data)
    
    # Extract REF and ALT alleles (vectorized bytes decoding)
    ref_data = vcf_data["variants/REF"]
    if isinstance(ref_data[0], (bytes, np.bytes_)):
        ref_data = np.char.decode(ref_data.astype("S"), "utf-8")
    else:
        ref_data = ref_data.astype(str)
    
    alt_data = vcf_data["variants/ALT"]
    # ALT can be 2D array (multiple ALT alleles per variant)
    if alt_data.ndim == 2:
        # Take first ALT allele
        alt_data = alt_data[:, 0]
    if isinstance(alt_data[0], (bytes, np.bytes_)):
        alt_data = np.char.decode(alt_data.astype("S"), "utf-8")
    else:
        alt_data = alt_data.astype(str)
    
    # Extract EAF from INFO field (AF) if available, otherwise compute from genotypes
    # Note: scikit-allel automatically extracts standard INFO fields like AF as variants/AF,
    # not as variants/INFO["AF"]. Check both locations for compatibility.
    eaf_data = np.full(n_variants, np.nan, dtype=float)
    
    # First, try variants/AF (scikit-allel's automatic extraction)
    if "variants/AF" in vcf_data:
        af_info = vcf_data["variants/AF"]
        if af_info.ndim == 2:
            eaf_data = af_info[:, 0].astype(float)
        else:
            eaf_data = af_info.astype(float)
    # Fallback: try variants/INFO["AF"] (structured array access)
    elif "variants/INFO" in vcf_data and hasattr(vcf_data["variants/INFO"].dtype, 'names') and "AF" in vcf_data["variants/INFO"].dtype.names:
        af_info = vcf_data["variants/INFO"]["AF"]
        if af_info.ndim == 2:
            eaf_data = af_info[:, 0].astype(float)
        else:
            eaf_data = af_info.astype(float)
    
    # If EAF still missing, try to compute from genotypes (vectorized)
    if np.any(np.isnan(eaf_data)) and "calldata/GT" in vcf_data:
        # Compute AF from genotypes for missing values (vectorized)
        gt_array = GenotypeArray(vcf_data["calldata/GT"])
        n_alt = gt_array.to_n_alt()  # Integer dtype, -1 for missing
        # Convert to float and handle missing genotypes (scikit-allel uses -1 for missing, not NaN)
        n_alt = n_alt.astype(float)
        n_alt[n_alt < 0] = np.nan  # Map missing genotypes (-1) to NaN
        # Vectorized: compute allele frequency: mean of alt allele counts / 2 (NaN-aware)
        computed_eaf = np.nanmean(n_alt, axis=1) / 2.0
        # Vectorized: fill in missing values using boolean indexing
        nan_mask = np.isnan(eaf_data)
        eaf_data[nan_mask] = computed_eaf[nan_mask]
    
    # If still missing, use a default (0.5 for common variants) - vectorized
    nan_count = np.isnan(eaf_data).sum()
    if nan_count > 0:
        log.warning(f" -{nan_count} variants have missing EAF, using default 0.5", verbose=verbose)
        eaf_data[np.isnan(eaf_data)] = 0.5
    
    # Extract INFO score if available
    info_data = np.full(n_variants, np.nan, dtype=float)
    if "variants/INFO" in vcf_data and "INFO" in vcf_data["variants/INFO"].dtype.names:
        info_data = vcf_data["variants/INFO"]["INFO"].astype(float)
    
    # Create Polars DataFrame directly for better performance
    # Convert object arrays to proper types for Polars
    data_dict = {
        "CHR": [str(c) for c in chr_data],  # Ensure string type
        "POS": pos_data,
        "EA": [str(a) for a in alt_data],  # Ensure string type
        "NEA": [str(r) for r in ref_data],  # Ensure string type
        "EAF": eaf_data,
    }
    
    # Add SNPID if available from VCF
    if snpid_data is not None:
        data_dict["SNPID"] = [str(s) if s is not None else None for s in snpid_data]
    
    # Add INFO if available
    if not np.all(np.isnan(info_data)):
        data_dict["INFO"] = info_data
    
    # Create Polars DataFrame with explicit types using Series to avoid object dtype
    # This ensures proper types from the start, which is required for sorting
    df = pl.DataFrame({
        "CHR": pl.Series("CHR", data_dict["CHR"], dtype=pl.Utf8),
        "POS": pl.Series("POS", data_dict["POS"], dtype=pl.Int64),
        "EA": pl.Series("EA", data_dict["EA"], dtype=pl.Utf8),
        "NEA": pl.Series("NEA", data_dict["NEA"], dtype=pl.Utf8),
        "EAF": pl.Series("EAF", data_dict["EAF"], dtype=pl.Float64),
    })
    
    # Add optional columns with proper types if they exist
    if "SNPID" in data_dict:
        df = df.with_columns(pl.Series("SNPID", data_dict["SNPID"], dtype=pl.Utf8))
    if "INFO" in data_dict:
        df = df.with_columns(pl.Series("INFO", data_dict["INFO"], dtype=pl.Float64))
    
    # Sort by position (Polars is much faster)
    # All columns are now properly typed (Utf8, Int64, Float64), so sorting will work
    df = df.sort(["CHR", "POS"])
    
    log.write(f" -Loaded {df.height} variants", verbose=verbose)
    return df


# =============================================================================
# Helper functions for global simulation (efficient X^T X approach)
# =============================================================================

def _get_standardized_genotypes_for_block(
    vcf_path: str,
    variant_df: Union[pd.DataFrame, pl.DataFrame],
    start_idx: int,
    end_idx: int,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    vcf_chr_dict: Optional[dict] = None,
    tabix: Optional[bool] = None,
    mapper: Optional[ChromosomeMapper] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract and standardize genotypes for a block of variants.
    
    Returns standardized genotype matrix X where:
        X_ij = (G_ij - 2p_j) / sqrt(2p_j(1-p_j))
    
    Parameters
    ----------
    vcf_path : str
        Path to VCF file
    variant_df : pd.DataFrame
        DataFrame with variant information (must have CHR, POS, EA, NEA, EAF)
    start_idx : int
        Start index of block in variant_df
    end_idx : int
        End index of block in variant_df
    region : tuple, optional
        Original region specification
    vcf_chr_dict : dict, optional
        Chromosome mapping dictionary
    tabix : bool, optional
        Whether to use tabix indexing
    mapper : ChromosomeMapper, optional
        Chromosome mapper instance
    log : Log, optional
        Logging object
    verbose : bool
        Verbose flag
        
    Returns
    -------
    tuple
        (X, p, matched_indices) where:
        - X: standardized genotype matrix (n_samples x n_variants), standardized using sample mean and SD
        - p: effect allele frequencies (n_variants,) - returned for compatibility but not used in standardization
        - matched_indices: indices in variant_df that were matched (n_variants,)
    """
    if log is None:
        log = Log()
    
    if tabix is None:
        tabix = which("tabix")
    
    # Get or create mapper
    if mapper is None:
        mapper = ChromosomeMapper(log=log, verbose=verbose)
        mapper.detect_reference_format(vcf_path)
    
    # Get block variant data - extract numpy arrays directly (avoid Polars→Pandas conversion)
    if isinstance(variant_df, pl.DataFrame):
        block_slice = variant_df.slice(start_idx, end_idx - start_idx)
        if block_slice.height == 0:
            return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
        # Extract numpy arrays directly from Polars (much faster than converting to Pandas)
        block_pos = block_slice["POS"].to_numpy().astype(np.int64)
        block_nea = block_slice["NEA"].to_numpy()
        block_ea = block_slice["EA"].to_numpy()
        block_chr = block_slice["CHR"].to_numpy()
        block_eaf = block_slice["EAF"].to_numpy()
        # Need block_df for EAF lookup later, but only create if needed
        block_df = None  # Will create on-demand if needed
    else:
        block_df = variant_df.iloc[start_idx:end_idx].copy()
        if len(block_df) == 0:
            return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
        block_pos = block_df["POS"].values.astype(np.int64)
        block_nea = block_df["NEA"].values
        block_ea = block_df["EA"].values
        block_chr = block_df["CHR"].values
        block_eaf = block_df["EAF"].values
    
    # Determine region for this block
    chr_ = block_chr[0]
    start_pos = block_pos.min()
    end_pos = block_pos.max()
    
    # Convert chromosome to reference format
    region_chr_ref = mapper.sumstats_to_reference(chr_, reference_file=vcf_path, as_string=True)
    region_str = f"{region_chr_ref}:{start_pos}-{end_pos}"
    
    # Read VCF data for this block - minimal fields: POS, REF, ALT, GT
    vcf_data = read_vcf(
        vcf_path,
        region=region_str,
        tabix=tabix,
        fields=[
            "variants/POS",
            "variants/REF",
            "variants/ALT",
            "calldata/GT",
        ]
    )
    
    if vcf_data is None or len(vcf_data["variants/POS"]) == 0:
        log.warning(f" -No VCF data for block {start_idx}:{end_idx}", verbose=verbose)
        return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
    
    # Match variants (vectorized bytes decoding)
    vcf_pos = vcf_data["variants/POS"].astype(np.int64)
    vcf_ref = vcf_data["variants/REF"]
    if isinstance(vcf_ref[0], (bytes, np.bytes_)):
        vcf_ref = np.char.decode(vcf_ref.astype("S"), "utf-8")
    else:
        vcf_ref = vcf_ref.astype(str)
    
    vcf_alt = vcf_data["variants/ALT"]
    if vcf_alt.ndim == 2:
        vcf_alt = vcf_alt[:, 0]
    if isinstance(vcf_alt[0], (bytes, np.bytes_)):
        vcf_alt = np.char.decode(vcf_alt.astype("S"), "utf-8")
    else:
        vcf_alt = vcf_alt.astype(str)
    
    # Match variants by position only (vectorized using np.isin)
    # Since VCF ALT = EA and VCF REF = NEA (by definition), we only need to match positions
    # Vectorized matching: find positions in block_pos that exist in vcf_pos
    matched_mask = np.isin(block_pos, vcf_pos)
    matched_block_indices = np.where(matched_mask)[0]
    
    if len(matched_block_indices) == 0:
        log.warning(f" -No variants matched in VCF for block {start_idx}:{end_idx}", verbose=verbose)
        return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
    
    # Map block positions to VCF indices (take first match for each position)
    # Create sorted index for vcf_pos for efficient lookup
    vcf_pos_sorted_idx = np.argsort(vcf_pos)
    vcf_pos_sorted = vcf_pos[vcf_pos_sorted_idx]
    
    # Use searchsorted to find first match for each block position
    matched_positions = block_pos[matched_block_indices]
    vcf_match_idx = np.searchsorted(vcf_pos_sorted, matched_positions, side='left')
    # Ensure we don't go out of bounds
    vcf_match_idx = np.clip(vcf_match_idx, 0, len(vcf_pos_sorted) - 1)
    # Verify matches (positions must actually match)
    valid_match = vcf_pos_sorted[vcf_match_idx] == matched_positions
    matched_block_indices = matched_block_indices[valid_match]
    vcf_match_idx = vcf_match_idx[valid_match]
    # Map back to original vcf indices
    matched_indices = vcf_pos_sorted_idx[vcf_match_idx]
    
    if len(matched_indices) == 0:
        log.warning(f" -No variants matched in VCF for block {start_idx}:{end_idx}", verbose=verbose)
        return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
    
    matched_indices = np.array(matched_indices, dtype=int)
    matched_block_indices = np.array(matched_block_indices, dtype=int)
    
    # Get genotypes
    if "calldata/GT" not in vcf_data:
        log.warning(" -No genotype data in VCF", verbose=verbose)
        return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
    
    gt_array = GenotypeArray(vcf_data["calldata/GT"][matched_indices])
    n_alt = gt_array.to_n_alt()  # Shape: (n_variants, n_samples), integer dtype, -1 for missing
    # VCF ALT = EA, so n_alt directly represents EA dosage (0/1/2)
    
    # Convert to float and handle missing genotypes (scikit-allel uses -1 for missing, not NaN)
    n_alt = n_alt.astype(float)
    n_alt[n_alt < 0] = np.nan  # Map missing genotypes (-1) to NaN
    
    # Check for monomorphic variants (use nanstd to handle missing genotypes)
    variant_std = np.nanstd(n_alt, axis=1)
    polymorphic_mask = variant_std > 1e-10
    
    if polymorphic_mask.sum() == 0:
        log.warning(f" -All variants are monomorphic in block {start_idx}:{end_idx}", verbose=verbose)
        return np.array([]).reshape(0, 0), np.array([]), np.array([], dtype=int)
    
    # Filter to polymorphic variants
    n_alt_poly = n_alt[polymorphic_mask, :]
    matched_block_indices_poly = matched_block_indices[polymorphic_mask]
    
    # Get EAF from block data for matched variants
    # VCF ALT = EA, so EAF from VCF is already correct (no flipping needed)
    # Use pre-extracted block_eaf array (much faster than accessing block_df)
    if block_df is None:
        # Polars case: use pre-extracted array
        p = block_eaf[matched_block_indices_poly]
    else:
        # Pandas case: extract from DataFrame
        p = block_df.iloc[matched_block_indices_poly]["EAF"].values
    
    p = p.astype(float)
    
    # Standardize genotypes using sample statistics (not theoretical HWE values)
    # This ensures diag(R) = 1 exactly, preventing λGC inflation
    # 
    # Step 1: Transpose genotype matrix
    # n_alt is (n_variants, n_samples), we need (n_samples, n_variants) for standardization
    G = n_alt_poly.T.astype(np.float64, copy=False)  # Shape: (n_samples, n_variants)
    n_samples = G.shape[0]
    
    # Step 2: Compute sample-based standardization: X = (G - μ) / σ
    # where μ and σ are computed from the actual genotype matrix (not theoretical HWE)
    # This ensures each column has exactly mean=0 and var=1 in the sample
    # This is critical: theoretical standardization X = (G - 2p) / sqrt(2p(1-p)) 
    # only works in expectation under HWE and infinite samples. In finite samples,
    # it can lead to diag(R) ≠ 1 and λGC inflation.
    mu = np.nanmean(G, axis=0).astype(np.float64)  # Sample mean per variant (handles NaN)
    sd = np.nanstd(G, axis=0, ddof=0).astype(np.float64)  # Sample SD per variant (ddof=0 for population SD)
    
    # Step 3: Avoid division by zero for monomorphic variants
    # Should be filtered already, but safety check to prevent numerical errors
    sd = np.clip(sd, 1e-6, None).astype(np.float64)
    
    # Step 4: Standardize: X = (G - μ) / σ
    # Broadcast mu and sd across samples to standardize each variant
    X = ((G - mu[np.newaxis, :]) / sd[np.newaxis, :]).astype(np.float64)
    
    # Step 5: Handle missing data (NaN in genotypes)
    # Set to 0 (mean-centered) after standardization
    # After standardization, NaN values should be 0 since we subtract the mean
    X = np.where(np.isfinite(X), X, 0.0).astype(np.float64)
    
    return X, p, matched_block_indices_poly


# =============================================================================
# Helper functions for simulation workflow
# =============================================================================

def _filter_variants(
    df: pl.DataFrame,
    maf_min: float = 0.0,
    maf_max: float = 0.5,
    eaf_min: float = 0.0,
    eaf_max: float = 1.0,
    exclude_rare: bool = False,
    rare_maf_threshold: float = 0.01,
    log: Optional[Log] = None,
    verbose: bool = True
) -> pl.DataFrame:
    """
    Apply variant-level filters (MAF, EAF, rare variant exclusion) using Polars.
    
    Returns filtered Polars DataFrame.
    """
    if log is None:
        log = Log()
    
    m = df.height
    if m == 0:
        return df
    
    # Log filter parameters when non-default values are used
    if maf_min > 0.0 or maf_max < 0.5 or eaf_min > 0.0 or eaf_max < 1.0 or exclude_rare:
        log.write(f" -Applying filters: MAF=[{maf_min:.4f}, {maf_max:.4f}], EAF=[{eaf_min:.4f}, {eaf_max:.4f}], exclude_rare={exclude_rare}", verbose=verbose)
    
    # Step 1: Compute MAF from EAF
    # MAF (Minor Allele Frequency) = min(EAF, 1-EAF)
    # This ensures MAF is always between 0 and 0.5, regardless of which allele is the effect allele
    df = df.with_columns([
        (pl.min_horizontal(pl.col("EAF"), 1.0 - pl.col("EAF"))).alias("_MAF")
    ])
    
    # Step 2: Build filter conditions
    # We build conditions separately and combine them at the end
    conditions = []
    
    # Step 2a: MAF filter condition
    # Always apply MAF filter (even if maf_min=0.0 and maf_max=0.5, it's a no-op but ensures consistency)
    # Exclude NaN values from passing the filter (they fail the condition)
    maf_condition = (pl.col("_MAF").is_not_null()) & (pl.col("_MAF") >= maf_min) & (pl.col("_MAF") <= maf_max)
    conditions.append(maf_condition)
    # Count how many variants pass/fail for logging
    n_filtered_maf = df.filter(~maf_condition).height
    n_passed_maf = df.filter(maf_condition).height
    if maf_min > 0.0 or maf_max < 0.5:
        log.write(f" -MAF filter: {n_passed_maf} variants passed, {n_filtered_maf} filtered (MAF in [{maf_min:.4f}, {maf_max:.4f}])", verbose=verbose)
    elif n_filtered_maf > 0:
        # Even with default bounds, log if any were filtered (e.g., due to NaN)
        log.write(f" -MAF filter: {n_passed_maf} variants passed, {n_filtered_maf} filtered (MAF in [{maf_min:.4f}, {maf_max:.4f}])", verbose=verbose)
    
    # Step 2b: EAF filter condition
    # Always apply EAF filter (even if eaf_min=0.0 and eaf_max=1.0, it's a no-op but ensures consistency)
    # Exclude NaN values from passing the filter (they fail the condition)
    eaf_condition = (pl.col("EAF").is_not_null()) & (pl.col("EAF") >= eaf_min) & (pl.col("EAF") <= eaf_max)
    conditions.append(eaf_condition)
    # Count how many variants pass/fail for logging
    n_filtered_eaf = df.filter(~eaf_condition).height
    n_passed_eaf = df.filter(eaf_condition).height
    if eaf_min > 0.0 or eaf_max < 1.0:
        log.write(f" -EAF filter: {n_passed_eaf} variants passed, {n_filtered_eaf} filtered (EAF in [{eaf_min:.4f}, {eaf_max:.4f}])", verbose=verbose)
    elif n_filtered_eaf > 0:
        # Even with default bounds, log if any were filtered (e.g., due to NaN)
        log.write(f" -EAF filter: {n_passed_eaf} variants passed, {n_filtered_eaf} filtered (EAF in [{eaf_min:.4f}, {eaf_max:.4f}])", verbose=verbose)
    
    # Step 2c: Rare variant exclusion (optional)
    # If exclude_rare=True, exclude variants with MAF below the threshold
    if exclude_rare:
        rare_condition = pl.col("_MAF") >= rare_maf_threshold
        conditions.append(rare_condition)
        n_rare = df.filter(~rare_condition).height
        if n_rare > 0:
            log.write(f" -Excluded {n_rare} rare variants (MAF < {rare_maf_threshold:.4f})", verbose=verbose)
    
    # Step 3: Apply all filter conditions
    # Combine all conditions with AND logic (all must be satisfied)
    if conditions:
        combined_condition = conditions[0]
        for cond in conditions[1:]:
            combined_condition = combined_condition & cond
        df = df.filter(combined_condition)
    
    # Step 4: Clean up temporary MAF column
    df = df.drop("_MAF")
    
    # Step 5: Log final filter results
    n_before_filter = m
    m_filtered = df.height
    n_filtered = n_before_filter - m_filtered
    
    if n_filtered > 0:
        log.write(f" -After filtering: {m_filtered} variants remaining ({n_filtered} excluded)", verbose=verbose)
    
    if m_filtered == 0:
        raise ValueError("No variants remaining after applying filters. Please adjust filter criteria.")
    
    return df


def _setup_sample_sizes(
    m: int,
    trait: str,
    n: Optional[int] = None,
    n_case: Optional[int] = None,
    n_ctrl: Optional[int] = None,
    n_drop_rate: float = 0.0,
    n_drop_min: float = 0.5,
    rng: Optional[np.random.Generator] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Setup sample sizes (N) and effective sample sizes (N_eff) for variants.
    
    Returns
    -------
    tuple
        (N, N_eff) arrays of length m
    """
    if rng is None:
        rng = np.random.default_rng()
    
    # Step 1: Set baseline N and N_eff based on trait type
    if trait == "binary":
        if n_case is None or n_ctrl is None:
            raise ValueError("For trait='binary', you must provide n_case and n_ctrl.")
        # For binary traits:
        # - Total N = n_case + n_ctrl
        # - Effective N accounts for unequal case/control sizes: N_eff ≈ 4 / (1/N_case + 1/N_ctrl)
        #   This formula comes from the fact that variance in case-control studies depends on both sample sizes
        neff_base = 4.0 / (1.0 / float(n_case) + 1.0 / float(n_ctrl))
        N = np.full(m, float(n_case + n_ctrl), dtype=float)
        N_eff = np.full(m, neff_base, dtype=float)
    else:
        if n is None:
            raise ValueError("For trait='quant', you must provide n.")
        # For quantitative traits:
        # - N = n (user-specified sample size)
        # - N_eff = N (no adjustment needed for quantitative traits)
        N = np.full(m, float(n), dtype=float)
        N_eff = N.copy()
    
    # Step 2: Apply per-SNP N variation to simulate missingness/QC issues
    # In real GWAS, not all SNPs have the same sample size due to:
    # - Missing genotypes
    # - Quality control filtering
    # - Different imputation quality
    if n_drop_rate > 0:
        # Randomly select n_drop_rate fraction of SNPs to have reduced N
        drop = rng.random(m) < n_drop_rate
        n_dropped = drop.sum()
        if n_dropped > 0:
            # Multiply N by a random factor in [n_drop_min, 1.0]
            mult = rng.uniform(n_drop_min, 1.0, size=n_dropped)
            N[drop] *= mult
            # For quantitative traits, N_eff follows N; for binary we keep neff_base constant
            if trait == "quant":
                N_eff[drop] *= mult
    
    return N, N_eff


def _setup_info_scores(
    df: pl.DataFrame,
    use_info: bool = True,
    info_mean: float = 0.95,
    info_sd: float = 0.05,
    info_min: float = 0.3,
    info_max: float = 1.0,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Setup INFO scores for variants (from VCF or simulated).
    
    Returns
    -------
    np.ndarray
        INFO scores array of length m
    """
    if rng is None:
        rng = np.random.default_rng()
    
    m = df.height
    
    if use_info and "INFO" in df.columns:
        info = np.clip(df["INFO"].to_numpy().astype(float), info_min, info_max)
    else:
        info = rng.normal(info_mean, info_sd, size=m)
        info = np.clip(info, info_min, info_max)
    
    return info


def _select_causal_variants(
    df: pl.DataFrame,
    mode: str,
    pi: float,
    n_causal: int,
    alpha: float,
    effect_sd: float,
    causal_maf_min: Optional[float] = None,
    causal_maf_max: Optional[float] = None,
    causal_eaf_min: Optional[float] = None,
    causal_eaf_max: Optional[float] = None,
    causal_idx: Optional[np.ndarray] = None,
    causal_beta: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Select causal variants and assign effect sizes.
    
    Returns
    -------
    tuple
        (beta_true, is_causal) arrays of length m
    """
    if log is None:
        log = Log()
    if rng is None:
        rng = np.random.default_rng()
    
    m = df.height
    beta_true = np.zeros(m, dtype=float)  # True causal effect sizes (0 for non-causal)
    is_causal = np.zeros(m, dtype=bool)    # Boolean mask indicating which variants are causal
    
    # Step 1: Compute MAF and EAF for causal variant filtering
    # Extract EAF values and compute MAF = min(EAF, 1-EAF)
    eaf_values = df["EAF"].to_numpy().astype(float)
    maf_values = np.minimum(eaf_values, 1.0 - eaf_values)
    
    # Step 2: Build causal variant candidate mask
    # Apply causal_maf_min/max and causal_eaf_min/max filters to restrict which variants
    # are eligible to be selected as causal (separate from variant-level filters)
    causal_candidate_mask = np.ones(m, dtype=bool)
    
    # Apply causal variant filters if specified
    if causal_maf_min is not None or causal_maf_max is not None:
        causal_maf_mask = np.ones(m, dtype=bool)
        if causal_maf_min is not None:
            causal_maf_mask = causal_maf_mask & (maf_values >= causal_maf_min)
        if causal_maf_max is not None:
            causal_maf_mask = causal_maf_mask & (maf_values <= causal_maf_max)
        causal_candidate_mask = causal_candidate_mask & causal_maf_mask
        n_excluded_maf = (~causal_maf_mask).sum()
        if n_excluded_maf > 0:
            log.write(f" -Excluded {n_excluded_maf} variants from causal selection by MAF filter", verbose=verbose)
    
    if causal_eaf_min is not None or causal_eaf_max is not None:
        causal_eaf_mask = np.ones(m, dtype=bool)
        if causal_eaf_min is not None:
            causal_eaf_mask = causal_eaf_mask & (eaf_values >= causal_eaf_min)
        if causal_eaf_max is not None:
            causal_eaf_mask = causal_eaf_mask & (eaf_values <= causal_eaf_max)
        causal_candidate_mask = causal_candidate_mask & causal_eaf_mask
        n_excluded_eaf = (~causal_eaf_mask).sum()
        if n_excluded_eaf > 0:
            log.write(f" -Excluded {n_excluded_eaf} variants from causal selection by EAF filter", verbose=verbose)
    
    # Get candidate indices
    causal_candidate_indices = np.where(causal_candidate_mask)[0]
    n_candidates = len(causal_candidate_indices)
    
    if n_candidates == 0:
        raise ValueError("No variants eligible for causal selection after applying causal filters. Please adjust filter criteria.")
    
    # Step 3: Select causal variants and assign effect sizes
    if causal_idx is not None and causal_beta is not None:
        # Option A: Use explicit causal variants provided by user
        # Validate inputs
        causal_idx = np.asarray(causal_idx, dtype=int)
        causal_beta = np.asarray(causal_beta, dtype=float)
        if causal_idx.ndim != 1 or causal_beta.ndim != 1 or len(causal_idx) != len(causal_beta):
            raise ValueError("causal_idx and causal_beta must be 1D and same length.")
        if np.any((causal_idx < 0) | (causal_idx >= m)):
            raise ValueError("causal_idx contains indices outside the variant range.")
        
        # Check if explicit causal indices pass the causal filters (warn if not)
        if not np.all(causal_candidate_mask[causal_idx]):
            invalid_causal = causal_idx[~causal_candidate_mask[causal_idx]]
            log.warning(f" -Warning: {len(invalid_causal)} explicit causal variants do not pass causal filters", verbose=verbose)
        
        # Assign explicit effect sizes
        beta_true[causal_idx] = causal_beta
        is_causal[causal_idx] = True
    else:
        # Option B: Select causals based on mode (polygenic or sparse)
        if mode == "sparse":
            # Sparse mode: Select exactly n_causal variants uniformly at random
            k = min(int(n_causal), n_candidates)
            if k < int(n_causal):
                log.warning(f" -Requested {n_causal} causal variants but only {n_candidates} candidates available, selecting {k}", verbose=verbose)
            selected_candidate_idx = rng.choice(n_candidates, size=k, replace=False)
            idx = causal_candidate_indices[selected_candidate_idx]
        elif mode == "polygenic":
            # Polygenic mode: Select m * pi variants (each variant is causal with probability pi)
            k = max(1, int(round(n_candidates * float(pi))))
            if k > n_candidates:
                k = n_candidates
                log.warning(f" -Requested {int(round(n_candidates * float(pi)))} causal variants but only {n_candidates} candidates available, selecting {k}", verbose=verbose)
            selected_candidate_idx = rng.choice(n_candidates, size=k, replace=False)
            idx = causal_candidate_indices[selected_candidate_idx]
        else:
            raise ValueError("mode must be 'polygenic' or 'sparse'")
        
        # Step 4: Sample causal effect sizes with MAF-dependent architecture
        if alpha > 0.0:
            # MAF-dependent architecture: β_j^raw ~ N(0, [2p_j(1-p_j)]^(-alpha))
            # where p_j is the effect allele frequency (EAF)
            # This means rarer variants have larger effect sizes on average
            # Use EAF (p_j) to match the standardization formula
            p_causal = eaf_values[idx]
            # Variance formula: [2p_j(1-p_j)]^(-alpha) where p_j is EAF
            # For alpha > 0, rarer variants (smaller p) have larger variance
            variance_scale = np.power(2.0 * p_causal * (1.0 - p_causal), -alpha)
            variance_scale = np.clip(variance_scale, 1e-10, 1e10)  # Avoid extreme values (numerical stability)
            beta_true[idx] = rng.normal(0.0, np.sqrt(variance_scale), size=len(idx))
        else:
            # MAF-independent architecture: β_j ~ N(0, effect_sd^2)
            # All causal variants have the same effect size variance regardless of frequency
            beta_true[idx] = rng.normal(0.0, float(effect_sd), size=len(idx))
        
        is_causal[idx] = True
    
    log.write(f" -Selected {is_causal.sum()} causal variants", verbose=verbose)
    
    return beta_true, is_causal


def _convert_z_to_sumstats(
    z: np.ndarray,
    N_eff: np.ndarray,
    log: Optional[Log] = None,
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert Z-scores to summary statistics (SE, BETA, P, MLOG10P).
    
    Returns
    -------
    tuple
        (se, beta_hat, p, mlog10p) arrays
    """
    if log is None:
        log = Log()
    
    m = len(z)
    
    # Validate N_eff
    N_eff = np.clip(N_eff, 1.0, None)
    
    # Standard error: SE = 1 / sqrt(N_eff)
    se = 1.0 / np.sqrt(N_eff)
    
    # Validate Z scores
    z_nonzero_count = np.sum(z != 0.0)
    z_finite_count = np.sum(np.isfinite(z))
    z_abs_max = np.abs(z[np.isfinite(z)]).max() if z_finite_count > 0 else 0.0
    z_abs_mean = np.abs(z[np.isfinite(z)]).mean() if z_finite_count > 0 else 0.0
    
    log.write(f" -Z-score statistics: {z_nonzero_count}/{m} non-zero, {z_finite_count}/{m} finite, max|Z|={z_abs_max:.6f}, mean|Z|={z_abs_mean:.6f}", verbose=verbose)
    
    if z_nonzero_count == 0:
        log.warning(" -WARNING: All Z scores are zero! This will result in all betas being zero.", verbose=verbose)
        log.warning(" -Possible causes:", verbose=verbose)
        log.warning("   1. All causal effects are zero (check causal variant selection)", verbose=verbose)
        log.warning("   2. LD matrices are all identity (variants not matching in VCF)", verbose=verbose)
        log.warning("   3. Numerical issues causing all Z to be replaced with 0", verbose=verbose)
    
    # Replace NaN/Inf Z scores with 0 (no effect)
    z_valid = np.where(np.isfinite(z), z, 0.0)
    
    # Effect estimate: BETA = Z * SE
    beta_hat = z_valid * se
    
    # Calculate MLOG10P directly from Z using log-space for numerical precision
    mlog10p = _z_to_mlog10p(z_valid)
    mlog10p = np.where(np.isfinite(mlog10p), mlog10p, 0.0)
    
    # Calculate P from MLOG10P: P = 10^(-MLOG10P)
    p = np.power(10.0, -mlog10p)
    p = np.where(np.isfinite(p), p, 0.5)
    
    # Log warning if there were invalid Z scores
    invalid_z_count = np.sum(~np.isfinite(z))
    if invalid_z_count > 0:
        log.warning(f" -Found {invalid_z_count} invalid Z scores (NaN/Inf), replaced with 0", verbose=verbose)
    
    return se, beta_hat, p, mlog10p


def _extract_causal_snp_ids(
    df: pl.DataFrame,
    is_causal: np.ndarray,
    log: Optional[Log] = None,
    verbose: bool = True
) -> List[str]:
    """
    Extract causal SNP IDs from Polars DataFrame.
    
    Uses VCF ID field if available, otherwise constructs as "CHR:POS:EA:NEA".
    """
    if log is None:
        log = Log()
    
    causal_indices = np.where(is_causal)[0]
    if len(causal_indices) == 0:
        return []
    
    # Filter to causal variants only (Polars is fast)
    # Use row indices to filter - create a boolean mask
    # Polars int_range().is_in() works with lists or Series
    is_causal_mask = pl.int_range(0, df.height).is_in(causal_indices)
    causal_df = df.filter(is_causal_mask)
    
    # Extract SNP IDs
    causal_snp_ids = []
    
    # Check if SNPID column exists (from VCF)
    if "SNPID" in df.columns:
        snpid_col = causal_df["SNPID"].to_list()
        chr_col = causal_df["CHR"].to_list()
        pos_col = causal_df["POS"].to_list()
        ea_col = causal_df["EA"].to_list()
        nea_col = causal_df["NEA"].to_list()
        
        for snp_id, chr_val, pos_val, ea_val, nea_val in zip(snpid_col, chr_col, pos_col, ea_col, nea_col):
            # Use VCF ID if available and not missing
            if snp_id is not None and str(snp_id) != 'None' and str(snp_id) != '.' and str(snp_id) != '':
                causal_snp_ids.append(str(snp_id))
            else:
                # Fallback: construct ID from CHR:POS:EA:NEA
                snp_id = f"{chr_val}:{pos_val}:{ea_val}:{nea_val}"
                causal_snp_ids.append(snp_id)
    else:
        # No SNPID column: construct from CHR:POS:EA:NEA
        chr_col = causal_df["CHR"].to_list()
        pos_col = causal_df["POS"].to_list()
        ea_col = causal_df["EA"].to_list()
        nea_col = causal_df["NEA"].to_list()
        
        for chr_val, pos_val, ea_val, nea_val in zip(chr_col, pos_col, ea_col, nea_col):
            snp_id = f"{chr_val}:{pos_val}:{ea_val}:{nea_val}"
            causal_snp_ids.append(snp_id)
    
    log.write(f" -Extracted {len(causal_snp_ids)} causal SNP IDs", verbose=verbose)
    
    return causal_snp_ids


# =============================================================================
# Main functions
# =============================================================================

def simulate_sumstats_region(
    vcf_path: str,
    *,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    # trait and sample size
    trait: str = "quant",                    # "quant" or "binary"
    n: Optional[int] = 300_000,              # required if trait="quant"
    n_case: Optional[int] = None,            # required if trait="binary"
    n_ctrl: Optional[int] = None,            # required if trait="binary"
    # genetic architecture
    mode: str = "sparse",                 # "polygenic" or "sparse"
    pi: float = 2e-3,                        # causal fraction if polygenic
    n_causal: int =1,                      # number of causals if sparse
    effect_sd: float = 0.05,                 # SD of causal effects (standardized units, used if alpha=0)
    alpha: float = 0,                       # MAF dependence parameter (0 = MAF-independent, uses effect_sd)
    # variant filtering
    maf_min: float = 0.01,                    # minimum MAF for variants (0.0 = no filter)
    maf_max: float = 0.5,                    # maximum MAF for variants (0.5 = no filter)
    eaf_min: float = 0.0,                    # minimum EAF for variants (0.0 = no filter)
    eaf_max: float = 1.0,                    # maximum EAF for variants (1.0 = no filter)
    exclude_rare: bool = False,              # if True, exclude rare variants (MAF < 0.01)
    rare_maf_threshold: float = 0.01,        # MAF threshold for rare variants
    thin: Optional[float] = 0.5,            # fraction of variants to keep after filtering (None = keep all, simulates incomplete coverage)
    # causal variant selection filters
    causal_maf_min: Optional[float] = None,  # minimum MAF for causal variants (None = no filter)
    causal_maf_max: Optional[float] = None,  # maximum MAF for causal variants (None = no filter)
    causal_eaf_min: Optional[float] = None,  # minimum EAF for causal variants (None = no filter)
    causal_eaf_max: Optional[float] = None,  # maximum EAF for causal variants (None = no filter)
    # realism knobs
    n_drop_rate: float = 0.15,               # fraction of SNPs with reduced N
    n_drop_min: float = 0.5,                 # minimum N multiplier for dropped SNPs
    use_info: bool = True,                   # if INFO column exists, use it
    info_mean: float = 0.95,                 # if simulating INFO
    info_sd: float = 0.05,
    info_min: float = 0.3,
    info_max: float = 1.0,
    # confounding / inflation parameters
    lambda_gc: float = 1.01,              # cryptic relatedness / global inflation (λ)
    sigma_strat: float = 0.05,            # population stratification (σ_strat)
    # deterministic
    seed: int = 1,
    # explicit causals (optional)
    causal_idx: Optional[np.ndarray] = None,
    causal_beta: Optional[np.ndarray] = None,
    # VCF reading options
    vcf_chr_dict: Optional[dict] = None,
    tabix: Optional[bool] = None,
    verbose: bool = True,
    log: Optional[Log] = None,
    # Sumstats object options
    study: str = "Simulated_Study",
    trait_name: str = "Simulated_Trait",
    build: str = "38",
) -> Tuple[Sumstats, List[str]]:
    """
    Simulate GWAS summary statistics (BETA/SE/Z/P) for a genomic region using a reference panel LD model.
    
    This function uses an efficient X^T X approach to simulate Z-scores without explicitly constructing
    large LD matrices. It supports MAF-dependent effect size architectures and processes the entire
    region at once for efficient simulation.

    Parameters
    ----------
    vcf_path : str
        Path to VCF/BCF file containing reference panel genotypes. The VCF should contain:
        - Genotype data (calldata/GT) for LD computation
        - AF field in INFO (or variants/AF) for effect allele frequency
        - Optional: INFO field for imputation quality scores
        
    region : tuple, optional
        Optional tuple (chr, start, end) specifying genomic region to simulate.
        If provided, only variants in this region are used.
        Example: region=("1", 100_000, 2_000_000)
        If None, uses all variants in the VCF.
        
    trait : str, default="quant"
        Trait type: "quant" for quantitative trait or "binary" for case-control
        
    n : int, optional, default=300000
        Sample size for quantitative trait (required if trait="quant")
        
    n_case : int, optional
        Number of cases for binary trait (required if trait="binary")
        
    n_ctrl : int, optional
        Number of controls for binary trait (required if trait="binary")
        
    mode : str, default="sparse"
        Genetic architecture mode:
        - "polygenic": Causal variants selected with probability pi
        - "sparse": Exactly n_causal causal variants selected
        
    pi : float, default=2e-3
        Causal fraction if mode="polygenic" (probability each variant is causal)
        
    n_causal : int, default=1
        Number of causal variants if mode="sparse"
        
    effect_sd : float, default=0.05
        Standard deviation of causal effects (standardized units).
        Only used when alpha=0 (MAF-independent effects).
        When alpha>0, effect size variance scales as [2p(1-p)]^(-alpha) and effect_sd is ignored.
        
    alpha : float, default=0
        MAF dependence parameter for effect size architecture.
        Effect size variance scales as [2p_j(1-p_j)]^(-alpha) where p_j is the effect
        allele frequency (EAF). 
        - alpha=0: MAF-independent effects (uses effect_sd)
        - alpha>0: Rarer variants have larger effect sizes on average
        Typical values: 0.0 (no dependence), 0.2 (moderate), 0.5 (strong)
        
    maf_min : float, default=0.01
        Minimum MAF for variants to include in simulation.
        Variants with MAF < maf_min are excluded. Set to 0.0 for no filter.
        
    maf_max : float, default=0.5
        Maximum MAF for variants to include in simulation.
        Variants with MAF > maf_max are excluded. Set to 0.5 for no filter.
        
    eaf_min : float, default=0.0
        Minimum EAF for variants to include in simulation.
        Variants with EAF < eaf_min are excluded. Set to 0.0 for no filter.
        
    eaf_max : float, default=1.0
        Maximum EAF for variants to include in simulation.
        Variants with EAF > eaf_max are excluded. Set to 1.0 for no filter.
        
    exclude_rare : bool, default=False
        If True, exclude rare variants with MAF < rare_maf_threshold
        
    rare_maf_threshold : float, default=0.01
        MAF threshold for rare variants when exclude_rare=True
        
    thin : float, optional
        Fraction of variants to randomly keep after filtering (between 0 and 1, exclusive).
        If None (default), keeps all variants after filtering.
        If specified (e.g., 0.1), randomly samples that fraction of variants to simulate
        incomplete GWAS coverage (real GWAS don't cover all sites).
        Useful for simulating genotyping arrays or imputation panels that only cover
        a subset of variants.
        
    causal_maf_min : float, optional
        Minimum MAF for variants eligible to be selected as causal (None = no filter).
        This is separate from maf_min which filters all variants.
        
    causal_maf_max : float, optional
        Maximum MAF for variants eligible to be selected as causal (None = no filter)
        
    causal_eaf_min : float, optional
        Minimum EAF for variants eligible to be selected as causal (None = no filter)
        
    causal_eaf_max : float, optional
        Maximum EAF for variants eligible to be selected as causal (None = no filter)
        
    n_drop_rate : float, default=0.15
        Fraction of SNPs with reduced sample size (simulates missingness/QC issues)
        
    n_drop_min : float, default=0.5
        Minimum N multiplier for dropped SNPs (N is multiplied by random factor in [n_drop_min, 1.0])
        
    use_info : bool, default=True
        If True and INFO column exists in VCF, use it for N_eff attenuation.
        Otherwise, simulates INFO scores.
        
    info_mean : float, default=0.95
        Mean INFO score if simulating (when INFO not in VCF or use_info=False)
        
    info_sd : float, default=0.05
        Standard deviation of INFO score if simulating
        
    info_min : float, default=0.3
        Minimum INFO score (clipped to this value)
        
    info_max : float, default=1.0
        Maximum INFO score (clipped to this value)
        
    lambda_gc : float, default=1.01
        Cryptic relatedness / global inflation parameter (λ).
        Scales the noise term: ε = λ * ε_0 where ε_0 ~ N(0, R).
        - lambda_gc = 1.0: No inflation (default)
        - lambda_gc > 1.0: Genomic control inflation (e.g., 1.0-1.1 mild, 1.1-1.3 noticeable)
        - Captures global effects of relatedness/unmodeled covariance in GWAS residuals
        
    sigma_strat : float, default=0.05
        Population stratification parameter (σ_strat).
        Adds LD-correlated bias: b_strat ~ N(0, σ_strat^2 R).
        - sigma_strat = 0.0: No stratification (default)
        - sigma_strat > 0.0: Creates LD-shaped structured confounding
        - Typical values: 0.0-0.2 (small), 0.2-0.5 (moderate), 0.5-1.0 (strong)
        - Creates genome-wide shifts and correlated "hills" of signal across LD blocks
        
    seed : int, default=1
        Random seed for reproducibility
        
    causal_idx : np.ndarray, optional
        Explicit indices of causal variants (0-based, relative to filtered variant order).
        If provided, overrides mode/pi/n_causal selection.
        
    causal_beta : np.ndarray, optional
        Explicit causal effect sizes (must match causal_idx length).
        If provided, overrides effect_sd/alpha effect size sampling.
        
    vcf_chr_dict : dict, optional
        Chromosome mapping dictionary for VCF. If None, auto-detected.
        Format: {sumstats_chr_format: vcf_chr_format}, e.g., {"1": "chr1"}
        
    tabix : bool, optional
        Whether to use tabix indexing for VCF. If None, auto-detected based on file extension.
        
    verbose : bool, default=True
        Verbose output with progress messages
        
    log : Log, optional
        Logging object. If None, creates a new one.
        
    study : str, default="Simulated_Study"
        Study name for Sumstats object metadata
        
    trait_name : str, default="Simulated_Trait"
        Trait name for Sumstats object metadata
        
    build : str, default="38"
        Genome build version for Sumstats object ("37" or "38")

    Returns
    -------
    tuple
        A tuple containing:
        - Sumstats: GWASLab Sumstats object with simulated summary statistics.
          Contains columns: CHR, POS, EA, NEA, EAF, BETA, SE, Z, P, MLOG10P, N, N_EFF, INFO,
          IS_CAUSAL, BETA_TRUE
        - List[str]: List of causal SNP IDs. Uses VCF ID field if available,
          otherwise constructed as "CHR:POS:EA:NEA"
        
    Notes
    -----
    **Workflow Overview**
    
    The function implements an efficient workflow for simulating GWAS summary statistics:
    
    1. **Variant Loading and Filtering**:
       - Reads variants from VCF/BCF file using Polars for fast DataFrame operations
       - Extracts: CHR, POS, EA (effect allele = VCF ALT), NEA (non-effect allele = VCF REF), 
         EAF (effect allele frequency from INFO/AF field)
       - If EAF missing from INFO, attempts to compute from genotypes; if still missing, 
         defaults to 0.5 (with warning)
       - Applies variant-level filters (always applied, even with default bounds):
         * MAF filters: excludes variants with MAF outside [maf_min, maf_max]
         * EAF filters: excludes variants with EAF outside [eaf_min, eaf_max]
         * Rare variant exclusion: if exclude_rare=True, removes variants with MAF < rare_maf_threshold
       - Computes MAF from EAF: MAF = min(EAF, 1-EAF)
       - Filters out variants with missing/null EAF values
       - If thin is specified: randomly samples a fraction of variants to simulate incomplete
         GWAS coverage (real GWAS don't cover all sites, only genotyped/imputed variants)
    
    2. **Sample Size and Effective Sample Size Setup**:
       - For quantitative traits: N = n (user-specified), N_eff = N
       - For binary traits: N = n_case + n_ctrl, N_eff = 4 / (1/n_case + 1/n_ctrl)
       - Applies per-SNP N variation to simulate missingness/QC:
         * Randomly selects n_drop_rate fraction of SNPs
         * Multiplies their N by a random factor in [n_drop_min, 1.0]
         * For quantitative traits, N_eff follows N; for binary, N_eff remains constant
    
    3. **INFO Score and Imputation Quality**:
       - If INFO column exists in VCF and use_info=True: uses actual INFO scores
       - Otherwise: simulates INFO scores from Normal(info_mean, info_sd)
       - Clips INFO to [info_min, info_max] range
       - Attenuates N_eff: N_eff <- N_eff * INFO
       - Lower INFO (imputation uncertainty) → larger SE → smaller Z-scores
    
    4. **Causal Variant Selection**:
       - If causal_idx provided: uses explicit causal variants
       - Otherwise, selects based on mode:
         * "polygenic": each variant is causal with probability pi
         * "sparse": exactly n_causal variants selected uniformly at random
       - Applies causal_maf_min/max and causal_eaf_min/max filters to candidate pool
       - If causal_beta provided: uses explicit effect sizes
       - Otherwise, samples effect sizes:
         * If alpha=0: beta_j ~ Normal(0, effect_sd^2)
         * If alpha>0: beta_j ~ Normal(0, [2p_j(1-p_j)]^(-alpha)) where p_j is EAF
    
    5. **Z-Score Generation (X^T X Method)**:
       - Processes entire region at once (no block-wise division)
       - Loads genotypes from VCF and standardizes: X_ij = (G_ij - μ_j) / σ_j
         where μ_j and σ_j are sample mean and standard deviation
       - Computes causal mean: μ_causal = sqrt(N_eff) * (1/n) X^T (X β)
       - Generates standard LD-noise: ε_0 = (1/sqrt(n)) X^T u, where u ~ Normal(0, I_n)
       - Applies cryptic relatedness inflation: ε = λ * ε_0 (scales noise by lambda_gc)
       - Generates population stratification bias: b_strat = σ_strat * (1/sqrt(n)) X^T v, where v ~ Normal(0, I_n)
       - Z-scores: z = μ_causal + b_strat + ε
       - This avoids explicit LD matrix construction (O(k²) memory) and Cholesky decomposition (O(k³) time)
    
    6. **Summary Statistics Conversion**:
       - Converts Z-scores to BETA, SE, P, MLOG10P using standard formulas:
         * SE = 1 / sqrt(N_eff)  [standard error]
         * BETA = Z * SE  [effect estimate, ensures Z = BETA/SE]
         * MLOG10P = -log10(P)  [negative log10 p-value, computed in log-space for precision]
         * P = 10^(-MLOG10P)  [two-sided p-value, derived from MLOG10P]
       - Handles invalid Z-scores (NaN/Inf) by replacing with 0
    
    7. **Output Assembly**:
       - Creates Polars DataFrame with all summary statistics
       - Adds columns: INFO, N, N_EFF, IS_CAUSAL, BETA_TRUE, Z, SE, BETA, P, MLOG10P
       - For binary traits, also adds N_CASE and N_CONTROL
       - Converts to Pandas for Sumstats object creation
       - Returns Sumstats object and list of causal SNP IDs
    
    Examples
    --------
    >>> from gwaslab import simulate_sumstats_region
    >>> 
    >>> # Simulate a 2Mb region with default settings
    >>> sumstats, causals = simulate_sumstats_region(
    ...     vcf_path="reference_panel.vcf.gz",
    ...     region=("1", 100000, 2100000),
    ...     n=100000,
    ...     n_causal=10,
    ...     maf_min=0.01
    ... )
    >>> 
    >>> # Simulate with MAF-dependent effects
    >>> sumstats, causals = simulate_sumstats_region(
    ...     vcf_path="reference_panel.vcf.gz",
    ...     region=("1", 100000, 2100000),
    ...     n=100000,
    ...     mode="polygenic",
    ...     pi=1e-3,
    ...     alpha=0.2,  # Rarer variants have larger effects
    ...     maf_min=0.01,
    ...     maf_max=0.5
    ... )
    
    See Also
    --------
    simulate_sumstats_global : Simulate summary statistics genome-wide with global heritability calibration
    
    Notes
    -----
    **Mathematical Model**
    
    The simulation uses a reference-panel LD model to generate Z-scores:
    
    z = μ_causal + b_strat + ε
    
    where:
    - μ_causal = sqrt(N_eff) * (1/n) X^T (X β)  [causal mean: effects spread via LD]
    - ε = λ * (1/sqrt(n)) X^T u, u ~ N(0, I_n)  [noise: scaled by λ for cryptic relatedness]
    - b_strat = σ_strat * (1/sqrt(n)) X^T v, v ~ N(0, I_n)  [bias: population stratification]
    - X: standardized genotype matrix, X_ij = (G_ij - μ_j) / σ_j
    - β: true causal effect sizes (mostly zeros, non-zero for causal variants)
    - N_eff: effective sample size (adjusted for trait type and INFO score)
    - λ: cryptic relatedness / global inflation parameter (lambda_gc)
    - σ_strat: population stratification parameter (sigma_strat)
    
    The final distribution is: z ~ N(μ_causal, (λ + σ_strat^2) R) where R = (1/n) X^T X.
    This avoids explicit construction of R (O(k²) memory) and Cholesky decomposition (O(k³) time).
    
    **Key Features**
    
    1. **Sample-based standardization**: Uses sample mean and SD (not theoretical HWE)
       to ensure diag(R) = 1 exactly, preventing λGC inflation.
    
    2. **MAF-dependent effects**: When alpha > 0, effect size variance scales as
       [2p(1-p)]^(-alpha), meaning rarer variants have larger effects on average.
    
    3. **Realistic sample size variation**: Simulates missingness/QC by randomly
       reducing N for a fraction of SNPs.
    
    4. **Imputation quality**: Attenuates N_eff by INFO score to simulate
       imputation uncertainty (lower INFO → larger SE → smaller Z).
    
    **Performance**
    
    - Uses Polars for fast DataFrame operations
    - Processes entire region at once (no block-wise division)
    - Vectorized NumPy operations throughout
    - Supports tabix-indexed VCF files for fast region queries
    """
    if log is None:
        log = Log()
    
    rng = np.random.default_rng(seed)
    
    log.write("Starting GWAS summary statistics simulation...", verbose=verbose)
    log.write(f" -VCF path: {vcf_path}", verbose=verbose)
    if region:
        log.write(f" -Region: chr{region[0]}:{region[1]}-{region[2]}", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # Step 0: Load variants from VCF file
    # -------------------------------------------------------------------------
    # Initialize chromosome mapper to handle different chromosome naming conventions
    # (e.g., "1" vs "chr1", "chr1" vs "1" in VCF vs sumstats format)
    log.write("Loading variants from VCF...", verbose=verbose)
    mapper = ChromosomeMapper(log=log, verbose=verbose)
    mapper.detect_reference_format(vcf_path)
    
    # Read variant metadata from VCF (CHR, POS, EA, NEA, EAF, INFO)
    # Only loads variant-level data, not genotypes (genotypes loaded later when needed)
    df = _read_vcf_variants(
        vcf_path=vcf_path,
        region=region,
        vcf_chr_dict=vcf_chr_dict,
        tabix=tabix,
        mapper=mapper,
        log=log,
        verbose=verbose
    )
    
    # Validate that required columns are present
    for col in ("CHR", "POS", "EA", "NEA", "EAF"):
        if col not in df.columns:
            raise ValueError(f"VCF must provide column '{col}'")
    
    # Extract position array and count variants
    # DataFrame is already sorted by POS from _read_vcf_variants
    pos = df["POS"].to_numpy().astype(np.int64)
    m = df.height
    
    if m == 0:
        raise ValueError("No variants returned for the requested region.")
    
    log.write(f" -Loaded {m} variants", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # Step 0.5: Apply variant-level filters (MAF/EAF thresholds, exclude rare)
    # -------------------------------------------------------------------------
    # Filter variants based on minor allele frequency (MAF) and effect allele frequency (EAF)
    # This reduces the number of variants to process and ensures only valid variants are included
    df = _filter_variants(
        df=df,
        maf_min=maf_min,
        maf_max=maf_max,
        eaf_min=eaf_min,
        eaf_max=eaf_max,
        exclude_rare=exclude_rare,
        rare_maf_threshold=rare_maf_threshold,
        log=log,
        verbose=verbose
    )
    m = df.height
    
    # -------------------------------------------------------------------------
    # Step 0.6: Thin variants to simulate incomplete GWAS coverage
    # -------------------------------------------------------------------------
    # Real GWAS don't cover all sites - they may only have genotyped/imputed variants
    # at a subset of positions. This step randomly samples a fraction of variants
    # to simulate this incomplete coverage.
    if thin is not None:
        if thin <= 0.0 or thin >= 1.0:
            raise ValueError(f"thin must be between 0 and 1 (exclusive), got {thin}")
        n_keep = max(1, int(m * thin))  # Keep at least 1 variant
        log.write(f" -Thinning variants: keeping {n_keep} out of {m} variants (thin={thin:.4f})", verbose=verbose)
        # Randomly sample indices without replacement
        keep_indices = rng.choice(m, size=n_keep, replace=False)
        keep_indices = np.sort(keep_indices)  # Sort to maintain position order
        df = df[keep_indices]
        m = df.height
        log.write(f" -After thinning: {m} variants remaining", verbose=verbose)
    
    pos = df["POS"].to_numpy().astype(np.int64)
    
    # -------------------------------------------------------------------------
    # Step 1: Setup sample sizes (N and N_eff) for each variant
    # -------------------------------------------------------------------------
    # For quantitative traits: N = n (user-specified), N_eff = N
    # For binary traits: N = n_case + n_ctrl, N_eff = 4 / (1/n_case + 1/n_ctrl)
    # Also applies per-SNP N variation to simulate missingness/QC issues:
    # Randomly selects n_drop_rate fraction of SNPs and reduces their N
    N, N_eff = _setup_sample_sizes(
        m=m,
        trait=trait,
        n=n,
        n_case=n_case,
        n_ctrl=n_ctrl,
        n_drop_rate=n_drop_rate,
        n_drop_min=n_drop_min,
        rng=rng
    )
    
    # -------------------------------------------------------------------------
    # Step 2: Setup INFO scores and apply attenuation to N_eff
    # -------------------------------------------------------------------------
    # INFO score represents imputation quality (0-1 scale)
    # If INFO column exists in VCF and use_info=True: uses actual INFO scores
    # Otherwise: simulates INFO scores from Normal(info_mean, info_sd)
    info = _setup_info_scores(
        df=df,
        use_info=use_info,
        info_mean=info_mean,
        info_sd=info_sd,
        info_min=info_min,
        info_max=info_max,
        rng=rng
    )
    
    # Attenuate N_eff by INFO: N_eff <- N_eff * INFO
    # Lower INFO (poor imputation) => less information => larger SE => smaller Z-scores
    # Ensure N_eff is always positive and valid (minimum of 1.0) to prevent division issues
    N_eff = np.clip(N_eff * info, 1.0, None)
    
    # -------------------------------------------------------------------------
    # Step 3: Select causal variants and assign effect sizes
    # -------------------------------------------------------------------------
    # Selects which variants are causal based on mode (polygenic or sparse)
    # Assigns effect sizes with MAF-dependent architecture if alpha > 0
    # Effect size variance scales as [2p(1-p)]^(-alpha) where p is EAF
    if alpha > 0.0:
        log.write(f" -Sampling causal effect sizes: β_j^raw ~ N(0, [2p_j(1-p_j)]^(-α)) where α = {alpha:.4f}", verbose=verbose)
    else:
        log.write(f" -Sampling causal effect sizes: β_j^raw ~ N(0, effect_sd²) where effect_sd = {effect_sd:.4f}", verbose=verbose)
    beta_true, is_causal = _select_causal_variants(
        df=df,
        mode=mode,
        pi=pi,
        n_causal=n_causal,
        alpha=alpha,
        effect_sd=effect_sd,
        causal_maf_min=causal_maf_min,
        causal_maf_max=causal_maf_max,
        causal_eaf_min=causal_eaf_min,
        causal_eaf_max=causal_eaf_max,
        causal_idx=causal_idx,
        causal_beta=causal_beta,
        rng=rng,
        log=log,
        verbose=verbose
    )
    
    # -------------------------------------------------------------------------
    # Step 4: Simulate Z-scores using efficient X^T X approach
    # -------------------------------------------------------------------------
    # This method avoids explicit LD matrix construction (O(k²) memory) and 
    # Cholesky decomposition (O(k³) time) by using matrix-vector products.
    # Instead of computing R = (1/n) X^T X and then sampling z ~ N(μ, R),
    # we directly compute z = μ_causal + b_strat + ε where:
    #   - μ_causal = sqrt(N_eff) * (1/n) X^T (X β)  [mean term: causal effects spread via LD]
    #   - ε = λ * (1/sqrt(n)) X^T u, u ~ N(0, I_n)  [noise term: scaled by λ for cryptic relatedness]
    #   - b_strat = σ_strat * (1/sqrt(n)) X^T v, v ~ N(0, I_n)  [bias term: population stratification]
    # This is mathematically equivalent but much more efficient for large k.
    z = np.empty(m, dtype=np.float64)
    
    log.write(" -Simulating Z-scores using X^T X method", verbose=verbose)
    
    # Load genotypes from VCF and standardize them
    # Standardization: X_ij = (G_ij - μ_j) / σ_j where μ_j and σ_j are sample mean and SD
    # This ensures each variant has mean=0 and var=1 in the sample
    X, p_block, matched_indices = _get_standardized_genotypes_for_block(
        vcf_path=vcf_path,
        variant_df=df,
        start_idx=0,
        end_idx=m,
        region=region,
        vcf_chr_dict=vcf_chr_dict,
        tabix=tabix,
        mapper=mapper,
        log=log,
        verbose=verbose
    )
    
    if X.shape[0] == 0 or X.shape[1] == 0:
        # No matched variants found in VCF, set all Z-scores to 0
        z[:] = 0.0
    else:
        n_samples = X.shape[0]  # Number of samples in reference panel
        n_variants_matched = X.shape[1]  # Number of variants matched between sumstats and VCF
        
        # matched_indices are indices within [0, m) indicating which variants were matched
        # Extract beta and N_eff only for matched variants
        beta_matched = beta_true[matched_indices]
        N_eff_matched = N_eff[matched_indices]
        N_eff_matched = np.clip(N_eff_matched, 1.0, None)
        
        # Compute mean term: μ_causal = sqrt(N_eff) * (1/n) X^T (X β)
        log.write(" -Computing causal mean: μ_causal = sqrt(N_eff) * (1/n) X^T (X β)", verbose=verbose)
        # Step 1: Compute X β (effect of causal variants on phenotypes)
        X_beta = X @ beta_matched
        # Step 2: Compute (1/n) X^T (X β) ≈ R β (LD-weighted sum of effects)
        R_beta = (X.T @ X_beta) / float(n_samples)
        # Step 3: Scale by sqrt(N_eff) to account for sample size
        sqrt_N_eff = np.sqrt(N_eff_matched)
        mu_causal = sqrt_N_eff * R_beta
        
        # Generate standard LD-noise: ε_0 ~ N(0, R)
        # NOTE: Scaling must be 1/sqrt(n), NOT 1/n, so that Var(ε_0) = (1/n) X^T X = R
        # This ensures the noise has the correct correlation structure matching the LD matrix
        log.write(" -Generating standard LD-noise: ε_0 = (1/sqrt(n)) X^T u, u ~ N(0, I_n)", verbose=verbose)
        u = rng.standard_normal(n_samples)
        epsilon_0 = (X.T @ u) / np.sqrt(float(n_samples))
        
        # Apply cryptic relatedness / global inflation: ε = λ * ε_0
        # This scales the noise variance: Cov(ε) = λ * R
        # If λ > 1, creates genomic control inflation (QQ plot inflates)
        log.write(f" -Applying cryptic relatedness inflation: ε = λ * ε_0 (λ = {lambda_gc:.4f})", verbose=verbose)
        epsilon = lambda_gc * epsilon_0
        
        # Generate population stratification bias: b_strat ~ N(0, σ_strat^2 R)
        # This creates LD-correlated bias that mimics population structure
        # Creates genome-wide shifts and correlated "hills" of signal across LD blocks
        b_strat = np.zeros(n_variants_matched, dtype=np.float64)
        if sigma_strat > 0.0:
            log.write(f" -Generating population stratification bias: b_strat = σ_strat * (1/sqrt(n)) X^T v (σ_strat = {sigma_strat:.4f})", verbose=verbose)
            v = rng.standard_normal(n_samples)
            b_strat = sigma_strat * (X.T @ v) / np.sqrt(float(n_samples))
        else:
            log.write(" -Skipping population stratification bias (σ_strat = 0.0)", verbose=verbose)
        
        # Z-scores: z = μ_causal + b_strat + ε
        # μ_causal: signal from causal variants (spread via LD)
        # b_strat: LD-correlated bias from population stratification
        # ε: scaled noise (correlated due to LD, inflated by λ)
        log.write(" -Computing final Z-scores: z = μ_causal + b_strat + ε", verbose=verbose)
        z_matched = mu_causal + b_strat + epsilon
        
        # Validate: replace any NaN/Inf values with 0 (shouldn't happen, but safety check)
        z_matched = np.where(np.isfinite(z_matched), z_matched, 0.0).astype(np.float64)
        
        # Store Z-scores in global array (only for matched variants)
        z[matched_indices] = z_matched
        
        # Set unmatched variants to 0 (variants in sumstats but not found in VCF)
        all_indices = np.arange(m)
        unmatched = np.setdiff1d(all_indices, matched_indices)
        if len(unmatched) > 0:
            z[unmatched] = 0.0
    
    # -------------------------------------------------------------------------
    # Step 5: Convert Z-scores to summary statistics (SE, BETA, P, MLOG10P)
    # -------------------------------------------------------------------------
    # Convert Z-scores to other summary statistics:
    #   - SE = 1 / sqrt(N_eff)  [standard error]
    #   - BETA = Z * SE  [effect estimate, ensures Z = BETA/SE]
    #   - P = 2 * P(Z > |z|)  [two-sided p-value]
    #   - MLOG10P = -log10(P)  [negative log10 p-value, computed in log-space for precision]
    se, beta_hat, p, mlog10p = _convert_z_to_sumstats(
        z=z,
        N_eff=N_eff,
        log=log,
        verbose=verbose
    )
    
    # -------------------------------------------------------------------------
    # Step 6: Assemble output DataFrame (using Polars for performance)
    # -------------------------------------------------------------------------
    # Combine all variant metadata with simulated summary statistics
    # Use Polars for fast DataFrame operations, then convert to Pandas for Sumstats object
    # Create Series first with explicit types to ensure proper dtypes
    out = df.with_columns([
        pl.Series("INFO", info, dtype=pl.Float64),
        pl.Series("N", np.round(N).astype(int), dtype=pl.Int64),
        pl.Series("N_EFF", N_eff, dtype=pl.Float64),
        pl.Series("IS_CAUSAL", is_causal, dtype=pl.Boolean),
        pl.Series("BETA_TRUE", beta_true, dtype=pl.Float64),
        pl.Series("Z", z, dtype=pl.Float64),
        pl.Series("SE", se, dtype=pl.Float64),
        pl.Series("BETA", beta_hat, dtype=pl.Float64),
        pl.Series("P", p, dtype=pl.Float64),
        pl.Series("MLOG10P", mlog10p, dtype=pl.Float64),
    ])
    
    # Add binary trait specific columns (N_CASE and N_CONTROL)
    if trait == "binary":
        out = out.with_columns([
            pl.Series("N_CASE", np.full(m, int(n_case), dtype=int), dtype=pl.Int64),
            pl.Series("N_CONTROL", np.full(m, int(n_ctrl), dtype=int), dtype=pl.Int64),
        ])
    
    # Convert to Pandas for Sumstats object (required by GWASLab)
    out = out.to_pandas()
    
    log.write("Simulation completed successfully!", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # Step 7: Extract causal SNP IDs
    # -------------------------------------------------------------------------
    # Extract IDs of causal variants for return value
    # Uses VCF ID field if available, otherwise constructs as "CHR:POS:EA:NEA"
    causal_snp_ids = _extract_causal_snp_ids(
        df=df,
        is_causal=is_causal,
        log=log,
        verbose=verbose
    )
    
    # -------------------------------------------------------------------------
    # Step 8: Create and return Sumstats object
    # -------------------------------------------------------------------------
    # Wrap the DataFrame in a GWASLab Sumstats object with proper metadata
    # This enables all GWASLab functionality (QC, harmonization, visualization, etc.)
    log.write("Creating Sumstats object...", verbose=verbose)
    sumstats_obj = Sumstats(
        sumstats=out,
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        eaf="EAF",
        n="N",
        beta="BETA",
        se="SE",
        z="Z",
        p="P",
        info="INFO",
        neff="N_EFF",
        ncase="N_CASE" if trait == "binary" else None,
        ncontrol="N_CONTROL" if trait == "binary" else None,
        study=study,
        trait=trait_name,
        build=build,
        verbose=verbose
    )
    
    return sumstats_obj, causal_snp_ids


def simulate_sumstats_global(
    vcf_path: str,
    *,
    chromosomes: Optional[List[Union[int, str]]] = None,  # If None, use all chromosomes
    # trait and sample size
    trait: str = "quant",                    # "quant" or "binary"
    n: Optional[int] = 300_000,              # required if trait="quant"
    n_case: Optional[int] = None,            # required if trait="binary"
    n_ctrl: Optional[int] = None,            # required if trait="binary"
    # genetic architecture
    mode: str = "polygenic",                 # "polygenic" or "sparse"
    pi: float = 2e-3,                        # causal fraction if polygenic
    n_causal: Optional[int] = None,          # number of causals if sparse (per chromosome if None)
    h2: float = 0.1,                          # target heritability
    alpha: float = 0.2,                       # MAF dependence parameter (0 = MAF-independent)
    # variant filtering
    maf_min: float = 0.01,                    # minimum MAF for variants
    maf_max: float = 0.5,                    # maximum MAF for variants
    eaf_min: float = 0.0,                    # minimum EAF for variants
    eaf_max: float = 1.0,                    # maximum EAF for variants
    exclude_rare: bool = False,              # if True, exclude rare variants (MAF < 0.01)
    rare_maf_threshold: float = 0.01,        # MAF threshold for rare variants
    # causal variant selection filters
    causal_maf_min: Optional[float] = None,  # minimum MAF for causal variants
    causal_maf_max: Optional[float] = None,  # maximum MAF for causal variants
    causal_eaf_min: Optional[float] = None,  # minimum EAF for causal variants
    causal_eaf_max: Optional[float] = None,  # maximum EAF for causal variants
    # block-wise processing
    window_bp: int = 1_000_000,              # window size in basepairs for block-wise simulation (1Mb default)
    # realism knobs
    n_drop_rate: float = 0.15,               # fraction of SNPs with reduced N
    n_drop_min: float = 0.5,                 # minimum N multiplier for dropped SNPs
    use_info: bool = True,                   # if INFO column exists, use it
    info_mean: float = 0.95,                 # if simulating INFO
    info_sd: float = 0.05,
    info_min: float = 0.3,
    info_max: float = 1.0,
    # deterministic
    seed: int = 1,
    # VCF reading options
    vcf_chr_dict: Optional[dict] = None,
    tabix: Optional[bool] = None,
    verbose: bool = True,
    log: Optional[Log] = None,
    # Sumstats object options
    study: str = "Simulated_Study",
    trait_name: str = "Simulated_Trait",
    build: str = "19",
) -> Tuple[Sumstats, List[str]]:
    """
    Simulate genome-wide GWAS summary statistics using efficient reference-panel LD model.
    
    This function implements the model described in the documentation, using the efficient
    X^T X approach to avoid constructing explicit LD matrices. It processes chromosomes
    sequentially and applies global heritability calibration.
    
    Parameters
    ----------
    vcf_path : str
        Path to VCF/BCF file containing reference panel genotypes
        
    chromosomes : list, optional
        List of chromosomes to simulate (e.g., [1, 2, 3] or ["1", "2", "3"]).
        If None, simulates all chromosomes found in VCF.
        
    trait : str, default="quant"
        "quant" for quantitative trait or "binary" for case-control
        
    n : int, optional, default=300000
        Sample size for quantitative trait (required if trait="quant")
        
    n_case : int, optional
        Number of cases for binary trait (required if trait="binary")
        
    n_ctrl : int, optional
        Number of controls for binary trait (required if trait="binary")
        
    window_bp : int, default=1000000
        Window size in basepairs for block-wise simulation (1Mb default)
        
    mode : str, default="polygenic"
        "polygenic" or "sparse" genetic architecture
        
    pi : float, default=2e-3
        Causal fraction if mode="polygenic"
        
    n_causal : int, optional
        Number of causal variants if mode="sparse". If None and mode="sparse",
        uses a default based on chromosome size.
        
    h2 : float, default=0.1
        Target heritability (genetic variance on standardized genotype scale)
        
    alpha : float, default=0.0
        MAF dependence parameter. Effect size variance scales as [2p(1-p)]^(-alpha).
        alpha=0: MAF-independent effects
        alpha>0: rarer variants have larger effect sizes on average
        
    maf_min : float, default=0.01
        Minimum MAF for variants to include
        
    maf_max : float, default=0.5
        Maximum MAF for variants to include
        
    eaf_min : float, default=0.0
        Minimum EAF for variants to include
        
    eaf_max : float, default=1.0
        Maximum EAF for variants to include
        
    exclude_rare : bool, default=False
        If True, exclude rare variants (MAF < rare_maf_threshold)
        
    rare_maf_threshold : float, default=0.01
        MAF threshold for rare variants when exclude_rare=True
        
    causal_maf_min : float, optional
        Minimum MAF for variants eligible to be selected as causal
        
    causal_maf_max : float, optional
        Maximum MAF for variants eligible to be selected as causal
        
    causal_eaf_min : float, optional
        Minimum EAF for variants eligible to be selected as causal
        
    causal_eaf_max : float, optional
        Maximum EAF for variants eligible to be selected as causal
        
    n_drop_rate : float, default=0.15
        Fraction of SNPs with reduced sample size
        
    n_drop_min : float, default=0.5
        Minimum N multiplier for dropped SNPs
        
    use_info : bool, default=True
        If True and INFO column exists in VCF, use it for N_eff attenuation
        
    info_mean : float, default=0.95
        Mean INFO score if simulating
        
    info_sd : float, default=0.05
        Standard deviation of INFO score if simulating
        
    info_min : float, default=0.3
        Minimum INFO score
        
    info_max : float, default=1.0
        Maximum INFO score
        
    seed : int, default=1
        Random seed for reproducibility
        
    vcf_chr_dict : dict, optional
        Chromosome mapping dictionary for VCF. If None, auto-detected.
        
    tabix : bool, optional
        Whether to use tabix indexing. If None, auto-detected.
        
    verbose : bool, default=True
        Verbose output
        
    log : Log, optional
        Logging object. If None, creates a new one.
        
    study : str, default="Simulated_Study"
        Study name for Sumstats object
        
    trait_name : str, default="Simulated_Trait"
        Trait name for Sumstats object
        
    build : str, default="19"
        Genome build version for Sumstats object

    Returns
    -------
    tuple
        A tuple containing:
        - Sumstats: GWASLab Sumstats object with simulated summary statistics
        - List[str]: List of causal SNP IDs
        
    Notes
    -----
    The function implements the reference-panel LD model:
    
    1. Standardizes genotypes using sample statistics: X_ij = (G_ij - μ_j) / σ_j
       where μ_j and σ_j are the sample mean and SD computed from the genotype matrix.
       This ensures diag(R) = 1 exactly, preventing λGC inflation.
    2. Samples causal effects with MAF-dependent architecture:
       β_j^raw ~ N(0, [2p_j(1-p_j)]^(-alpha)) for causal j
    3. Computes raw genetic variance: V_g^raw = β^T R β ≈ (1/n) (Xβ)^T (Xβ)
    4. Rescales effects globally: β = β^raw * sqrt(h2 / V_g^raw)
    5. Generates Z-scores: z = sqrt(N_eff) * (1/n) X^T (Xβ) + (1/sqrt(n)) X^T u
       where u ~ N(0, I_n)
       NOTE: noise scaling must be 1/sqrt(n), NOT 1/n, so that Var(eps) = (1/n) X^T X = R
    
    This approach scales as O(nm) per block and avoids constructing explicit
    LD matrices, making it suitable for whole-genome simulation.
    """
    if log is None:
        log = Log()
    
    rng = np.random.default_rng(seed)
    
    log.write("Starting genome-wide GWAS summary statistics simulation...", verbose=verbose)
    log.write(f" -VCF path: {vcf_path}", verbose=verbose)
    log.write(f" -Target heritability: {h2:.4f}", verbose=verbose)
    log.write(f" -MAF dependence (alpha): {alpha:.4f}", verbose=verbose)
    
    # Initialize mapper
    mapper = ChromosomeMapper(log=log, verbose=verbose)
    mapper.detect_reference_format(vcf_path)
    
    # -------------------------------------------------------------------------
    # Phase 1: Load variants and select causals
    # -------------------------------------------------------------------------
    log.write("Phase 1: Loading variants and selecting causal variants...", verbose=verbose)
    
    # Determine chromosomes to process
    if chromosomes is None:
        # Auto-detect chromosomes: need to load all variants first
        log.write(" -Auto-detecting chromosomes: loading all variants from VCF...", verbose=verbose)
        df = _read_vcf_variants(
            vcf_path=vcf_path,
            region=None,  # Read all variants
            vcf_chr_dict=vcf_chr_dict,
            tabix=tabix,
            mapper=mapper,
            log=log,
            verbose=verbose
        )
        
        if df.height == 0:
            raise ValueError("No variants found in VCF")
        
        # Auto-detect chromosomes from loaded variants
        unique_chrs = df["CHR"].unique().to_list()
        # Sort chromosomes naturally
        def chr_sort_key(x):
            x_str = str(x).lstrip('chrCHR')
            try:
                return (0, int(x_str))  # Numeric chromosomes first
            except ValueError:
                return (1, x_str)  # Non-numeric (e.g., X, Y, MT) after
        
        chromosomes = sorted(unique_chrs, key=chr_sort_key)
        log.write(f" -Found {len(chromosomes)} chromosomes: {chromosomes[:5]}..." if len(chromosomes) > 5 else f" -Found {len(chromosomes)} chromosomes: {chromosomes}", verbose=verbose)
    else:
        # Load only specified chromosomes (more efficient)
        log.write(f" -Loading variants for {len(chromosomes)} specified chromosomes...", verbose=verbose)
        chr_str_list = [str(c) for c in chromosomes]
        
        # Load each chromosome separately and combine
        # For tabix-indexed VCFs, this is efficient as it only loads the specified chromosome
        all_dfs = []
        for chr_val in chromosomes:
            log.write(f"   -Loading chromosome {chr_val}...", verbose=verbose)
            # Load chromosome by using a very large end position
            # For tabix-indexed VCFs, this efficiently loads only the specified chromosome
            # For non-indexed VCFs, it will still scan but is better than loading everything
            # Using 500Mb as upper bound (largest human chromosome is ~250Mb)
            df_chr = _read_vcf_variants(
                vcf_path=vcf_path,
                region=(chr_val, 1, 500_000_000),  # Load whole chromosome (large upper bound)
                vcf_chr_dict=vcf_chr_dict,
                tabix=tabix,
                mapper=mapper,
                log=log,
                verbose=False  # Less verbose for each chromosome
            )
            
            # Filter to exact chromosome match (handle string/numeric differences)
            # This ensures we only get variants from the requested chromosome
            if df_chr.height > 0:
                df_chr = df_chr.filter(pl.col("CHR").cast(pl.Utf8) == str(chr_val))
                if df_chr.height > 0:
                    all_dfs.append(df_chr)
                    log.write(f"     -Loaded {df_chr.height} variants from chromosome {chr_val}", verbose=verbose)
        
        if len(all_dfs) == 0:
            raise ValueError(f"No variants found for specified chromosomes: {chromosomes}")
        
        # Combine all chromosomes using Polars (much faster)
        df = pl.concat(all_dfs)
        log.write(f" -Loaded variants from {len(chromosomes)} chromosomes", verbose=verbose)
    
    # Sort by CHR and POS (Polars is much faster)
    df = df.sort(["CHR", "POS"])
    
    m = df.height
    log.write(f" -Total variants loaded: {m}", verbose=verbose)
    
    # Apply variant filters
    df = _filter_variants(
        df=df,
        maf_min=maf_min,
        maf_max=maf_max,
        eaf_min=eaf_min,
        eaf_max=eaf_max,
        exclude_rare=exclude_rare,
        rare_maf_threshold=rare_maf_threshold,
        log=log,
        verbose=verbose
    )
    m = df.height
    
    # Setup N and N_eff
    N, N_eff = _setup_sample_sizes(
        m=m,
        trait=trait,
        n=n,
        n_case=n_case,
        n_ctrl=n_ctrl,
        n_drop_rate=n_drop_rate,
        n_drop_min=n_drop_min,
        rng=rng
    )
    
    # Setup INFO scores and apply attenuation
    info = _setup_info_scores(
        df=df,
        use_info=use_info,
        info_mean=info_mean,
        info_sd=info_sd,
        info_min=info_min,
        info_max=info_max,
        rng=rng
    )
    N_eff = np.clip(N_eff * info, 1.0, None)
    
    # Select causal variants (for global, we use raw effects that will be rescaled)
    # Note: For global simulation, we need to handle n_causal=None differently
    if mode == "sparse" and n_causal is None:
        # Default: ~10 per chromosome (rough estimate)
        n_causal_effective = max(1, int(round(m * 0.001)))
    else:
        n_causal_effective = n_causal if n_causal is not None else 3
    
    beta_raw, is_causal = _select_causal_variants(
        df=df,
        mode=mode,
        pi=pi,
        n_causal=n_causal_effective,
        alpha=alpha,
        effect_sd=1.0,  # Will be rescaled globally, so use 1.0 as placeholder
        causal_maf_min=causal_maf_min,
        causal_maf_max=causal_maf_max,
        causal_eaf_min=causal_eaf_min,
        causal_eaf_max=causal_eaf_max,
        causal_idx=None,  # Not supported in global mode
        causal_beta=None,  # Not supported in global mode
        rng=rng,
        log=log,
        verbose=verbose
    )
    
    # -------------------------------------------------------------------------
    # Phase 2: Compute raw genetic variance across all blocks
    # -------------------------------------------------------------------------
    log.write("Phase 2: Computing raw genetic variance for heritability calibration...", verbose=verbose)
    
    pos = df["POS"].to_numpy().astype(np.int64)
    chr_array = df["CHR"].to_numpy()
    
    # Group variants by chromosome
    chr_to_indices = {}
    for i, chr_val in enumerate(chr_array):
        if chr_val not in chr_to_indices:
            chr_to_indices[chr_val] = []
        chr_to_indices[chr_val].append(i)
    
    V_g_raw_total = 0.0
    n_blocks_processed = 0
    
    for chr_val in chromosomes:
        if chr_val not in chr_to_indices:
            continue
        
        chr_indices = np.array(chr_to_indices[chr_val])
        chr_pos = pos[chr_indices]
        
        # Create blocks for this chromosome
        blocks = _make_blocks_by_bp(chr_pos, int(window_bp))
    
    for block_idx, (start_i, end_i) in enumerate(blocks):
            # Map block indices to global indices
            # Note: end_i can be equal to len(chr_pos) (half-open interval [start, end))
            global_start = chr_indices[start_i]
            if end_i < len(chr_indices):
                global_end = chr_indices[end_i]
            else:
                # end_i equals length, meaning we want all variants from start_i to end of chromosome
                # Use one past the last variant index for this chromosome (exclusive end)
                # Cap at len(df) to avoid out-of-bounds
                global_end = min(chr_indices[-1] + 1 if len(chr_indices) > 0 else global_start, df.height)
            
            # Get standardized genotypes for this block
            X, p_block, matched_indices = _get_standardized_genotypes_for_block(
            vcf_path=vcf_path,
            variant_df=df,
                start_idx=global_start,
                end_idx=global_end,
                region=None,
            vcf_chr_dict=vcf_chr_dict,
            tabix=tabix,
            mapper=mapper,
            log=log,
                verbose=False
            )
            
            if X.shape[0] == 0 or X.shape[1] == 0:
                continue
            
            n_samples = X.shape[0]
            n_variants_block = X.shape[1]
            
            # Get beta_raw for matched variants
            # Map matched_indices (local to block, 0-based within block slice) to global indices
            # matched_indices are indices within the block slice [global_start, global_end)
            global_matched = global_start + matched_indices
            beta_block = beta_raw[global_matched]
            
            # Compute V_g^raw = (1/n) (Xβ)^T (Xβ)
            X_beta = X @ beta_block
            V_g_block = np.dot(X_beta, X_beta) / float(n_samples)
            
            V_g_raw_total += V_g_block
            n_blocks_processed += 1
    
    if n_blocks_processed == 0:
        raise ValueError("No blocks could be processed for genetic variance calculation")
    
    log.write(f" -Raw genetic variance: {V_g_raw_total:.6f}", verbose=verbose)
    log.write(f" -Processed {n_blocks_processed} blocks", verbose=verbose)
    
    # Compute rescaling factor
    if V_g_raw_total <= 0:
        log.warning(" -Raw genetic variance is non-positive, using identity scaling", verbose=verbose)
        scale_factor = 1.0
    else:
        scale_factor = np.sqrt(h2 / V_g_raw_total)
    
    log.write(f" -Rescaling factor: {scale_factor:.6f}", verbose=verbose)
    
    # Rescale effects
    beta_true = beta_raw * scale_factor
    
    # -------------------------------------------------------------------------
    # Phase 3: Generate Z-scores using efficient X^T X approach
    # -------------------------------------------------------------------------
    log.write("Phase 3: Generating Z-scores...", verbose=verbose)
    
    z = np.zeros(m, dtype=float)
    
    for chr_val in chromosomes:
        if chr_val not in chr_to_indices:
            continue
        
        chr_indices = np.array(chr_to_indices[chr_val])
        chr_pos = pos[chr_indices]
        
        blocks = _make_blocks_by_bp(chr_pos, int(window_bp))
        
        log.write(f" -Processing chromosome {chr_val}: {len(blocks)} blocks", verbose=verbose)
        
        for block_idx, (start_i, end_i) in enumerate(blocks):
            # Map block indices to global indices
            # Note: end_i can be equal to len(chr_pos) (half-open interval [start, end))
            global_start = chr_indices[start_i]
            if end_i < len(chr_indices):
                global_end = chr_indices[end_i]
            else:
                # end_i equals length, meaning we want all variants from start_i to end of chromosome
                # Use one past the last variant index for this chromosome (exclusive end)
                # Cap at len(df) to avoid out-of-bounds
                global_end = min(chr_indices[-1] + 1 if len(chr_indices) > 0 else global_start, df.height)
            
            # Get standardized genotypes
            X, p_block, matched_indices = _get_standardized_genotypes_for_block(
                vcf_path=vcf_path,
                variant_df=df,
                start_idx=global_start,
                end_idx=global_end,
                region=None,
                vcf_chr_dict=vcf_chr_dict,
                tabix=tabix,
                mapper=mapper,
                log=log,
                verbose=False
            )
            
            if X.shape[0] == 0 or X.shape[1] == 0:
                continue
            
            n_samples = X.shape[0]
            n_variants_block = X.shape[1]
            
            # Map matched indices to global
            # matched_indices are indices within the block slice [global_start, global_end)
            global_matched = global_start + matched_indices
            
            # Get rescaled effects and N_eff for this block
            beta_block = beta_true[global_matched]
            N_eff_block = N_eff[global_matched]
            N_eff_block = np.clip(N_eff_block, 1.0, None)
            
            # Compute mean: μ = sqrt(N_eff) * (1/n) X^T (Xβ)
            X_beta = X @ beta_block
            R_beta = (X.T @ X_beta) / float(n_samples)
            sqrt_N_eff = np.sqrt(N_eff_block)
            mu = sqrt_N_eff * R_beta
            
            # Generate noise: ε = (1/sqrt(n)) X^T u, where u ~ N(0, I_n)
            # NOTE: scaling must be 1/sqrt(n), NOT 1/n, so that Var(ε) = (1/n) X^T X = R
            u = rng.standard_normal(n_samples)
            epsilon = (X.T @ u) / np.sqrt(float(n_samples))
            
            # Z-scores: z = μ + ε
            z_block = mu + epsilon
            
            # Validate
            z_block = np.where(np.isfinite(z_block), z_block, 0.0)
            
            # Store in global array
            z[global_matched] = z_block
    
    # -------------------------------------------------------------------------
    # Phase 4: Convert to summary statistics
    # -------------------------------------------------------------------------
    log.write("Phase 4: Converting to summary statistics...", verbose=verbose)
    se, beta_hat, p, mlog10p = _convert_z_to_sumstats(
        z=z,
        N_eff=N_eff,
        log=log,
        verbose=verbose
    )
    
    # -------------------------------------------------------------------------
    # Phase 5: Assemble output (using Polars)
    # -------------------------------------------------------------------------
    log.write("Phase 5: Assembling output...", verbose=verbose)
    
    # Add all columns using Polars (much faster)
    # Create Series with explicit dtype to avoid casting issues
    out = df.with_columns([
        pl.Series("INFO", info, dtype=pl.Float64),
        pl.Series("N", np.round(N).astype(int), dtype=pl.Int64),
        pl.Series("N_EFF", N_eff, dtype=pl.Float64),
        pl.Series("IS_CAUSAL", is_causal, dtype=pl.Boolean),
        pl.Series("BETA_TRUE", beta_true, dtype=pl.Float64),
        pl.Series("Z", z, dtype=pl.Float64),
        pl.Series("SE", se, dtype=pl.Float64),
        pl.Series("BETA", beta_hat, dtype=pl.Float64),
        pl.Series("P", p, dtype=pl.Float64),
        pl.Series("MLOG10P", mlog10p, dtype=pl.Float64),
    ])
    
    if trait == "binary":
        out = out.with_columns([
            pl.Series("N_CASE", np.full(m, int(n_case), dtype=int), dtype=pl.Int64),
            pl.Series("N_CONTROL", np.full(m, int(n_ctrl), dtype=int), dtype=pl.Int64),
        ])
    
    # Convert to Pandas for Sumstats object (required)
    out = out.to_pandas()
    
    # Extract causal SNP IDs
    causal_snp_ids = _extract_causal_snp_ids(
        df=df,
        is_causal=is_causal,
        log=log,
        verbose=verbose
    )
    
    # Create Sumstats object
    sumstats_obj = Sumstats(
        sumstats=out,
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        eaf="EAF",
        n="N",
        beta="BETA",
        se="SE",
        z="Z",
        p="P",
        info="INFO",
        neff="N_EFF",
        ncase="N_CASE" if trait == "binary" else None,
        ncontrol="N_CONTROL" if trait == "binary" else None,
        study=study,
        trait=trait_name,
        build=build,
        verbose=verbose
    )
    
    log.write("Genome-wide simulation completed successfully!", verbose=verbose)
    
    return sumstats_obj, causal_snp_ids
