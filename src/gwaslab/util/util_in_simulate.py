from __future__ import annotations
import numpy as np
import pandas as pd
from math import sqrt
from typing import Callable, Optional, Tuple, List, Any, Union
from pysam import VariantFile
from allel import read_vcf, GenotypeArray, rogers_huff_r_between
from scipy.special import erfc
import scipy.stats as ss
from gwaslab.g_Sumstats import Sumstats
from gwaslab.io.io_vcf import auto_check_vcf_chr_dict
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log
from shutil import which


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


def _ensure_psd_corr(R: np.ndarray, jitter: float = 1e-8, max_tries: int = 6) -> np.ndarray:
    """
    Ensure R is symmetric, correlation-like (diag=1), and numerically PSD so that
    Cholesky decomposition works.

    Why:
    - LD matrices estimated from finite samples/rounding can be slightly non-PSD.
    - MVN sampling uses Cholesky: R = L L^T (requires PSD).
    """
    R = (R + R.T) / 2.0
    np.fill_diagonal(R, 1.0)

    for t in range(max_tries):
        eps = jitter * (10 ** t)
        try:
            np.linalg.cholesky(R + np.eye(R.shape[0]) * eps)
            return R + np.eye(R.shape[0]) * eps
        except np.linalg.LinAlgError:
            continue

    # Fallback: eigenvalue clipping
    w, V = np.linalg.eigh(R)
    w = np.clip(w, 1e-10, None)
    R2 = (V * w) @ V.T
    R2 = (R2 + R2.T) / 2.0
    np.fill_diagonal(R2, 1.0)
    return R2


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
) -> pd.DataFrame:
    """
    Read variants from VCF file and return DataFrame with CHR, POS, EA, NEA, EAF, INFO.
    
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
    pd.DataFrame
        DataFrame with columns: CHR, POS, EA, NEA, EAF, INFO (if available)
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
    if region is not None:
        chr_, start, end = region
        # Convert chromosome to reference format
        region_chr_ref = mapper.sumstats_to_reference(chr_, reference_file=vcf_path, as_string=True)
        region_str = f"{region_chr_ref}:{start}-{end}"
        log.write(f" -Loading VCF region: {region_str}", verbose=verbose)
        vcf_data = read_vcf(vcf_path, region=region_str, tabix=tabix)
    else:
        log.write(" -Loading all variants from VCF", verbose=verbose)
        vcf_data = read_vcf(vcf_path, tabix=tabix)
    
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
    
    # Extract ID (SNP ID) from VCF if available
    snpid_data = None
    if "variants/ID" in vcf_data:
        snpid_data = vcf_data["variants/ID"]
        # Handle bytes encoding
        if isinstance(snpid_data[0], bytes):
            snpid_data = np.array([s.decode() if isinstance(s, bytes) else str(s) for s in snpid_data], dtype=object)
        else:
            snpid_data = np.array([str(s) if s != b'.' and s != '.' else None for s in snpid_data], dtype=object)
        # Replace missing IDs ('.' or None) with None
        snpid_data = np.where((snpid_data == '.') | (snpid_data == 'None') | (snpid_data == ''), None, snpid_data)
    
    # Extract REF and ALT alleles (vectorized)
    ref_data = vcf_data["variants/REF"]
    if isinstance(ref_data[0], bytes):
        # Vectorized bytes decoding
        ref_data = np.array([r.decode() if isinstance(r, bytes) else str(r) for r in ref_data], dtype=object)
    else:
        ref_data = np.array([str(r) for r in ref_data], dtype=object)
    
    alt_data = vcf_data["variants/ALT"]
    # ALT can be 2D array (multiple ALT alleles per variant)
    if alt_data.ndim == 2:
        # Take first ALT allele (vectorized)
        alt_data = alt_data[:, 0]
    if isinstance(alt_data[0], bytes):
        # Vectorized bytes decoding
        alt_data = np.array([a.decode() if isinstance(a, bytes) else str(a) for a in alt_data], dtype=object)
    else:
        alt_data = np.array([str(a) for a in alt_data], dtype=object)
    
    # Extract EAF from INFO field (AF) if available, otherwise compute from genotypes
    eaf_data = np.full(n_variants, np.nan, dtype=float)
    if "variants/INFO" in vcf_data and "AF" in vcf_data["variants/INFO"].dtype.names:
        af_info = vcf_data["variants/INFO"]["AF"]
        if af_info.ndim == 2:
            eaf_data = af_info[:, 0].astype(float)
        else:
            eaf_data = af_info.astype(float)
    
    # If EAF still missing, try to compute from genotypes (vectorized)
    if np.any(np.isnan(eaf_data)) and "calldata/GT" in vcf_data:
        # Compute AF from genotypes for missing values (vectorized)
        gt_array = GenotypeArray(vcf_data["calldata/GT"])
        n_alt = gt_array.to_n_alt()
        # Vectorized: compute allele frequency: mean of alt allele counts / 2
        computed_eaf = np.mean(n_alt, axis=1) / 2.0
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
    
    # Create DataFrame
    df = pd.DataFrame({
        "CHR": chr_data,
        "POS": pos_data,
        "EA": alt_data,  # ALT is typically the effect allele
        "NEA": ref_data,  # REF is typically the non-effect allele
        "EAF": eaf_data,
    })
    
    # Add SNPID if available from VCF
    if snpid_data is not None:
        df["SNPID"] = snpid_data
    
    # Add INFO if available
    if not np.all(np.isnan(info_data)):
        df["INFO"] = info_data
    
    # Sort by position
    df = df.sort_values(["CHR", "POS"]).reset_index(drop=True)
    
    log.write(f" -Loaded {len(df)} variants", verbose=verbose)
    return df


def _get_ld_matrix_for_block(
    vcf_path: str,
    variant_df: pd.DataFrame,
    start_idx: int,
    end_idx: int,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    vcf_chr_dict: Optional[dict] = None,
    tabix: Optional[bool] = None,
    mapper: Optional[ChromosomeMapper] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> np.ndarray:
    """
    Get LD correlation matrix (r, not r^2) for a block of variants.
    
    Parameters
    ----------
    vcf_path : str
        Path to VCF file
    variant_df : pd.DataFrame
        DataFrame with variant information (must have CHR, POS, EA, NEA)
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
    np.ndarray
        LD correlation matrix (r values, signed)
    """
    if log is None:
        log = Log()
    
    if tabix is None:
        tabix = which("tabix")
    
    # Get or create mapper
    if mapper is None:
        mapper = ChromosomeMapper(log=log, verbose=verbose)
        mapper.detect_reference_format(vcf_path)
    
    # Get block variants
    block_df = variant_df.iloc[start_idx:end_idx].copy()
    if len(block_df) == 0:
        return np.array([]).reshape(0, 0)
    
    # Determine region for this block
    chr_ = block_df["CHR"].iloc[0]
    start_pos = block_df["POS"].min()
    end_pos = block_df["POS"].max()
    
    # Convert chromosome to reference format
    region_chr_ref = mapper.sumstats_to_reference(chr_, reference_file=vcf_path, as_string=True)
    region_str = f"{region_chr_ref}:{start_pos}-{end_pos}"
    
    # Read VCF data for this block
    vcf_data = read_vcf(vcf_path, region=region_str, tabix=tabix)
    
    if vcf_data is None or len(vcf_data["variants/POS"]) == 0:
        log.warning(f" -No VCF data for block {start_idx}:{end_idx}, returning identity matrix", verbose=verbose)
        k = end_idx - start_idx
        return np.eye(k)
    
    # Match variants in block to VCF data (vectorized)
    vcf_pos = vcf_data["variants/POS"].astype(np.int64)
    vcf_ref = vcf_data["variants/REF"]
    if isinstance(vcf_ref[0], bytes):
        vcf_ref = np.array([r.decode() if isinstance(r, bytes) else str(r) for r in vcf_ref], dtype=object)
    else:
        vcf_ref = np.array([str(r) for r in vcf_ref], dtype=object)
    
    vcf_alt = vcf_data["variants/ALT"]
    if vcf_alt.ndim == 2:
        vcf_alt = vcf_alt[:, 0]
    if isinstance(vcf_alt[0], bytes):
        vcf_alt = np.array([a.decode() if isinstance(a, bytes) else str(a) for a in vcf_alt], dtype=object)
    else:
        vcf_alt = np.array([str(a) for a in vcf_alt], dtype=object)
    
    # Vectorized variant matching: match block variants to VCF variants
    # Following the approach in io_vcf.py for robust matching
    block_pos = block_df["POS"].values.astype(np.int64)
    block_nea = block_df["NEA"].values
    block_ea = block_df["EA"].values
    
    # Create position lookup dictionary for O(1) lookup
    pos_to_indices = {}
    for idx, pos in enumerate(vcf_pos):
        if pos not in pos_to_indices:
            pos_to_indices[pos] = []
        pos_to_indices[pos].append(idx)
    
    # Match variants: similar to match_variant function in io_vcf.py
    matched_indices = []
    matched_block_indices = []  # Track which block indices were matched
    
    for block_idx, (pos, nea, ea) in enumerate(zip(block_pos, block_nea, block_ea)):
        if pos not in pos_to_indices:
            continue
        
        # Get all VCF variants at this position
        candidate_indices = pos_to_indices[pos]
        
        if len(candidate_indices) == 1:
            # Single match - check if alleles match
            vcf_idx = candidate_indices[0]
            # Check forward match: NEA==REF and EA==ALT
            if vcf_ref[vcf_idx] == nea and vcf_alt[vcf_idx] == ea:
                matched_indices.append(vcf_idx)
                matched_block_indices.append(block_idx)
            # Check reverse match: NEA==ALT and EA==REF (flipped)
            elif vcf_ref[vcf_idx] == ea and vcf_alt[vcf_idx] == nea:
                matched_indices.append(vcf_idx)
                matched_block_indices.append(block_idx)
        else:
            # Multiple matches at same position - check each one
            matched = False
            for vcf_idx in candidate_indices:
                # Check forward match
                if vcf_ref[vcf_idx] == nea:
                    # Check if EA is in ALT (can be multi-allelic)
                    if isinstance(vcf_alt[vcf_idx], (list, np.ndarray)):
                        if ea in vcf_alt[vcf_idx]:
                            matched_indices.append(vcf_idx)
                            matched_block_indices.append(block_idx)
                            matched = True
                            break
                    elif vcf_alt[vcf_idx] == ea:
                        matched_indices.append(vcf_idx)
                        matched_block_indices.append(block_idx)
                        matched = True
                        break
                # Check reverse match (flipped)
                elif vcf_ref[vcf_idx] == ea:
                    if isinstance(vcf_alt[vcf_idx], (list, np.ndarray)):
                        if nea in vcf_alt[vcf_idx]:
                            matched_indices.append(vcf_idx)
                            matched_block_indices.append(block_idx)
                            matched = True
                            break
                    elif vcf_alt[vcf_idx] == nea:
                        matched_indices.append(vcf_idx)
                        matched_block_indices.append(block_idx)
                        matched = True
                        break
                # Check if NEA is in ALT and EA is REF (reverse)
                elif isinstance(vcf_alt[vcf_idx], (list, np.ndarray)):
                    if nea in vcf_alt[vcf_idx] and vcf_ref[vcf_idx] == ea:
                        matched_indices.append(vcf_idx)
                        matched_block_indices.append(block_idx)
                        matched = True
                        break
    
    if len(matched_indices) < len(block_df):
        log.warning(f" -Only {len(matched_indices)}/{len(block_df)} variants matched in VCF for LD calculation", verbose=verbose)
    
    if len(matched_indices) == 0:
        log.warning(f" -No variants matched in VCF for block {start_idx}:{end_idx}, returning identity matrix", verbose=verbose)
        k = end_idx - start_idx
        return np.eye(k)
    
    # Convert to numpy arrays for indexing
    matched_indices = np.array(matched_indices, dtype=int)
    matched_block_indices = np.array(matched_block_indices, dtype=int)
    
    # Get genotypes for matched variants
    if "calldata/GT" not in vcf_data:
        log.warning(" -No genotype data in VCF, returning identity matrix", verbose=verbose)
        k = end_idx - start_idx
        return np.eye(k)
    
    # Get genotype array and convert to alternate allele counts
    gt_array = GenotypeArray(vcf_data["calldata/GT"][matched_indices])
    n_alt = gt_array.to_n_alt()
    
    # Filter out monomorphic variants (no variation) - these cause NaN in LD calculation
    # A variant is monomorphic if all samples have the same genotype
    n_variants = n_alt.shape[0]
    if n_variants == 0:
        k = end_idx - start_idx
        return np.eye(k)
    
    # Check for monomorphic variants: std == 0 means no variation
    variant_std = np.std(n_alt, axis=1)
    polymorphic_mask = variant_std > 1e-10  # Small threshold to handle numerical precision
    
    n_polymorphic = polymorphic_mask.sum()
    n_monomorphic = n_variants - n_polymorphic
    
    if n_monomorphic > 0:
        log.write(f" -Filtering out {n_monomorphic} monomorphic variants (no variation)", verbose=verbose)
    
    if n_polymorphic == 0:
        log.warning(f" -All {n_variants} variants are monomorphic, returning identity matrix", verbose=verbose)
        k = end_idx - start_idx
        return np.eye(k)
    
    # Filter to only polymorphic variants
    n_alt_poly = n_alt[polymorphic_mask, :]
    matched_indices_poly = matched_indices[polymorphic_mask]
    
    # Check for missing data (NaN in genotypes)
    missing_mask = np.isnan(n_alt_poly).any(axis=1)
    n_with_missing = missing_mask.sum()
    if n_with_missing > 0:
        log.write(f" -Found {n_with_missing} variants with missing genotypes", verbose=verbose)
        # For variants with missing data, we'll let rogers_huff_r handle it
        # (it typically handles missing data by using pairwise complete observations)
    
    # Calculate LD correlation matrix (r, signed) for polymorphic variants
    log.write(f" -Computing LD matrix for {n_polymorphic} polymorphic variants...", verbose=verbose)
    r_matrix_poly = rogers_huff_r_between(n_alt_poly, n_alt_poly)
    
    # Validate r_matrix - check for NaN/Inf values
    if np.any(~np.isfinite(r_matrix_poly)):
        invalid_count = np.sum(~np.isfinite(r_matrix_poly))
        total_elements = r_matrix_poly.size
        invalid_fraction = invalid_count / total_elements if total_elements > 0 else 0.0
        
        log.warning(f" -LD matrix has {invalid_count}/{total_elements} invalid values (NaN/Inf, {invalid_fraction*100:.1f}%)", verbose=verbose)
        
        # Replace invalid values: NaN/Inf -> 0 (no correlation)
        # This is more conservative than identity - only diagonal should be 1
        r_matrix_poly = np.where(np.isfinite(r_matrix_poly), r_matrix_poly, 0.0)
        # Ensure diagonal is 1 (perfect correlation with self)
        np.fill_diagonal(r_matrix_poly, 1.0)
    
    # If we filtered out monomorphic variants, we need to expand the matrix back
    if n_polymorphic < n_variants:
        # Create full-size matrix with identity for monomorphic variants
        r_matrix = np.eye(n_variants)
        # Map polymorphic variants back to their original positions
        poly_indices_in_full = np.where(polymorphic_mask)[0]
        for i, orig_idx_i in enumerate(poly_indices_in_full):
            for j, orig_idx_j in enumerate(poly_indices_in_full):
                r_matrix[orig_idx_i, orig_idx_j] = r_matrix_poly[i, j]
    else:
        r_matrix = r_matrix_poly
    
    # Map LD matrix back to block variant positions
    # If we have fewer matched variants than block size, we need to map correctly
    k = len(block_df)  # Block size
    if len(matched_indices) < k:
        # Create full-size matrix with identity (no LD for unmatched variants)
        R_full = np.eye(k)
        # Map LD values to correct positions based on matched_block_indices
        # matched_block_indices tells us which block positions have LD data
        if len(matched_block_indices) > 0:
            # Create mapping: block position -> matched position in r_matrix
            for i, block_i in enumerate(matched_block_indices):
                for j, block_j in enumerate(matched_block_indices):
                    if block_i < k and block_j < k:
                        R_full[block_i, block_j] = r_matrix[i, j]
        return R_full
    elif len(matched_indices) == k:
        # Perfect match - r_matrix should already be k x k
        if r_matrix.shape[0] != k:
            # If shape doesn't match (due to filtering), pad/truncate
            R_full = np.eye(k)
            min_size = min(r_matrix.shape[0], k)
            R_full[:min_size, :min_size] = r_matrix[:min_size, :min_size]
            return R_full
        return r_matrix
    else:
        # More matched than block size (shouldn't happen, but handle it)
        return r_matrix[:k, :k]


# =============================================================================
# Main function
# =============================================================================

def simulate_sumstats(
    vcf_path: str,
    *,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    # trait and sample size
    trait: str = "quant",                    # "quant" or "binary"
    n: Optional[int] = 300_000,              # required if trait="quant"
    n_case: Optional[int] = None,            # required if trait="binary"
    n_ctrl: Optional[int] = None,            # required if trait="binary"
    # block simulation settings
    window_bp: int = 1_000_000,              # simulate in 1Mb windows
    # genetic architecture
    mode: str = "polygenic",                 # "polygenic" or "sparse"
    pi: float = 2e-3,                        # causal fraction if polygenic
    n_causal: int =3,                      # number of causals if sparse
    effect_sd: float = 0.05,                 # SD of causal effects (standardized units)
    # variant filtering
    maf_min: float = 0.01,                    # minimum MAF for variants (0.0 = no filter)
    maf_max: float = 0.5,                    # maximum MAF for variants (0.5 = no filter)
    eaf_min: float = 0.0,                    # minimum EAF for variants (0.0 = no filter)
    eaf_max: float = 1.0,                    # maximum EAF for variants (1.0 = no filter)
    exclude_rare: bool = False,              # if True, exclude rare variants (MAF < 0.01)
    rare_maf_threshold: float = 0.01,        # MAF threshold for rare variants
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
    ld_shrink: float = 0.02,                 # shrink R toward I (panel mismatch)
    jitter: float = 1e-8,                    # numerical stability
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
    Simulate GWAS summary statistics (BETA/SE/Z/P) for a region using an LD panel from VCF.

    Parameters
    ----------
    vcf_path : str
        Path to VCF/BCF file containing reference panel genotypes
        
    region : tuple, optional
        Optional tuple (chr, start, end). If provided, only variants in this region are used.
        Example: region=("1", 100_000, 2_000_000)
        
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
        
    n_causal : int, default=10
        Number of causal variants if mode="sparse"
        
    effect_sd : float, default=0.02
        Standard deviation of causal effects (standardized units)
        
    n_drop_rate : float, default=0.15
        Fraction of SNPs with reduced sample size
        
    n_drop_min : float, default=0.5
        Minimum N multiplier for dropped SNPs
        
    use_info : bool, default=True
        If True and INFO column exists in VCF, use it for N_eff attenuation
        
    info_mean : float, default=0.95
        Mean INFO score if simulating (when INFO not in VCF)
        
    info_sd : float, default=0.05
        Standard deviation of INFO score if simulating
        
    info_min : float, default=0.3
        Minimum INFO score
        
    info_max : float, default=1.0
        Maximum INFO score
        
    ld_shrink : float, default=0.02
        Shrinkage factor for LD matrix toward identity (panel mismatch)
        
    jitter : float, default=1e-8
        Numerical stability jitter for PSD correction
        
    seed : int, default=1
        Random seed for reproducibility
        
    causal_idx : np.ndarray, optional
        Explicit indices of causal variants (0-based, relative to variant order)
        
    causal_beta : np.ndarray, optional
        Explicit causal effect sizes (must match causal_idx length)
        
    maf_min : float, default=0.0
        Minimum MAF for variants to include in simulation (0.0 = no filter)
        
    maf_max : float, default=0.5
        Maximum MAF for variants to include in simulation (0.5 = no filter)
        
    eaf_min : float, default=0.0
        Minimum EAF for variants to include in simulation (0.0 = no filter)
        
    eaf_max : float, default=1.0
        Maximum EAF for variants to include in simulation (1.0 = no filter)
        
    exclude_rare : bool, default=False
        If True, exclude rare variants (MAF < rare_maf_threshold)
        
    rare_maf_threshold : float, default=0.01
        MAF threshold for rare variants when exclude_rare=True
        
    causal_maf_min : float, optional
        Minimum MAF for variants eligible to be selected as causal (None = no filter)
        
    causal_maf_max : float, optional
        Maximum MAF for variants eligible to be selected as causal (None = no filter)
        
    causal_eaf_min : float, optional
        Minimum EAF for variants eligible to be selected as causal (None = no filter)
        
    causal_eaf_max : float, optional
        Maximum EAF for variants eligible to be selected as causal (None = no filter)
        
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
        
    build : str, default="38"
        Genome build version for Sumstats object

    Returns
    -------
    tuple
        A tuple containing:
        - Sumstats: GWASLab Sumstats object with simulated summary statistics
        - List[str]: List of causal SNP IDs. Uses VCF ID field if available,
          otherwise constructed as "CHR:POS:EA:NEA"
        
    Notes
    -----
    The function implements a comprehensive workflow for simulating GWAS summary statistics:
    
    1. **Variant Loading and Filtering**:
       - Reads variants from VCF/BCF file (optionally restricted to a genomic region)
       - Extracts required columns: CHR, POS, EA (effect allele), NEA (non-effect allele), EAF (effect allele frequency)
       - Applies variant-level filters:
         * MAF filters: excludes variants outside [maf_min, maf_max] range
         * EAF filters: excludes variants outside [eaf_min, eaf_max] range
         * Rare variant exclusion: if exclude_rare=True, removes variants with MAF < rare_maf_threshold
       - Computes MAF from EAF: MAF = min(EAF, 1-EAF)
    
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
         * Lower INFO (imputation uncertainty) → larger SE → smaller Z-scores
    
    4. **Causal Variant Selection**:
       - Applies causal variant filters (if specified):
         * causal_maf_min/max: restricts candidates by MAF
         * causal_eaf_min/max: restricts candidates by EAF
       - If explicit causals provided (causal_idx, causal_beta): uses those directly
       - Otherwise, selects causals based on mode:
         * "sparse": randomly selects n_causal variants from candidates
         * "polygenic": randomly selects m * pi variants from candidates
       - Assigns causal effect sizes: beta_true ~ N(0, effect_sd^2) for selected variants
    
    5. **Block-wise LD Matrix Computation**:
       - Partitions variants into blocks of ~window_bp basepairs (default 1Mb)
       - For each block:
         * Matches variants between summary stats and VCF by position and alleles
         * Extracts genotype data from VCF for matched variants
         * Computes signed LD correlation matrix R using Rogers-Huff correlation
         * Validates R: replaces NaN/Inf with identity matrix (no LD)
         * Applies LD shrinkage: R <- (1 - ld_shrink) * R + ld_shrink * I
           (simulates panel/ancestry mismatch)
         * Ensures R is numerically positive semi-definite (PSD) for Cholesky decomposition:
           - Symmetrizes: R <- (R + R^T) / 2
           - Sets diagonal to 1.0
           - Adds jitter if needed: R <- R + eps * I
           - If still fails, uses eigenvalue clipping
    
    6. **Z-score Simulation via Multivariate Normal (MVN)**:
       - For each block with k variants:
         * Computes mean vector: mu = sqrt(N_eff) * (R @ beta_true)
           - R @ beta_true spreads causal effects to tagged SNPs via LD
           - sqrt(N_eff) scales with sample size (Z ~ sqrt(N))
         * Samples MVN: z ~ N(mu, R)
           - Uses Cholesky decomposition: R = L L^T
           - Generates epsilon ~ N(0, I)
           - Computes: z = mu + L @ epsilon
         * Validates Z-scores: replaces NaN/Inf with 0.0 (no effect)
    
    7. **Summary Statistics Conversion**:
       - Standard error: SE = 1 / sqrt(N_eff)
       - Effect estimate: BETA = Z * SE (ensures Z = BETA/SE)
       - MLOG10P: -log10(P) computed directly from Z using extreme mode (log-space):
         * Uses log-space for numerical precision: log_pvalue = log(2) + norm.logsf(|Z|)
         * Converts to -log10: mlog10p = -log_pvalue / log(10)
         * Handles extreme P-values without numerical underflow (same as util_in_fill_data.py)
         * This is the primary P-value representation (no intermediate P calculation)
       - P-value: P = 10^(-MLOG10P) (computed from MLOG10P for compatibility)
         * Derived from MLOG10P, not calculated directly from Z
         * Allows natural P-value range without clipping
    
    8. **Output Assembly**:
       - Creates DataFrame with all original VCF columns plus:
         * SNPID: variant ID from VCF (if available in VCF ID field)
         * INFO, N, N_EFF: sample size and effective sample size
         * IS_CAUSAL, BETA_TRUE: causal variant indicators and true effects
         * Z, SE, BETA, P, MLOG10P: simulated summary statistics
       - Wraps in GWASLab Sumstats object with metadata:
         * study, trait_name, build (genome build version)
       - Returns tuple: (Sumstats object, list of causal SNP IDs)
         * Causal SNP IDs use VCF ID if available, otherwise constructed as "CHR:POS:EA:NEA"
    
    Mathematical Model
    --------------------
    The simulation approximates marginal Z-scores as multivariate normal:
    
    z ~ N(mu, R)
    
    where:
    - mu = sqrt(N_eff) * (R @ beta_true)
      * R @ beta_true: causal effects spread to tagged SNPs via LD
      * sqrt(N_eff): Z-scores scale with square root of effective sample size
    - R: signed LD correlation matrix (computed from VCF genotypes)
    - beta_true: true causal effect sizes (mostly zeros, non-zero for causal variants)
    
    The model accounts for:
    - Linkage disequilibrium (LD) between variants
    - Imputation quality (via INFO attenuation of N_eff)
    - Sample size variation across SNPs
    - Panel/ancestry mismatch (via LD shrinkage)
    
    Performance
    -----------
    - Vectorized operations for maximum performance
    - Block-wise processing to handle large regions efficiently
    - LD matrices computed on-the-fly from VCF (no pre-computed LD required)
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
    # 0) Load variants from VCF
    # -------------------------------------------------------------------------
    log.write("Loading variants from VCF...", verbose=verbose)
    mapper = ChromosomeMapper(log=log, verbose=verbose)
    mapper.detect_reference_format(vcf_path)
    
    df = _read_vcf_variants(
        vcf_path=vcf_path,
        region=region,
        vcf_chr_dict=vcf_chr_dict,
        tabix=tabix,
        mapper=mapper,
        log=log,
        verbose=verbose
    )
    
    for col in ("CHR", "POS", "EA", "NEA", "EAF"):
        if col not in df.columns:
            raise ValueError(f"VCF must provide column '{col}'")
    
    # We assume df is already sorted by POS
    pos = df["POS"].to_numpy(dtype=np.int64)
    m = len(df)
    
    if m == 0:
        raise ValueError("No variants returned for the requested region.")
    
    log.write(f" -Loaded {m} variants", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # 0.5) Apply variant filters (MAF/EAF thresholds, exclude rare)
    # -------------------------------------------------------------------------
    # Compute MAF from EAF: MAF = min(EAF, 1-EAF)
    eaf_values = df["EAF"].to_numpy(dtype=float)
    maf_values = np.minimum(eaf_values, 1.0 - eaf_values)
    
    # Build filter mask
    variant_filter_mask = np.ones(m, dtype=bool)
    
    # MAF filters
    if maf_min > 0.0 or maf_max < 0.5:
        maf_mask = (maf_values >= maf_min) & (maf_values <= maf_max)
        variant_filter_mask = variant_filter_mask & maf_mask
        n_filtered_maf = (~maf_mask).sum()
        if n_filtered_maf > 0:
            log.write(f" -Filtered {n_filtered_maf} variants by MAF (keeping MAF in [{maf_min:.4f}, {maf_max:.4f}])", verbose=verbose)
    
    # EAF filters
    if eaf_min > 0.0 or eaf_max < 1.0:
        eaf_mask = (eaf_values >= eaf_min) & (eaf_values <= eaf_max)
        variant_filter_mask = variant_filter_mask & eaf_mask
        n_filtered_eaf = (~eaf_mask).sum()
        if n_filtered_eaf > 0:
            log.write(f" -Filtered {n_filtered_eaf} variants by EAF (keeping EAF in [{eaf_min:.4f}, {eaf_max:.4f}])", verbose=verbose)
    
    # Exclude rare variants
    if exclude_rare:
        rare_mask = maf_values >= rare_maf_threshold
        variant_filter_mask = variant_filter_mask & rare_mask
        n_rare = (~rare_mask).sum()
        if n_rare > 0:
            log.write(f" -Excluded {n_rare} rare variants (MAF < {rare_maf_threshold:.4f})", verbose=verbose)
    
    # Apply filters
    n_before_filter = m
    df = df.loc[variant_filter_mask].reset_index(drop=True)
    m = len(df)
    n_filtered = n_before_filter - m
    
    if n_filtered > 0:
        log.write(f" -After filtering: {m} variants remaining ({n_filtered} excluded)", verbose=verbose)
    
    if m == 0:
        raise ValueError("No variants remaining after applying filters. Please adjust filter criteria.")
    
    # Update position array after filtering
    pos = df["POS"].to_numpy(dtype=np.int64)
    
    # -------------------------------------------------------------------------
    # 1) Baseline N and N_eff
    # -------------------------------------------------------------------------
    if trait == "binary":
        if n_case is None or n_ctrl is None:
            raise ValueError("For trait='binary', you must provide n_case and n_ctrl.")
        # Effective N for case-control:
        #   N_eff ≈ 4 / (1/N_case + 1/N_ctrl)
        neff_base = 4.0 / (1.0 / float(n_case) + 1.0 / float(n_ctrl))
        # Total N is common metadata in summary stats
        N = np.full(m, float(n_case + n_ctrl), dtype=float)
        # Use N_eff for Z scaling
        N_eff = np.full(m, neff_base, dtype=float)
    else:
        if n is None:
            raise ValueError("For trait='quant', you must provide n.")
        N = np.full(m, float(n), dtype=float)
        N_eff = N.copy()
    
    # -------------------------------------------------------------------------
    # 2) Per-SNP N variation (missingness / QC / availability)
    # -------------------------------------------------------------------------
    if n_drop_rate > 0:
        # Vectorized: generate random values and boolean mask
        drop = rng.random(m) < n_drop_rate
        n_dropped = drop.sum()
        if n_dropped > 0:
            # Vectorized: generate multipliers only for dropped SNPs
            mult = rng.uniform(n_drop_min, 1.0, size=n_dropped)
            # Vectorized assignment using boolean indexing
            N[drop] *= mult
            # For quant traits, N_eff follows N; for binary we keep neff_base constant (simple).
            if trait == "quant":
                N_eff[drop] *= mult
    
    # -------------------------------------------------------------------------
    # 3) INFO and INFO attenuation: N_eff <- N_eff * INFO
    # -------------------------------------------------------------------------
    if use_info and "INFO" in df.columns:
        info = np.clip(df["INFO"].to_numpy(dtype=float), info_min, info_max)
    else:
        info = rng.normal(info_mean, info_sd, size=m)
        info = np.clip(info, info_min, info_max)
    
    # Lower INFO => less information => larger SE => smaller Z, via N_eff
    # Ensure N_eff is always positive and valid BEFORE using it in calculations
    N_eff = np.clip(N_eff * info, 1.0, None)  # Minimum of 1.0 to prevent division issues
    
    # -------------------------------------------------------------------------
    # 4) Choose causal effects beta_true
    # -------------------------------------------------------------------------
    beta_true = np.zeros(m, dtype=float)
    is_causal = np.zeros(m, dtype=bool)
    
    # Compute MAF for causal variant filtering
    eaf_values = df["EAF"].to_numpy(dtype=float)
    maf_values = np.minimum(eaf_values, 1.0 - eaf_values)
    
    # Build causal variant candidate mask (variants eligible to be causal)
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
    
    if causal_idx is not None and causal_beta is not None:
        causal_idx = np.asarray(causal_idx, dtype=int)
        causal_beta = np.asarray(causal_beta, dtype=float)
        if causal_idx.ndim != 1 or causal_beta.ndim != 1 or len(causal_idx) != len(causal_beta):
            raise ValueError("causal_idx and causal_beta must be 1D and same length.")
        if np.any((causal_idx < 0) | (causal_idx >= m)):
            raise ValueError("causal_idx contains indices outside the region variant range.")
        
        # Check if explicit causal indices pass the causal filters
        if not np.all(causal_candidate_mask[causal_idx]):
            invalid_causal = causal_idx[~causal_candidate_mask[causal_idx]]
            log.warning(f" -Warning: {len(invalid_causal)} explicit causal variants do not pass causal filters", verbose=verbose)
        
        beta_true[causal_idx] = causal_beta
        is_causal[causal_idx] = True
    else:
        # Select causal variants from candidates only
        if mode == "sparse":
            k = min(int(n_causal), n_candidates)
            if k < int(n_causal):
                log.warning(f" -Requested {n_causal} causal variants but only {n_candidates} candidates available, selecting {k}", verbose=verbose)
            # Select from candidate indices
            selected_candidate_idx = rng.choice(n_candidates, size=k, replace=False)
            idx = causal_candidate_indices[selected_candidate_idx]
        elif mode == "polygenic":
            k = max(1, int(round(n_candidates * float(pi))))
            if k > n_candidates:
                k = n_candidates
                log.warning(f" -Requested {int(round(n_candidates * float(pi)))} causal variants but only {n_candidates} candidates available, selecting {k}", verbose=verbose)
            # Select from candidate indices
            selected_candidate_idx = rng.choice(n_candidates, size=k, replace=False)
            idx = causal_candidate_indices[selected_candidate_idx]
        else:
            raise ValueError("mode must be 'polygenic' or 'sparse'")
        
        beta_true[idx] = rng.normal(0.0, float(effect_sd), size=len(idx))
        is_causal[idx] = True
    
    log.write(f" -Selected {is_causal.sum()} causal variants", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # 5) Simulate Z-scores block-wise: z ~ N(mu, R)
    # -------------------------------------------------------------------------
    blocks = _make_blocks_by_bp(pos, int(window_bp))
    z = np.empty(m, dtype=float)
    
    log.write(f" -Simulating in {len(blocks)} blocks", verbose=verbose)
    
    for block_idx, (start_i, end_i) in enumerate(blocks):
        k = end_i - start_i
        log.write(f" -Processing block {block_idx + 1}/{len(blocks)}: variants {start_i}-{end_i} (size: {k})", verbose=verbose)
        
        # 5.1) Fetch signed LD correlation matrix R for this block
        R = _get_ld_matrix_for_block(
            vcf_path=vcf_path,
            variant_df=df,
            start_idx=start_i,
            end_idx=end_i,
            region=region,
            vcf_chr_dict=vcf_chr_dict,
            tabix=tabix,
            mapper=mapper,
            log=log,
            verbose=verbose
        )
        
        if R.shape != (k, k):
            # Pad or truncate to match block size
            if R.shape[0] < k:
                R_padded = np.eye(k)
                R_padded[:R.shape[0], :R.shape[1]] = R
                R = R_padded
            elif R.shape[0] > k:
                R = R[:k, :k]
        
        # 5.2) LD mismatch realism: shrink R toward identity
        if ld_shrink > 0:
            R = (1.0 - ld_shrink) * R + ld_shrink * np.eye(k)
        
        # 5.3) Make R numerically PSD for Cholesky
        R = _ensure_psd_corr(R, jitter=jitter)
        
        # 5.4) Mean vector:
        # mu = sqrt(N_eff) * (R @ beta_true)
        # This is the expected Z-score under the causal model
        beta_b = beta_true[start_i:end_i]
        
        # Ensure N_eff values are valid for this block
        N_eff_block = N_eff[start_i:end_i]
        N_eff_block = np.clip(N_eff_block, 1.0, None)  # Ensure minimum of 1.0
        
        # Compute mean: mu = sqrt(N_eff) * (R @ beta_b)
        # R @ beta_b spreads causal effects to tagged SNPs via LD
        # Validate R before matrix multiplication
        if np.any(~np.isfinite(R)):
            R = np.where(np.isfinite(R), R, np.eye(k))
        
        R_beta = R @ beta_b
        
        # Validate R_beta
        if np.any(~np.isfinite(R_beta)):
            R_beta = np.where(np.isfinite(R_beta), R_beta, 0.0)
        
        # Validate N_eff_block before sqrt
        if np.any(~np.isfinite(N_eff_block)) or np.any(N_eff_block <= 0):
            N_eff_block = np.where((np.isfinite(N_eff_block)) & (N_eff_block > 0), N_eff_block, 1.0)
        
        sqrt_N_eff = np.sqrt(N_eff_block)
        mu = sqrt_N_eff * R_beta
        
        # Validate mu before sampling
        if np.any(~np.isfinite(mu)):
            mu = np.where(np.isfinite(mu), mu, 0.0)
        
        # 5.5) Sample MVN: z ~ N(mu, R)
        # z = mu + L @ epsilon, where epsilon ~ N(0, I) and R = L L^T
        try:
            L = np.linalg.cholesky(R)
        except np.linalg.LinAlgError:
            L = np.eye(k)
        
        # Generate random noise
        epsilon = rng.standard_normal(k)
        z_block = mu + (L @ epsilon)
        
        # Validate Z scores: replace NaN/Inf with 0
        invalid_z_mask = ~np.isfinite(z_block)
        z_block = np.where(invalid_z_mask, 0.0, z_block)
        
        z[start_i:end_i] = z_block
    
    # -------------------------------------------------------------------------
    # 6) Convert Z -> SE, BETA, P
    # -------------------------------------------------------------------------
    # Validate N_eff to prevent division by zero or negative values
    # (Already validated above, but ensure again for safety)
    N_eff = np.clip(N_eff, 1.0, None)  # Ensure minimum of 1.0
    
    se = 1.0 / np.sqrt(N_eff)
    
    # Validate Z scores before computing betas and P values
    # Check Z score statistics
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
    beta_hat = z_valid * se
    
    # Calculate MLOG10P directly from Z using extreme mode (log-space for numerical precision)
    # This is the same approach as util_in_fill_data.py fill_extreme_mlog10p
    # Uses: log_pvalue = log(2) + norm.logsf(|Z|), then mlog10p = -log_pvalue / log(10)
    mlog10p = _z_to_mlog10p(z_valid)
    # Replace NaN/Inf MLOG10P with 0 (non-significant)
    mlog10p = np.where(np.isfinite(mlog10p), mlog10p, 0.0)
    
    # Calculate P from MLOG10P if needed (for compatibility, but MLOG10P is the primary output)
    # P = 10^(-MLOG10P)
    p = np.power(10.0, -mlog10p)
    # Replace NaN/Inf P values with 0.5 (non-significant)
    p = np.where(np.isfinite(p), p, 0.5)
    
    # Log warning if there were invalid Z scores
    invalid_z_count = np.sum(~np.isfinite(z))
    if invalid_z_count > 0:
        log.warning(f" -Found {invalid_z_count} invalid Z scores (NaN/Inf), replaced with 0", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # 7) Assemble output DataFrame
    # -------------------------------------------------------------------------
    out = df.copy()
    out["INFO"] = info
    out["N"] = N
    out["N_EFF"] = N_eff
    out["IS_CAUSAL"] = is_causal
    out["BETA_TRUE"] = beta_true
    out["Z"] = z
    out["SE"] = se
    out["BETA"] = beta_hat
    out["P"] = p
    out["MLOG10P"] = mlog10p
    
    # Add binary trait specific columns
    if trait == "binary":
        out["N_CASE"] = np.full(m, float(n_case), dtype=float)
        out["N_CONTROL"] = np.full(m, float(n_ctrl), dtype=float)
    
    log.write("Simulation completed successfully!", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # 8) Extract causal SNP IDs
    # -------------------------------------------------------------------------
    causal_indices = np.where(is_causal)[0]
    causal_snp_ids = []
    
    # Check if SNPID column exists (from VCF)
    if "SNPID" in df.columns:
        for idx in causal_indices:
            snp_id = df.iloc[idx]["SNPID"]
            # Use VCF ID if available and not missing
            if snp_id is not None and str(snp_id) != 'None' and str(snp_id) != '.' and str(snp_id) != '':
                causal_snp_ids.append(str(snp_id))
            else:
                # Fallback: construct ID from CHR:POS:EA:NEA
                chr_val = df.iloc[idx]["CHR"]
                pos_val = df.iloc[idx]["POS"]
                ea_val = df.iloc[idx]["EA"]
                nea_val = df.iloc[idx]["NEA"]
                snp_id = f"{chr_val}:{pos_val}:{ea_val}:{nea_val}"
                causal_snp_ids.append(snp_id)
    else:
        # No SNPID column: construct from CHR:POS:EA:NEA
        for idx in causal_indices:
            chr_val = df.iloc[idx]["CHR"]
            pos_val = df.iloc[idx]["POS"]
            ea_val = df.iloc[idx]["EA"]
            nea_val = df.iloc[idx]["NEA"]
            snp_id = f"{chr_val}:{pos_val}:{ea_val}:{nea_val}"
            causal_snp_ids.append(snp_id)
    
    log.write(f" -Extracted {len(causal_snp_ids)} causal SNP IDs", verbose=verbose)
    
    # -------------------------------------------------------------------------
    # 9) Create and return Sumstats object
    # -------------------------------------------------------------------------
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
