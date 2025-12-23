"""
Simplified SSF format validator using only pandas/numpy (no external dependencies).
This module replicates the original gwas_sumstats_tools validation behavior as closely as possible.
"""
from pathlib import Path
from typing import Tuple, List, Optional, Dict
import pandas as pd
import numpy as np
import re
import gzip

from gwaslab.info.g_Log import Log


# SSF format requirements - matching original schema
SSF_REQUIRED_COLUMNS = [
    "chromosome",
    "base_pair_location", 
    "effect_allele",
    "other_allele",
    "standard_error",
    "effect_allele_frequency",
    "p_value"
]

SSF_EFFECT_FIELDS = ["beta", "odds_ratio", "hazard_ratio"]
SSF_PVALUE_FIELDS = ["p_value", "neg_log_10_p_value"]

VALID_FILE_EXTENSIONS = {".tsv", ".tsv.gz"}


def validate_ssf_file(
    filename: Path,
    log: Log = Log(),
    verbose: bool = True,
    minimum_rows: int = 100_000,
    pval_zero: bool = False,
    sample_size: int = 100_000,
    chunksize: int = 1_000_000
) -> Tuple[bool, str, Optional[List[str]]]:
    """
    Validate an SSF format file - replicates original validator behavior.
    
    Parameters
    ----------
    filename : Path
        Path to the SSF file to validate
    log : Log, optional
        Logger instance for messages
    verbose : bool, default True
        Whether to print validation messages
    minimum_rows : int, default 100_000
        Minimum number of rows required
    pval_zero : bool, default False
        Whether to allow p-values of zero
    sample_size : int, default 100_000
        Number of rows to sample for initial validation
    chunksize : int, default 1_000_000
        Chunk size for reading large files
        
    Returns
    -------
    Tuple[bool, str, Optional[List[str]]]
        (is_valid, message, error_list)
    """
    errors = []
    primary_error_type = None
    
    # Step 1: Validate file extension
    if verbose:
        log.write("Validating file extension...", verbose=verbose)
    file_ext = "".join(Path(filename).suffixes)
    if not any(file_ext.endswith(ext) for ext in VALID_FILE_EXTENSIONS):
        primary_error_type = "file_ext"
        errors.append(f"Extension '{file_ext}' not in valid set: {VALID_FILE_EXTENSIONS}")
        return False, f"Extension '{file_ext}' not in valid set: {VALID_FILE_EXTENSIONS}", errors
    
    if verbose:
        log.write("✓ File extension OK", verbose=verbose)
    
    # Step 2: Read header
    try:
        if str(filename).endswith('.gz'):
            with gzip.open(filename, 'rt') as f:
                header_line = f.readline().strip()
        else:
            with open(filename, 'r') as f:
                header_line = f.readline().strip()
        
        header = header_line.split('\t')
    except Exception as e:
        errors.append(f"Failed to read file header: {e}")
        return False, f"Failed to read file: {e}", errors
    
    # Step 3: Validate required columns exist
    missing_required = [col for col in SSF_REQUIRED_COLUMNS if col not in header]
    if missing_required:
        primary_error_type = "headers"
        errors.append(f"Missing required columns: {missing_required}")
        return False, f"Missing required columns: {missing_required}", errors
    
    # Step 4: Detect effect and p-value fields (matching original logic)
    effect_field = _detect_effect_field(header)
    p_value_field = _detect_p_value_field(header)
    
    if effect_field is None:
        primary_error_type = "headers"
        errors.append("Missing effect field: must have at least one of beta, odds_ratio, or hazard_ratio")
        return False, "Missing effect field", errors
    
    if p_value_field is None:
        primary_error_type = "headers"
        errors.append("Missing p-value field: must have p_value or neg_log_10_p_value")
        return False, "Missing p-value field", errors
    
    # Step 5: Validate column order (STRICT - matching original)
    if verbose:
        log.write("Validating column order...", verbose=verbose)
    
    required_order = _get_required_field_order(effect_field, p_value_field)
    actual_order = header[:len(required_order)]
    
    if actual_order != required_order:
        primary_error_type = "field order"
        error_msg = (f"Fields not in the required order:\n"
                    f"Headers given:    {list(actual_order)}\n"
                    f"Headers required: {list(required_order)}")
        errors.append(error_msg)
        return False, error_msg, errors
    
    if verbose:
        log.write("✓ Column order OK", verbose=verbose)
    
    # Step 6: Validate chromosomes (checking for missing autosomes)
    if verbose:
        log.write("Validating chromosomes...", verbose=verbose)
    
    chr_valid, chr_message = _validate_chromosomes(filename, header, log, verbose)
    if not chr_valid:
        primary_error_type = "missing_chromsomes"
        errors.append(chr_message)
        return False, chr_message, errors
    
    if verbose:
        log.write(f"✓ Chromosomes OK ({chr_message})", verbose=verbose)
    
    # Step 7: Read and validate data
    try:
        # Read sample first
        nrows = max(sample_size, minimum_rows)
        if verbose:
            log.write(f"Validating minimum row count...", verbose=verbose)
        
        sample_df = _read_ssf_sample(filename, nrows)
        
        # Check minimum rows
        if len(sample_df) < minimum_rows:
            primary_error_type = "minrows"
            error_msg = f"The file has fewer than the minimum rows required: {len(sample_df)} < {minimum_rows}"
            errors.append(error_msg)
            return False, error_msg, errors
        
        if verbose:
            log.write("✓ Minimum row count OK", verbose=verbose)
        
        # Validate sample data
        if verbose:
            log.write(f"Validating the first {nrows} rows...", verbose=verbose)
        
        sample_errors = _validate_dataframe(
            sample_df, header, effect_field, p_value_field, pval_zero, log, verbose
        )
        
        if sample_errors:
            primary_error_type = sample_errors[0].get('error_type', 'data')
            error_list = [e['message'] if isinstance(e, dict) else str(e) for e in sample_errors]
            errors.extend(error_list)
            return False, f"Data validation failed on sample: {len(sample_errors)} errors found", errors
        
        if verbose:
            log.write("✓ Sample validation OK", verbose=verbose)
        
        # If file is larger, validate remaining chunks
        total_rows = _count_rows(filename)
        if total_rows > nrows:
            if verbose:
                log.write(f"Validating the rest of the file ({total_rows - nrows} rows)...", verbose=verbose)
            
            chunk_errors = _validate_file_chunks(
                filename, header, effect_field, p_value_field, nrows, chunksize, pval_zero, log, verbose
            )
            
            if chunk_errors:
                primary_error_type = chunk_errors[0].get('error_type', 'data')
                error_list = [e['message'] if isinstance(e, dict) else str(e) for e in chunk_errors]
                errors.extend(error_list)
                return False, f"Data validation failed on remaining data: {len(chunk_errors)} errors found", errors
        
        if verbose:
            log.write(f"✓ SSF validation successful: {total_rows} rows validated", verbose=verbose)
        
        return True, f"SSF validation successful: {total_rows} rows validated", None
        
    except Exception as e:
        errors.append(f"Error during data validation: {e}")
        return False, f"Validation error: {e}", errors


def _detect_effect_field(header: List[str]) -> Optional[str]:
    """Detect effect field from header - matches original logic (checks index 4)."""
    if len(header) > 4:
        field_4 = header[4]
        if field_4 in SSF_EFFECT_FIELDS:
            return field_4
    
    # Fallback: check if any effect field exists
    for field in SSF_EFFECT_FIELDS:
        if field in header:
            return field
    
    return None


def _detect_p_value_field(header: List[str]) -> Optional[str]:
    """Detect p-value field from header - matches original logic (checks index 7)."""
    if len(header) > 7:
        field_7 = header[7]
        if field_7 in SSF_PVALUE_FIELDS:
            return field_7
    
    # Fallback: check if any p-value field exists
    for field in SSF_PVALUE_FIELDS:
        if field in header:
            return field
    
    return None


def _get_required_field_order(effect_field: str, p_value_field: str) -> List[str]:
    """Get required field order matching original schema.field_order()."""
    # Base required fields (excluding effect and p-value which vary)
    base_fields = [
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        effect_field,  # Insert effect field at position 4
        "standard_error",
        "effect_allele_frequency",
        p_value_field  # Insert p-value field at position 7
    ]
    return base_fields


def _validate_chromosomes(
    filename: Path,
    header: List[str],
    log: Log,
    verbose: bool
) -> Tuple[bool, str]:
    """
    Validate chromosomes - matches original logic:
    - Requires all autosomes (1-22)
    - Treats 23-25 (X, Y, MT) as optional
    """
    if "chromosome" not in header:
        return False, "Chromosome column is missing from the input file."
    
    # Read chromosome column
    try:
        if str(filename).endswith('.gz'):
            chr_df = pd.read_csv(filename, sep='\t', usecols=['chromosome'], 
                               compression='gzip', low_memory=False, dtype=str)
        else:
            chr_df = pd.read_csv(filename, sep='\t', usecols=['chromosome'], 
                               low_memory=False, dtype=str)
        
        unique_chr = set(chr_df['chromosome'].astype(str).unique())
        autosomes_chromosomes = set(map(str, range(1, 23)))
        optional_chromosomes = set(map(str, range(23, 26)))
        
        missing_autosomes = sorted(autosomes_chromosomes - unique_chr, key=int)
        missing_optional = sorted(optional_chromosomes - unique_chr, key=int)
        
        # Special case: only chromosome X
        if unique_chr == {"23"}:
            return True, "This file only contains chromosome X."
        
        # Fail if any autosomes are missing
        if missing_autosomes:
            return False, f"Chromosome column missing values: {missing_autosomes}"
        
        # Warn about optional chromosomes but don't fail
        if missing_optional:
            return True, f"All autosomes exist. Optional chromosomes {missing_optional} do not exist."
        
        return True, "All chromosomes, including X, Y, and MT, exist."
        
    except Exception as e:
        return False, f"Error validating chromosomes: {e}"


def _read_ssf_sample(filename: Path, nrows: int) -> pd.DataFrame:
    """Read a sample of rows from SSF file with type coercion."""
    na_values = ["", "#NA", "NA", "N/A", "NaN", "NR"]
    
    if str(filename).endswith('.gz'):
        df = pd.read_csv(filename, sep='\t', nrows=nrows, compression='gzip', 
                        low_memory=False, na_values=na_values, dtype=str)
    else:
        df = pd.read_csv(filename, sep='\t', nrows=nrows, low_memory=False, 
                        na_values=na_values, dtype=str)
    
    # Coerce numeric columns (matching Pandera's coerce=True behavior)
    return _coerce_numeric_types(df)


def _coerce_numeric_types(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce string numbers to appropriate numeric types."""
    df = df.copy()
    
    # Numeric columns
    numeric_cols = ['chromosome', 'base_pair_location', 'beta', 'odds_ratio', 
                    'hazard_ratio', 'standard_error', 'effect_allele_frequency',
                    'p_value', 'neg_log_10_p_value', 'ci_upper', 'ci_lower', 
                    'info', 'n']
    
    for col in numeric_cols:
        if col in df.columns:
            try:
                if col in ['chromosome', 'base_pair_location', 'n']:
                    # Integer columns
                    df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
                else:
                    # Float columns
                    df[col] = pd.to_numeric(df[col], errors='coerce').astype(float)
            except:
                pass  # Keep as string if conversion fails
    
    return df


def _count_rows(filename: Path) -> int:
    """Count total rows in file (excluding header)."""
    if str(filename).endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return sum(1 for _ in f) - 1  # Subtract header
    else:
        with open(filename, 'r') as f:
            return sum(1 for _ in f) - 1  # Subtract header


def _validate_file_chunks(
    filename: Path,
    header: List[str],
    effect_field: str,
    p_value_field: str,
    skiprows: int,
    chunksize: int,
    pval_zero: bool,
    log: Log,
    verbose: bool
) -> List[Dict[str, str]]:
    """Validate file in chunks."""
    errors = []
    na_values = ["", "#NA", "NA", "N/A", "NaN", "NR"]
    
    try:
        if str(filename).endswith('.gz'):
            chunk_iter = pd.read_csv(filename, sep='\t', chunksize=chunksize, 
                                    skiprows=range(1, skiprows + 1),
                                    compression='gzip', low_memory=False,
                                    na_values=na_values, dtype=str)
        else:
            chunk_iter = pd.read_csv(filename, sep='\t', chunksize=chunksize, 
                                    skiprows=range(1, skiprows + 1),
                                    low_memory=False, na_values=na_values, dtype=str)
        
        for i, chunk_df in enumerate(chunk_iter):
            chunk_df = _coerce_numeric_types(chunk_df)
            chunk_errors = _validate_dataframe(
                chunk_df, header, effect_field, p_value_field, pval_zero, log, False
            )
            errors.extend(chunk_errors)
            if len(errors) > 10:  # Limit error reporting
                break
                
    except Exception as e:
        errors.append({'error_type': 'data', 'message': f"Error reading chunks: {e}"})
    
    return errors


def _validate_dataframe(
    df: pd.DataFrame,
    header: List[str],
    effect_field: str,
    p_value_field: str,
    pval_zero: bool,
    log: Log,
    verbose: bool
) -> List[Dict[str, str]]:
    """
    Validate a pandas DataFrame against SSF requirements.
    Matches original Pandera validation logic.
    """
    errors = []
    
    # Split p-values into mantissa and exponent for precision (CRITICAL)
    df = _pval_to_mantissa_and_exponent(df, p_value_field)
    
    # Validate chromosome (1-25)
    if "chromosome" in df.columns:
        invalid_chr = df[~df["chromosome"].between(1, 25, inclusive='both')]
        if len(invalid_chr) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid chromosome values: {len(invalid_chr)} rows with chromosome not in [1-25]"
            })
    
    # Validate base_pair_location (>= 0)
    if "base_pair_location" in df.columns:
        invalid_pos = df[df["base_pair_location"] < 0]
        if len(invalid_pos) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid positions: {len(invalid_pos)} rows with negative base_pair_location"
            })
    
    # Validate effect_allele and other_allele (nucleotides)
    if "effect_allele" in df.columns:
        invalid_ea = df[~df["effect_allele"].astype(str).str.match(r'^(LONG_STRING|[ACTGactg]+)$', na=False)]
        if len(invalid_ea) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid effect_allele: {len(invalid_ea)} rows with non-nucleotide sequences"
            })
    
    if "other_allele" in df.columns:
        invalid_oa = df[~df["other_allele"].astype(str).str.match(r'^(LONG_STRING|[ACTGactg]+)$', na=False)]
        if len(invalid_oa) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid other_allele: {len(invalid_oa)} rows with non-nucleotide sequences"
            })
    
    # Validate effect field (only the detected one)
    if effect_field in df.columns:
        if effect_field in ["odds_ratio", "hazard_ratio"]:
            invalid = df[df[effect_field] < 0]
            if len(invalid) > 0:
                errors.append({
                    'error_type': 'data',
                    'message': f"Invalid {effect_field}: {len(invalid)} rows with negative values"
                })
    
    # Validate effect_allele_frequency (0-1)
    if "effect_allele_frequency" in df.columns:
        invalid_eaf = df[~df["effect_allele_frequency"].between(0, 1, inclusive='both')]
        if len(invalid_eaf) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid effect_allele_frequency: {len(invalid_eaf)} rows outside [0, 1]"
            })
    
    # Validate p-value (0-1)
    if p_value_field in df.columns:
        invalid_p = df[~df[p_value_field].between(0, 1, inclusive='both')]
        if len(invalid_p) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid {p_value_field}: {len(invalid_p)} rows outside [0, 1]"
            })
    
    # Validate p-value mantissa (CRITICAL for very small p-values)
    if "_p_value_mantissa" in df.columns:
        if pval_zero:
            invalid_mantissa = df[df["_p_value_mantissa"] < 0]
        else:
            invalid_mantissa = df[df["_p_value_mantissa"] <= 0]
        
        if len(invalid_mantissa) > 0:
            error_type = 'p_val' if not pval_zero else 'data'
            errors.append({
                'error_type': error_type,
                'message': f"Invalid p-value mantissa: {len(invalid_mantissa)} rows with mantissa {'<=' if not pval_zero else '<'} 0"
            })
    
    # Validate neg_log_10_p_value (>= 0) if present
    if "neg_log_10_p_value" in df.columns:
        invalid_mlog = df[df["neg_log_10_p_value"] < 0]
        if len(invalid_mlog) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid neg_log_10_p_value: {len(invalid_mlog)} rows with negative values"
            })
    
    # Validate info (0-1) if present
    if "info" in df.columns:
        invalid_info = df[~df["info"].between(0, 1, inclusive='both')]
        if len(invalid_info) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid info: {len(invalid_info)} rows outside [0, 1]"
            })
    
    # Validate n (>= 0) if present
    if "n" in df.columns:
        invalid_n = df[df["n"] < 0]
        if len(invalid_n) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid n: {len(invalid_n)} rows with negative values"
            })
    
    # Validate rsid pattern if present (strict - no empty strings)
    if "rsid" in df.columns:
        invalid_rsid = df[~df["rsid"].astype(str).str.match(r'^rs[0-9]+$', na=False)]
        if len(invalid_rsid) > 0:
            errors.append({
                'error_type': 'data',
                'message': f"Invalid rsid format: {len(invalid_rsid)} rows with non-standard rsID format"
            })
    
    return errors


def _pval_to_mantissa_and_exponent(df: pd.DataFrame, p_value_field: str) -> pd.DataFrame:
    """
    Split p-values into mantissa and exponent for precision validation.
    This is CRITICAL for handling very small p-values that would round to 0 in float64.
    Matches original pval_to_mantissa_and_exponent() logic.
    """
    df = df.copy()
    
    if p_value_field not in df.columns:
        return df
    
    # Convert p-value to string for splitting
    pval_str = df[p_value_field].astype(str)
    
    # Split on 'e' or 'E' (scientific notation)
    split_result = pval_str.str.split(r'e|E', n=1, expand=True, regex=True)
    
    if len(split_result.columns) >= 2:
        df["_p_value_mantissa"] = pd.to_numeric(split_result[0], errors='coerce')
        df["_p_value_exponent"] = pd.to_numeric(split_result[1], errors='coerce').astype('Int64')
    else:
        # No exponent found, mantissa is the whole value
        df["_p_value_mantissa"] = pd.to_numeric(pval_str, errors='coerce')
        df["_p_value_exponent"] = pd.NA
    
    return df
