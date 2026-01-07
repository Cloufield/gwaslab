"""
UCSC BEDPE format I/O utilities.

This module provides functions for reading and processing BEDPE (Browser Extensible Data Paired-End)
format files according to the UCSC Genome Browser specification.

BEDPE format specification: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

BEDPE format uses:
- 0-based, half-open intervals [chromStart, chromEnd) for both ends
- chromStart: 0-based (first base in chromosome is numbered 0)
- chromEnd: exclusive (not included in the feature display)
- Required fields: chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2
- Optional fields: name, score, strand1, strand2
"""

from typing import TYPE_CHECKING, Union, Optional, Tuple, List
import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    pass

# BEDPE file suffix definitions
BEDPE_SUFFIXES = ('.bedpe', '.bedpe.gz')

# BEDPE field names (standard BEDPE10 format)
BEDPE_REQUIRED_FIELDS = ['chrom1', 'chromStart1', 'chromEnd1', 'chrom2', 'chromStart2', 'chromEnd2']
BEDPE_OPTIONAL_FIELDS = ['name', 'score', 'strand1', 'strand2']
BEDPE_ALL_FIELDS = BEDPE_REQUIRED_FIELDS + BEDPE_OPTIONAL_FIELDS


def read_bedpe(
    bedpe_path: str,
    usecols: Optional[List[int]] = None,
    header: Optional[bool] = None,
    comment: Optional[str] = None,
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Read a BEDPE file into a pandas DataFrame.
    
    BEDPE files use 0-based, half-open intervals [chromStart, chromEnd) for both ends.
    The chromEnd position is exclusive (not included in the feature).
    
    Parameters
    ----------
    bedpe_path : str
        Path to BEDPE file. Supports uncompressed (.bedpe) and gzipped (.bedpe.gz) files.
    usecols : list of int, optional
        Column indices to read (0-based). If None, reads all columns.
        Default columns: chrom1 (0), chromStart1 (1), chromEnd1 (2), chrom2 (3), 
        chromStart2 (4), chromEnd2 (5), and optional fields.
    header : bool, optional
        Whether the file has a header line. BEDPE files typically don't have headers,
        but custom tracks may have "browser" or "track" lines.
        If None, auto-detects by checking for "browser" or "track" lines.
    comment : str, optional
        Character(s) that indicate comment lines to skip. Default is "#".
        BEDPE custom tracks may have "browser" or "track" header lines.
    verbose : bool, default=False
        If True, prints information about the file being read.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with BEDPE data. Column names depend on number of columns:
        - BEDPE6: chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2
        - BEDPE7+: chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2, name, ...
        Columns are named according to BEDPE specification.
    
    Notes
    -----
    BEDPE format specification:
    - chrom1/chrom2: Chromosome names (e.g., chr1, chrX, 1, X)
    - chromStart1/chromStart2: Starting positions (0-based)
    - chromEnd1/chromEnd2: Ending positions (1-based, exclusive)
    - Optional fields: name, score, strand1, strand2
    
    Coordinate system:
    - BEDPE uses 0-based, half-open intervals [chromStart, chromEnd)
    - For example, chromStart=0, chromEnd=100 spans bases 0-99
    - In 1-based coordinates, this corresponds to positions 1-100
    
    References
    ----------
    UCSC BEDPE format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """
    if not os.path.exists(bedpe_path):
        raise FileNotFoundError(f"BEDPE file not found: {bedpe_path}")
    
    if log and verbose:
        log.write(f" -Reading BEDPE file: {bedpe_path}", verbose=verbose)
    
    # Determine if file is compressed
    is_gzipped = bedpe_path.endswith('.gz')
    
    # Auto-detect header lines (browser or track lines)
    if header is None:
        header = _detect_bedpe_header(bedpe_path, is_gzipped)
    
    # Read BEDPE file
    if is_gzipped:
        bedpe_df = pd.read_csv(
            bedpe_path,
            sep=r"\s+",
            header=None,
            comment=comment if comment is not None else "#",
            compression='gzip' if is_gzipped else None,
            dtype={0: "string", 1: "Int64", 2: "Int64", 3: "string", 4: "Int64", 5: "Int64"}
        )
    else:
        bedpe_df = pd.read_csv(
            bedpe_path,
            sep=r"\s+",
            header=None,
            comment=comment if comment is not None else "#",
            dtype={0: "string", 1: "Int64", 2: "Int64", 3: "string", 4: "Int64", 5: "Int64"}
        )
    
    if len(bedpe_df) == 0:
        if log:
            log.write(" -Warning: BEDPE file is empty.", verbose=verbose)
        return pd.DataFrame(columns=BEDPE_REQUIRED_FIELDS)
    
    # Select columns if specified
    if usecols is not None:
        bedpe_df = bedpe_df.iloc[:, usecols]
    
    # Name columns based on BEDPE format
    n_cols = bedpe_df.shape[1]
    if n_cols >= 6:
        # At least BEDPE6 format
        bedpe_df.columns = BEDPE_ALL_FIELDS[:n_cols]
    else:
        # Invalid BEDPE format (less than 6 columns)
        raise ValueError(f"BEDPE file must have at least 6 columns, found {n_cols}")
    
    # Validate BEDPE coordinates
    if not _validate_bedpe_coordinates(bedpe_df):
        if log:
            log.warning(" -Warning: BEDPE file contains invalid coordinates (chromStart >= chromEnd).")
    
    return bedpe_df


def _detect_bedpe_header(bedpe_path: str, is_gzipped: bool) -> bool:
    """
    Detect if BEDPE file has header lines (browser or track lines).
    
    Parameters
    ----------
    bedpe_path : str
        Path to BEDPE file
    is_gzipped : bool
        Whether the file is gzipped
    
    Returns
    -------
    bool
        True if header lines are detected, False otherwise
    """
    try:
        if is_gzipped:
            f = gzip.open(bedpe_path, 'rt')
        else:
            f = open(bedpe_path, 'r')
        
        # Check first few lines for "browser" or "track"
        for i, line in enumerate(f):
            if i >= 10:  # Check first 10 lines
                break
            line_lower = line.strip().lower()
            if line_lower.startswith('browser') or line_lower.startswith('track'):
                f.close()
                return True
        
        f.close()
        return False
    except Exception:
        return False


def _validate_bedpe_coordinates(bedpe_df: pd.DataFrame) -> bool:
    """
    Validate BEDPE coordinates (chromStart < chromEnd for both ends).
    
    Parameters
    ----------
    bedpe_df : pd.DataFrame
        BEDPE DataFrame with chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2 columns
    
    Returns
    -------
    bool
        True if all coordinates are valid, False otherwise
    """
    required_cols = ['chromStart1', 'chromEnd1', 'chromStart2', 'chromEnd2']
    if not all(col in bedpe_df.columns for col in required_cols):
        return False
    
    # Check that chromStart < chromEnd for both ends
    valid1 = (bedpe_df['chromStart1'] < bedpe_df['chromEnd1']).all()
    valid2 = (bedpe_df['chromStart2'] < bedpe_df['chromEnd2']).all()
    return valid1 and valid2


def bedpe_to_sumstats_coordinates(
    bedpe_df: pd.DataFrame,
    mapper: Optional[ChromosomeMapper] = None,
    chrom1_col: str = 'chrom1',
    start1_col: str = 'chromStart1',
    end1_col: str = 'chromEnd1',
    chrom2_col: str = 'chrom2',
    start2_col: str = 'chromStart2',
    end2_col: str = 'chromEnd2',
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Convert BEDPE coordinates to sumstats-compatible format.
    
    Converts:
    - Chromosome notation: BEDPE chrom1/chrom2 → sumstats CHR format
    - Coordinates: BEDPE 0-based [start, end) → 1-based [start+1, end] for matching
    
    Parameters
    ----------
    bedpe_df : pd.DataFrame
        BEDPE DataFrame with chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2 columns
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome conversion.
        If None, creates a new mapper.
    chrom1_col : str, default='chrom1'
        Column name for chromosome 1 in BEDPE file
    start1_col : str, default='chromStart1'
        Column name for start position 1 in BEDPE file
    end1_col : str, default='chromEnd1'
        Column name for end position 1 in BEDPE file
    chrom2_col : str, default='chrom2'
        Column name for chromosome 2 in BEDPE file
    start2_col : str, default='chromStart2'
        Column name for start position 2 in BEDPE file
    end2_col : str, default='chromEnd2'
        Column name for end position 2 in BEDPE file
    verbose : bool, default=False
        If True, prints conversion information
    log : Log, optional
        Logger instance for messages
    
    Returns
    -------
    pd.DataFrame
        DataFrame with converted coordinates:
        - CHR1: Chromosome 1 in sumstats format
        - START1: Start position 1 for matching (1-based, inclusive)
        - END1: End position 1 for matching (1-based, inclusive)
        - CHR2: Chromosome 2 in sumstats format
        - START2: Start position 2 for matching (1-based, inclusive)
        - END2: End position 2 for matching (1-based, inclusive)
        Original BEDPE columns are preserved with '_bedpe' suffix.
    
    Notes
    -----
    BEDPE uses 0-based, half-open intervals [chromStart, chromEnd).
    For matching with 1-based sumstats POS:
    - BEDPE [start, end) in 0-based = [start+1, end] in 1-based
    - We convert to: START = start + 1, END = end (both inclusive for matching)
    """
    if mapper is None:
        mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    
    # Create a copy to avoid modifying original
    result_df = bedpe_df.copy()
    
    # Convert chromosome notation to numeric for chrom1
    bedpe_chr1_numeric = result_df[chrom1_col].apply(
        lambda x: mapper.sumstats_to_number(str(x).strip().lstrip("chr").lstrip("CHR"))
        if pd.notna(x) else None
    )
    
    # Convert chromosome notation to numeric for chrom2
    bedpe_chr2_numeric = result_df[chrom2_col].apply(
        lambda x: mapper.sumstats_to_number(str(x).strip().lstrip("chr").lstrip("CHR"))
        if pd.notna(x) else None
    )
    
    # Filter out unconvertible chromosomes
    valid_mask = bedpe_chr1_numeric.notna() & bedpe_chr2_numeric.notna()
    n_invalid = (~valid_mask).sum()
    
    if n_invalid > 0 and log:
        log.write(f" -Warning: {n_invalid} BEDPE regions have unconvertible chromosome notation.", verbose=verbose)
    
    result_df = result_df[valid_mask].copy()
    bedpe_chr1_numeric = bedpe_chr1_numeric[valid_mask]
    bedpe_chr2_numeric = bedpe_chr2_numeric[valid_mask]
    
    # Convert numeric chromosomes to sumstats format
    result_df['CHR1'] = bedpe_chr1_numeric.apply(
        lambda x: mapper.number_to_sumstats(x) if pd.notna(x) else None
    )
    result_df['CHR2'] = bedpe_chr2_numeric.apply(
        lambda x: mapper.number_to_sumstats(x) if pd.notna(x) else None
    )
    
    # Convert BEDPE coordinates to 1-based for matching
    # BEDPE [start, end) in 0-based = [start+1, end] in 1-based
    # For matching, we use: START = start + 1, END = end (both inclusive)
    result_df['START1'] = result_df[start1_col] + 1  # Convert 0-based to 1-based
    result_df['END1'] = result_df[end1_col]  # End is already 1-based (exclusive in BEDPE = inclusive for matching)
    result_df['START2'] = result_df[start2_col] + 1  # Convert 0-based to 1-based
    result_df['END2'] = result_df[end2_col]  # End is already 1-based (exclusive in BEDPE = inclusive for matching)
    
    # Rename original columns with _bedpe suffix
    rename_dict = {
        chrom1_col: f'{chrom1_col}_bedpe',
        start1_col: f'{start1_col}_bedpe',
        end1_col: f'{end1_col}_bedpe',
        chrom2_col: f'{chrom2_col}_bedpe',
        start2_col: f'{start2_col}_bedpe',
        end2_col: f'{end2_col}_bedpe'
    }
    result_df = result_df.rename(columns=rename_dict)
    
    return result_df


def bedpe_contains_position(
    bedpe_df: pd.DataFrame,
    chrom: Union[str, int],
    pos: int,
    chrom1_col: str = 'CHR1',
    start1_col: str = 'START1',
    end1_col: str = 'END1',
    chrom2_col: str = 'CHR2',
    start2_col: str = 'START2',
    end2_col: str = 'END2',
    mapper: Optional[ChromosomeMapper] = None
) -> bool:
    """
    Check if a position (chrom, pos) is contained in any BEDPE region (either end).
    
    Parameters
    ----------
    bedpe_df : pd.DataFrame
        BEDPE DataFrame with converted coordinates (from bedpe_to_sumstats_coordinates)
    chrom : str or int
        Chromosome identifier (in sumstats format)
    pos : int
        Position (1-based)
    chrom1_col : str, default='CHR1'
        Column name for chromosome 1 in converted BEDPE DataFrame
    start1_col : str, default='START1'
        Column name for start position 1 in converted BEDPE DataFrame
    end1_col : str, default='END1'
        Column name for end position 1 in converted BEDPE DataFrame
    chrom2_col : str, default='CHR2'
        Column name for chromosome 2 in converted BEDPE DataFrame
    start2_col : str, default='START2'
        Column name for start position 2 in converted BEDPE DataFrame
    end2_col : str, default='END2'
        Column name for end position 2 in converted BEDPE DataFrame
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome conversion.
        If None, assumes chrom is already in the same format as bedpe_df[chrom1_col].
    
    Returns
    -------
    bool
        True if position is contained in any BEDPE region (either end), False otherwise
    """
    if len(bedpe_df) == 0:
        return False
    
    # Convert chrom to match BEDPE format if mapper is provided
    if mapper is not None:
        chrom_num = mapper.sumstats_to_number(chrom)
        if pd.isna(chrom_num):
            return False
        bedpe_chr1_num = bedpe_df[chrom1_col].apply(
            lambda x: mapper.sumstats_to_number(x) if pd.notna(x) else None
        )
        bedpe_chr2_num = bedpe_df[chrom2_col].apply(
            lambda x: mapper.sumstats_to_number(x) if pd.notna(x) else None
        )
        matches1 = (bedpe_chr1_num == chrom_num) & (bedpe_df[start1_col] <= pos) & (pos <= bedpe_df[end1_col])
        matches2 = (bedpe_chr2_num == chrom_num) & (bedpe_df[start2_col] <= pos) & (pos <= bedpe_df[end2_col])
    else:
        # Direct comparison (assumes same format)
        matches1 = (bedpe_df[chrom1_col] == chrom) & (bedpe_df[start1_col] <= pos) & (pos <= bedpe_df[end1_col])
        matches2 = (bedpe_df[chrom2_col] == chrom) & (bedpe_df[start2_col] <= pos) & (pos <= bedpe_df[end2_col])
    
    return (matches1 | matches2).any()

