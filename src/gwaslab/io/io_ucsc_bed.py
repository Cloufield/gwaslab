"""
UCSC BED format I/O utilities.

This module provides functions for reading and processing BED (Browser Extensible Data)
format files according to the UCSC Genome Browser specification.

BED format specification: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

BED format uses:
- 0-based, half-open intervals [chromStart, chromEnd)
- chromStart: 0-based (first base in chromosome is numbered 0)
- chromEnd: exclusive (not included in the feature display)
- Required fields: chrom, chromStart, chromEnd
- Optional fields: name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
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

# BED file suffix definitions
BED_SUFFIXES = ('.bed', '.bed.gz')

# BED field names (standard BED12 format)
BED_REQUIRED_FIELDS = ['chrom', 'chromStart', 'chromEnd']
BED_OPTIONAL_FIELDS = ['name', 'score', 'strand', 'thickStart', 'thickEnd', 
                       'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
BED_ALL_FIELDS = BED_REQUIRED_FIELDS + BED_OPTIONAL_FIELDS


def read_bed(
    bed_path: str,
    usecols: Optional[List[int]] = None,
    header: Optional[bool] = None,
    comment: Optional[str] = None,
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Read a BED file into a pandas DataFrame.
    
    BED files use 0-based, half-open intervals [chromStart, chromEnd).
    The chromEnd position is exclusive (not included in the feature).
    
    Parameters
    ----------
    bed_path : str
        Path to BED file. Supports uncompressed (.bed) and gzipped (.bed.gz) files.
    usecols : list of int, optional
        Column indices to read (0-based). If None, reads all columns.
        Default columns: chrom (0), chromStart (1), chromEnd (2), and optional fields.
    header : bool, optional
        Whether the file has a header line. BED files typically don't have headers,
        but custom tracks may have "browser" or "track" lines.
        If None, auto-detects by checking for "browser" or "track" lines.
    comment : str, optional
        Character(s) that indicate comment lines to skip. Default is "#".
        BED custom tracks may have "browser" or "track" header lines.
    verbose : bool, default=False
        If True, prints information about the file being read.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with BED data. Column names depend on number of columns:
        - BED3: chrom, chromStart, chromEnd
        - BED4+: chrom, chromStart, chromEnd, name, ...
        Columns are named according to BED specification.
    
    Notes
    -----
    BED format specification:
    - chrom: Chromosome name (e.g., chr1, chrX, 1, X)
    - chromStart: Starting position (0-based)
    - chromEnd: Ending position (1-based, exclusive)
    - Optional fields follow the standard BED12 format
    
    Coordinate system:
    - BED uses 0-based, half-open intervals [chromStart, chromEnd)
    - For example, chromStart=0, chromEnd=100 spans bases 0-99
    - In 1-based coordinates, this corresponds to positions 1-100
    
    References
    ----------
    UCSC BED format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    
    if log and verbose:
        log.write(f" -Reading BED file: {bed_path}", verbose=verbose)
    
    # Determine if file is compressed
    is_gzipped = bed_path.endswith('.gz')
    
    # Auto-detect header lines (browser or track lines)
    if header is None:
        header = _detect_bed_header(bed_path, is_gzipped)
    
    # Read BED file
    if is_gzipped:
        bed_df = pd.read_csv(
            bed_path,
            sep=r"\s+",
            header=None,
            comment=comment if comment is not None else "#",
            compression='gzip' if is_gzipped else None,
            dtype={0: "string", 1: "Int64", 2: "Int64"}
        )
    else:
        bed_df = pd.read_csv(
            bed_path,
            sep=r"\s+",
            header=None,
            comment=comment if comment is not None else "#",
            dtype={0: "string", 1: "Int64", 2: "Int64"}
        )
    
    if len(bed_df) == 0:
        if log:
            log.write(" -Warning: BED file is empty.", verbose=verbose)
        return pd.DataFrame(columns=BED_REQUIRED_FIELDS)
    
    # Select columns if specified
    if usecols is not None:
        bed_df = bed_df.iloc[:, usecols]
    
    # Name columns based on BED format
    n_cols = bed_df.shape[1]
    if n_cols >= 3:
        # At least BED3 format
        bed_df.columns = BED_ALL_FIELDS[:n_cols]
    else:
        # Invalid BED format (less than 3 columns)
        raise ValueError(f"BED file must have at least 3 columns, found {n_cols}")
    
    # Validate BED coordinates
    if not _validate_bed_coordinates(bed_df):
        if log:
            log.warning(" -Warning: BED file contains invalid coordinates (chromStart >= chromEnd).")
    
    return bed_df


def _detect_bed_header(bed_path: str, is_gzipped: bool) -> bool:
    """
    Detect if BED file has header lines (browser or track lines).
    
    Parameters
    ----------
    bed_path : str
        Path to BED file
    is_gzipped : bool
        Whether the file is gzipped
    
    Returns
    -------
    bool
        True if header lines are detected, False otherwise
    """
    try:
        if is_gzipped:
            f = gzip.open(bed_path, 'rt')
        else:
            f = open(bed_path, 'r')
        
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


def _validate_bed_coordinates(bed_df: pd.DataFrame) -> bool:
    """
    Validate BED coordinates (chromStart < chromEnd).
    
    Parameters
    ----------
    bed_df : pd.DataFrame
        BED DataFrame with chrom, chromStart, chromEnd columns
    
    Returns
    -------
    bool
        True if all coordinates are valid, False otherwise
    """
    if 'chromStart' not in bed_df.columns or 'chromEnd' not in bed_df.columns:
        return False
    
    # Check that chromStart < chromEnd for all rows
    valid = (bed_df['chromStart'] < bed_df['chromEnd']).all()
    return valid


def bed_to_sumstats_coordinates(
    bed_df: pd.DataFrame,
    mapper: Optional[ChromosomeMapper] = None,
    chrom_col: str = 'chrom',
    start_col: str = 'chromStart',
    end_col: str = 'chromEnd',
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Convert BED coordinates to sumstats-compatible format.
    
    Converts:
    - Chromosome notation: BED chrom → sumstats CHR format
    - Coordinates: BED 0-based [start, end) → 1-based [start+1, end] for matching
    
    Parameters
    ----------
    bed_df : pd.DataFrame
        BED DataFrame with chrom, chromStart, chromEnd columns
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome conversion.
        If None, creates a new mapper.
    chrom_col : str, default='chrom'
        Column name for chromosome in BED file
    start_col : str, default='chromStart'
        Column name for start position in BED file
    end_col : str, default='chromEnd'
        Column name for end position in BED file
    verbose : bool, default=False
        If True, prints conversion information
    log : Log, optional
        Logger instance for messages
    
    Returns
    -------
    pd.DataFrame
        DataFrame with converted coordinates:
        - CHR: Chromosome in sumstats format
        - START: Start position for matching (1-based, inclusive)
        - END: End position for matching (1-based, inclusive)
        Original BED columns are preserved with '_bed' suffix.
    
    Notes
    -----
    BED uses 0-based, half-open intervals [chromStart, chromEnd).
    For matching with 1-based sumstats POS:
    - BED [start, end) in 0-based = [start+1, end] in 1-based
    - We convert to: START = start + 1, END = end (both inclusive for matching)
    """
    if mapper is None:
        mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    
    # Create a copy to avoid modifying original
    result_df = bed_df.copy()
    
    # Convert chromosome notation to numeric (middle layer)
    bed_chr_numeric = result_df[chrom_col].apply(
        lambda x: mapper.sumstats_to_number(str(x).strip().lstrip("chr").lstrip("CHR"))
        if pd.notna(x) else None
    )
    
    # Filter out unconvertible chromosomes
    valid_mask = bed_chr_numeric.notna()
    n_invalid = (~valid_mask).sum()
    
    if n_invalid > 0 and log:
        log.write(f" -Warning: {n_invalid} BED regions have unconvertible chromosome notation.", verbose=verbose)
    
    result_df = result_df[valid_mask].copy()
    bed_chr_numeric = bed_chr_numeric[valid_mask]
    
    # Convert numeric chromosomes to sumstats format
    result_df['CHR'] = bed_chr_numeric.apply(
        lambda x: mapper.number_to_sumstats(x) if pd.notna(x) else None
    )
    
    # Convert BED coordinates to 1-based for matching
    # BED [start, end) in 0-based = [start+1, end] in 1-based
    # For matching, we use: START = start + 1, END = end (both inclusive)
    result_df['START'] = result_df[start_col] + 1  # Convert 0-based to 1-based
    result_df['END'] = result_df[end_col]  # End is already 1-based (exclusive in BED = inclusive for matching)
    
    # Rename original columns with _bed suffix
    rename_dict = {chrom_col: f'{chrom_col}_bed', start_col: f'{start_col}_bed', end_col: f'{end_col}_bed'}
    result_df = result_df.rename(columns=rename_dict)
    
    return result_df


def bed_contains_position(
    bed_df: pd.DataFrame,
    chrom: Union[str, int],
    pos: int,
    chrom_col: str = 'CHR',
    start_col: str = 'START',
    end_col: str = 'END',
    mapper: Optional[ChromosomeMapper] = None
) -> bool:
    """
    Check if a position (chrom, pos) is contained in any BED region.
    
    Parameters
    ----------
    bed_df : pd.DataFrame
        BED DataFrame with converted coordinates (from bed_to_sumstats_coordinates)
    chrom : str or int
        Chromosome identifier (in sumstats format)
    pos : int
        Position (1-based)
    chrom_col : str, default='CHR'
        Column name for chromosome in converted BED DataFrame
    start_col : str, default='START'
        Column name for start position in converted BED DataFrame
    end_col : str, default='END'
        Column name for end position in converted BED DataFrame
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome conversion.
        If None, assumes chrom is already in the same format as bed_df[chrom_col].
    
    Returns
    -------
    bool
        True if position is contained in any BED region, False otherwise
    """
    if len(bed_df) == 0:
        return False
    
    # Convert chrom to match BED format if mapper is provided
    if mapper is not None:
        chrom_num = mapper.sumstats_to_number(chrom)
        if pd.isna(chrom_num):
            return False
        bed_chr_num = bed_df[chrom_col].apply(
            lambda x: mapper.sumstats_to_number(x) if pd.notna(x) else None
        )
        matches = (bed_chr_num == chrom_num) & (bed_df[start_col] <= pos) & (pos <= bed_df[end_col])
    else:
        # Direct comparison (assumes same format)
        matches = (bed_df[chrom_col] == chrom) & (bed_df[start_col] <= pos) & (pos <= bed_df[end_col])
    
    return matches.any()

