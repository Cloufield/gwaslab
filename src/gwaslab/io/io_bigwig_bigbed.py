"""
bigWig and bigBed I/O utilities for gwaslab.

This module provides functions for reading and writing bigWig and bigBed format files
using the pyBigWig library. bigWig files store continuous-valued data (e.g., signal tracks),
while bigBed files store discrete intervals with associated data (e.g., annotations).

bigWig/bigBed format specification: https://genome.ucsc.edu/FAQ/FAQformat.html#format6.1

Coordinate system:
- bigWig and bigBed files use 0-based, half-open intervals [start, end)
- start: 0-based (first base in chromosome is numbered 0)
- end: exclusive (not included in the feature)
"""

from typing import TYPE_CHECKING, Union, Optional, List, Tuple, Dict, Any
import pandas as pd
import numpy as np
import os
from pathlib import Path
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

if TYPE_CHECKING:
    pass

# File suffix definitions
BIGWIG_SUFFIXES = ('.bw', '.bigwig', '.bigWig')
BIGBED_SUFFIXES = ('.bb', '.bigbed', '.bigBed')

# Try to import pyBigWig
try:
    import pyBigWig
    PYBigWig_AVAILABLE = True
except ImportError:
    PYBigWig_AVAILABLE = False
    import warnings
    warnings.warn(
        "pyBigWig not available. Install pyBigWig for bigWig/bigBed support: "
        "pip install pybigwig",
        UserWarning
    )


def _find_matching_chromosome(
    chrom: str,
    available_chroms: Dict[str, int],
    mapper: Optional[ChromosomeMapper] = None,
    log: Optional[Log] = None,
    verbose: bool = False
) -> Optional[str]:
    """
    Find matching chromosome name in available chromosomes using ChromosomeMapper.
    
    This function uses the same chromosome matching logic as BED files:
    - Strips "chr" or "CHR" prefixes
    - Uses ChromosomeMapper to convert to numeric and back
    - Handles case-insensitive matching
    
    Parameters
    ----------
    chrom : str
        Chromosome name to match (e.g., 'chr1', '1', 'Chr1')
    available_chroms : dict
        Dictionary of available chromosome names and sizes from bigWig/bigBed file
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance. If None, creates a new one.
    log : Log, optional
        Logger instance for messages
    verbose : bool, default=False
        Whether to show verbose messages
    
    Returns
    -------
    str or None
        Matching chromosome name from available_chroms, or None if not found
    """
    if mapper is None:
        mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    
    # First, try direct match
    if chrom in available_chroms:
        return chrom
    
    # Try case-insensitive direct match
    chrom_lower = chrom.lower()
    for c in available_chroms.keys():
        if c.lower() == chrom_lower:
            if log and verbose:
                log.write(f"  -Matched chromosome '{chrom}' to '{c}' (case-insensitive)", verbose=verbose)
            return c
    
    # Use ChromosomeMapper: convert input chrom to numeric, then try to match
    # Strip "chr" prefix and convert to numeric
    chrom_stripped = str(chrom).strip().lstrip("chr").lstrip("CHR")
    chrom_num = mapper.sumstats_to_number(chrom_stripped)
    
    if chrom_num is not None:
        # Try to find a chromosome in available_chroms that maps to the same number
        for c in available_chroms.keys():
            c_stripped = str(c).strip().lstrip("chr").lstrip("CHR")
            c_num = mapper.sumstats_to_number(c_stripped)
            if c_num == chrom_num:
                if log and verbose:
                    log.write(f"  -Matched chromosome '{chrom}' to '{c}' (via ChromosomeMapper)", verbose=verbose)
                return c
    
    # If still not found, return None
    return None


def read_bigwig(
    bw_path: str,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    as_dataframe: bool = True,
    numpy: bool = False,
    verbose: bool = False,
    log: Optional[Log] = None
) -> Union[pd.DataFrame, List[float], np.ndarray]:
    """
    Read values from a bigWig file.
    
    Parameters
    ----------
    bw_path : str
        Path to bigWig file (.bw, .bigwig, .bigWig)
    chrom : str, optional
        Chromosome name (e.g., 'chr1', '1'). If None, returns header information.
    start : int, optional
        Start position (0-based). Required if chrom is specified.
    end : int, optional
        End position (0-based, exclusive). Required if chrom is specified.
    as_dataframe : bool, default=True
        If True, returns a pandas DataFrame with columns: chrom, start, end, value.
        If False, returns a list of values (one per base position).
    numpy : bool, default=False
        If True and as_dataframe=False, returns a numpy array instead of list.
        Only works if pyBigWig was compiled with numpy support.
    verbose : bool, default=False
        If True, prints information about the operation.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    pd.DataFrame, list, or np.ndarray
        - If chrom/start/end specified and as_dataframe=True: DataFrame with columns
          (chrom, start, end, value) for each interval
        - If chrom/start/end specified and as_dataframe=False: List or array of values
          (one per base position, NaN for uncovered positions)
        - If chrom/start/end not specified: DataFrame with header information
          (chromosomes and their sizes)
    
    Notes
    -----
    bigWig files use 0-based, half-open intervals [start, end).
    The values() function returns one value per base position in the range.
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    if not os.path.exists(bw_path):
        raise FileNotFoundError(f"bigWig file not found: {bw_path}")
    
    if log and verbose:
        log.write(f" -Reading bigWig file: {bw_path}", verbose=verbose)
    
    bw = pyBigWig.open(bw_path)
    
    try:
        # If no chrom specified, return header information
        if chrom is None:
            if log and verbose:
                log.write(" -No chromosome specified, returning header information", verbose=verbose)
            
            chroms = bw.chroms()
            header_df = pd.DataFrame([
                {"chrom": chrom_name, "size": size}
                for chrom_name, size in chroms.items()
            ])
            return header_df
        
        # Validate start and end
        if start is None or end is None:
            raise ValueError("Both start and end must be specified when chrom is provided")
        
        if start >= end:
            raise ValueError(f"start ({start}) must be less than end ({end})")
        
        # Get chromosome size and validate/clamp coordinates
        chroms = bw.chroms()
        matching_chrom = _find_matching_chromosome(chrom, chroms, log=log, verbose=verbose)
        
        if matching_chrom is None:
            available_chroms = list(chroms.keys())[:10]  # Show first 10
            raise ValueError(
                f"Chromosome '{chrom}' not found in bigWig file. "
                f"Available chromosomes (showing first 10): {available_chroms}. "
                f"Total chromosomes: {len(chroms)}"
            )
        
        chrom = matching_chrom
        
        chrom_size = chroms[chrom]
        
        # Clamp coordinates to valid range
        start = max(0, start)
        end = min(chrom_size, end)
        
        if start >= end:
            if log:
                log.warning(f"  -Adjusted coordinates are invalid (start={start}, end={end}, chrom_size={chrom_size}). Returning empty result.")
            if as_dataframe:
                return pd.DataFrame(columns=["chrom", "start", "end", "value"])
            else:
                return [] if not numpy else np.array([])
        
        # Get values for the specified region
        values = bw.values(chrom, start, end, numpy=numpy)
        
        if as_dataframe:
            # Convert to DataFrame with intervals
            # Get intervals instead of per-base values for more efficient representation
            intervals = bw.intervals(chrom, start, end)
            
            if intervals is None or len(intervals) == 0:
                # No data in this region
                return pd.DataFrame(columns=["chrom", "start", "end", "value"])
            
            # Convert intervals to DataFrame
            # intervals is a list of tuples: [(start, end, value), ...]
            df = pd.DataFrame(intervals, columns=["start", "end", "value"])
            df.insert(0, "chrom", chrom)
            return df
        else:
            # Return raw values (list or numpy array)
            return values
    
    finally:
        bw.close()


def read_bigwig_intervals(
    bw_path: str,
    chrom: str,
    start: int,
    end: int,
    with_string: bool = True,
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Read intervals from a bigWig file for a specific region.
    
    Parameters
    ----------
    bw_path : str
        Path to bigWig file
    chrom : str
        Chromosome name
    start : int
        Start position (0-based)
    end : int
        End position (0-based, exclusive)
    with_string : bool, default=True
        If True, includes string values in intervals (for bigBed compatibility).
        For bigWig, this is typically not used.
    verbose : bool, default=False
        If True, prints information about the operation.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    if not os.path.exists(bw_path):
        raise FileNotFoundError(f"bigWig file not found: {bw_path}")
    
    if log and verbose:
        log.write(f" -Reading bigWig intervals from: {bw_path}", verbose=verbose)
    
    bw = pyBigWig.open(bw_path)
    
    try:
        intervals = bw.intervals(chrom, start, end)
        
        if intervals is None or len(intervals) == 0:
            return pd.DataFrame(columns=["chrom", "start", "end", "value"])
        
        df = pd.DataFrame(intervals, columns=["start", "end", "value"])
        df.insert(0, "chrom", chrom)
        return df
    
    finally:
        bw.close()


def read_bigwig_stats(
    bw_path: str,
    chrom: str,
    start: int,
    end: int,
    stat: str = "mean",
    n_bins: int = 1,
    verbose: bool = False,
    log: Optional[Log] = None
) -> Union[float, List[float]]:
    """
    Get statistics from a bigWig file for a region.
    
    Parameters
    ----------
    bw_path : str
        Path to bigWig file
    chrom : str
        Chromosome name
    start : int
        Start position (0-based)
    end : int
        End position (0-based, exclusive)
    stat : str, default='mean'
        Statistic to compute. Options: 'mean', 'std', 'min', 'max', 'coverage', 'sum'
    n_bins : int, default=1
        Number of bins to divide the region into. If 1, returns a single value.
        If > 1, returns a list of values (one per bin).
    verbose : bool, default=False
        If True, prints information about the operation.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    float or list of float
        Statistic value(s) for the region
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    if not os.path.exists(bw_path):
        raise FileNotFoundError(f"bigWig file not found: {bw_path}")
    
    if log and verbose:
        log.write(f" -Computing {stat} from bigWig file: {bw_path}", verbose=verbose)
    
    bw = pyBigWig.open(bw_path)
    
    try:
        stats = bw.stats(chrom, start, end, type=stat, nBins=n_bins)
        
        if stats is None:
            return None if n_bins == 1 else []
        
        if n_bins == 1:
            return stats[0] if len(stats) > 0 else None
        else:
            return stats
    
    finally:
        bw.close()


def write_bigwig(
    output_path: str,
    intervals_df: Optional[pd.DataFrame] = None,
    chroms: Optional[List[Tuple[str, int]]] = None,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    value_col: str = "value",
    max_zooms: int = 10,
    validate: bool = True,
    verbose: bool = False,
    log: Optional[Log] = None
) -> str:
    """
    Write a bigWig file from a DataFrame of intervals.
    
    Parameters
    ----------
    output_path : str
        Path to output bigWig file (.bw, .bigwig, .bigWig)
    intervals_df : pd.DataFrame, optional
        DataFrame with columns: chrom, start, end, value.
        If None, only creates header (requires chroms parameter).
    chroms : list of tuples, optional
        List of (chromosome_name, size) tuples for header.
        Required if intervals_df doesn't contain all chromosomes.
        Example: [("chr1", 1000000), ("chr2", 1500000)]
    chrom_col : str, default='chrom'
        Column name for chromosome in intervals_df
    start_col : str, default='start'
        Column name for start position in intervals_df
    end_col : str, default='end'
        Column name for end position in intervals_df
    value_col : str, default='value'
        Column name for value in intervals_df
    max_zooms : int, default=10
        Maximum number of zoom levels to create. Set to 0 to disable zoom levels
        (not recommended as IGV and other tools require at least one zoom level).
    validate : bool, default=True
        If True, validates that entries are added in order.
        Set to False for faster writing if you're sure entries are ordered.
    verbose : bool, default=False
        If True, prints information about the operation.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    str
        Path to the created bigWig file
    
    Notes
    -----
    Entries must be added in order (sorted by chromosome, then by start position).
    The function will automatically sort if validate=True, but this adds overhead.
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    if log and verbose:
        log.write(f" -Writing bigWig file: {output_path}", verbose=verbose)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    bw = pyBigWig.open(output_path, "w")
    
    try:
        # Determine chromosomes and sizes
        if intervals_df is not None and len(intervals_df) > 0:
            # Get unique chromosomes and their max positions
            chrom_sizes = {}
            for chrom in intervals_df[chrom_col].unique():
                chrom_data = intervals_df[intervals_df[chrom_col] == chrom]
                max_end = chrom_data[end_col].max()
                chrom_sizes[chrom] = int(max_end)
            
            # Merge with provided chroms if any
            if chroms is not None:
                for chrom_name, size in chroms:
                    if chrom_name in chrom_sizes:
                        chrom_sizes[chrom_name] = max(chrom_sizes[chrom_name], size)
                    else:
                        chrom_sizes[chrom_name] = size
            
            header = [(chrom, size) for chrom, size in sorted(chrom_sizes.items())]
        elif chroms is not None:
            header = sorted(chroms)
        else:
            raise ValueError("Either intervals_df or chroms must be provided")
        
        # Add header
        bw.addHeader(header, maxZooms=max_zooms)
        
        if log and verbose:
            log.write(f" -Added header with {len(header)} chromosomes", verbose=verbose)
        
        # Add entries if provided
        if intervals_df is not None and len(intervals_df) > 0:
            # Sort by chromosome and start position
            if validate:
                intervals_df = intervals_df.sort_values([chrom_col, start_col]).copy()
            
            # Group by chromosome and add entries
            for chrom_name in sorted(intervals_df[chrom_col].unique()):
                chrom_data = intervals_df[intervals_df[chrom_col] == chrom_name].copy()
                
                chroms_list = chrom_data[chrom_col].tolist()
                starts = chrom_data[start_col].astype(int).tolist()
                ends = chrom_data[end_col].astype(int).tolist()
                values = chrom_data[value_col].astype(float).tolist()
                
                bw.addEntries(chroms_list, starts, ends=ends, values=values, validate=validate)
            
            if log and verbose:
                log.write(f" -Added {len(intervals_df)} intervals", verbose=verbose)
    
    finally:
        bw.close()
    
    if log and verbose:
        log.write(f" -Successfully created bigWig file: {output_path}", verbose=verbose)
    
    return output_path


def read_bigbed(
    bb_path: str,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    with_string: bool = True,
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Read entries from a bigBed file.
    
    Parameters
    ----------
    bb_path : str
        Path to bigBed file (.bb, .bigbed, .bigBed)
    chrom : str, optional
        Chromosome name. If None, returns header information.
    start : int, optional
        Start position (0-based). Required if chrom is specified.
    end : int, optional
        End position (0-based, exclusive). Required if chrom is specified.
    with_string : bool, default=True
        If True, includes string values associated with each entry.
        If False, returns only coordinates (saves memory).
    verbose : bool, default=False
        If True, prints information about the operation.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    pd.DataFrame
        - If chrom/start/end specified: DataFrame with columns:
          - chrom, start, end (always present)
          - rest (string values) if with_string=True
        - If chrom/start/end not specified: DataFrame with header information
          (chromosomes and their sizes)
    
    Notes
    -----
    bigBed files use 0-based, half-open intervals [start, end).
    The 'rest' column contains the remaining fields from the BED entry.
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    if not os.path.exists(bb_path):
        raise FileNotFoundError(f"bigBed file not found: {bb_path}")
    
    if log and verbose:
        log.write(f" -Reading bigBed file: {bb_path}", verbose=verbose)
    
    bb = pyBigWig.open(bb_path)
    
    try:
        # If no chrom specified, return header information
        if chrom is None:
            if log and verbose:
                log.write(" -No chromosome specified, returning header information", verbose=verbose)
            
            chroms = bb.chroms()
            header_df = pd.DataFrame([
                {"chrom": chrom_name, "size": size}
                for chrom_name, size in chroms.items()
            ])
            return header_df
        
        # Validate start and end
        if start is None or end is None:
            raise ValueError("Both start and end must be specified when chrom is provided")
        
        if start >= end:
            raise ValueError(f"start ({start}) must be less than end ({end})")
        
        # Get chromosome size and validate/clamp coordinates
        chroms = bb.chroms()
        matching_chrom = _find_matching_chromosome(chrom, chroms, log=log, verbose=verbose)
        
        if matching_chrom is None:
            available_chroms = list(chroms.keys())[:10]  # Show first 10
            raise ValueError(
                f"Chromosome '{chrom}' not found in bigBed file. "
                f"Available chromosomes (showing first 10): {available_chroms}. "
                f"Total chromosomes: {len(chroms)}"
            )
        
        chrom = matching_chrom
        
        chrom_size = chroms[chrom]
        
        # Clamp coordinates to valid range
        start = max(0, start)
        end = min(chrom_size, end)
        
        if start >= end:
            if log:
                log.warning(f"  -Adjusted coordinates are invalid (start={start}, end={end}, chrom_size={chrom_size}). Returning empty DataFrame.")
            return pd.DataFrame(columns=["chrom", "start", "end"])
        
        # Get entries for the specified region
        entries = bb.entries(chrom, start, end, withString=with_string)
        
        if entries is None or len(entries) == 0:
            return pd.DataFrame(columns=["chrom", "start", "end"])
        
        # Convert entries to DataFrame
        # entries is a list of tuples: [(start, end, rest), ...] or [(start, end), ...]
        if with_string and len(entries[0]) == 3:
            df = pd.DataFrame(entries, columns=["start", "end", "rest"])
        else:
            df = pd.DataFrame(entries, columns=["start", "end"])
        
        df.insert(0, "chrom", chrom)
        return df
    
    finally:
        bb.close()


def write_bigbed(
    output_path: str,
    intervals_df: pd.DataFrame,
    chroms: List[Tuple[str, int]],
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    rest_col: Optional[str] = None,
    max_zooms: int = 10,
    validate: bool = True,
    verbose: bool = False,
    log: Optional[Log] = None
) -> str:
    """
    Write a bigBed file from a DataFrame of intervals.
    
    Parameters
    ----------
    output_path : str
        Path to output bigBed file (.bb, .bigbed, .bigBed)
    intervals_df : pd.DataFrame
        DataFrame with columns: chrom, start, end, and optionally rest (string values).
        The 'rest' column should contain the remaining BED fields as strings.
    chroms : list of tuples
        List of (chromosome_name, size) tuples for header.
        Example: [("chr1", 1000000), ("chr2", 1500000)]
    chrom_col : str, default='chrom'
        Column name for chromosome in intervals_df
    start_col : str, default='start'
        Column name for start position in intervals_df
    end_col : str, default='end'
        Column name for end position in intervals_df
    rest_col : str, optional
        Column name for rest (string values) in intervals_df.
        If None, no string values are included.
    max_zooms : int, default=10
        Maximum number of zoom levels to create.
    validate : bool, default=True
        If True, validates that entries are added in order.
    verbose : bool, default=False
        If True, prints information about the operation.
    log : Log, optional
        Logger instance for messages.
    
    Returns
    -------
    str
        Path to the created bigBed file
    
    Notes
    -----
    Entries must be added in order (sorted by chromosome, then by start position).
    """
    if not PYBigWig_AVAILABLE:
        raise ImportError(
            "pyBigWig is not installed. Install it with: pip install pybigwig"
        )
    
    if log and verbose:
        log.write(f" -Writing bigBed file: {output_path}", verbose=verbose)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    bb = pyBigWig.open(output_path, "w")
    
    try:
        # Add header
        header = sorted(chroms)
        bb.addHeader(header, maxZooms=max_zooms)
        
        if log and verbose:
            log.write(f" -Added header with {len(header)} chromosomes", verbose=verbose)
        
        # Sort by chromosome and start position
        if validate:
            intervals_df = intervals_df.sort_values([chrom_col, start_col]).copy()
        
        # Group by chromosome and add entries
        for chrom_name in sorted(intervals_df[chrom_col].unique()):
            chrom_data = intervals_df[intervals_df[chrom_col] == chrom_name].copy()
            
            chroms_list = chrom_data[chrom_col].tolist()
            starts = chrom_data[start_col].astype(int).tolist()
            ends = chrom_data[end_col].astype(int).tolist()
            
            if rest_col is not None and rest_col in chrom_data.columns:
                rests = chrom_data[rest_col].astype(str).tolist()
                bb.addEntries(chroms_list, starts, ends=ends, rest=rests, validate=validate)
            else:
                bb.addEntries(chroms_list, starts, ends=ends, validate=validate)
        
        if log and verbose:
            log.write(f" -Added {len(intervals_df)} entries", verbose=verbose)
    
    finally:
        bb.close()
    
    if log and verbose:
        log.write(f" -Successfully created bigBed file: {output_path}", verbose=verbose)
    
    return output_path


def bigwig_to_sumstats_coordinates(
    bw_df: pd.DataFrame,
    mapper: Optional[ChromosomeMapper] = None,
    chrom_col: str = 'chrom',
    start_col: str = 'start',
    end_col: str = 'end',
    value_col: str = 'value',
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Convert bigWig coordinates to sumstats-compatible format.
    
    Converts:
    - Chromosome notation: bigWig chrom → sumstats CHR format
    - Coordinates: bigWig 0-based [start, end) → 1-based [start+1, end] for matching
    
    Parameters
    ----------
    bw_df : pd.DataFrame
        bigWig DataFrame with chrom, start, end, value columns
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome conversion.
        If None, creates a new mapper.
    chrom_col : str, default='chrom'
        Column name for chromosome in bigWig DataFrame
    start_col : str, default='start'
        Column name for start position in bigWig DataFrame
    end_col : str, default='end'
        Column name for end position in bigWig DataFrame
    value_col : str, default='value'
        Column name for value in bigWig DataFrame
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
        - VALUE: Value from bigWig (preserved)
        Original bigWig columns are preserved with '_bw' suffix.
    
    Notes
    -----
    bigWig uses 0-based, half-open intervals [start, end).
    For matching with 1-based sumstats POS:
    - bigWig [start, end) in 0-based = [start+1, end] in 1-based
    - We convert to: START = start + 1, END = end (both inclusive for matching)
    """
    if mapper is None:
        mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    
    # Create a copy to avoid modifying original
    result_df = bw_df.copy()
    
    # Convert chromosome notation to numeric (middle layer)
    bw_chr_numeric = result_df[chrom_col].apply(
        lambda x: mapper.sumstats_to_number(str(x).strip().lstrip("chr").lstrip("CHR"))
        if pd.notna(x) else None
    )
    
    # Filter out unconvertible chromosomes
    valid_mask = bw_chr_numeric.notna()
    n_invalid = (~valid_mask).sum()
    
    if n_invalid > 0 and log:
        log.write(f" -Warning: {n_invalid} bigWig intervals have unconvertible chromosome notation.", verbose=verbose)
    
    result_df = result_df[valid_mask].copy()
    bw_chr_numeric = bw_chr_numeric[valid_mask]
    
    # Convert numeric chromosomes to sumstats format
    result_df['CHR'] = bw_chr_numeric.apply(
        lambda x: mapper.number_to_sumstats(x) if pd.notna(x) else None
    )
    
    # Convert bigWig coordinates to 1-based for matching
    # bigWig [start, end) in 0-based = [start+1, end] in 1-based
    # For matching, we use: START = start + 1, END = end (both inclusive)
    result_df['START'] = result_df[start_col] + 1  # Convert 0-based to 1-based
    result_df['END'] = result_df[end_col]  # End is already 1-based (exclusive in bigWig = inclusive for matching)
    
    # Preserve value column
    if value_col in result_df.columns:
        result_df['VALUE'] = result_df[value_col]
    
    # Rename original columns with _bw suffix
    rename_dict = {chrom_col: f'{chrom_col}_bw', start_col: f'{start_col}_bw', end_col: f'{end_col}_bw'}
    if value_col in result_df.columns and value_col != 'VALUE':
        rename_dict[value_col] = f'{value_col}_bw'
    result_df = result_df.rename(columns=rename_dict)
    
    return result_df


def bigbed_to_sumstats_coordinates(
    bb_df: pd.DataFrame,
    mapper: Optional[ChromosomeMapper] = None,
    chrom_col: str = 'chrom',
    start_col: str = 'start',
    end_col: str = 'end',
    verbose: bool = False,
    log: Optional[Log] = None
) -> pd.DataFrame:
    """
    Convert bigBed coordinates to sumstats-compatible format.
    
    Converts:
    - Chromosome notation: bigBed chrom → sumstats CHR format
    - Coordinates: bigBed 0-based [start, end) → 1-based [start+1, end] for matching
    
    Parameters
    ----------
    bb_df : pd.DataFrame
        bigBed DataFrame with chrom, start, end columns (and optionally rest)
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome conversion.
        If None, creates a new mapper.
    chrom_col : str, default='chrom'
        Column name for chromosome in bigBed DataFrame
    start_col : str, default='start'
        Column name for start position in bigBed DataFrame
    end_col : str, default='end'
        Column name for end position in bigBed DataFrame
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
        - REST: Rest column from bigBed (if present, preserved)
        Original bigBed columns are preserved with '_bb' suffix.
    
    Notes
    -----
    bigBed uses 0-based, half-open intervals [start, end).
    For matching with 1-based sumstats POS:
    - bigBed [start, end) in 0-based = [start+1, end] in 1-based
    - We convert to: START = start + 1, END = end (both inclusive for matching)
    """
    if mapper is None:
        mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    
    # Create a copy to avoid modifying original
    result_df = bb_df.copy()
    
    # Convert chromosome notation to numeric (middle layer)
    bb_chr_numeric = result_df[chrom_col].apply(
        lambda x: mapper.sumstats_to_number(str(x).strip().lstrip("chr").lstrip("CHR"))
        if pd.notna(x) else None
    )
    
    # Filter out unconvertible chromosomes
    valid_mask = bb_chr_numeric.notna()
    n_invalid = (~valid_mask).sum()
    
    if n_invalid > 0 and log:
        log.write(f" -Warning: {n_invalid} bigBed entries have unconvertible chromosome notation.", verbose=verbose)
    
    result_df = result_df[valid_mask].copy()
    bb_chr_numeric = bb_chr_numeric[valid_mask]
    
    # Convert numeric chromosomes to sumstats format
    result_df['CHR'] = bb_chr_numeric.apply(
        lambda x: mapper.number_to_sumstats(x) if pd.notna(x) else None
    )
    
    # Convert bigBed coordinates to 1-based for matching
    # bigBed [start, end) in 0-based = [start+1, end] in 1-based
    # For matching, we use: START = start + 1, END = end (both inclusive)
    result_df['START'] = result_df[start_col] + 1  # Convert 0-based to 1-based
    result_df['END'] = result_df[end_col]  # End is already 1-based (exclusive in bigBed = inclusive for matching)
    
    # Preserve rest column if present
    if 'rest' in result_df.columns:
        result_df['REST'] = result_df['rest']
    
    # Rename original columns with _bb suffix
    rename_dict = {chrom_col: f'{chrom_col}_bb', start_col: f'{start_col}_bb', end_col: f'{end_col}_bb'}
    if 'rest' in result_df.columns and 'rest' != 'REST':
        rename_dict['rest'] = 'rest_bb'
    result_df = result_df.rename(columns=rename_dict)
    
    return result_df
