"""
Simplified FASTA I/O module for gwaslab.

This module provides simplified functions for reading and writing FASTA format files.
Uses pysam FastxFile for fast FASTA reading (3-4x faster than previous implementation).
"""

import gzip
import numpy as np
import pandas as pd
from typing import Iterator, TextIO, Union, Dict, Tuple, Optional
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_common_data import _maketrans
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

# FASTA file suffix definitions
FASTA_SUFFIXES = ('.fa.gz', '.fasta.gz', '.fa.bgz', '.fasta.bgz', '.fa', '.fasta')

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    import warnings
    warnings.warn(
        "pysam not available. Falling back to slower implementation. "
        "Install pysam for better performance: pip install pysam",
        UserWarning
    )

# Constants for FASTA sequence translation
# chr(0) should not be used in the mapping dict because it's a reserved value.
# Instead of starting from chr(1), we start from chr(2) because this could be useful in the future
# to compute the complementary allele with a simple XOR operation (e.g. 2 ^ 1 = 3, 3 ^ 1 = 2, 4 ^ 1 = 5, 5 ^ 1 = 4, ...)
_FASTA_MAPPING = {
    "A": chr(2),
    "T": chr(3),
    "C": chr(4),
    "G": chr(5),
    "N": chr(6),
}
assert all(value != chr(0) for value in _FASTA_MAPPING.values()), "Mapping in the dictionary should not be equal to chr(0). This is a reserved value"
_FASTA_TRANSLATE_TABLE = _maketrans(_FASTA_MAPPING)


class FastaRecord:
    """
    Simple FASTA record class.
    
    This class provides compatibility with code that expects FASTA record objects
    with .id and .seq._data attributes.
    """
    def __init__(self, id: str, seq: str):
        """
        Initialize a FastaRecord.
        
        Parameters
        ----------
        id : str
            Record identifier (title line without '>')
        seq : str
            Sequence string
        """
        self.id = id
        self.seq = _Sequence(seq)
    
    def __repr__(self):
        return f"FastaRecord(id='{self.id}', length={len(self.seq._data)})"


class _Sequence:
    """
    Internal class to provide sequence with _data attribute.
    """
    def __init__(self, seq: str):
        # Store as bytes for compatibility with translate() operations
        if isinstance(seq, str):
            self._data = seq.encode('ascii')
        else:
            self._data = seq
    
    def __len__(self):
        return len(self._data)


def _open_fasta_handle(path: str) -> TextIO:
    """
    Open a FASTA file handle, automatically detecting compression.
    
    Parameters
    ----------
    path : str
        Path to FASTA file. Supports .fa, .fasta, .fa.gz, .fasta.gz
        Note: .bgz files are treated as .gz files (using gzip decompression)
        
    Returns
    -------
    TextIO
        File handle opened in text mode
        
    Raises
    ------
    ValueError
        If file extension is not recognized
    """
    path_lower = path.lower()
    
    # BGZF (bgzip) detection - treat as gzip
    # Note: BGZF is a variant of gzip, so gzip can often read it
    if path_lower.endswith((".bgz", ".bgzf")):
        import warnings
        warnings.warn(
            "BGZF file detected. Using gzip decompression. "
            "Some BGZF-specific features may not be supported.",
            UserWarning
        )
        return gzip.open(path, "rt")
    
    # Standard gzipped FASTA
    elif path_lower.endswith(".gz"):
        return gzip.open(path, "rt")
    
    # Plain FASTA
    elif path_lower.endswith((".fa", ".fasta")):
        return open(path, "r")
    
    else:
        raise ValueError(
            f"Unrecognized FASTA file extension: {path}. "
            "Supported extensions: .fa, .fasta, .fa.gz, .fasta.gz, .fa.bgz, .fasta.bgz"
        )


def parse_fasta_simple(handle: TextIO) -> Iterator[Tuple[str, str]]:
    """
    Iterate over FASTA records as string tuples (title, sequence).
    
    This is a simplified FASTA parser (fallback implementation when pysam is not available).
    For each record, returns a tuple of (title, sequence) where:
    - title: The FASTA title line without the leading '>' character
    - sequence: The sequence with whitespace removed
    
    Parameters
    ----------
    handle : TextIO
        Input stream opened in text mode
        
    Yields
    ------
    tuple[str, str]
        Tuple of (title, sequence) for each FASTA record
        
    Examples
    --------
    >>> with open("example.fasta") as handle:
    ...     for title, seq in parse_fasta_simple(handle):
    ...         print(f"{title}: {seq[:10]}...")
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    title = None
    for line in handle:
        if line.startswith(">"):
            title = line[1:].rstrip()
            break
    else:
        # No break encountered - probably an empty file
        return
    
    # Main logic
    # Note: remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files)
    lines = []
    for line in handle:
        if line.startswith(">"):
            sequence = "".join(lines).replace(" ", "").replace("\r", "")
            yield title, sequence
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    
    # Yield the last record
    if title is not None:
        sequence = "".join(lines).replace(" ", "").replace("\r", "")
        yield title, sequence


def parse_fasta(path: str, as_dict: bool = True) -> Union[Dict[str, str], Iterator[Tuple[str, str]]]:
    """
    Parse a FASTA file and return records.
    
    Uses pysam FastxFile for fast reading (3-4x faster than previous implementation).
    Falls back to slower implementation if pysam is not available.
    
    Parameters
    ----------
    path : str
        Path to FASTA file. Supports plain and gzipped formats.
    as_dict : bool, optional
        If True (default), return a dictionary mapping titles to sequences.
        If False, return an iterator of (title, sequence) tuples.
        Note: When as_dict=False, the file handle is kept open until iteration completes.
        
    Returns
    -------
    dict[str, str] or Iterator[tuple[str, str]]
        If as_dict=True: Dictionary mapping record titles to sequences.
        If as_dict=False: Iterator yielding (title, sequence) tuples.
        
    Examples
    --------
    >>> # Load as dictionary
    >>> records = parse_fasta("reference.fasta")
    >>> print(records["chr1"][:10])
    
    >>> # Load as iterator (memory efficient for large files)
    >>> for title, seq in parse_fasta("reference.fasta", as_dict=False):
    ...     print(f"{title}: {len(seq)} bp")
    """
    if PYSAM_AVAILABLE:
        # Use pysam FastxFile for fast reading
        if as_dict:
            records = {}
            try:
                with pysam.FastxFile(path) as fastx:
                    for entry in fastx:
                        records[entry.name] = entry.sequence
                return records
            except Exception as e:
                # Fallback to old implementation if pysam fails
                import warnings
                warnings.warn(
                    f"pysam FastxFile failed, falling back to slower implementation: {e}",
                    UserWarning
                )
                with _open_fasta_handle(path) as handle:
                    return dict(parse_fasta_simple(handle))
        else:
            # Iterator mode with pysam
            def _parse_generator():
                try:
                    with pysam.FastxFile(path) as fastx:
                        for entry in fastx:
                            yield entry.name, entry.sequence
                except Exception as e:
                    # Fallback to old implementation if pysam fails
                    import warnings
                    warnings.warn(
                        f"pysam FastxFile failed, falling back to slower implementation: {e}",
                        UserWarning
                    )
                    handle = _open_fasta_handle(path)
                    try:
                        yield from parse_fasta_simple(handle)
                    finally:
                        handle.close()
            
            return _parse_generator()
    else:
        # Fallback to old implementation
        if as_dict:
            with _open_fasta_handle(path) as handle:
                return dict(parse_fasta_simple(handle))
        else:
            def _parse_generator():
                handle = _open_fasta_handle(path)
                try:
                    yield from parse_fasta_simple(handle)
                finally:
                    handle.close()
            
            return _parse_generator()


def load_fasta_auto(path: str, as_seqrecord: bool = True):
    """
    Automatically load FASTA or gzipped FASTA files.
    
    This function is compatible with the existing load_fasta_auto in hm_harmonize_sumstats.py.
    Returns an iterator of FastaRecord objects or (title, sequence) tuples.
    Uses pysam FastxFile for fast reading (3-4x faster than previous implementation).
    
    Note: The file handle is kept open until iteration completes. For proper resource
    management, ensure you consume the entire iterator or use it within a context manager.
    
    Parameters
    ----------
    path : str
        Path to FASTA file. Supports: .fa, .fasta, .fa.gz, .fasta.gz, .fa.bgz, .fasta.bgz
    as_seqrecord : bool, optional
        If True (default), return FastaRecord objects with .id and .seq._data attributes.
        If False, return (title, sequence) tuples.
        
    Returns
    -------
    Iterator[FastaRecord] or Iterator[tuple[str, str]]
        Iterator yielding FastaRecord objects or (title, sequence) tuples.
        
    Examples
    --------
    >>> # Get FastaRecord objects
    >>> for record in load_fasta_auto("reference.fasta.gz"):
    ...     print(f"{record.id}: {len(record.seq._data)} bp")
    
    >>> # Get simple tuples
    >>> for title, seq in load_fasta_auto("reference.fasta.gz", as_seqrecord=False):
    ...     print(f"{title}: {len(seq)} bp")
    """
    def _load_generator():
        if PYSAM_AVAILABLE:
            try:
                with pysam.FastxFile(path) as fastx:
                    for entry in fastx:
                        if as_seqrecord:
                            yield FastaRecord(entry.name, entry.sequence)
                        else:
                            yield entry.name, entry.sequence
            except Exception as e:
                # Fallback to old implementation if pysam fails
                import warnings
                warnings.warn(
                    f"pysam FastxFile failed, falling back to slower implementation: {e}",
                    UserWarning
                )
                handle = _open_fasta_handle(path)
                try:
                    for title, sequence in parse_fasta_simple(handle):
                        if as_seqrecord:
                            yield FastaRecord(title, sequence)
                        else:
                            yield title, sequence
                finally:
                    handle.close()
        else:
            # Fallback to old implementation
            handle = _open_fasta_handle(path)
            try:
                for title, sequence in parse_fasta_simple(handle):
                    if as_seqrecord:
                        yield FastaRecord(title, sequence)
                    else:
                        yield title, sequence
            finally:
                handle.close()
    
    return _load_generator()


def load_fasta_filtered(
    path: str,
    chromlist_set: set,
    chroms_in_sumstats_set: set,
    mapper: Optional[ChromosomeMapper] = None,
    log: Log = Log(),
    verbose: bool = True
) -> Dict[str, FastaRecord]:
    """
    Load and filter FASTA records in a single pass for better performance.
    
    This function combines loading and filtering to avoid creating FastaRecord objects
    for records that will be filtered out. Only creates records for chromosomes that
    are needed.
    Uses pysam FastxFile for fast reading (3-4x faster than previous implementation).
    
    Parameters
    ----------
    path : str
        Path to FASTA file. Supports: .fa, .fasta, .fa.gz, .fasta.gz, .fa.bgz, .fasta.bgz
    chromlist_set : set
        Set of valid chromosome identifiers
    chroms_in_sumstats_set : set
        Set of chromosomes present in the summary statistics
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided, creates a default mapper.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object
    verbose : bool, default=True
        If True, print progress messages
        
    Returns
    -------
    dict[str, FastaRecord]
        Dictionary mapping chromosome identifiers to FastaRecord objects
    """
    # Create mapper if not provided
    if mapper is None:
        mapper = ChromosomeMapper(log=log, verbose=verbose)
    all_records_dict = {}
    
    if PYSAM_AVAILABLE:
        try:
            with pysam.FastxFile(path) as fastx:
                for entry in fastx:
                    title = entry.name
                    sequence = entry.sequence
                    
                    # Filter during loading - avoid creating FastaRecord if not needed
                    record_chr = title.strip()
                    
                    # Use mapper to convert to numeric format
                    try:
                        if mapper._sumstats_format is None:
                            mapper.detect_sumstats_format(pd.Series([record_chr]))
                        i = mapper.to_numeric(record_chr)
                        if isinstance(i, str):
                            try:
                                i = int(i)
                            except ValueError:
                                pass
                    except (KeyError, ValueError, AttributeError):
                        # Fallback: try stripping chr prefix
                        record_chr_stripped = record_chr.strip("chrCHR").upper()
                        try:
                            if mapper._sumstats_format is None:
                                mapper.detect_sumstats_format(pd.Series([record_chr_stripped]))
                            i = mapper.to_numeric(record_chr_stripped)
                            if isinstance(i, str):
                                try:
                                    i = int(i)
                                except ValueError:
                                    pass
                        except (KeyError, ValueError, AttributeError):
                            i = record_chr_stripped
                    
                    # Only create FastaRecord if it passes filters
                    if (i in chromlist_set) and (i in chroms_in_sumstats_set):
                        log.write(record_chr, " ", end="", show_time=False, verbose=verbose)
                        all_records_dict[i] = FastaRecord(title, sequence)
            return all_records_dict
        except Exception as e:
            # Fallback to old implementation if pysam fails
            import warnings
            warnings.warn(
                f"pysam FastxFile failed, falling back to slower implementation: {e}",
                UserWarning
            )
    
    # Fallback to old implementation
    handle = _open_fasta_handle(path)
    try:
        # Skip any text before the first record
        title = None
        for line in handle:
            if line.startswith(">"):
                title = line[1:].rstrip()
                break
        else:
            # No records found
            return all_records_dict
        
        # Main parsing and filtering loop
        lines = []
        for line in handle:
            if line.startswith(">"):
                # Process the previous record before starting a new one
                if title is not None:
                    sequence = "".join(lines).replace(" ", "").replace("\r", "")
                    # Filter during loading - avoid creating FastaRecord if not needed
                    record_chr = title.strip()
                    
                    # Use mapper to convert to numeric format
                    try:
                        if mapper._sumstats_format is None:
                            mapper.detect_sumstats_format(pd.Series([record_chr]))
                        i = mapper.to_numeric(record_chr)
                        if isinstance(i, str):
                            try:
                                i = int(i)
                            except ValueError:
                                pass
                    except (KeyError, ValueError, AttributeError):
                        # Fallback: try stripping chr prefix
                        record_chr_stripped = record_chr.strip("chrCHR").upper()
                        try:
                            if mapper._sumstats_format is None:
                                mapper.detect_sumstats_format(pd.Series([record_chr_stripped]))
                            i = mapper.to_numeric(record_chr_stripped)
                            if isinstance(i, str):
                                try:
                                    i = int(i)
                                except ValueError:
                                    pass
                        except (KeyError, ValueError, AttributeError):
                            i = record_chr_stripped
                    
                    # Only create FastaRecord if it passes filters
                    if (i in chromlist_set) and (i in chroms_in_sumstats_set):
                        log.write(record_chr, " ", end="", show_time=False, verbose=verbose)
                        all_records_dict[i] = FastaRecord(title, sequence)
                
                # Start new record
                lines = []
                title = line[1:].rstrip()
                continue
            lines.append(line.rstrip())
        
        # Process the last record
        if title is not None:
            sequence = "".join(lines).replace(" ", "").replace("\r", "")
            record_chr = title.strip()
            
            # Use mapper to convert to numeric format
            try:
                if mapper._sumstats_format is None:
                    mapper.detect_sumstats_format(pd.Series([record_chr]))
                i = mapper.to_numeric(record_chr)
                if isinstance(i, str):
                    try:
                        i = int(i)
                    except ValueError:
                        pass
            except (KeyError, ValueError, AttributeError):
                # Fallback: try stripping chr prefix
                record_chr_stripped = record_chr.strip("chrCHR").upper()
                try:
                    if mapper._sumstats_format is None:
                        mapper.detect_sumstats_format(pd.Series([record_chr_stripped]))
                    i = mapper.to_numeric(record_chr_stripped)
                    if isinstance(i, str):
                        try:
                            i = int(i)
                        except ValueError:
                            pass
                except (KeyError, ValueError, AttributeError):
                    i = record_chr_stripped
            
            if (i in chromlist_set) and (i in chroms_in_sumstats_set):
                log.write(record_chr, " ", end="", show_time=False, verbose=verbose)
                all_records_dict[i] = FastaRecord(title, sequence)
    
    finally:
        handle.close()
    
    return all_records_dict


def write_fasta(
    records: Union[Dict[str, str], Iterator[Tuple[str, str]]],
    path: str,
    wrap: int = 60,
    mode: str = "w"
) -> None:
    """
    Write FASTA records to a file.
    
    Parameters
    ----------
    records : dict[str, str] or Iterator[tuple[str, str]]
        FASTA records to write. Can be:
        - Dictionary mapping titles to sequences
        - Iterator of (title, sequence) tuples
    path : str
        Output file path. Supports plain and gzipped formats based on extension.
    wrap : int, optional
        Line length for sequence wrapping. Default is 60. Use 0 or None for no wrapping.
    mode : str, optional
        File open mode. Default is "w" (write). Use "a" for append.
        
    Examples
    --------
    >>> # Write from dictionary
    >>> records = {"chr1": "ATCGATCG", "chr2": "GCTAGCTA"}
    >>> write_fasta(records, "output.fasta")
    
    >>> # Write from iterator
    >>> records = [("chr1", "ATCGATCG"), ("chr2", "GCTAGCTA")]
    >>> write_fasta(records, "output.fasta.gz", wrap=80)
    """
    # Determine if we should compress based on extension
    path_lower = path.lower()
    if path_lower.endswith((".bgz", ".bgzf")):
        # Treat BGZF as gzip for writing
        import warnings
        warnings.warn(
            "BGZF file extension detected. Using gzip compression. "
            "Output will be standard gzip format, not BGZF.",
            UserWarning
        )
        handle = gzip.open(path, f"{mode}t")
    elif path_lower.endswith(".gz"):
        handle = gzip.open(path, f"{mode}t")
    else:
        handle = open(path, mode)
    
    try:
        # Handle dictionary input
        if isinstance(records, dict):
            records = records.items()
        
        # Write each record
        for title, sequence in records:
            # Clean title (remove newlines)
            title = title.replace("\n", "").replace("\r", "")
            
            # Write title line
            handle.write(f">{title}\n")
            
            # Write sequence with optional wrapping
            if wrap and wrap > 0:
                for i in range(0, len(sequence), wrap):
                    handle.write(sequence[i:i + wrap] + "\n")
            else:
                handle.write(sequence + "\n")
    
    finally:
        handle.close()


def get_fasta_record(path: str, title: str) -> Optional[str]:
    """
    Get a specific FASTA record by title.
    
    Uses pysam FastxFile for fast reading. For indexed files, could use pysam FastaFile
    for even faster random access, but FastxFile is sufficient for single record lookup.
    
    Parameters
    ----------
    path : str
        Path to FASTA file
    title : str
        Title of the record to retrieve (without the '>' character)
        
    Returns
    -------
    str or None
        Sequence for the requested record, or None if not found.
        
    Examples
    --------
    >>> seq = get_fasta_record("reference.fasta", "chr1")
    >>> if seq:
    ...     print(f"chr1 length: {len(seq)}")
    """
    if PYSAM_AVAILABLE:
        try:
            with pysam.FastxFile(path) as fastx:
                for entry in fastx:
                    if entry.name == title or entry.name.split()[0] == title:
                        return entry.sequence
            return None
        except Exception as e:
            # Fallback to old implementation if pysam fails
            import warnings
            warnings.warn(
                f"pysam FastxFile failed, falling back to slower implementation: {e}",
                UserWarning
            )
    
    # Fallback to old implementation
    with _open_fasta_handle(path) as handle:
        for record_title, sequence in parse_fasta_simple(handle):
            if record_title == title or record_title.split()[0] == title:
                return sequence
    return None


def load_and_build_fasta_records(
    path: str,
    chromlist_set: set,
    chroms_in_sumstats_set: set,
    mapper: Optional[ChromosomeMapper] = None,
    pos_as_dict: bool = True,
    log: Log = Log(),
    verbose: bool = True
):
    """
    Load, filter, and build numpy fasta records in a single pass for maximum performance.
    
    This function combines loading, filtering, and numpy array conversion in one pass,
    avoiding the creation of intermediate FastaRecord objects. This is significantly
    faster than loading and then building records separately.
    
    Parameters
    ----------
    path : str
        Path to FASTA file. Supports: .fa, .fasta, .fa.gz, .fasta.gz, .fa.bgz, .fasta.bgz
    chromlist_set : set
        Set of valid chromosome identifiers
    chroms_in_sumstats_set : set
        Set of chromosomes present in the summary statistics
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided, creates a default mapper.
    pos_as_dict : bool, default=True
        If True, return starting_positions and records_len as dictionaries
    log : gwaslab.g_Log.Log, default=Log()
        Logging object
    verbose : bool, default=True
        If True, print progress messages
        
    Returns
    -------
    tuple
        (record, starting_positions, records_len_dict) where:
        - record: concatenated numpy array of uint8 integers
        - starting_positions: dict or array of starting positions for each chromosome
        - records_len_dict: dict or array of lengths for each chromosome
    """
    log.write("   -Loading and building numpy fasta records:", end="", verbose=verbose)
    
    # Create mapper if not provided (for backward compatibility)
    if mapper is None:
        # Default mapper - will auto-detect format when needed
        mapper = ChromosomeMapper(
            log=log,
            verbose=verbose
        )
    
    all_r = []
    records_len = []
    chrom_keys = []  # Track chromosome keys in order
    
    pysam_success = False
    if PYSAM_AVAILABLE:
        try:
            with pysam.FastxFile(path) as fastx:
                for entry in fastx:
                    title = entry.name
                    sequence = entry.sequence
                    
                    # Get chromosome name from FASTA title (may be in various formats)
                    # Try to extract chromosome identifier from title
                    record_chr = title.strip()
                    
                    # Use mapper to convert to numeric format (sumstats format)
                    try:
                        # Detect format if not already detected
                        if mapper._sumstats_format is None:
                            mapper.detect_sumstats_format(pd.Series([record_chr]))
                        i = mapper.to_numeric(record_chr)
                        # Ensure it's numeric for comparison
                        if isinstance(i, str):
                            # Try to convert to numeric if it's a string number
                            try:
                                i = int(i)
                            except ValueError:
                                pass
                    except (KeyError, ValueError, AttributeError):
                        # Fallback: try stripping chr prefix and mapping
                        record_chr_stripped = record_chr.strip("chrCHR").upper()
                        try:
                            if mapper._sumstats_format is None:
                                mapper.detect_sumstats_format(pd.Series([record_chr_stripped]))
                            i = mapper.to_numeric(record_chr_stripped)
                            if isinstance(i, str):
                                try:
                                    i = int(i)
                                except ValueError:
                                    pass
                        except (KeyError, ValueError, AttributeError):
                            # Last resort: use as-is
                            i = record_chr_stripped
                    
                    # Only process if it passes filters
                    if (i in chromlist_set) and (i in chroms_in_sumstats_set):
                        log.write(record_chr, " ", end="", show_time=False, verbose=verbose)
                        # Convert directly to numpy array without creating FastaRecord
                        # Translate bytes using the translation table
                        r = sequence.encode('ascii').translate(_FASTA_TRANSLATE_TABLE)
                        r_len = len(r)
                        r = np.array([r], dtype=f'<U{r_len}').view('<u4').astype(np.uint8)
                        all_r.append(r)
                        records_len.append(r_len)
                        chrom_keys.append(i)
            pysam_success = True  # Successfully processed with pysam, skip fallback
        except Exception as e:
            # Fallback to old implementation if pysam fails
            import warnings
            warnings.warn(
                f"pysam FastxFile failed, falling back to slower implementation: {e}",
                UserWarning
            )
    
    # Use fallback only if pysam was not successful
    if not pysam_success:
        # Fallback to old implementation
        handle = _open_fasta_handle(path)
        try:
            # Skip any text before the first record
            title = None
            for line in handle:
                if line.startswith(">"):
                    title = line[1:].rstrip()
                    break
            else:
                # No records found - return empty arrays
                if pos_as_dict:
                    return np.array([], dtype=np.uint8), {}, {}
                else:
                    return np.array([], dtype=np.uint8), np.array([], dtype=np.int64), np.array([], dtype=np.int64)
            
            # Main parsing, filtering, and conversion loop
            lines = []
            for line in handle:
                if line.startswith(">"):
                    # Process the previous record before starting a new one
                    if title is not None:
                        sequence = "".join(lines).replace(" ", "").replace("\r", "")
                        # Get chromosome name from FASTA title
                        record_chr = title.strip()
                        
                        # Use mapper to convert to numeric format
                        try:
                            if mapper._sumstats_format is None:
                                mapper.detect_sumstats_format(pd.Series([record_chr]))
                            i = mapper.to_numeric(record_chr)
                            if isinstance(i, str):
                                try:
                                    i = int(i)
                                except ValueError:
                                    pass
                        except (KeyError, ValueError, AttributeError):
                            # Fallback: try stripping chr prefix
                            record_chr_stripped = record_chr.strip("chrCHR").upper()
                            try:
                                if mapper._sumstats_format is None:
                                    mapper.detect_sumstats_format(pd.Series([record_chr_stripped]))
                                i = mapper.to_numeric(record_chr_stripped)
                                if isinstance(i, str):
                                    try:
                                        i = int(i)
                                    except ValueError:
                                        pass
                            except (KeyError, ValueError, AttributeError):
                                i = record_chr_stripped
                        
                        # Only process if it passes filters
                        if (i in chromlist_set) and (i in chroms_in_sumstats_set):
                            log.write(record_chr, " ", end="", show_time=False, verbose=verbose)
                            # Convert directly to numpy array without creating FastaRecord
                            # Translate bytes using the translation table
                            r = sequence.encode('ascii').translate(_FASTA_TRANSLATE_TABLE)
                            r_len = len(r)
                            r = np.array([r], dtype=f'<U{r_len}').view('<u4').astype(np.uint8)
                            all_r.append(r)
                            records_len.append(r_len)
                            chrom_keys.append(i)
                    
                    # Start new record
                    lines = []
                    title = line[1:].rstrip()
                    continue
                lines.append(line.rstrip())
            
            # Process the last record
            if title is not None:
                sequence = "".join(lines).replace(" ", "").replace("\r", "")
                record_chr = title.strip()
                
                # Use mapper to convert to numeric format
                try:
                    if mapper._sumstats_format is None:
                        mapper.detect_sumstats_format(pd.Series([record_chr]))
                    i = mapper.to_numeric(record_chr)
                    if isinstance(i, str):
                        try:
                            i = int(i)
                        except ValueError:
                            pass
                except (KeyError, ValueError, AttributeError):
                    # Fallback: try stripping chr prefix
                    record_chr_stripped = record_chr.strip("chrCHR").upper()
                    try:
                        if mapper._sumstats_format is None:
                            mapper.detect_sumstats_format(pd.Series([record_chr_stripped]))
                        i = mapper.to_numeric(record_chr_stripped)
                        if isinstance(i, str):
                            try:
                                i = int(i)
                            except ValueError:
                                pass
                    except (KeyError, ValueError, AttributeError):
                        i = record_chr_stripped
                
                if (i in chromlist_set) and (i in chroms_in_sumstats_set):
                    log.write(record_chr, " ", end="", show_time=False, verbose=verbose)
                    # Convert directly to numpy array
                    r = sequence.encode('ascii').translate(_FASTA_TRANSLATE_TABLE)
                    r_len = len(r)
                    r = np.array([r], dtype=f'<U{r_len}').view('<u4').astype(np.uint8)
                    all_r.append(r)
                    records_len.append(r_len)
                    chrom_keys.append(i)
        
        finally:
            handle.close()
    
    log.write("", show_time=False, verbose=verbose)
    
    if len(all_r) == 0:
        if pos_as_dict:
            return np.array([], dtype=np.uint8), {}, {}
        else:
            return np.array([], dtype=np.uint8), np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    
    # Convert lengths to array and compute starting positions
    records_len = np.array(records_len, dtype=np.int64)
    starting_positions = np.cumsum(records_len) - records_len
    
    if pos_as_dict:
        # Use dict() constructor with zip for better performance
        starting_positions = dict(zip(chrom_keys, starting_positions))
        records_len_dict = dict(zip(chrom_keys, records_len))
    else:
        records_len_dict = records_len
    
    # Concatenate all arrays
    record = np.concatenate(all_r)
    del all_r  # free memory
    
    return record, starting_positions, records_len_dict


def build_fasta_records(fasta_records_dict, pos_as_dict=True, log=Log(), verbose=True):
    """
    Build numpy fasta records from a dictionary of FastaRecord objects.
    
    This function converts FASTA records to a single numpy array of integers for fast lookup.
    It uses a translation table to map nucleotides to integer codes.
    
    Parameters
    ----------
    fasta_records_dict : dict
        Dictionary mapping chromosome names to FastaRecord objects
    pos_as_dict : bool, default=True
        If True, return starting_positions and records_len as dictionaries
    log : gwaslab.g_Log.Log, default=Log()
        Logging object
    verbose : bool, default=True
        If True, print progress messages
        
    Returns
    -------
    tuple
        (record, starting_positions, records_len_dict) where:
        - record: concatenated numpy array of uint8 integers
        - starting_positions: dict or array of starting positions for each chromosome
        - records_len_dict: dict or array of lengths for each chromosome
    """
    log.write("   -Building numpy fasta records from dict", verbose=verbose)

    # Convert the fasta record to a numpy array of integers in a very fast way.
    # fasta_record.seq._data is a byte-string, so we can use the bytes.maketrans to apply a translation.
    # Here we map the bytes to the unicode character representing the desired integer as defined in the mapping dict
    # (i.e. b'A' -> '\x02', b'T' -> '\x03', b'C' -> '\x04', b'G' -> '\x05', b'N' -> '\x06')
    # Then, using np.array(... dtype=<U..) we convert the string to a numpy array of unicode characters.
    # Then, we do a magic with view('<u4') to convert the unicode characters to 4-byte integers, so we obtain the actual integer representation of the characters
    # Lastly, we cast the array to np.uint8 to convert the 4-byte integers to 1-byte integers to save memory
    # Full example:
    # fasta_record.seq._data = b'ACTGN' -> b'\x02\x04\x03\x05\x06' -> np.array(['\x02\x04\x03\x05\x06'], dtype='<U5') -> np.array([2, 4, 3, 5, 6], dtype=uint32) -> np.array([2, 4, 3, 5, 6], dtype=uint8)
    all_r = []
    records_len = []  # Build lengths list as we go to avoid second pass
    for r in fasta_records_dict.values():
        r = r.seq._data.translate(_FASTA_TRANSLATE_TABLE)
        r_len = len(r)
        r = np.array([r], dtype=f'<U{r_len}').view('<u4').astype(np.uint8)
        all_r.append(r)
        records_len.append(r_len)
    
    # We've just created a list of numpy arrays, so we can concatenate them to obtain a single numpy array
    # Then we keep track of the starting position of each record in the concatenated array. This will be useful later
    # to index the record array depending on the position of the variant and the chromosome
    records_len = np.array(records_len, dtype=np.int64)  # Convert to array once

    starting_positions = np.cumsum(records_len) - records_len

    
    if pos_as_dict:
        # Use dict() constructor with zip for better performance than dict comprehension
        keys = fasta_records_dict.keys()
        starting_positions = dict(zip(keys, starting_positions))
        records_len_dict = dict(zip(keys, records_len))
    else:
        records_len_dict = records_len
    record = np.concatenate(all_r)
    del all_r  # free memory

    return record, starting_positions, records_len_dict

