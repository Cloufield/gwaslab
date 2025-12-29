# FASTA I/O in GWASLab

GWASLab provides efficient functions for reading and writing FASTA format files, with support for compressed formats and optimized memory usage.

!!! note "Chromosome Conversion for High-Level Functions"
    For high-level functions like `check_ref()`, chromosome format conversion is handled automatically by `ChromosomeMapper`. No manual `chr_dict` parameter is needed. The `chr_dict` parameters shown in this documentation are for low-level utility functions that may still use them for backward compatibility.

## Supported File Formats

GWASLab supports the following FASTA file extensions:
- `.fa`, `.fasta` - Plain text FASTA files
- `.fa.gz`, `.fasta.gz` - Gzip-compressed FASTA files
- `.fa.bgz`, `.fasta.bgz` - BGZF-compressed FASTA files (treated as gzip)

## Reading FASTA Files

### `parse_fasta()`

Parse a FASTA file and return records as a dictionary or iterator.

```python
from gwaslab.io.io_fasta import parse_fasta

# Load as dictionary (default)
records = parse_fasta("reference.fasta")
print(records["chr1"][:10])  # First 10 bases of chr1

# Load as iterator (memory efficient for large files)
for title, seq in parse_fasta("reference.fasta", as_dict=False):
    print(f"{title}: {len(seq)} bp")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `path` | `string` | Path to FASTA file. Supports plain and gzipped formats | Required |
| `as_dict` | `boolean` | If True, return a dictionary mapping titles to sequences. If False, return an iterator of (title, sequence) tuples | `True` |

**Returns:**
- If `as_dict=True`: `dict[str, str]` - Dictionary mapping record titles to sequences
- If `as_dict=False`: `Iterator[tuple[str, str]]` - Iterator yielding (title, sequence) tuples

!!! note "Memory Efficiency"
    When `as_dict=False`, the file handle is kept open until iteration completes. This is more memory-efficient for large files.

---

### `load_fasta_auto()`

Automatically load FASTA or gzipped FASTA files, returning FastaRecord objects or tuples.

```python
from gwaslab.io.io_fasta import load_fasta_auto

# Get FastaRecord objects (default)
for record in load_fasta_auto("reference.fasta.gz"):
    print(f"{record.id}: {len(record.seq._data)} bp")

# Get simple tuples
for title, seq in load_fasta_auto("reference.fasta.gz", as_seqrecord=False):
    print(f"{title}: {len(seq)} bp")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `path` | `string` | Path to FASTA file. Supports: .fa, .fasta, .fa.gz, .fasta.gz, .fa.bgz, .fasta.bgz | Required |
| `as_seqrecord` | `boolean` | If True, return FastaRecord objects with .id and .seq._data attributes. If False, return (title, sequence) tuples | `True` |

**Returns:**
- `Iterator[FastaRecord]` or `Iterator[tuple[str, str]]` - Iterator yielding FastaRecord objects or (title, sequence) tuples

!!! warning "File Handle Management"
    The file handle is kept open until iteration completes. For proper resource management, ensure you consume the entire iterator or use it within a context manager.

---

### `load_fasta_filtered()`

Load and filter FASTA records in a single pass for better performance. Only creates records for chromosomes that are needed.

```python
from gwaslab.io.io_fasta import load_fasta_filtered
from gwaslab.info.g_Log import Log

chromlist_set = {"1", "2", "3"}
chroms_in_sumstats_set = {"1", "2"}
chr_dict = {"1": "1", "2": "2", "3": "3"}
chr_dict_keys = set(chr_dict.keys())
log = Log()

records = load_fasta_filtered(
    path="reference.fasta.gz",
    chromlist_set=chromlist_set,
    chroms_in_sumstats_set=chroms_in_sumstats_set,
    chr_dict=chr_dict,
    chr_dict_keys=chr_dict_keys,
    log=log,
    verbose=True
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `path` | `string` | Path to FASTA file | Required |
| `chromlist_set` | `set` | Set of valid chromosome identifiers | Required |
| `chroms_in_sumstats_set` | `set` | Set of chromosomes present in the summary statistics | Required |
| `chr_dict` | `dict` | Dictionary mapping chromosome names to standardized format | Required |
| `chr_dict_keys` | `set` | Set of keys in chr_dict for fast lookup | Required |
| `log` | `gwaslab.g_Log.Log` | Logging object | `Log()` |
| `verbose` | `boolean` | If True, print progress messages | `True` |

**Returns:**
- `dict[str, FastaRecord]` - Dictionary mapping chromosome identifiers to FastaRecord objects

---

### `get_fasta_record()`

Get a specific FASTA record by title.

```python
from gwaslab.io.io_fasta import get_fasta_record

seq = get_fasta_record("reference.fasta", "chr1")
if seq:
    print(f"chr1 length: {len(seq)}")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `path` | `string` | Path to FASTA file | Required |
| `title` | `string` | Title of the record to retrieve (without the '>' character) | Required |

**Returns:**
- `str` or `None` - Sequence for the requested record, or None if not found

---

## Writing FASTA Files

### `write_fasta()`

Write FASTA records to a file.

```python
from gwaslab.io.io_fasta import write_fasta

# Write from dictionary
records = {"chr1": "ATCGATCG", "chr2": "GCTAGCTA"}
write_fasta(records, "output.fasta")

# Write from iterator
records = [("chr1", "ATCGATCG"), ("chr2", "GCTAGCTA")]
write_fasta(records, "output.fasta.gz", wrap=80)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `records` | `dict[str, str]` or `Iterator[tuple[str, str]]` | FASTA records to write. Can be a dictionary mapping titles to sequences or an iterator of (title, sequence) tuples | Required |
| `path` | `string` | Output file path. Supports plain and gzipped formats based on extension | Required |
| `wrap` | `int` | Line length for sequence wrapping. Use 0 or None for no wrapping | `60` |
| `mode` | `string` | File open mode. Use "w" for write or "a" for append | `"w"` |

**Supported output formats:**
- Plain text: `.fa`, `.fasta`
- Gzip-compressed: `.fa.gz`, `.fasta.gz`
- BGZF: `.fa.bgz`, `.fasta.bgz` (output will be standard gzip format)

---

## Advanced Functions

### `load_and_build_fasta_records()`

Load, filter, and build numpy FASTA records in a single pass for maximum performance. This function combines loading, filtering, and numpy array conversion in one pass, avoiding the creation of intermediate FastaRecord objects.

```python
from gwaslab.io.io_fasta import load_and_build_fasta_records
from gwaslab.info.g_Log import Log

chromlist_set = {"1", "2"}
chroms_in_sumstats_set = {"1", "2"}
chr_dict = {"1": "1", "2": "2"}
chr_dict_keys = set(chr_dict.keys())
log = Log()

record, starting_positions, records_len_dict = load_and_build_fasta_records(
    path="reference.fasta.gz",
    chromlist_set=chromlist_set,
    chroms_in_sumstats_set=chroms_in_sumstats_set,
    chr_dict=chr_dict,
    chr_dict_keys=chr_dict_keys,
    pos_as_dict=True,
    log=log,
    verbose=True
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `path` | `string` | Path to FASTA file | Required |
| `chromlist_set` | `set` | Set of valid chromosome identifiers | Required |
| `chroms_in_sumstats_set` | `set` | Set of chromosomes present in the summary statistics | Required |
| `chr_dict` | `dict` | Dictionary mapping chromosome names to standardized format | Required |
| `chr_dict_keys` | `set` | Set of keys in chr_dict for fast lookup | Required |
| `pos_as_dict` | `boolean` | If True, return starting_positions and records_len as dictionaries | `True` |
| `log` | `gwaslab.g_Log.Log` | Logging object | `Log()` |
| `verbose` | `boolean` | If True, print progress messages | `True` |

**Returns:**
- `tuple` - (record, starting_positions, records_len_dict) where:
  - `record`: concatenated numpy array of uint8 integers
  - `starting_positions`: dict or array of starting positions for each chromosome
  - `records_len_dict`: dict or array of lengths for each chromosome

---

### `build_fasta_records()`

Build numpy FASTA records from a dictionary of FastaRecord objects. This function converts FASTA records to a single numpy array of integers for fast lookup.

```python
from gwaslab.io.io_fasta import build_fasta_records, load_fasta_filtered
from gwaslab.info.g_Log import Log

# First load records
fasta_records_dict = load_fasta_filtered(...)

# Then build numpy records
record, starting_positions, records_len_dict = build_fasta_records(
    fasta_records_dict,
    pos_as_dict=True,
    log=Log(),
    verbose=True
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `fasta_records_dict` | `dict` | Dictionary mapping chromosome names to FastaRecord objects | Required |
| `pos_as_dict` | `boolean` | If True, return starting_positions and records_len as dictionaries | `True` |
| `log` | `gwaslab.g_Log.Log` | Logging object | `Log()` |
| `verbose` | `boolean` | If True, print progress messages | `True` |

**Returns:**
- `tuple` - (record, starting_positions, records_len_dict) where:
  - `record`: concatenated numpy array of uint8 integers
  - `starting_positions`: dict or array of starting positions for each chromosome
  - `records_len_dict`: dict or array of lengths for each chromosome

---

## FastaRecord Class

The `FastaRecord` class provides compatibility with code that expects FASTA record objects with `.id` and `.seq._data` attributes.

```python
from gwaslab.io.io_fasta import FastaRecord

record = FastaRecord(id="chr1", seq="ATCGATCG")
print(record.id)  # "chr1"
print(len(record.seq._data))  # 8
```

**Attributes:**
- `id` (str): Record identifier (title line without '>')
- `seq._data` (bytes): Sequence data as bytes

---

## Examples

!!! example "Basic FASTA reading"
    ```python
    from gwaslab.io.io_fasta import parse_fasta
    
    # Load entire file as dictionary
    records = parse_fasta("reference.fasta.gz")
    chr1_seq = records["chr1"]
    print(f"chr1 length: {len(chr1_seq)} bp")
    ```

!!! example "Memory-efficient iteration"
    ```python
    from gwaslab.io.io_fasta import parse_fasta
    
    # Process large file without loading everything into memory
    total_length = 0
    for title, seq in parse_fasta("large_reference.fasta.gz", as_dict=False):
        total_length += len(seq)
        print(f"Processed {title}: {len(seq)} bp")
    print(f"Total length: {total_length} bp")
    ```

!!! example "Writing FASTA files"
    ```python
    from gwaslab.io.io_fasta import write_fasta
    
    # Write from dictionary
    sequences = {
        "chr1": "ATCGATCGATCG",
        "chr2": "GCTAGCTAGCTA"
    }
    write_fasta(sequences, "output.fasta", wrap=60)
    
    # Write compressed file
    write_fasta(sequences, "output.fasta.gz", wrap=80)
    ```

!!! example "Filtered loading for harmonization"
    ```python
    from gwaslab.io.io_fasta import load_fasta_filtered
    from gwaslab.info.g_Log import Log
    
    # Only load chromosomes present in your sumstats
    chroms_in_sumstats = {"1", "2", "3"}
    valid_chroms = {"1", "2", "3", "4", "5"}
    chr_dict = {str(i): str(i) for i in range(1, 23)}
    chr_dict_keys = set(chr_dict.keys())
    
    records = load_fasta_filtered(
        path="reference.fasta.gz",
        chromlist_set=valid_chroms,
        chroms_in_sumstats_set=chroms_in_sumstats,
        chr_dict=chr_dict,
        chr_dict_keys=chr_dict_keys,
        log=Log(),
        verbose=True
    )
    ```

!!! example "High-performance numpy array conversion"
    ```python
    from gwaslab.io.io_fasta import load_and_build_fasta_records
    from gwaslab.info.g_Log import Log
    
    # Load and convert to numpy arrays in one pass
    record, positions, lengths = load_and_build_fasta_records(
        path="reference.fasta.gz",
        chromlist_set={"1", "2"},
        chroms_in_sumstats_set={"1", "2"},
        chr_dict={"1": "1", "2": "2"},
        chr_dict_keys={"1", "2"},
        pos_as_dict=True,
        log=Log(),
        verbose=True
    )
    
    # Access sequence data
    chr1_start = positions["1"]
    chr1_length = lengths["1"]
    chr1_seq = record[chr1_start:chr1_start + chr1_length]
    ```

---

## Performance Tips

1. **Use iterators for large files**: When processing large FASTA files, use `as_dict=False` or `load_fasta_auto()` to avoid loading everything into memory.

2. **Filter during loading**: Use `load_fasta_filtered()` to only load chromosomes you need, avoiding unnecessary memory usage.

3. **Use numpy arrays for fast lookup**: For repeated sequence lookups, use `load_and_build_fasta_records()` or `build_fasta_records()` to convert to numpy arrays.

4. **Compressed files**: GWASLab automatically handles gzip and BGZF compression, so you can work directly with compressed files without manual decompression.

