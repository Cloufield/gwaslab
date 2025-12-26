# GSF Format Specification

**GSF (GWASLab Standard Format)** is an efficient, columnar storage format specifically designed for GWAS summary statistics. **GSF is based on Apache Parquet format**, providing optimal compression, fast I/O performance, and support for predicate pushdown filtering.

## Overview

GSF format is GWASLab's native storage format based on **Apache Parquet**. It is optimized for:
- **Space efficiency**: Typically smaller than CSV.gz files (Parquet's columnar compression)
- **Fast I/O**: Parquet's columnar storage enables rapid loading and filtering
- **Type safety**: Parquet preserves data types and metadata
- **Filtering**: Supports efficient predicate pushdown during loading (Parquet feature)
- **Partitioning**: Optional partitioning by chromosome or other columns (Parquet partitioning)

## File Extension

GSF files use the `.gsf` extension. **Internally, GSF files are standard Apache Parquet files**. The `.gsf` extension indicates they follow GWASLab's specific conventions and optimizations, but they can be read by any Parquet-compatible tool by simply renaming the extension to `.parquet` or using Parquet readers directly.

## Format Specification

### Parquet-Based Format

**GSF is based on Apache Parquet format**. GSF files are standard Apache Parquet files with GWASLab-specific conventions and optimizations. The `.gsf` extension indicates these conventions, but the files can be read by any Parquet-compatible tool.

!!! important "Parquet Format"
    GSF files are **standard Parquet files**. They follow the Apache Parquet specification and can be read by any Parquet-compatible tool (PyArrow, pandas, R arrow, DuckDB, etc.) without modification. The `.gsf` extension is purely a convention to indicate GWASLab-specific optimizations.

Parquet provides:
- **Columnar storage** for efficient compression
- **Schema preservation** with data types
- **Metadata embedding** in file structure
- **Cross-platform compatibility** (Linux, macOS, Windows, cloud storage)
- **Predicate pushdown** for efficient filtering
- **Column pruning** for reading only needed columns

### Default Compression

GSF uses **Zstandard (zstd)** compression by default, which provides:
- Excellent compression ratios (comparable to gzip)
- Fast compression and decompression speeds
- Wide compatibility

Other supported compression codecs:
- `snappy` - Fast compression, moderate ratio
- `gzip` - High compression ratio, slower
- `brotli` - Very high compression ratio
- `lz4` - Very fast compression

## Column Organization

GSF files organize columns in an optimal order for GWAS data:

1. **Information columns**: SNPID, CHR, POS, EA, NEA, STATUS, rsID
2. **Statistics columns**: EAF, BETA, SE, P, OR, Z, etc.
3. **Other columns**: Any additional columns in the dataset

This ordering improves compression efficiency by grouping similar data types together.

## Data Type Optimization

GSF automatically optimizes data types for storage:

| GWASLab Type | GSF Storage Type | Notes |
|-------------|------------------|-------|
| `category` | `string` | Categories converted to strings for Parquet compatibility |
| `Int64` | `int64` | Preserved as integer |
| `float64` | `float64` | Preserved as float |
| `float32` | `float32` | Preserved as float |
| `string` | `string` | Preserved as string |

## Metadata Storage

GSF files embed comprehensive metadata in the Parquet file schema:

### Embedded Metadata Fields

- `gwaslab_format`: Format identifier ("gsf")
- `gwaslab_version`: Format version ("1.0")
- `gwaslab_study_name`: Study name (if available)
- `gwaslab_genome_build`: Genome build (e.g., "19", "38")
- `gwaslab_species`: Species (e.g., "human")
- `gwaslab_full_meta`: Complete metadata dictionary as JSON

The full metadata dictionary includes:
- GWASLab version information
- Study metadata
- Processing history
- Quality control information
- Any custom metadata added by the user

## Usage

### Saving to GSF Format

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("sumstats.txt.gz", ...)

# Save to GSF format
mysumstats.to_gsf("sumstats.gsf")

# With custom options
mysumstats.to_gsf(
    "sumstats.gsf",
    compression="zstd",  # or "snappy", "gzip", "brotli", "lz4"
    verbose=True
)
```

### Loading from GSF Format

```python
# Basic loading
mysumstats = gl.load_gsf("sumstats.gsf")

# With column selection
mysumstats = gl.load_gsf("sumstats.gsf", columns=["CHR", "POS", "P", "BETA"])

# With filtering (predicate pushdown)
mysumstats = gl.load_gsf("sumstats.gsf", filters="CHR == 7")
mysumstats = gl.load_gsf("sumstats.gsf", filters="P < 5e-8")
mysumstats = gl.load_gsf("sumstats.gsf", filters="CHR == 7 & P < 5e-8")
mysumstats = gl.load_gsf("sumstats.gsf", filters="CHR in [1, 2, 3]")
```

### Filter Syntax

GSF supports efficient filtering using predicate pushdown. Filter expressions support:

**Comparison operators:**
- `==` - Equality
- `!=` - Inequality
- `<` - Less than
- `<=` - Less than or equal
- `>` - Greater than
- `>=` - Greater than or equal

**Logical operators:**
- `&` - AND
- `|` - OR

**Special operators:**
- `in` - Membership test (e.g., `CHR in [1, 2, 3]`)

**Examples:**
```python
# Single condition
filters="CHR == 7"
filters="P < 5e-8"
filters="BETA > 0"

# Multiple conditions (AND)
filters="CHR == 7 & P < 5e-8"
filters="P < 5e-8 & BETA > 0"

# Multiple conditions (OR)
filters="CHR == 1 | CHR == 2"

# IN operator
filters="CHR in [1, 2, 3, 7, 22]"

# Complex expressions
filters="CHR == 7 & P < 5e-8 & BETA > 0"
```

**Note:** Filtering uses predicate pushdown, meaning filters are applied at the storage level before loading data into memory. This provides significant performance benefits for large datasets.

## Partitioned GSF Files

GSF supports partitioning for very large datasets. Partitioned files are stored as directories containing multiple Parquet files.

### Creating Partitioned Files

```python
# Partition by chromosome
mysumstats.to_gsf("sumstats_partitioned", partition_cols=["CHR"])

# This creates a directory structure:
# sumstats_partitioned/
#   CHR=1/
#     part-0.parquet
#   CHR=2/
#     part-0.parquet
#   ...
```

### Loading Partitioned Files

```python
# Load partitioned file (same API)
mysumstats = gl.load_gsf("sumstats_partitioned")

# Filtering on partition columns is especially efficient
mysumstats = gl.load_gsf("sumstats_partitioned", filters="CHR == 7")
```

**Benefits of partitioning:**
- Faster queries on partition columns
- Can process individual partitions independently
- Better compression (data grouped by partition value)

### Reading with Other Tools

**GSF files are standard Parquet files**, so they can be read by any Parquet-compatible tool without modification:

**Python:**
```python
import pyarrow.parquet as pq
table = pq.read_table("sumstats.gsf")
df = table.to_pandas()
```

## Examples

### Complete Workflow

```python
import gwaslab as gl

# Load data
mysumstats = gl.Sumstats("input.txt.gz", ...)

# Process and QC
mysumstats.basic_check()
mysumstats.harmonize(...)

# Save to GSF
mysumstats.to_gsf("processed_sumstats.gsf", compression="zstd")

# Later: Load with filtering
significant_variants = gl.load_gsf(
    "processed_sumstats.gsf",
    filters="P < 5e-8 & CHR in [1, 2, 3, 7, 22]",
    columns=["CHR", "POS", "EA", "NEA", "BETA", "SE", "P"]
)

# Work with filtered data
print(f"Found {len(significant_variants.data)} significant variants")
```

### Partitioned Workflow

```python
# Save partitioned by chromosome
mysumstats.to_gsf(
    "sumstats_partitioned",
    partition_cols=["CHR"],
    compression="zstd"
)

# Load specific chromosome
chr7_data = gl.load_gsf("sumstats_partitioned", filters="CHR == 7")

# Process chromosome by chromosome
for chr_num in range(1, 23):
    chr_data = gl.load_gsf("sumstats_partitioned", filters=f"CHR == {chr_num}")
    # Process chromosome data
    ...
```

## Technical Details

### Internal Structure

**GSF files are standard Apache Parquet files** with GWASLab-specific optimizations:
- **Parquet row groups**: Data organized into row groups for efficient access
- **Parquet column chunks**: Each column stored separately (Parquet's columnar format)
- **Parquet metadata**: Embedded in Parquet file metadata structure
- **Parquet schema**: Preserves data types and column information using Parquet's type system
- **GWASLab metadata**: Custom metadata stored in Parquet's key-value metadata fields

### Type Casting for Filters

GSF automatically handles type casting when filtering:
- String columns compared to numeric values are automatically cast
- Ensures filters work correctly regardless of storage type
- Maintains performance through efficient casting

### Metadata Preservation

All GWASLab metadata is preserved:
- Study information
- Processing history
- Quality control flags
- Custom user metadata
- Genome build information

## Limitations

1. **Schema evolution**: Adding/removing columns requires rewriting the file

## Migration Guide

### From CSV.gz to GSF

```python
# Load CSV
mysumstats = gl.Sumstats("sumstats.csv.gz", ...)

# Save as GSF
mysumstats.to_gsf("sumstats.gsf")

# Verify
loaded = gl.load_gsf("sumstats.gsf")
assert len(loaded.data) == len(mysumstats.data)
```

### From Pickle to GSF

```python
# Load Pickle
mysumstats = gl.load_pickle("sumstats.pickle")

# Save as GSF
mysumstats.to_gsf("sumstats.gsf")
```

## Version History

- **v1.0** (2024): Initial GSF format specification
  - Parquet-based storage
  - Zstd compression default
  - Metadata embedding
  - Filtering support
  - Partitioning support

## References

- [Apache Parquet Format](https://parquet.apache.org/)
- [PyArrow Documentation](https://arrow.apache.org/docs/python/)
- [Zstandard Compression](https://facebook.github.io/zstd/)
- [GWASLab Documentation](https://cloufield.github.io/gwaslab/)

## See Also

- [Format Load Save Documentation](format_load_save.md) - General I/O operations
- [Format Documentation](Format.md) - Other supported formats
- [Sumstats Object Documentation](SumstatsObject.md) - GWASLab Sumstats class

