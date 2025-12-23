# Logging Strategy in GWASLab

## Overview

GWASLab uses a custom logging system based on a `Log` class (`gwaslab.g_Log.Log`) rather than Python's standard `logging` module. This system provides both console output and persistent log storage.

## Architecture

### Log Class (`src/gwaslab/g_Log.py`)

The `Log` class maintains:
- **`log_text`**: A string buffer that accumulates all log messages with timestamps
- **Core Methods**:
  - `write()`: Primary logging method (also aliased as `log()`)
  - `warning()`: For warning messages (prefixed with "#WARNING!")
  - `show()`: Display the entire log
  - `save(path)`: Write log to file
  - `get_log_for_last_operation()`: Extract logs for the most recent operation
  - `combine()`: Merge logs from multiple sources

- **Standardized Logging Methods** (for unified formatting):
  - All methods support `indent` parameter (default 0) for hierarchical logging
  - `log_variants_filtered(count, reason=None, indent=0)`: Log filtered variants count
  - `log_variants_removed(count, reason=None, indent=0)`: Log removed variants count
  - `log_variants_kept(count, reason=None, indent=0)`: Log kept variants count
  - `log_column_added(column_name, indent=0)`: Log addition of a column
  - `log_column_dropped(column_name, reason=None, indent=0)`: Log removal of a column
  - `log_shape_change(before_shape, after_shape, indent=0)`: Log DataFrame shape changes
  - `log_dataframe_shape(sumstats, indent=0)`: Log current DataFrame shape and memory usage
  - `log_operation_start(operation_name, version=None, indent=0)`: Log operation start
  - `log_operation_finish(operation_name, indent=0)`: Log operation completion
  - `log_filtering_condition(column, operator, threshold, count, action="Removing", indent=0)`: Log filtering conditions
  - `log_operation(message, prefix=" -", indent=0)`: Log general operation messages
  - `log_status_change(digit, before, after, count=None, reason=None, indent=0)`: Log STATUS column digit changes
  - `log_datatype_change(column, from_dtype, to_dtype, success=True, indent=0)`: Log datatype conversions
  - `log_datatype_attempt(column, from_dtype, to_dtype, indent=0)`: Log datatype conversion attempts
  - `log_formula(target_column, formula, source_columns=None, indent=0)`: Log formulas/calculations used
  - `log_reference_path(ref_type, path, indent=0)`: Log reference file paths (VCF, FASTA, TSV, etc.)
  - `log_threads(n_cores, indent=0)`: Log number of threads/cores being used

### Key Characteristics

1. **Dual Output**: Every log message is:
   - Conditionally printed to stdout (controlled by `verbose` parameter)
   - Always stored in `log_text` string buffer

2. **Timestamps**: All messages include timestamps in format `YYYY/MM/DD HH:MM:SS`

3. **No Log Levels**: Unlike standard logging, there's no concept of DEBUG/INFO/WARNING/ERROR levels - only a binary `verbose` flag

## Usage Patterns

### Initialization

Log objects are created per `Sumstats` instance:
```
self.log = Log()  # Created in Sumstats.__init__()
```

### Verbose Parameter

- **Default**: `verbose=True` (most functions default to verbose output)
- **Control**: Passed through function calls to control console output
- **Storage**: Messages are stored regardless of `verbose` setting

### Logging Decorators

Two main decorators add standardized logging:

1. **`with_logging()`** (`qc/qc_decorator.py`):
   - Logs operation start/finish
   - Logs reference files (VCF, FASTA, TSV)
   - Logs thread/core counts
   - Validates required columns and arguments
   - Optionally checks DataFrame shapes and dtypes

2. **`with_logging_filter()`** (`util/util_in_filter_value.py`):
   - Similar to `with_logging()` but specialized for filtering operations
   - Logs variant counts before/after filtering

### Typical Log Message Flow

```
Start to [operation] ...(version)
 -Reference VCF: [path]
 -Number of threads/cores to use: [n]
 -[DataFrame shape info]
Finished [operation].
```

## Current Issues (Verbosity Concerns)

### 1. **Default Verbosity**
- Most functions default to `verbose=True`
- Users must explicitly set `verbose=False` to reduce output
- No global verbosity control

### 2. **No Log Levels**
- All messages treated equally (no distinction between critical errors and routine info)
- Cannot filter by importance/severity
- Warnings use string prefix but no structural separation

### 3. **Extensive Logging**
- Every operation logs start/finish messages
- Shape checks, column validations, and dtype checks all produce log messages
- Reference file paths, thread counts logged for every operation
- Estimated 3668+ log calls across 197 files

### 4. **Always-On Storage**
- Even with `verbose=False`, all messages are stored in memory
- Can lead to large `log_text` strings for long-running operations

### 5. **No Selective Filtering**
- Cannot enable verbose logging for specific operations while keeping others quiet
- Binary on/off control only

## Statistics

- **Log calls**: 3668+ matches across 197 files
- **Default verbose=True**: 997+ instances
- **Main decorators**: 2 (`with_logging`, `with_logging_filter`)

## Standardized Logging Methods

To unify logging format across the codebase, the `Log` class provides standardized methods for common operations:

### Variant Operations
```
# Filtering variants
log.log_variants_filtered(150, reason="with P > 0.05")
# Output: " -Filtered out 150 variants with P > 0.05"

# Removing variants
log.log_variants_removed(25, reason="with nan in CHR column")
# Output: " -Removed 25 variants with nan in CHR column"

# Keeping variants
log.log_variants_kept(1000, reason="in specified regions")
# Output: " -Keeping 1000 variants in specified regions"
```

### Column Operations
```
# Adding columns
log.log_column_added("MLOG10P")
# Output: " -Added column: MLOG10P"

# Dropping columns
log.log_column_dropped("INFO", reason="all values are invalid")
# Output: " -Dropped column INFO (all values are invalid)"
```

### Shape and Data Operations
```
# Shape changes
log.log_shape_change((1000, 20), (850, 21))
# Output: " -Shape changed: 1000 x 20 -> 850 x 21 (-150 rows, +1 columns)"

# Current shape
log.log_dataframe_shape(sumstats)
# Output: " -Current Dataframe shape : 1000 x 20 ; Memory usage: 19.95 MB"

# Operation lifecycle
log.log_operation_start("filter variants by condition", version="v3.4.38")
log.log_operation_finish("filtering variants")
```

### Filtering Conditions
```
# Filtering conditions
log.log_filtering_condition("P", ">", 0.05, 150, action="Removing")
# Output: " -Removing 150 variants with P > 0.05 ..."
```

### General Operations
```
# General operation logging
log.log_operation("Sorting variants by chromosome and position")
# Output: " -Sorting variants by chromosome and position"

log.log_operation("Processing batch 1 of 10", prefix="   ")
# Output: "   Processing batch 1 of 10"
```

### Indentation for Nested Operations

All standardized logging methods support an `indent` parameter to show hierarchical structure:

```
# Main operation (indent=0, default)
log.log_operation_start("process data", version="v3.4.38")
# Output: "Start to process data ...(v3.4.38)"

# Inner operation (indent=1, 2 spaces)
log.log_operation_start("filter variants", indent=1)
# Output: "  Start to filter variants ..."

# Nested inner operation (indent=2, 4 spaces)
log.log_variants_filtered(150, reason="with P > 0.05", indent=2)
# Output: "    -Filtered out 150 variants with P > 0.05"

log.log_variants_removed(25, reason="duplicate variants", indent=2)
# Output: "    -Removed 25 variants duplicate variants"

# Back to inner operation level
log.log_operation_finish("filtering variants", indent=1)
# Output: "  Finished filtering variants."

# Back to main operation level
log.log_operation_finish("processing data")
# Output: "Finished processing data."
```

**Example with nested function calls:**
```
def process_data(log, verbose=True):
    log.log_operation_start("process data", verbose=verbose)
    log.log_dataframe_shape(sumstats, verbose=verbose)
    
    # Call inner function
    filter_variants(sumstats, log, verbose=verbose, indent=1)
    
    log.log_operation_finish("processing data", verbose=verbose)

def filter_variants(sumstats, log, verbose=True, indent=0):
    log.log_operation_start("filter variants", indent=indent, verbose=verbose)
    
    # Nested operations
    log.log_variants_filtered(150, reason="with P > 0.05", indent=indent+1, verbose=verbose)
    log.log_variants_removed(25, reason="duplicates", indent=indent+1, verbose=verbose)
    
    log.log_operation_finish("filtering variants", indent=indent, verbose=verbose)
```

**Output:**
```
Start to process data ...
 -Current Dataframe shape : 1000 x 20 ; Memory usage: 19.95 MB
  Start to filter variants ...
    -Filtered out 150 variants with P > 0.05
    -Removed 25 variants duplicates
  Finished filtering variants.
Finished processing data.
```

### Status Changes
```
# STATUS column digit changes
log.log_status_change(digit=1, before="9", after="1", count=5000, reason="genome build update")
# Output: " -Updated STATUS digit 1: 9 -> 1 (5000 variants) - genome build update"

log.log_status_change(digit=6, before=["4", "5"], after=["1", "2"])
# Output: " -Updated STATUS digit 6: 4/5 -> 1/2"
```

### Datatype Changes
```
# Datatype conversions
log.log_datatype_change("CHR", "object", "Int64", success=True)
# Output: " -Converted datatype for CHR: object -> Int64"

log.log_datatype_change("POS", "float64", "Int64", success=False)
# Output: " -Failed to convert datatype for POS: float64 -> Int64"

# For attempts where result is logged separately
log.log_datatype_attempt("CHR", "object", "Int64")
# Output: " -Trying to convert datatype for CHR: object -> Int64..." (no newline)
```

### Formulas and Calculations
```
# Formula/calculation logging
log.log_formula("MLOG10P", "from P")
# Output: "    Filling MLOG10P from P..."

log.log_formula("Z", "from BETA/SE", source_columns=["BETA", "SE"])
# Output: "    Filling Z from BETA/SE (using BETA, SE)..."

log.log_formula("OR", "from BETA")
# Output: "    Filling OR from BETA..."
```

### Reference Paths and Resources
```
# Reference file paths
log.log_reference_path("VCF", "/path/to/reference.vcf.gz")
# Output: " -Reference VCF: /path/to/reference.vcf.gz"

log.log_reference_path("FASTA", "/path/to/genome.fa")
# Output: " -Reference FASTA: /path/to/genome.fa"

log.log_reference_path("TSV", "/path/to/annotation.tsv")
# Output: " -Reference TSV: /path/to/annotation.tsv"

# Thread/core usage
log.log_threads(8)
# Output: " -Number of threads/cores to use: 8"
```

**Benefits of Standardized Methods:**
- Consistent formatting across the codebase
- Easier to parse and analyze logs programmatically
- Reduced code duplication
- Easier to update log format globally
- Better readability and maintainability

## Potential Improvements

1. **Add Log Levels**: Implement DEBUG/INFO/WARNING/ERROR levels
2. **Global Verbosity Control**: Centralized verbosity setting
3. **Selective Verbosity**: Per-module or per-operation verbosity control
4. **Lazy Evaluation**: Only store messages when needed (e.g., when saving log)
5. **Structured Logging**: Use structured format (JSON) for better parsing
6. **Progress Indicators**: Replace frequent messages with progress bars for long operations
7. **Log Rotation**: Limit memory usage for very long operations
8. **Migrate to Standardized Methods**: Gradually replace ad-hoc `log.write()` calls with standardized methods

