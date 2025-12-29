# Chromosome Handling in GWASLab

GWASLab provides a flexible system for handling chromosome identifiers across different formats, species, and reference files using a **numeric middle layer architecture**.

## Architecture

The chromosome mapping system uses three layers:

1. **Middle Layer (Numeric)**: Species-specific numeric identifiers (e.g., human: 1-22 autosomes, 23=X, 24=Y, 25=MT)
2. **Sumstats Layer**: Maps sumstats chromosome notation ↔ numbers ↔ sumstats notation
3. **Reference Layer**: Maps reference file chromosome notation ↔ numbers ↔ reference notation

**Conversion Flow:**
```
Sumstats Format → Numeric → Reference Format
     (e.g., "chr1")  →  (1)  →  ("NC_000001.11")
```

This design enables seamless matching between sumstats and reference files, even when they use different notations.

## Supported Formats

GWASLab automatically detects and handles four chromosome formats:

| Format | Examples | Use Case |
|--------|----------|----------|
| **Numeric** | `1`, `2`, `23`, `24`, `25` | Common in GWAS files (23=X, 24=Y, 25=MT) |
| **String** | `"1"`, `"2"`, `"X"`, `"Y"`, `"MT"` | Text-based formats (case-insensitive) |
| **Chr-prefixed** | `"chr1"`, `"chr2"`, `"chrX"` | UCSC-style, many VCF files (case variations supported) |
| **NCBI RefSeq** | `"NC_000001.11"`, `"NC_000023.11"` | Reference genome files (build-specific) |

## Species Support

GWASLab supports 12+ species with species-specific chromosome structures:

| Species | Autosomes | Sex Chromosomes | Examples |
|---------|-----------|----------------|----------|
| Human | 1-22 | X, Y | Default species |
| Mouse | 1-19 | X, Y | `species="mus musculus"` |
| Chicken | 1-28 | Z, W | `species="gallus gallus"` |
| Zebrafish | 1-25 | - | `species="danio rerio"` |

Other supported species: Rat, Fruit fly, Pig, Cattle, Dog, Horse, Rice, Arabidopsis.

## Build Detection

**For Reference Files:**
- Build is **automatically detected** from NCBI RefSeq notation (NC_*)
- Works with VCF, FASTA, GTF, and Chain files

**For Sumstats Data:**
- Build parameter **must be explicitly specified**: `ChromosomeMapper(build="38")`
- Supported builds: "19" (hg19/GRCh37), "38" (hg38/GRCh38), and species-specific builds

## Usage Examples

### Basic Usage

```python
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
import pandas as pd

# Initialize mapper (defaults to human)
mapper = ChromosomeMapper(species="homo sapiens")

# Detect sumstats format
mapper.detect_sumstats_format(pd.Series([1, 2, 23, "X", "MT"]))
# Returns: 'numeric'

# Convert sumstats to reference format
mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
# Returns: 'chr1' (if VCF uses chr-prefixed format)
```

### Format Conversions

```python
mapper = ChromosomeMapper()

# Numeric format
mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
mapper.sumstats_to_number(1)  # Returns: 1

# String format
mapper.detect_sumstats_format(pd.Series(["1", "2", "X"]))
mapper.sumstats_to_number("X")  # Returns: 23

# Chr-prefixed format
mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))
mapper.sumstats_to_number("chr1")  # Returns: 1

# NCBI RefSeq (requires build for sumstats)
mapper = ChromosomeMapper(build="38")
mapper.detect_sumstats_format(pd.Series(["NC_000001.11", "NC_000023.11"]))
mapper.sumstats_to_number("NC_000001.11")  # Returns: 1

# NCBI RefSeq in reference files (auto-detects build)
mapper = ChromosomeMapper()  # build not specified
mapper.detect_reference_format("reference.vcf.gz")  # Contains NC_000001.11
# Build "38" is automatically detected
```

### Matching Sumstats with Reference Files

```python
mapper = ChromosomeMapper()
mapper.detect_sumstats_format(pd.Series([1, 2, 23]))

# Convert sumstats (numeric) to reference format (chr-prefixed)
mapper.sumstats_to_reference(1, reference_file="reference.vcf.gz")
# Returns: 'chr1'

# Reverse conversion
mapper.reference_to_sumstats("chr1", reference_file="reference.vcf.gz")
# Returns: 1
```

### Convenience Methods

```python
mapper = ChromosomeMapper()
mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))

mapper.to_numeric("chr1")  # Returns: 1
mapper.to_string(1)  # Returns: '1'
mapper.to_chr(1)  # Returns: 'chr1'

# With build specified
mapper = ChromosomeMapper(build="38")
mapper.to_nc(1)  # Returns: 'NC_000001.11'
```

## Integration with Sumstats Object

The `Sumstats` object automatically creates and manages a `ChromosomeMapper`:

```python
from gwaslab.g_Sumstats import Sumstats

mysumstats = Sumstats(data, verbose=False)
# Mapper auto-detects format from data
# Automatically handles chromosome conversion with reference files
```

## Reference File Support

GWASLab automatically detects chromosome format from:
- **VCF files**: `.vcf`, `.vcf.gz`, `.bcf`
- **FASTA files**: `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`, `.fa.bgz`, `.fasta.bgz`
- **GTF/GFF files**: `.gtf`, `.gff`, `.gtf.gz`, `.gff.gz`
- **Chain files**: `.chain`, `.chain.gz`

## Handling Unconvertible Chromosomes

GWASLab gracefully handles chromosomes that cannot be converted:
- **Alternative loci**: `"1_KI270766v1_alt"`, `"KI270766.1"`
- **Unplaced sequences**: `"GL000195.1"`, `"GL000009.2"`
- **Chromosome arms**: `"1p"`, `"1q"`, `"chr1p"`
- **Malformed identifiers**: Invalid or unexpected formats

**Behavior**: Returns `pd.NA` instead of raising an error, allowing processing to continue.

```python
mapper = ChromosomeMapper()
mapper.sumstats_to_number("1_KI270766v1_alt")  # Returns: <NA>
mapper.sumstats_to_number("1p")  # Returns: <NA>
```

## Best Practices

1. **Always detect format**: Call `detect_sumstats_format()` for correct handling
2. **Specify build for NC notation**: Use `ChromosomeMapper(build="38")` when using NCBI RefSeq IDs in sumstats
3. **Let reference files auto-detect**: Reference files automatically detect format and build
4. **Handle unconvertible chromosomes**: Filter or handle `pd.NA` values as needed
5. **Use Sumstats object**: Automatically manages chromosome mapping

## Summary

GWASLab's chromosome handling provides:
- ✅ Flexible format support (numeric, string, chr-prefixed, NCBI RefSeq)
- ✅ Multi-species support (12+ species)
- ✅ Automatic format detection from sumstats and reference files
- ✅ Build auto-detection for reference files with NCBI RefSeq notation
- ✅ Graceful error handling (returns `pd.NA` for unconvertible chromosomes)
- ✅ Seamless integration with `Sumstats` object
