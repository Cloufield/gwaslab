# Command Line Interface (CLI)

GWASLab provides a unified command-line interface for processing GWAS summary statistics. The CLI supports quality control (QC), harmonization, format conversion, and various output formatting options.

## Basic Usage

The CLI follows a unified interface pattern:

```bash
gwaslab --input <file> --fmt <format> [--options] --to-fmt <format> --out <file>
```

### Quick Examples

```bash
# Show version
gwaslab version

# Basic QC and output
gwaslab --input sumstats.tsv --fmt auto --qc --out cleaned.tsv --to-fmt gwaslab

# Harmonization with reference
gwaslab --input sumstats.tsv --fmt auto --ref-seq ref.fasta --harmonize --out harmonized.tsv --to-fmt gwaslab

# Format conversion only
gwaslab --input sumstats.tsv --fmt gwaslab --out sumstats.ldsc --to-fmt ldsc
```

## Command Structure

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--input` | Input sumstats file path | `--input data/sumstats.tsv` |
| `--fmt` | Input format (default: `auto`) | `--fmt gwaslab` or `--fmt auto` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--out` | Output file path (prefix) | None |
| `--to-fmt` | Output format | `gwaslab` |
| `--tab-fmt` | Tabular format (`tsv`, `csv`, `parquet`) | `tsv` |
| `--nrows` | Number of rows to read (for testing) | None |
| `--quiet` | Suppress output messages | False |
| `--threads` | Number of threads for parallel processing | 1 |

## Processing Options

### Quality Control (QC)

Perform quality control on sumstats using `basic_check()`:

```bash
# Basic QC
gwaslab --input sumstats.tsv --fmt auto --qc --out cleaned.tsv --to-fmt gwaslab

# QC with remove bad variants
gwaslab --input sumstats.tsv --fmt auto --qc --remove --out cleaned.tsv --to-fmt gwaslab

# QC with remove duplicates
gwaslab --input sumstats.tsv --fmt auto --qc --remove-dup --out cleaned.tsv --to-fmt gwaslab

# QC with normalize indels
gwaslab --input sumstats.tsv --fmt auto --qc --normalize --out cleaned.tsv --to-fmt gwaslab

# QC with all options
gwaslab --input sumstats.tsv --fmt auto --qc --remove --remove-dup --normalize --threads 4 --out cleaned.tsv --to-fmt gwaslab
```

**QC Options:**

| Option | Description |
|--------|-------------|
| `--qc` | Perform quality control (basic_check) |
| `--remove` | Remove bad quality variants detected during QC |
| `--remove-dup` | Remove duplicated or multi-allelic variants |
| `--normalize` | Normalize indels (e.g., ATA:AA -> AT:A) |

### Harmonization

Harmonize sumstats with reference data:

```bash
# Basic harmonization (without reference files)
gwaslab --input sumstats.tsv --fmt auto --harmonize --out harmonized.tsv --to-fmt gwaslab

# Harmonization with reference sequence for allele flipping
gwaslab --input sumstats.tsv --fmt auto --harmonize --ref-seq /path/to/reference.fasta --out harmonized.tsv --to-fmt gwaslab

# Harmonization with rsID assignment
gwaslab --input sumstats.tsv --fmt auto --harmonize --ref-rsid-vcf /path/to/reference.vcf.gz --out harmonized.tsv --to-fmt gwaslab

# Full harmonization pipeline
gwaslab --input sumstats.tsv --fmt auto \
  --harmonize \
  --ref-seq /path/to/reference.fasta \
  --ref-rsid-vcf /path/to/rsid.vcf.gz \
  --ref-infer /path/to/inference.vcf.gz \
  --ref-alt-freq AF \
  --maf-threshold 0.40 \
  --ref-maf-threshold 0.5 \
  --ref-seq-mode v \
  --sweep-mode \
  --threads 8 \
  --out harmonized.tsv \
  --to-fmt gwaslab
```

**Harmonization Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--harmonize` | Perform harmonization | False |
| `--basic-check` | Run basic QC in harmonization | True |
| `--no-basic-check` | Skip basic QC in harmonization | - |
| `--ref-seq` | Reference sequence file (FASTA) for allele flipping | None |
| `--ref-rsid-tsv` | Reference rsID TSV file | None |
| `--ref-rsid-vcf` | Reference rsID VCF/BCF file | None |
| `--ref-infer` | Reference VCF/BCF file for strand inference | None |
| `--ref-alt-freq` | Allele frequency field name in VCF INFO | `AF` |
| `--ref-maf-threshold` | MAF threshold for reference | 0.5 |
| `--maf-threshold` | MAF threshold for sumstats | 0.40 |
| `--ref-seq-mode` | Reference sequence mode (`v`=vectorized, `s`=sequential) | `v` |
| `--sweep-mode` | Use sweep mode for large datasets | False |

### Assign rsID

Assign rsID to variants using reference data:

```bash
# Assign rsID from TSV file
gwaslab --input sumstats.tsv --fmt auto --assign-rsid --ref-rsid-tsv /path/to/rsid.tsv --out output.tsv --to-fmt gwaslab

# Assign rsID from VCF file
gwaslab --input sumstats.tsv --fmt auto --assign-rsid --ref-rsid-vcf /path/to/rsid.vcf.gz --overwrite empty --out output.tsv --to-fmt gwaslab
```

**Assign rsID Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--assign-rsid` | Assign rsID to variants | False |
| `--overwrite` | Overwrite mode (`all`, `invalid`, `empty`) | `empty` |

### rsID to CHR:POS

Convert rsID to CHR:POS coordinates:

```bash
# Convert rsID to CHR:POS
gwaslab --input sumstats.tsv --fmt auto --rsid-to-chrpos --ref-tsv /path/to/chrpos.tsv --build 19 --out output.tsv --to-fmt gwaslab
```

**rsID to CHR:POS Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--rsid-to-chrpos` | Convert rsID to CHR:POS | False |
| `--build` | Genome build version | `19` |
| `--overwrite-rtc` | Overwrite existing CHR:POS | False |
| `--chunksize` | Chunk size for processing | 5000000 |

## Output Formatting Options

### Basic Formatting

```bash
# Output in gwaslab format (default)
gwaslab --input sumstats.tsv --fmt auto --out output --to-fmt gwaslab

# Output in LDSC format
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt ldsc

# Output in PLINK format
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt plink

# Output as CSV instead of TSV
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --tab-fmt csv

# Output without compression
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --no-gzip
```

### Advanced Formatting

```bash
# Output with bgzip compression and tabix index
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --bgzip --tabix

# Extract only HapMap3 variants
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --hapmap3

# Exclude HLA region
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --exclude-hla

# Exclude HLA with custom range (in Mbp)
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --exclude-hla --hla-lower 20 --hla-upper 30

# Add chromosome prefix
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --chr-prefix chr

# Use numeric notation for X, Y, MT (23, 24, 25)
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --xymt-number

# Add N column with specified value
gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab --n 10000
```

**Output Formatting Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--to-fmt` | Output format | `gwaslab` |
| `--tab-fmt` | Tabular format (`tsv`, `csv`, `parquet`) | `tsv` |
| `--no-gzip` | Disable gzip compression | False (gzip enabled) |
| `--bgzip` | Use bgzip compression | False |
| `--tabix` | Create tabix index (requires bgzip) | False |
| `--hapmap3` | Extract HapMap3 variants only | False |
| `--exclude-hla` | Exclude HLA region | False |
| `--hla-lower` | HLA region lower bound (Mbp) | 25 |
| `--hla-upper` | HLA region upper bound (Mbp) | 34 |
| `--n` | Add N column with specified value | None |
| `--chr-prefix` | Prefix for chromosome column | `""` |
| `--xymt-number` | Use numeric notation for X, Y, MT | False |

## Complete Workflow Examples

### Example 1: Basic QC Pipeline

```bash
# Load, QC, and save sumstats
gwaslab --input raw_sumstats.tsv \
  --fmt auto \
  --qc \
  --remove \
  --remove-dup \
  --normalize \
  --threads 4 \
  --out cleaned_sumstats \
  --to-fmt gwaslab \
  --quiet
```

This will:
1. Load sumstats from `raw_sumstats.tsv` with auto-detection
2. Perform QC (fix IDs, chromosomes, positions, alleles)
3. Remove bad quality variants
4. Remove duplicates
5. Normalize indels
6. Save to `cleaned_sumstats.gwaslab.tsv.gz`

### Example 2: Harmonization Pipeline

```bash
# Full harmonization with all reference files
gwaslab --input sumstats.tsv \
  --fmt auto \
  --harmonize \
  --ref-seq /data/reference/hg19.fasta.gz \
  --ref-rsid-vcf /data/reference/dbsnp.vcf.gz \
  --ref-infer /data/reference/1kg.vcf.gz \
  --ref-alt-freq AF \
  --maf-threshold 0.40 \
  --ref-maf-threshold 0.5 \
  --ref-seq-mode v \
  --sweep-mode \
  --threads 8 \
  --out harmonized_sumstats \
  --to-fmt gwaslab \
  --tab-fmt tsv
```

This will:
1. Load sumstats
2. Run basic QC
3. Check alignment with reference sequence
4. Flip alleles if needed
5. Assign rsIDs from VCF
6. Infer strand for palindromic SNPs
7. Save harmonized sumstats

### Example 3: Format Conversion for LDSC

```bash
# Convert to LDSC format for heritability analysis
gwaslab --input sumstats.tsv \
  --fmt gwaslab \
  --out sumstats_ldsc \
  --to-fmt ldsc \
  --no-gzip
```

### Example 4: Prepare Data for Meta-analysis

```bash
# QC, harmonize, and format for meta-analysis
gwaslab --input study1.tsv \
  --fmt auto \
  --qc \
  --remove-dup \
  --harmonize \
  --ref-seq /data/reference/hg19.fasta.gz \
  --ref-rsid-vcf /data/reference/dbsnp.vcf.gz \
  --out study1_ready \
  --to-fmt gwaslab \
  --exclude-hla \
  --n 50000
```

### Example 5: Extract HapMap3 Variants

```bash
# Extract only HapMap3 variants for replication
gwaslab --input discovery_sumstats.tsv \
  --fmt gwaslab \
  --out hapmap3_sumstats \
  --to-fmt gwaslab \
  --hapmap3 \
  --build 19
```

### Example 6: Multi-step Processing

```bash
# Step 1: QC
gwaslab --input raw.tsv --fmt auto --qc --remove-dup --out step1_qc --to-fmt gwaslab

# Step 2: Harmonize
gwaslab --input step1_qc.gwaslab.tsv.gz --fmt gwaslab --harmonize \
  --ref-seq ref.fasta --out step2_harmonized --to-fmt gwaslab

# Step 3: Format for specific tool
gwaslab --input step2_harmonized.gwaslab.tsv.gz --fmt gwaslab \
  --out final --to-fmt plink --tab-fmt tsv
```

## Output File Naming

The CLI follows a consistent naming pattern for output files:

- **Basic format**: `{output_path}.{format}.{tab_fmt}.gz`
  - Example: `output.gwaslab.tsv.gz`

- **With filters**: `{output_path}.{filter}.{format}.{tab_fmt}.gz`
  - Example: `output.hapmap3.gwaslab.tsv.gz` (with `--hapmap3`)
  - Example: `output.noMHC.gwaslab.tsv.gz` (with `--exclude-hla`)

- **Without compression**: `{output_path}.{format}.{tab_fmt}`
  - Example: `output.gwaslab.tsv` (with `--no-gzip`)

- **Log file**: `{output_path}.{format}.log`
  - Example: `output.gwaslab.log`

## Tips and Best Practices

### 1. Input Format Detection

Use `--fmt auto` to let GWASLab automatically detect the input format:

```bash
gwaslab --input sumstats.tsv --fmt auto --qc --out output --to-fmt gwaslab
```

### 2. Parallel Processing

Use `--threads` to speed up processing for large files:

```bash
gwaslab --input large_sumstats.tsv --fmt auto --qc --threads 8 --out output --to-fmt gwaslab
```

### 3. Quiet Mode

Use `--quiet` to suppress verbose output in scripts:

```bash
gwaslab --input sumstats.tsv --fmt auto --qc --out output --to-fmt gwaslab --quiet
```

### 4. Testing with Small Samples

Use `--nrows` to test commands on a subset of data:

```bash
gwaslab --input large_sumstats.tsv --fmt auto --nrows 1000 --qc --out test_output --to-fmt gwaslab
```

### 5. Combining Operations

You can combine multiple processing steps in a single command:

```bash
# QC + Harmonization + Format conversion
gwaslab --input sumstats.tsv --fmt auto \
  --qc --remove-dup \
  --harmonize --ref-seq ref.fasta \
  --out output --to-fmt gwaslab
```

### 6. Reference Files

For harmonization, you typically need:
- **Reference sequence (FASTA)**: For allele flipping
- **rsID reference (VCF/TSV)**: For rsID assignment
- **Inference reference (VCF)**: For strand inference of palindromic SNPs

Download reference files from:
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
- [1000 Genomes](https://www.internationalgenome.org/)
- [UCSC Genome Browser](https://genome.ucsc.edu/)

### 7. Memory Considerations

For very large files:
- Use `--sweep-mode` for harmonization (faster for large datasets)
- Process in chunks if memory is limited
- Consider using `--nrows` for testing first

## Common Use Cases

### Use Case 1: Quick Format Check

```bash
# Just load and check the file (no output)
gwaslab --input sumstats.tsv --fmt auto
```

### Use Case 2: Standard QC Workflow

```bash
gwaslab --input sumstats.tsv --fmt auto \
  --qc --remove --remove-dup --normalize \
  --out qc_sumstats --to-fmt gwaslab
```

### Use Case 3: Prepare for LDSC

```bash
gwaslab --input sumstats.tsv --fmt gwaslab \
  --out ldsc_input --to-fmt ldsc --no-gzip
```

### Use Case 4: Prepare for Meta-analysis

```bash
gwaslab --input sumstats.tsv --fmt auto \
  --qc --remove-dup \
  --harmonize --ref-seq ref.fasta --ref-rsid-vcf dbsnp.vcf.gz \
  --out meta_ready --to-fmt gwaslab --exclude-hla
```

### Use Case 5: Extract Replication Set

```bash
gwaslab --input discovery.tsv --fmt gwaslab \
  --out replication --to-fmt gwaslab --hapmap3 --build 19
```

## Troubleshooting

### Issue: File not found

**Error**: `FileNotFoundError: [Errno 2] No such file or directory`

**Solution**: Check that the input file path is correct and accessible.

### Issue: Format detection fails

**Error**: Format auto-detection doesn't work

**Solution**: Specify the format explicitly with `--fmt`:
```bash
gwaslab --input sumstats.tsv --fmt gwaslab --qc --out output --to-fmt gwaslab
```

### Issue: Memory errors with large files

**Error**: Out of memory errors

**Solution**: 
- Use `--threads 1` to reduce memory usage
- Process in smaller chunks
- Use `--sweep-mode` for harmonization

### Issue: Reference file errors

**Error**: Reference file not found or invalid

**Solution**: 
- Verify reference file paths are correct
- Check that reference files are properly formatted
- Ensure VCF files are indexed if using VCF references

## Getting Help

For more information:

```bash
# Show help message
gwaslab --help

# Show version
gwaslab version
```

For detailed documentation on specific functions, see:
- [QC & Filtering](QC&Filtering.md)
- [Harmonization](Harmonization.md)
- [Format](Format.md)
- [Standardization](Standardization.md)

## Supported Formats

GWASLab supports many input and output formats through the [formatbook](https://github.com/Cloufield/formatbook) repository. Common formats include:

- **Input**: `auto`, `gwaslab`, `plink`, `ldsc`, `vcf`, and many more
- **Output**: `gwaslab`, `ldsc`, `plink`, `plink2`, `saige`, `fastgwa`, `regenie`, `vcf`, and many more

Check the formatbook repository for the complete list of supported formats.

