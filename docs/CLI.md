# Command Line Interface (CLI)

GWASLab provides a unified command-line interface for processing GWAS summary statistics. The CLI supports quality control (QC), harmonization, format conversion, and various output formatting options.

> **Warning — heavy development:** The CLI is changing quickly (flags, defaults, and behavior may shift between releases). Pin a GWASLab version for reproducible pipelines, check `gwaslab --help` after upgrading, and expect occasional breaking changes until the interface stabilizes.

## Basic Usage

The CLI follows a unified interface pattern:

```bash
gwaslab --input <file> --fmt <format> [--options] --to-fmt <format> --out <file>
```

### Quick Examples

```bash
# Show version
gwaslab version

# Show default and current config paths
gwaslab config

# Show one configured path by key
gwaslab config show config
gwaslab config show reference

# Resolve built-in path key
gwaslab path config

# List all formats in formatbook
gwaslab formatbook list

# Show one format mapping
gwaslab formatbook show metal

# Reference catalog + downloads (see "Download" and "list ref" sections)
gwaslab list ref --available
gwaslab download ref 1kg_eas_hg19
# Sumstats download: output dir is -o / --output-dir / -d / --directory (all equivalent)
gwaslab download sumstats GCST90270926 --directory downloads
gwaslab download-sumstats GCST90270926 --output-dir downloads

# Basic QC and output
gwaslab --input sumstats.tsv --fmt auto --qc --out cleaned.tsv --to-fmt gwaslab

# Harmonization with reference
gwaslab --input sumstats.tsv --fmt auto --ref-seq ref.fasta --harmonize --out harmonized.tsv --to-fmt gwaslab

# Format conversion only
gwaslab --input sumstats.tsv --fmt gwaslab --out sumstats.ldsc --to-fmt ldsc

# Plot (CLI runs fix_chr + fix_pos if basic_check was not run)
gwaslab --input sumstats.tsv --plot manhattan --out manhattan.png

# Assign rsID (CLI runs fix_chr + fix_pos if basic_check was not run)
gwaslab --input sumstats.tsv --fmt auto --assign-rsid --ref-rsid-vcf /path/to/rsid.vcf.gz --out output.tsv --to-fmt gwaslab

# Liftover (CLI runs fix_chr + fix_pos if basic_check was not run)
gwaslab --input sumstats.tsv --liftover 19 38 --out lifted_hg38.tsv

# Infer build (hg19/hg38) from HapMap3 coordinates
gwaslab --input sumstats.tsv --infer-build --out inferred.tsv

# Extract lead signals
gwaslab --input sumstats.tsv --get lead --out lead.tsv
```

## Utility Subcommands

### config

Inspect GWASLab path configuration and query a single configured path.

```bash
# Show default + current path config
gwaslab config

# Show one configured path by keyword
gwaslab config show config
gwaslab config show reference

# Persist a path override (written to ~/.gwaslab/settings.json)
gwaslab config set data_directory /data/refs
gwaslab config set config /data/refs/config.json

# JSON output (easy to parse)
gwaslab config --json
```

**Environment variables** (override defaults without editing settings):

| Variable | Purpose |
|----------|---------|
| `GWASLAB_DATA_DIR` | Default download / data root |
| `GWASLAB_CONFIG` | Path to local registry `config.json` |

**Options:**

| Option | Description |
|--------|-------------|
| `--json` | Print as JSON |
| `show <keyword>` | Show JSON content for `config`/`reference`/`formatbook`; otherwise print resolved path |
| `set <key> <path>` | Persist `data_directory` or `config` (`key` must be one of those two) |

### path

Resolve a local path by built-in key or downloaded reference keyword.

```bash
# Built-in keys
gwaslab path config
gwaslab path reference
gwaslab path formatbook
gwaslab path data_directory

# Downloaded reference keyword
gwaslab path <downloaded_keyword>
```

**Options:**

| Option | Description |
|--------|-------------|
| `keyword` | Built-in key or downloaded reference keyword |

### formatbook

Inspect and update format definitions in the formatbook (e.g., `saige`, `metal`).

```bash
# List all available formats
gwaslab formatbook list

# Show mapping for one format
gwaslab formatbook show saige

# JSON output
gwaslab formatbook list --json

# Update local formatbook from remote repository
gwaslab formatbook update
```

**Actions and options:**

| Command | Description |
|--------|-------------|
| `formatbook list` | List available formats in formatbook |
| `formatbook show <format>` | Show header mapping for one format |
| `formatbook update` | Update formatbook from remote source |
| `--json` | Print output in JSON format (`list` only) |

### Download (references and GWAS Catalog sumstats)

GWASLab keeps **flat legacy commands** (`download-ref`, `download-sumstats`) and adds a **grouped form** so documentation and flags line up:

| Goal | Grouped command | Legacy (same behavior) |
|------|-------------------|-------------------------|
| Fetch a packaged reference by keyword | `gwaslab download ref KEY` | `gwaslab download-ref KEY` |
| Fetch GWAS Catalog sumstats by GCST | `gwaslab download sumstats GCST…` | `gwaslab download-sumstats GCST…` |

**Unified output directory flag for sumstats:** any of `-o`, `--output-dir`, `-d`, or `--directory` sets the download folder (same underlying option).

```bash
# References
gwaslab download ref 1kg_eas_hg19
gwaslab download ref 1kg_eas_hg19 --directory ~/.gwaslab --overwrite
gwaslab download-ref 1kg_eas_hg19 --directory ~/.gwaslab

# GWAS Catalog sumstats
gwaslab download sumstats GCST90270926
gwaslab download sumstats GCST90270926 --directory ./downloads
gwaslab download-sumstats GCST90270926 -o ./downloads
```

Reference downloads also accept `--local-filename` and `--overwrite` (see `gwaslab download ref --help`).

### list ref

List **available** reference keywords (from the bundled/updated reference catalog) and/or **downloaded** entries registered in your local config. With no scope flags, both sections are shown.

```bash
gwaslab list ref
gwaslab list ref --available
gwaslab list ref --downloaded
gwaslab list ref --downloaded --source catalog
gwaslab list ref --downloaded --kind ref
gwaslab list ref --downloaded --kind sumstats
gwaslab list ref --available --downloaded --json
```

| Option | Description |
|--------|-------------|
| `--available` | Only keywords you can install with `gwaslab download ref …` |
| `--downloaded` | Only keywords already recorded under `downloaded` in config |
| `--source` | Filter `--downloaded` by registry `source` (`catalog`, `gwas_catalog`, or `local`) |
| `--kind` | Filter `--downloaded` by registry `kind` (`ref` or `sumstats`) |
| `--json` | Machine-readable output |
| `-q` / `--quiet` | Less library logging |

GWAS Catalog sumstats do not have a browseable `list` command in the CLI (you supply a `GCST…` ID).

### ref

Manage local registry entries without downloading from the catalog.

```bash
gwaslab ref add 1kg_eas_hg19 /data/refs/EAS.vcf.gz --tbi /data/refs/EAS.vcf.gz.tbi
gwaslab ref remove 1kg_eas_hg19
gwaslab ref remove 1kg_eas_hg19 --delete-file
```

| Subcommand | Description |
|------------|-------------|
| `add KEYWORD PATH` | Register a local file (`gl.add_local_data`) |
| `remove KEYWORD` | Drop registry entry; add `--delete-file` to remove blobs |

### init

Prepare the local reference registry: create `~/.gwaslab/`, migrate legacy package-local `config.json` if needed, optionally set `data_directory`, then scan for files to register.

```bash
gwaslab init
gwaslab init --directory /data/refs
gwaslab init -d /data/refs --recursive --quiet
```

| Option | Description |
|--------|-------------|
| `-d` / `--directory` | Set and persist `data_directory`, then scan this folder |
| `-r` / `--recursive` | Scan subdirectories when matching files to catalog keywords |
| `-q` / `--quiet` | Less library logging |

## Command Structure

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--input` | Input sumstats file path (required for main processing mode) | `--input data/sumstats.tsv` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--out` | Output file path (prefix) | None |
| `--to-fmt` | Output format | `gwaslab` |
| `--tab-fmt` | Tabular format (`tsv`, `csv`, `parquet`) | `tsv` |
| `--nrows` | Number of rows to read (for testing) | None |
| `--quiet` | Suppress output messages | False |
| `--threads` | Number of threads for parallel processing | 1 |

Many other flags (`--fix-chr`, `--fix-chr-pos`, variant filters, `--get`, `--plot-chr`, harmonization, liftover, etc.) are summarized in **Processing Options** below and in `gwaslab --help`.

### Plot / Get shared optional arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--sig-level` | Significance threshold used in Manhattan/MQQ/Regional plotting | `5e-8` |
| `--ylim <min> <max>` | Y-axis limits for plotting | None |
| `--highlight <id ...>` | Variant IDs (e.g. rsID/SNPID) to highlight in Manhattan/MQQ/Regional plots | None |
| `--sig-level-extract` | P-value threshold for `--get` operations | `5e-8` |
| `--windowsizekb` | Window size (kb) for `--get lead` | `500` |

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

### Statistics fill (`--fill`)

Derive missing statistical columns using [`fill_data()`](Conversion.md) — the same conversions available in the Python API. Runs **after QC** and **before** region/variant filters.

```bash
# MLOG10P → P before LDSC export
gwaslab --input sumstats.tsv --fill P --out out.ldsc --to-fmt ldsc

# EAF → MAF before MAF filtering
gwaslab --input sumstats.tsv --qc --fill MAF --maf 0.01 --out filtered.tsv --to-fmt gwaslab

# Extreme MLOG10P from Z/BETA/SE (P < 1e-300)
gwaslab --input sumstats.tsv --fill MLOG10P --fill-extreme --out out.tsv --to-fmt gwaslab

# Overwrite existing P values
gwaslab --input sumstats.tsv --fill P --fill-overwrite --out out.tsv
```

**Fill options:**

| Option | Description |
|--------|-------------|
| `--fill COL [COL ...]` | Target column(s) to derive. Valid: `OR`, `OR_95L`, `OR_95U`, `BETA`, `SE`, `P`, `Z`, `CHISQ`, `MLOG10P`, `MAF`, `SIG` (case-insensitive) |
| `--fill-overwrite` | Overwrite existing values in target columns (default: skip columns that already exist) |
| `--fill-extreme` | Use log-space MLOG10P calculation for extreme P values |
| `--fill-only-sig` | Only fill significant variants (uses `--sig-level`) |
| `--fill-df COL` | Column with degrees of freedom when filling `CHISQ` |

Note: `SIG` creates a `SIGNIFICANT` boolean column. See [Statistics conversion](Conversion.md) for conversion priority and formulas.

### Optional coordinate and ID fixes

Run individual `Sumstats.fix_*()` steps without full `--qc`. These run **before** QC when combined in one command. Use them to normalize columns before export or downstream steps.

| Option | Description |
|--------|-------------|
| `--fix-chr` | `fix_chr()` only (chromosome notation) |
| `--fix-pos` | `fix_pos()` only (position dtype / range) |
| `--fix-chr-pos` | `fix_chr()` then `fix_pos()` (same as `--fix-chr --fix-pos`) |
| `--fix-chr-pos-allele` | `fix_chr()`, `fix_pos()`, and `fix_allele()` |
| `--fix-allele` | `fix_allele()` only (allele notation) |
| `--fix-id` | `fix_id()` only (SNPID / rsID column) |

If both `fix_chr` and `fix_pos` run (via `--fix-chr-pos` or the pair `--fix-chr --fix-pos`), the CLI treats coordinates as ready for steps that normally auto-run `fix_chr` + `fix_pos` (e.g. liftover, `--assign-rsid`).

### Variant filters (PLINK-style)

Subset rows **after** `--filter-region` (if any) and **before** harmonization, assign-rsid, liftover, plotting, or `--get`. List-based filters require a **`SNPID`** or **`rsID`** column. BED and chromosome filters run `fix_chr` + `fix_pos` first if coordinates are not already normalized.

See also `examples/10_cli/09_variant_filters.sh` for runnable demos (uses `../../src` when run from a git checkout).

```bash
# Keep only IDs in a file (one per line; # comments and extra columns ignored)
gwaslab --input sumstats.tsv --extract keep.txt --out subset.tsv --to-fmt gwaslab

# Drop IDs from a file
gwaslab --input sumstats.tsv --exclude drop.txt --out pruned.tsv --to-fmt gwaslab

# Autosomes only (example)
gwaslab --input sumstats.tsv --chr 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 --out autosomes.tsv --to-fmt gwaslab

# MAF / MAC (uses MAF, or EAF/FRQ as min(f, 1−f); MAC from column or 2×N×MAF if N present)
gwaslab --input sumstats.tsv --maf 0.01 --max-maf 0.5 --mac 20 --out filtered.tsv --to-fmt gwaslab

# SNPs only; imputation quality
gwaslab --input sumstats.tsv --snps-only --min-info 0.8 --out qc.tsv --to-fmt gwaslab
```

**Variant filter options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--extract` | Path to file of variant IDs to **keep** (first column per line) | None |
| `--exclude` | Path to file of variant IDs to **remove** | None |
| `--extract-bed` | BED path: keep variants overlapping intervals (0-based half-open; uses `--build`) | None |
| `--exclude-bed` | BED path: remove variants overlapping intervals | None |
| `--chr` | Keep only listed chromosomes (space-separated; PLINK 2 `--chr` analog). For **`--plot regional`**, a **single** `--chr` together with `--start` and `--end` is interpreted as the plot chromosome instead of a filter | None |
| `--maf` | Minimum MAF | None |
| `--max-maf` | Maximum MAF | None |
| `--mac` | Minimum MAC (uses `MAC` column, or `2 × N × MAF` when `N` and a frequency column exist) | None |
| `--snps-only` | Keep rows where `EA` and `NEA` are single-nucleotide | False |
| `--min-info` | Minimum **INFO** (requires a column named `INFO`, case-insensitive) | None |

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
  --maf-threshold 0.40 \
  --ref-maf-threshold 0.4 \
  --sweep-mode \
  --threads 8 \
  --out harmonized.tsv \
  --to-fmt gwaslab
```

**Harmonization Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--harmonize` | Perform harmonization | False |
| `--basic-check` | Force basic QC in harmonization | `True` when explicitly set |
| `--no-basic-check` | Skip basic QC in harmonization | - |
| `--ref-seq` | Reference sequence file (FASTA) for allele flipping | None |
| `--ref-rsid-tsv` | Reference rsID HDF5 file (legacy name, accepts HDF5 path) | None |
| `--ref-rsid-vcf` | Reference rsID VCF/BCF file | None |
| `--ref-infer` | Reference VCF/BCF file for strand inference | None |
| `--ref-alt-freq` | INFO field name for ALT allele frequency when using `--ref-infer` | `AF` |
| `--ref-maf-threshold` | MAF threshold for reference | 0.4 |
| `--maf-threshold` | MAF threshold for sumstats | 0.40 |
| `--sweep-mode` | Use sweep mode for large datasets | False |

Default behavior when neither `--basic-check` nor `--no-basic-check` is provided: harmonization uses `basic_check = not --qc` (i.e., it runs basic check by default unless `--qc` was already run in the same command).

### Infer Build

Infer genome build (`hg19`/`hg38`) from HapMap3 SNP coordinates:

```bash
gwaslab --input sumstats.tsv --infer-build --out inferred.tsv --to-fmt gwaslab
```

**Infer Build Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--infer-build` | Infer genome build from HapMap3 coordinates | False |

### Assign rsID

Assign rsID to variants using reference data:

```bash
# Assign rsID from HDF5 file
gwaslab --input sumstats.tsv --fmt auto --assign-rsid --ref-rsid-tsv /path/to/rsid.hdf5 --out output.tsv --to-fmt gwaslab

# Assign rsID from VCF file
gwaslab --input sumstats.tsv --fmt auto --assign-rsid --ref-rsid-vcf /path/to/rsid.vcf.gz --overwrite empty --out output.tsv --to-fmt gwaslab
```

**Assign rsID Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--assign-rsid` | Assign rsID to variants (auto runs `fix_chr` + `fix_pos` if basic_check not run) | False |
| `--ref-rsid-tsv` | Reference rsID HDF5 file (legacy name, accepts HDF5 path) | None |
| `--ref-rsid-vcf` | Reference rsID VCF/BCF file | None |
| `--overwrite` | Overwrite mode (`all`, `invalid`, `empty`) | `empty` |
| `--threads` | Number of threads for parallel processing | 1 |

### rsID to CHR:POS

Convert rsID to CHR:POS coordinates:

```bash
# Convert rsID to CHR:POS using VCF (auto-generates HDF5)
gwaslab --input sumstats.tsv --fmt auto --rsid-to-chrpos --ref-rsid-vcf /path/to/reference.vcf.gz --build 19 --out output.tsv --to-fmt gwaslab

# Convert rsID to CHR:POS using existing HDF5 file
gwaslab --input sumstats.tsv --fmt auto --rsid-to-chrpos --ref-rsid-tsv /path/to/reference.hdf5 --build 19 --out output.tsv --to-fmt gwaslab
```

**rsID to CHR:POS Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--rsid-to-chrpos` | Convert rsID to CHR:POS | False |
| `--ref-rsid-vcf` | Reference VCF file for rsID to CHR:POS conversion (auto-generates HDF5) | None |
| `--ref-rsid-tsv` | Reference HDF5 file path for rsID to CHR:POS conversion | None |
| `--build` | Genome build version | `19` |
| `--threads` | Number of threads for parallel processing | 4 (when using rsid-to-chrpos) |

### Liftover

Convert coordinates between genome builds:

```bash
# Liftover from hg19 to hg38
gwaslab --input sumstats.tsv --fmt auto --liftover 19 38 --out lifted_hg38.tsv --to-fmt gwaslab
```

**Liftover Notes:**

- CLI auto-runs `fix_chr` + `fix_pos` before liftover if `basic_check` was not run.
- You can still run `--qc` earlier in the same command when full QC is preferred.

### Plotting

Generate plots from one input sumstats file:

```bash
# Manhattan
gwaslab --input sumstats.tsv --plot manhattan --out manhattan.png

# QQ
gwaslab --input sumstats.tsv --plot qq --out qq.png

# Combined Manhattan+QQ
gwaslab --input sumstats.tsv --plot mqq --out mqq.png

# Regional (explicit chromosome flag)
gwaslab --input sumstats.tsv --plot regional --plot-chr 6 --start 26000000 --end 34000000 --out region.png

# Regional (legacy: one --chr with --start/--end is treated as the plot region, not a chromosome filter)
gwaslab --input sumstats.tsv --plot regional --chr 6 --start 26000000 --end 34000000 --out region.png
```

Forest plots are not supported on the CLI; use the Python API (`gl.plot_forest()`, see [Forest plot](ForestPlot.md)).

**Plot Options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--plot` | Plot type: `manhattan`, `qq`, `mqq`, `regional`, `miami` | None |
| `--sig-level` | Significance threshold for Manhattan/MQQ/Regional plots | `5e-8` |
| `--ylim` | Y-axis range for Manhattan/MQQ/Regional plots (`min max`) | None |
| `--highlight` | Variant IDs to highlight in Manhattan/MQQ/Regional plots | None |
| `--plot-chr` | Chromosome for `--plot regional` (use with `--start`, `--end`) | None |
| `--start`, `--end` | Genomic interval for `--plot regional` (1-based positions as in CLI) | None |
| `--chr` | With **`--plot regional`**: if exactly **one** value is given **and** `--start`/`--end` are set, that value is the regional chromosome (same role as `--plot-chr`). With **multiple** values, or without regional plotting, `--chr` is a **variant filter** (see [Variant filters](#variant-filters-plink-style)) | None |

Notes:
- `--plot miami` is currently not available from single-input CLI mode and exits with guidance to use Python API.
- As with liftover/assign-rsid, CLI runs `fix_chr` + `fix_pos` before plotting if `basic_check` was not run.

### Get / variant lists

**`--get`** writes lead, novel, or proxy results to `--out` / `--output` and exits (does not run the same path as `--extract` file filtering). **`--extract`** is only a variant-ID list for in-pipeline filtering (see [Variant filters](#variant-filters-plink-style)).

Extract lead or novel variants:

```bash
# Lead variants
gwaslab --input sumstats.tsv --get lead --out lead.tsv

# Novel variants with GWAS Catalog EFO trait(s)
gwaslab --input sumstats.tsv --get novel --efo EFO_0004340 --out novel.tsv
```

**`--get` options:**

| Option | Description | Default |
|--------|-------------|---------|
| `--get` | `lead`, `novel`, or `proxy` — write results to `--out` / `--output`, then exit (`proxy` not yet implemented) | None |
| `--sig-level-extract` | P-value threshold for `--get` | `5e-8` |
| `--windowsizekb` | Lead-variant window size (kb) | `500` |
| `--efo` | One or more EFO IDs for novel extraction | None |
| `--only-novel` | Return only truly novel hits in novel extraction | False |

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

For complete workflow examples, see [CLI Workflow Examples](CLIWorkflowExamples.md).

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

## Processing Options

### Execution Order

When multiple flags are provided in a single command, GWASLab executes operations in the following order:

1. **Optional fixes**: `--fix-chr`, `--fix-pos`, `--fix-chr-pos`, `--fix-chr-pos-allele`, `--fix-allele`, `--fix-id`
2. **QC**: `--qc`, `--remove`, `--remove-dup`, `--normalize`
3. **Fill statistics**: `--fill`, `--fill-overwrite`, `--fill-extreme`, `--fill-only-sig`, `--fill-df`
4. **Region filter**: `--filter-region`
5. **Variant filters**: `--extract`, `--exclude`, `--extract-bed`, `--exclude-bed`, `--chr`, `--maf`, `--max-maf`, `--mac`, `--snps-only`, `--min-info`
6. **Harmonization**: `--harmonize` (and related `--ref-*` options)
7. **rsID assignment**: `--assign-rsid`
8. **rsID to CHR:POS conversion**: `--rsid-to-chrpos`
9. **Build inference**: `--infer-build`
10. **Liftover**: `--liftover`
11. **Plotting**: `--plot`
12. **Association extraction**: `--get`
13. **Final export**: `to_format` output (`--out` / `--output` with `--to-fmt`)

### Early-Exit Rules

- `--plot` writes the requested figure and exits immediately.  
  If `--plot` is present, `--get` and final `to_format` export are not executed in that run.
- `--get` writes extraction results and exits immediately.  
  If `--get` is present, final `to_format` export is not executed in that run.

For exact runtime behavior, always verify with `gwaslab --help`.

## Supported Formats

GWASLab supports many input and output formats through the [formatbook](https://github.com/Cloufield/formatbook) repository. Common formats include:

- **Input**: `auto`, `gwaslab`, `plink`, `ldsc`, `vcf`, and many more
- **Output**: `gwaslab`, `ldsc`, `plink`, `plink2`, `saige`, `fastgwa`, `regenie`, `vcf`, and many more

Check the formatbook repository for the complete list of supported formats.

