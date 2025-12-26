# QC Report Generation

The QC report generation feature provides a comprehensive quality control report for GWAS summary statistics, combining basic QC checks, lead variant extraction, visualization, and summary statistics into a single HTML or PDF document.

## Overview

The `generate_qc_report()` function and `Sumstats.report()` method generate a complete QC report that includes:

- **Basic QC Summary**: Quality control statistics and checks
- **Lead Variants Table**: Genome-wide significant variants identified using sliding window approach
- **MQQ Plot**: Manhattan and QQ plots for genome-wide visualization
- **Regional Plots**: Locus-specific plots for each lead variant (optional, requires VCF reference)
- **Summary Statistics**: Comprehensive dataset summary including chromosome distribution, missing values, MAF distribution, p-value statistics, and variant status
- **Processing Log**: Complete log of all operations performed

## Basic Usage

### Using the Sumstats.report() method

```python
import gwaslab as gl

# Load your sumstats
mysumstats = gl.Sumstats(
    "path/to/sumstats.txt.gz",
    snpid="SNP",
    chrom="CHR",
    pos="POS",
    ea="ALT",
    nea="REF",
    neaf="Frq",
    beta="BETA",
    se="SE",
    p="P",
    build="19",
    n="N"
)

# Generate HTML report
mysumstats.report("qc_report.html")
```

### Using the standalone function

```python
import gwaslab as gl

mysumstats = gl.Sumstats("path/to/sumstats.txt.gz", ...)

# Generate report using the function
gl.generate_qc_report(
    mysumstats,
    output_path="qc_report.html"
)
```

## Report Formats

The report format is determined by the file extension of `output_path`:

- **HTML**: `.html` extension (default)
- **PDF**: `.pdf` extension (requires `weasyprint`)

```python
# Generate HTML report
mysumstats.report("report.html")

# Generate PDF report (requires weasyprint)
mysumstats.report("report.pdf")
```

### PDF Library Requirements

For PDF generation, `weasyprint` must be installed:

```bash
pip install weasyprint
```

On some systems, you may also need to install system dependencies:
```bash
# Ubuntu/Debian:
sudo apt-get install python3-weasyprint
```

If `weasyprint` is not available, the report will be generated as HTML instead.

## Customizing the Report

### Basic QC Options

```python
mysumstats.report(
    "report.html",
    basic_check_kwargs={
        "remove": False,          # Don't remove variants during QC
        "remove_dup": False,      # Don't remove duplicates
        "normalize": True,        # Normalize variants
        "verbose": False
    }
)
```

### Lead Variant Extraction Options

```python
mysumstats.report(
    "report.html",
    get_lead_kwargs={
        "sig_level": 5e-8,       # Significance threshold
        "windowsizekb": 500,     # Window size in kb for clumping
        "anno": False,           # Don't annotate variants
        "verbose": False
    }
)
```

### MQQ Plot Options

```python
mysumstats.report(
    "report.html",
    mqq_plot_kwargs={
        "save": True,
        "dpi": 300,
        "verbose": False
    }
)
```

### Regional Plot Options

Regional plots require a VCF reference file for LD information:

```python
mysumstats.report(
    "report.html",
    regional_plot_kwargs={
        "vcf_path": "/path/to/reference.vcf.gz",
        "region_size": 500,      # Region size in kb
        "save": True,
        "dpi": 300,
        "verbose": False
    }
)
```

### Custom Report Title

```python
mysumstats.report(
    "report.html",
    report_title="My Custom QC Report"
)
```

## Complete Example

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats(
    "examples/0_sample_data/t2d_bbj.txt.gz",
    snpid="SNP",
    chrom="CHR",
    pos="POS",
    ea="ALT",
    nea="REF",
    neaf="Frq",
    beta="BETA",
    se="SE",
    p="P",
    build="19",
    n="N"
)

# Generate comprehensive QC report with regional plots
report_path = mysumstats.report(
    output_path="t2d_qc_report.pdf",
    report_title="T2D GWAS QC Report",
    basic_check_kwargs={
        "remove": False,
        "remove_dup": False,
        "normalize": True,
        "verbose": False
    },
    get_lead_kwargs={
        "sig_level": 5e-8,
        "windowsizekb": 500,
        "anno": False,
        "verbose": False
    },
    mqq_plot_kwargs={
        "save": True,
        "dpi": 300,
        "verbose": False
    },
    regional_plot_kwargs={
        "vcf_path": "/path/to/reference.vcf.gz",
        "region_size": 500,
        "save": True,
        "dpi": 300,
        "verbose": False
    }
)

print(f"Report generated: {report_path}")
```

## Report Contents

### 1. Basic QC Summary

The report includes a summary of the basic QC checks performed, including:
- Total variants before and after QC
- Variants removed and reasons
- Data quality metrics

### 2. Lead Variants Table

A table of all genome-wide significant lead variants, including:
- SNPID/rsID
- Chromosome and position
- Effect and non-effect alleles
- P-value and -log10(P)
- Effect size (BETA) and standard error (SE)

### 3. MQQ Plot

A combined Manhattan and QQ plot showing:
- Genome-wide association signals
- Expected vs observed p-values
- Lambda GC statistic

### 4. Regional Plots

For each lead variant, a regional plot showing:
- Association signals in the region
- LD information (if VCF reference provided)
- Gene annotations
- Recombination rate

### 5. Summary Statistics

Comprehensive dataset summary including:
- **Overview**: Row/column counts, column names, QC status
- **Chromosome Distribution**: Variant counts per chromosome
- **Missing Values**: Missing data counts by column
- **MAF Distribution**: Minor allele frequency categories
- **P-value Statistics**: Minimum p-value, counts at significance thresholds
- **Variant Status**: QC status code distribution
- **Variant Metadata**: Total variants, minimum p-value, minimum MAF
- **Sample Metadata**: Sample size statistics

### 6. Processing Log

Complete log of all operations performed, including:
- Timestamps for each operation
- Data shape changes
- QC steps performed
- Warnings and errors

## Method Chaining

The `report()` method supports method chaining:

```python
mysumstats = gl.Sumstats("data.txt.gz", ...)

result = (mysumstats
    .basic_check(remove=False, normalize=True)
    .get_lead(sig_level=5e-8)
    .report("report.html")
)
```

## Output Files

The report generation creates:

1. **Main report file**: HTML or PDF report
2. **Plots directory**: `{report_name}_plots/` containing:
   - `mqq_plot.png`: Manhattan and QQ plot
   - `regional_plot_chr{chr}_{pos}.png`: Regional plots for each lead variant

## Performance Considerations

- **Caching**: GTF and recombination rate data are cached in memory for faster repeated access
- **Regional Plots**: Generating regional plots for many lead variants can be time-consuming
- **PDF Generation**: PDF conversion may take additional time for large reports

## Notes

- If no lead variants are found at the specified significance threshold, the report will still be generated with available information
- Regional plots are only generated if a VCF reference path is provided
- The report includes all available information even if some steps fail (e.g., if regional plots cannot be generated)

