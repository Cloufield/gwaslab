# Inferring and Checking Allele Frequencies

This document describes the functions in GWASLab for inferring and checking allele frequencies (AF) in summary statistics using reference VCF files.

## Overview

GWASLab provides several functions to work with allele frequencies:

1. **Infer EAF from reference VCF**: `infer_af()` / `infer_af2()`
2. **Check EAF against reference**: `check_af()` / `check_af2()`
3. **Infer EAF from MAF using reference**: `infer_eaf_from_maf()` / `infer_eaf_from_maf2()`

Each function has two versions:
- **Version 1 (normal mode)**: Per-variant VCF queries using tabix
- **Version 2 (sweep mode)**: Bulk lookup table extraction + in-memory matching

Both versions produce **identical results** but differ in performance characteristics.

## Table of Contents

- [Infer EAF from Reference VCF](#infer-eaf-from-reference-vcf)
- [Check EAF Against Reference](#check-eaf-against-reference)
- [Infer EAF from MAF Using Reference](#infer-eaf-from-maf-using-reference)
- [How Allele Orientation is Handled](#how-allele-orientation-is-handled)
- [Version Comparison](#version-comparison)
- [Best Practices](#best-practices)

---

## Infer EAF from Reference VCF

### `infer_af()` / `infer_af2()`

Infers Effect Allele Frequency (EAF) in summary statistics by matching variants with a reference VCF file and extracting the alternative allele frequency (ALT_AF).

**When to use:**
- When sumstats are missing EAF values but have valid CHR, POS, EA, NEA
- When you want to fill in EAF from a reference population (e.g., 1000 Genomes, gnomAD)
- As a preprocessing step before strand inference or allele frequency comparison

**Basic usage:**

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("sumstats.txt.gz")

# Infer EAF using reference VCF (normal mode)
mysumstats.infer_af(
    ref_infer="reference.vcf.gz",
    ref_alt_freq="AF"  # INFO field name in VCF
)

# Infer EAF using reference VCF (sweep mode - faster for large datasets)
mysumstats.infer_af2(
    vcf_path="reference.vcf.gz",
    ref_alt_freq="AF"
)
```

**How it works:**

1. **Matching variants**: Matches variants in sumstats with reference VCF by CHR:POS:EA:NEA
2. **Allele orientation detection**: 
   - Checks if EA matches ALT in VCF → uses ALT_AF directly
   - Checks if EA matches REF in VCF → uses 1-ALT_AF (flipped)
3. **EAF assignment**: Updates EAF column with inferred values

**Example:**

```python
# Sumstats variant: CHR=7, POS=123456, EA=C, NEA=A, EAF=NaN
# VCF record: CHR=7, POS=123456, REF=A, ALT=C, AF=0.3

# Result: EAF=0.3 (EA=C matches ALT=C, so use ALT_AF directly)

# If VCF record was: REF=C, ALT=A, AF=0.7
# Result: EAF=0.3 (EA=C matches REF=C, so use 1-ALT_AF = 1-0.7 = 0.3)
```

**Parameters:**

- `ref_infer` / `vcf_path`: Path to reference VCF/BCF file (must be tabix-indexed for ver1)
- `ref_alt_freq`: Field name for alternative allele frequency in VCF INFO section (e.g., "AF", "AF_popmax", "gnomAD_AF")
- `force`: If True, infer EAF for all variants regardless of STATUS codes

!!! note "Chromosome Conversion"
    Chromosome format conversion is handled automatically by `ChromosomeMapper`. No manual `chr_dict` parameter is needed for high-level functions like `infer_af()`, `check_af()`, etc.

**Returns:**

- Updates the EAF column in the sumstats DataFrame
- Reports number of variants with EAF successfully inferred
- Reports number of variants still missing EAF (not found in reference)

---

## Check EAF Against Reference

### `check_af()` / `check_af2()`

Calculates the difference between EAF in sumstats and ALT_AF in reference VCF (DAF = Difference in Allele Frequency).

**When to use:**
- Quality control: Identify variants with large EAF discrepancies
- Population differences: Compare sumstats EAF with reference population
- After EAF inference: Verify that inferred EAF values are reasonable

**Basic usage:**

```python
# Check EAF difference (normal mode)
mysumstats.check_af(
    ref_infer="reference.vcf.gz",
    ref_alt_freq="AF"
)

# Check EAF difference (sweep mode)
mysumstats.check_af2(
    vcf_path="reference.vcf.gz",
    ref_alt_freq="AF"
)
```

**How it works:**

1. **Matching variants**: Matches variants in sumstats with reference VCF by CHR:POS:EA:NEA
2. **Allele orientation adjustment**: 
   - If EA matches ALT: uses ALT_AF directly
   - If EA matches REF: uses 1-ALT_AF (adjusted to EA frequency)
3. **DAF calculation**: `DAF = EAF (sumstats) - EA_frequency (reference)`

**Output:**

- Creates a `DAF` column with the difference
- Positive DAF: sumstats EAF > reference EA frequency
- Negative DAF: sumstats EAF < reference EA frequency
- Large absolute DAF values may indicate:
  - Population differences
  - Data quality issues
  - Allele mismatches

**Example:**

```python
# Sumstats: EAF=0.4
# Reference: ALT_AF=0.3, EA matches ALT
# DAF = 0.4 - 0.3 = 0.1

# If EA matches REF instead:
# Reference: ALT_AF=0.7, EA matches REF → EA_freq = 1-0.7 = 0.3
# DAF = 0.4 - 0.3 = 0.1 (same result after adjustment)
```

**Note:** DAF here is **not** the derived allele frequency. It's simply the difference between sumstats EAF and reference EA frequency.

---

## Infer EAF from MAF Using Reference

### `infer_eaf_from_maf()` / `infer_eaf_from_maf2()`

Infers EAF from Minor Allele Frequency (MAF) in sumstats using reference VCF ALT frequency to determine allele orientation.

**When to use:**
- When sumstats have MAF values but are missing EAF values
- When you want to infer EAF from MAF using a reference population
- As an alternative to `infer_af()` when you have MAF but not direct EAF information

**Basic usage:**

```python
# Infer EAF from MAF (normal mode)
mysumstats.infer_eaf_from_maf(
    ref_infer="reference.vcf.gz",
    ref_alt_freq="AF"
)

# Infer EAF from MAF (sweep mode)
mysumstats.infer_eaf_from_maf2(
    vcf_path="reference.vcf.gz",
    ref_alt_freq="AF"
)
```

**How it works:**

1. **Extract reference AF**: Gets ALT_AF from reference VCF for each variant
2. **Adjust for allele orientation**: 
   - If EA matches ALT: uses ALT_AF directly
   - If EA matches REF: uses 1-ALT_AF (adjusted to EA frequency)
3. **Compare with MAF**: 
   - If reference EA frequency and MAF suggest different major/minor alleles → flip (use 1-MAF)
   - Otherwise → use MAF directly

**Example:**

```python
# Sumstats: MAF=0.2
# Reference: ALT_AF=0.8, EA matches ALT → EA_freq = 0.8
# Since EA_freq=0.8 (major) and MAF=0.2 (minor), they agree
# → EAF = MAF = 0.2 (EA is minor)

# If reference: ALT_AF=0.2, EA matches ALT → EA_freq = 0.2
# Since EA_freq=0.2 (minor) and MAF=0.2 (minor), they agree
# → EAF = MAF = 0.2 (EA is minor)

# If reference: ALT_AF=0.8, EA matches REF → EA_freq = 0.2
# Since EA_freq=0.2 (minor) and MAF=0.2 (minor), they agree
# → EAF = MAF = 0.2 (EA is minor)
```

**Parameters:**

- `ref_infer` / `vcf_path`: Path to reference VCF/BCF file
- `ref_alt_freq`: Field name for alternative allele frequency in VCF INFO section
- `maf`: Column name for minor allele frequency in sumstats (default: "MAF")
- `force`: If True, infer EAF for all variants regardless of STATUS codes

---

## How Allele Orientation is Handled

The functions automatically detect and handle cases where alleles in sumstats are in opposite orientation compared to the reference VCF:

- **When EA matches ALT**: Uses ALT_AF directly for EA frequency
- **When EA matches REF**: Uses 1-ALT_AF for EA frequency (since ALT_AF is the frequency of NEA)

**Example:**
```python
Sumstats: EA=C, NEA=A
VCF:      REF=A, ALT=C, ALT_AF=0.3
→ EA=C matches ALT=C
→ EAF = ALT_AF = 0.3

Sumstats: EA=C, NEA=A
VCF:      REF=C, ALT=A, ALT_AF=0.7
→ EA=C matches REF=C
→ ALT_AF=0.7 is the frequency of ALT=A (which is NEA)
→ EAF = 1 - ALT_AF = 1 - 0.7 = 0.3 (frequency of EA=C)
```

This ensures that **EAF always represents EA frequency** regardless of how alleles are encoded in the reference VCF.

---

## Version Comparison

Both versions produce **identical results** but differ in performance:

| Aspect | Version 1 (Normal Mode) | Version 2 (Sweep Mode) |
|--------|------------------------|------------------------|
| **Method** | Per-variant VCF queries | Bulk lookup table + in-memory matching |
| **Memory** | Low | Higher (loads lookup table) |
| **Speed** | Slower for large datasets | Faster for large datasets |
| **Allele orientation handling** | Implicit | Automatic |
| **Best for** | Small datasets (< 10K variants) | Large datasets (> 10K variants) |

**Recommendation:** Use Version 2 (`*_af2()`) for datasets with > 10,000 variants. Both versions are equivalent in results (100% match rate tested).

---

## Best Practices

1. **Preprocessing**: Run `harmonize()` before inferring AF to ensure alleles are standardized
2. **Reference selection**: Use population-matched reference VCF (e.g., EAS for East Asian studies)
3. **Version choice**: Use Version 2 (`*_af2()`) for datasets with > 10,000 variants
4. **Quality control**: After inference, use `check_af()` to identify variants with large DAF values
5. **Workflow**: Harmonize → Infer EAF → Check DAF

**Example workflow:**

```python
import gwaslab as gl

mysumstats = gl.Sumstats("sumstats.txt.gz")
mysumstats.harmonize(basic_check=True)
mysumstats.infer_af2(vcf_path="reference.vcf.gz", ref_alt_freq="AF")
mysumstats.check_af2(vcf_path="reference.vcf.gz", ref_alt_freq="AF")
```

---

## Common Issues and Solutions

### Issue: No EAF values inferred

**Possible causes:**
- Variants not in reference VCF
- Chromosome naming mismatch (ChromosomeMapper should handle this automatically)
- Alleles don't match (check harmonization)
- `ref_alt_freq` field missing in VCF

**Solutions:**
- Check if variants are in reference VCF region
- Verify that ChromosomeMapper correctly detected chromosome format (check `mysumstats.mapper._reference_format`)
- Run `harmonize()` to standardize alleles
- Verify `ref_alt_freq` field exists in VCF INFO

### Issue: EAF values seem inconsistent with reference

**This is normal!** Allele orientation differences between sumstats and reference VCF are automatically detected and handled. The functions ensure EAF always represents EA frequency regardless of orientation.

### Issue: Large DAF values

**Possible causes:**
- Population differences (expected)
- Data quality issues
- Allele mismatches

**Solutions:**
- Check if large DAF is consistent across variants
- Verify harmonization status
- Compare with known population frequencies
- Investigate specific variants with very large DAF

---

## See Also

- [Harmonization Documentation](Harmonize.md): Standardizing alleles before AF inference
- [Strand Inference](InferStrand.md): Inferring strand orientation
- [Quality Control](QC.md): General quality control procedures

---

## Function Reference

### `infer_af()` / `infer_af2()`

**Parameters:**
- `ref_infer` / `vcf_path`: Path to reference VCF/BCF file
- `ref_alt_freq`: INFO field name for ALT frequency (default: "AF")
- `force`: Process all variants regardless of STATUS (default: False)

!!! note "Chromosome Conversion"
    Chromosome format conversion is handled automatically by `ChromosomeMapper`. No manual `chr_dict` parameter is needed.

**Returns:**
- Updates EAF column in sumstats
- Reports inference statistics

### `check_af()` / `check_af2()`

**Parameters:**
- `ref_infer` / `vcf_path`: Path to reference VCF/BCF file
- `ref_alt_freq`: INFO field name for ALT frequency (default: "AF")
- `column_name`: Name for DAF column (default: "DAF")
- `force`: Process all variants regardless of STATUS (default: False)

**Returns:**
- Creates DAF column with EAF differences
- Reports DAF statistics (mean, max, min, etc.)

### `infer_eaf_from_maf()` / `infer_eaf_from_maf2()`

**Parameters:**
- `ref_infer` / `vcf_path`: Path to reference VCF/BCF file
- `ref_alt_freq`: INFO field name for ALT frequency (default: "AF")
- `maf`: Column name for MAF in sumstats (default: "MAF")
- `force`: Process all variants regardless of STATUS (default: False)

**Returns:**
- Updates EAF column in sumstats
- Reports inference statistics

