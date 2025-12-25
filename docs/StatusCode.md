# Status Code

## Overview

GWASLab uses a **7-digit status code** system to track the quality and processing status of each variant in your summary statistics. This code provides a comprehensive, traceable record of variant validation, standardization, and harmonization.

### Design Principles

- **Traceable**: Every digit encodes specific information about variant processing
- **Higher value → higher uncertainty**: Lower digit values generally indicate better quality/confidence
- **Comprehensive**: Covers genome build, IDs, coordinates, alleles, alignment, and variant type

### Status Code Structure

Each status code is a 7-digit integer where:

| Digit Position | Description | Values |
|----------------|-------------|--------|
| **1-2** | Genome build | 13, 19, 38, 97, 98, 99 |
| **3** | rsID & SNPID validation | 0-9 |
| **4** | CHR & POS validation | 0, 2-9 |
| **5** | EA & NEA standardization | 0-9 |
| **6** | REF-NEA alignment | 0-9 |
| **7** | Palindromic SNPs & Indels | 0-9 |

## Digit Meanings

### Digits 1-2: Genome Build

| Code | Build | Description |
|------|-------|-------------|
| `13` | CHM13 | CHM13/T2T reference genome |
| `19` | hg19 | GRCh37/hg19 reference genome |
| `38` | hg38 | GRCh38/hg38 reference genome |
| `97` | Unmapped | Variant could not be mapped to any build |
| `98` | Unknown | Genome build is unknown |
| `99` | Unchecked | Genome build has not been checked |

### Digit 3: rsID & SNPID Validation

| Code | Description |
|------|-------------|
| `0` | rsid valid & SNPID valid |
| `1` | rsid valid & SNPID invalid |
| `2` | rsid invalid & SNPID valid |
| `3` | rsid invalid & SNPID invalid |
| `4` | rsid valid & SNPID valid |
| `5` | rsid valid & SNPID unknown |
| `6` | rsid unknown & SNPID valid |
| `7` | rsid invalid & SNPID unknown |
| `8` | rsid unknown & SNPID invalid |
| `9` | Unchecked |

### Digit 4: CHR & POS Validation

| Code | Description |
|------|-------------|
| `0` | CHR valid & POS valid |
| `2` | CHR invalid & POS invalid |
| `3` | CHR invalid & POS valid |
| `4` | CHR valid & POS invalid |
| `5` | CHR valid & POS unknown |
| `6` | CHR unknown & POS valid |
| `7` | CHR invalid & POS unknown |
| `8` | CHR unknown & POS invalid |
| `9` | Unchecked |

**Note:** Code `1` is not used in this digit.

### Digit 5: EA & NEA Standardization

| Code | Description |
|------|-------------|
| `0` | standardized SNP |
| `1` | standardized & normalized insertion |
| `2` | standardized & normalized deletion |
| `3` | standardized & normalized indel |
| `4` | standardized indel |
| `5` | indistinguishable or not normalized allele |
| `6` | invalid allele notation |
| `7` | Unknown |
| `9` | Unchecked |

### Digit 6: REF-NEA Alignment

| Code | Description |
|------|-------------|
| `0` | Match: NEA=REF |
| `1` | Flipped_fixed |
| `2` | Reverse_complementary_fixed |
| `3` | Flipped |
| `4` | Reverse_complementary |
| `5` | Reverse_complementary+Flipped |
| `6` | Both_alleles_on_ref+indistinguishable |
| `8` | Not_on_reference_genome |
| `9` | Unchecked |

**Note:** Code `7` is not used in this digit.

### Digit 7: Palindromic SNPs & Indels

| Code | Description |
|------|-------------|
| `0` | Not_palindromic_SNPs |
| `1` | Palindromic+strand |
| `2` | Palindromic-strand_fixed |
| `3` | Indel_match |
| `4` | Indel_flipped_fixed |
| `5` | Palindromic-strand |
| `6` | Indel_flipped |
| `7` | Indistinguishable |
| `8` | No_matching_or_no_info |
| `9` | Unchecked |

## Usage

### View Status Summary

The `summary()` method provides an overview of status codes in your dataset:

```python
mysumstats.summary()
```

This displays:
- Total variant counts
- Status code distribution
- Percentage breakdown by status category

![image](https://user-images.githubusercontent.com/40289485/211861728-bc56389d-84a9-4929-95ab-0517ac063dfd.png)

### Look Up Status Codes

Use `lookup_status()` to decode and analyze all status codes in your dataset:

```python
status_df = mysumstats.lookup_status()
print(status_df)
```

**Returns:** `pandas.DataFrame`

The method returns a pandas DataFrame with columns:
- **Genome_Build**: Reference genome build information (CHM13, hg19, hg38, Unmapped, Unknown, or Unchecked)
- **rsID&SNPID**: Validation status of ID fields
- **CHR&POS**: Chromosome/position validation status
- **Stadardize&Normalize**: Standardization and normalization status
- **Align**: Reference alignment status
- **Panlidromic_SNP&Indel**: Variant type characteristics
- **Count**: Absolute count of each status code
- **Percentage(%)**: Relative percentage of total variants



## Notes

- Status codes are automatically assigned during data loading and processing
- Codes are updated during harmonization, standardization, and quality control steps
- The STATUS column is a reserved header in GWASLab and it's automatically created and maintained
- Status codes are stored as integers but can be converted to strings for digit extraction
- Higher digit values generally indicate lower confidence or more uncertainty about that aspect of the variant

