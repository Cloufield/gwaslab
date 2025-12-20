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
| `0` | rsID valid & SNPID valid |
| `1` | rsID valid & SNPID invalid |
| `2` | rsID invalid & SNPID valid |
| `3` | rsID invalid & SNPID invalid |
| `4` | rsID valid & SNPID valid (alternative) |
| `5` | rsID valid & SNPID unknown |
| `6` | rsID unknown & SNPID valid |
| `7` | rsID invalid & SNPID unknown |
| `8` | rsID unknown & SNPID invalid |
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

### Digit 5: EA & NEA Standardization

| Code | Description |
|------|-------------|
| `0` | Standardized SNP |
| `1` | Standardized & normalized insertion |
| `2` | Standardized & normalized deletion |
| `3` | Standardized & normalized indel |
| `4` | Standardized indel |
| `5` | Indistinguishable or not normalized allele |
| `6` | Invalid allele notation |
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

This returns a DataFrame with columns:
- **Genome_Build**: Reference genome build information
- **rsID&SNPID**: Validation status of ID fields
- **CHR&POS**: Chromosome/position validation status
- **Stadardize&Normalize**: Standardization status
- **Align**: Reference alignment status
- **Panlidromic_SNP&Indel**: Variant type characteristics
- **Count**: Absolute count of each status code
- **Percentage(%)**: Relative percentage of total variants

![image](https://user-images.githubusercontent.com/40289485/211861846-8309be9f-05ea-456e-ad8a-3265872826f9.png)


## Reference Table

For a complete visual reference of all status code combinations:

![image](https://user-images.githubusercontent.com/40289485/196681586-eb79707c-d866-4393-a0f5-b825686c9b04.png)

## Notes

- Status codes are automatically assigned during data loading and processing
- Codes are updated during harmonization, standardization, and quality control steps
- The STATUS column is a reserved header in GWASLab and it's automatically created and maintained
- Status codes are stored as integers but can be converted to strings for digit extraction
- Higher digit values generally indicate lower confidence or more uncertainty about that aspect of the variant

