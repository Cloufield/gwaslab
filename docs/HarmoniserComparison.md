# GWASLab vs OpenTargets Genetics-Sumstat-Harmoniser Comparison

## Overview

This document compares the harmonized output from **GWASLab** and the [OpenTargets genetics-sumstat-harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser) when processing the same input data.

**Note:** All input variants are assumed to be already aligned with the reference FASTA (NEA matches REF). The comparison uses **maf_threshold=0.42** and **ref_maf_threshold=0.42** for both tools (GWASLab default is 0.4, but 0.42 is used here to match the harmoniser's default threshold).

**Key difference in approach:** GWASLab performs a two-step check: (1) align NEA/EA against the reference FASTA sequence (`check_ref`, setting STATUS digit 6), then (2) use the reference VCF for strand inference of palindromic SNPs and indistinguishable indels (`infer_strand`, setting STATUS digit 7). The OpenTargets harmoniser only checks the reference VCF, assuming NEA already matches REF.

## Input Data

Test dataset with 13 variants covering different harmonization scenarios:

| Case | SNPID | CHR | POS | EA | NEA | EAF | BETA | SE | P | Description |
|------|-------|-----|-----|----|----|-----|------|----|---|-------------|
| 1 | 1:1000_A_G | 1 | 1000 | G | A | 0.30 | 0.1 | 0.01 | 1e-5 | Non-palindromic, correct |
| 2 | 1:2000_C_T | 1 | 2000 | T | C | 0.40 | 0.1 | 0.01 | 1e-5 | Non-palindromic, correct |
| 3 | 1:3000_G_A | 1 | 3000 | A | G | 0.80 | 0.1 | 0.01 | 1e-5 | Non-palindromic, correct |
| 4 | 1:4000_T_C | 1 | 4000 | C | T | 0.35 | 0.1 | 0.01 | 1e-5 | Non-palindromic, correct |
| 5 | 1:5000_A_T | 1 | 5000 | T | A | 0.25 | 0.1 | 0.01 | 1e-5 | Palindromic A/T, forward strand |
| 6 | 1:6000_G_C | 1 | 6000 | C | G | 0.28 | 0.1 | 0.01 | 1e-5 | Palindromic G/C, forward strand |
| 7 | 1:7000_A_T | 1 | 7000 | T | A | 0.25 | 0.1 | 0.01 | 1e-5 | Palindromic A/T, reverse strand |
| 8 | 1:8000_G_C | 1 | 8000 | C | G | 0.28 | 0.1 | 0.01 | 1e-5 | Palindromic G/C, reverse strand |
| 9 | 1:9000_A_T | 1 | 9000 | T | A | 0.55 | 0.1 | 0.01 | 1e-5 | Palindromic, high MAF |
| 10 | 1:10000_A_AT | 1 | 10000 | AT | A | 0.30 | 0.1 | 0.01 | 1e-5 | Indel, in VCF |
| 11 | 1:20000_A_G | 1 | 20000 | G | A | 0.30 | 0.1 | 0.01 | 1e-5 | **Non-palindromic, NOT in VCF** |
| 12 | 1:21000_A_T | 1 | 21000 | T | A | 0.25 | 0.1 | 0.01 | 1e-5 | **Palindromic A/T, NOT in VCF** |
| 13 | 1:22000_A_AT | 1 | 22000 | AT | A | 0.30 | 0.1 | 0.01 | 1e-5 | **Indel (both alleles on ref), NOT in VCF** |

### Reference VCF

The reference VCF contains 10 variants used for harmonization:

| CHROM | POS | ID | REF | ALT | AF_NFE |
|-------|-----|----|----|----|--------|
| 1 | 1000 | rs1 | A | G | 0.30 |
| 1 | 2000 | rs2 | C | T | 0.40 |
| 1 | 3000 | rs3 | G | A | 0.20 |
| 1 | 4000 | rs4 | T | C | 0.35 |
| 1 | 5000 | rs5 | A | T | 0.25 |
| 1 | 6000 | rs6 | G | C | 0.28 |
| 1 | 7000 | rs7 | A | T | 0.75 |
| 1 | 8000 | rs8 | G | C | 0.72 |
| 1 | 9000 | rs9 | A | T | 0.45 |
| 1 | 10000 | rs10 | A | AT | 0.30 |

**Note:** Variants at positions 20000, 21000, and 22000 are **not** in the reference VCF. The FASTA reference has `AT` at positions 22000-22001, so both alleles of the Case 13 indel (NEA=A and EA=AT) match the reference, making it an indistinguishable indel (digit 6 = 6).

## Harmonized Output Comparison

### Variants with Matching Results (9/13)

| SNPID | CHR | POS | GWASLab EA | GWASLab NEA | GWASLab EAF | GWASLab BETA | Harmoniser EA | Harmoniser NEA | Harmoniser EAF | Harmoniser BETA | Match |
|-------|-----|-----|------------|-------------|-------------|--------------|---------------|-----------------|----------------|-----------------|-------|
| 1:1000_A_G | 1 | 1000 | G | A | 0.30 | 0.1 | G | A | 0.30 | 0.1 | ✓ |
| 1:2000_C_T | 1 | 2000 | T | C | 0.40 | 0.1 | T | C | 0.40 | 0.1 | ✓ |
| 1:3000_G_A | 1 | 3000 | A | G | 0.80 | 0.1 | A | G | 0.80 | 0.1 | ✓ |
| 1:4000_T_C | 1 | 4000 | C | T | 0.35 | 0.1 | C | T | 0.35 | 0.1 | ✓ |
| 1:5000_A_T | 1 | 5000 | T | A | 0.25 | 0.1 | T | A | 0.25 | 0.1 | ✓ |
| 1:6000_G_C | 1 | 6000 | C | G | 0.28 | 0.1 | C | G | 0.28 | 0.1 | ✓ |
| 1:7000_A_T | 1 | 7000 | T | A | 0.75 | -0.1 | T | A | 0.75 | -0.1 | ✓ |
| 1:8000_G_C | 1 | 8000 | C | G | 0.72 | -0.1 | C | G | 0.72 | -0.1 | ✓ |
| 1:10000_A_AT | 1 | 10000 | AT | A | 0.30 | 0.1 | AT | A | 0.30 | 0.1 | ✓ |

### Variants with Different Results (4/13)

#### 1. Non-palindromic variant not in reference VCF (Case 11)

**Position:** 1:20000

| Tool | EA | NEA | EAF | BETA | STATUS | hm_code |
|------|----|----|-----|------|--------|---------|
| GWASLab | G | A | 0.30 | 0.1 | 9960000 (digit 6=0, digit 7=0) | — |
| Harmoniser | — | — | — | — | — | 15 (not in VCF) |

**Difference:** GWASLab preserves this variant but its STATUS code (`9960000`) is **indistinguishable from a valid matched variant** (e.g., Cases 1-4). Digit 6 = 0 because NEA matches the FASTA reference at that position, and digit 7 = 0 because non-palindromic SNPs skip the VCF lookup entirely. The harmoniser identifies it as absent from the reference VCF (hm_code=15) and drops it.

**This is a known gap:** non-palindromic variants not in the reference VCF cannot be filtered by STATUS code.

#### 2. Palindromic variant not in reference VCF (Case 12)

**Position:** 1:21000 (A/T palindromic, EAF=0.25)

| Tool | EA | NEA | EAF | BETA | STATUS | hm_code |
|------|----|----|-----|------|--------|---------|
| GWASLab | T | A | 0.25 | 0.1 | 9960008 (digit 6=0, digit 7=8) | — |
| Harmoniser | — | — | — | — | — | 15 (not in VCF) |

**Difference:** GWASLab flags this variant with digit 7 = 8 ("No_matching_or_no_info") because the VCF-based strand inference found no match. This variant **can** be filtered by STATUS code. The harmoniser drops it (hm_code=15).

#### 3. Indistinguishable indel not in reference VCF (Case 13)

**Position:** 1:22000 (indel, both alleles on reference genome, EAF=0.30)

| Tool | EA | NEA | EAF | BETA | STATUS | hm_code |
|------|----|----|-----|------|--------|---------|
| GWASLab | AT | A | 0.30 | 0.1 | 9960368 (digit 6=6, digit 7=8) | — |
| Harmoniser | — | — | — | — | — | 15 (not in VCF) |

**Difference:** GWASLab identifies this as an indistinguishable indel (digit 6 = 6, both alleles match the reference) and flags it with digit 7 = 8 after VCF lookup finds no match. This variant **can** be filtered by STATUS code. The harmoniser drops it (hm_code=15).

#### 4. Palindromic SNP with high MAF (Case 9)

**Position:** 1:9000 (A/T palindromic, EAF=0.55, MAF=0.45 > 0.42 threshold)

| Tool | EA | NEA | EAF | BETA | STATUS | hm_code |
|------|----|----|-----|------|--------|---------|
| GWASLab | T | A | 0.55 | 0.1 | 9960007 (digit 7=7) | — |
| Harmoniser | — | — | — | — | — | 18 (high MAF) |

**Difference:** Both tools agree this variant's strand cannot be inferred. GWASLab keeps it with digit 7 = 7 indicating MAF too high for inference. The harmoniser drops it (hm_code=18).

## Variants Not in Reference VCF — STATUS Code Behavior

A key finding from this comparison is how different variant types are handled when they are **absent from the reference VCF**:

| Case | Variant Type | GWASLab STATUS | Digit 6 | Digit 7 | VCF Queried? | Filterable? |
|------|-------------|----------------|---------|---------|--------------|-------------|
| 11 | Non-palindromic SNP | 9960000 | 0 (Match) | 0 (Not palindromic) | **No** | **No** |
| 12 | Palindromic SNP | 9960008 | 0 (Match) | 8 (No match/info) | Yes | Yes |
| 13 | Indistinguishable indel | 9960368 | 6 (Both on ref) | 8 (No match/info) | Yes | Yes |

- **Palindromic SNPs and indistinguishable indels** correctly get digit 7 = 8 because their code paths query the reference VCF.
- **Non-palindromic SNPs** get digit 7 = 0 unconditionally (without VCF lookup), making them indistinguishable from valid matched variants. This is a known gap.

## Key Differences

1. **Non-palindromic variants not in reference VCF:** GWASLab preserves them but STATUS code does not reflect absence from VCF (digit 7 hardcoded to 0). Harmoniser drops them (hm_code=15). 
2. **Palindromic variants not in reference VCF:** GWASLab correctly flags them with digit 7 = 8. Harmoniser drops them (hm_code=15).
3. **Indistinguishable indels not in reference VCF:** GWASLab correctly flags them with digit 6 = 6 and digit 7 = 8. Harmoniser drops them (hm_code=15).
4. **High MAF palindromic SNPs:** GWASLab preserves them with digit 7 = 7, harmoniser drops them (hm_code=18).
5. **Standard variants:** Both tools produce identical harmonized results.
6. **Low MAF palindromic SNPs:** Both tools produce identical harmonized results with correct strand inference.

## Test Output

Full comparison results are saved to: `test/output/gwaslab_harmoniser_toy_comparison.tsv`

To run the consistency test:

```bash
python test/test_gwaslab_harmoniser_consistency.py
```

## References

- **GWASLab:** [https://cloufield.github.io/gwaslab/](https://cloufield.github.io/gwaslab/)
- **OpenTargets Harmoniser:** [https://github.com/opentargets/genetics-sumstat-harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)
