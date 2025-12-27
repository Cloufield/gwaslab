# GWASLab vs OpenTargets Genetics-Sumstat-Harmoniser Comparison

## Overview

This document compares the harmonized output from **GWASLab** and the [OpenTargets genetics-sumstat-harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser) when processing the same input data.

**Note:** All input variants are assumed to be already aligned with the reference FASTA (NEA matches REF). The comparison uses **maf_threshold=0.42** and **ref_maf_threshold=0.42** for both tools (GWASLab default is 0.4, but 0.42 is used here to match the harmoniser's default threshold).

## Input Data

Test dataset with 11 variants covering different harmonization scenarios:

| SNPID | CHR | POS | EA | NEA | EAF | BETA | SE | P |
|-------|-----|-----|----|----|-----|------|----|---|
| 1:1000_A_G | 1 | 1000 | G | A | 0.30 | 0.1 | 0.01 | 1e-5 |
| 1:2000_C_T | 1 | 2000 | T | C | 0.40 | 0.1 | 0.01 | 1e-5 |
| 1:3000_G_A | 1 | 3000 | A | G | 0.80 | 0.1 | 0.01 | 1e-5 |
| 1:4000_T_C | 1 | 4000 | C | T | 0.35 | 0.1 | 0.01 | 1e-5 |
| 1:5000_A_T | 1 | 5000 | T | A | 0.25 | 0.1 | 0.01 | 1e-5 |
| 1:6000_G_C | 1 | 6000 | C | G | 0.28 | 0.1 | 0.01 | 1e-5 |
| 1:7000_A_T | 1 | 7000 | T | A | 0.75 | 0.1 | 0.01 | 1e-5 |
| 1:8000_G_C | 1 | 8000 | C | G | 0.72 | 0.1 | 0.01 | 1e-5 |
| 1:9000_A_T | 1 | 9000 | T | A | 0.55 | 0.1 | 0.01 | 1e-5 |
| 1:10000_A_AT | 1 | 10000 | AT | A | 0.30 | 0.1 | 0.01 | 1e-5 |
| 1:20000_A_G | 1 | 20000 | G | A | 0.30 | 0.1 | 0.01 | 1e-5 |

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

**Note:** Variant at position 20000 is not in the reference VCF.

## Harmonized Output Comparison

### Variants with Matching Results (9/11)

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

**Note:** Variants at positions 7000 and 8000 have negative BETA values after harmonization, indicating that alleles were flipped (BETA sign reversed and EAF flipped to 1-EAF).

### Variants with Different Results (2/11)

#### 1. Variant Not in Reference VCF

**Position:** 1:20000

| Tool | EA | NEA | EAF | BETA | Status |
|------|----|----|-----|------|--------|
| GWASLab | G | A | 0.30 | 0.1 | Kept with original data, STATUS code updated to reflect not found in reference |
| Harmoniser | - | - | - | - | Dropped (not found in reference) |

**Difference:** GWASLab preserves variants not found in the reference VCF (with STATUS code updated to indicate this status), while the harmoniser excludes them from the output.

#### 2. Palindromic SNP with High MAF

**Position:** 1:9000 (A/T palindromic, EAF=0.55, MAF=0.45 > 0.42 threshold)

| Tool | EA | NEA | EAF | BETA | Status |
|------|----|----|-----|------|--------|
| GWASLab | T | A | 0.55 | 0.1 | Kept with original data, STATUS code updated to indicate indistinguishable (high MAF) |
| Harmoniser | - | - | - | - | Dropped (high MAF, cannot infer strand) |

**Difference:** GWASLab keeps palindromic SNPs with high MAF (above the inference threshold) with STATUS code updated to reflect that strand cannot be inferred, while the harmoniser drops them because strand orientation cannot be reliably inferred.

## Summary Statistics

| Metric | Match Rate |
|--------|-----------|
| EA matches | 81.8% (9/11) |
| NEA matches | 81.8% (9/11) |
| EAF matches | 81.8% (9/11) |
| BETA matches | 81.8% (9/11) |
| All results match | 81.8% (9/11) |

## Key Differences

1. **Variants not in reference VCF:** GWASLab preserves them with STATUS code updated to reflect this status, harmoniser drops them
2. **High MAF palindromic SNPs:** GWASLab preserves them with STATUS code updated to indicate indistinguishable status, harmoniser drops them
3. **Standard variants:** Both tools produce identical harmonized results
4. **Low MAF palindromic SNPs:** Both tools produce identical harmonized results with correct strand inference

## Test Output

Full comparison results are saved to: `test/output/gwaslab_harmoniser_toy_comparison.tsv`

To run the consistency test:

```bash
python test/test_gwaslab_harmoniser_consistency.py
```

## References

- **GWASLab:** [https://cloufield.github.io/gwaslab/](https://cloufield.github.io/gwaslab/)
- **OpenTargets Harmoniser:** [https://github.com/opentargets/genetics-sumstat-harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)
