# HLA region

For post-GWAS analyses, we often exclude the HLA (Human Leukocyte Antigen) region due to its exceptional complexity of linkage disequilibrium (LD). The HLA region on chromosome 6 contains the most polymorphic genes in the human genome, with extensive LD patterns that can confound association analyses.

In GWASLab, `exclude_hla()` is designed to remove variants in the extended MHC (xMHC) region by default.

## Default exclusion range: chr6:25-34 Mb

By default, `exclude_hla()` removes variants in **chr6:25,000,000 ~ 34,000,000** (25-34 Mb) for both hg19 and hg38.

### Why 25-34 Mb instead of exact gene boundaries?

The actual gene boundaries for the extended MHC region are:

- **xMHC (HIST1H2AA ~ RPL12P1)**:
  - hg38: 25,726,063 ~ 33,400,644
  - hg19: 25,726,291 ~ 33,368,421

However, GWASLab uses the simplified **25-34 Mb range** for several practical reasons:

1. **Safety margin**: The rounded boundaries ensure complete coverage of the complex LD region, accounting for:
   - Variants near gene boundaries that may still be in strong LD
   - Coordinate differences between genome builds
   - Edge cases where variants might be slightly outside exact gene boundaries

2. **Consistency**: The same range (25-34 Mb) works for both hg19 and hg38, simplifying usage and avoiding build-specific confusion

3. **Standard practice**: Many GWAS analyses and tools use rounded boundaries (e.g., 25-34 Mb) for HLA exclusion, making results more comparable across studies

4. **Practical simplicity**: Round numbers are easier to remember, communicate, and use in analyses

5. **Comprehensive coverage**: The 25-34 Mb range fully encompasses the extended MHC region while providing a conservative buffer zone

## Region definitions

### Extended MHC (xMHC) region

The extended MHC region spans from HIST1H2AA to RPL12P1, covering approximately 7.6 Mb.

**Reference**: Horton, R., Wilming, L., Rand, V., Lovering, R. C., Bruford, E. A., Khodiyar, V. K., ... & Beck, S. (2004). Gene map of the extended human MHC. Nature Reviews Genetics, 5(12), 889-899.

**Exact gene boundaries**:
- Start: [HIST1H2AA](https://www.ncbi.nlm.nih.gov/gene/221613)
- End: [RPL12P1](https://www.ncbi.nlm.nih.gov/gene/729727)
- hg38: 25,726,063 ~ 33,400,644
- hg19: 25,726,291 ~ 33,368,421
- T2T-CHM13v2.0: 25,591,841 ~ 33,222,006

**GWASLab default**: chr6:25,000,000 ~ 34,000,000 (mode="xmhc")

### Classical HLA region

The classical HLA region spans from GABBR1 to KIFC1, covering approximately 3.78 Mb.

**Reference**: Shiina, T., Hosomichi, K., Inoko, H., & Kulski, J. K. (2009). The HLA genomic loci map: expression, interaction, diversity and disease. Journal of human genetics, 54(1), 15-39.

**Exact gene boundaries**:
- Start: [GABBR1](https://www.ncbi.nlm.nih.gov/gene/2550)
- End: [KIFC1](https://www.ncbi.nlm.nih.gov/gene/3833)
- hg38: 29,602,238 ~ 33,409,896
- hg19: 29,570,015 ~ 33,377,673
- T2T-CHM13v2.0: 29,476,949 ~ 33,231,258

**GWASLab default**: chr6:29,500,000 ~ 33,500,000 (mode="hla" or mode="mhc")

## Usage

```
# Exclude extended MHC region (default, 25-34 Mb)
mysumstats.exclude_hla()

# Exclude classical HLA region only (29.5-33.5 Mb)
mysumstats.exclude_hla(mode="hla")

# Custom range
mysumstats.exclude_hla(lower=25000000, upper=34000000)

# Specify genome build (uses build-specific ranges)
mysumstats.exclude_hla(build="38", mode="xmhc")
```

!!! note "Mode options"
    - `mode="xmhc"` (default): Extended MHC region, chr6:25-34 Mb
    - `mode="hla"` or `mode="mhc"`: Classical HLA region, chr6:29.5-33.5 Mb