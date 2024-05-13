
# HLA region

(20240513)

For post-GWAS analyses, we often exclude HLA region due to its complexity of LD.

In gwaslab, `exclude_hla` is designed to remove variants in the extended MHC (xMHC) region.

From a practical point of view, by default, `exclude_hla` will remove variants in chr6:25mb~34mb (for both hg38 and hg19).

## xMHC region : HIST1H2AA ~ 7.6mb ~ RPL12P1

- reference: Horton, R., Wilming, L., Rand, V., Lovering, R. C., Bruford, E. A., Khodiyar, V. K., ... & Beck, S. (2004). Gene map of the extended human MHC. Nature Reviews Genetics, 5(12), 889-899.
- range: start of [HIST1H2AA](https://www.ncbi.nlm.nih.gov/gene/221613) ~ end of [RPL12P1](https://www.ncbi.nlm.nih.gov/gene/729727)
    - hg38:  25,726,063 ~ 33,400,644
    - hg19 : 25,726,291 ~ 33,368,421
    - T2T-CHM13v2.0 : 25,591,841 ~ 33,222,006

## HLA region : GABBR1 ~ 3.78mb ~ KIFC1

- reference: Shiina, T., Hosomichi, K., Inoko, H., & Kulski, J. K. (2009). The HLA genomic loci map: expression, interaction, diversity and disease. Journal of human genetics, 54(1), 15-39. 
- range: start of [GABBR1](https://www.ncbi.nlm.nih.gov/gene/2550) ~ end of [KIFC1](https://www.ncbi.nlm.nih.gov/gene/3833)
    - hg38:  29,602,238 ~ 33,409,896
    - hg19:  29,570,015 ~ 33,377,673
    - T2T-CHM13v2.0 : 29,476,949 ~ 33,231,258