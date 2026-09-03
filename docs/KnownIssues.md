# Known issues

## P value conversion during extracting lead variants 

- gwaslab<=3.4.18 : In older versions, gwaslab converts MLOG10P to P (if P column is not available) and then uses P to extract lead variants using `.get_lead()`. But during conversion, since MLOG10P was set as float32, and then P was accordingly set as float32. The minimum of float32 is around 1e-38, below which P values will become 0. This could cause error when extracting lead variants. 
- Solution: simply use scaled=True to directly use MLOG10P instead of P values.

## Skipped lead variant in extracting lead variants 

- gwaslab<=3.4.22 : When extracting lead variants using `.get_lead()`, if the lead variant is the last variant in all significant variants, it will be mistakenly skipped.
- Solution: Update to new version of gwaslab.

## Regional plots : color issue

- gwaslab<=3.4.39 : the color assigned to each variant is actually the color for the lower LD r2 category. For example, variants with LD>0.8 will be colored with the color for 0.8>LD>0.6.
- Solution: Update to new version of gwaslab.

## Alignment for unsorted sumstats

-  3.4.43 <= gwaslab <= 3.5.7 : If the coordinates in sumstats are not unsorted, the alignment may fail. 
- Solution: Update to new version of gwaslab.

## Incomplete VCF tabix/CSI index can silently miss rsIDs

- A `.tbi` or `.csi` that was not built from the local VCF can be inconsistent with it. `bcftools view -r 14` may start mid-chromosome (for example around 50 Mb) even though a direct interval query such as `14:30000400` still returns the variant. Sweep-mode `.assign_rsid2()` used whole-chromosome `-r CHR` with `-T` targets, so those prefix variants looked absent from the reference VCF. Downloaded dbSNP v151 VCF/index pairs are one known example.
- This is not a GWASLab bug in the VCF itself. Rebuild the index from the local VCF (`bcftools index -t -f <vcf>`), then re-run assignment. Current GWASLab queries `-r CHR:minPOS-maxPOS` and warns under these conditions (see [Assigning rsID](AssignrsID.md) and [Download reference data](Download.md)):
    - **Lookup stream** (sweep extract, `n_targets > 0`): no records returned for that chromosome, **or** first returned POS is after the smallest target POS and **> 5%** of that chromosome’s targets sit before the first returned record.
    - **Assignment-rate outlier** (after VCF/sweep fill): among chromosomes with **≥ 500** eligible variants (need **≥ 5** such chromosomes), fill rate is **≥ 0.15** below the median (`median_rate − rate ≥ 0.15`).