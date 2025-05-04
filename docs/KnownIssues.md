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