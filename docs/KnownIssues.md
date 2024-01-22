# Known issues

## P value conversion during extracting lead variants 
- gwaslab<=3.4.18 : In older versions, gwaslab will convert MLOG10P to P (if P column is not avaiable) and then using P to extract lead variants using `.get_lead()`. But during conversion, since MLOG10P was set as float32, and then P was accordinglly set as float32. The minimum of float32 is around 1e-38, below which P values will become 0 and cause gwaslab to extract wrong lead variants. 
- Solution: simply use scaled=True to directly use MLOG10P instead of P values.

