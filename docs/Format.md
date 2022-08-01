# Formatting sumstats

```
sumstats.format()
```

Format the sumstats to the formats that were accepted by commonly used tools including LDSC, MAGMA, FUMA , METAL, (bed ,vcf,...)

`path`: output file path

`format`: Currently, GWASlab support ldsc, (vcf...)

`extract`: a list of SNPIDs.

`exclude`: a list of SNPIDs.

`exclude_hla`: True or False. If True, exclude HLA region when exporting sumstats.

`hapmap3` : True or False. If True,, only exporting Hapmap3 SNPs.

`to_csvargs`: arguments for `pandas.to_csv()` function


```sumstats.to_fmt()```