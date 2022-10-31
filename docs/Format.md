# Output sumstats in certain formats

```
to_format(
          path="./sumstats",
          fmt="ldsc",   
          extract=None,
          exclude=None,
          id_use="rsID",
          hapmap3=False,
          exclude_hla=False,  
          build="19", 
          verbose=True,
          output_log=True,
          to_csvargs={},
          float_formats={},
          xymt_number=False,
          xymt=["X","Y","MT"],
          chr_prefix=None,
          bgzip=False,
          tabix=False
          )
```
- `path` : `string`, the path for the output, only prefix is needed.
- `fmt`="ldsc": output format for sumstats. Currently support `plink` ,`plink2`, `ldsc`, `saige`, `fastgwa`, `regenie` . For details , please check [ https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook).
- `extract`:`list`, a list of variants to include.
- `exclude`:`list`, a list of variants to exclude.
- `id_use`:`string`, specify which ID to use when merging with Hapmap3 SNPs.
- `hapmap3`:`boolean` , if True, only output Hapmap3 SNPs.
- `exclude_hla`:`boolean` , if True, exclude variants in MHC region from output.
- `build` : `string`, reference genome build. 
- `xymt_number` : if True, output chrX/Y/MT as 23/24/25.
- `xymt` : `list`, descript how to convert chromosome 23,24,25.
- `chr_prefix` : `string`, add a prefix to chromosomes when output chr. 6 -> Chr6.
- `bgzip` : `boolean`, if True, bgzip the output file. Only works for bed format.
- `tabix` : `boolean`, if True, use tabix to index the bgzipped output file. Only works for bed format.
- `to_csvargs` : `dict` , extra parameters for pd.to_csv()
- `float_formats` : `dict`, a dictionary to specify the float format for each column.
- `verbose` : `boolean`, if True, print logs.
- `output_log` : `boolean`, if True, save log to a file.

## Example 1:
```
import gwaslab as gl
# load your raw sumstats
mysumstats = gl.Sumstats(...)
# basic QC
mysumstats.basic_check()
# output metal format
mysumstats.to_format("./test",fmt="metal")
```

log :
```
Tue Sep 13 18:00:41 2022 Start to format the output sumstats in:  metal  format
Tue Sep 13 18:00:41 2022  -Formatting statistics ...
Tue Sep 13 18:00:41 2022  - Float statistics formats:
Tue Sep 13 18:00:41 2022   - Columns: ['EAF', 'BETA', 'SE', 'P']
Tue Sep 13 18:00:41 2022   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
Tue Sep 13 18:00:41 2022  - Start outputting sumstats in metal format...
Tue Sep 13 18:00:41 2022  -metal format will be loaded...
Tue Sep 13 18:00:41 2022  -metal format meta info:
Tue Sep 13 18:00:41 2022   - format_name  :  metal
Tue Sep 13 18:00:41 2022   - format_source  :  https://genome.sph.umich.edu/wiki/METAL_Documentation
Tue Sep 13 18:00:41 2022   - format_version  :  20220726
Tue Sep 13 18:00:41 2022  -gwaslab to metal format dictionary:
Tue Sep 13 18:00:41 2022   - gwaslab keys: ['SNPID', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'DIRECTION']
Tue Sep 13 18:00:41 2022   - metal values: ['MarkerName', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'Direction']
Tue Sep 13 18:00:41 2022  -Output columns: Index(['MarkerName', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr',
       'P-value'],
      dtype='object')
Tue Sep 13 18:00:41 2022  -Output path: ./test.metal.tsv.gz
Tue Sep 13 18:00:41 2022  -Saving log file: ./test.metal.log
Tue Sep 13 18:00:41 2022 Finished outputting successfully!
```

## Example 2: LDSC format, extract hapmap3 SNPs and exclude SNPs in HLA region
```
## format the sumstats to ldsc format
## extract only hapmap3 SNPs
## exclude SNPs in HLA region
mysumstats.to_format("./test",fmt="ldsc", hapmap3=True, exclude_hla=False, build="19")
```

## Example 3: bed-like format
```
# output 1-based bed-like files for vep
mysumstats.to_format("./test",fmt="vep",xymt_number=True,chr_prefix="Chr")

# output 0-based bed-like file, and then bgzip and index the file.
mysumstats.to_format("./test",fmt="bed",bgzip=True,tabix=True)
```

## Example 4: vcf format
```
# output vcf file, and then bgzip and index the file.
mysumstats.to_format("./test",fmt="vcf",bgzip=True,tabix=True)
```

## Example 5: GWAS-ssf
```
# output  GWAS-ssf format
mysumstats.to_format("./test",fmt="ssf")
```

For sample codes, please check [https://github.com/Cloufield/gwaslab/blob/main/examples/format_load_save.ipynb](https://github.com/Cloufield/gwaslab/blob/main/examples/format_load_save.ipynb)
