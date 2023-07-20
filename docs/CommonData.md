# Commonly used data in GWASLab

GWASLab integrates a few pre-defined datasets for quick conversion and access, which are used in the functions of GWASLab.

## Chromosme notation conversion dictionary

### Full chromosome list

`gl.get_chr_list()`: Get a full list of chromosomes (string datatype):

```
gl.get_chr_list()

['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','X','Y','M','MT']


gl.get_chr_list(n=2)

['1', '2', 'X', 'Y', 'M', 'MT']
```

### CHR - NC conversion

`gl.get_chr_to_NC(build)`

```
gl.get_chr_to_NC(build="19")

{'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9', '6': 'NC_000006.11', '7': 'NC_000007.13', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10', '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9', '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10', '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}
```

`gl.get_number_to_NC(build)`

```
gl.get_number_to_NC(build="19")

{1: 'NC_000001.10', 2: 'NC_000002.11', 3: 'NC_000003.11', 4: 'NC_000004.11', 5: 'NC_000005.9', 6: 'NC_000006.11', 7: 'NC_000007.13', 8: 'NC_000008.10', 9: 'NC_000009.11', 10: 'NC_000010.10', 11: 'NC_000011.9', 12: 'NC_000012.11', 13: 'NC_000013.10', 14: 'NC_000014.8', 15: 'NC_000015.9', 16: 'NC_000016.9', 17: 'NC_000017.10', 18: 'NC_000018.9', 19: 'NC_000019.9', 20: 'NC_000020.10', 21: 'NC_000021.8', 22: 'NC_000022.10', 23: 'NC_000023.10', 24: 'NC_000024.9', 25: 'NC_012920.1'}
```

`gl.get_NC_to_chr(build)`

```
gl.get_NC_to_chr(build="19")

{'NC_000001.10': '1', 'NC_000002.11': '2', 'NC_000003.11': '3', 'NC_000004.11': '4', 'NC_000005.9': '5', 'NC_000006.11': '6', 'NC_000007.13': '7', 'NC_000008.10': '8', 'NC_000009.11': '9', 'NC_000010.10': '10', 'NC_000011.9': '11', 'NC_000012.11': '12', 'NC_000013.10': '13', 'NC_000014.8': '14', 'NC_000015.9': '15', 'NC_000016.9': '16', 'NC_000017.10': '17', 'NC_000018.9': '18', 'NC_000019.9': '19', 'NC_000020.10': '20', 'NC_000021.8': '21', 'NC_000022.10': '22', 'NC_000023.10': 'X', 'NC_000024.9': 'Y', 'NC_012920.1': 'MT'}
```

`gl.get_NC_to_number(build)`

```
gl.get_NC_to_number(build="19")

{'NC_000001.10': 1, 'NC_000002.11': 2, 'NC_000003.11': 3, 'NC_000004.11': 4, 'NC_000005.9': 5, 'NC_000006.11': 6, 'NC_000007.13': 7, 'NC_000008.10': 8, 'NC_000009.11': 9, 'NC_000010.10': 10, 'NC_000011.9': 11, 'NC_000012.11': 12, 'NC_000013.10': 13, 'NC_000014.8': 14, 'NC_000015.9': 15, 'NC_000016.9': 16, 'NC_000017.10': 17, 'NC_000018.9': 18, 'NC_000019.9': 19, 'NC_000020.10': 20, 'NC_000021.8': 21, 'NC_000022.10': 22, 'NC_000023.10': 23, 'NC_000024.9': 24, 'NC_012920.1': 25}
```

### CHR Datatype conversion (string <-> int) 

`gl.get_chr_to_number()`: `string` to `int`

`gl.get_number_to_chr()`: `int` to `string`


## High-LD region bed

`gl.get_high_ld(build)`: get the path for the BED file of high-LD regions

```
gl.get_high_ld("19")

'/home/yunye/anaconda3/envs/gwaslab/lib/python3.8/site-packages/gwaslab/data/high_ld/high_ld_hla_hg19.bed.gz'
```

## Sumstats Format list

`gl.get_format_dict(fmt)`: check the details of a format.

```
gl.get_format_dict(fmt="gwaslab")

({'format_name': 'ldsc',
  'format_source': 'https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format',
  'format_source2': 'https://github.com/bulik/ldsc/blob/master/munge_sumstats.py',
  'format_version': 20150306},
 {'SNP': 'rsID',
  'A2': 'NEA',
  'A1': 'EA',
  'Frq': 'EAF',
  'N': 'N',
  'Beta': 'BETA',
  'P': 'P',
  'Z': 'Z',
  'INFO': 'INFO',
  'OR': 'OR',
  'CHR': 'CHR',
  'POS': 'POS'})
```

`gl.get_formats_list()`: show all available format GWASLab supports.

```
['auto',
 'bolt_lmm',
 'fastgwa',
 'gwascatalog',
 'gwascatalog_hm',
 'gwaslab',
 'ldsc',
 'metal',
 'mrmega',
 'mtag',
 'pgscatalog',
 'pgscatalog_hm',
 'pheweb',
 'plink',
 'plink2',
 'regenie',
 'saige',
 'ssf',
 'template',
 'vcf']
```
