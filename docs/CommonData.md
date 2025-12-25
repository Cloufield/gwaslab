# Commonly used data in GWASLab

GWASLab integrates a few pre-defined datasets for quick conversion and access, which are used in the functions of GWASLab.

## Chromosome notation conversion dictionary

### Full chromosome list

`gl.get_chr_list()`: Get a full list of chromosomes (string datatype):

**Parameters:**

- `add_number=False`: If True, include both string and numeric representations
- `n=25`: Maximum chromosome number to include
- `only_number=False`: If True, return only numeric chromosome numbers

**Examples:**
```python
gl.get_chr_list()

['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','X','Y','M','MT']

gl.get_chr_list(n=2)

['1', '2', 'X', 'Y', 'M', 'MT']

gl.get_chr_list(only_number=True)

[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
```

### CHR - NC conversion

### CHR - NC conversion

These functions convert between chromosome identifiers and NCBI accession IDs.

**`gl.get_chr_to_NC(build, inverse=False)`**: Convert chromosome string to NCBI accession ID

**Parameters:**

- `build`: Genome build version ('19' for hg19/GRCh37 or '38' for hg38/GRCh38)
- `inverse`: If True, return NCBI ID to chromosome mapping (equivalent to `get_NC_to_chr`)

**Example:**
```python
gl.get_chr_to_NC(build="19")

{'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9', '6': 'NC_000006.11', '7': 'NC_000007.13', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10', '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9', '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10', '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}
```

**`gl.get_number_to_NC(build, inverse=False)`**: Convert chromosome number (int) to NCBI accession ID

**Parameters:**

- `build`: Genome build version ('19' or '38')
- `inverse`: If True, return NCBI ID to chromosome number mapping (equivalent to `get_NC_to_number`)

**Example:**
```python
gl.get_number_to_NC(build="19")

{1: 'NC_000001.10', 2: 'NC_000002.11', 3: 'NC_000003.11', 4: 'NC_000004.11', 5: 'NC_000005.9', 6: 'NC_000006.11', 7: 'NC_000007.13', 8: 'NC_000008.10', 9: 'NC_000009.11', 10: 'NC_000010.10', 11: 'NC_000011.9', 12: 'NC_000012.11', 13: 'NC_000013.10', 14: 'NC_000014.8', 15: 'NC_000015.9', 16: 'NC_000016.9', 17: 'NC_000017.10', 18: 'NC_000018.9', 19: 'NC_000019.9', 20: 'NC_000020.10', 21: 'NC_000021.8', 22: 'NC_000022.10', 23: 'NC_000023.10', 24: 'NC_000024.9', 25: 'NC_012920.1'}
```

**`gl.get_NC_to_chr(build)`**: Convert NCBI accession ID to chromosome string

**Parameters:**

- `build`: Genome build version ('19' or '38')

**Example:**
```python
gl.get_NC_to_chr(build="19")

{'NC_000001.10': '1', 'NC_000002.11': '2', 'NC_000003.11': '3', 'NC_000004.11': '4', 'NC_000005.9': '5', 'NC_000006.11': '6', 'NC_000007.13': '7', 'NC_000008.10': '8', 'NC_000009.11': '9', 'NC_000010.10': '10', 'NC_000011.9': '11', 'NC_000012.11': '12', 'NC_000013.10': '13', 'NC_000014.8': '14', 'NC_000015.9': '15', 'NC_000016.9': '16', 'NC_000017.10': '17', 'NC_000018.9': '18', 'NC_000019.9': '19', 'NC_000020.10': '20', 'NC_000021.8': '21', 'NC_000022.10': '22', 'NC_000023.10': 'X', 'NC_000024.9': 'Y', 'NC_012920.1': 'MT'}
```

**`gl.get_NC_to_number(build)`**: Convert NCBI accession ID to chromosome number (int)

**Parameters:**

- `build`: Genome build version ('19' or '38')

**Example:**
```python
gl.get_NC_to_number(build="19")

{'NC_000001.10': 1, 'NC_000002.11': 2, 'NC_000003.11': 3, 'NC_000004.11': 4, 'NC_000005.9': 5, 'NC_000006.11': 6, 'NC_000007.13': 7, 'NC_000008.10': 8, 'NC_000009.11': 9, 'NC_000010.10': 10, 'NC_000011.9': 11, 'NC_000012.11': 12, 'NC_000013.10': 13, 'NC_000014.8': 14, 'NC_000015.9': 15, 'NC_000016.9': 16, 'NC_000017.10': 17, 'NC_000018.9': 18, 'NC_000019.9': 19, 'NC_000020.10': 20, 'NC_000021.8': 21, 'NC_000022.10': 22, 'NC_000023.10': 23, 'NC_000024.9': 24, 'NC_012920.1': 25}
```

### CHR Datatype conversion (string <-> int) 

**`gl.get_chr_to_number(out_chr=False, xymt=["X","Y","MT"], xymt_num=[23,24,25])`**: Convert chromosome string to number

**Parameters:**

- `out_chr`: If True, returns dictionary with string keys and string values
- `xymt`: List of non-numeric chromosome identifiers (default: ["X","Y","MT"])
- `xymt_num`: Corresponding numeric values for xymt (default: [23,24,25])

**Example:**
```python
gl.get_chr_to_number()

{'1': 1, '2': 2, ..., 'X': 23, 'Y': 24, 'MT': 25, 'M': 25}
```

**`gl.get_number_to_chr(in_chr=False, xymt=["X","Y","MT"], xymt_num=[23,24,25], prefix="")`**: Convert chromosome number to string

**Parameters:**

- `in_chr`: If True, returns dictionary with string keys and values
- `xymt`: List of non-numeric chromosome identifiers (default: ["X","Y","MT"])
- `xymt_num`: Corresponding numeric values for xymt (default: [23,24,25])
- `prefix`: Optional prefix for chromosome identifiers

**Example:**
```python
gl.get_number_to_chr()

{1: '1', 2: '2', ..., 23: 'X', 24: 'Y', 25: 'MT'}
```


## High-LD region bed

**`gl.get_high_ld(build="19")`**: Get the path to the BED file of high-LD regions

**Parameters:**

- `build`: Genome build version ('19' for hg19/GRCh37 or '38' for hg38/GRCh38)

**Returns:**
- `str`: Path to the high LD region BED file

**Example:**
```python
gl.get_high_ld("19")

'/home/yunye/anaconda3/envs/gwaslab/lib/python3.8/site-packages/gwaslab/data/high_ld/high_ld_hla_hg19.bed.gz'
```

## Sumstats Format list

**`gl.get_format_dict(fmt, inverse=False)`**: Get format dictionary and metadata for a specified format

**Parameters:**

- `fmt`: Format name to look up in the format book
- `inverse`: If True, return inverted dictionary with value-key mapping

**Returns:**
- `tuple`: (metadata_dict, format_dict) - Metadata and format dictionary mapping fields

**Example:**
```python
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

**`gl.get_formats_list()`**: Get a list of all available formats that GWASLab supports

**Returns:**
- `list`: Format names available in the format book

**Example:**
```python
gl.get_formats_list()

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
