# Commonly used data in gwaslab
## Chromosme notation conversion dictionary
- `get_chr_list()` : get a full list of chromosomes
- chr/ number <-> NC dictionary
    - `get_chr_to_NC(build)` : 
    - `get_number_to_NC(build)`
    - `get_NC_to_chr(build)`
    - `get_NC_to_number(build)`
- chr <-> number 
  - `get_chr_to_number()`
  - `get_number_to_chr()`
## High-ld region bed
- `get_high_ld(build)`

## Sumstats Format list
- `get_format_dict(fmt)`
- `get_formats_list()`

## Gene gtf
- `get_gtf(chrom, build="19",source="ensembl")`

## gwaslab information
- `gwaslab_info()`
