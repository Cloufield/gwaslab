# Download reference data

## Downloading and file management system

|functions|options|description|
|-|-|-|
|`gl.check_available_ref()`|| list current available reference data GWASLab can use |
|`gl.check_downloaded_ref()`||list downloaded reference data |
|`gl.download_ref(keyword)`|keyword|download the refernce data|
|`gl.remove_file(keyword)`|keyword|remove the downloaded data|
|`gl.get_path(keyword)`|keyword|get the path to the refernce data|

Currently, you can download the following files using GWASLab with the keywords:

!!! error "Processed datasets are currently hosted on Dropbox which may not be accessible for users in certain regions."

|Keyword|Description|Note|
|-|-|-|
|'1kg_eas_hg19','1kg_eas_hg19_tbi'|1KGp3v5 low-coverage EAS VCF and index (hg19)|Processed, Dropbox|
|'1kg_eur_hg19','1kg_eur_hg19_tbi'|1KGp3v5 low-coverage EUR VCF and index (hg19)|Processed, Dropbox|
|'1kg_eas_hg38','1kg_eas_hg38_tbi'|1KGp3v5 30x EAS VCF and index (hg38)|Processed, Dropbox|
|'1kg_eur_hg38','1kg_eur_hg38_tbi'|1KGp3v5 30x EUR VCF and index (hg38)|Processed, Dropbox|
|'dbsnp_v151_hg19'|dbSNP151 (hg19, !!very large)|Original source|
|'dbsnp_v151_hg38'|dbSNP151 (hg38, !!very large)|Original source|
|'ucsc_genome_hg19'|UCSC human reference genome (hg19)|Original source|
|'ucsc_genome_hg38'|UCSC human reference genome (hg38)|Original source|
|'1kg_dbsnp151_hg19_auto'|1KGp3v5 variants SNPID-rsID conversion table (hg19)|Processed, Dropbox|
|'1kg_dbsnp151_hg38_auto'|1KGp3v5 variants SNPID-rsID conversion table (hg38)|Processed, Dropbox|
|'recombination_hg19'|Recombination rate reference files from Hapmap Project (hg19)|Processed, Dropbox|
|'recombination_hg38'|Recombination rate reference files from Hapmap Project (hg38)|Processed, Dropbox|
|'ensembl_hg19_gtf'|GTF file for genes from Ensembl (hg19)|Original source|
|'ensembl_hg38_gtf'|GTF file for genes from Ensembl (hg38)|Original source|
|'refseq_hg19_gtf'|GTF file for genes from Refseq (hg19)|Original source|
|'refseq_hg38_gtf'|GTF file for genes from Refseq (hg19)|Original source|

## Configurations

GWASLab uses 3 files and a default path for reference management:

- `config` : a dictionary of `keyword`:`local path` for local file management.
- `reference` : a dictionary of `keyword`:`url` for automatically downloading reference files.
- `formatbook` : a dictionary used for format header conversions. 
- `data_directory`: the path for downloaded reference file. default: (`~/.gwaslab`)

You can use `gl.options.paths` to check the paths of the three files.

```
gl.options.paths
{'config': '/Users/he/work/gwaslab/src/gwaslab/data/config.json',
 'reference': '/Users/he/work/gwaslab/src/gwaslab/data/reference.json',
 'formatbook': '/Users/he/work/gwaslab/src/gwaslab/data/formatbook.json',
 'data_directory': '/Users/he/.gwaslab/'}
```

Sometimes you might need to use your own files, which can be done using `gl.options.set_option(key, newpath)` (simply run this after loadnig the gwaslab package):

```
# change the path for formatbook

gl.options.set_option("formatbook","/newpath/formatbook.json")
gl.options.paths
{'config': '/Users/he/work/gwaslab/src/gwaslab/data/config.json',
 'reference': '/Users/he/work/gwaslab/src/gwaslab/data/reference.json',
 'formatbook': '/newpath/formatbook.json',
 'data_directory': '/Users/he/.gwaslab/'}
```
 


