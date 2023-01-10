# Download reference data

## Downloading and file management system

|functions|options|description|
|-|-|-|
|gl.check_available_ref()|| list current available reference data gwaslab can use |
|gl.check_downloaded_ref()||list downloaded reference data |
|gl.download_ref(name)|name|download the refernce data|
|gl.remove_file(name)|name|remove the downloaded data|
|gl.get_path(name)|name|get the path to the refernce data|

Currently, you can download the following files using gwaslab:

- '1kg_eas_hg19'
- '1kg_eas_hg19_tbi'
- '1kg_eur_hg19'
- '1kg_eur_hg19_tbi'
- '1kg_eas_hg38'
- '1kg_eas_hg38_tbi'
- '1kg_eur_hg38'
- '1kg_eur_hg38_tbi'
- 'dbsnp_v151_hg19'
- 'dbsnp_v151_hg38'
- 'ucsc_genome_hg19'
- 'ucsc_genome_hg38'
- '1kg_dbsnp151_hg19_auto'
- '1kg_dbsnp151_hg38_auto'
- 'recombination_hg19'
- 'ensembl_hg19_gtf'
- 'ensembl_hg19_gtf_protein_coding'
- 'ensembl_hg38_gtf'
- 'ensembl_hg38_gtf_protein_coding'
- 'refseq_hg19_gtf'
- 'refseq_hg19_gtf_protein_coding'
- 'refseq_hg38_gtf'
- 'refseq_hg38_gtf_protein_coding'

## Configurations

gwaslab uses 3 files and a default path for reference management:

- `config` : a dictionary of `keyword` : `local path` for local file management. 
- `reference` : a dictionary of `keyword` : `url` for automatically downloading reference files.
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
 


