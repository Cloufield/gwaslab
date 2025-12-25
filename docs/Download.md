# Download reference data

## Downloading and file management system

| Functions                   | Options | Datatype | Description                                       |
|-----------------------------|---------|----------|---------------------------------------------------|
| `gl.check_available_ref()`  | `show_all=False`, `verbose=True` | `boolean` | list available reference datasets GWASLab can use |
| `gl.check_downloaded_ref()` |         |          | list downloaded reference dataset                 |
| `gl.download_ref(keyword)`  | `keyword`, `directory=None`, `local_filename=None`, `overwrite=False` | `string` | download the reference dataset                    |
| `gl.remove_file(keyword)`   | keyword | `string` | remove the downloaded dataset                     |
| `gl.get_path(keyword)`      | keyword, `verbose=True` | `string` | get the path to the reference dataset              |
| `gl.set_default_directory(path)` | path | `string` | set the default directory for downloaded files |
| `gl.get_default_directory()` | | | get the default directory path |
| `gl.scan_downloaded_files()` | `verbose=True` | `boolean` | scan data directory for unregistered files and match with available references |
| `gl.add_local_data(keyword, local_path, ...)` | keyword, local_path, format, description, md5sum, suggested_use, tbi, csi | | register a local file without downloading |
| `gl.remove_local_record(keyword)` | keyword | `string` | remove a local data record from configuration (does not delete the file) |

You can download the following files using GWASLab with the keywords:

!!! error "Processed datasets are currently hosted on Dropbox which may not be accessible for users in certain regions."

Datasets you need to download explicitly if needed.

- **`1kg_eas_hg19`, `1kg_eas_hg19_tbi`**: Autosomes; 1KGp3v5 low-coverage EAS VCF and index (hg19) - *Processed, Dropbox*
- **`1kg_eur_hg19`, `1kg_eur_hg19_tbi`**: Autosomes; 1KGp3v5 low-coverage EUR VCF and index (hg19) - *Processed, Dropbox*
- **`1kg_amr_hg19`, `1kg_amr_hg19_tbi`**: Autosomes; 1KGp3v5 low-coverage AMR VCF and index (hg19) - *Processed, Dropbox*
- **`1kg_afr_hg19`, `1kg_afr_hg19_tbi`**: Autosomes; 1KGp3v5 low-coverage AFR VCF and index (hg19) - *Processed, Dropbox*
- **`1kg_sas_hg19`, `1kg_sas_hg19_tbi`**: Autosomes; 1KGp3v5 low-coverage SAS VCF and index (hg19) - *Processed, Dropbox*
- **`1kg_pan_hg19`, `1kg_pan_hg19_tbi`**: Autosomes; 1KGp3v5 low-coverage PAN VCF and index (hg19); all ancestries - *Processed, Dropbox*
- **`1kg_eas_hg38`, `1kg_eas_hg38_tbi`**: Autosomes; 1KGp3v5 30x EAS VCF and index (hg38) - *Processed, Dropbox*
- **`1kg_eur_hg38`, `1kg_eur_hg38_tbi`**: Autosomes; 1KGp3v5 30x EUR VCF and index (hg38) - *Processed, Dropbox*
- **`1kg_afr_hg38`, `1kg_afr_hg38_tbi`**: Autosomes; 1KGp3v5 30x AFR VCF and index (hg38) - *Processed, Dropbox*
- **`1kg_amr_hg38`, `1kg_amr_hg38_tbi`**: Autosomes; 1KGp3v5 30x AMR VCF and index (hg38) - *Processed, Dropbox*
- **`1kg_sas_hg38`, `1kg_sas_hg38_tbi`**: Autosomes; 1KGp3v5 30x SAS VCF and index (hg38) - *Processed, Dropbox*
- **`1kg_pan_hg38`, `1kg_pan_hg38_tbi`**: Autosomes; 1KGp3v5 30x PAN VCF and index (hg38); all ancestries - *Processed, Dropbox*
- **`dbsnp_v151_hg19`**: dbSNP v151 (hg19, !!very large) - *Original source*
- **`dbsnp_v151_hg38`**: dbSNP v151 (hg38, !!very large) - *Original source*
- **`dbsnp_v157_hg19`**: dbSNP v157 (hg19, !!very large) - *Original source*
- **`dbsnp_v157_hg38`**: dbSNP v157 (hg38, !!very large) - *Original source*
- **`ucsc_genome_hg19`**: UCSC human reference genome (hg19) - *Original source*
- **`ucsc_genome_hg38`**: UCSC human reference genome (hg38) - *Original source*
- **`1kg_dbsnp151_hg19_auto`**: Autosomes; 1KGp3v5 variants SNPID-rsID conversion table (hg19) - *Processed, Dropbox*
- **`1kg_dbsnp151_hg38_auto`**: Autosomes; 1KGp3v5 variants SNPID-rsID conversion table (hg38) - *Processed, Dropbox*
- **`1kg_hm3_hg19_eaf`**: 1KG HapMap3 hg19 EAF (GRCh37) provides allele frequency estimates for the HapMap3 populations on the hg19 reference genome, derived from 1000 Genomes Project data - *Processed, Dropbox*
- **`1kg_hm3_hg38_eaf`**: 1KG HapMap3 hg38 EAF (GRCh38) provides allele frequency estimates for the HapMap3 populations on the hg38 reference genome, derived from 1000 Genomes Project data - *Processed, Dropbox*

!!! note "tbi index file will be automatically downloaded when you download VCF files"

!!! example "Download and check the local path of `1kg_eas_hg19`"
    ```python
    gl.download_ref("1kg_eas_hg19")

    gl.get_path("1kg_eas_hg19")
    '/Users/he/work/gwaslab/src/gwaslab/data/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz'
    ```

!!! example "Download to a custom directory"
    ```python
    gl.set_default_directory("/path/to/custom/directory")
    gl.download_ref("1kg_eas_hg19")
    ```

!!! example "Register a local file"
    ```python
    gl.add_local_data(
        keyword="my_custom_ref",
        local_path="/path/to/my/file.vcf.gz",
        format="vcf",
        description="My custom reference file",
        suggested_use="Harmonization"
    )
    ```

Other reference datasets GWASLab uses to create regional plot and annotate variants (GWASLab will automatically download and process these files when needed): 

| Keyword              | Description                                                   | Note               |
|----------------------|---------------------------------------------------------------|--------------------|
| `recombination_hg19` | Recombination rate reference files from Hapmap Project (hg19) | Processed, Dropbox |
| `recombination_hg38` | Recombination rate reference files from Hapmap Project (hg38) | Processed, Dropbox |
| `ensembl_hg19_gtf`   | GTF file for genes from Ensembl (hg19)                        | Original source    |
| `ensembl_hg38_gtf`   | GTF file for genes from Ensembl (hg38)                        | Original source    |
| `refseq_hg19_gtf`    | GTF file for genes from Refseq (hg19)                         | Original source    |
| `refseq_hg38_gtf`    | GTF file for genes from Refseq (hg38)                         | Original source    |


## Configurations

GWASLab uses 3 files and a default path for reference management:

- `config` : a dictionary of `keyword`:`local path` for local file management.
- `reference` : a dictionary of `keyword`:`url` for automatically downloading reference files.
- `formatbook` : a dictionary used for format header conversions. 
- `data_directory`: the path for downloaded reference file. default: (`~/.gwaslab`)

You can use `gl.options.paths` to check the paths of the three files.

```python
gl.options.paths
{'config': '/Users/he/work/gwaslab/src/gwaslab/data/config.json',
 'reference': '/Users/he/work/gwaslab/src/gwaslab/data/reference.json',
 'formatbook': '/Users/he/work/gwaslab/src/gwaslab/data/formatbook.json',
 'data_directory': '/Users/he/.gwaslab/'}
```

Sometimes you might need to use your own files, which can be done using `gl.options.set_option(key, newpath)` (simply run this after loading the gwaslab package):

```python
# change the path for formatbook

gl.options.set_option("formatbook","/newpath/formatbook.json")
gl.options.paths
{'config': '/Users/he/work/gwaslab/src/gwaslab/data/config.json',
 'reference': '/Users/he/work/gwaslab/src/gwaslab/data/reference.json',
 'formatbook': '/newpath/formatbook.json',
 'data_directory': '/Users/he/.gwaslab/'}
```
 


