# Download reference data

!!! example
    ```python
    import gwaslab as gl
    ```

## Configuration and default paths

Downloaded **files** go under `~/.gwaslab/` (`data_directory`). The local **registry** (`config.json`, mapping keywords → paths) also defaults to `~/.gwaslab/config.json`. On first import, GWASLab migrates a legacy registry from `{package}/data/config.json` if present.

| Key | Default | Override |
|-----|---------|----------|
| `data_directory` | `~/.gwaslab/` | `gl.set_default_directory()`, `gwaslab init --directory DIR`, `gwaslab config set data_directory PATH`, or `GWASLAB_DATA_DIR` |
| `config` | `~/.gwaslab/config.json` | `gl.options.set_option("config", path)`, `gwaslab config set config PATH`, or `GWASLAB_CONFIG` |
| `reference` | package `data/reference.json` | `gl.update_available_ref()` to refresh from GitHub |

User overrides for `data_directory` and `config` are persisted in `~/.gwaslab/settings.json`.

!!! example
    ```python
    gl.set_default_directory("/data/refs/")
    gl.options.paths["config"]  # ~/.gwaslab/config.json by default
    gl.get_path("1kg_eas_hg19")  # resolves from registry; does not download
    ```

CLI:

```bash
gwaslab init                          # create dirs, migrate config, scan default data_directory
gwaslab init --directory /data/refs   # set data_directory (persisted) and scan
gwaslab config set data_directory /data/refs
gwaslab path config                   # show registry JSON path
```

### Cache and storage layout

| Location | Purpose |
|----------|---------|
| `~/.gwaslab/` (`data_directory`) | Downloaded genomic refs, sumstats subdirs (`GCST…/`), recombination extracts |
| `~/.gwaslab/config.json` | Registry: keyword → `{local_path, kind, source, …}` |
| `~/.gwaslab/settings.json` | Persisted path overrides |
| `~/.gwaslab/lookup/` | Harmonize rsID sweep cache (per workflow) |
| `~/.gwaslab/recombination/hg19|hg38/` | Lazy-downloaded recombination maps |
| Package `data/` | Read-only catalog (`reference.json`), built-in HapMap3/chains |
| `platformdirs` user cache (`gwaspipe`) | infer_strand HDF5 cache (outside `data_directory`) |

Registry entries include **`kind`**: `ref` (catalog/local genomic refs) or `sumstats` (GWAS Catalog GCST). Filter with `gwaslab list ref --downloaded --kind ref`.

## Check available reference data

Processed files are hosted on Dropbox. Other files will be downloaded from their original source.

The catalog is nested JSON in `reference.json` (one keyword per entry with `url`, `description`, `suggested_use`, optional `tbi`, `md5sum`).

!!! example
    ```python
    gl.check_available_ref()
    ```

**stdout (abbreviated):**
```
Start to check available reference files...
 - Available keywords: 1kg_eas_hg19 1kg_eur_hg19 ucsc_genome_hg19 ...
Finished checking available reference files...
```

Each keyword maps to metadata such as:

```json
{
  "url": "https://…/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz?dl=1",
  "description": "1000 Genomes Project East Asian (1KG EAS) VCF on hg19 …",
  "suggested_use": "LD reference panel for region plot; infer strand for EAS",
  "tbi": {"url": "https://…/EAS.….vcf.gz.tbi?dl=1"},
  "md5sum": "6162a93cb80935168c0bfa519748b054"
}
```

## Download reference data

GWASLab default directory for saving reference data is `~/.gwaslab`

!!! example
    ```python
    gl.download_ref("testlink")
    ```

**stdout:**
```
Sat Feb  3 13:45:00 2024 Start to download  testlink  ...
Sat Feb  3 13:45:00 2024  -Downloading to: /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz
Sat Feb  3 13:46:24 2024  -Updating record in config file...
Sat Feb  3 13:46:24 2024  -File /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz.tbi exists.
Sat Feb  3 13:46:24 2024  -Updating record in config file...
Sat Feb  3 13:46:24 2024  -Downloading to: /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz.tbi
Sat Feb  3 13:46:24 2024 Downloaded  testlink  successfully!
```

## Check downloaded reference data 

!!! example
    ```python
    gl.check_downloaded_ref()
    ```

**stdout:**
```
Sat Feb  3 13:46:24 2024 Start to check downloaded reference files...
Sat Feb  3 13:46:24 2024  -Checking the config file:/home/yunye/work/gwaslab/src/gwaslab/data/config.json
Sat Feb  3 13:46:24 2024  -Config file exists.
Sat Feb  3 13:46:24 2024  -Updating config.json...
Sat Feb  3 13:46:24 2024   - ensembl_hg19_gtf  :  /home/yunye/.gwaslab/Homo_sapiens.GRCh37.87.chr.gtf.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_hg19  :  /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_hg19_tbi  :  /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - recombination_hg38  :  /home/yunye/.gwaslab/recombination/hg38/recombination_hg38.tar.gz
Sat Feb  3 13:46:24 2024   - ensembl_hg38_gtf  :  /home/yunye/.gwaslab/Homo_sapiens.GRCh38.109.chr.gtf.gz
Sat Feb  3 13:46:24 2024   - ucsc_genome_hg19  :  /home/yunye/.gwaslab/hg19.fa
Sat Feb  3 13:46:24 2024   - ucsc_genome_hg38  :  /home/yunye/.gwaslab/hg38.fa
Sat Feb  3 13:46:24 2024   - refseq_hg19_gtf  :  /home/yunye/.gwaslab/GRCh37_latest_genomic.gtf.gz
Sat Feb  3 13:46:24 2024   - refseq_hg38_gtf  :  /home/yunye/.gwaslab/GRCh38_latest_genomic.gtf.gz
Sat Feb  3 13:46:24 2024   - 1kg_dbsnp151_hg19_auto  :  /home/yunye/.gwaslab/1kg_dbsnp151_hg19_auto.txt.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_x_hg19  :  /home/yunye/.gwaslab/EAS.chrX.split_norm_af.1kgp3v5.vcf.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_x_hg19_tbi  :  /home/yunye/.gwaslab/EAS.chrX.split_norm_af.1kgp3v5.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - 1kg_afr_hg19  :  /home/yunye/.gwaslab/AFR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
Sat Feb  3 13:46:24 2024   - 1kg_afr_hg19_tbi  :  /home/yunye/.gwaslab/AFR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - testlink_tbi  :  /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - testlink  :  /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz
```

```
{'ensembl_hg19_gtf': '/home/yunye/.gwaslab/Homo_sapiens.GRCh37.87.chr.gtf.gz',
 '1kg_eas_hg19': '/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz',
 '1kg_eas_hg19_tbi': '/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz.tbi',
 'recombination_hg38': '/home/yunye/.gwaslab/recombination/hg38/recombination_hg38.tar.gz',
 'ensembl_hg38_gtf': '/home/yunye/.gwaslab/Homo_sapiens.GRCh38.109.chr.gtf.gz',
 'ucsc_genome_hg19': '/home/yunye/.gwaslab/hg19.fa',
 'ucsc_genome_hg38': '/home/yunye/.gwaslab/hg38.fa',
 'refseq_hg19_gtf': '/home/yunye/.gwaslab/GRCh37_latest_genomic.gtf.gz',
 'refseq_hg38_gtf': '/home/yunye/.gwaslab/GRCh38_latest_genomic.gtf.gz',
 '1kg_dbsnp151_hg19_auto': '/home/yunye/.gwaslab/1kg_dbsnp151_hg19_auto.txt.gz',
 '1kg_eas_x_hg19': '/home/yunye/.gwaslab/EAS.chrX.split_norm_af.1kgp3v5.vcf.gz',
 '1kg_eas_x_hg19_tbi': '/home/yunye/.gwaslab/EAS.chrX.split_norm_af.1kgp3v5.vcf.gz.tbi',
 '1kg_afr_hg19': '/home/yunye/.gwaslab/AFR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz',
 '1kg_afr_hg19_tbi': '/home/yunye/.gwaslab/AFR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz.tbi',
 'testlink_tbi': '/home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz.tbi',
 'testlink': '/home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz'}
```

## Remove downloaded reference data

!!! example
    ```python
    gl.remove_file("testlink")
    ```

**stdout:**
```
Sat Feb  3 13:46:24 2024 Start to remove  testlink  ...
Sat Feb  3 13:46:24 2024 Removed : /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz
Sat Feb  3 13:46:24 2024 Start to check downloaded reference files...
Sat Feb  3 13:46:24 2024  -Checking the config file:/home/yunye/work/gwaslab/src/gwaslab/data/config.json
Sat Feb  3 13:46:24 2024  -Config file exists.
Sat Feb  3 13:46:24 2024  -Updating config.json...
Sat Feb  3 13:46:24 2024   - ensembl_hg19_gtf  :  /home/yunye/.gwaslab/Homo_sapiens.GRCh37.87.chr.gtf.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_hg19  :  /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_hg19_tbi  :  /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - recombination_hg38  :  /home/yunye/.gwaslab/recombination/hg38/recombination_hg38.tar.gz
Sat Feb  3 13:46:24 2024   - ensembl_hg38_gtf  :  /home/yunye/.gwaslab/Homo_sapiens.GRCh38.109.chr.gtf.gz
Sat Feb  3 13:46:24 2024   - ucsc_genome_hg19  :  /home/yunye/.gwaslab/hg19.fa
Sat Feb  3 13:46:24 2024   - ucsc_genome_hg38  :  /home/yunye/.gwaslab/hg38.fa
Sat Feb  3 13:46:24 2024   - refseq_hg19_gtf  :  /home/yunye/.gwaslab/GRCh37_latest_genomic.gtf.gz
Sat Feb  3 13:46:24 2024   - refseq_hg38_gtf  :  /home/yunye/.gwaslab/GRCh38_latest_genomic.gtf.gz
Sat Feb  3 13:46:24 2024   - 1kg_dbsnp151_hg19_auto  :  /home/yunye/.gwaslab/1kg_dbsnp151_hg19_auto.txt.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_x_hg19  :  /home/yunye/.gwaslab/EAS.chrX.split_norm_af.1kgp3v5.vcf.gz
Sat Feb  3 13:46:24 2024   - 1kg_eas_x_hg19_tbi  :  /home/yunye/.gwaslab/EAS.chrX.split_norm_af.1kgp3v5.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - 1kg_afr_hg19  :  /home/yunye/.gwaslab/AFR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
Sat Feb  3 13:46:24 2024   - 1kg_afr_hg19_tbi  :  /home/yunye/.gwaslab/AFR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz.tbi
Sat Feb  3 13:46:24 2024   - testlink_tbi  :  /home/yunye/.gwaslab/EAS.chr22.split_norm_af.1kgp3v5.vcf.gz.tbi
```

or you can simply delete files in `~/.gwaslab`

## Get the path of reference data

Get the path using keywords. The path can be passed to other functions.

!!! example
    ```python
    gl.get_path("1kg_eas_hg19")
    ```

```
'/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz'
```

If you haven't downloaded it. It will return False.

!!! example
    ```python
    gl.get_path("1kg_eur_hg19")
    ```

**stdout:**
```
Sat Feb  3 13:46:24 2024 No records in config file. Please download first.
```

```
False
```

## Update available reference list

!!! example
    ```python
    gl.update_available_ref()
    ```

**stdout:**
```
Sat Feb  3 13:46:24 2024 Updating available_ref list from: https://raw.github.com/Cloufield/gwaslab/main/src/gwaslab/data/reference.json
Sat Feb  3 13:46:25 2024 Available_ref list has been updated!
```
