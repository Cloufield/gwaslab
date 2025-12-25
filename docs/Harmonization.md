# Harmonization

GWASLab provides reference-dependent harmonization functions.

## Methods summary

| Sumstats Methods       | Options                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | Description                                                                                                                       |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `.check_ref()`         | `ref_path`,<br /> `chr_dict=get_chr_to_number()`                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | Check alignment with a reference genome FASTA file                                                                                |
| `.assign_rsid()`       | `ref_rsid_tsv`,<br /> `ref_rsid_vcf`,<br /> `threads=1`, <br />`chunksize=5000000`, <br />`chr_dict=get_number_to_chr()`, <br />`overwrite="empty"`                                                                                                                                                                                                                                                                                                                                                                       | Annotate **rsID** using a reference tabular file or VCF/BCF file                                                                      |
| `.infer_strand()`      | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`remove_snp=""`,<br />`daf_tolerance=0.2`, ,<br />`mode="pi"`,<br />`threads=1`,<br />`remove_indel=""`                                                                                                                                                                                                                                                                                                                                            | Infer the strand of palindromic variants/indistinguishable indels using reference VCF/BCF files based on allele frequency in **INFO** |
| `.check_af()`         | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`                                                                                                                                                                                                                                                                                                                                                                                                                                        | Calculate difference in allele frequencies (**DAF**) between sumstats **EAF** and reference VCF **ALT** frequency                                                         |
| `.flip_allele_stats()` |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | After alignment and inferring, flip the alleles and allele-specific statistics to harmonize the variants.                         |
| `.harmonize()`         | `basic_check=True`, <br /> `ref_seq=None`,<br />`ref_rsid_tsv=None`,<br />`ref_rsid_vcf=None`,<br />`ref_infer=None`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`threads=1`,<br />`remove=False`,<br />`sweep_mode=False`,<br />`checkref_args={}`,<br />`removedup_args={}`,<br />`assignrsid_args={}`,<br />`inferstrand_args={}`,<br />`flipallelestats_args={}`,<br />`fixid_args={}`,<br />`fixchr_agrs={}`,<br />`fixpos_args={}`,<br />`fixallele_args={}`,<br />`sanitycheckstats_args={}`,<br />`normalizeallele_args={}` | All-in-one function for harmonization                                                                                             |

## Align NEA with REF in the reference genome

`.check_ref()`:  Check if **NEA** is aligned with the reference sequence. After checking, the tracking status code will be changed accordingly.

| `.check_ref()` options                  | DataType   | Description                                                                               | Default                            |
|-----------------------------------------|------------|-------------------------------------------------------------------------------------------|------------------------------------|
| `ref_path`                              | `string`   | path to the reference genome FASTA file.                                                  |                                    |
| `ref_seq_mode` (available since v3.4.42) | `v` or `s` | `v` for vectorized implementation (faster); `s` for single row iteration mode                     | `v` since v3.4.42 (except v3.4.43) |
| `chr_dict`                              | `dict`     | a conversion dictionary for chromosome notations in reference FASTA and those in sumstats | `gl.get_chr_to_number()`           |

!!! example
    ```python
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    ```
!!! note
    `check_ref()` only change the status code. Use [flip function](#flip-allele-specific-statistics) `.flip_allele_stats()` to flip the allele-specific stats.

## Assign rsID according to CHR, POS, REF/ALT

`.assign_rsid()` : Annotated variants with **rsID** using a reference tsv file (1KG variants) and reference vcf file (tabix indexed, entire dbSNP).

See [https://cloufield.github.io/gwaslab/AssignrsID/](https://cloufield.github.io/gwaslab/AssignrsID/)

!!! example
    ```python
    mysumstats.assign_rsid(ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg19_auto"), 
                           ref_rsid_vcf = "/home/yunye/mydata/d_disk/dbsnp/GCF_000001405.25.vcf.gz",
                           chr_dict = gl.get_number_to_NC(build="19") 
                           threads=1)
    ```

- For TSV file, variants will be matched using **SNPID** (**CHR**:**POS**:**NEA**:**EA**) for quick assigning.
- For VCF file, GWASLab will first extract all variants in the reference file with matching **CHR** and **POS**. And then compare **EA**/**NEA** in sumstats with **REF**/**ALT** in reference vcf. When matching, it will annotate the variant in sumstats with the matching **rsID** in reference vcf.  

## Check palindromic SNPs or indistinguishable Indels

`.infer_strand()` and `.infer_strand2()`:

- Infer the strand for palindromic SNPs (A/T, T/A, G/C, C/G) with **MAF** ≤ `maf_threshold`. 
- Check the alignment status of indels with the **REF** allele in a reference VCF/BCF file and verify if the allele frequencies are consistent (**DAF** ≤ `daf_tolerance`).
- Update **STATUS** code 7th digit to reflect strand inference results (0=non-palindromic, 1=forward, 2=reverse, 3=indel forward, 4=ambiguous indel, 6=indel reverse, 7=**MAF** too high, 8=not found, 9=unchecked).

!!! info "DAF in GWASLab"
    **DAF** in GWASLab: Difference between effect allele frequency (**EAF**) in sumstats and **ALT** frequency in reference VCF/BCF file

!!! note
    This **DAF** in GWASLab is not the derived allele frequency in evolutionary genetics.

| `.infer_strand()` options | DataType          | Description                                                                                                                 | Default                  |
|---------------------------|-------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------------------|
| `ref_infer`               | `string`          | path to the reference VCF/BCF file (must be tabix-indexed for efficient querying).                                                                |                          |
| `ref_alt_freq`            | `string`          | the field for alternative allele frequency in VCF INFO (e.g., "AF", "AF_popmax", "gnomAD_AF"). If None, auto-detection is attempted.                                                                          |                          |
| `chr_dict`                | `dict`            | a conversion dictionary for chromosome notations in sumstats and those in reference VCF/BCF                                 | `gl.get_number_to_chr()` |
| `maf_threshold`           | `float`          | maximum minor allele frequency threshold. Only palindromic SNPs with **MAF** ≤ `maf_threshold` will be inferred                                                          | `0.4`                    |
| `ref_maf_threshold`       | `float`          | maximum minor allele frequency threshold for reference **RAF**. Used as additional filter for both palindromic SNPs and indels. | `0.4`                    |
| `daf_tolerance`           | `float`          | difference in allele frequency tolerance for indels. Only indels with |**EAF** - **RAF**| or |**EAF** - (1-**RAF**)| ≤ `daf_tolerance` will be inferred                       | `0.2`                    |
| `remove_snp`              | ` `, `7` or `8`   | `7` remove palindromic SNPs with **MAF** unable to infer. `8`: remove palindromic SNPs with No information in reference VCF/BCF | ``                       |
| `remove_indel`            | ` ` or `8`        | `8`: indistinguishable indels with No information in reference VCF/BCF                                                      | ``                       |
| `mode`                    | `p`, `i`, or `pi` | `p`: infer palindromic SNPs only. `i`: infer indels only. `pi`: infer both (default).                                                                            | `pi`                     |
| `threads`                 | `int`             | number of CPU threads to use                                                                                                | `1`                      |
| `cache_options` (available since v3.4.42)           | `dict`            | options for using cache to speed up this step                                                                                                     | `None`                   |

| `.infer_strand2()` options | DataType          | Description                                                                                                                 | Default                  |
|----------------------------|-------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------------------|
| `path` / `vcf_path` / `tsv_path` | `string` | path to reference file (VCF/BCF or TSV). `vcf_path` overrides `path`, `tsv_path` overrides both. If TSV not provided, generated from VCF. |                          |
| `assign_cols`              | `tuple` or `str`  | column names to extract from reference (default: ("AF",)). First column will be renamed to RAF.                            | `("AF",)`                |
| `maf_threshold`           | `float`          | maximum minor allele frequency threshold for sumstats EAF                                                          | `0.4`                    |
| `ref_maf_threshold`       | `float`          | maximum minor allele frequency threshold for reference RAF                                 | `0.4`                    |
| `daf_tolerance`           | `float`          | difference in allele frequency tolerance for indels                       | `0.2`                    |
| `threads`                 | `int`             | number of CPU threads for processing                                                                                                | `6`                      |
| `reuse_lookup`            | `bool`            | if True, reuse existing lookup table TSV file if found                                                                     | `True`                   |
| `convert_to_bcf`          | `bool`            | if True, convert VCF to BCF format for faster processing                                                                   | `False`                  |
| `strip_info`              | `bool`            | if True, strip INFO fields from VCF during lookup table generation                                                          | `True`                   |
| `chr_dict`                | `dict`            | chromosome name mapping dictionary                                                                                           | `None`                   |


| `cache_options` | DataType                                  | Description                                                                                                                                                                                                                                             | Default |
|-----------------|-------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `cache_manager` | CacheManager object or `None`               | If any between cache_loader and cache_process is not None, or use_cache is True, a CacheManager object will be created automatically.                                                                                                                   |         |
| `trust_cache`   | `bool`                                      | Whether to completely trust the cache or not. Trusting the cache means that any key not found inside the cache will be considered as a missing value even in the VCF file.                                                                              | True    |
| `cache_loader`  | Object with a get_cache() method or `None`  | Object with a get_cache() method or None.                                                                                                                                                                                                               |         |
| `cache_process` | Object with an apply_fn() method or `None` | Object with an apply_fn() method or None.                                                                                                                                                                                                               |         |
| `use_cache`     | `bool`                                      | If any of the cache_manager, cache_loader or cache_process is not None, this will be set to True automatically. If set to True and all between cache_manager, cache_loader and cache_process are None, the cache will be loaded (or built) on the spot. | False   |

!!! example
    ```python
    # Old method: per-variant VCF queries (good for small datasets)
    mysumstats.infer_strand(ref_infer=gl.get_path("1kg_eas_hg19"),
                            ref_alt_freq="AF",
                            maf_threshold=0.40,
                            daf_tolerance=0.20,
                            threads=4)
    mysumstats.flip_allele_stats()

    # New optimized method: bulk lookup table (faster for large datasets)
    mysumstats.infer_strand2(vcf_path=gl.get_path("1kg_eas_hg19"),
                             maf_threshold=0.40,
                             ref_maf_threshold=0.40,
                             daf_tolerance=0.20,
                             threads=6)
    mysumstats.flip_allele_stats()

    # Using cache to speed up old method
    mysumstats.infer_strand(ref_infer=gl.get_path("1kg_eas_hg19"),
                            cache_options={"use_cache":True})
    mysumstats.flip_allele_stats()
    ```
    
!!! note
    `infer_strand()` and `infer_strand2()` only change the status code. Use [flip function](#flip-allele-specific-statistics) `.flip_allele_stats()` to flip the allele-specific stats.
    
!!! tip "When to use infer_strand2()"
    Use `.infer_strand2()` for:
    - Large datasets (millions of variants)
    - When you have `bcftools` installed
    - When you want faster processing with reduced I/O operations
    - When you plan to reuse lookup tables for multiple runs
    
    Use `.infer_strand()` for:
    - Small to medium datasets (< 1 million variants)
    - When `bcftools` is not available
    - When you need maximum parallelization across CPU cores


## Check the difference in allele frequency

`.check_af()` : Check the allele frequency discrepancy with a reference VCF/BCF file. This function calculates the difference (**DAF**) between the effect allele frequency (**EAF**) in your sumstats and the alternative allele frequency (**ALT_AF**) in the reference VCF.

**Purpose:**

- Quality control: Identify variants with large differences in allele frequency between sumstats and reference
- Validation: Verify that **EAF** values in sumstats are consistent with reference population frequencies
- Detect potential issues: Large |**DAF**| values (> 0.2) may indicate allele mismatches, strand flips, or data quality issues

**Important:** Please make sure your sumstats are already harmonized, and the variants in reference VCF are also aligned. GWASLab will retrieve information only for matched variants (**CHR**, **POS**, **EA**-**ALT**, and **NEA**-**REF**).

| `.check_af()` options | DataType | Description                                                                                 | Default                  |
|------------------------|----------|---------------------------------------------------------------------------------------------|--------------------------|
| `ref_infer`            | `string` | path to the reference VCF/BCF file (must be tabix-indexed for efficient querying).                                |                          |
| `ref_alt_freq`         | `string` | the field for alternative allele frequency in VCF INFO (e.g., "AF", "AF_popmax", "gnomAD_AF"). If None, auto-detection is attempted.                                          |                          |
| `maf_threshold`        | `float`  | minor allele frequency threshold for filtering variants. Only variants with MAF ≤ threshold in both sumstats and reference are included in DAF statistics. | `0.4`                    |
| `column_name`           | `string` | name of the column to store DAF values                                                      | `"DAF"`                  |
| `suffix`               | `string` | suffix to append to column name (useful when comparing multiple reference populations)      | `""`                     |
| `chr_dict`             | `dict`   | a conversion dictionary for chromosome notations in sumstats and those in reference VCF/BCF | `gl.get_number_to_chr()` |
| `n_cores`              | `int`    | number of CPU threads to use                                                                | `1`                      |
| `force`                | `bool`   | if True, check all variants regardless of STATUS codes                                       | `False`                  |

!!! example
    ```python
    mysumstats.check_af(ref_infer=gl.get_path("1kg_eas_hg19"), 
                        ref_alt_freq="AF",
                        n_cores=2)
    ```

**Terminology:**

- **DAF**: Difference in Allele Frequency = **EAF** (sumstats) - **ALT_AF** (reference VCF)
- **EAF**: Effect allele frequency (from your sumstats)
- **RAF**: Reference **ALT** allele frequency (from reference VCF/BCF file)

**Interpretation:**

- Positive **DAF**: **EAF** in sumstats is higher than reference
- Negative **DAF**: **EAF** in sumstats is lower than reference
- Large |**DAF**| values (> 0.2) may indicate:
  - Population differences (expected for population-specific variants)
  - Allele mismatches (check **EA**/**NEA** alignment)
  - Strand flips (use `infer_strand()` to resolve)
  - Data quality issues (verify EAF calculation in sumstats)

You may want to check the allele frequency discrepancy with a reference VCF. Just specify the path and the right allele frequency for your target ancestry in INFO field.

## Allele frequency correlation plot

GWASLab will simply calculate DAF = EAF (sumstats) - frequency in VCF file, and store the results in DAF column. 

DAF can then be used for plotting (`.plot_daf()`) or filter variants.

!!! example
    ```python
    mysumstats.plot_daf(threshold=0.12)
    ```
![image](https://github.com/Cloufield/gwaslab/assets/40289485/0c607470-bbb6-4f11-93fe-038a53f6eebb)

## Flip allele-specific statistics

`.flip_allele_stats()` :  Flip allele-specific statistics to harmonize the variants based on the tracking status code in STATUS. 

!!! example
    ```python 
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    
    mysumstats.infer_strand()
    mysumstats.flip_allele_stats()
    ```

## Sweep Mode for Large Datasets

!!! info "Available since v3.4.42"

Sweep mode is an optimized processing mode designed for large datasets that significantly improves performance when harmonizing summary statistics with reference VCF/BCF files by reducing I/O operations. Sweep mode uses **bcftools** for efficient VCF/BCF processing.

!!! warning "bcftools Required"
    Sweep mode requires **bcftools** to be installed and available in your PATH. If bcftools is not found, GWASLab will raise a `RuntimeError`.

### How Sweep Mode Works

When `sweep_mode=True`, GWASLab uses a different processing strategy that creates lookup tables instead of querying the reference file per-variant:

- **Normal mode (sweep_mode=False)**: 
  - For **strand inference** (`ref_infer`): Uses `_parallelize_infer_strand` which queries the reference VCF/BCF file using tabix for each variant (or chunk of variants). Each query requires a separate I/O operation to the reference file.
  - For **rsID assignment** (`ref_rsid_vcf`): Uses `_parallelize_assign_rsid` which also performs per-variant or per-chunk tabix queries against the reference VCF/BCF file.

- **Sweep mode (sweep_mode=True)**: 
  - For **strand inference** (`ref_infer`): Uses `_infer_strand_with_annotation` which:
    1. First extracts all needed variants from the reference VCF/BCF in **one pass** using **bcftools** (`bcftools view`, `bcftools query`) to create a lookup table (TSV)
    2. Then annotates all variants in your sumstats using this lookup table (fast in-memory operations)
    3. Finally infers strand orientation using the annotated allele frequencies
  - For **rsID assignment** (`ref_rsid_vcf`): Uses `_assign_rsid` which:
    1. First extracts a lookup table from the reference VCF/BCF in **one pass** using **bcftools** (`bcftools view`, `bcftools query`, `bcftools index`) via `_extract_lookup_table_from_vcf_bcf`
    2. Then assigns rsIDs to all variants using `_assign_from_lookup` (fast in-memory operations)

### Key Advantages of Sweep Mode

- **Reduced I/O operations**: Instead of many tabix queries (one per variant/chunk), sweep mode makes a single pass through the reference file to extract all needed data
- **Better performance for large datasets**: The lookup table approach scales better when processing millions of variants
- **Reusable lookup tables**: The extracted lookup tables can be saved and reused in subsequent runs (if `lookup_path` is specified)
- **Lower overhead**: Fewer file system operations and less network I/O (if using remote files)

### When to Use Sweep Mode

Use sweep mode when:
- Processing **large datasets** (millions of variants)
- Using **large reference VCF/BCF files** (e.g., full dbSNP, 1000 Genomes)
- You want to **minimize I/O operations** and improve throughput
- You plan to **reuse lookup tables** for multiple harmonization runs
- Processing variants from **remote/network-mounted** reference files

Use normal mode when:
- Processing **small to medium datasets** (< 1 million variants)
- You need **maximum parallelization** across CPU cores (normal mode parallelizes tabix queries)
- The reference file is **small or already indexed** in a way that makes per-variant queries efficient
- You prefer **lower memory usage** (sweep mode loads lookup tables into memory)

### Usage Example

!!! example
    ```python
    # Harmonization with sweep mode for large dataset
    mysumstats.harmonize(
        ref_seq="/path/to/reference.fasta",
        ref_rsid_vcf="/path/to/dbsnp.vcf.gz",
        ref_infer="/path/to/reference.vcf.gz",
        ref_alt_freq="AF",
        sweep_mode=True,  # Enable sweep mode
        threads=8
    )
    
    ```

!!! note
    - Sweep mode only affects processing with `ref_rsid_vcf` (VCF/BCF files) and `ref_infer` (for strand inference)
    - TSV-based rsID assignment (`ref_rsid_tsv`) always uses the same efficient method regardless of sweep mode
    - Sweep mode uses **bcftools** for VCF/BCF processing - ensure bcftools is installed and in your PATH
    - Sweep mode creates temporary lookup TSV files that can be reused if you specify `lookup_path` in the kwargs
    - The lookup table extraction step may take some time initially, but subsequent annotation is very fast
    - bcftools operations are parallelized using the `threads` parameter

## Assign CHR and POS according to rsID and reference data

```python
mysumstats.rsid_to_chrpos()  

mysumstats.rsid_to_chrpos2()  
```

