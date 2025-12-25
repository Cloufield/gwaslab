# Scatter & Distribution plot : allele frequency comparison

!! info "Available from v3.4.15"


## .check_af()

Calculate the difference in allele frequency (**DAF**) between the effect allele frequency (**EAF**) in your sumstats and the alternative allele frequency (**ALT_AF**) in a reference VCF/BCF file.

**Purpose:**

- Quality control: Identify variants with large differences in allele frequency
- Validation: Verify **EAF** values are consistent with reference population frequencies
- Detect issues: Large |**DAF**| values (> 0.2) may indicate allele mismatches, strand flips, or data quality issues

**DAF Calculation:**

- **DAF** = **EAF** (sumstats) - **ALT_AF** (reference VCF)
- Positive **DAF**: **EAF** in sumstats is higher than reference
- Negative **DAF**: **EAF** in sumstats is lower than reference

```python
# Check the difference between the EAF in the sumstats and the allele frequency in VCF files
mysumstats.check_af(ref_infer="path/to/reference.vcf.gz",
                  ref_alt_freq="AF",
                  maf_threshold=0.40,
                  n_cores=2)
```

**Parameters:**

- `ref_infer`: Path to reference VCF/BCF file (must be tabix-indexed)
- `ref_alt_freq`: Field name for ALT frequency in VCF INFO (e.g., "AF", "AF_popmax", "gnomAD_AF")
- `maf_threshold`: **MAF** threshold for filtering variants (default: 0.40)
- `column_name`: Column name to store **DAF** values (default: "DAF")
- `n_cores`: Number of CPU threads (default: 1)

For more details, see [Harmonization documentation](https://cloufield.github.io/gwaslab/Harmonization/#check-the-difference-in-allele-frequency).

## .plot_daf()

```python
#allele frequnecy correlation plot
mysumstats.plot_daf()
```

You need to run 'check_af()' first before plotting. For check_af(), see [here](https://cloufield.github.io/gwaslab/Harmonization/#check-the-difference-in-allele-frequency).

**Options for `plot_daf`:**
`threshold`: `float`, the threshold used to determine outliers.

## Examples

!!! example
    ```python
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
                 snpid="SNP",
                 chrom="CHR",
                 pos="POS",
                 ea="ALT",
                 nea="REF",
                 neaf="Frq",
                 beta="BETA",
                 se="SE",
                 p="P",
                 direction="Dir",
                 n="N",nrows=10000)
    
    # harmonize
    mysumstats.harmonize(basic_check = True, 
                         ref_seq=gl.get_path("ucsc_genome_hg19"))
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/dc3716d6-1b85-42dd-99ca-d0a38a3ad6a1">

    ```python
    # check the difference in allele frequency with reference vcf
    mysumstats.check_af(ref_infer=gl.get_path("1kg_eas_hg19"), 
                        ref_alt_freq="AF",
                        threads=2)
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/8cdc53b4-f661-40d6-a61b-ab5f23af6a4c">

    ```python
    plot and get the outliers
    outliers = mysumstats.plot_daf(threshold=0.12, 
                                    save="af_correlation.png",
                                    save_args={"dpi":300})
    ```
    <img width=600 src="https://github.com/Cloufield/gwaslab/assets/40289485/0e9c92d4-3bea-4734-8ba6-76f7fc6e329b">

    ```python
    outliers[1]
    ```

