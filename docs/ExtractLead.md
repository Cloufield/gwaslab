# Extract Lead Variants

GWASLab can extract the lead variants based on MLOG10P values or P values (by default, GWASLab will use MLOG10P first since v3.4.37) from identified significant loci using a sliding window, and return the result as a pandas.DataFrame or gl.Sumstats Object.

## .get_lead()

```
mysumstats.get_lead(
    scaled=False,
    use_p=False,
    windowsizekb=500,
    sig_level=5e-8,
    anno=False,
    build="19",
    source="ensembl",
    gtf_path=None,
    wc_correction=False,
    verbose=True,
    gls=False
)
```

## Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `scaled` | `boolean` | (deprecated since v3.4.37) If True, use MLOG10P for extraction instead of P values | `False` |
| `use_p` | `boolean` | (available since v3.4.37) If True, use P for extraction instead of MLOG10P values | `False` |
| `windowsizekb` | `int` | Specify the sliding window size in **kb** | `500` |
| `sig_level` | `float` | Specify the P value threshold | `5e-8` |
| `anno` | `boolean` | If True, annotate the lead variants with nearest gene names | `False` |
| `source` | `"ensembl"` or `"refseq"` | When `anno=True`, annotate variants using GTF files from `ensembl` or `refseq` | `"ensembl"` |
| `gtf_path` | `string` | Path to custom GTF file for annotation. If provided, overrides `source` | `None` |
| `build` | `"19"` or `"38"` | Genome build version "19" (GRCh37/hg19) or "38" (GRCh38/hg38) | `"19"` |
| `wc_correction` | `boolean` | If True, apply Winner's Curse correction to effect sizes | `False` |
| `verbose` | `boolean` | If True, print logs | `True` |
| `gls` | `boolean` | If True, return a new gl.Sumstats Object instead of pandas DataFrame | `False` |

!!! note "Lead Variant Extraction"
    `.get_lead()` simply extracts the lead variants of each significant loci. It is different from clumping.

!!! note "Prerequisites"
    Please try running `.basic_check()` to standardize the sumstats if there are any errors.

!!! quote "GBMI Definition"
    GWASLab adopted the definition for novel loci from Global Biobank Meta-analysis Initiative flagship paper. 
    
    **"We defined genome-wide significant loci by iteratively spanning the ±500 kb region around the most significant variant and merging overlapping regions until no genome-wide significant variants were detected within ±1 Mb."** 
    
    (Details are described in Zhou, W., Kanai, M., Wu, K. H. H., Rasheed, H., Tsuo, K., Hirbo, J. B., ... & Study, C. O. H. (2022). Global Biobank Meta-analysis Initiative: Powering genetic discovery across human disease. Cell Genomics, 2(10), 100192.)
    
    GWASLab currently iteratively extends ± `windowsizekb` kb region around the most significant variant and merges overlapping regions until no genome-wide significant variants were detected within ± `windowsizekb`. (slightly different from the GBMI paper. When `windowsizekb=1000`, it is equivalent to GBMI's definition.)

## Examples

!!! example "Basic lead variant extraction"
    ```
    # Extract lead variants
    lead_variants = mysumstats.get_lead(
        windowsizekb=500,
        sig_level=5e-8
    )
    ```

!!! example "Extract lead variants with gene annotation"
    ```
    # Extract lead variants with gene annotation
    lead_variants = mysumstats.get_lead(
        windowsizekb=500,
        sig_level=5e-8,
        anno=True,
        build="38",
        source="ensembl"
    )
    ```

!!! example "Extract lead variants with Winner's Curse correction"
    ```
    # Extract lead variants with Winner's Curse correction
    lead_variants = mysumstats.get_lead(
        windowsizekb=500,
        sig_level=5e-8,
        wc_correction=True
    )
    ```

!!! example "Return as Sumstats object"
    ```
    # Return as Sumstats object instead of DataFrame
    lead_sumstats = mysumstats.get_lead(
        windowsizekb=500,
        sig_level=5e-8,
        gls=True
    )
    ```

!!! example "Using custom GTF file"
    ```
    # Use custom GTF file for annotation
    lead_variants = mysumstats.get_lead(
        windowsizekb=500,
        sig_level=5e-8,
        anno=True,
        gtf_path="/path/to/custom.gtf.gz"
    )
    ```

## Notes

- By default, GWASLab uses MLOG10P for extraction (since v3.4.37). Use `use_p=True` to use P values instead.
- The `scaled` parameter is deprecated but still supported for backward compatibility.
- Gene annotation requires internet connection to download GTF files (unless `gtf_path` is provided).
- Winner's Curse correction adjusts effect sizes for variants that are significant due to sampling variability.
