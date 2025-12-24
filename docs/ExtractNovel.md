# Checking if Lead Variants are Novel

GWASLab can check if the lead variants of your summary statistics overlap with reported variants or not based on the physical distance.

## .get_novel()

```
mysumstats.get_novel(
    known=None,
    efo=None,
    only_novel=False,
    windowsizekb_for_novel=1000,
    windowsizekb=500,
    sig_level=5e-8,
    if_get_lead=True,
    group_key=None,
    use_p=False,
    anno=False,
    wc_correction=False,
    use_cache=True,
    cache_dir="./",
    build="19",
    source="ensembl",
    gwascatalog_source="NCBI",
    output_known=False,
    verbose=True
)
```

GWASLab checks overlap with a local file of variants or records in GWAS Catalog.

## Required Parameters

**Either `known` or `efo` must be provided:**

- **`known`**: `string` or `pandas.DataFrame`, path to the local file of reported variants or a DataFrame containing known variants with **CHR**/**POS** columns

- **`efo`**: `string` or `list`, EFO ID, MONDO ID, or trait name for the target trait, which is used for querying the GWAS Catalog API v2.
  
  Examples:
  - EFO ID: `"EFO_0001360"` (type 2 diabetes)
  - MONDO ID: `"MONDO_0005148"` (type 2 diabetes)
  - Trait name: `"type 2 diabetes mellitus"`
  - Multiple traits: `["EFO_0001360", "EFO_0001361"]`

## Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `known` | `string` or `DataFrame` | Path to local file of reported variants or DataFrame with **CHR**/**POS** columns | `None` |
| `efo` | `string` or `list` | EFO ID, MONDO ID, or trait name(s) for querying GWAS Catalog | `None` |
| `only_novel` | `boolean` | If True, output only novel variants | `False` |
| `windowsizekb_for_novel` | `int` | Window size (kb) for determining if lead variants overlap with reported variants in GWAS Catalog | `1000` |
| `windowsizekb` | `int` | Window size (kb) for lead variant extraction | `500` |
| `sig_level` | `float` | P value threshold for lead variant extraction | `5e-8` |
| `if_get_lead` | `boolean` | If True, first extract lead variants using `.get_lead()` | `True` |
| `group_key` | `string` | Column name for grouping variants (e.g., trait/phenotype ID) | `None` |
| `use_p` | `boolean` | If True, use P values instead of MLOG10P for extraction | `False` |
| `anno` | `boolean` | If True, annotate variants with gene information | `False` |
| `wc_correction` | `boolean` | If True, apply Winner's Curse correction to effect sizes | `False` |
| `use_cache` | `boolean` | If True, use cached GWAS catalog data | `True` |
| `cache_dir` | `string` | Directory for caching downloaded GWAS catalog data | `"./"` |
| `build` | `"19"` or `"38"` | Genome build version "19" (GRCh37/hg19) or "38" (GRCh38/hg38) | `"19"` |
| `source` | `"ensembl"` or `"refseq"` | Database source for gene annotation | `"ensembl"` |
| `gwascatalog_source` | `"NCBI"` or `"EBI"` | Source for GWAS catalog data | `"NCBI"` |
| `output_known` | `boolean` | If True, additionally output the reported variants | `False` |
| `verbose` | `boolean` | If True, print logs | `True` |

## Return Value

Returns a pandas.DataFrame containing variants with novelty status and metadata:

- **`NOVEL`**: Boolean indicating novelty status
- **`DISTANCE_TO_KNOWN`**: Distance to nearest known variant (in base pairs)
- **`LOCATION_OF_KNOWN`**: Relative position to known variant
- **`KNOWN_ID`**: ID of matching known variant
- Additional metadata from known variants (if available)

!!! info "EFO ID"
    You can find the EFO ID by simply searching in GWAS Catalog. For example, the EFO ID for T2D can be obtained:
    
    <img width="700" alt="Screenshot 2023-02-03 at 13 36 16" src="https://user-images.githubusercontent.com/40289485/216513724-27f4e742-a03c-4346-a6bc-6b1002e22847.png">

!!! info "EFO ID, MONDO ID, or Trait Name"
    The `efo` parameter now supports multiple formats:
    - **EFO ID**: `"EFO_0001360"` (e.g., for type 2 diabetes)
    - **MONDO ID**: `"MONDO_0005148"` (e.g., for type 2 diabetes)
    - **Trait name**: `"type 2 diabetes mellitus"` (e.g., for type 2 diabetes)
    - **Multiple traits**: `["EFO_0001360", "EFO_0001361"]` (list of trait identifiers)
    
    The function automatically handles MONDO to EFO conversion and trait name lookups. If a MONDO ID is provided, it will attempt to find the corresponding EFO ID. If that fails, it will try using the trait name.

!!! note "Genome Build"
    When using GWAS Catalog (`efo` parameter), ensure your sumstats are based on the correct genome build. GWAS Catalog data is available for both GRCh37 (build="19") and GRCh38 (build="38"). Make sure to specify the correct `build` parameter.

!!! note "GWAS Catalog Trait Associations"
    Only associations with the specified EFO trait will be obtained. This does not include associations with child traits.

## Examples

!!! example "Check novelty using GWAS Catalog (EFO ID)"
    ```
    # Check if lead variants are novel for type 2 diabetes
    novel_variants = mysumstats.get_novel(
        efo="EFO_0001360",
        build="38",
        windowsizekb_for_novel=1000
    )
    ```

!!! example "Check novelty using GWAS Catalog (trait name)"
    ```
    # Check if lead variants are novel using trait name
    novel_variants = mysumstats.get_novel(
        efo="type 2 diabetes mellitus",
        build="38"
    )
    ```

!!! example "Check novelty using local file"
    ```
    # Check if lead variants are novel using local file
    novel_variants = mysumstats.get_novel(
        known="/path/to/known_variants.txt",
        windowsizekb_for_novel=1000
    )
    ```

!!! example "Check novelty using DataFrame"
    ```
    # Check if lead variants are novel using DataFrame
    known_df = pd.read_csv("/path/to/known_variants.txt")
    novel_variants = mysumstats.get_novel(
        known=known_df,
        windowsizekb_for_novel=1000
    )
    ```

!!! example "Get only novel variants"
    ```
    # Return only novel variants
    novel_only = mysumstats.get_novel(
        efo="EFO_0001360",
        only_novel=True,
        build="38"
    )
    ```

!!! example "Output known variants as well"
    ```
    # Also output the known variants for comparison
    result = mysumstats.get_novel(
        efo="EFO_0001360",
        output_known=True,
        build="38"
    )
    ```

!!! example "Skip lead variant extraction"
    ```
    # If you already have lead variants extracted
    novel_variants = mysumstats.get_novel(
        efo="EFO_0001360",
        if_get_lead=False,  # Skip lead variant extraction
        build="38"
    )
    ```

!!! example "Multiple traits"
    ```
    # Check novelty against multiple traits
    novel_variants = mysumstats.get_novel(
        efo=["EFO_0001360", "EFO_0001361"],  # Multiple traits
        build="38"
    )
    ```

!!! example "Custom cache directory"
    ```
    # Use custom cache directory for GWAS catalog data
    novel_variants = mysumstats.get_novel(
        efo="EFO_0001360",
        cache_dir="/path/to/cache",
        build="38"
    )
    ```

## Notes

- The function first extracts lead variants (unless `if_get_lead=False`) and then checks their novelty.
- Novelty is determined based on physical distance: variants within `windowsizekb_for_novel` kb of a known variant are considered "known".
- GWAS Catalog queries require internet connection and may take some time depending on the trait.
- Cached GWAS Catalog data is stored in `cache_dir` and reused in subsequent runs (if `use_cache=True`).
- The function supports both GRCh37 (build="19") and GRCh38 (build="38") genome builds.
- When using `group_key`, variants are compared within groups, useful for multi-trait analyses.

## Reference

- Buniello, A., MacArthur, J. A. L., Cerezo, M., Harris, L. W., Hayhurst, J., Malangone, C., ... & Parkinson, H. (2019). The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. Nucleic acids research, 47(D1), D1005-D1012.
