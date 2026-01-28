# Checking if Lead Variants are Novel

GWASLab can check if the lead variants of your summary statistics overlap with reported variants or not based on the physical distance.

## .get_novel()

```python
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
    show_child_traits=True,
    verbose=True
)
```

GWASLab checks overlap with a local file of variants or records in GWAS Catalog.

When **`build="19"`** (hg19/GRCh37), coordinates are first lifted over to hg38, then the same steps run. GWAS catalog queries and novelty checks use hg38; returned coordinates are in hg38 in that case.

## Required Parameters

**Either `known` or `efo` must be provided:**

- **`known`**: `string` or `pandas.DataFrame`, path to the local file of reported variants or a DataFrame containing known variants with **CHR**/**POS** columns

- **`efo`**: `string` or `list`, EFO ID, MONDO ID, or trait name for the target trait(s), used for querying the GWAS Catalog API v2. Supports IDs and trait names; a list may mix formats.
  
  Examples:
  - EFO ID: `"EFO_0001360"` (type 2 diabetes)
  - MONDO ID: `"MONDO_0005148"` (type 2 diabetes)
  - Trait name: `"type 2 diabetes mellitus"`
  - Multiple traits (mixed): `efo=['type 2 diabetes mellitus', 'MONDO_0004247', 'EFO_0004330']`

## Parameters

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `known` | `string` or `DataFrame` | Path to local file of reported variants or DataFrame with **CHR**/**POS** columns | `None` |
| `efo` | `string` or `list` | EFO ID(s), MONDO ID(s), or trait name(s); a list may mix formats, e.g. `['coffee consumption','MONDO_0004247','EFO_0004330']` | `None` |
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
| `build` | `"19"` or `"38"` | Genome build of sumstats. Use "38" (GRCh38/hg38) or "19" (GRCh37/hg19); when "19", coordinates are lifted to hg38 first and results are in hg38 | `"19"` |
| `source` | `"ensembl"` or `"refseq"` | Database source for gene annotation | `"ensembl"` |
| `gwascatalog_source` | `"NCBI"` or `"EBI"` | Source for GWAS catalog data | `"NCBI"` |
| `output_known` | `boolean` | If True, additionally output the reported variants | `False` |
| `show_child_traits` | `boolean` | If True, include child traits in GWAS Catalog results when querying by efo | `True` |
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
    The `efo` parameter supports EFO IDs, MONDO IDs, and trait names. A list may mix any of these:
    - **EFO ID**: `"EFO_0001360"` (e.g., for type 2 diabetes)
    - **MONDO ID**: `"MONDO_0005148"` (e.g., for type 2 diabetes)
    - **Trait name**: `"type 2 diabetes mellitus"` or `"coffee consumption"`
    - **Mixed list for `get_novel`**: `efo=['coffee consumption', 'MONDO_0004247', 'EFO_0004330']`
    
    The function automatically handles MONDO to EFO conversion and trait name lookups. If a MONDO ID is provided, it will attempt to find the corresponding EFO ID. If that fails, it will try using the trait name.

!!! note "Genome Build"
    When using GWAS Catalog (`efo` parameter), specify your sumstats build via `build`. If `build="19"` (hg19/GRCh37), coordinates are automatically lifted to hg38 before querying the catalog and comparing; returned coordinates are in hg38. If `build="38"`, no liftover is applied.

!!! note "GWAS Catalog Trait Associations"
    By default (`show_child_traits=True`), associations include the specified EFO trait(s) and their child traits. Set `show_child_traits=False` to restrict results to the specified trait(s) only.

## Examples

!!! example "Check novelty using GWAS Catalog (EFO ID)"
    ```python
    # Check if lead variants are novel for type 2 diabetes
    novel_variants = mysumstats.get_novel(
        efo="EFO_0001360",
        build="38",
        windowsizekb_for_novel=1000
    )
    ```

!!! example "Check novelty using GWAS Catalog (trait name)"
    ```python
    # Check if lead variants are novel using trait name
    novel_variants = mysumstats.get_novel(
        efo="type 2 diabetes mellitus",
        build="38"
    )
    ```

!!! example "Check novelty using local file"
    ```python
    # Check if lead variants are novel using local file
    novel_variants = mysumstats.get_novel(
        known="/path/to/known_variants.txt",
        windowsizekb_for_novel=1000
    )
    ```

!!! example "Check novelty using DataFrame"
    ```python
    # Check if lead variants are novel using DataFrame
    known_df = pd.read_csv("/path/to/known_variants.txt")
    novel_variants = mysumstats.get_novel(
        known=known_df,
        windowsizekb_for_novel=1000
    )
    ```

!!! example "Get only novel variants"
    ```python
    # Return only novel variants
    novel_only = mysumstats.get_novel(
        efo="EFO_0001360",
        only_novel=True,
        build="38"
    )
    ```

!!! example "Output known variants as well"
    ```python
    # Also output the known variants for comparison
    result = mysumstats.get_novel(
        efo="EFO_0001360",
        output_known=True,
        build="38"
    )
    ```

!!! example "Skip lead variant extraction"
    ```python
    # If you already have lead variants extracted
    novel_variants = mysumstats.get_novel(
        efo="EFO_0001360",
        if_get_lead=False,  # Skip lead variant extraction
        build="38"
    )
    ```

!!! example "Multiple traits (EFO ID, MONDO ID, or trait name)"
    ```python
    # efo supports a mix of trait names, MONDO IDs, and EFO IDs
    novel_variants = mysumstats.get_novel(
        efo=['coffee consumption', 'MONDO_0004247', 'EFO_0004330'],
        build="38"
    )
    ```

!!! example "Custom cache directory"
    ```python
    # Use custom cache directory for GWAS catalog data
    novel_variants = mysumstats.get_novel(
        efo="EFO_0001360",
        cache_dir="/path/to/cache",
        build="38"
    )
    ```

## Caching (when using `efo`)

When you use the `efo` parameter to query the GWAS Catalog, `use_cache` and `cache_dir` control where results are stored and reused:

- **`use_cache`** (default `True`): Use cached data when available.
- **`cache_dir`** (default `"./"`): Directory for cache files.

**Cache file naming:**

- **API v2**: JSON files named  
  `GWASCatalog_v2_{efo}_associationsByTraitSummary_text_{date}_{suffix}.json`  
  (the suffix encodes significance level and `show_child_traits`).
- **Legacy fallback** (e.g. when the v2 API fails or for some MONDO traits): JSON files named  
  `GWASCatalog_{efo}_associationsByTraitSummary_text_{date}_childTrue.json`  
  or `..._childFalse.json` in the same `cache_dir`.

Caching applies only when known variants are fetched via `efo`; it is not used when you pass a local `known` file or DataFrame.

## Notes

- The function first extracts lead variants (unless `if_get_lead=False`) and then checks their novelty.
- When `build="19"`, sumstats are lifted to hg38 before lead extraction and novelty checks; returned coordinates are in hg38.
- Novelty is determined based on physical distance: variants within `windowsizekb_for_novel` kb of a known variant are considered "known".
- GWAS Catalog queries require internet connection and may take some time depending on the trait.
- The function supports both GRCh37 (build="19") and GRCh38 (build="38"); with build="19", liftover to hg38 is applied automatically.
- When using `group_key`, variants are compared within groups, useful for multi-trait analyses.

## Reference

- Buniello, A., MacArthur, J. A. L., Cerezo, M., Harris, L. W., Hayhurst, J., Malangone, C., ... & Parkinson, H. (2019). The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. Nucleic acids research, 47(D1), D1005-D1012.
