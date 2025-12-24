# GWAS Catalog API

GWASLab provides integration with the GWAS Catalog REST API v2, allowing you to query and retrieve association data, studies, variants, and traits directly from the GWAS Catalog.

## Overview

The GWAS Catalog API integration in GWASLab provides two main ways to access association data:

1. **Direct API queries**: Use `gl.get_associations()` to query the GWAS Catalog API directly
2. **Sumstats integration**: Use `mysumstats.get_associations()` to fetch associations for variants in your summary statistics

## Direct API Access

### `gl.get_associations()`

Query the GWAS Catalog API v2 directly to retrieve associations.

```
import gwaslab as gl

# Get associations for a specific variant
associations = gl.get_associations(rs_id="rs1050316")

# Get associations for a specific trait
associations = gl.get_associations(efo_trait="type 2 diabetes mellitus")

# Get associations for a specific study
associations = gl.get_associations(accession_id="GCST000001")
```

#### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `efo_trait` | `str`, optional | EFO trait name to filter associations (e.g., "type 2 diabetes mellitus") | `None` |
| `rs_id` | `str`, optional | Variant rsID to filter associations (e.g., "rs1050316") | `None` |
| `accession_id` | `str`, optional | Study accession ID to filter associations (e.g., "GCST000001") | `None` |
| `sort` | `str`, optional | Field to sort by (e.g., "p_value", "or_value") | `None` |
| `direction` | `str`, optional | Sort direction: "asc" or "desc" | `None` |
| `show_child_traits` | `bool`, optional | If True, include child traits in results. If False, only return data annotated directly with the query term | `None` |
| `extended_geneset` | `bool`, optional | If True, use extended gene set (all Ensembl/RefSeq genes within 50kb). If False, use standard gene set | `False` |
| `page` | `int`, optional | Page number for pagination | `None` |
| `size` | `int`, optional | Page size for pagination | `None` |
| `get_all` | `bool` | If True, retrieve all pages automatically | `False` |
| `verbose` | `bool` | If True, print log messages | `True` |

#### Returns

- `pandas.DataFrame`: If multiple associations are returned
- `dict` or `list`: If a single association or raw API response is returned

#### Examples

!!! example "Get associations for a specific variant"
    ```
    import gwaslab as gl
    
    # Get all associations for rs1050316
    associations = gl.get_associations(rs_id="rs1050316", get_all=True)
    print(associations.head())
    ```

!!! example "Get associations for a trait"
    ```
    # Get associations for type 2 diabetes
    associations = gl.get_associations(
        efo_trait="type 2 diabetes mellitus",
        size=10,
        sort="p_value",
        direction="asc"
    )
    ```

!!! example "Get all associations for a study"
    ```
    # Get all associations from a specific study
    associations = gl.get_associations(
        accession_id="GCST000001",
        get_all=True
    )
    ```

### Using GWASCatalogClient

For more advanced usage, you can use the `GWASCatalogClient` class directly:

```
from gwaslab.extension.gwascatalog import GWASCatalogClient

# Create a client
client = GWASCatalogClient(verbose=True)

# Get associations
associations = client.get_associations(rs_id="rs1050316", get_all=True)

# Get studies
studies = client.get_studies(efo_trait="EFO_0001360")

# Get variants
variants = client.get_variants(rs_id="rs1050316")
```

#### `get_known_variants_for_trait()`

A specialized method for retrieving known variants from GWAS Catalog for use with `get_novel()`. This method handles MONDO to EFO conversion, extracts CHR/POS from locations, and returns data in the format expected by `get_novel()`.

```
# Get known variants for a trait (supports EFO ID, MONDO ID, or trait name)
known_variants = client.get_known_variants_for_trait(
    efo="EFO_0001360",  # or "MONDO_0005148" or "type 2 diabetes mellitus"
    sig_level=5e-8,
    verbose=True
)
```

**Parameters:**

- `efo` (str): EFO trait ID (e.g., "EFO_0001360"), MONDO ID (e.g., "MONDO_0005148"), or trait name (e.g., "type 2 diabetes mellitus")
- `sig_level` (float): P-value threshold for filtering associations (default: 5e-8)
- `verbose` (bool): Whether to print log messages (default: True)

**Returns:**

- `pandas.DataFrame`: DataFrame with columns: SNPID, CHR, POS, P, BETA, SE, OR, TRAIT, REPORT_GENENAME, PUBMEDID, AUTHOR, STUDY

**Features:**
- Automatically converts MONDO IDs to EFO IDs when possible
- Falls back to trait names if EFO ID is not available
- Extracts CHR and POS from the `locations` field in API responses
- Handles pagination automatically with reduced logging verbosity
- Falls back to old API v1 for MONDO IDs if API v2 doesn't return results

## Sumstats Integration

### `mysumstats.get_associations()`

Extract GWAS Catalog associations for variants in your summary statistics. This method queries the GWAS Catalog API for each unique variant (rsID) in your sumstats and returns a summary of associations.

!!! warning "Input Limit"
    The input is limited to **100 unique variants**. If your sumstats contains more than 100 unique variants, only the first 100 will be processed and a warning will be issued. This limit prevents excessive API calls and ensures reasonable processing times.

```
import gwaslab as gl

# Load your sumstats
mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")

# Get associations for variants in your sumstats
associations_summary = mysumstats.get_associations()
```

#### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `rsid` | `str` | Name of the rsID column in sumstats | `"rsID"` |
| `fetch_metadata` | `bool` | If True, fetch additional metadata (traits, studies, variants). If False, only fetch associations (faster) | `False` |
| `fetch_traits` | `bool`, optional | If True, fetch traits. If None, uses `fetch_metadata` value | `None` |
| `fetch_studies` | `bool`, optional | If True, fetch studies. If None, uses `fetch_metadata` value | `None` |
| `fetch_variants` | `bool`, optional | If True, fetch variants. If None, uses `fetch_metadata` value | `None` |
| `use_gcv2_format` | `bool` | If True, transform associations to GWASLab GCV2 format and merge with sumstats | `True` |
| `gcv2_columns` | `list`, optional | List of columns to include in GCV2 format output | `None` |
| `verbose` | `bool` | If True, print log messages | `True` |

#### Returns

- `pandas.DataFrame`: Summary of associations with key information (association ID, trait, p-value, beta, etc.)

#### Notes

- **Input limit**: The method processes a maximum of 100 unique variants. If more variants are provided, only the first 100 will be processed and a warning will be issued.
- The full associations are stored in `mysumstats.associations` after calling this method
- If `use_gcv2_format=True`, the associations are transformed to GWASLab's GCV2 format and can be merged with your sumstats
- The method automatically handles API rate limiting and pagination
- API responses are cached to avoid redundant requests

#### Examples

!!! example "Basic usage - get associations only"
    ```
    import gwaslab as gl
    
    mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
    
    # Get associations (fast, no metadata)
    associations = mysumstats.get_associations(fetch_metadata=False)
    
    # View summary
    print(associations.head())
    
    # Access full associations
    print(mysumstats.associations.head())
    ```

!!! example "Get associations with metadata"
    ```
    # Get associations with traits, studies, and variants
    associations = mysumstats.get_associations(
        fetch_metadata=True,
        fetch_traits=True,
        fetch_studies=True,
        fetch_variants=True
    )
    ```

!!! example "Use GCV2 format (default)"
    ```
    # Get associations in GCV2 format (merged with sumstats)
    associations = mysumstats.get_associations(
        use_gcv2_format=True,
        fetch_metadata=False
    )
    
    # The associations are automatically aligned with your sumstats
    # Beta values are aligned to match your effect alleles
    ```

!!! example "Plot associations"
    ```
    # After getting associations, you can plot them
    mysumstats.get_associations()
    mysumstats.plot_associations()
    ```

## Output Format

### Association DataFrame Columns

The association DataFrame typically contains the following columns:

- `associationId`: Unique association ID from GWAS Catalog
- `rsID`: Variant rsID
- `trait` or `GWASCATALOG_TRAIT`: Trait name(s)
- `shortForm`: EFO trait short form
- `pvalue` or `P-value`: P-value
- `betaNum` or `Beta`: Beta value (numeric)
- `betaUnit` or `Unit`: Beta unit
- `riskFrequency` or `RAF`: Risk allele frequency
- `snp_effect_allele`: Effect allele
- `mapped_genes`: Mapped gene names
- Study and variant metadata (if `fetch_metadata=True`)

### GCV2 Format

When `use_gcv2_format=True`, the associations are transformed to GWASLab's GCV2 format, which includes:

- Standardized column names
- Beta values aligned with your sumstats effect alleles
- Merged with your sumstats DataFrame
- Additional GWAS Catalog metadata columns

## Performance Tips

1. **Input limit**: `mysumstats.get_associations()` is limited to 100 unique variants. For larger datasets, filter your sumstats to the most important variants first (e.g., significant variants)
2. **Use `fetch_metadata=False`** for faster queries when you only need basic association information
3. **Use `get_all=True`** when querying the API directly to retrieve all pages automatically
4. **API responses are cached** - repeated queries for the same variants will use cached data
5. **For large datasets**, consider processing variants in batches to avoid API rate limits
6. **Reduced logging**: When processing multiple pages, detailed request logging is suppressed for pages 2 onwards to reduce log verbosity. Progress is shown every 20 pages.

## Related Functions

- `gl.get_studies()`: Query GWAS Catalog studies
- `gl.get_variants()`: Query GWAS Catalog variants
- `mysumstats.plot_associations()`: Plot associations after retrieval
- `mysumstats.get_novel()`: Extract novel loci using GWAS Catalog data (uses `get_known_variants_for_trait()` internally)
- `GWASCatalogClient.get_known_variants_for_trait()`: Retrieve known variants for a trait (supports EFO IDs, MONDO IDs, and trait names)

## Technical Details

### API v2 Integration

GWASLab uses the GWAS Catalog REST API v2, which provides:
- Improved data structure and metadata
- Better support for pagination
- Enhanced trait and variant information

### CHR/POS Extraction

The API v2 returns genomic coordinates in a `locations` field (e.g., `['12:111803962']`). GWASLab automatically extracts chromosome and position from this field, handling various formats:
- Standard format: `"12:111803962"`
- With chr prefix: `"chr12:111803962"`
- Position ranges: `"12:111803962-111803963"` (uses first position)

### MONDO ID Support

The API v2 primarily uses EFO IDs for traits. GWASLab automatically handles MONDO IDs by:
1. Looking up the trait information to find the corresponding EFO ID
2. Using the trait name as a fallback if EFO ID is not available
3. Falling back to the old API v1 if API v2 doesn't return results for MONDO IDs

### Logging and Performance

- **Reduced verbosity**: When processing multiple pages, detailed request logging is suppressed for pages 2 onwards. Progress is shown every 20 pages to reduce log noise.
- **Automatic pagination**: The `get_all=True` parameter automatically retrieves all pages of results.
- **Rate limiting**: Built-in rate limiting respects the GWAS Catalog API rate limits (15 requests per second).

## See Also

- [GWAS Catalog API Documentation](https://www.ebi.ac.uk/gwas/rest/docs/api)
- [Extract Novel Loci](https://cloufield.github.io/gwaslab/ExtractNovel/)

