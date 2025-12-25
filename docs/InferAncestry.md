# Infer Ancestry

GWASLab infers the genetic ancestry of your summary statistics by comparing effective allele frequencies (EAF) with 1000 Genomes Project reference data using Fst (fixation index) calculations.

## How It Works

The function compares the effective allele frequencies from your sumstats against allele frequencies from 26 populations in the 1000 Genomes Project. It calculates Fst values for each population and identifies the population with the minimum average Fst, which indicates the closest genetic ancestry match.

1. **Loads 1000 Genomes reference data** (Hapmap3 SNPs) for the specified genome build
2. **Matches CHR:POS coordinates** between your sumstats and reference data
3. **Aligns alleles** (EA/NEA) with reference ALT/REF alleles
4. **Calculates Fst values** for each variant against all 26 populations
5. **Identifies closest ancestry** based on minimum average Fst
6. **Updates metadata** in the Sumstats object

## .infer_ancestry()

```
mysumstats.infer_ancestry(build="19", verbose=True)
```

## Parameters

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `build` | `str` | Genome build version. Must be "19" (hg19/GRCh37) or "38" (hg38/GRCh38). Required. | Uses `mysumstats.build` |
| `ancestry_af` | `str` | Path to custom allele frequency file. If None, uses built-in 1kg_hm3_hg19_eaf or 1kg_hm3_hg38_eaf based on build. | `None` |
| `verbose` | `boolean` | If True, print detailed logs including Fst values for all populations | `True` |

!!! note "Required Columns"
    The function requires the following columns in your sumstats:
    - `CHR`: Chromosome
    - `POS`: Base pair position
    - `EA`: Effect allele
    - `NEA`: Non-effect allele
    - `EAF`: Effect allele frequency

!!! note "Build Parameter"
    The `build` parameter is required and must match your sumstats genome build. If not specified, the function will use `mysumstats.build`.

## Return Value

The function returns a string indicating the closest ancestry (e.g., "EUR", "EAS", "AFR") and updates the Sumstats object metadata:
- Updates `mysumstats.meta["gwaslab"]["inferred_ancestry"]` with the inferred ancestry code

## Reference Data

The function uses built-in reference datasets derived from the 1000 Genomes Project:

### Default Reference Files

- **`1kg_hm3_hg19_eaf`** (for build="19"): HapMap3 allele frequency file for hg19/GRCh37
  - File: `PAN.hapmap3.hg19.EAF.tsv.gz`
  - Contains allele frequencies for HapMap3 SNPs from 1000 Genomes Project Phase 3 data
  - Provides population-specific allele frequencies for 26 populations across 5 super-populations
  - Genome build: GRCh37/hg19

- **`1kg_hm3_hg38_eaf`** (for build="38"): HapMap3 allele frequency file for hg38/GRCh38
  - File: `PAN.hapmap3.hg38.EAF.tsv.gz`
  - Contains allele frequencies for HapMap3 SNPs from 1000 Genomes Project 30x data
  - Provides population-specific allele frequencies for 26 populations across 5 super-populations
  - Genome build: GRCh38/hg38

### Reference Data Characteristics

- **SNP Set**: HapMap3 SNPs (~1.2 million variants)
  - High-quality, well-characterized variants commonly used in GWAS
  - Provides good coverage across the genome
  - Standard reference panel for population genetics analyses

- **Source**: 1000 Genomes Project
  - Phase 3 (hg19): 2,504 individuals from 26 populations
  - 30x (hg38): High-coverage sequencing data
  - Publicly available reference dataset for human genetic variation

- **Data Format**: Tab-separated values (TSV) with columns:
  - `CHR`: Chromosome
  - `POS`: Base pair position
  - `REF`: Reference allele
  - `ALT`: Alternative allele
  - Population-specific allele frequency columns (26 populations)

- **Automatic Download**: Reference files are automatically downloaded on first use via `get_path()`
  - Files are cached locally for subsequent use
  - Download location: `~/.gwaslab/` directory

### Custom Reference Files

You can provide a custom allele frequency file using the `ancestry_af` parameter. The custom file should:
- Be in TSV format (tab-separated)
- Contain columns: `CHR`, `POS`, `REF`, `ALT`, and population-specific allele frequency columns
- Match the genome build of your sumstats

## Supported Ancestry Populations

The function compares against 26 populations from the 1000 Genomes Project:

**European (EUR)**: GBR, FIN, IBS, CEU, TSI  
**East Asian (EAS)**: CHS, CDX, CHB, JPT, KHV  
**African (AFR)**: YRI, LWK, GWD, ESN, MSL, ACB, ASW  
**South Asian (SAS)**: GIH, PJL, BEB, STU, ITU  
**American (AMR)**: MXL, PUR, CLM, PEL  

## Examples

!!! example "Basic usage"
    ```
    # Infer ancestry (build will be taken from mysumstats.build)
    mysumstats.infer_ancestry()
    
    # Check the inferred ancestry
    print(mysumstats.meta["gwaslab"]["inferred_ancestry"])  # e.g., "EUR"
    ```

!!! example "Specify build explicitly"
    ```
    # Infer ancestry for hg19 data
    mysumstats.infer_ancestry(build="19")
    
    # Infer ancestry for hg38 data
    mysumstats.infer_ancestry(build="38")
    ```

!!! example "Quiet mode"
    ```
    # Infer ancestry without verbose output
    mysumstats.infer_ancestry(build="19", verbose=False)
    ```

!!! example "Custom reference file"
    ```
    # Use a custom allele frequency file
    mysumstats.infer_ancestry(build="19", ancestry_af="/path/to/custom_af.tsv")
    ```

## Output

The function provides detailed logging output showing Fst values for all populations:

```
Start to infer ancestry based on Fst ...(version)
 -Estimating Fst using 12345 variants...
 -FST_GBR : 0.001234
 -FST_FIN : 0.001456
 -FST_CHS : 0.002345
 ...
 -FST_EUR : 0.001234
 -Closest Ancestry: EUR
Finished inferring ancestry.
```

The populations are sorted by Fst value (lowest to highest), with the closest match (lowest Fst) indicating the inferred ancestry.

## Notes

- **Fst Interpretation**: Lower Fst values indicate greater genetic similarity. The population with the minimum average Fst is identified as the closest ancestry match.
- **EAF Accuracy**: The accuracy of ancestry inference depends on the accuracy of EAF values in your sumstats. Inconsistent or mislabeled EAF values may lead to incorrect ancestry inference.
- **Reference Data**: Uses HapMap3 SNPs from the 1000 Genomes Project (see [Reference Data](#reference-data) section above for details). The reference files are automatically downloaded on first use and cached locally.
- **Allele Alignment**: The function automatically handles allele flipping to ensure proper comparison between sumstats and reference data.
- **Metadata Storage**: The inferred ancestry is stored in the Sumstats object metadata for future reference.
- **File Format**: The reference files are compressed TSV files (`.tsv.gz`) that are automatically decompressed during loading.

## When to Use

- When the genetic ancestry of your sumstats is unknown
- To verify the reported ancestry of your dataset
- As part of quality control workflows to detect potential data issues
- Before performing population-specific analyses
- To ensure proper matching with reference datasets

## Related Functions

- `.infer_build()`: Infer the genome build of your sumstats
- `.check_af()`: Compare allele frequencies with reference data
- `.harmonize()`: Harmonize alleles with reference data (may be needed before ancestry inference)

