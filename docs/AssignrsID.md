# Assigning rsID

GWASLab uses a two-step strategy (both steps are optional) to assign rsIDs to variants in your summary statistics.

- **Step 1 (TSV annotation)**: For quick annotation, GWASLab iterates over a **SNPID**-**rsID** table and assigns **rsID** by joining on **SNPID** (**CHR**:**POS**:**REF**:**ALT**) with sumstats. GWASLab provides curated tables (1KG autosome variants).

- **Step 2 (VCF annotation)**: For full annotation, GWASLab will query a large reference VCF file (dbSNP for example, >20GB) by **CHR**, **POS**, **NEA**, **EA**. It will assign the ID in VCF file to sumstats if the **CHR**, **POS** and **EA**/**NEA** match.

!!! info "New in v4.0.0"
    The rsID assignment process has been optimized for better performance and includes improved error handling. The function now supports better STATUS code filtering to ensure only properly standardized variants receive rsID assignments.

## Quick Start

```
# 1. Download reference data (one-time setup)
gl.download_ref("1kg_dbsnp151_hg19_auto")

# 2. Load and prepare your sumstats
mysumstats = gl.Sumstats("your_sumstats.txt.gz", ...)
mysumstats.basic_check()  # Always run this first!

# 3. Assign rsIDs (simple case)
mysumstats.assign_rsid(
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto")
)
```

## Reference Data

Before assigning rsIDs, you need reference data. GWASLab supports two types:

### SNPID-rsID Table (Recommended for Common Variants)

GWASLab provides curated tables containing ~80M 1KG variants that can be downloaded automatically:

- **hg19 (GRCh37)**: `gl.download_ref("1kg_dbsnp151_hg19_auto")`
- **hg38 (GRCh38)**: `gl.download_ref("1kg_dbsnp151_hg38_auto")`

!!! note "1kg_dbsnp151_hg19_auto format"
    ```
    ~/.gwaslab$ zcat 1kg_dbsnp151_hg19_auto.txt.gz |head
    **SNPID**   **rsID**    **CHR**     **POS**     **NEA**     **EA**
    1:10177:A:AC    rs367896724     1       10177   A       AC
    1:10235:T:TA    rs540431307     1       10235   T       TA
    1:10352:T:TA    rs555500075     1       10352   T       TA
    1:10505:A:T     rs548419688     1       10505   A       T
    1:10511:G:A     rs534229142     1       10511   G       A
    1:10539:C:A     rs537182016     1       10539   C       A
    1:10542:C:T     rs572818783     1       10542   C       T
    1:10579:C:A     rs538322974     1       10579   C       A
    1:10616:CCGCCGTTGCAAAGGCGCGCCG:C        rs376342519     1       10616   CCGCCGTTGCAAAGGCGCGCCG  C
    ```

### VCF/BCF Files (For Rare Variants)

For comprehensive annotation including rare variants, you can download dbSNP VCF files:

**hg19 (GRCh37)** - As of 20240205:
- VCF: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
- Index (tbi): https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

**hg38 (GRCh38)** - As of 20240205:
- VCF: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz 
- Index (tbi): https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

!!! note "VCF file from dbSNP"
    ```
     zcat GCF_000001405.25.vcf.gz | head -100 | tail -10
    NC_000001.10    10059   rs1570391745    C       G       .       .       RS=1570391745;dbSNPBuildID=154;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9997,0.0003425|dbGaP_PopFreq:1,0
    NC_000001.10    10060   rs1639544146    C       CT      .       .       RS=1639544146;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10060   rs1639544159    CT      C       .       .       RS=1639544159;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=DEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10063   rs1010989343    A       C,G     .       .       RS=1010989343;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9928,0.004112,0.003084|Siberian:0.5,0.5,.|dbGaP_PopFreq:1,.,0
    NC_000001.10    10067   rs1489251879    T       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC      .       .       RS=1489251879;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=GnomAD:1,1.789e-05
    NC_000001.10    10067   rs1639545042    T       C       .       .       RS=1639545042;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10067   rs1639545104    TA      T       .       .       RS=1639545104;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10068   rs1639545079    A       T       .       .       RS=1639545079;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10069   rs1570391755    A       C,G     .       .       RS=1570391755;dbSNPBuildID=154;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9966,.,0.003425|dbGaP_PopFreq:1,0,0
    NC_000001.10    10069   rs1639545200    A       AC      .       .       RS=1639545200;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    ```

## Methods Overview

GWASLab provides two methods for assigning rsIDs:

| Method | Best For | Reference Files | Processing Mode |
|--------|---------|----------------|-----------------|
| **`.assign_rsid()`** | Small to medium datasets (< 1M variants) | TSV and/or VCF/BCF | Per-variant/chunk lookup |
| **`.assign_rsid2()`** | Large datasets (millions of variants) | VCF/BCF files | Sweep mode (one-pass extraction) |

!!! info "Sweep Mode vs Standard Mode"
    **Sweep Mode (`.assign_rsid2()`)**: 
    - Extracts all needed variants from VCF/BCF in **one pass** using bcftools
    - Creates a lookup table (TSV) that can be reused
    - Faster for **large datasets** (millions of variants)
    - Better for **large reference VCF/BCF files** (e.g., full dbSNP)
    - Reduces I/O operations significantly
    - Requires **bcftools** to be installed
    
    **Standard Mode (`.assign_rsid()`)**: 
    - Uses per-variant or per-chunk tabix queries
    - Good for **small to medium datasets** (< 1 million variants)
    - Better parallelization across CPU cores
    - Lower memory usage
    - Works with both TSV and VCF/BCF files

!!! warning "Prerequisites"
    - Always run `.basic_check()` first to ensure proper standardization and normalization
    - For sweep mode (`.assign_rsid2()`), **bcftools** must be installed and in your PATH
    - For VCF files, tabix/csi indexing is recommended for optimal performance

## Method 1: .assign_rsid() - Standard Mode

### Basic Usage

```
mysumstats.basic_check()
mysumstats.assign_rsid( 
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    ref_rsid_vcf="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    chr_dict=gl.get_number_to_NC(build="19"),
    threads=2
)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `ref_rsid_tsv` | `string` | TSV file path for annotation of commonly used variants using SNPID (like 1:725932:G:A) as key. This is the first step and is faster for common variants. | `None` |
| `ref_rsid_vcf` | `string` | VCF/BCF file path for annotation of variants with rsID not assigned. .tbi/.csi file is required for indexed VCF files. This is the second step for rare variants not in the TSV file. | `None` |
| `chr_dict` | `dict` | A dictionary for converting chromosome numbers (1-25) to chromosome names in the VCF files. For example, the notation in dbSNP vcf file is based on RefSeq (like NC_000001.10). `gwaslab` provides built-in conversion dictionaries: `gl.get_number_to_NC(build="19")` for hg19 and `gl.get_number_to_NC(build="38")` for hg38. | `None` |
| `threads` | `int` | Number of threads to use for parallel processing. More threads can speed up processing for large VCF files. | `1` |
| `overwrite` | `string` | Overwrite mode for rsID assignment. Options: `"all"`, `"invalid"`, or `"empty"`. See [Overwrite Modes](#overwrite-modes) for details. | `"empty"` |
| `chunksize` | `int` | Size of chunks for processing large reference TSV files. Larger chunks use more memory but may be faster. | `5000000` |

## Method 2: .assign_rsid2() - Sweep Mode

### Basic Usage

```
mysumstats.basic_check()
mysumstats.assign_rsid2(
    vcf_path="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    chr_dict=gl.get_number_to_NC(build="19"),
    threads=6,
    overwrite="empty"
)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `path` | `string` | Path to reference file (VCF/BCF or TSV). If both `path` and `vcf_path`/`tsv_path` are provided, `vcf_path` takes precedence, then `path`, then `tsv_path`. | `None` |
| `vcf_path` | `string` | VCF/BCF file path. Overrides `path` and `tsv_path`. For sweep mode, this is the recommended way to specify VCF/BCF files. | `None` |
| `tsv_path` | `string` | Precomputed lookup TSV file path. If provided, will use this directly instead of extracting from VCF. | `None` |
| `lookup_path` | `string` | Path to save/load the extracted lookup table. If the file exists and is valid, it will be reused. Useful for caching lookup tables between runs. | `None` |
| `chr_dict` | `dict` | A dictionary for converting chromosome numbers (1-25) to chromosome names in the VCF files. For example, the notation in dbSNP vcf file is based on RefSeq (like NC_000001.10). `gwaslab` provides built-in conversion dictionaries: `gl.get_number_to_NC(build="19")` for hg19 and `gl.get_number_to_NC(build="38")` for hg38. | `None` |
| `threads` | `int` | Number of threads for bcftools operations and parallel processing. More threads can significantly speed up lookup table extraction. | `6` |
| `overwrite` | `string` | Overwrite mode for rsID assignment. Options: `"all"`, `"invalid"`, or `"empty"`. See [Overwrite Modes](#overwrite-modes) for details. | `"empty"` |
| `convert_to_bcf` | `bool` | If True, convert VCF to BCF before processing. BCF format is more efficient for large files. | `False` |
| `strip_info` | `bool` | If True, strip INFO fields when converting VCF to BCF. Reduces file size and speeds up processing. | `True` |

!!! note "Sweep Mode Requirements"
    - **bcftools** must be installed and available in your PATH
    - For VCF files, tabix/csi indexing is recommended for optimal performance
    - Sweep mode creates temporary lookup TSV files that can be cached using `lookup_path`

## Chromosome Dictionary

For VCF files from dbSNP, you need to convert chromosome notation from numbers (1-25) to RefSeq names (NC_000001.10, etc.). GWASLab provides built-in conversion dictionaries:

```
# For hg19 (GRCh37)
gl.get_number_to_NC(build="19")
# {1: 'NC_000001.10', 2: 'NC_000002.11', 3: 'NC_000003.11', ...}

# For hg38 (GRCh38)
gl.get_number_to_NC(build="38")
# {1: 'NC_000001.11', 2: 'NC_000002.12', 3: 'NC_000003.12', ...}
```

These dictionaries map chromosome numbers (1-25) to RefSeq chromosome names used in dbSNP VCF files.

## Overwrite Modes

The `overwrite` parameter controls which existing rsID values should be replaced:

- **`"empty"`** (default): Only assign rsID for variants with missing/NA rsID values. This is the safest option and preserves existing rsID assignments.
- **`"invalid"`**: Assign rsID for variants with invalid rsID format (not matching the pattern `rs[0-9]+`). Useful for fixing incorrectly formatted rsIDs.
- **`"all"`**: Overwrite all rsIDs for eligible variants, regardless of existing values. Use with caution as this will replace all existing rsID assignments.

!!! note "STATUS Code Filtering"
    Only variants with proper STATUS codes (digit 4 = 0, digit 5 = 0-4) are eligible for rsID assignment. This ensures that only standardized and normalized variants receive rsID assignments. Run `.basic_check()` first to ensure proper STATUS codes.

## Examples

### Example 1: Basic rsID Assignment with TSV File

```
# Download reference SNPID-rsID table first
gl.download_ref("1kg_dbsnp151_hg19_auto") 

# Load sumstats
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
             n="N")

# Run basic_check first to standardize and normalize
mysumstats.basic_check() 

# If your SNPID is like 1:725932_G_A, you can use fix_id to fix the separator
mysumstats.fix_id(fixsep=True)

# rsID annotation using TSV file (fast, covers common variants)
mysumstats.assign_rsid(
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    threads=2
)
```

![image](https://user-images.githubusercontent.com/40289485/211799713-032b3571-ae01-4307-909b-973cdf043e17.png)

![image](https://user-images.githubusercontent.com/40289485/211799561-dddf7649-09fc-4eb5-b8e3-d7301a8944d1.png)

### Example 2: Complete rsID Assignment with Both TSV and VCF Files

```
# rsID annotation using both TSV and VCF files
# This maximizes annotation coverage
mysumstats.assign_rsid(
    threads=2,
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    ref_rsid_vcf="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    chr_dict=gl.get_number_to_NC(build="19"),
    overwrite="empty"  # Only fill missing rsIDs
)
```

The log output shows the annotation process:

```
Start to assign rsID using reference file...
 -Current Dataframe shape : 10000 x 12
 -SNPID-rsID text file: /home/yunye/.gwaslab/1kg_dbsnp151_hg19_auto.txt.gz
 -10000 rsID could be possibly fixed...
 -Setting block size: 5000000
 -Loading block: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
 -rsID Annotation for 58 need to be fixed!
 -Annotated 9942 rsID successfully!
Start to assign rsID using reference file...
 -Current Dataframe shape : 10000 x 13
 -CPU Cores to use : 2
 -Reference VCF file: /path/to/dbsnp/GCF_000001405.25.vcf.gz
 -Assigning rsID based on CHR:POS and REF:ALT/ALT:REF...
 -rsID Annotation for 1 need to be fixed!
 -Annotated 57 rsID successfully!
```

As you can see, the SNPID-rsID table (`1kg_dbsnp151_hg19_auto`) annotated 9942 rsIDs, and the large reference VCF file (from dbSNP) annotated an additional 57 rare rsIDs that were not in the TSV file.

![image](https://user-images.githubusercontent.com/40289485/211800319-68f33eaa-4c48-4ba4-aa52-afb9ac145dee.png)

### Example 3: Using Overwrite Modes

```
# Only fill missing rsIDs (default, safest)
mysumstats.assign_rsid(
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    overwrite="empty"
)

# Fix invalid rsID formats (e.g., "rs123" -> "rs123456")
mysumstats.assign_rsid(
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    overwrite="invalid"
)

# Overwrite all rsIDs (use with caution!)
mysumstats.assign_rsid(
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    overwrite="all"
)
```

### Example 4: Using .assign_rsid2() with Sweep Mode

```
# For large datasets with VCF/BCF files, use sweep mode
mysumstats.basic_check()

# Option 1: Using vcf_path (recommended)
mysumstats.assign_rsid2(
    vcf_path="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    chr_dict=gl.get_number_to_NC(build="19"),
    threads=8,
    overwrite="empty"
)

# Option 2: Using path parameter
mysumstats.assign_rsid2(
    path="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    chr_dict=gl.get_number_to_NC(build="19"),
    threads=8
)

# Option 3: Using pre-computed lookup table (fastest for repeated runs)
mysumstats.assign_rsid2(
    tsv_path="/path/to/cached_lookup.tsv.gz",
    overwrite="empty"
)

# Option 4: Save lookup table for future use
mysumstats.assign_rsid2(
    vcf_path="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    lookup_path="/path/to/cached_lookup.tsv.gz",  # Save for reuse
    chr_dict=gl.get_number_to_NC(build="19"),
    threads=8
)
```

### Example 5: Using .assign_rsid2() with BCF Conversion

```
# Convert VCF to BCF first for better performance
mysumstats.basic_check()
mysumstats.assign_rsid2(
    vcf_path="/path/to/dbsnp/GCF_000001405.25.vcf.gz",
    convert_to_bcf=True,  # Convert to BCF format
    strip_info=True,       # Strip INFO fields to reduce size
    chr_dict=gl.get_number_to_NC(build="19"),
    threads=8
)
```

## Performance Tips

1. **Use TSV first**: The TSV file is much faster and covers most common variants. Only use VCF for rare variants.
2. **Choose the right method**: 
   - Use `.assign_rsid()` for small to medium datasets (< 1 million variants)
   - Use `.assign_rsid2()` (sweep mode) for large datasets (millions of variants) with VCF/BCF files
3. **Parallel processing**: Increase `threads` for faster processing, especially with large VCF files. Sweep mode benefits significantly from more threads.
4. **Indexed VCF files**: Ensure your VCF file has a `.tbi` (tabix) or `.csi` (csi) index for faster queries.
5. **Chunk size**: For `.assign_rsid()`, adjust `chunksize` based on available memory. Larger chunks may be faster but use more memory.
6. **Cache lookup tables**: When using `.assign_rsid2()`, save the lookup table using `lookup_path` to reuse it in subsequent runs.
7. **BCF format**: For very large VCF files, consider using `convert_to_bcf=True` to convert to BCF format, which is more efficient.

## Important Notes

- Always run `.basic_check()` before assigning rsIDs to ensure proper standardization and normalization
- Both functions only assign rsIDs to variants with proper STATUS codes (standardized and normalized)
- For VCF files, tabix/csi indexing is highly recommended for performance
- The TSV file approach is much faster and should be used first for common variants
- **Sweep mode (`.assign_rsid2()`) requires bcftools** - ensure it's installed and in your PATH
- Sweep mode creates temporary lookup tables that can be cached and reused using `lookup_path`
- For very large datasets, sweep mode can be significantly faster than standard mode
