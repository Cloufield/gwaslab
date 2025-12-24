# Assigning CHR/POS

GWASLab can update **CHR**/**POS** using a pre-processed **rsID**-to-**CHR**:**POS** mapping table stored in HDF5 format. The HDF5-based approach provides optimized performance for processing large datasets with parallel processing capabilities.

!!! note "Note on .rsid_to_chrpos()"
    The old TSV-based implementation of `.rsid_to_chrpos()` has been removed due to poor performance. Please use `.rsid_to_chrpos2()` with HDF5 files instead, which provides significantly better performance for large datasets.

## .rsid_to_chrpos2()

**Available since v3.4.31**

`.rsid_to_chrpos2()` is a high-performance function to assign **CHR** and **POS** based on **rsID** using a HDF5 file derived from dbSNP reference VCF files. This method is optimized for large datasets and uses parallel processing for fast lookups.

### How it works

The function uses pre-processed HDF5 files (one per chromosome) that contain **rsID**-to-**POS** mappings organized into 10 groups based on modulo 10 grouping (`rsID % 10`). This grouping strategy:

- **Creates approximately equal-sized groups** for balanced data distribution
- **Enables efficient parallel lookups** by processing chromosome-group combinations independently
- **Minimizes memory usage** by loading only relevant groups into memory
- **Supports fast random access** to specific rsID ranges using pandas index operations
- **Separate files per chromosome** allow true parallel writes without locking overhead
- **CHR not stored** - assigned from chromosome number (extracted from filename or provided parameter)
- **Index-based matching** - uses pandas index intersection for fast O(1) average-case lookups
- **Smart processing** - if **CHR** column is available, only processes relevant chromosome files (much faster)

### Reference VCF from dbSNP

First, download reference VCF file from dbSNP ftp site.

```
# For example, dbSNP v155 hg19
GCF_000001405.25.gz (24G)
GCF_000001405.25.gz.tbi (2.9M)
```

**Requirements:**

- The VCF file must be indexed (`.tbi` for bgzip-compressed VCF, or `.csi` for BCF format)
- **`bcftools` must be installed and available in your PATH** - the function uses `bcftools query` to extract data from VCF files
  - Install bcftools: `conda install -c bioconda bcftools` or download from https://www.htslib.org/download/
  - Verify installation: `bcftools --version`

### gl.process_vcf_to_hfd5()

**From gwaslab v3.4.42,  `process_ref_vcf()` was renamed to `process_vcf_to_hfd5()`**

```
gl.process_vcf_to_hfd5()
```

**Prerequisites:**

- **`bcftools` must be installed** - this function uses `bcftools query` to extract POS and rsID data from VCF files
  - Install: `conda install -c bioconda bcftools` or from https://www.htslib.org/download/
  - Verify: `bcftools --version`

Process the VCF file and convert it to HDF5 files (one per chromosome) using `.process_vcf_to_hfd5()`. This step may take 1-3 hours for a 30GB VCF file, depending on your system. Each chromosome is written to a separate HDF5 file for maximum parallel processing speed.

| Option      | DataType | Description                                                 | Default              |
|-------------|----------|-------------------------------------------------------------|----------------------|
| `vcf`       | `string` | the path to dbSNP VCF file (must be indexed with .tbi/.csi) | -                    |
| `directory` | `string` | the directory where you want output the converted HDF5 file | the same as VCF file |
| `chr_dict`  | `dict`   | dictionary for chromosome mapping (e.g., NC_000001.10 -> 1) | `None`               |
| `complevel` | `int`    | compression level for HDF5 (0-9, higher = more compression) | `3`                  |
| `threads`   | `int`    | number of threads for parallel processing                   | `1`                  |
| `chr_list`  | `list`   | list of chromosomes to process (None = all 1-25)            | `None`               |
| `overwrite` | `bool`   | if True, overwrite existing HDF5 files; if False, skip existing files | `False`              |

**Example:**

```
directory="/home/yunye/work/gwaslab/examples/vcf_hd5/"
vcf = "/home/yunye/CommonData/Reference/ncbi_dbsnp/ncbi_dbsnp/db155/GCF_000001405.25.gz"

gl.process_vcf_to_hfd5(vcf=vcf,
                   directory=directory,
                   chr_dict=gl.get_NC_to_number(build="19"),
                   complevel=3,
                   threads=6)

# To overwrite existing HDF5 files, set overwrite=True:
# gl.process_vcf_to_hfd5(vcf=vcf,
#                    directory=directory,
#                    chr_dict=gl.get_NC_to_number(build="19"),
#                    complevel=3,
#                    threads=6,
#                    overwrite=True)

2025/12/24 14:35:36 Start to process VCF file to HDF5 using bcftools:
2025/12/24 14:35:36 -Reference VCF path: /home/yunye/CommonData/dbsnp/GCF_000001405.25.gz
2025/12/24 14:35:36 -Grouping method: modulo 10 (creating 10 groups per chromosome)
2025/12/24 14:35:36 -Output format: separate HDF5 file per chromosome (no locking, true parallel writes)
2025/12/24 14:35:36 -Compression level: 3
2025/12/24 14:35:36 -Number of threads: 6
2025/12/24 14:35:36 -Overwrite existing files: False
2025/12/24 14:35:36 -HDF5 Output directory: /home/yunye/CommonData/dbsnp
2025/12/24 14:35:36 -HDF5 File pattern: GCF_000001405.25.chr{chr_num}.rsID_CHR_POS_mod10.h5
2025/12/24 14:35:36 -Processing chromosomes in parallel...
```

### HDF5 File Structure

The generated HDF5 files are organized as **one file per chromosome**, each containing 10 groups (group_0 through group_9):

**File Organization:**

- One HDF5 file per chromosome: `{vcf_file_name}.chr{chr_num}.rsID_CHR_POS_mod10.h5`
- Example: `GCF_000001405.25.chr1.rsID_CHR_POS_mod10.h5`, `GCF_000001405.25.chr2.rsID_CHR_POS_mod10.h5`, etc.
- CHR is **not stored** in the HDF5 file - it's encoded in the filename (e.g., `chr1` → CHR = 1)

**Group Organization:**

- Groups are named `group_0`, `group_1`, ..., `group_9`
- Each rsID is assigned to a group based on: `group_id = rsID_number % 10`
- For example: `rs123456789` → `789 % 10 = 9` → stored in `group_9`

**Data Structure within each group:**

Each group contains a DataFrame with **rsn as index** and **POS as column** (CHR is not stored):

| Structure | Data Type | Description                                    | Example      |
|-----------|-----------|------------------------------------------------|--------------|
| Index: `rsn` | `int64`   | rsID number (without "rs" prefix) - used as index for fast matching | `123456789`  |
| Column: `POS` | `int32`   | Base pair position                             | `10177`      |

**Note**: The `rsn` column is stored as the DataFrame index, enabling fast index-based matching using pandas index intersection operations.

**Example of data storage:**

```
HDF5 Directory: /path/to/dbsnp/
├── GCF_000001405.25.chr1.rsID_CHR_POS_mod10.h5
│   ├── group_0
│   │   └── DataFrame: rsn (int64) as index, POS (int32) as column
│   │       Example rows:
│   │              POS
│   │       rsn        
│   │       100    10177
│   │       200    10235
│   │       300    10352
│   │       ...
│   │
│   ├── group_1
│   │   └── DataFrame: rsn (int64) as index, POS (int32) as column
│   │       ...
│   │
│   └── group_9
│       └── DataFrame: rsn (int64) as index, POS (int32) as column
│           Example rows:
│                  POS
│           rsn        
│           99     10616
│           199    10789
│           ...
│
├── GCF_000001405.25.chr2.rsID_CHR_POS_mod10.h5
│   └── (same structure with chromosome 2 data)
│
└── ... (one file per chromosome)
```

**Accessing HDF5 data directly:**

You can inspect the HDF5 file structure using pandas:

```
import pandas as pd
import glob
import os

# Find all chromosome HDF5 files
h5_dir = "/path/to/dbsnp/"
h5_files = glob.glob(os.path.join(h5_dir, "*.chr*.rsID_CHR_POS_mod10.h5"))

# Access a specific chromosome file (e.g., chromosome 1)
chr1_file = [f for f in h5_files if '.chr1.' in f][0]

with pd.HDFStore(chr1_file, mode='r') as store:
    # List all groups
    print("Available groups:", store.keys())
    
    # Access a specific group
    group_0 = store['group_0']
    print("\nGroup 0 structure:")
    print(group_0.head())
    print(f"\nGroup 0 shape: {group_0.shape}")
    print(f"Group 0 index name: {group_0.index.name}")
    print(f"Group 0 dtypes:\n{group_0.dtypes}")
    print("\nNote: CHR is not stored - extract from filename")
    print("Note: rsn is stored as index for fast matching")
    # Extract CHR from filename
    import re
    match = re.search(r'\.chr(\d+)\.', chr1_file)
    if match:
        chr_num = int(match.group(1))
        print(f"CHR from filename: {chr_num}")
    
    # Search for a specific rsID (e.g., rs123456789) using index
    rsn_to_find = 123456789
    group_id = rsn_to_find % 10  # = 9
    group_9 = store[f'group_{group_id}']
    # Use index-based lookup (faster than column filtering)
    if rsn_to_find in group_9.index:
        result = group_9.loc[[rsn_to_find]]
        print(f"\nFound rs{rsn_to_find} in chromosome {chr_num}:")
        print(result)
        print(f"CHR: {chr_num} (from filename)")
        print(f"POS: {result['POS'].iloc[0]}")
```

**Storage optimization:**

The HDF5 format uses optimized data types to minimize storage:
- `int32` for POS (4 bytes per value, range: up to ~2.1 billion)
- `int64` for rsn (8 bytes per value, sufficient for all rsID numbers)
- **CHR not stored** - encoded in filename (saves 1 byte per variant)

This results in approximately **12 bytes per variant** (8 + 4), plus HDF5 overhead and compression. With compression level 3, a typical dbSNP v155 file (~1 billion variants) may compress to ~9-14 GB. The separate files per chromosome also enable better parallel processing and reduce file size per chromosome.

### .rsid_to_chrpos2()

```
.rsid_to_chrpos2()
```

### Options
| Option                    | DataType | Description                                                          | Default |
|---------------------------|----------|----------------------------------------------------------------------|---------|
| `path`                    | `string` | the path to the HDF5 directory (containing chromosome-specific files) | -       |
| `ref_rsid_to_chrpos_hdf5`| `string` | Path to HDF5 directory/file (takes precedence over VCF path)        | `None`  |
| `ref_rsid_to_chrpos_vcf` | `string` | Path to VCF file (auto-generates HDF5 path from VCF location)      | `None`  |
| `threads`                 | `int`    | number of threads to use for parallel processing                     | `4`     |
| `build`                   | `string` | genome build version for **CHR** and **POS** (e.g., "19", "38")              | `"99"`  |
| `rsid`                    | `string` | column name containing **rsID** values                                  | `"rsID"`|
| `chrom`                   | `string` | column name for chromosome values to be updated                      | `"CHR"` |
| `pos`                     | `string` | column name for position values to be updated                        | `"POS"` |
| `status`                  | `string` | column name for status codes                                         | `"STATUS"`|
| `verbose`                 | `bool`   | if True, print progress messages                                    | `True`  |

**Processing Strategy:**

The function automatically chooses the optimal processing strategy based on your data:

- **If **CHR** column has values**: Processes by chromosome first (much faster)

  - Groups data by chromosome, then by group (modulo 10)
  - Only processes relevant chromosome HDF5 files
  - Example: If data has **CHR**=1, only checks `*.chr1.rsID_CHR_POS_mod10.h5`
  
- **If **CHR** column is missing/empty**: Searches across all chromosomes
  - Groups by group (modulo 10), then tries all chromosome files
  - Automatically deduplicates results (keeps best match per variant)
  - Slower but necessary when CHR is unknown

### Example

```
# Load summary statistics
mysumstats = gl.Sumstats("my_sumstats.txt")

# Before: missing CHR and POS
mysumstats.data[['rsID', 'CHR', 'POS']].head()

# Output:
#        rsID   CHR   POS
# 0  rs123456   NaN   NaN
# 1  rs234567   NaN   NaN
# 2  rs345678   NaN   NaN
# 3  rs456789   NaN   NaN
# 4  rs567890   NaN   NaN

# Option 1: Direct path to HDF5 directory (contains chromosome-specific files)
mysumstats.rsid_to_chrpos2(
    path="/home/yunye/work/gwaslab/examples/vcf_hd5/",
    threads=6
)

# Option 2: Using ref_rsid_to_chrpos_hdf5 parameter
mysumstats.rsid_to_chrpos2(
    ref_rsid_to_chrpos_hdf5="/home/yunye/work/gwaslab/examples/vcf_hd5/",
    threads=6
)

# Option 3: Using VCF path (auto-generates HDF5 path from VCF location)
mysumstats.rsid_to_chrpos2(
    ref_rsid_to_chrpos_vcf="/path/to/dbsnp/GCF_000001405.25.gz",
    threads=6
)

# After: CHR and POS assigned
mysumstats.data[['rsID', 'CHR', 'POS']].head()

# Output:
#        rsID  CHR      POS
# 0  rs123456    1   10177
# 1  rs234567    1   10235
# 2  rs345678    1   10352
# 3  rs456789    1   10505
# 4  rs567890    1   10511
```

**Processing workflow:**

1. **rsID preprocessing**: Extracts numeric **rsID** values (strips first 2 characters) and identifies valid **rsID**s
2. **HDF5 file discovery**: Finds all chromosome-specific HDF5 files in the directory (e.g., `*.chr*.**rsID**_**CHR**_**POS**_mod10.h5`)
3. **Processing strategy**:

   - **If **CHR** column available**: Groups data by chromosome first, then by group (modulo 10), processes each chromosome-group combination
   - **If **CHR** column not available**: Groups by group, searches across all chromosome files, then deduplicates results
4. **Group assignment**: Each **rsID** is assigned to one of 10 groups based on `rsID % 10`
5. **Parallel lookup**: Multiple threads process tasks in parallel (each task = chromosome + group combination)
6. **Index-based matching**: Uses pandas index intersection for fast matching (rsn as index in both sumstats and reference data)
7. **CHR assignment**: **CHR** is assigned from the chromosome number (extracted from filename or provided parameter)
8. **POS assignment**: **POS** is extracted from the matched reference data
9. **Deduplication**: When searching across all chromosomes, keeps the first match per variant (prioritizes matches with both **CHR** and **POS**)

**Performance characteristics:**

- **Speed**: Processing 10,000 variants typically takes 1-3 seconds with 4-6 threads
- **Memory**: Only loads relevant HDF5 groups into memory (one group at a time per thread)
- **Scalability**: Can handle datasets with millions of variants efficiently
- **Parallelization**: Uses multiprocessing for concurrent task processing (chromosome × group combinations)
- **Index-based matching**: Uses pandas index intersection for O(1) average-case lookups
- **Optimized design**: Separate files per chromosome enable true parallel writes and efficient lookups
- **Smart processing**: When CHR is available, only processes relevant chromosome files (much faster)

**Log output example:**

```
2025/12/24 17:30:44 Start to assign CHR and POS using rsIDs ...(v4.0.0)
2025/12/24 17:30:44  -Current Dataframe shape : 10000 x 11 ; Memory usage: 0.64 MB
2025/12/24 17:30:44  -Source hdf5 file: /path/to/dbsnp/
2025/12/24 17:30:44  -Threads to use: 6
2025/12/24 17:30:44  -Grouping method: modulo 10 (rsID % 10)
2025/12/24 17:30:44  -Non-Valid rsIDs: 218
2025/12/24 17:30:44  -Valid rsIDs: 9782
2025/12/24 17:30:44  -Found HDF5 files for 25 chromosome(s)
2025/12/24 17:30:44  -Initiating CHR ... 
2025/12/24 17:30:44  -Initiating POS ... 
2025/12/24 17:30:44  -Searching across all chromosomes (CHR column not available or single-file format)
2025/12/24 17:30:44  -Created 250 processing tasks: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
2025/12/24 17:32:44  -Variants matched: 1922 / 9782
2025/12/24 17:32:44  -Merging group data... 
2025/12/24 17:32:44  -Append data... 
2025/12/24 17:32:44 Start to fix chromosome notation (CHR) ...(v4.0.0)
2025/12/24 17:32:44  -Variants with standardized chromosome notation: 1922
2025/12/24 17:32:44  -Variants with NA chromosome notations: 8078
2025/12/24 17:32:44 Finished fixing chromosome notation (CHR).
2025/12/24 17:32:44 Start to fix basepair positions (POS) ...(v4.0.0)
2025/12/24 17:32:44  -Position bound:(0 , 250,000,000)
2025/12/24 17:32:44  -Removed variants outliers: 0
2025/12/24 17:32:44 Finished fixing basepair positions (POS).
2025/12/24 17:32:44 Finished assigning CHR and POS using rsIDs.
```

### Tips for optimal performance

1. **Use appropriate thread count**: Set `threads` to match your CPU cores (typically 4-8 for most systems)
2. **Pre-generate HDF5 files**: Process VCF files once and reuse the HDF5 directory for multiple datasets
3. **Compression level**: Use `complevel=3` (default) for a good balance between file size and processing speed
4. **Memory considerations**: The function loads one HDF5 group per thread, so memory usage scales with thread count
5. **Large datasets**: For datasets with >10M variants, consider processing in batches or using more threads
6. **Separate files per chromosome**: The new format allows true parallel writes during processing, significantly faster than single-file approach
7. **CHR column advantage**: If your data already has CHR values, the function will process much faster by only checking relevant chromosome files
8. **Index-based matching**: Uses pandas index operations for efficient matching (faster than binary search for this use case)
9. **Task-based progress**: Task numbers are printed during processing to show progress (e.g., "1 2 3 4 5...")
10. **Automatic deduplication**: When searching across all chromosomes, the function automatically keeps the best match per variant