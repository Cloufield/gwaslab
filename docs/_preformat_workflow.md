# `_preformat` Workflow and Priority

## Overview
The `_preformat` function loads and standardizes summary statistics data into GWASLab format. It follows a strict 9-step workflow where each step builds upon the previous ones.

## Priority Order and Workflow

### **Step 1: Initialize Parameters** (Priority: Foundation)
**Purpose**: Set up basic data structures and validate inputs

**Actions**:
- Merge `readargs` and `kwreadargs` into a single dictionary
- Initialize `log` object if not provided
- Initialize empty lists for `other`, `exclude`, `include` if None

**Why First**: All subsequent steps depend on these initialized structures.

---

### **Step 2: Handle Parquet Format** (Priority: Early Exit)
**Purpose**: Special handling for parquet files that must load before other processing

**Actions**:
- If `tab_fmt == "parquet"`, load the parquet file immediately
- Convert to DataFrame for subsequent processing

**Why Second**: Parquet files need to be loaded before format detection and column mapping can occur. This is a special case that bypasses normal file reading.

---

### **Step 3: Load Format Configuration** (Priority: Configuration)
**Purpose**: Load predefined format mappings from formatbook if `fmt` is specified

**Actions**:
- If `fmt` is provided, load format metadata and rename dictionary from formatbook
- Apply format-specific settings:
  - Separator (`format_separator`)
  - NA values (`format_na`)
  - Comment characters (`format_comment`)
  - Additional columns (`format_other_cols`)
- Default separator to `"\t"` if not specified

**Priority Logic**:
1. **Formatbook** (`fmt` parameter) - Highest priority for column mappings
2. **User-specified headers** - Applied in Step 5
3. **Include/Exclude filters** - Applied in Step 7

**Why Third**: Format configuration provides the base rename dictionary that user parameters will extend/modify.

---

### **Step 4: Check Path and Header** (Priority: Discovery)
**Purpose**: Discover available columns in the input data

**Actions**:
- If `sumstats` is a string path:
  - Detect chromosome-separated files (pattern with `@`)
  - Handle compressed/uncompressed file variants
  - Read header to get `raw_cols` (available columns)
- If `sumstats` is a DataFrame:
  - Extract `raw_cols` directly from DataFrame columns
- Update `usecols` and `dtype_dictionary` based on format mappings found in `rename_dictionary`

**Why Fourth**: We need to know what columns are actually available before we can build mappings or load data.

---

### **Step 5: Build Column Mappings** (Priority: User Override)
**Purpose**: Map user-specified column names to GWASLab standard names

**Actions**:
- Process standard column mappings (SNPID, CHR, POS, EA, NEA, BETA, SE, P, etc.)
- Handle special cases:
  - EAF/NEAF: If `eaf` provided, use it; else if `neaf` provided, map to EAF (will be converted later)
  - Numeric columns (n, ncase, ncontrol, neff): Only add if string (int values handled later)
- Add `other` columns to usecols
- Check for reserved keyword overlaps

**Priority Logic**:
1. **Formatbook mappings** (from Step 3) - Base mappings
2. **User-specified mappings** (Step 5) - Override/extend formatbook
3. **Include/Exclude** (Step 7) - Final filter

**Why Fifth**: User parameters should override formatbook defaults, but we need formatbook loaded first to know the base state.

---

### **Step 6: Handle VCF Format Special Case** (Priority: Format-Specific)
**Purpose**: VCF format requires special column handling

**Actions**:
- If `fmt == "vcf"`:
  - Save current `usecols` as `vcf_usecols` (needed for post-processing)
  - Replace `usecols` with VCF fixed columns + study column
  - Study column comes from `study` parameter or defaults to column 9

**Why Sixth**: VCF format has a unique structure that requires different columns during loading vs. after parsing.

---

### **Step 7: Apply Include/Exclude Filters** (Priority: Final Selection)
**Purpose**: Apply user-specified column inclusion/exclusion filters

**Actions**:
- If `include` is specified:
  - Create reverse mapping (GWASLab name → original name)
  - Filter `usecols` to only include specified columns
- If `exclude` is specified:
  - Create reverse mapping
  - Remove excluded columns from `usecols`

**Priority Logic**:
1. Formatbook + User mappings → Full column set
2. Include filter → Subset to only included columns
3. Exclude filter → Remove excluded columns

**Why Seventh**: Filters must be applied after all mappings are built, but before data loading to avoid loading unnecessary columns.

---

### **Step 8: Load Data** (Priority: Data Acquisition)
**Purpose**: Actually load the data from file(s) or DataFrame

**Actions**:
- **If path string**:
  - **Multi-chromosome files** (`@` in path):
    - Load each chromosome file separately
    - Merge all chromosomes into single DataFrame
  - **Single file**:
    - If `chrom_pat`: Load with chromosome filtering (chunked)
    - If `snpid_pat`: Load with SNP ID pattern filtering (chunked)
    - Otherwise: Load full file
- **If DataFrame**:
  - Copy DataFrame
  - Apply dtype conversions based on `dtype_dictionary`

**Why Eighth**: Data loading happens after we know exactly which columns to load and how to interpret them.

---

### **Step 9: Post-Process Data** (Priority: Standardization)
**Purpose**: Transform loaded data into final GWASLab format

**Sub-steps** (in order):

1. **VCF Parsing** (if `fmt == "vcf"`):
   - Parse FORMAT column into individual columns
   - Restore `vcf_usecols` for final column selection

2. **Column Renaming**:
   - Apply `rename_dictionary` to standardize column names
   - Log column mapping information

3. **Add Constant Values**:
   - If `n`, `ncase`, `ncontrol` are integers, add as constant columns

4. **Process Status**:
   - Create STATUS column if missing (based on `build`)
   - Convert STATUS to integer type

5. **Process Alleles**:
   - Handle EA/NEA/REF/ALT relationships
   - Convert to categorical types for memory efficiency

6. **Process NEAF to EAF**:
   - Convert NEAF to EAF if needed (1 - NEAF)
   - Filter invalid frequency values

7. **Sort Columns**:
   - Reorder columns to GWASLab standard order

8. **Convert Data Types**:
   - Optimize data types for memory efficiency

9. **Ensure SNPID**:
   - Create SNPID from CHR:POS[:NEA:EA] if both rsID and SNPID missing

10. **Final Validation**:
    - Check data types
    - Check memory usage
    - Garbage collect

**Why Ninth**: All transformations happen after data is loaded, ensuring we work with actual data rather than metadata.

---

## Priority Summary

### Configuration Priority (Steps 1-3)
1. Initialize → 2. Parquet (if needed) → 3. Format config

### Mapping Priority (Steps 4-7)
4. Discover columns → 5. User mappings → 6. VCF special case → 7. Include/Exclude filters

### Data Priority (Steps 8-9)
8. Load data → 9. Post-process

## Key Design Principles

1. **Configuration Before Data**: All configuration (format, mappings, filters) happens before data loading
2. **Formatbook Override**: User parameters override formatbook defaults
3. **Include Before Exclude**: Include filter creates subset, exclude filter removes from that subset
4. **Transform After Load**: All data transformations happen after data is loaded
5. **Early Validation**: Column existence and type validation happens during discovery phase

## Example Flow

```
Input: sumstats="data.tsv", fmt="gwaslab", chrom="CHR", beta="BETA"

Step 1: Initialize → readargs={}, log=Log(), other=[], exclude=[], include=[]
Step 2: Skip (not parquet)
Step 3: Load formatbook → rename_dictionary={"CHR": "CHR", "BETA": "BETA", ...}
Step 4: Check path → raw_cols=["CHR", "POS", "BETA", "SE", "P", ...]
Step 5: Build mappings → rename_dictionary["CHR"]="CHR", rename_dictionary["BETA"]="BETA"
Step 6: Skip (not VCF)
Step 7: Apply filters → usecols=["CHR", "POS", "BETA", "SE", "P"]
Step 8: Load data → pd.read_table("data.tsv", usecols=usecols)
Step 9: Post-process → rename columns, add STATUS, process alleles, etc.
```

