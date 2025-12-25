# GTF I/O in GWASLab

GWASLab provides fast and efficient functions for reading GTF (Gene Transfer Format) files, with support for filtering, attribute expansion, and automatic chromosome handling.

## Supported File Formats

GWASLab supports:
- Plain text GTF files (`.gtf`)
- Gzip-compressed GTF files (`.gtf.gz`)

## Reading GTF Files

### `read_gtf()`

Fast GTF file reader using Polars for improved performance. Returns a pandas DataFrame for compatibility.

```python
from gwaslab.io.io_gtf import read_gtf

# Read entire GTF file
gtf = read_gtf("annotation.gtf.gz")

# Read with column selection
gtf = read_gtf(
    "annotation.gtf.gz",
    usecols=["seqname", "start", "end", "strand", "feature", "gene_id", "gene_name"]
)

# Filter by chromosome early for speed
gtf = read_gtf("annotation.gtf.gz", chrom="1")

# Filter by feature type
gtf = read_gtf("annotation.gtf.gz", features={"gene", "transcript", "exon"})
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `filepath_or_buffer` | `str` or buffer | Path to GTF file (may be gzip compressed) or buffer object | Required |
| `usecols` | `list[str]` or `None` | Restrict which columns are loaded. If None, load all columns. Common columns: seqname, start, end, strand, feature, gene_biotype, gene_id, gene_name | `None` |
| `features` | `set[str]` or `None` | Drop rows which aren't one of the features in the supplied set (e.g., {'gene', 'transcript', 'exon'}) | `None` |
| `chrom` | `str` or `None` | Filter by chromosome/seqname early for speed. If None, load all chromosomes. Can be chromosome number (e.g., "1", "23" for X) or name (e.g., "X", "chr1"). For X chromosome, can use "X", "chrX", or "23" | `None` |
| `expand_attribute_column` | `bool` | Expand the 'attribute' column into separate columns (default: True) | `True` |
| `infer_biotype_column` | `bool` | Infer biotype from 'source' column if gene_biotype/transcript_biotype missing | `False` |

**Returns:**
- `pandas.DataFrame` - DataFrame containing parsed GTF data

**Standard GTF Columns:**
- `seqname` - Sequence name (chromosome)
- `source` - Source of annotation
- `feature` - Feature type (gene, transcript, exon, etc.)
- `start` - Start position (1-based)
- `end` - End position (1-based)
- `score` - Score
- `strand` - Strand (+ or -)
- `frame` - Reading frame (0, 1, 2, or .)
- `attribute` - Additional attributes (expanded by default)

**Expanded Attribute Columns:**
When `expand_attribute_column=True`, common attributes from the attribute column are automatically expanded into separate columns, such as:
- `gene_id`
- `gene_name`
- `gene_biotype`
- `transcript_id`
- `transcript_name`
- `transcript_biotype`
- And other attributes present in the file

!!! tip "Performance Optimization"
    Filtering by chromosome (`chrom`) or feature type (`features`) early in the reading process significantly improves performance for large GTF files, as it reduces the amount of data that needs to be processed.

---

### `read_gtf_file()`

Convenience function to read a GTF file with commonly used columns.

```python
from gwaslab.io.io_gtf import read_gtf_file

gtf = read_gtf_file("annotation.gtf.gz")
```

This is equivalent to:
```python
read_gtf(
    gtf_path,
    usecols=[
        "seqname",
        "start",
        "end",
        "strand",
        "feature",
        "gene_biotype",
        "gene_id",
        "gene_name",
    ],
)
```

---

### `get_gtf()`

Get GTF data for a specific chromosome from built-in reference datasets. Optimized to filter by chromosome early during file reading for speed.

```python
from gwaslab.io.io_gtf import get_gtf

# Get Ensembl hg19 GTF for chromosome 1
gtf = get_gtf(chrom="1", build="19", source="ensembl")

# Get Ensembl hg38 GTF for chromosome X
gtf = get_gtf(chrom="X", build="38", source="ensembl")

# Get RefSeq hg38 GTF for chromosome 1
gtf = get_gtf(chrom="1", build="38", source="refseq")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `chrom` | `str` or `int` | Chromosome number or name (e.g., "1", "X", 23) | Required |
| `build` | `str` | Genome build ("19" or "38") | `"19"` |
| `source` | `str` | Data source ("ensembl" or "refseq") | `"ensembl"` |

**Returns:**
- `pandas.DataFrame` - GTF data for the specified chromosome

**Supported builds and sources:**
- `build="19"`, `source="ensembl"` - Ensembl hg19 GTF
- `build="38"`, `source="ensembl"` - Ensembl hg38 GTF
- `build="19"`, `source="refseq"` - RefSeq hg19 GTF
- `build="38"`, `source="refseq"` - RefSeq hg38 GTF

!!! note "Automatic Download"
    If the reference GTF file is not available locally, GWASLab will automatically download it on first use.

---

## Processing GTF Files

### `gtf_to_protein_coding()`

Extract protein-coding genes from a GTF file and save to a new file.

```python
from gwaslab.io.io_gtf import gtf_to_protein_coding
from gwaslab.info.g_Log import Log

protein_coding_path = gtf_to_protein_coding(
    "annotation.gtf.gz",
    log=Log(),
    verbose=True
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `gtfpath` | `str` | Path to input GTF file | Required |
| `log` | `gwaslab.g_Log.Log` | Logging object | `Log()` |
| `verbose` | `bool` | If True, print progress messages | `True` |

**Returns:**
- `str` - Path to the output protein-coding GTF file (`.protein_coding.gtf.gz`)

**Output:**
- Creates a new file with suffix `.protein_coding.gtf.gz` containing only protein-coding gene records
- If the output file already exists, returns the path without re-processing

---

### `gtf_to_all_gene()`

Extract all gene records from a GTF file and save to a new file.

```python
from gwaslab.io.io_gtf import gtf_to_all_gene
from gwaslab.info.g_Log import Log

all_gene_path = gtf_to_all_gene(
    "annotation.gtf.gz",
    log=Log(),
    verbose=True
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `gtfpath` | `str` | Path to input GTF file | Required |
| `log` | `gwaslab.g_Log.Log` | Logging object | `Log()` |
| `verbose` | `bool` | If True, print progress messages | `True` |

**Returns:**
- `str` - Path to the output all-gene GTF file (`.all_genes.gtf.gz`)

**Output:**
- Creates a new file with suffix `.all_genes.gtf.gz` containing all gene records
- If the output file already exists, returns the path without re-processing

---

## Chromosome Handling

GWASLab automatically handles different chromosome naming conventions:

- **Numeric chromosomes**: "1", "2", ..., "22"
- **Sex chromosomes**: "X", "Y" are automatically converted to "23", "24" internally
- **Prefixed chromosomes**: "chr1", "chrX" are handled automatically
- **RefSeq IDs**: For RefSeq GTF files, chromosome names are automatically converted to standard format

When filtering by chromosome, you can use any of these formats:
```python
# All of these work for chromosome X:
gtf = read_gtf("annotation.gtf.gz", chrom="X")
gtf = read_gtf("annotation.gtf.gz", chrom="chrX")
gtf = read_gtf("annotation.gtf.gz", chrom="23")
```

---

## Examples

!!! example "Basic GTF reading"
    ```python
    from gwaslab.io.io_gtf import read_gtf
    
    # Read entire GTF file
    gtf = read_gtf("annotation.gtf.gz")
    print(f"Total records: {len(gtf)}")
    print(f"Columns: {gtf.columns.tolist()}")
    ```

!!! example "Filter by chromosome and feature"
    ```python
    from gwaslab.io.io_gtf import read_gtf
    
    # Read only genes from chromosome 1
    genes = read_gtf(
        "annotation.gtf.gz",
        chrom="1",
        features={"gene"}
    )
    print(f"Genes on chr1: {len(genes)}")
    ```

!!! example "Select specific columns"
    ```python
    from gwaslab.io.io_gtf import read_gtf
    
    # Read only essential columns
    gtf = read_gtf(
        "annotation.gtf.gz",
        usecols=["seqname", "start", "end", "gene_id", "gene_name", "gene_biotype"]
    )
    ```

!!! example "Get reference GTF data"
    ```python
    from gwaslab.io.io_gtf import get_gtf
    
    # Get Ensembl hg38 GTF for chromosome 1
    gtf = get_gtf(chrom="1", build="38", source="ensembl")
    
    # Filter for protein-coding genes
    protein_coding = gtf[gtf["gene_biotype"] == "protein_coding"]
    print(f"Protein-coding genes on chr1: {len(protein_coding)}")
    ```

!!! example "Extract protein-coding genes"
    ```python
    from gwaslab.io.io_gtf import gtf_to_protein_coding
    from gwaslab.info.g_Log import Log
    
    # Extract and save protein-coding genes
    protein_coding_path = gtf_to_protein_coding(
        "annotation.gtf.gz",
        log=Log(),
        verbose=True
    )
    
    # Read the extracted file
    from gwaslab.io.io_gtf import read_gtf
    protein_coding = read_gtf(protein_coding_path)
    ```

!!! example "Gene annotation workflow"
    ```python
    from gwaslab.io.io_gtf import read_gtf
    
    # Read genes with essential information
    gtf = read_gtf(
        "annotation.gtf.gz",
        features={"gene"},
        usecols=["seqname", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype"]
    )
    
    # Filter for protein-coding genes
    protein_coding = gtf[gtf["gene_biotype"] == "protein_coding"]
    
    # Get gene coordinates
    for _, gene in protein_coding.iterrows():
        print(f"{gene['gene_name']}: {gene['seqname']}:{gene['start']}-{gene['end']}")
    ```

---

## Performance Tips

1. **Filter early**: Use `chrom` and `features` parameters to filter data during reading, not after loading.

2. **Select columns**: Use `usecols` to load only the columns you need, reducing memory usage.

3. **Use built-in references**: For common reference genomes, use `get_gtf()` which automatically handles downloading and caching.

4. **Compressed files**: GWASLab automatically handles gzip-compressed GTF files, so you can work directly with `.gtf.gz` files.

5. **Polars backend**: The GTF reader uses Polars internally for fast file reading, then converts to pandas for compatibility.

---

## GTF Format Specification

GTF (Gene Transfer Format) is a tab-delimited text format for storing gene annotations. Each line represents one feature (gene, transcript, exon, etc.) with the following structure:

```python
seqname    source    feature    start    end    score    strand    frame    attribute
```

- **seqname**: Chromosome or contig name
- **source**: Annotation source (e.g., "ensembl", "refseq")
- **feature**: Feature type (gene, transcript, exon, CDS, etc.)
- **start**: Start position (1-based, inclusive)
- **end**: End position (1-based, inclusive)
- **score**: Score (often ".")
- **strand**: Strand (+ or -)
- **frame**: Reading frame (0, 1, 2, or ".")
- **attribute**: Semicolon-separated key-value pairs with additional information

For more details, see the [GTF format specification](https://www.ensembl.org/info/website/upload/gff.html).

