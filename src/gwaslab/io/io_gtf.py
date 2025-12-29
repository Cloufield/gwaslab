import pandas as pd
import polars as pl
from os import path
from os.path import exists
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_download import check_and_download
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

# GTF/GFF file suffix definitions
GTF_GFF_SUFFIXES = ('.gtf.gz', '.gff.gz', '.gtf', '.gff')

# GTF required columns
REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]


def _parse_frame(value):
    """Parse frame column, converting '.' to 0."""
    if value == "." or value is None:
        return 0
    try:
        return int(value)
    except (ValueError, TypeError):
        return 0

def read_gtf(
    filepath_or_buffer,
    usecols=None,
    features=None,
    chrom=None,
    expand_attribute_column=True,
    infer_biotype_column=False,
):
    """
    Fast GTF file reader using Polars.

    Returns a pandas DataFrame for compatibility.
    
    Parameters
    ----------
    filepath_or_buffer : str or buffer
        Path to GTF file (may be gzip compressed) or buffer object
    
    usecols : list of str or None
        Restrict which columns are loaded. If None, load all columns.
        Common columns: seqname, start, end, strand, feature, 
        gene_biotype, gene_id, gene_name
    
    features : set of str or None
        Drop rows which aren't one of the features in the supplied set
        (e.g., {'gene', 'transcript', 'exon'})
    
    chrom : str or None
        Filter by chromosome/seqname early for speed. If None, load all chromosomes.
        Can be chromosome number (e.g., "1", "23" for X) or name (e.g., "X", "chr1").
        For X chromosome, can use "X", "chrX", or "23".
        For Y chromosome, can use "Y", "chrY", or "24".
        For MT chromosome, can use "MT", "chrMT", "M", "chrM", or "25".
    
    expand_attribute_column : bool
        Expand the 'attribute' column into separate columns (default: True)
    
    infer_biotype_column : bool
        Infer biotype from 'source' column if gene_biotype/transcript_biotype missing
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed GTF data
    """
    import re
    
    if isinstance(filepath_or_buffer, str) and not exists(filepath_or_buffer):
        raise ValueError("GTF file does not exist: %s" % filepath_or_buffer)
    
    if features is not None:
        features = set(features)
    
    # Determine which attribute columns we need to extract
    restrict_attribute_columns = None
    if usecols is not None:
        standard_cols = set(REQUIRED_COLUMNS) - {"attribute"}
        restrict_attribute_columns = [c for c in usecols if c not in standard_cols]
        if not restrict_attribute_columns:
            restrict_attribute_columns = None
    
    # Read GTF file with Polars (much faster than pandas)
    # Read all columns as strings first to handle invalid values
    # Use infer_schema_length=0 to prevent type inference and read everything as string
    df = pl.read_csv(
        filepath_or_buffer,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=REQUIRED_COLUMNS,
        null_values=".",  # Only '.' is null, 'X' is valid chromosome name
        try_parse_dates=False,
        infer_schema_length=0,  # Don't infer types, read everything as string
    )
    
    if len(df) == 0:
        if usecols is not None:
            return pd.DataFrame(columns=usecols)
        return pd.DataFrame(columns=REQUIRED_COLUMNS)
    
    # Convert sex chromosomes to numeric in seqname column using ChromosomeMapper
    # Convert Polars DataFrame to pandas for ChromosomeMapper, then back
    df_pd_temp = df.select("seqname").to_pandas()
    mapper = ChromosomeMapper()
    mapper.detect_sumstats_format(df_pd_temp["seqname"])
    df_pd_temp["seqname"] = mapper.to_numeric(df_pd_temp["seqname"])
    # Update the Polars DataFrame
    df = df.with_columns(
        pl.Series("seqname", df_pd_temp["seqname"].astype(str))
    )
    
    # Early chromosome filtering for speed (before processing attributes)
    if chrom is not None:
        # Use ChromosomeMapper to normalize chromosome identifier
        # Convert to numeric format (which is what seqname column now contains)
        mapper = ChromosomeMapper()
        mapper.detect_sumstats_format(pd.Series([chrom]))
        chrom_numeric = mapper.to_numeric(chrom)
        chrom_numeric_str = str(chrom_numeric)
        
        # Also get string format for matching
        chrom_str = str(mapper.to_string(chrom_numeric))
        
        # Build list of possible values to match
        # After conversion above, seqname contains numeric values (23, 24, 25 for X, Y, MT)
        chrom_values = [chrom_numeric_str, chrom_str, str(chrom)]
        
        # Add common variations
        if chrom_str.upper() == "X":
            chrom_values.extend(["23", "chrX", "X"])
        elif chrom_str.upper() == "Y":
            chrom_values.extend(["24", "chrY", "Y"])
        elif chrom_str.upper() in ["MT", "M"]:
            chrom_values.extend(["25", "chrMT", "chrM", "MT", "M"])
        
        df = df.filter(pl.col("seqname").is_in(chrom_values))
        if len(df) == 0:
            if usecols is not None:
                return pd.DataFrame(columns=usecols)
            return pd.DataFrame(columns=REQUIRED_COLUMNS)
    
    # Parse frame column (convert '.' to 0, handle nulls)
    df = df.with_columns(
        pl.col("frame")
        .fill_null(".")
        .map_elements(_parse_frame, return_dtype=pl.Int8)
        .alias("frame")
    )
    
    # Convert numeric columns, handling nulls and invalid values
    # Replace invalid values (like 'X') with null before casting
    df = df.with_columns([
        # Start: replace invalid values with null, then cast to int64
        pl.when(pl.col("start").is_in(["X", ".", ""]) | pl.col("start").is_null())
        .then(None)
        .otherwise(pl.col("start"))
        .cast(pl.Int64)
        .alias("start"),
        # End: replace invalid values with null, then cast to int64
        pl.when(pl.col("end").is_in(["X", ".", ""]) | pl.col("end").is_null())
        .then(None)
        .otherwise(pl.col("end"))
        .cast(pl.Int64)
        .alias("end"),
        # Score: replace invalid values with null, then cast to float32
        pl.when(pl.col("score").is_in(["X", ".", ""]) | pl.col("score").is_null())
        .then(None)
        .otherwise(pl.col("score"))
        .cast(pl.Float32)
        .alias("score"),
    ])
    
    # Early feature filtering
    if features is not None:
        df = df.filter(pl.col("feature").is_in(list(features)))
        if len(df) == 0:
            if usecols is not None:
                return pd.DataFrame(columns=usecols)
            return pd.DataFrame(columns=REQUIRED_COLUMNS)
    
    # Expand attributes if needed
    if expand_attribute_column and "attribute" in df.columns:
        # Fix broken quotes first (handles some Ensembl GTF issues)
        df = df.with_columns(
            pl.col("attribute")
            .str.replace_all(';"', '"')
            .str.replace_all(";-", "-")
            .alias("attribute")
        )
        
        # Convert to pandas temporarily for attribute parsing (Polars string ops are limited)
        # This is still faster overall because Polars reads the file much faster
        df_pd = df.to_pandas()
        
        # Parse attributes using string operations (faster than regex)
        expanded_dict = {}
        attr_series = df_pd["attribute"]
        
        # Determine which attributes to extract
        if restrict_attribute_columns is None:
            # Find all unique keys in attributes
            all_keys = set()
            for attr_str in attr_series:
                if attr_str:
                    pairs = attr_str.split(';')
                    for pair in pairs:
                        pair = pair.strip()
                        if pair:
                            space_idx = pair.find(' ')
                            if space_idx > 0:
                                key = pair[:space_idx].strip()
                                all_keys.add(key)
            keys_to_extract = list(all_keys)
        else:
            keys_to_extract = list(restrict_attribute_columns)
        
        # Extract each attribute
        for key in keys_to_extract:
            values = []
            for attr_str in attr_series:
                if not attr_str:
                    values.append(None)
                    continue
                
                # Find all occurrences of this key
                found_values = []
                pairs = attr_str.split(';')
                for pair in pairs:
                    pair = pair.strip()
                    if pair.startswith(key + ' '):
                        # Extract value
                        value = pair[len(key):].strip()
                        # Remove quotes
                        if value.startswith('"') and value.endswith('"'):
                            value = value[1:-1]
                        elif value.startswith('"'):
                            value = value[1:]
                        found_values.append(value)
                
                if found_values:
                    values.append(','.join(found_values))
                else:
                    values.append(None)
            
            expanded_dict[key] = values
        
        # Drop attribute column and add expanded columns
        df_pd = df_pd.drop("attribute", axis=1)
        for column_name, values in expanded_dict.items():
            df_pd[column_name] = values
        
        # Convert back to Polars
        df = pl.from_pandas(df_pd)
    
    # Infer biotype column if requested
    if infer_biotype_column:
        if "protein_coding" in df["source"].unique().to_list():
            if "gene_biotype" not in df.columns:
                df = df.with_columns(pl.col("source").alias("gene_biotype"))
            if "transcript_biotype" not in df.columns:
                df = df.with_columns(pl.col("source").alias("transcript_biotype"))
    
    # Select only requested columns
    if usecols is not None:
        available_cols = [c for c in usecols if c in df.columns]
        df = df.select(available_cols)
    
    # Convert to pandas for compatibility
    return df.to_pandas()


def read_gtf_file(gtf_path):
    return read_gtf(
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

# Module-level cache for GTF data
_GTF_CACHE = {}

def get_gtf(chrom, build="19", source="ensembl"):
    """
    Get GTF data for a specific chromosome.
    
    Optimized to filter by chromosome early during file reading for speed.
    Results are cached for faster subsequent access.
    
    Parameters
    ----------
    chrom : str or int
        Chromosome number or name (e.g., "1", "X", 23)
    build : str
        Genome build ("19" or "38")
    source : str
        Data source ("ensembl" or "refseq")
    
    Returns
    -------
    pandas.DataFrame
        GTF data for the specified chromosome
    """
    # Create cache key
    cache_key = f"{build}_{source}_{chrom}"
    
    # Check cache first
    if cache_key in _GTF_CACHE:
        return _GTF_CACHE[cache_key].copy()
    
    gtf = None
    if source == "ensembl":
        if build == "19":
            data_path = check_and_download("ensembl_hg19_gtf")
            # Filter by chromosome early for speed
            gtf = read_gtf(
                data_path,
                chrom=str(chrom),  # Early filtering by chromosome
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
        if build == "38":
            data_path = check_and_download("ensembl_hg38_gtf")
            # Filter by chromosome early for speed
            gtf = read_gtf(
                data_path,
                chrom=str(chrom),  # Early filtering by chromosome
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
    if source == "refseq":
        # Use ChromosomeMapper to convert chromosome to NCBI RefSeq ID format
        mapper = ChromosomeMapper(species="homo sapiens", build=build)
        mapper.detect_sumstats_format(pd.Series([chrom]))
        chrom_NC = mapper.to_nc(chrom)
        
        if build == "19":
            data_path = check_and_download("refseq_hg19_gtf")
            # Filter by chromosome early for speed
            gtf = read_gtf(
                data_path,
                chrom=chrom_NC,  # Early filtering by chromosome (NC format)
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
            # Convert seqname back to chromosome number using ChromosomeMapper
            mapper_back = ChromosomeMapper(species="homo sapiens", build=build)
            mapper_back.detect_sumstats_format(gtf["seqname"])
            gtf["seqname"] = mapper_back.to_string(gtf["seqname"]).astype(str)
        if build == "38":
            data_path = check_and_download("refseq_hg38_gtf")
            # Filter by chromosome early for speed
            gtf = read_gtf(
                data_path,
                chrom=chrom_NC,  # Early filtering by chromosome (NC format)
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
            # Convert seqname back to chromosome number using ChromosomeMapper
            mapper_back = ChromosomeMapper(species="homo sapiens", build=build)
            mapper_back.detect_sumstats_format(gtf["seqname"])
            gtf["seqname"] = mapper_back.to_string(gtf["seqname"]).astype(str)
    if gtf is None:
        gtf = pd.DataFrame(
            columns=[
                "seqname",
                "start",
                "end",
                "strand",
                "feature",
                "gene_biotype",
                "gene_id",
                "gene_name",
            ]
        )
    
    # Cache the result for future use
    _GTF_CACHE[cache_key] = gtf.copy()
    
    return gtf

def gtf_to_protein_coding(gtfpath, log=Log(), verbose=True):
    protein_coding_path = gtfpath[:-6] + "protein_coding.gtf.gz"
    if not path.isfile(protein_coding_path):
        log.write(
            " - Extracting protein_coding genes from {}".format(gtfpath),
            verbose=verbose,
        )
        gtf = read_gtf(
            gtfpath,
            usecols=["feature", "gene_biotype", "gene_id", "gene_name"],
        )
        gene_list = (
            gtf.loc[
                (gtf["feature"] == "gene") & (gtf["gene_biotype"] == "protein_coding"),
                "gene_id",
            ]
            .values
        )
        log.write(
            " - Loaded {} protein_coding genes.".format(len(gene_list)),
            verbose=verbose,
        )
        gtf_raw = pd.read_csv(gtfpath, sep="\t", header=None, comment="#", dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[gtf_raw["_gene_id"].isin(gene_list), :]
        gtf_raw = gtf_raw.drop("_gene_id", axis=1)
        log.write(
            " - Extracted records are saved to : {} ".format(protein_coding_path),
            verbose=verbose,
        )
        gtf_raw.to_csv(protein_coding_path, header=None, index=None, sep="\t")
    return protein_coding_path

def gtf_to_all_gene(gtfpath, log=Log(), verbose=True):
    all_gene_path = gtfpath[:-6] + "all_genes.gtf.gz"
    if not path.isfile(all_gene_path):
        log.write(" - Extracting genes from {}".format(gtfpath), verbose=verbose)
        gtf = read_gtf(
            gtfpath,
            usecols=["feature", "gene_biotype", "gene_id", "gene_name"],
        )
        gene_list = gtf.loc[gtf["feature"] == "gene", "gene_id"].values
        log.write(" - Loaded {} genes.".format(len(gene_list)), verbose=verbose)
        gtf_raw = pd.read_csv(gtfpath, sep="\t", header=None, comment="#", dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[gtf_raw["_gene_id"].isin(gene_list), :]
        gtf_raw = gtf_raw.drop("_gene_id", axis=1)
        log.write(
            " - Extracted records are saved to : {} ".format(all_gene_path),
            verbose=verbose,
        )
        gtf_raw.to_csv(all_gene_path, header=None, index=None, sep="\t")
    return all_gene_path
