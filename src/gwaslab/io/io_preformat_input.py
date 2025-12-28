import pandas as pd
import numpy as np
import scipy.stats as ss
import gzip
import os
import re
import gc
from gwaslab.bd.bd_common_data import get_format_dict
from gwaslab.qc.qc_fix_sumstats import _sort_column
from gwaslab.qc.qc_fix_sumstats import _process_build
from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import quick_convert_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_memory_usage
from gwaslab.qc.qc_reserved_headers import _check_overlap_with_reserved_keys
from gwaslab.info.g_vchange_status import STATUS_CATEGORIES
from gwaslab.info.g_Log import Log


#### Helper functions for _preformat #######################################################################

def _initialize_preformat_parameters(*, readargs=None, kwreadargs=None, log=None, other=None, exclude=None, include=None):
    """Initialize and validate parameters for preformatting."""
    if readargs is None:
        readargs = dict()
    if kwreadargs is None:
        kwreadargs = dict()
    readargs = readargs | kwreadargs
    
    if log is None:
        log = Log()
    
    if other is None:
        other = list()
    if exclude is None:
        exclude = list()
    if include is None:
        include = list()
    
    return readargs, log, other, exclude, include


def _load_format_config(*, fmt=None, readargs=None, other=None, log=None, verbose=False):
    """Load format configuration from formatbook if format is specified."""
    meta_data = None
    rename_dictionary = {}
    
    if fmt is not None:
        log.write("Start to load format from formatbook....", verbose=verbose)
        meta_data, rename_dictionary = get_format_dict(fmt)
        
        _print_format_info(fmt=fmt, meta_data=meta_data, rename_dictionary=rename_dictionary, 
                         verbose=verbose, log=log)
        
        # Apply format-specific readargs
        if "format_separator" in meta_data:
            if "sep" not in readargs:
                readargs["sep"] = meta_data["format_separator"]
            elif readargs["sep"] != meta_data["format_separator"]:
                log.write('  - format_separator will be changed to: "{}"'.format(readargs["sep"]), verbose=verbose)
        
        if "format_na" in meta_data:
            readargs["na_values"] = meta_data["format_na"]
        
        if "format_comment" in meta_data:
            readargs["comment"] = meta_data["format_comment"]
        
        if "format_other_cols" in meta_data:
            other += meta_data["format_other_cols"]
        
        if "sep" not in readargs:
            readargs["sep"] = "\t"
    
    return meta_data, rename_dictionary, readargs, other


def _build_column_mappings(*, column_params=None, rename_dictionary=None, usecols=None, dtype_dictionary=None, log=None, verbose=False):
    """
    Build column mappings from user-specified parameters.
    
    User-specified mappings override formatbook mappings. If a user mapping targets
    a standard name (e.g., "EA"), any existing mappings to that same target are removed
    to prevent duplicate column names.
    
    Parameters
    ----------
    column_params : dict
        Dictionary of column parameter names and values (e.g., {'snpid': 'SNP', 'chrom': 'CHR'})
    rename_dictionary : dict
        Dictionary to update with column mappings (may contain formatbook mappings)
    usecols : list
        List to update with columns to use
    dtype_dictionary : dict
        Dictionary to update with dtype specifications
    log : Log
        Logging object
    verbose : bool
        Verbose flag
    """
    def _remove_conflicting_mappings(target_name, rename_dictionary, usecols):
        """
        Remove any existing mappings that target the same standard name.
        This ensures user mappings override formatbook mappings.
        """
        # Find all keys that map to the target_name
        conflicting_keys = [k for k, v in rename_dictionary.items() if v == target_name]
        # Remove conflicting mappings
        for key in conflicting_keys:
            del rename_dictionary[key]
            # Also remove from usecols if present
            if key in usecols:
                usecols.remove(key)
    
    # Standard column mappings
    column_mappings = [
        ('snpid', "SNPID", None),
        ('rsid', "rsID", None),
        ('chrom', "CHR", "string"),
        ('pos', "POS", None),
        ('ea', "EA", "string"),
        ('nea', "NEA", "string"),
        ('ref', "REF", "string"),
        ('alt', "ALT", "string"),
        ('maf', "MAF", None),
        ('beta', "BETA", None),
        ('beta_95L', "BETA_95L", None),
        ('beta_95U', "BETA_95U", None),
        ('se', "SE", None),
        ('chisq', "CHISQ", None),
        ('z', "Z", None),
        ('q', "Q", None),
        ('p', "P", None),
        ('t', "T", None),
        ('f', "F", None),
        ('mlog10p', "MLOG10P", None),
        ('test', "TEST", None),
        ('info', "INFO", None),
        ('OR', "OR", None),
        ('OR_95L', "OR_95L", None),
        ('OR_95U', "OR_95U", None),
        ('HR', "HR", None),
        ('HR_95L', "HR_95L", None),
        ('HR_95U', "HR_95U", None),
        ('phet', "P_HET", None),
        ('i2', "I2", None),
        ('snpr2', "SNPR2", None),
        ('dof', "DOF", None),
        ('direction', "DIRECTION", None),
        ('status', "STATUS", "Int64"),
    ]
    
    for param_name, target_name, dtype in column_mappings:
        col_name = column_params.get(param_name)
        if col_name:
            # Remove any existing mappings to the same target (user overrides formatbook)
            _remove_conflicting_mappings(target_name, rename_dictionary, usecols)
            usecols.append(col_name)
            rename_dictionary[col_name] = target_name
            if dtype:
                dtype_dictionary[col_name] = dtype
    
    # Handle EAF/NEAF (special case: neaf can be converted to eaf)
    eaf = column_params.get('eaf')
    neaf = column_params.get('neaf')
    if eaf:
        # Remove any existing mappings to EAF
        _remove_conflicting_mappings("EAF", rename_dictionary, usecols)
        usecols.append(eaf)
        rename_dictionary[eaf] = "EAF"
    elif neaf:
        # Remove any existing mappings to EAF
        _remove_conflicting_mappings("EAF", rename_dictionary, usecols)
        usecols.append(neaf)
        rename_dictionary[neaf] = "EAF"
    
    # Handle numeric columns that can be str or int
    numeric_columns = [
        ('n', "N"),
        ('ncase', "N_CASE"),
        ('ncontrol', "N_CONTROL"),
        ('neff', "N_EFF"),
    ]
    
    for param_name, target_name in numeric_columns:
        col_name = column_params.get(param_name)
        if col_name and isinstance(col_name, str):
            # Remove any existing mappings to the same target
            _remove_conflicting_mappings(target_name, rename_dictionary, usecols)
            usecols.append(col_name)
            rename_dictionary[col_name] = target_name
    
    # Handle other columns
    other = column_params.get('other', [])
    if other:
        overlapped = _check_overlap_with_reserved_keys(other)
        if overlapped:
            log.warning("Columns with headers overlapping with GWASLab reserved keywords:{}".format(overlapped), verbose=verbose)
        usecols.extend(other)
        for col in other:
            rename_dictionary[col] = col


def _apply_column_filters(*, include=None, exclude=None, usecols=None, rename_dictionary=None, log=None, verbose=False):
    """Apply include/exclude filters to column list."""
    if include:
        # Create reverse mapping for faster lookup
        reverse_rename = {v: k for k, v in rename_dictionary.items()}
        usecols_new = [reverse_rename.get(i) for i in include if i in reverse_rename]
        # Use set intersection for faster filtering
        usecols_set = set(usecols)
        usecols_valid = [i for i in usecols_new if i in usecols_set]
        log.write(f' -Include columns :{",".join(usecols_valid)}', verbose=verbose)
        usecols = usecols_valid
    
    if exclude:
        # Create reverse mapping for faster lookup
        reverse_rename = {v: k for k, v in rename_dictionary.items()}
        exclude_cols = [reverse_rename.get(i) for i in exclude if i in reverse_rename]
        log.write(f' -Exclude columns :{",".join(exclude_cols)}', verbose=verbose)
        # Use set for faster removal
        exclude_set = set(exclude_cols)
        usecols = [i for i in usecols if i not in exclude_set]
    
    return usecols


def _load_sumstats_data(*, sumstats=None, inpath=None, inpath_chr_list=None, inpath_chr_num_list=None, 
                        usecols=None, dtype_dictionary=None, readargs=None, rename_dictionary=None,
                        chrom_pat=None, snpid_pat=None, log=None, verbose=False):
    """Load sumstats data from file path or DataFrame."""
    if isinstance(sumstats, str):
        return _load_sumstats_from_path(
            inpath=inpath, inpath_chr_list=inpath_chr_list, inpath_chr_num_list=inpath_chr_num_list,
            usecols=usecols, dtype_dictionary=dtype_dictionary, readargs=readargs, rename_dictionary=rename_dictionary,
            chrom_pat=chrom_pat, snpid_pat=snpid_pat, log=log, verbose=verbose
        )
    elif isinstance(sumstats, pd.DataFrame):
        return _load_sumstats_from_dataframe(
            sumstats=sumstats, dtype_dictionary=dtype_dictionary, rename_dictionary=rename_dictionary, 
            log=log, verbose=verbose
        )
    else:
        raise ValueError("Please input a path or a pd.DataFrame, and make sure it contain the columns.")


def _load_sumstats_from_path(*, inpath=None, inpath_chr_list=None, inpath_chr_num_list=None,
                             usecols=None, dtype_dictionary=None, readargs=None, rename_dictionary=None,
                             chrom_pat=None, snpid_pat=None, log=None, verbose=False):
    """Load sumstats from file path(s)."""
    if "@" in inpath:
        # Load multiple chromosome files
        log.write("Start to initialize gl.Sumstats from files with pattern :" + inpath, verbose=verbose)
        sumstats_chr_list = []
        for i in inpath_chr_list:
            log.write(" -Loading:" + i)
            skip_rows = _get_skip_rows(i)
            readargs_copy = readargs.copy()
            readargs_copy["skiprows"] = skip_rows
            explicit = {"usecols", "dtype_dictionary"}
            readargs_copy = {k: v for k, v in readargs_copy.items() if k not in explicit}
            sumstats_chr = pd.read_table(i, usecols=set(usecols), dtype=dtype_dictionary, **readargs_copy)
            sumstats_chr_list.append(sumstats_chr)
        log.write(" -Merging sumstats for chromosomes:", ",".join(inpath_chr_num_list), verbose=verbose)
        sumstats = pd.concat(sumstats_chr_list, axis=0, ignore_index=True)
        del sumstats_chr_list
        gc.collect()
        return sumstats
    else:
        # Load single file
        skip_rows = _get_skip_rows(inpath)
        readargs["skiprows"] = skip_rows
        log.write("Start to initialize gl.Sumstats from file :" + inpath, verbose=verbose)
        
        if chrom_pat is not None:
            return _load_single_chr(
                inpath=inpath, usecols=usecols, dtype_dictionary=dtype_dictionary, readargs=readargs,
                rename_dictionary=rename_dictionary, chrom_pat=chrom_pat, log=log, verbose=verbose
            )
        elif snpid_pat is not None:
            return _load_variants_with_pattern(
                inpath=inpath, usecols=usecols, dtype_dictionary=dtype_dictionary, readargs=readargs,
                rename_dictionary=rename_dictionary, snpid_pat=snpid_pat, log=log, verbose=verbose
            )
        else:
            explicit = {"usecols", "dtype_dictionary"}
            readargs = {k: v for k, v in readargs.items() if k not in explicit}
            return pd.read_table(inpath, usecols=set(usecols), dtype=dtype_dictionary, **readargs)


def _load_sumstats_from_dataframe(*, sumstats=None, dtype_dictionary=None, rename_dictionary=None, log=None, verbose=False):
    """Load and prepare sumstats from existing DataFrame."""
    log.write("Start to initialize gl.Sumstats from pandas DataFrame ...", verbose=verbose)
    sumstats = sumstats.copy()
    # Cache column set for faster lookups
    sumstats_cols = set(sumstats.columns)
    for key, value in dtype_dictionary.items():
        if key in sumstats_cols:
            astype = value
            if key in rename_dictionary and rename_dictionary[key] == "CHR":
                astype = "Int64"
            try:
                sumstats[key] = sumstats[key].astype(astype)
            except:
                sumstats[key] = sumstats[key].astype("string")
    return sumstats


def _postprocess_sumstats(*, sumstats=None, fmt=None, format_cols=None, study=None, vcf_usecols=None, usecols=None,
                         rename_dictionary=None, n=None, ncase=None, ncontrol=None, build=None, status=None, neaf=None, log=None, verbose=False):
    """Apply post-processing steps to loaded sumstats."""
    # Handle VCF format
    if fmt == "vcf":
        sumstats = _parse_vcf_study(sumstats, format_cols, study, vcf_usecols, log=log, verbose=verbose)
        usecols = vcf_usecols
    
    # Log column information
    usecols_set = set(usecols)
    converted_columns = [rename_dictionary[x] for x in usecols_set]
    log.write(" -Reading columns          :", ",".join(usecols_set), verbose=verbose)
    log.write(" -Renaming columns to      :", ",".join(converted_columns), verbose=verbose)
    log.write(" -Current Dataframe shape :", len(sumstats), " x ", len(sumstats.columns), verbose=verbose)
    
    # Rename columns
    sumstats = sumstats.rename(columns=rename_dictionary)
    
    # Add constant numeric values if provided
    if isinstance(n, int):
        sumstats["N"] = n
    if isinstance(ncase, int):
        sumstats["N_CASE"] = ncase
    if isinstance(ncontrol, int):
        sumstats["N_CONTROL"] = ncontrol
    
    # Process status
    sumstats = _process_status(sumstats=sumstats, build=build, status=status, log=log, verbose=verbose)
    
    # Process alleles
    sumstats = _process_allele(sumstats=sumstats, log=log, rename_dictionary=rename_dictionary, verbose=verbose)
    
    # Process NEAF to EAF
    if neaf is not None or ("NEAF" in sumstats.columns and "EAF" not in sumstats.columns):
        sumstats = _process_neaf(sumstats=sumstats, log=log, verbose=verbose)
    
    # Sort columns and convert datatypes
    sumstats = _sort_column(sumstats_obj=sumstats, log=log, verbose=verbose)
    sumstats = quick_convert_datatype(sumstats, log=log, verbose=verbose)
    
    # Create SNPID if missing
    sumstats = _ensure_snpid_column(sumstats=sumstats, log=log, verbose=verbose)
    
    # Final checks
    check_datatype(sumstats, log=log, verbose=verbose)
    gc.collect()
    check_dataframe_memory_usage(sumstats, log=log, verbose=verbose)
    
    return sumstats


def _ensure_snpid_column(*, sumstats=None, log=None, verbose=False):
    """Create SNPID column if both rsID and SNPID are absent."""
    sumstats_cols = sumstats.columns
    if "rsID" not in sumstats_cols and "SNPID" not in sumstats_cols:
        if "CHR" in sumstats_cols and "POS" in sumstats_cols:
            # Optimize string concatenation by converting once and using vectorized operations
            chr_str = sumstats["CHR"].astype("string")
            pos_str = sumstats["POS"].astype("string")
            if "EA" in sumstats_cols and "NEA" in sumstats_cols:
                sumstats["SNPID"] = chr_str + ":" + pos_str + ":" + sumstats["NEA"].astype("string") + ":" + sumstats["EA"].astype("string")
            else:
                sumstats["SNPID"] = chr_str + ":" + pos_str
            log.write(" -No rsID/SNPID found; created SNPID from CHR:POS[:NEA:EA]", verbose=verbose)
        else:
            sumstats["SNPID"] = pd.Series([None] * len(sumstats), dtype="string")
            log.warning(" -No rsID/SNPID and missing CHR/POS; created empty SNPID", verbose=verbose)
    return sumstats


#20221030
def _preformat(sumstats,
          fmt=None,
          tab_fmt="tsv",
          snpid=None,
          rsid=None,
          chrom=None,
          pos=None,
          ea=None,
          nea=None,
          ref=None,
          alt=None,
          eaf=None,
          neaf=None,
          maf=None,
          n=None,
          beta=None,
          se=None,
          chisq=None,
          z=None,
          f=None,
          t=None,
          p=None,
          q=None,
          mlog10p=None,
          test=None,
          info=None,
          OR=None,
          OR_95L=None,
          OR_95U=None,
          beta_95L=None,
          beta_95U=None,
          HR=None,
          HR_95L=None,
          HR_95U=None,
          i2=None,
          snpr2=None,
          phet=None,
          dof=None,
          ncase=None,
          ncontrol=None,
          neff=None,
          direction=None,
          status=None,
          study=None,
          trait=None,
          build=None,
          other=None,
          exclude=None,
          include=None,
          chrom_pat=None,
          snpid_pat=None,
          verbose=False,
          log=None,
          readargs=None,
          **kwreadargs):
    """
    Load and preformat summary statistics data into standardized GWASLab format.
    
    Workflow and Priority
    ---------------------
    The function follows a strict 9-step workflow where each step builds upon previous ones:
    
    **Phase 1: Configuration (Steps 1-3)**
    1. Initialize parameters - Set up basic data structures and validate inputs
    2. Handle parquet format - Special handling for parquet files (early exit)
    3. Load format configuration - Load predefined format mappings from formatbook if fmt specified
    
    **Phase 2: Mapping (Steps 4-7)**
    4. Check path and header - Discover available columns in input data
    5. Build column mappings - Map user-specified column names to GWASLab standard names
    6. Handle VCF format special case - VCF requires special column handling
    7. Apply include/exclude filters - Apply user-specified column inclusion/exclusion filters
    
    **Phase 3: Data (Steps 8-9)**
    8. Load data - Actually load data from file(s) or DataFrame
    9. Post-process data - Transform loaded data into final GWASLab format
    
    **Priority Order for Column Mappings:**
    1. Formatbook (fmt parameter) - Base mappings (Step 3)
    2. User-specified headers - Override/extend formatbook (Step 5)
    3. Include filter - Subset to specified columns (Step 7)
    4. Exclude filter - Remove excluded columns (Step 7)
    
    **Design Principles:**
    - Configuration before data: All setup happens before data loading
    - User overrides formatbook: User parameters take precedence over formatbook defaults
    - Transform after load: All data transformations happen after data is loaded
    - Early validation: Column existence and type validation happens during discovery phase

    Parameters
    ----------
    sumstats : str or pandas.DataFrame
        Input summary statistics, provided either as a file path or a DataFrame.
        When summary statistics are split by chromosome, a single path pattern
        may be supplied using the `@` symbol as a placeholder for the chromosome
        number. For example: "gwas/chr@.sumstats.gz" will load
        "gwas/chr1.sumstats.gz", "gwas/chr2.sumstats.gz", ... automatically.
    fmt : str, optional
        Format name to get predefined mapping if provided. (e.g., 'gwaslab', 'vcf').
    tab_fmt : str, default: 'tsv'
        Table format ('tsv', 'parquet').
    snpid : str, optional
        Column name for SNP identifiers in the input data. Expected format is CHR:POS:NEA:EA (e.g., 1:123:A:G), although the delimiter may vary depending on the source.
    rsid : str, optional
        Column name for rsID in the input data. Values should follow the standard "rs" prefix (e.g., rs12345).
    chrom : str, optional
        Column name for chromosome in input data.
    pos : str, optional
        Column name for position in input data.
    ea : str, optional
        Column name for effect allele in input data. (assuming alternative allele)
    nea : str, optional
        Column name for non-effect allele in input data. (assuming reference allele)
    eaf : str, optional
        Column name for effect allele frequency in input data.
    neaf : str, optional
        Column name for non-effect allele frequency in input data.
    maf : str, optional
        Column name for minor allele frequency in input data.
    n : str or int, optional
        Column name or constant value for sample size.
    beta : str, optional
        Column name for beta in input data.
    se : str, optional
        Column name for standard error in input data.
    chisq : str, optional
        Column name for chi-square in input data.
    z : str, optional
        Column name for z-score in input data.
    f : str, optional
        Column name for F-statistic in input data.
    t : str, optional
        Column name for T-statistic in input data.
    p : str, optional
        Column name for p-value in input data.
    q : str, optional
        Column name for Q-statistic in input data.
    mlog10p : str, optional
        Column name for -log10(p) in input data.
    test : str, optional
        Column name for test type in input data.
    info : str, optional
        Column name for imputation info in input data.
    OR : str, optional
        Column name for odds ratio in input data.
    OR_95L : str, optional
        Column name for lower 95% CI of OR in input data.
    OR_95U : str, optional
        Column name for upper 95% CI of OR in input data.
    beta_95L : str, optional
        Column name for lower 95% CI of beta in input data.
    beta_95U : str, optional
        Column name for upper 95% CI of beta in input data.
    HR : str, optional
        Column name for hazard ratio in input data.
    HR_95L : str, optional
        Column name for lower 95% CI of HR in input data.
    HR_95U : str, optional
        Column name for upper 95% CI of HR in input data.
    i2 : str, optional
        Column name for I2 statistic in input data.
    snpr2 : str, optional
        Column name for SNP R2 in input data.
    phet : str, optional
        Column name for p-heterogeneity in input data.
    dof : str, optional
        Column name for degrees of freedom in input data.
    ncase : str or int, optional
        Column name or constant value for case count.
    ncontrol : str or int, optional
        Column name or constant value for control count.
    neff : str or int, optional
        Column name or constant value for effective sample size.
    direction : str, optional
        Column name for meta-analysis effect-direction strings, where each character
        ("+", "-", "0") represents the direction of effect for one cohort (e.g., "++-0+").
    status : str, optional
        Column name for status in input data.
    study : str, optional, default="Study_1"
        Column name for study ID in input data.
    trait : str, optional, default="Trait_1"
        Column name for trait in input data.
    build : str, optional
        Genome build version (e.g., '19' and '38').
    species : str, default="homo sapiens"
        species
    other : list, optional
        Additional columns in the raw file to load.
    chrom_pat : str, optional
        Regex pattern to filter chromosomes like'chrX'
    snpid_pat : str, optional
        Regex pattern to filter SNPs based on snpid like'chrX:'.
    verbose : bool, default: False
        Enable verbose output.
    readargs : dict, optional
        Additional arguments for reading files using pd.read_csv() like `nrows`, `comment`. 
        Example, {"nrows": 1000} means to load first 1000 rows.

    Returns
    -------
    pandas.DataFrame
        Formatted summary statistics with standardized column names.
    dict
        Updated readargs dictionary.

    Raises
    ------
    ValueError
        If input is not a path or DataFrame, or if columns are missing.

    Less used parameters
    -------------------------
    exclude : list, optional
        Columns to exclude explicitly. Columns should be passed as GWASLab built-in HEADER keywords in uppercase like BETA, DIRECTION. Not original headers.
    include : list, optional
        Columns to include explicitly. Columns should be passed as GWASLab built-in HEADER keywords in uppercase like BETA, DIRECTION. Not original headers.
    """
    # Step 1: Initialize parameters
    readargs, log, other, exclude, include = _initialize_preformat_parameters(
        readargs=readargs, kwreadargs=kwreadargs, log=log, other=other, exclude=exclude, include=include
    )
    
    # Step 2: Handle parquet format (must be done before other processing)
    if tab_fmt == "parquet":
        if isinstance(sumstats, str):
            log.write("Start to load data from parquet file....", verbose=verbose)
            log.write(" -path: {}".format(sumstats), verbose=verbose)
            sumstats = pd.read_parquet(sumstats, **readargs)
            log.write("Finished loading parquet file into pd.DataFrame....", verbose=verbose)
        else:
            raise ValueError("Please input a path for parquet file.")
    
    # Step 3: Load format configuration
    meta_data, rename_dictionary, readargs, other = _load_format_config(
        fmt=fmt, readargs=readargs, other=other, log=log, verbose=verbose
    )
    
    # Initialize data structures
    usecols = []
    dtype_dictionary = {}
    
    # Step 4: Check path and header, get available columns
    inpath, inpath_chr_list, inpath_chr_num_list, format_cols, raw_cols, usecols, dtype_dictionary = _check_path_and_header(
        sumstats, fmt, meta_data, readargs, usecols, dtype_dictionary, rename_dictionary, log, verbose
    )
    
    # Step 5: Build column mappings from user parameters
    column_params = {
        'snpid': snpid, 'rsid': rsid, 'chrom': chrom, 'pos': pos, 'ea': ea, 'nea': nea,
        'ref': ref, 'alt': alt, 'eaf': eaf, 'neaf': neaf, 'maf': maf, 'n': n,
        'ncase': ncase, 'ncontrol': ncontrol, 'neff': neff, 'beta': beta,
        'beta_95L': beta_95L, 'beta_95U': beta_95U, 'se': se, 'chisq': chisq,
        'z': z, 'q': q, 'p': p, 't': t, 'f': f, 'mlog10p': mlog10p, 'test': test,
        'info': info, 'OR': OR, 'OR_95L': OR_95L, 'OR_95U': OR_95U, 'HR': HR,
        'HR_95L': HR_95L, 'HR_95U': HR_95U, 'i2': i2, 'snpr2': snpr2, 'phet': phet,
        'dof': dof, 'direction': direction, 'status': status, 'other': other
    }
    _build_column_mappings(
        column_params=column_params, rename_dictionary=rename_dictionary, 
        usecols=usecols, dtype_dictionary=dtype_dictionary, log=log, verbose=verbose
    )
    
    # Step 6: Handle VCF format special case
    vcf_usecols = None
    if fmt == "vcf":
        vcf_usecols = usecols.copy()
        usecols = meta_data["format_fixed"]
        if study is not None:
            usecols = usecols + [study]
        else:
            study = raw_cols[9]
            usecols = usecols + [study]
    
    # Step 7: Apply include/exclude filters
    usecols = _apply_column_filters(
        include=include, exclude=exclude, usecols=usecols, 
        rename_dictionary=rename_dictionary, log=log, verbose=verbose
    )
    
    # Step 8: Load data
    sumstats = _load_sumstats_data(
        sumstats=sumstats, inpath=inpath, inpath_chr_list=inpath_chr_list, inpath_chr_num_list=inpath_chr_num_list,
        usecols=usecols, dtype_dictionary=dtype_dictionary, readargs=readargs, rename_dictionary=rename_dictionary,
        chrom_pat=chrom_pat, snpid_pat=snpid_pat, log=log, verbose=verbose
    )
    
    # Step 9: Post-process data (rename, transform, validate)
    sumstats = _postprocess_sumstats(
        sumstats=sumstats, fmt=fmt, format_cols=format_cols, study=study, vcf_usecols=vcf_usecols, usecols=usecols,
        rename_dictionary=rename_dictionary, n=n, ncase=ncase, ncontrol=ncontrol, 
        build=build, status=status, neaf=neaf, log=log, verbose=verbose
    )
    
    log.write("Finished loading data successfully!", verbose=verbose)
    return sumstats


#### helper #######################################################################
def _isfile_casesensitive(path):
    if not os.path.isfile(path): 
        return False   # exit early
    directory, filename = os.path.split(path)
    return filename in os.listdir(directory)

def _get_readargs_header(inpath,readargs):
    if "vcf.gz" in inpath:
        with gzip.open(inpath,'r') as file:      
            skip=0
            for line in file:        
                if line.decode('utf-8').startswith('##'):
                    skip+=1
                else:
                    readargs["skiprows"]=skip
                    readargs["sep"]="\t"
                    break
    readargs_header = readargs.copy()
    readargs_header["nrows"]=1
    readargs_header["dtype"]="string"
    return readargs_header

def _get_skip_rows(inpath):
    if "vcf.gz" in inpath:
        with gzip.open(inpath,'r') as file:      
            skip=0
            for line in file:        
                if line.decode('utf-8').startswith('##'):
                    skip+=1
                else:
                    return skip
    else:
        return 0

def _parse_vcf_study(sumstats,format_cols,study,vcf_usecols,log,verbose=True):
    log.write(" -Parsing based on FORMAT: ", format_cols,verbose=verbose)
    log.write(" -Parsing vcf study : ", study,verbose=verbose)
    sumstats[format_cols] = sumstats[study].str.split(":",expand=True).values
    sumstats = sumstats.drop(["FORMAT",study],axis=1)
    sumstats = sumstats[ vcf_usecols]
    gc.collect()
    return sumstats

def _print_format_info(fmt,meta_data, rename_dictionary, verbose, log,output=False, skip_meta_records=None):
    log.write(" -"+fmt+" format meta info:",verbose=verbose)   
    if skip_meta_records is None:
        skip_meta_records =[]
    for key,value in meta_data.items():
        if key in skip_meta_records:
            continue
        if value is None:
            continue
        if type(value) is str:
            if "\n" in value:
                value_first_line=value.split("\n")[0]
                log.write("  -",key," : "+value_first_line.strip()+"...",verbose=verbose)
            elif value==" ":
                log.write('  -',key,' : \\s ',verbose=verbose)     
            elif value=="\t":
                log.write('  -',key,' : \\t',verbose=verbose)    
            else:
                log.write("  -",key," : "+value.strip(),verbose=verbose)  
        elif type(value) is list:
            log.write("  -",key," : "+','.join(value),verbose=verbose)  
        else:
            log.write("  -",key," : ",value,verbose=verbose)  
    keys=[]
    values=[]
    for key,value in rename_dictionary.items():
        keys.append(key)
        values.append(value)
    if fmt!="gwaslab":
        if output == False:
            if fmt!="auto":
                log.write(" -"+fmt+" to gwaslab format dictionary:",verbose=verbose)  
                log.write("  - "+fmt+" keys:",",".join(keys),verbose=verbose) 
                log.write("  - gwaslab values:",",".join(values),verbose=verbose) 
            else:
                log.write("  - Auto-detection mode. Note: auto-detection assumes A1=EA; Alt=EA and Frq=EAF...",verbose=verbose)
                log.write("  - Header conversion source: https://github.com/Cloufield/formatbook/blob/main/formats/auto.json",verbose=verbose)  
        else:
            log.write(" -gwaslab to "+fmt+" format dictionary:",verbose=verbose)  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            log.write("  - gwaslab keys:",  ','.join(keys),verbose=verbose) 
            log.write("  - "+fmt+" values:"  , ','.join(values),verbose=verbose)       

def _process_neaf(sumstats,log,verbose):
    log.write(" -NEAF is specified...",verbose=verbose) 
    pre_number=len(sumstats)
    log.write(" -Checking if 0<= NEAF <=1 ...",verbose=verbose) 
    if "NEAF" in sumstats.columns:
        sumstats["NEAF"] = pd.to_numeric(sumstats["NEAF"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["NEAF"]>=0) & (sumstats["NEAF"]<=1),:].copy()
        sumstats.loc[:, "EAF"] = 1 - sumstats["NEAF"]
        sumstats.drop(columns=["NEAF"], inplace=True)
    else:
        sumstats["EAF"] = pd.to_numeric(sumstats["EAF"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["EAF"]>=0) & (sumstats["EAF"]<=1),:].copy()
        sumstats.loc[:, "EAF"] = 1 - sumstats["EAF"]
    log.write(" -Converted NEAF to EAF.",verbose=verbose) 
    after_number=len(sumstats)
    log.write(" -Removed "+str(pre_number - after_number)+" variants with bad NEAF.",verbose=verbose) 
    return sumstats

def _process_allele(sumstats,
                   log,
                   rename_dictionary,
                   verbose):
    
    if "EA" in sumstats.columns:

        if "REF" in sumstats.columns and "ALT" in sumstats.columns:

            if "NEA" not in sumstats.columns:
                log.write(" NEA not available: assigning REF to NEA...",verbose=verbose) 
                sumstats["NEA"]=sumstats["REF"]    
            
            log.write(" -EA,REF and ALT columns are available: assigning NEA...",verbose=verbose) 
            ea_alt = sumstats["EA"]==sumstats["ALT"]
            
            log.write(" -For variants with EA == ALT : assigning REF to NEA ...",verbose=verbose) 
            sumstats.loc[ea_alt,"NEA"] = sumstats.loc[ea_alt,"REF"]
            
            ea_not_alt = sumstats["EA"]!=sumstats["ALT"]
            log.write(" -For variants with EA != ALT : assigning ALT to NEA ...",verbose=verbose) 
            sumstats.loc[ea_not_alt,"NEA"] = sumstats.loc[ea_not_alt,"ALT"]

            #sumstats = sumstats.drop(labels=["REF","ALT"],axis=1)
            sumstats["REF"]=sumstats["REF"].astype("category") 
            sumstats["ALT"]=sumstats["ALT"].astype("category") 
            
        sumstats["EA"]=sumstats["EA"].astype("category")     
    if "NEA" in sumstats.columns:
        sumstats["NEA"]=sumstats["NEA"].astype("category")  

    if "NEA" not in sumstats.columns and "EA" not in sumstats.columns:
        if "REF" in sumstats.columns and "ALT" in sumstats.columns:
            log.write(" -Converting REF and ALT to NEA and EA ...")
            sumstats["REF"]=sumstats["REF"].astype("category") 
            sumstats["ALT"]=sumstats["ALT"].astype("category") 
            sumstats["NEA"]=sumstats["REF"].astype("category") 
            sumstats["EA"]=sumstats["ALT"].astype("category") 
            sumstats = sumstats.drop(columns=["ALT","REF"])
    return sumstats

def _process_status(sumstats,build,status, log,verbose):
    if status is None:
        log.write(" -Initiating a status column: STATUS ...",verbose=verbose)
        #sumstats["STATUS"] = int(build)*(10**5) +99999
        build = _process_build(build,log,verbose)
        # Create integer status code: build (2 digits) + 99999 (5 digits) = 7 digits
        sumstats["STATUS"] = int(build + "99999")

    # Convert STATUS to integer if it's string (for backward compatibility)
    if sumstats["STATUS"].dtype == 'object' or sumstats["STATUS"].dtype.name == 'category':
        sumstats["STATUS"] = sumstats["STATUS"].astype(str).astype(int)
    elif sumstats["STATUS"].dtype not in ['int64', 'Int64', 'int32', 'Int32']:
        sumstats["STATUS"] = sumstats["STATUS"].astype(int)
    
    return sumstats


def _load_single_chr(*, inpath=None, usecols=None, dtype_dictionary=None, readargs=None, rename_dictionary=None, chrom_pat=None, log=None, verbose=False):
    explicit = {"usecols","dtype_dictionary","chunksize","iterator"}
    readargs = {k: v for k, v in readargs.items() if k not in explicit}
    sumstats_iter = pd.read_table(inpath,
                usecols=set(usecols),
                dtype=dtype_dictionary, 
                iterator=True, 
                chunksize=500000,
                **readargs)
    # get chr 
    for k,v in rename_dictionary.items():
        if v=="CHR":
            if k in usecols:
                log.write(" -Columns used to filter variants: {}".format(k),verbose=verbose)
                chunk_chrom = k
                break

    log.write(" -Loading only variants on chromosome with pattern : {} ...".format(chrom_pat),verbose=verbose)
    sumstats_filtered = pd.concat([chunk[chunk[chunk_chrom].str.match(chrom_pat, case=False,na=False) ] for chunk in sumstats_iter])
    log.write(" -Loaded {} variants on chromosome with pattern :{} ...".format(len(sumstats_filtered), chrom_pat),verbose=verbose)
    return sumstats_filtered

def _load_variants_with_pattern(*, inpath=None, usecols=None, dtype_dictionary=None, readargs=None, rename_dictionary=None, snpid_pat=None, log=None, verbose=False):
    explicit = {"usecols","dtype_dictionary","chunksize","iterator"}
    readargs = {k: v for k, v in readargs.items() if k not in explicit}
    sumstats_iter = pd.read_table(inpath,
                usecols=set(usecols),
                dtype=dtype_dictionary, 
                iterator=True, 
                chunksize=500000,
                **readargs)
    # get chr 
    for k,v in rename_dictionary.items():
        if v=="SNPID":
            if k in usecols:
                log.write(" -Columns used to filter variants: {}".format(k),verbose=verbose)
                chunk_snpid = k
                break

    log.write(" -Loading only variants with pattern :  {} ...".format(snpid_pat),verbose=verbose)
    sumstats_filtered = pd.concat([chunk[chunk[chunk_snpid].str.match(snpid_pat, case=False,na=False) ] for chunk in sumstats_iter])
    log.write(" -Loaded {} variants with pattern : {} ...".format(len(sumstats_filtered), snpid_pat),verbose=verbose)
    return sumstats_filtered


def _check_path_and_header(sumstats=None, 
                          fmt=None, 
                          meta_data=None, 
                          readargs=None, 
                          usecols=None, 
                          dtype_dictionary=None, 
                          rename_dictionary=None, 
                          log=None, 
                          verbose=None):
    

    if isinstance(sumstats, str):
        ## loading data from path #################################################
        inpath = sumstats
        
        try:
            format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list = _process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose)
        
        except (FileNotFoundError, IndexError):
            log.warning("Loading {} failed...Tesing if compressed/uncompressed...".format(inpath),verbose=verbose)
            try:
                if inpath[-3:]==".gz":
                    inpath = inpath[:-3]
                    log.write(" -Trying to load {}...".format(inpath),verbose=verbose)
                    format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list = _process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose)
                else:
                    inpath = inpath+".gz"
                    log.write(" -Trying to load {}...".format(inpath),verbose=verbose)
                    format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list = _process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose)
            except:
                raise ValueError("Please input a valid path, and make sure the separator is correct and the columns you specified are in the file.") 

        ###################################################################################### 
    elif isinstance(sumstats, pd.DataFrame):
        inpath = None
        format_cols = None
        inpath_chr_list = None
        inpath_chr_num_list = None
        ## loading data from dataframe
        raw_cols = sumstats.columns

    ################################################
    for key,value in rename_dictionary.items():
        # check available keys  key->raw header
        # usecols : a list of raw headers to load from file/DataFrame 
        if key in raw_cols:
            usecols.append(key)
        if value in ["EA","NEA"]:
            dtype_dictionary[key]="category"
        if value in ["STATUS"]:
            dtype_dictionary[key]="string"     
        if value in ["CHR"]:
            dtype_dictionary[key]="string"  

    return inpath, inpath_chr_list, inpath_chr_num_list, format_cols, raw_cols, usecols, dtype_dictionary

def _process_inpath_and_load_header(inpath, fmt, meta_data,  readargs, log, verbose):
    
    format_cols = None
    inpath_chr_list = None
    inpath_chr_num_list = None

    if "@" in inpath:
        log.write(" -Detected @ in path: load sumstats by each chromosome...",verbose=verbose)
        inpath_chr_list=[]
        inpath_chr_num_list=[]
        
        # create a regex pattern for matching
        pat = os.path.basename(inpath).replace("@", r"(\w+)")
        
        # get dir
        dirname = os.path.dirname(inpath)
        
        # all files in the directory
        files = os.listdir(dirname)
        
        files.sort()

        for file in files:
            # match
            result = re.match(pat, file)
            if result:
                # get chr
                chr_matched = str(result.group(1))
                inpath_chr_num_list.append(chr_matched)
                inpath_chr_list.append(inpath.replace("@",str(chr_matched))  )
        
        log.write(" -Chromosomes detected:",",".join(inpath_chr_num_list),verbose=verbose)

        #if inpath_chr_list is empty-> IndexError
        readargs_header = _get_readargs_header(inpath = inpath_chr_list[0], readargs = readargs)
        row_one = pd.read_table(inpath_chr_list[0],**readargs_header)
        # columns in the sumstats
        raw_cols = row_one.columns
    else:
    ##### loading data from tabular file#################################################
    #if file not found, FileNotFoundError
        readargs_header = _get_readargs_header(inpath = inpath, readargs = readargs)
        row_one = pd.read_table(inpath,**readargs_header)
        raw_cols = row_one.columns

    if fmt=="vcf":
        # expanded
        format_cols = list(row_one["FORMAT"].str.split(":"))[0]
        # fixed + study1 + expanded
        raw_cols = meta_data["format_fixed"] + [raw_cols[9]] + format_cols

    return format_cols, raw_cols, inpath_chr_list, inpath_chr_num_list
