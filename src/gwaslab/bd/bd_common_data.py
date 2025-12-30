from typing import Optional, Dict, List, Union, Tuple
import os
import json
import tarfile
import requests
import pandas as pd
from os import path
from pathlib import Path
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_download import download_ref
from gwaslab.bd.bd_download import check_and_download
from gwaslab.bd.bd_download import update_formatbook
from gwaslab.bd.bd_config import options
from gwaslab.bd.bd_sex_chromosomes import Chromosomes
from gwaslab.qc.qc_build import _process_build

# Species-specific NCBI RefSeq accession ID mappings
# Format: {(species, build_code): {chromosome: ncbi_accession_id, ...}}
# Loaded from JSON file: src/gwaslab/data/chromosomes/chromosomes_nc.json
_NCBI_ACCESSION_IDS = None

def _load_ncbi_accession_ids() -> Dict[Tuple[str, str], Dict[str, str]]:
    """Load NCBI accession IDs from JSON file."""
    global _NCBI_ACCESSION_IDS
    if _NCBI_ACCESSION_IDS is not None:
        return _NCBI_ACCESSION_IDS
    
    # Get the path to the JSON file
    json_path = Path(__file__).parent.parent / "data" / "chromosomes" / "chromosomes_nc.json"
    
    try:
        with open(json_path, 'r') as f:
            json_data = json.load(f)
        
        # Convert JSON format {species: {build: {chromosome: accession}}}
        # to code format {(species, build): {chromosome: accession}}
        _NCBI_ACCESSION_IDS = {}
        for species, builds in json_data.items():
            for build, chromosomes in builds.items():
                key = (species.lower(), str(build))
                _NCBI_ACCESSION_IDS[key] = chromosomes.copy()
        
        return _NCBI_ACCESSION_IDS
    except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
        # Fallback to empty dict if file not found or invalid
        _NCBI_ACCESSION_IDS = {}
        return _NCBI_ACCESSION_IDS

def _get_ncbi_accession_mapping(species: Optional[str], build: str, log: Optional[Log] = None, verbose: bool = True) -> Optional[Dict[str, str]]:
    """
    Get NCBI RefSeq accession ID mapping for a species and build.
    
    Parameters
    ----------
    species : str
        Species name (case-insensitive)
    build : str
        Genome build identifier (will be processed to normalized build code)
    log : Log, optional
        Logging object
    verbose : bool, default True
        Whether to print warnings
    
    Returns
    -------
    dict or None
        Dictionary mapping chromosome identifiers to NCBI accession IDs,
        or None if not available
    """
    # Load NCBI accession IDs from JSON if not already loaded
    ncbi_ids = _load_ncbi_accession_ids()
    
    if species is None:
        species = "homo sapiens"
    
    species_lower = species.lower()
    
    # Process build to get normalized build code
    try:
        processed_build = _process_build(build, log=log or Log(), verbose=False, species=species)
    except (ValueError, KeyError):
        # If build processing fails, try using the original build string
        processed_build = str(build)
    
    # Check if mapping exists
    key = (species_lower, processed_build)
    if key in ncbi_ids:
        return ncbi_ids[key].copy()
    
    # Check for aliases (try both directions)
    species_aliases = {
        "homo sapiens": "human",
        "human": "homo sapiens",
        "mus musculus": "mouse",
        "mouse": "mus musculus",
        "rattus norvegicus": "rat",
        "rat": "rattus norvegicus",
        "gallus gallus": "chicken",
        "chicken": "gallus gallus",
        "danio rerio": "zebrafish",
        "zebrafish": "danio rerio",
        "drosophila melanogaster": "fruit fly",
        "fruit fly": "drosophila melanogaster",
        "sus scrofa": "pig",
        "pig": "sus scrofa",
        "bos taurus": "cattle",
        "cattle": "bos taurus",
        "cow": "bos taurus",
        "canis lupus familiaris": "dog",
        "dog": "canis lupus familiaris",
        "equus caballus": "horse",
        "horse": "equus caballus",
        "oryza sativa": "rice",
        "rice": "oryza sativa",
        "arabidopsis thaliana": "arabidopsis",
        "arabidopsis": "arabidopsis thaliana",
    }
    
    # Try alias
    alias = species_aliases.get(species_lower)
    if alias:
        key = (alias.lower(), processed_build)
        if key in ncbi_ids:
            return ncbi_ids[key].copy()
    
    # Not found - return None (caller should handle this)
    return None

#hard-coded data
def get_chr_to_NC(build: str, inverse: bool = False, species: str = "homo sapiens", log: Optional[Log] = None, verbose: bool = True) -> Dict[str, str]:
    """
    Create a dictionary mapping chromosome identifiers to NCBI RefSeq accession IDs.
    
    Parameters
    ----------
    build : str
        Genome build version (e.g., '19', '38', 'mm10', 'rn6'). Will be normalized based on species.
    inverse : bool, default False
        If True, return NCBI ID to chromosome mapping
    species : str, default "homo sapiens"
        Species name (case-insensitive). Currently only human mappings are fully supported.
    log : Log, optional
        Logging object for warnings
    verbose : bool, default True
        Whether to print warnings
    
    Returns
    -------
    dict
        Dictionary mapping chromosome identifiers (string) to NCBI accession IDs (string).
        Returns empty dict for unsupported species/build combinations.
    """
    if log is None:
        from gwaslab.info.g_Log import Log
        log = Log()
    
    # Get NCBI accession mapping
    dic = _get_ncbi_accession_mapping(species, build, log=log, verbose=verbose)
    
    if dic is None:
        # No mapping available for this species/build
        dic = {}
    
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic

def get_NC_to_chr(build: str, species: str = "homo sapiens", log: Optional[Log] = None, verbose: bool = True) -> Dict[str, str]:
    """
    Create a dictionary mapping NCBI RefSeq accession IDs to chromosome identifiers.
    
    Parameters
    ----------
    build : str
        Genome build version (e.g., '19', '38', 'mm10', 'rn6'). Will be normalized based on species.
    species : str, default "homo sapiens"
        Species name (case-insensitive). Currently only human mappings are fully supported.
    log : Log, optional
        Logging object for warnings
    verbose : bool, default True
        Whether to print warnings
    
    Returns
    -------
    dict
        Dictionary mapping NCBI accession IDs (string) to chromosome identifiers (string)
    """
    return get_chr_to_NC(build=build, inverse=True, species=species, log=log, verbose=verbose)


def get_number_to_NC(build: str, inverse: bool = False, species: str = "homo sapiens", log: Optional[Log] = None, verbose: bool = True) -> Dict[Union[int, str], str]:
    """
    Create a dictionary mapping chromosome numbers (int) to NCBI RefSeq accession IDs (string).
    
    Parameters
    ----------
    build : str
        Genome build version (e.g., '19', '38', 'mm10', 'rn6'). Will be normalized based on species.
    inverse : bool, default False
        If True, return NCBI ID to chromosome number mapping
    species : str, default "homo sapiens"
        Species name (case-insensitive). Currently only human mappings are fully supported.
    log : Log, optional
        Logging object for warnings
    verbose : bool, default True
        Whether to print warnings
    
    Returns
    -------
    dict
        Dictionary mapping chromosome numbers (int) to NCBI accession IDs (string).
        Returns empty dict for unsupported species/build combinations.
    """
    if log is None:
        from gwaslab.info.g_Log import Log
        log = Log()
    
    # Get string-based mapping first
    chr_to_nc = get_chr_to_NC(build=build, inverse=False, species=species, log=log, verbose=verbose)
    
    if not chr_to_nc:
        # Return empty dict if no mapping available
        return {}
    
    # Convert to number-based mapping
    # Get chromosome to number mapping for this species
    chromosomes_obj = Chromosomes(species=species)
    chr_to_num = chromosomes_obj.get_chr_to_number_dict(out_chr=False, xymt_num=[23, 24, 25])
    
    # Build number to NC mapping
    dic = {}
    for chrom_str, nc_id in chr_to_nc.items():
        # Get numeric value for this chromosome
        if chrom_str in chr_to_num:
            chrom_num = chr_to_num[chrom_str]
            dic[chrom_num] = nc_id
    
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic


def get_NC_to_number(build: str, species: str = "homo sapiens", log: Optional[Log] = None, verbose: bool = True) -> Dict[str, Union[int, str]]:
    """
    Create a dictionary mapping NCBI RefSeq accession IDs to chromosome numbers (int).
    
    Parameters
    ----------
    build : str
        Genome build version (e.g., '19', '38', 'mm10', 'rn6'). Will be normalized based on species.
    species : str, default "homo sapiens"
        Species name (case-insensitive). Currently only human mappings are fully supported.
    log : Log, optional
        Logging object for warnings
    verbose : bool, default True
        Whether to print warnings
    
    Returns
    -------
    dict
        Dictionary mapping NCBI accession IDs to chromosome numbers (int)
    """
    return get_number_to_NC(build=build, inverse=True, species=species, log=log, verbose=verbose)

def get_chr_list(add_number: bool = False, n: int = 25, only_number: bool = False, species: str = "homo sapiens") -> List[Union[str, int]]:
    """
    Generate a list of chromosome identifiers.
    
    Parameters:
    -----------
    add_number : bool, default=False
        If True, include both string and numeric representations
    n : int, default=25
        Maximum chromosome number to include (deprecated, now uses species-specific chromosomes)
    only_number : bool, default=False
        If True, return only numeric chromosome numbers
    species : str, default="homo sapiens"
        Species name (case-insensitive). If provided, uses species-specific chromosomes.
        Otherwise falls back to legacy behavior with n parameter.
    
    Returns:
    --------
    list
        List of chromosome identifiers in string format by default, 
        or numeric format if specified
    """
    # Use Chromosomes class if species is provided (or default to human)
    if species is not None:
        chromosomes_obj = Chromosomes(species=species)
        return chromosomes_obj.get_chr_list(add_number=add_number, only_number=only_number)
    
    # Legacy behavior: fallback to old implementation if species is explicitly None
    chrom_list = [str(i) for i in range(1, n+1)] + ["X", "Y", "M", "MT"]
    
    if add_number:
        chrom_list = [str(i) for i in range(1, n+1)] + ["X", "Y", "M", "MT"] + [i for i in range(1, n+1)]

    if only_number:
        chrom_list = [i for i in range(1, n+1)]
    return chrom_list
def get_chr_to_number(out_chr: bool = False, xymt: Optional[List[str]] = None, xymt_num: Optional[List[int]] = None, species: str = "homo sapiens", max_chr: int = 200) -> Dict[str, Union[int, str]]:
    """
    Create a dictionary mapping chromosome identifiers to numeric representations.
    
    Parameters:
    -----------
    out_chr : bool, default=False
        If True, returns dictionary with string keys and values
    xymt : list, default=["X","Y","MT"]
        List of non-numeric chromosome identifiers (deprecated, now uses species-specific chromosomes)
    xymt_num : list, default=[23,24,25]
        Corresponding numeric values for xymt (deprecated, now uses species-specific mappings)
    species : str, default="homo sapiens"
        Species name (case-insensitive). If provided, uses species-specific chromosomes.
        Otherwise falls back to legacy behavior with xymt/xymt_num parameters.
    max_chr : int, default=200
        Maximum chromosome number to include in dictionary
    
    Returns:
    --------
    dict
        Dictionary mapping chromosome identifiers to numeric values or strings
        depending on the out_chr parameter
    """
    if xymt is None:
        xymt = ["X","Y","MT"]
    if xymt_num is None:
        xymt_num = [23,24,25]
    # Use Chromosomes class if species is provided (or default to human)
    if species is not None:
        chromosomes_obj = Chromosomes(species=species)
        return chromosomes_obj.get_chr_to_number_dict(out_chr=out_chr, xymt_num=xymt_num, max_chr=max_chr)
    
    # Legacy behavior: fallback to old implementation if species is explicitly None
    if out_chr:
        dic = {str(i): str(i) for i in range(1, max_chr + 1)}
        if len(xymt) > 0:
            dic[xymt[0]] = str(xymt_num[0])
        if len(xymt) > 1:
            dic[xymt[1]] = str(xymt_num[1])
        if len(xymt) > 2:
            dic[xymt[2]] = str(xymt_num[2])
    else:
        dic = {str(i): i for i in range(1, max_chr + 1)}
        if len(xymt) > 0:
            dic[xymt[0]] = xymt_num[0]
        if len(xymt) > 1:
            dic[xymt[1]] = xymt_num[1]
        if len(xymt) > 2:
            dic[xymt[2]] = xymt_num[2]
            dic["M"] = xymt_num[2]
    return dic
def get_number_to_chr(in_chr: bool = False, xymt: Optional[List[str]] = None, xymt_num: Optional[List[int]] = None, prefix: str = "", species: str = "homo sapiens", max_chr: int = 200) -> Dict[Union[int, str], str]:
    """
    Create a dictionary mapping chromosome numbers to string representations.
    
    Parameters:
    -----------
    in_chr : bool, default=False
        If True, returns dictionary with string keys and values
    xymt : list, default=["X","Y","MT"]
        List of non-numeric chromosome identifiers (deprecated, now uses species-specific chromosomes)
    xymt_num : list, default=[23,24,25]
        Corresponding numeric values for xymt (deprecated, now uses species-specific mappings)
    prefix : str, default=""
        Optional prefix for chromosome identifiers
    species : str, default="homo sapiens"
        Species name (case-insensitive). If provided, uses species-specific chromosomes.
        Otherwise falls back to legacy behavior with xymt/xymt_num parameters.
    max_chr : int, default=200
        Maximum chromosome number to include in dictionary
    
    Returns:
    --------
    dict
        Dictionary mapping chromosome numbers to string representations
    """
    if xymt is None:
        xymt = ["X","Y","MT"]
    if xymt_num is None:
        xymt_num = [23,24,25]
    # Use Chromosomes class if species is provided (or default to human)
    if species is not None:
        chromosomes_obj = Chromosomes(species=species)
        return chromosomes_obj.get_number_to_chr_dict(in_chr=in_chr, xymt_num=xymt_num, prefix=prefix, max_chr=max_chr)
    
    # Legacy behavior: fallback to old implementation if species is explicitly None
    if in_chr:
        dic = {str(i): prefix + str(i) for i in range(1, max_chr + 1)}
        if len(xymt) > 0 and len(xymt_num) > 0:
            dic[str(xymt_num[0])] = prefix + xymt[0]
        if len(xymt) > 1 and len(xymt_num) > 1:
            dic[str(xymt_num[1])] = prefix + xymt[1]
        if len(xymt) > 2 and len(xymt_num) > 2:
            dic[str(xymt_num[2])] = prefix + xymt[2]
    else:
        dic = {i: prefix + str(i) for i in range(1, max_chr + 1)}
        if len(xymt) > 0 and len(xymt_num) > 0:
            dic[xymt_num[0]] = prefix + xymt[0]
        if len(xymt) > 1 and len(xymt_num) > 1:
            dic[xymt_num[1]] = prefix + xymt[1]
        if len(xymt) > 2 and len(xymt_num) > 2:
            dic[xymt_num[2]] = prefix + xymt[2]
    return dic
# reading from files    
###################################################################################################################    
def get_high_ld(build: str = "19") -> str:
    """
    Get the path to the high LD region file for the specified genome build.
    
    Parameters:
    build (str): Genome build version ('19' or '38') indicating which reference genome to use
    
    Returns:
    str: Path to the high LD region BED file for the specified genome build
    """
    if build=="19":
        #data_path =  path.dirname(__file__) + '/data/high_ld/high_ld_hla_hg19.bed.gz'
        data_path = path.join( Path(__file__).parents[1], "data","high_ld","high_ld_hla_hg19.bed.gz")
    elif build=="38":
        #data_path =  path.dirname(__file__) + '/data/high_ld/high_ld_hla_hg38.bed.gz'
        data_path = path.join( Path(__file__).parents[1], "data","high_ld","high_ld_hla_hg38.bed.gz")
    return data_path
def get_format_dict(fmt: str, inverse: bool = False) -> tuple:
    """
    Retrieve format dictionary and metadata for a specified format.
    
    Parameters:
    fmt (str): Format name to look up in the format book
    inverse (bool): If True, return inverted dictionary with value-key mapping
    
    Returns:
    tuple: 
    - dict: Metadata associated with the format
    - dict: Format dictionary mapping fields, or inverted mapping if specified
    """
    #data_path =  path.dirname(__file__) + '/data/formatbook.json'
    data_path = options.paths["formatbook"]
    if not path.exists(data_path):
        update_formatbook()
    dicts = json.load(open(data_path))
    dic_meta = dicts[fmt]["meta_data"]
    dic_dict = dicts[fmt]["format_dict"]
    if inverse is True:
        inv_dic = {v: k for k, v in dic_dict.items()}
        return dic_meta,inv_dic
    return dic_meta, dic_dict
def get_formats_list() -> List[str]:
    """
    Retrieve a list of available format names from the format book.
    
    Returns:
    list: Format names available in the format book
    """
    #data_path =  path.dirname(__file__) + '/data/formatbook.json'
    data_path = options.paths["formatbook"]
    dicts = json.load(open(data_path))
    format_list = list(dicts.keys())
    return format_list
# Module-level cache for recombination rate data
_RECOMBINATION_RATE_CACHE = {}

def get_recombination_rate(chrom: Union[str, int], build: str = "19") -> pd.DataFrame:
    """
    Retrieve recombination rate data for a specific chromosome and genome build.
    
    Results are cached for faster subsequent access.
    
    Parameters:
    chrom (str or int): Chromosome number to retrieve recombination data for
    build (str): Genome build version ('19' or '38') to use
    
    Returns:
    pandas.DataFrame: Recombination rate data with columns 'Rate(cM/Mb)' and 'Position(bp)'
                     Returns empty DataFrame if build is not supported or data is unavailable
    """
    # Create cache key
    cache_key = f"{build}_{chrom}"
    
    # Check cache first
    if cache_key in _RECOMBINATION_RATE_CACHE:
        return _RECOMBINATION_RATE_CACHE[cache_key].copy()
    
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110106_recombination_hotspots/
    if build=="19":
        data_path =  options.paths["data_directory"] + 'recombination/hg19/genetic_map_GRCh37_chr'+str(chrom)+'.txt.gz'
        if path.exists(data_path):
            recombination_rate = pd.read_csv(data_path,sep="\t")
        else:
            if not path.exists(options.paths["data_directory"] + 'recombination/hg19/'):
                os.makedirs(options.paths["data_directory"] + 'recombination/hg19/')
            download_ref("recombination_hg19",directory = options.paths["data_directory"] + 'recombination/hg19/')
            file = tarfile.open(options.paths["data_directory"]+'recombination/hg19/recombination_hg19.tar.gz')
            file.extractall(options.paths["data_directory"]+'recombination/hg19/')
            file.close()
            recombination_rate = pd.read_csv(data_path,sep="\t")
    elif build=="38":
        data_path =  options.paths["data_directory"] + 'recombination/hg38/genetic_map_GRCh38_chr'+str(chrom)+'.txt.gz'
        if path.exists(data_path):
            recombination_rate = pd.read_csv(data_path,sep="\t")
        else:
            if not path.exists( options.paths["data_directory"] + 'recombination/hg38/'):
                os.makedirs( options.paths["data_directory"] + 'recombination/hg38/')
            download_ref("recombination_hg38",directory = options.paths["data_directory"]+'recombination/hg38/')
            file = tarfile.open(options.paths["data_directory"]+'recombination/hg38/recombination_hg38.tar.gz')
            file.extractall(options.paths["data_directory"]+'recombination/hg38/')
            file.close()
            recombination_rate = pd.read_csv(data_path,sep="\t")
    else:
        recombination_rate = pd.DataFrame(columns=["Rate(cM/Mb)","Position(bp)"])
    
    # Cache the result for future use
    _RECOMBINATION_RATE_CACHE[cache_key] = recombination_rate.copy()
    
    return recombination_rate
####################################################################################################################
def get_chain(from_build: str = "19", to_build: str = "38") -> str:    
    """
    Get the path to a chain file for liftover between genome builds.
    
    First checks for built-in chain files in the package data directory.
    If not found, falls back to downloading from reference.json.
    
    Parameters
    ----------
    from_build : str
        Source genome build (e.g., "19" for hg19/GRCh37)
    to_build : str
        Target genome build (e.g., "38" for hg38/GRCh38)
    
    Returns
    -------
    str
        Path to the chain file
    """
    # Map build numbers to chain file names
    chain_files = {
        ("19", "38"): "hg19ToHg38.over.chain.gz",
        ("38", "19"): "hg38ToHg19.over.chain.gz",
    }
    
    # Check for built-in chain file first
    chain_filename = chain_files.get((str(from_build), str(to_build)))
    if chain_filename:
        builtin_chain_path = path.join(Path(__file__).parents[1], "data", "chains", chain_filename)
        if path.exists(builtin_chain_path):
            return builtin_chain_path
    
    # Fall back to download mechanism
    chain_path = check_and_download("{}to{}".format(from_build, to_build))
    return chain_path
####################################################################################################################
from gwaslab.io.io_gtf import gtf_to_protein_coding
from gwaslab.io.io_gtf import gtf_to_all_gene

####################################################################################################################   
# From BioPython: https://github.com/biopython/biopython/blob/c5a6b1374267d769b19c1022b4b45472316e78b4/Bio/Seq.py#L36
def _maketrans(complement_mapping: Dict[str, str]) -> bytes:
    """Make a python string translation table.

    Arguments:
     - complement_mapping - a dictionary.

    Returns a translation table (a bytes object of length 256) for use with
    the python string's translate method.

    Compatible with lower case and upper case sequences.
    """
    keys = "".join(complement_mapping.keys()).encode("ASCII")
    values = "".join(complement_mapping.values()).encode("ASCII")

    return bytes.maketrans(keys + keys.lower(), values + values.lower())
        
####################################################################################################################   
        
def _inch_to_point(inch: float) -> float:
    #dpi: Dots per Inch
    #points: 1/72 inch
    return inch*72
        
NA_STRINGS=["na","NA","Na","Nan","NaN","<NA>","null","NULL","#N/A","#VALUE!","N/A","n/a","missing",""]
