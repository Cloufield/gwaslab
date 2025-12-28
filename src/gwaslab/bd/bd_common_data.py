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

#hard-coded data
def get_chr_to_NC(build, inverse=False):
    """
    Create a dictionary mapping chromosome numbers to NCBI accession IDs.
    
    Parameters:
    build (str): Genome build version ('19' or '38')
    inverse (bool): If True, return NCBI ID to chromosome mapping
    
    Returns:
    dict: Dictionary mapping chromosome identifiers (string) to NCBI accession IDs (string)
    """
    #https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
    if build =="19":
        dic={
        "1":"NC_000001.10",
        "2":"NC_000002.11",
        "3":"NC_000003.11",
        "4":"NC_000004.11",
        "5":"NC_000005.9",
        "6":"NC_000006.11",
        "7":"NC_000007.13",
        "8":"NC_000008.10",
        "9":"NC_000009.11",
        "10":"NC_000010.10",
        "11":"NC_000011.9",
        "12":"NC_000012.11",
        "13":"NC_000013.10",
        "14":"NC_000014.8",
        "15":"NC_000015.9",
        "16":"NC_000016.9",
        "17":"NC_000017.10",
        "18":"NC_000018.9",
        "19":"NC_000019.9",
        "20":"NC_000020.10",
        "21":"NC_000021.8",
        "22":"NC_000022.10",
        "X":"NC_000023.10",
        "Y":"NC_000024.9",
        "MT":"NC_012920.1"}
    elif build=="38":
        dic={
        "1":"NC_000001.11",
        "2":"NC_000002.12",
        "3":"NC_000003.12",
        "4":"NC_000004.12",
        "5":"NC_000005.10",
        "6":"NC_000006.12",
        "7":"NC_000007.14",
        "8":"NC_000008.11",
        "9":"NC_000009.12",
        "10":"NC_000010.11",
        "11":"NC_000011.10",
        "12":"NC_000012.12",
        "13":"NC_000013.11",
        "14":"NC_000014.9",
        "15":"NC_000015.10",
        "16":"NC_000016.10",
        "17":"NC_000017.11",
        "18":"NC_000018.10",
        "19":"NC_000019.10",
        "20":"NC_000020.11",
        "21":"NC_000021.9",
        "22":"NC_000022.11",
        "X":"NC_000023.11",
        "Y":"NC_000024.1",
        "MT":"NC_012920.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic

def get_NC_to_chr(build):
    """
    Create a dictionary mapping NCBI accession IDs to chromosome numbers.
    
    Parameters:
    build (str): Genome build version ('19' or '38')
    
    Returns:
    dict: Dictionary mapping NCBI accession IDs (string) to chromosome numbers (string)
    """
    return get_chr_to_NC(build=build,inverse=True)


def get_number_to_NC(build, inverse=False):
    """
    Create a dictionary mapping chromosome numbers (int) to NCBI accession IDs (string).
    
    Parameters:
    build (str): Genome build version ('19' or '38')
    inverse (bool): If True, return NCBI ID to chromosome mapping
    
    Returns:
    dict: Dictionary mapping chromosome identifiers (int) to NCBI accession IDs
    """
    #https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
    if build =="19":
        dic={
        1:"NC_000001.10",
        2:"NC_000002.11",
        3:"NC_000003.11",
        4:"NC_000004.11",
        5:"NC_000005.9",
        6:"NC_000006.11",
        7:"NC_000007.13",
        8:"NC_000008.10",
        9:"NC_000009.11",
        10:"NC_000010.10",
        11:"NC_000011.9",
        12:"NC_000012.11",
        13:"NC_000013.10",
        14:"NC_000014.8",
        15:"NC_000015.9",
        16:"NC_000016.9",
        17:"NC_000017.10",
        18:"NC_000018.9",
        19:"NC_000019.9",
        20:"NC_000020.10",
        21:"NC_000021.8",
        22:"NC_000022.10",
        23:"NC_000023.10",
        24:"NC_000024.9",
        25:"NC_012920.1"}
    elif build=="38":
        dic={
        1:"NC_000001.11",
        2:"NC_000002.12",
        3:"NC_000003.12",
        4:"NC_000004.12",
        5:"NC_000005.10",
        6:"NC_000006.12",
        7:"NC_000007.14",
        8:"NC_000008.11",
        9:"NC_000009.12",
        10:"NC_000010.11",
        11:"NC_000011.10",
        12:"NC_000012.12",
        13:"NC_000013.11",
        14:"NC_000014.9",
        15:"NC_000015.10",
        16:"NC_000016.10",
        17:"NC_000017.11",
        18:"NC_000018.10",
        19:"NC_000019.10",
        20:"NC_000020.11",
        21:"NC_000021.9",
        22:"NC_000022.11",
        23:"NC_000023.11",
        24:"NC_000024.1",
        25:"NC_012920.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic


def get_NC_to_number(build):
    """
    Create a dictionary mapping NCBI accession IDs to chromosome numbers (int).
    
    Parameters:
    build (str): Genome build version ('19' or '38')
    
    Returns:
    dict: Dictionary mapping NCBI accession IDs to chromosome numbers (int)
    """
    return get_number_to_NC(build=build,inverse=True)

def get_chr_list(add_number=False, n=25, only_number=False, species="homo sapiens"):
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
def get_chr_to_number(out_chr=False, xymt=["X","Y","MT"], xymt_num=[23,24,25], species="homo sapiens", max_chr=200):
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
def get_number_to_chr(in_chr=False, xymt=["X","Y","MT"], xymt_num=[23,24,25], prefix="", species="homo sapiens", max_chr=200):
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
def get_high_ld(build="19"):
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
def get_format_dict(fmt, inverse=False):
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
def get_formats_list():
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

def get_recombination_rate(chrom, build="19"):
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
def get_chain(from_build="19", to_build="38"):    
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
def _maketrans(complement_mapping):
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
        
def _inch_to_point(inch):
    #dpi: Dots per Inch
    #points: 1/72 inch
    return inch*72
        
NA_STRINGS=["na","NA","Na","Nan","NaN","<NA>","null","NULL","#N/A","#VALUE!","N/A","n/a","missing",""]
