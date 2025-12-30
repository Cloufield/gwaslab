from typing import Optional, Dict, Any, Union, List, Tuple
import gzip
import os
from os import path
import json
import requests
import shutil
import hashlib
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_config import options
import re
from copy import deepcopy
from pathlib import Path
from gwaslab.info.g_version import _get_version

"""
Update configuration by deep-merging metadata into config["downloaded"][key].

FUNCTION CORRELATION:
add_local_data() → update_description() → _deep_merge() → config.json
download_ref() → update_record() → config.json
remove_local_record() → direct config deletion

Key workflow:
1. add_local_data() prepares metadata → passes to update_description()
2. update_description() uses _deep_merge() to safely update config
3. download_ref() uses update_record() for atomic path registration
4. remove_local_record() directly manipulates config structure
5. All operations maintain consistent config.json structure

This creates a cohesive reference management system where:
- Configuration updates flow through standardized pathways
- Path handling is consistent across all operations
- Metadata integrity is maintained through deep merging
- Local and downloaded resources share the same configuration schema
"""
#### config ##############################################################################################
# config.json
# {
#  "downloaded":{}
#}
############################################################################################################

def initiate_config(log: Log = Log()) -> None:
    '''
    Create an empty configuration file for managing downloaded reference data records.

    This function initializes a new configuration file with an empty "downloaded" dictionary
    structure to track reference data records. The configuration file is created at the
    default configuration path specified in the options.
    '''
    config_path = options.paths["config"]
    log.write(" -Creating an empty config file... ")
    with open(config_path, 'w') as f:
        dict={'downloaded':{}}
        json.dump(dict,f,indent=4) 
        log.write(" -Config file path:",config_path)

def update_config(log: Log = Log(), verbose: bool = True, show_all: bool = False) -> Dict[str, Any]:
    '''
    Update the configuration file or create a new one if missing.

    This function ensures the configuration file exists and contains valid data,
    creating a new one if necessary. It also cleans up invalid entries by removing
    references to non-existent files.

    Returns
    -------
    dict
        Dictionary of currently recorded downloaded files with their metadata. The keys
        are reference identifiers and values contain metadata including local paths.
    '''
    config_path = options.paths["config"]
    log.write(" -Updating config.json...")
    try:
        # if config exists
        dicts = json.load(open(config_path))
    except:
        # if config not exists
        if not path.exists(config_path):
            initiate_config()
            dicts = json.load(open(config_path))
    
    # check if the ref file exists. If not, remove it from dicts.
    to_remove=[]

    key_list = []
    filtered_dicts = {}

    for key,value in dicts["downloaded"].items():
        if type(value) is not dict: 
            to_remove.append(key)
            continue
        if not path.exists(value["local_path"]):
            to_remove.append(key)
        else:
            key_list.append(key)
            filtered_value = {}
            for value_key, value_value in value.items():
                if value_key in ["local_path","description","suggested_use"]:
                    filtered_value[value_key] = value_value
            filtered_dicts[key] = filtered_value

    log.write(" - Local file keywords: "," ".join(key_list), verbose=verbose)

    # remove from dict
    for i in to_remove:
        dicts["downloaded"].pop(i, None)
    
    # write the dicts
    with open(config_path, 'w') as f:
        json.dump(dicts,f,indent=4)
    # also return the dics
    return filtered_dicts if not show_all else dicts["downloaded"]

def set_default_directory(path: str) -> None:
    '''
    Set a temporary default directory for reference data downloads.

    Parameters
    ----------
    path : str
        Directory path to use for storing downloaded reference data. This overrides
        the default data directory specified in the configuration.
    '''
    options.set_option("data_directory", path)

def get_default_directory() -> str:
    '''
    Get the configured default directory for reference data storage.

    Returns
    -------
    str
        Path to the default data directory from configuration
    '''
    return options.paths["data_directory"]

##################################################################################
def check_available_ref(log: Log = Log(), 
                        show_all: bool = False,
                        verbose: bool = True) -> Dict[str, Any]:
    '''
    Load and return the list of available reference files from configuration for downloading.

    Returns
    -------
    dict
        Dictionary mapping reference names to URLs
    '''
    log.write("Start to check available reference files...", verbose=verbose)
    ref_path = options.paths["reference"]
    if not path.exists(ref_path):
        # if reference.json not exists
        update_available_ref()
    dicts = json.load(open(ref_path))
    if dicts is not None:
        key_list = []
        filtered_dicts = {}
        for key, value in dicts.items():
            key_list.append(key)
            filtered_value = {}
            for value_key, value_value in value.items():
                if value_key in ["description","suggested_use"]:
                    filtered_value[value_key] = value_value
            
            #log.write(" -",key,":",filtered_value, verbose=verbose, show_time=False)
            filtered_dicts[key] = filtered_value

        log.write(" - Available keywords: "," ".join(key_list), verbose=verbose)
        return filtered_dicts if not show_all else dicts
    else:
        log.write(" -No available reference files.", verbose=verbose)
    log.write("Finished checking available reference files...", verbose=verbose)
    return {}

def update_available_ref(log: Log = Log()) -> None:
    '''
    Download and update the reference file dictionary from GitHub.
    '''
    url = 'https://raw.github.com/Cloufield/gwaslab/main/src/gwaslab/data/reference.json'
    log.write("Updating available_ref list from:",url)
    r = requests.get(url, allow_redirects=True)
    data_path = options.paths["reference"]
    with open(data_path, 'wb') as file:
        file.write(r.content)
    log.write("Available_ref list has been updated!")

##################################################################################

def check_downloaded_ref(log: Log = Log()) -> Dict[str, Any]:
    '''
    Verify and return records of local reference files that are already downloaded or manually added by user.

    Returns
    -------
    dict
        Dictionary mapping reference keywords to local file paths
    '''
    log.write("Start to check downloaded reference files...")
    config_path = options.paths["config"]
    log.write(" -Checking the config file:{}".format(config_path))
    if not path.exists(config_path):
        log.write(" -Config file is missing.")
        initiate_config()      
    else:
        log.write(" -Config file exists.")
        dicts = update_config()
        return dicts
    log.write("Finished checking downloaded reference files...")

##################################################################################

def get_path(name: str, log: Log = Log(), verbose: bool = True) -> Union[str, bool]:
    '''
    Retrieve the local file path for a specified reference file using keywords. 

    Keywords can be found using:
    - check_downloaded_ref for already downloaded files.
    - check_available_ref for available files for downloading.

    Parameters
    ----------
    name : str
        Reference file identifier
    verbose : bool, optional
        Whether to show detailed logging output

    Returns
    -------
    str or bool
        File path if found, False otherwise
    '''
    config_path = options.paths["config"]
    if not path.exists(config_path):
        log.write("Config file not exists...", verbose=verbose)
        log.write("Created new config file...", verbose=verbose)
        initiate_config()      
    else:
        try:
            dicts = json.load(open(config_path))["downloaded"]
            if path.exists(dicts[name]["local_path"]):
                return dicts[name]["local_path"]
            else:
                log.write("File not exist.", verbose=verbose)
        except:
            log.write("No records in config file. Please download first.", verbose=verbose)
    return False

##################################################################################
def download_ref(
    name: str,
    directory: Optional[str] = None,
    local_filename: Optional[str] = None,
    overwrite: bool = False,
    log: Log = Log()
) -> None:
    '''
    Download a reference file based on its identifier from reference.json.

    Parameters
    ----------
    name : str
        Reference file identifier

    Returns
    -------
    None
        File is saved to disk with path recorded in config
    '''
    from_dropbox=0
    dicts = check_available_ref(log,verbose=False,show_all=True)
    
    if name in dicts.keys():
        # get url for name
        url = dicts[name]["url"]
        
        log.write("Start to download ",name," ...")

        # get file local path
        if directory is None:
            directory = options.paths["data_directory"]
            if not path.exists(directory):
                os.makedirs(directory)
        
        local_filename, from_dropbox = url_to_local_file_name(local_filename, url, from_dropbox)

        local_path = os.path.join(directory, local_filename)
        log.write(" -Downloading to:",local_path)
        
        # if existing in default path
        if search_local(local_path) == True:
            log.write(" -File {} exists.".format(local_path))
            if overwrite == True:
                log.write(" -Overwriting the existing file.")
                download_file(url,local_path)
        else:
            download_file(url,local_path)
        
        # update record in config json
        if "md5" in  dicts[name].keys():
            file_status = check_file_integrity(local_path=local_path, md5sum=dicts[name]["md5"],log=log)
            if file_status==0:
                log.write("Md5sum verification of ",name," failed! Please check again.")
        update_record("downloaded", name,"local_path",value=local_path)
        
        # if vcf.gz -> check tbi
        if local_path.endswith("vcf.gz"):
                if "tbi" in  dicts[name].keys():
                    tbi_url = dicts[name]["tbi"]["url"]
                try:
                    log.write(" -Downloading to:",local_path+".tbi")
                    tbi_path = local_path + ".tbi"
                    if search_local(tbi_path) == True:
                        log.write(" -File {} exists.".format(tbi_path))
                        if overwrite == True:
                            log.write(" -Overwriting the existing file.")
                            download_file(tbi_url, tbi_path)
                    else:
                        download_file(tbi_url, tbi_path)

                    update_record("downloaded", name, "tbi", "local_path", value=tbi_path)
                except:
                    pass
        
        # if fasta.gz or fa.gz, decompress
        if local_path.endswith("fa.gz") or local_path.endswith("fasta.gz"):
            log.write(" -gunzip :",local_path)
            try:
                with gzip.open(local_path, 'rb') as f_in:
                    decompressed_path = local_path[:-3]
                    with open(decompressed_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        update_record("downloaded", name, "local_path", value=decompressed_path)
            except:
                pass
        update_description(name, dicts[name])
        log.write("Downloaded ",name," successfully!")
    else:
        log.write(name," is not available. Please use check_available_ref() to check available reference files.")

#### helper #############################################################################################

def check_file_integrity(local_path: str, md5sum: str, log: Log) -> int:
    '''
    Calculate and verify the MD5 checksum of a file.

    Parameters
    ----------
    local_path : str
        Path to the file to check
    md5sum : str
        Expected MD5 checksum
    log : gwaslab.g_Log.Log
        Logger object for recording progress and status messages

    Returns
    -------
    int
        1 if MD5 matches, 0 if verification fails
    '''
    md5_hash = hashlib.md5()
    with open(local_path,"rb") as f:
        # Read and update hash in chunks
        for byte_block in iter(lambda: f.read(4096*1000),b""):
            md5_hash.update(byte_block)
    log.write(" -File path: {}".format(local_path))
    log.write(" -MD5 check: {}".format(str(md5_hash.hexdigest())))
    if str(md5_hash.hexdigest()) == md5sum:
        log.write(" -MD5 verified.")
        return 1
    else:
        log.warning("-MD5 VERIFICATION FAILED!")
        return 0

def remove_file(name: str, log: Log = Log()) -> None:
    '''
    Remove a reference file and its record from the configuration.

    Parameters
    ----------
    name : str
        Reference file identifier
    '''
    log.write("Start to remove ",name," ...")
    config_path = options.paths["config"]
    if not path.exists(config_path):
        log.write("Config file not exists...")
    else:
        try:
            dicts = json.load(open(config_path))["downloaded"]
            if path.exists(dicts[name]):
                os.remove(dicts[name])
                log.write("Removed :" , dicts[name])
                check_downloaded_ref()
            else:
                log.write("File not exist.")
        except:
            log.write("No records in config file. Please download first.")

#### helper #############################################################################################
def update_record(*keys: str, value: Any, log: Log = Log()) -> None:
    """
    Update configuration with a value at a nested key path.
    Automatically creates intermediate levels if they do not exist.
    
    Example 1:
    {
    "downloaded": {
        "<resource_key>": {
        "local_path": "<absolute/path/to/local/file>",
        "url": "<download_URL or empty>",
        "description": "<short human description>",
        "suggested_use": "<what this file is typically used for>",
        "related_files": "<comma-separated keys or empty>",
        "tbi": {
            "local_path": "<absolute/path/to/local/file.tbi>",
            "url": "<download_url_for_tbi_or_empty>"
        },
        "md5sum": "<expected_md5_or_empty>"
        }
    }
    }
    
    Example 2:
    {
    "downloaded": {
        "<resource_key>": {
        "local_path": "<absolute/path/to/local/file>"}
        }
    }

    Parameters
    ----------
    *keys : str
        One or more hierarchical keys. For example:
        update_record("1kg_eas_hg19", "tbi", value="/path/to/file")
    value : any
        Value to be stored at that nested path.
    """
    if log is not None:
        log.write(" - Updating record in config file...")

    config_path = Path(options.paths["config"])

    # Load config or initialize
    if config_path.exists():
        with open(config_path, "r", encoding="utf-8") as f:
            try:
                cfg = json.load(f)
            except Exception:
                cfg = {}
    else:
        cfg = {}

    # Navigate and build nested dict
    node = cfg
    for k in keys[:-1]:
        if k not in node or not isinstance(node[k], dict):
            node[k] = {}
        node = node[k]

    # Assign the final value
    node[keys[-1]] = value

    # Write back
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=4, ensure_ascii=False)

def add_local_data(
    keyword: str,
    local_path: str,
    format: Optional[str] = None,
    description: Optional[str] = None,
    md5sum: Optional[str] = None,
    suggested_use: Optional[str] = None,
    tbi: Optional[str] = None,
    csi: Optional[str] = None,
    log: Log = Log()
) -> bool:
    """
    Add or update local data file to configuration without downloading.
    
    Parameters
    ----------
    keyword : str
        Reference identifier for the local data
    local_path : str
        Absolute path to the local data file
    format : str
        Data format identifier (e.g., 'vcf', 'bcf', 'plink', 'gtf', 'fasta','tsv')
    description : str
        Human-readable description of the data
    md5sum : str
        MD5 checksum for file integrity verification
    suggested_use : str
        Description of typical use cases for this data
    tbi : str, optional
        Path to corresponding .tbi index file (if applicable)
    csi : str, optional
        Path to corresponding .csi index file (if applicable)
    log : gwaslab.g_Log.Log, optional
        Logger instance for recording operations
    
    Returns
    -------
    bool
        True if successfully added, False otherwise
    """
    # Verify file exists
    if not os.path.exists(local_path):
        log.write(f"Error: Local file {local_path} does not exist.")
        return False
    
    # Prepare metadata structure
    meta = {
        "local_path": local_path,
        "description": description,
        "suggested_use": suggested_use,
        "md5sum": md5sum,
        "format": format
    }
    
    # Add TBI index if provided
    if tbi and os.path.exists(tbi):
        meta["tbi"] = {"local_path": tbi}
    if csi and os.path.exists(tbi):
        meta["csi"] = {"local_path": csi}

    # Update configuration with metadata
    update_description(keyword, meta, log=log)
    log.write(f"Successfully registered local data '{keyword}' in configuration")
    return True

def remove_local_record(keyword: str, log: Log = Log()) -> bool:
    """
    Remove a local data record from configuration without deleting the physical file.
    
    Parameters
    ----------
    keyword : str
        Reference identifier for the data to remove from configuration
    log : gwaslab.g_Log.Log, optional
        Logger instance for recording operations
    
    Returns
    -------
    bool
        True if successfully removed from configuration, False otherwise
    """
    config_path = options.paths["config"]
    
    # Load config
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            cfg = json.load(f)
    except Exception:
        log.write("Error: Could not load config file.")
        return False
    
    # Check if the keyword exists in downloaded records
    if "downloaded" not in cfg or keyword not in cfg["downloaded"]:
        log.write(f"Error: No record found for '{keyword}' in configuration.")
        return False
    
    # Remove the record
    del cfg["downloaded"][keyword]
    
    # Write back the updated config
    try:
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=4, ensure_ascii=False)
        log.write(f"Successfully removed record '{keyword}' from configuration")
        return True
    except Exception as e:
        log.write(f"Error: Failed to update config file - {str(e)}")
        return False

def scan_downloaded_files(log=Log(), verbose=True):
    """
    Scan data directory for files not in config and match with available references.

    Parameters
    ----------
    verbose : bool, optional
        Whether to show detailed logging output. Default is True.

    Returns
    -------
    bool
        True if successful
    """
    log.write("Starting to scan data directory for unregistered files...", verbose=verbose)
    
    # Get directory paths
    data_dir = options.paths["data_directory"]
    log.write(f" -Scanning directory: {data_dir}", verbose=verbose)
    
    # Get list of files in data directory
    try:
        files_in_dir = [f for f in os.listdir(data_dir) 
                      if os.path.isfile(os.path.join(data_dir, f))]
    except Exception as e:
        log.write(f" -Error reading directory: {str(e)}")
        return False
    
    # Get downloaded records from config
    downloaded = check_downloaded_ref(log)
    
    # Get available references
    available = check_available_ref(log, verbose=False,show_all=True)
    
    # Process each file in directory
    for filename in files_in_dir:
        # Skip if already in config
        if filename in downloaded:
            continue
            
        log.write(f" -Found unregistered file: {filename}", verbose=verbose)
        
        # Try to match with available references
        matched_ref = None
        for ref_name, ref_info in available.items():
            # Check if filename matches URL's last component
            url_filename = os.path.basename(ref_info['url'])
            if filename == url_filename:
                matched_ref = ref_name
                log.write(f"  -Matched with reference: {ref_name}", verbose=verbose)
                break
        
        if not matched_ref:
            log.write(f"  -No matching reference found for {filename}", verbose=verbose)
            continue
        
        # Update config with matched reference
        local_path = os.path.join(data_dir, filename)
        try:
            update_record("downloaded", matched_ref, "local_path", value=local_path)
            log.write(f"  -Updated config for {matched_ref} with path: {local_path}", verbose=verbose)
        except Exception as e:
            log.write(f"  -Error updating config for {filename}: {str(e)}")
    
    log.write("Completed scanning data directory", verbose=verbose)
    return True

def _deep_merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merge dictionary b into dictionary a.

    Parameters
    ----------
    a : dict
        The base dictionary to merge into
    b : dict
        The dictionary to merge into a

    Returns
    -------
    dict
        The merged dictionary
    """
    # If a is not a dict, b replaces it entirely
    if not isinstance(a, dict):
        return b

    # Both a and b must be dicts to merge further
    for k, v in b.items():
        if isinstance(v, dict):
            a[k] = _deep_merge(a.get(k), v)
        else:
            a[k] = v

    return a

def update_description(key: str, dicts: Dict[str, Any], log: Optional[Log] = None) -> None:
    """
    Update configuration by deep-merging metadata into config["downloaded"][key].

    Parameters
    ----------
    key : str
        Reference identifier.
    dicts : dict
        Possibly nested dictionary to merge into the entry.
    """
    if log is not None:
        log.write(" - Updating description records in config file...")

    config_path = options.paths["config"]

    # Load or initialize config
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            cfg = json.load(f)
    except Exception:
        cfg = {}

    # Ensure nested structure
    cfg.setdefault("downloaded", {})
    cfg["downloaded"].setdefault(key, {})

    # Deep merge
    cfg["downloaded"][key] = _deep_merge(cfg["downloaded"][key], dicts)

    # Write back
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=4, ensure_ascii=False)

    if log is not None:
        log.write(f"   -> Deep-merged into downloaded/{key}")

def download_file(url: str, file_path: Optional[str] = None) -> Optional[str]:
    '''
    Low-level file download utility with streaming support.

    Parameters
    ----------
    url : str
        Source URL to download from
    file_path : str, optional
        Destination path to save file

    Returns
    -------
    str
        Path where file was saved
    '''
    # download file from url to file_path
    if file_path is not None:
        with requests.get(url, stream=True,timeout=(20, 20)) as r:
            # checking status
            r.raise_for_status()
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)     
        return file_path

def url_to_local_file_name(
    local_filename: Optional[str],
    url: str,
    from_dropbox: int
) -> Tuple[str, int]:
    '''
    Convert URL to valid local filename, handling Dropbox special cases.

    Parameters
    ----------
    local_filename : str, optional
        Optional custom filename
    url : str
        Source URL
    from_dropbox : int
        Dropbox indicator flag

    Returns
    -------
    tuple
        (processed filename, dropbox flag)
    '''
    if local_filename is None:
        # if local name not provided, grab it from url
        local_filename = url.split('/')[-1]
        
    if local_filename.endswith("dl=1"):
        # if file are downloaded form dropbox
        # set from_dropbox indicator to 1
        from_dropbox=1
        # remove "?dl=1" suffix
        local_filename = re.match(r'([^\?]+)(\?rlkey=[\w]+)?[&\?]dl=1$', local_filename)
        local_filename = local_filename.group(1)
    return local_filename, from_dropbox

##########################################################################################################

def check_and_download(name: str) -> str:
    '''
    Ensure a reference file exists by downloading if necessary.

    Parameters
    ----------
    name : str
        Reference file identifier

    Returns
    -------
    str
        Path to the reference file (existing or newly downloaded)
    '''
    # if file exsits , return file path
    if get_path(name) is not False:
        data_path = get_path(name)
        if path.exists(data_path):
            return data_path    
    # if file not exsits, download and return new path
    dir_path = get_default_directory()
    if not path.exists(dir_path):
        os.makedirs(dir_path)
    download_ref(name,directory = dir_path)
    data_path = get_path(name)
    return data_path

def search_local(file_path: str) -> bool:
    return path.exists(file_path)

##### format book ###################################################################################################

def update_formatbook(log: Log = Log()) -> None:
    '''
    Download and update the formatbook dictionary from GitHub.
    '''
    url = 'https://raw.github.com/Cloufield/formatbook/main/formatbook.json'
    log.write("Updating formatbook from:",url)
    r = requests.get(url, allow_redirects=True)
    data_path = options.paths["formatbook"]
    log.write("Overwrite formatbook to : ",data_path)
    with open(data_path, 'wb') as file:
        file.write(r.content)
    book=json.load(open(data_path))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))
    log.write("Formatbook has been updated!")

def list_formats(log: Log = Log()) -> List[str]:
    '''
    Display all available formats in the formatbook for GWASLab.
    '''
    data_path = options.paths["formatbook"]
    book=json.load(open(data_path))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))    
    return available_formats

def check_format(fmt: str, log: Log = Log()) -> Dict[str, str]:
    '''
    Check the header conversion dictionary between a given format and GWASLab format.

    Parameters
    ----------
    fmt : str
        Format name to check
    log : Log, optional
        Logger instance
    
    Returns
    -------
    Dict[str, str]
        Dictionary mapping format headers to GWASLab headers
    '''
    data_path = options.paths["formatbook"]
    book=json.load(open(data_path))
    log.write("Available formats:",end="")
    for i in book[fmt].keys():
        log.write(i,end="")
    log.write("") 
    for i in book[fmt].values():
        log.write(i,end="")
    return book[fmt]
