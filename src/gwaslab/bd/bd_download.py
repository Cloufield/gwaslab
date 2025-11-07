import gzip
import os
from os import path
import json
import requests
import shutil
import hashlib
from gwaslab.g_Log import Log
from gwaslab.bd.bd_config import options
import re
from copy import deepcopy
from pathlib import Path
from gwaslab.g_version import _get_version
#### config ##############################################################################################
# config.json
# {
#  "downloaded":{}
#}
############################################################################################################

def initiate_config(log=Log()):
    '''
    Create an empty configuration file for managing downloaded reference data records.
    
    Args:
        log (Log): Logging object for tracking operations. Defaults to new Log instance.
    '''
    #config_path=path.dirname(__file__) + '/data/config.json'
    config_path = options.paths["config"]
    log.write(" -Creating an empty config file... ")
    with open(config_path, 'w') as f:
        dict={'downloaded':{}}
        json.dump(dict,f,indent=4) 
        log.write(" -Config file path:",config_path)

def update_config(log=Log()):
    '''
    Update the configuration file or create a new one if missing.
    
    Args:
        log (Log): Logging object for tracking operations. Defaults to new Log instance.
    
    Returns:
        dict: Dictionary of currently recorded downloaded files
    '''
    #config_path=path.dirname(__file__) + '/data/config.json',
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
    for key,value in dicts["downloaded"].items():
        if type(value) is not dict: 
            to_remove.append(key)
            continue
        if not path.exists(value["local_path"]):
            to_remove.append(key)
        else:
            log.write("  -",key,":", dicts["downloaded"][key])

    # remove from dict
    for i in to_remove:
        dicts["downloaded"].pop(i, None)
    
    # write the dicts
    with open(config_path, 'w') as f:
        json.dump(dicts,f,indent=4)
    # also return the dics
    return dicts["downloaded"]

def set_default_directory(path):
    '''
    Set a temporary default directory for reference data downloads.
    
    Args:
        path (str): Directory path to use for storing downloaded reference data
    '''

    options.set_option("data_directory", path)


def get_default_directory():
    '''
    Get the configured default directory for reference data storage.
    
    Returns:
        str: Path to the default data directory from configuration
    '''
    return options.paths["data_directory"]



##################################################################################
def check_available_ref(log=Log(),verbose=True):
    '''
    Load and return the list of available reference files from configuration for downloading.
    
    Args:
        log (Log): Logging object for tracking operations
        verbose (bool): Whether to show detailed logging output
    
    Returns:
        dict: Dictionary mapping reference names to URLs
    '''
    log.write("Start to check available reference files...", verbose=verbose)
    #ref_path = path.dirname(__file__) + '/data/reference.json'
    ref_path = options.paths["reference"]
    if not path.exists(ref_path):
        # if reference.json not exists
        update_available_ref()
    dicts = json.load(open(ref_path))
    if dicts is not None:
        for key,value in dicts.items():
            log.write(" -",key," : ",value, verbose=verbose)
        return dicts
    else:
        log.write(" -No available reference files.", verbose=verbose)
    log.write("Finished checking available reference files...", verbose=verbose)
    return {}

def update_available_ref(log=Log()):
    '''
    Download and update the reference file dictionary from GitHub.
    
    Args:
        log (Log): Logging object for tracking operations
    '''
    url = 'https://raw.github.com/Cloufield/gwaslab/main/src/gwaslab/data/reference.json'
    log.write("Updating available_ref list from:",url)
    r = requests.get(url, allow_redirects=True)
    #data_path =  path.dirname(__file__) + '/data/reference.json'
    data_path = options.paths["reference"]
    with open(data_path, 'wb') as file:
        file.write(r.content)
    log.write("Available_ref list has been updated!")

##################################################################################

def check_downloaded_ref(log=Log()):
    '''
    Verify and return records of already downloaded reference files.
    
    Args:
        log (Log): Logging object for tracking operations
    
    Returns:
        dict: Dictionary mapping reference keywords to local file paths
    '''
    log.write("Start to check downloaded reference files...")
    #config_path =  path.dirname(__file__) + '/data/config.json'
    config_path = options.paths["config"]
    log.write(" -Checking the config file:{}".format(config_path))
    if not path.exists(config_path):
        log.write(" -Config file is missing.")
        initiate_config()      
    else:
        log.write(" -Config file exists.")
        #try:
        dicts = update_config()
        return dicts
        #except:
        #    log.write(" -No records in config file.")
    log.write("Finished checking downloaded reference files...")

##################################################################################

def get_path(name,log=Log(),verbose=True):
    '''
    Retrieve the local file path for a specified reference file using keywords. 
    keywords can be found using:
      check_downloaded_ref for already downloaded files.
      check_available_ref for avaiable files for downloading.
    
    Args:
        name (str): Reference file identifier
        log (Log): Logging object for tracking operations
        verbose (bool): Whether to show detailed logging output
    
    Returns:
        str|bool: File path if found, False otherwise
    '''
    #config_path =  path.dirname(__file__) + '/data/config.json'
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
def download_ref(name,
                directory=None,
                local_filename=None,
                overwrite=False,
                log=Log()):
    '''
    Download a reference file based on its identifier from reference.json.
    
    Args:
        name (str): Reference file identifier
        log (Log): Logging object for tracking operations
    
    Returns:
        None: File is saved to disk with path recorded in config
    '''
    from_dropbox=0
    dicts = check_available_ref(log,verbose=False)
    
    if name in dicts.keys():
        # get url for name
        
        url = dicts[name]["url"]
        
        log.write("Start to download ",name," ...")

        # get file local path
        if directory is None:
            #directory = path.dirname(__file__) + '/data/'
            directory = options.paths["data_directory"]
            if not path.exists(directory):
                os.makedirs(directory)
        
        local_filename, from_dropbox = url_to_local_file_name(local_filename, url, from_dropbox)

        local_path = directory + local_filename
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
                    if search_local(local_path+".tbi") == True:
                        log.write(" -File {} exists.".format(local_path+".tbi"))
                        if overwrite == True:
                            log.write(" -Overwriting the existing file.")
                            download_file(tbi_url,local_path+".tbi")
                    else:
                        download_file(tbi_url,local_path+".tbi")

                    update_record("downloaded", name, "tbi", "local_path", value=local_path+ ".tbi")
                except:
                    pass
        
        # if fasta.gz or fa.gz, decompress
        if local_path.endswith("fa.gz") or local_path.endswith("fasta.gz"):
            log.write(" -gunzip :",local_path)
            try:
                with gzip.open(local_path, 'rb') as f_in:
                    with open(local_path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        update_record("downloaded",name, "local_path", value=local_path[:-3])
            except:
                    pass
        update_description(name, dicts[name])
        log.write("Downloaded ",name," successfully!")
    else:
        log.write(name," is not available. Please use check_available_ref() to check available reference files.")

#### helper #############################################################################################

def check_file_integrity(local_path, md5sum,log):
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

def remove_file(name,log=Log()):
    log.write("Start to remove ",name," ...")
    #config_path =  path.dirname(__file__) + '/data/config.json'
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
def update_record(*keys, value, log=Log()):
    """
    Update configuration with a value at a nested key path.
    Automatically creates intermediate levels if they do not exist.

    Parameters
    ----------
    *keys : str
        One or more hierarchical keys. For example:
        update_record("1kg_eas_hg19", "tbi", value="/path/to/file")
    value : Any
        Value to be stored at that nested path.
    log : Log, optional
        Logging object for tracking operations.
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

def scan_downloaded_files(log=Log(), verbose=True):
    """
    Scan data directory for files not in config and match with available references.
    
    Args:
        log (Log): Logging object for tracking operations
        verbose (bool): Whether to show detailed logging output
    
    Returns:
        bool: True if successful
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
    available = check_available_ref(log, verbose=False)
    
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

    #if log is not None:
    #    log.write(f"   -> Updated path: {'/'.join(keys)}")

def _deep_merge(a, b):
    """
    Recursively merge dictionary b into dictionary a.
    If 'a' is not a dict, return 'b' (full replacement).
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

def update_description(key, dicts, log=None):
    """
    Update configuration by deep-merging metadata into config["downloaded"][key].

    Parameters
    ----------
    key : str
        Reference identifier.
    dicts : dict
        Possibly nested dictionary to merge into the entry.
    log : Log, optional
        Logging object for tracking operations.
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


def download_file(url, file_path=None):
    '''
    Low-level file download utility with streaming support.
    
    Args:
        url (str): Source URL to download from
        file_path (str): Destination path to save file
    
    Returns:
        str: Path where file was saved
    '''
    # download file from url to file_path
    if file_path is not None:
        with requests.get(url, stream=True,timeout=(20, 20)) as r:
            # checking status
            r.raise_for_status()
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)     
        return file_path

def url_to_local_file_name(local_filename, url, from_dropbox):
    '''
    Convert URL to valid local filename, handling Dropbox special cases.
    
    Args:
        local_filename (str): Optional custom filename
        url (str): Source URL
        from_dropbox (int): Dropbox indicator flag
    
    Returns:
        tuple: (processed filename, dropbox flag) 
    '''
    if local_filename is None:
        # if local name not provided, grab it from url
        local_filename = url.split('/')[-1]
        
    if local_filename.endswith("dl=1"):
        # if file are downloaded form dropbox
        # set from_dropbox indicator to 1
        from_dropbox=1
        # remove "?dl=1" suffix
        #local_filename = local_filename[:-5]
        local_filename = re.match(r'([^\?]+)(\?rlkey=[\w]+)?[&\?]dl=1$', local_filename)
        local_filename = local_filename.group(1)
    return local_filename, from_dropbox

##########################################################################################################

def check_and_download(name):
    '''
    Ensure a reference file exists by downloading if necessary.
    
    Args:
        name (str): Reference file identifier
    
    Returns:
        str: Path to the reference file (existing or newly downloaded)
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

def search_local(file_path):
    return path.exists(file_path)
##### format book ###################################################################################################

def update_formatbook(log=Log()):
    '''
    Download and update the formatbook dictionary from GitHub.
    
    Args:
        log (Log): Logging object for tracking operations
    '''
    url = 'https://raw.github.com/Cloufield/formatbook/main/formatbook.json'
    log.write("Updating formatbook from:",url)
    r = requests.get(url, allow_redirects=True)
    #data_path =  path.dirname(__file__) + '/data/formatbook.json'
    data_path = options.paths["formatbook"]
    log.write("Overwrite formatbook to : ",data_path)
    with open(data_path, 'wb') as file:
        file.write(r.content)
    book=json.load(open(data_path))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))
    log.write("Formatbook has been updated!")
         
def list_formats(log=Log()):
    '''
    Display all available formats in the formatbook.
    
    Args:
        log (Log): Logging object for tracking operations
    '''
    #data_path =  path.dirname(__file__) + '/data/formatbook.json'
    data_path = options.paths["formatbook"]
    book=json.load(open(data_path))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))    

def check_format(fmt,log=Log()):
    '''
    Check the header conversion dictionary between a given format and GWASLab format.
    
    Args:
        fmt (str): Format name to check
        log (Log): Logging object for tracking operations
    '''
    #data_path =  path.dirname(__file__) + '/data/formatbook.json'
    data_path = options.paths["formatbook"]
    book=json.load(open(data_path))
    log.write("Available formats:",end="")
    for i in book[fmt].keys():
        log.write(i,end="")
    log.write("") 
    for i in book[fmt].values():
        log.write(i,end="")
