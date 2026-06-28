from typing import Optional, Dict, Any, Union, List, Tuple
import gzip
import os
from os import path
import json
import requests
import shutil
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_config import options
import re
from copy import deepcopy
from pathlib import Path
from gwaslab.info.g_version import _get_version

_GCST_KEYWORD_RE = re.compile(r"^GCST\d+$", re.IGNORECASE)


def infer_registry_kind(keyword: str, meta: Optional[Dict[str, Any]] = None) -> str:
    """Classify a registry entry as catalog reference or GWAS Catalog sumstats.
"""
    if isinstance(meta, dict):
        kind = meta.get("kind")
        if kind in ("ref", "sumstats"):
            return kind
        if meta.get("source") == "gwas_catalog":
            return "sumstats"
    if _GCST_KEYWORD_RE.match(keyword or ""):
        return "sumstats"
    return "ref"


def backfill_registry_kinds(cfg: Dict[str, Any]) -> bool:
    """Ensure every downloaded entry has a ``kind`` field. Returns True if modified.
"""
    changed = False
    for key, meta in cfg.get("downloaded", {}).items():
        if isinstance(meta, dict) and not meta.get("kind"):
            meta["kind"] = infer_registry_kind(key, meta)
            changed = True
    return changed


def filter_downloaded_registry(
    entries: Dict[str, Any],
    *,
    source: Optional[str] = None,
    kind: Optional[str] = None,
) -> Dict[str, Any]:
    """Filter registry entries by ``source`` and/or ``kind`` metadata.
"""
    filtered: Dict[str, Any] = {}
    for key, meta in entries.items():
        if not isinstance(meta, dict):
            continue
        entry_kind = meta.get("kind") or infer_registry_kind(key, meta)
        if kind is not None and entry_kind != kind:
            continue
        if source is not None and meta.get("source") != source:
            continue
        filtered[key] = meta
    return filtered

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
    '''Create an empty configuration file for managing downloaded reference data records.

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
    '''Update the configuration file or create a new one if missing.

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
    log.write(" -Updating config.json...", verbose=verbose)
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

    for key, value in dicts["downloaded"].items():
        if type(value) is not dict:
            to_remove.append(key)
            continue
        local_path = value.get("local_path")
        if not local_path or not path.exists(local_path):
            to_remove.append(key)
        else:
            if not value.get("kind"):
                value["kind"] = infer_registry_kind(key, value)
            key_list.append(key)
            filtered_dicts[key] = value

    log.write(" - Local file keywords: ", " ".join(key_list), verbose=verbose)

    for i in to_remove:
        dicts["downloaded"].pop(i, None)

    backfill_registry_kinds(dicts)

    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(dicts, f, indent=4)
    return filtered_dicts if not show_all else dicts["downloaded"]

def set_default_directory(dir_path: str, persist: bool = True) -> None:
    '''Set the default directory for reference data downloads.

Parameters
----------
dir_path : str
    Directory path to use for storing downloaded reference data.
persist : bool, optional
    When True (default), save the choice to ~/.gwaslab/settings.json.
'''
    dir_path = path.expanduser(dir_path)
    if not dir_path.endswith(os.sep) and not dir_path.endswith("/"):
        dir_path = dir_path + "/"
    options.set_option("data_directory", dir_path, persist=persist)
    os.makedirs(dir_path, exist_ok=True)

def get_default_directory() -> str:
    '''Get the configured default directory for reference data storage.

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
    '''Load and return the list of available reference files from configuration for downloading.

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
    '''Download and update the reference file dictionary from GitHub.
'''
    url = 'https://raw.github.com/Cloufield/gwaslab/main/src/gwaslab/data/reference.json'
    log.write("Updating available_ref list from:",url)
    r = requests.get(url, allow_redirects=True)
    data_path = options.paths["reference"]
    with open(data_path, 'wb') as file:
        file.write(r.content)
    log.write("Available_ref list has been updated!")

##################################################################################

def check_downloaded_ref(log: Log = Log(), verbose: bool = True) -> Dict[str, Any]:
    '''Verify and return records of local reference files that are already downloaded or manually added by user.

Returns
-------
dict
    Dictionary mapping reference keywords to local file paths
'''
    log.write("Start to check downloaded reference files...", verbose=verbose)
    config_path = options.paths["config"]
    log.write(" -Checking the config file:{}".format(config_path), verbose=verbose)
    if not path.exists(config_path):
        log.write(" -Config file is missing.", verbose=verbose)
        initiate_config()
        dicts = update_config(log, verbose=verbose)
        log.write("Finished checking downloaded reference files...", verbose=verbose)
        return dicts
    log.write(" -Config file exists.", verbose=verbose)
    dicts = update_config(log, verbose=verbose)
    log.write("Finished checking downloaded reference files...", verbose=verbose)
    return dicts

##################################################################################

def get_path(
    name: str,
    log: Log = Log(),
    verbose: bool = True,
    raise_on_missing: bool = False,
) -> Union[str, bool]:
    '''Retrieve the local file path for a specified reference file using keywords.

    Returns the file path if found, otherwise False (or raises if raise_on_missing).
'''
    config_path = options.paths["config"]
    if not path.exists(config_path):
        log.write("Config file not exists...", verbose=verbose)
        log.write("Created new config file...", verbose=verbose)
        initiate_config()

    try:
        with open(config_path, "r", encoding="utf-8") as handle:
            registry = json.load(handle)
        entry = registry["downloaded"][name]
        if not isinstance(entry, dict):
            raise KeyError(name)
        local_path = entry.get("local_path")
        if not local_path:
            raise KeyError(f"missing local_path for '{name}'")
        if path.exists(local_path):
            return local_path
        log.write(f"File not found on disk: {local_path}", verbose=verbose)
    except KeyError as exc:
        log.write(
            f"No registry record for '{name}'. Please download or add_local_data first.",
            verbose=verbose,
        )
        if raise_on_missing:
            raise KeyError(str(exc)) from exc
    except (OSError, json.JSONDecodeError) as exc:
        log.write(f"Failed to read config: {exc}", verbose=verbose)
        if raise_on_missing:
            raise
    return False

##################################################################################
def download_ref(
    name: str,
    directory: Optional[str] = None,
    local_filename: Optional[str] = None,
    overwrite: bool = False,
    log: Log = Log()
) -> Optional[str]:
    '''Download a reference file based on its identifier from reference.json.

Parameters
----------
name : str
    Reference file identifier
Returns
-------
str or None
    Local path where the file was saved, or None if download failed
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
        expected_md5 = dicts[name].get("md5sum") or dicts[name].get("md5")
        if expected_md5:
            file_status = check_file_integrity(local_path=local_path, md5sum=expected_md5, log=log)
            if file_status == 0:
                log.write("Md5sum verification of ", name, " failed! Please check again.")
                return None
        update_record("downloaded", name, "local_path", value=local_path)
        
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
        update_description(name, {**dicts[name], "source": "catalog", "kind": "ref"})
        log.write("Downloaded ", name, " successfully!")
        return local_path
    else:
        log.write(name, " is not available. Please use check_available_ref() to check available reference files.")
        return None

#### helper #############################################################################################

def check_file_integrity(local_path: str, md5sum: str, log: Log) -> int:
    '''Calculate and verify the MD5 checksum of a file.

    Returns 1 if MD5 matches, 0 if verification fails.
'''
    from gwaslab.bd.bd_io import verify_md5
    return 1 if verify_md5(local_path, md5sum, log=log) else 0

def remove_file(name: str, log: Log = Log()) -> None:
    '''Remove a reference file and its registry record (alias for remove_local_record).
'''
    remove_local_record(name, delete_file=True, log=log)

#### helper #############################################################################################
def update_record(*keys: str, value: Any, log: Log = Log(), verbose: bool = True) -> None:
    """Update configuration with a value at a nested key path.
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
        log.write(" - Updating record in config file...", verbose=verbose)

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
    """Add or update local data file to configuration without downloading.
    
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
        "local_path": os.path.abspath(local_path),
        "description": description,
        "suggested_use": suggested_use,
        "md5sum": md5sum,
        "format": format,
        "source": "local",
        "kind": "ref",
    }
    
    # Add TBI index if provided
    if tbi and os.path.exists(tbi):
        meta["tbi"] = {"local_path": tbi}
    if csi and os.path.exists(csi):
        meta["csi"] = {"local_path": csi}

    # Update configuration with metadata
    update_description(keyword, meta, log=log)
    log.write(f"Successfully registered local data '{keyword}' in configuration")
    return True

def remove_local_record(
    keyword: str,
    delete_file: bool = False,
    log: Log = Log(),
) -> bool:
    """Remove a registry entry. Optionally delete the main file and .tbi sidecar.
"""
    config_path = options.paths["config"]

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            cfg = json.load(f)
    except Exception:
        log.write("Error: Could not load config file.")
        return False

    if "downloaded" not in cfg or keyword not in cfg["downloaded"]:
        log.write(f"Error: No record found for '{keyword}' in configuration.")
        return False

    meta = cfg["downloaded"][keyword]
    if delete_file and isinstance(meta, dict):
        local_path = meta.get("local_path")
        if local_path and os.path.exists(local_path):
            os.remove(local_path)
            log.write(f"Removed file: {local_path}")
        tbi_meta = meta.get("tbi")
        if isinstance(tbi_meta, dict):
            tbi_path = tbi_meta.get("local_path")
            if tbi_path and os.path.exists(tbi_path):
                os.remove(tbi_path)
                log.write(f"Removed index: {tbi_path}")

    del cfg["downloaded"][keyword]

    try:
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=4, ensure_ascii=False)
        log.write(f"Successfully removed record '{keyword}' from configuration")
        return True
    except Exception as e:
        log.write(f"Error: Failed to update config file - {str(e)}")
        return False

def _iter_scan_files(data_dir: str, recursive: bool) -> List[Tuple[str, str]]:
    """Return (basename, absolute_path) pairs under data_dir.
"""
    pairs: List[Tuple[str, str]] = []
    if not recursive:
        for name in os.listdir(data_dir):
            fp = os.path.join(data_dir, name)
            if os.path.isfile(fp):
                pairs.append((name, fp))
        return pairs
    for root, _dirs, files in os.walk(data_dir):
        for name in files:
            fp = os.path.join(root, name)
            pairs.append((name, fp))
    return pairs


def scan_downloaded_files(
    log=Log(),
    verbose=True,
    directory: Optional[str] = None,
    recursive: bool = False,
):
    """Scan data directory for files not in config and match with available references.
"""
    log.write("Starting to scan data directory for unregistered files...", verbose=verbose)

    if directory is None:
        data_dir = os.path.abspath(os.path.expanduser(options.paths["data_directory"]))
    else:
        data_dir = os.path.abspath(os.path.expanduser(directory))

    log.write(f" -Scanning directory: {data_dir}", verbose=verbose)
    if recursive:
        log.write(" -Recursive scan enabled", verbose=verbose)

    if not os.path.isdir(data_dir):
        log.write(f" -Directory does not exist: {data_dir}", verbose=verbose)
        return False

    try:
        file_pairs = _iter_scan_files(data_dir, recursive=recursive)
    except Exception as e:
        log.write(f" -Error reading directory: {str(e)}", verbose=verbose)
        return False

    downloaded = check_downloaded_ref(log, verbose=verbose) or {}

    registered_basenames = set()
    for _key, meta in downloaded.items():
        if isinstance(meta, dict) and meta.get("local_path"):
            registered_basenames.add(os.path.basename(meta["local_path"]))
        if isinstance(meta, dict) and isinstance(meta.get("tbi"), dict):
            tbi_path = meta["tbi"].get("local_path")
            if tbi_path:
                registered_basenames.add(os.path.basename(tbi_path))

    available = check_available_ref(log, verbose=False, show_all=True)

    for filename, local_path in file_pairs:
        if filename in registered_basenames:
            continue

        matched_ref = None
        matched_info = None
        for ref_name, ref_info in available.items():
            url_filename = os.path.basename(ref_info["url"])
            if filename == url_filename:
                matched_ref = ref_name
                matched_info = ref_info
                log.write(f"  -Matched with reference: {ref_name}", verbose=verbose)
                break

        if not matched_ref or matched_info is None:
            continue

        try:
            update_description(
                matched_ref,
                {
                    **matched_info,
                    "local_path": local_path,
                    "source": "catalog",
                    "kind": "ref",
                },
                log=log,
            )
            log.write(
                f"  -Registered {matched_ref} with path: {local_path}",
                verbose=verbose,
            )
        except Exception as e:
            log.write(f"  -Error updating config for {filename}: {str(e)}", verbose=verbose)

    log.write("Completed scanning data directory", verbose=verbose)
    return True

def _deep_merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge dictionary b into dictionary a.

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
    """Update configuration by deep-merging metadata into config["downloaded"][key].

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

def download_file(url: str, file_path: Optional[str] = None, overwrite: bool = True) -> Optional[str]:
    '''Low-level file download utility with streaming support.

    Returns the path where the file was saved.
'''
    from gwaslab.bd.bd_io import stream_download
    if file_path is not None:
        stream_download(url, file_path, timeout=20, overwrite=overwrite)
        return file_path
    return None

def url_to_local_file_name(
    local_filename: Optional[str],
    url: str,
    from_dropbox: int
) -> Tuple[str, int]:
    '''Convert URL to valid local filename, handling Dropbox special cases.

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

def check_and_download(name: str, log: Log = Log()) -> str:
    '''Ensure a reference file exists by downloading if necessary.

    Returns the local path. Raises RuntimeError when download fails.
'''
    data_path = get_path(name, log=log, verbose=False)
    if data_path is not False and path.exists(data_path):
        return data_path

    dir_path = get_default_directory()
    if not path.exists(dir_path):
        os.makedirs(dir_path)
    downloaded_path = download_ref(name, directory=dir_path, log=log)
    if downloaded_path is None:
        raise RuntimeError(f"Failed to download reference '{name}'.")
    return downloaded_path

def search_local(file_path: str) -> bool:
    return path.exists(file_path)

##### format book ###################################################################################################

def update_formatbook(log: Log = Log()) -> None:
    '''Download and update the formatbook dictionary from GitHub.
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

def list_formats_with_descriptions(log: Log = Log(), *, silent: bool = False) -> List[Tuple[str, str]]:
    """Return sorted (format keyword, description) entries from the formatbook.

    Description is taken from ``meta_data.format_name`` when present.

Parameters
----------
silent : bool, default False
    If True, do not write the "Available formats" line to the log (for CLI tables).
"""
    data_path = options.paths["formatbook"]
    with open(data_path, encoding="utf-8") as f:
        book = json.load(f)
    rows: List[Tuple[str, str]] = []
    for key in sorted(book.keys()):
        meta = book[key].get("meta_data") or {}
        raw = meta.get("format_name", "")
        if raw is None:
            desc = ""
        elif isinstance(raw, str):
            desc = raw
        else:
            desc = str(raw)
        rows.append((key, desc))
    if not silent:
        log.write("Available formats:", ",".join(k for k, _ in rows))
    return rows


def list_formats(log: Log = Log()) -> List[str]:
    '''Display all available formats in the formatbook for GWASLab.
'''
    return [k for k, _ in list_formats_with_descriptions(log)]

def check_format(fmt: str, log: Log = Log()) -> Dict[str, str]:
    '''Check the header conversion dictionary between a given format and GWASLab format.

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
