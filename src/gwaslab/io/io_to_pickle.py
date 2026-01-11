from typing import TYPE_CHECKING, Optional, List
import pickle
import os
import gc
from gwaslab.info.g_Log import Log 
import sys
from gwaslab import g_Sumstats
from gwaslab.info import g_Log
import pandas as pd

# Try to import ctypes for malloc_trim (Linux only)
try:
    import ctypes
    import platform
    MALLOC_TRIM_AVAILABLE = (platform.system() == 'Linux')
except ImportError:
    MALLOC_TRIM_AVAILABLE = False
    ctypes = None

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def dump_pickle(glsumstats: 'Sumstats', path: str = "~/mysumstats.pickle", overwrite: bool = False) -> None:
    glsumstats.log.write("Start to dump the Sumstats Object.")
    if overwrite==False and os.path.exists(path):
        glsumstats.log.write(" -File exists. Skipping. If you want to overwrite, please use overwrite=True.")
    else:
        with open(path, 'wb') as file:
            glsumstats.log.write(" -Dump the Sumstats Object to : ", path)
            pickle.dump(glsumstats, file)
    glsumstats.log.write("Finished dumping.")

def load_pickle(path: str) -> 'Sumstats':
    if os.path.exists(path):
        try:
            with open(path, 'rb') as file:
                glsumstats =  pickle.load(file)
                glsumstats.log.write("Loaded dumped Sumstats object created using gwaslab>=v3.4.32")    
                glsumstats.log.write("Loaded dumped Sumstats object from : ", path)
                return glsumstats
        except:
            sys.modules['gwaslab.Sumstats'] = g_Sumstats
            sys.modules['gwaslab.Log'] = g_Log
            
            with open(path, 'rb') as file:
                glsumstats =  pickle.load(file)
                glsumstats.log.write("Loaded dumped Sumstats object created using gwaslab<v3.4.32")    
                glsumstats.log.write("Loaded dumped Sumstats object from : ", path)
                return glsumstats
    else:
        Log().write("File not exists : ", path)

def load_data_from_pickle(path: str, usecols: Optional[List[str]] = None) -> pd.DataFrame:
    data = load_pickle(path).data
    existing_cols = []
    if usecols is not None:
        for i in usecols:
            if i in data.columns:
                existing_cols.append(i)
        data = data.loc[:,existing_cols]
        gc.collect()
    return data

def _force_memory_release():
    """
    Force Python's memory allocator to release memory back to OS.
    Works on Linux by calling malloc_trim() from glibc.
    
    Note: This only works on Linux with glibc. On other systems,
    Python's memory allocator may not release memory back to the OS
    even after objects are deleted and garbage collected.
    """
    if MALLOC_TRIM_AVAILABLE:
        try:
            # Call malloc_trim(0) to force glibc to return freed memory to OS
            # This is the key to actually seeing memory reduction in RSS
            libc = ctypes.CDLL("libc.so.6")
            result = libc.malloc_trim(0)
            # malloc_trim returns 1 if memory was released, 0 otherwise
            return result == 1
        except (OSError, AttributeError):
            # malloc_trim not available or failed
            return False
    return False


def _offload(df: pd.DataFrame, path: str, log: Log) -> None:
    # Get DataFrame size for logging
    df_size_mb = df.memory_usage(deep=True).sum() / (1024 * 1024) if hasattr(df, 'memory_usage') else 0
    
    with open(path, 'wb') as file:
        pickle.dump(df, file)
        log.write("Dumpping dataframe to : ", path)
        if df_size_mb > 0:
            log.write("  -DataFrame size: {:.2f} MB".format(df_size_mb))
    
    # Force memory release after pickling
    del df
    gc.collect()
    _force_memory_release()  # Force memory allocator to release memory to OS (Linux only)

def _reload(path: str, log: Log, delete_files: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Reload data from temporary pickle file.
    
    Parameters
    ----------
    path : str
        Path to the pickle file to reload
    log : Log
        Logger instance
    delete_files : list of str, optional
        Additional files to delete after successful reload
        
    Returns
    -------
    pd.DataFrame
        Reloaded dataframe
    """
    with open(path, 'rb') as file:
        df =  pickle.load(file)
        log.write("Loaded dataframe back from : ", path)
    try:
        os.remove(path)
    except:
        pass
    
    # Delete additional files if provided
    if delete_files is not None:
        n_deleted = 0
        for file_path in delete_files:
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                    n_deleted += 1
                except Exception as e:
                    log.write(" -Warning: Could not delete file {}: {}".format(file_path, str(e)))
        if n_deleted > 0:
            log.write(" -Cleaned up {} additional file(s) after successful reload...".format(n_deleted))
    
    return df