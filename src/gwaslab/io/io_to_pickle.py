from typing import TYPE_CHECKING, Optional, List, Union
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
    from gwaslab.g_SumstatsPair import SumstatsPair
    from gwaslab.g_SumstatsMulti import SumstatsMulti

def dump_pickle(glsumstats: 'Sumstats', path: str = "~/mysumstats.pickle", overwrite: bool = False) -> None:
    path = os.path.expanduser(path)
    glsumstats.log.write("Start to dump the Sumstats Object.")
    if overwrite==False and os.path.exists(path):
        glsumstats.log.write(" -File exists. Skipping. If you want to overwrite, please use overwrite=True.")
    else:
        with open(path, 'wb') as file:
            glsumstats.log.write(" -Dump the Sumstats Object to : ", path)
            pickle.dump(glsumstats, file)
    glsumstats.log.write("Finished dumping.")

def dump_pickle_pair(glsumpair: 'SumstatsPair', path: str = "~/mysumpair.pickle", overwrite: bool = False) -> None:
    path = os.path.expanduser(path)
    glsumpair.log.write("Start to dump the SumstatsPair Object.")
    if overwrite==False and os.path.exists(path):
        glsumpair.log.write(" -File exists. Skipping. If you want to overwrite, please use overwrite=True.")
    else:
        with open(path, 'wb') as file:
            glsumpair.log.write(" -Dump the SumstatsPair Object to : {}".format(path))
            pickle.dump(glsumpair, file)
    glsumpair.log.write("Finished dumping.")

def dump_pickle_multi(glsummulti: 'SumstatsMulti', path: str = "~/mysummulti.pickle", overwrite: bool = False) -> None:
    path = os.path.expanduser(path)
    glsummulti.log.write("Start to dump the SumstatsMulti Object.")
    if overwrite==False and os.path.exists(path):
        glsummulti.log.write(" -File exists. Skipping. If you want to overwrite, please use overwrite=True.")
    else:
        with open(path, 'wb') as file:
            glsummulti.log.write(" -Dump the SumstatsMulti Object to : {}".format(path))
            pickle.dump(glsummulti, file)
    glsummulti.log.write("Finished dumping.")

def load_pickle(path: str) -> Optional[Union['Sumstats', 'SumstatsPair', 'SumstatsMulti']]:
    """
    Load a previously saved GWASLab object from a pickle file.
    
    Automatically detects and loads Sumstats, SumstatsPair, or SumstatsMulti objects.
    
    Parameters
    ----------
    path : str
        File path to the pickle file. Supports `~` for home directory expansion.
    
    Returns
    -------
    Optional[Union[Sumstats, SumstatsPair, SumstatsMulti]]
        The loaded GWASLab object (type is automatically detected).
        Returns None if the file does not exist.
    """
    path = os.path.expanduser(path)
    
    if not os.path.exists(path):
        Log().write("File not exists : ", path)
        return None
    
    try:
        # Try loading with current module structure
        with open(path, 'rb') as file:
            obj = pickle.load(file)
            
            # Detect object type and log appropriately
            obj_type = type(obj).__name__
            if obj_type == 'SumstatsPair':
                obj.log.write("Loaded dumped SumstatsPair object created using gwaslab>=v3.4.32")
                obj.log.write("Loaded dumped SumstatsPair object from : {}".format(path))
            elif obj_type == 'SumstatsMulti':
                obj.log.write("Loaded dumped SumstatsMulti object created using gwaslab>=v3.4.32")
                obj.log.write("Loaded dumped SumstatsMulti object from : {}".format(path))
            else:
                # Default to Sumstats for backward compatibility
                obj.log.write("Loaded dumped Sumstats object created using gwaslab>=v3.4.32")
                obj.log.write("Loaded dumped Sumstats object from : {}".format(path))
            
            return obj
            
    except Exception as e:
        # Try with legacy module structure for backward compatibility
        try:
            sys.modules['gwaslab.Sumstats'] = g_Sumstats
            sys.modules['gwaslab.Log'] = g_Log
            
            with open(path, 'rb') as file:
                obj = pickle.load(file)
                
                # Detect object type and log appropriately
                obj_type = type(obj).__name__
                if obj_type == 'SumstatsPair':
                    obj.log.write("Loaded dumped SumstatsPair object created using gwaslab<v3.4.32")
                    obj.log.write("Loaded dumped SumstatsPair object from : {}".format(path))
                elif obj_type == 'SumstatsMulti':
                    obj.log.write("Loaded dumped SumstatsMulti object created using gwaslab<v3.4.32")
                    obj.log.write("Loaded dumped SumstatsMulti object from : {}".format(path))
                else:
                    # Default to Sumstats for backward compatibility
                    obj.log.write("Loaded dumped Sumstats object created using gwaslab<v3.4.32")
                    obj.log.write("Loaded dumped Sumstats object from : {}".format(path))
                
                return obj
        except Exception as e2:
            # If both attempts fail, raise the original error
            Log().write("Error loading pickle file {}: {}".format(path, str(e)))
            raise e

def load_data_from_pickle(path: str, usecols: Optional[List[str]] = None) -> pd.DataFrame:
    obj = load_pickle(path)
    if obj is None:
        raise FileNotFoundError("Pickle file not found: {}".format(path))
    data = obj.data
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