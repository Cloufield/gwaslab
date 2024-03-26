from pathlib import Path
import os
import pickle
import concurrent.futures
import threading
import multiprocessing
import time

from gwaslab.g_Log import Log

from platformdirs import user_cache_dir
from pysam import VariantFile

APPNAME = "gwaspipe"
APPAUTHOR = "ht_diva"

CACHE_EXT = '.cache'


def get_cache_path(base_path):
    cache_filename = str(Path(base_path).stem) + CACHE_EXT
    cache_path = os.path.join(os.path.dirname(base_path), cache_filename)
    if os.path.exists(cache_path):
        return cache_path
    else:
        cache_dir = user_cache_dir(APPNAME, APPAUTHOR)
        user_cache_path = os.path.join(cache_dir, cache_filename)
        if os.path.exists(user_cache_path):
            return user_cache_path

    return None

def get_write_path(base_path):
    cache_filename = str(Path(base_path).stem) + CACHE_EXT
    if os.access(os.path.dirname(base_path), os.W_OK):
        # if we have write access to the directory where the original input file is located
        return os.path.join(os.path.dirname(base_path), cache_filename)
    else:
        cache_dir = user_cache_dir(APPNAME, APPAUTHOR)
        if os.access(cache_dir, os.W_OK):
            # if we have write access to the user cache directory
            return os.path.join(cache_dir, cache_filename)
        
    raise Exception('No write access to any cache directory')

def build_cache(base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True, return_cache=False):
    cache_builder = CacheBuilder(base_path, ref_alt_freq=ref_alt_freq, n_cores=n_cores, log=log, verbose=verbose)
    cache_builder.start_building()
    cache_builder.save_cache(get_write_path(base_path)) # save_cache will wait for all threads to finish building cache
    if return_cache:
        return cache_builder.get_cache(complete=True) # wait for all threads to finish building cache
    

################################################# CACHE MANAGERs #################################################

class CacheMainManager:
    def __init__(self, base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.n_cores = n_cores
        self.log = log
        self.verbose = verbose

    def _get_cache_path(self):
        return get_cache_path(self.base_path)
    
    def _get_write_path(self):
        if self.base_path is not None:
            return get_write_path(self.base_path)
        else:
            raise Exception('base_path is None')
        
    @property
    def cache_len(self):
        return len(self.cache)

    @property
    def cache(self):
        if not hasattr(self, '_cache'):
            raise Exception('Cache not loaded')
        return self._cache

    def build_cache(self):
        self._cache = build_cache(self.base_path, ref_alt_freq=self.ref_alt_freq, n_cores=self.n_cores, log=self.log, verbose=self.verbose, return_cache=True)

    def load_cache(self, cache_path):
        if cache_path and os.path.exists(cache_path):
            with open(cache_path, 'rb') as f:
                self._cache = pickle.load(f)
        else:
            raise Exception('Cache file not found')


class CacheManager(CacheMainManager):
    def __init__(self, base_path=None, cache_loader=None, cache_process=None, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        none_value = sum([cache_loader is not None, cache_process is not None])
        assert none_value in [0, 1], 'Only one between cache_loader and cache_process should be provided'
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, n_cores=n_cores, log=log, verbose=verbose)
        if none_value == 1:
            self.base_path = None # unset base_path if cache_loader or cache_process is provided

        self.cache_loader = cache_loader
        self.cache_process = cache_process

        if cache_loader is not None:
            assert callable(getattr(cache_loader, 'get_cache', None)), 'cache_loader must have a get_cache method'
        elif cache_process is not None:
            assert isinstance(cache_process, CacheProcess), 'cache_process must be an instance of CacheProcess'
        else:
            cache_path = self._get_cache_path()
            if cache_path is not None:
                self.log.write(f'Start loading cache from {cache_path}...', verbose=self.verbose)
                self.load_cache(cache_path)
                self.log.write('Finshed loading cache.', verbose=self.verbose)
            else:
                self.log.write(f'Start building cache from {base_path}...', verbose=self.verbose)
                self.build_cache()
                self.log.write('Finished building cache.', verbose=self.verbose)

    @property
    def cache_len(self):
        if self.cache_process is not None:
            return self.cache_process.cache_len()
        else:
            return len(self.cache)

    @property
    def cache(self):
        if self.cache_loader is not None:
            return self.cache_loader.get_cache()
        else:
            if not hasattr(self, '_cache'):
                raise Exception('Cache not loaded or class not exposing cache')
            return self._cache
        
    def apply_fn(self, fn, *args, **kwargs):
        assert 'cache' not in kwargs, "'cache' can't be inside kwargs"
        if self.cache_process is not None:
            return self.cache_process.apply_fn(fn, *args, **kwargs)
        else:
            return fn(*args, cache=self.cache, **kwargs)

    def _get_cache_path(self):
        if self.cache_loader is None and self.cache_process is None:
            return super()._get_cache_path()
        return None


class CacheProcess(multiprocessing.Process):
    '''
    A class for managing a cache in a separate process. It is used to reduce memory consumption when the cache is very large.
    This class will load the cache in a separate process and provide methods to perform operations on the cache directly on the subprocess.
    In this way, the cache is not copied to the main process, but the operations are performed on the cache in the subprocess and only the
    input and output of the operations are communicated (i.e. copied) between the main and the subprocess.
    
    This is very useful when the cache is huge (e.g. 40GB in memory) and we want to perform operations on it based on a relatively small input
    (e.g. a "small" dataframe, where small is relative to the cache size) and the output is also relatively small.
    '''
    def __init__(self, base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        super().__init__()
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.n_cores = n_cores
        self.log = log
        self.verbose = verbose

        self.daemon = True # When parent process exits, it will attempt to terminate all of its daemonic child processes.

        self.manager = multiprocessing.Manager()
        self.input_queue = multiprocessing.Queue()  # Queue for communication between processes
        self.result_queue = multiprocessing.Queue()
        self.result_produced = multiprocessing.Value('b', True)

        cache_path = self._get_cache_path()
        if cache_path is None or not os.path.exists(cache_path):
            self.build_cache()
        else:
            if n_cores > 1:
                self.log.warning('[CacheProcess: since the cache already exists, the parameter n_cores could be set to 1 without any performance loss]', verbose=self.verbose)

    def _get_cache_path(self):
        return get_cache_path(self.base_path)
    
    def build_cache(self):
        build_cache(self.base_path, ref_alt_freq=self.ref_alt_freq, n_cores=self.n_cores, log=self.log, verbose=self.verbose, return_cache=False)

    def run(self):
        cache_path = self._get_cache_path()
        self.log.write(f'[CacheProcess: Start loading cache from {cache_path}...]', verbose=self.verbose)   
        with open(cache_path, 'rb') as handle:
            cache = pickle.load(handle)
        self.log.write('[CacheProcess: Finshed loading cache.]', verbose=self.verbose)       
        
        # Continuously listen for method calls
        while True:
            method, args, kwargs = self.input_queue.get()
            if method == 'get_from_cache':
                key = args[0]
                self.result_queue.put(cache[key])
                self.result_produced.value = True
            elif method == 'apply_fn':
                assert 'cache' not in kwargs, "'cache' can't be inside kwargs"
                fn, *args = args
                result = fn(*args, cache=cache, **kwargs)
                self.result_queue.put(result)
                self.result_produced.value = True
            elif method == 'cache_len':
                self.result_queue.put(len(cache))
                self.result_produced.value = True
            elif method == "terminate":
                self.result_produced.value = True
                break

    def _call_method(self, method, *args, **kwargs):
        self.result_produced.value = False
        self.input_queue.put((method, args, kwargs))
        
        # wait until the result is produced
        while not self.result_produced.value:
            pass
        
    def get_from_cache(self, key):
        self._call_method('get_from_cache', key)
        return self.result_queue.get()
    
    def apply_fn(self, fn, *args, **kwargs):
        ''' Apply an arbitrary function to the cache. The function should take the cache as an argument.'''
        self._call_method('apply_fn', fn, *args, **kwargs)
        return self.result_queue.get()
    
    def cache_len(self):
        self._call_method('cache_len')
        return self.result_queue.get()

    def terminate(self):
        self._call_method("terminate")


################################################# CACHE BUILDER #################################################

class CacheBuilder:
    def __init__(self, ref_infer, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        self.ref_infer = ref_infer
        self.ref_alt_freq = ref_alt_freq
        self.n_cores = n_cores
        self.log = log
        self.verbose = verbose

        self.cache = {}
        self.lock = threading.Lock()  # For thread-safe cache access
        self.cancelled = False  # Flag for cancelling the cache building process
        self.running = False
        self.executor = None  # Thread pool executor
        self.futures = None  # Stores Future objects

    def start_building(self):
        if self.running:
            print("Cache building is already running. If you want to restart, please stop the current process first.")
            return
        
        n_cores = self.n_cores
        chroms = [str(i) for i in range(1, 23, 1)]

        self.cancelled = False
        self.running = True

        self.log.write(f" -Building cache on {n_cores} cores...", verbose=self.verbose)
        self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=n_cores)
        self.futures = [self.executor.submit(self.build_cache, chrom) for chrom in chroms]

    def build_cache(self, chrom):
        vcf_reader = VariantFile(self.ref_infer, drop_samples=True)
        seq = vcf_reader.fetch(chrom)

        for record in seq:
            chrom = record.chrom
            start = record.pos - 1
            end = record.pos
            cache_key = f"{chrom}:{start}:{end}"
            to_add = [record.pos, record.ref, record.alts, record.info[self.ref_alt_freq][0]]
            self.add_to_cache(cache_key, to_add)

    def stop_building(self, wait=False, verbose=False):
        if self.futures:
            self.cancelled = True
            for future in self.futures:
                future.cancel()
            self.executor.shutdown(wait=wait)  # Whether to wait for threads to finish
            self.futures = None
            self.executor = None
            self.running = False

        if verbose:
            print(f"Cache contains {len(self.get_cache())} variants")

    def add_to_cache(self, key, value):
        self.lock.acquire()
        if key in self.cache:
            self.cache[key].append(value)
        else:
            self.cache[key] = [value]
        self.lock.release()

    def get_cache(self, complete=False):
        if complete:
            concurrent.futures.wait(self.futures)

        self.lock.acquire()
        cache = self.cache
        self.lock.release()
        return cache
    
    def reset_cache(self):
        self.lock.acquire()
        self.cache = {}
        self.lock.release()

    def save_cache(self, save_path):
        cache = self.get_cache(complete=True)
        self.log.write(f' -Saving cache to {save_path}', verbose=self.verbose)
        with open(save_path, 'wb') as f:
            pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)
        self.log.write(' -Cache saved', verbose=self.verbose)


################################################# CACHE LOADERs #################################################
# Classes for loading the cache in a separate thread or process in the background while the main process is running.
# However, right now, the most efficient way to load the cache and perform operations on it is to use the CacheProcess class.

class CacheLoader:
    def __new__(cls, *args, **kwargs):
        if cls is CacheLoader:
            raise TypeError(f"You are trying to instantiate an abstract class {cls.__name__}. Please use a concrete subclass.")
        return super().__new__(cls)
    
    def __init__(self, base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.n_cores = n_cores
        self.log = log
        self.verbose = verbose

    def _get_cache_path(self):
        return get_cache_path(self.base_path)
    
    def build_cache(self):
        self.cache = build_cache(self.base_path, ref_alt_freq=self.ref_alt_freq, n_cores=self.n_cores, log=self.log, verbose=self.verbose, return_cache=True)

    def add_to_cache(self, key, value):
        self.cache[key] = value

    def get_cache(self):
        return self.cache
    
    def reset_cache(self):
        self.cache = {}


class CacheLoaderThread(CacheLoader):
    '''
    A class for loading a cache in a separate thread. It is used to load the cache in the background while the main process is running.
    
    In theory, this should be the best and simplest approach to directly load the cache in the same process as the main process, without further 
    copying the cache to the main process. However, due to the GIL (Global Interpreter Lock) in Python, this approach is not efficient and
    it slows down the main process.
    '''
    def __init__(self, base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, n_cores=n_cores, log=log, verbose=verbose)

        self.cache = {}
        self.lock = threading.Lock()  # For thread-safe cache access
        self.running = False
        self.executor = None  # Thread pool executor
        self.future = None  # Stores Future objects

    def start_loading(self):
        if self.running:
            print("Cache loading is already running. If you want to restart, please stop the current process first.")
            return
        
        cache_path = self._get_cache_path()
        
        if cache_path and os.path.exists(cache_path) is False:
            self.log.write("No cache file found. Start building (and loading) cache...", verbose=self.verbose)
            self.build_cache() # this will also load the cache
        else:
            self.running = True
            self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
            self.future = self.executor.submit(self.load_cache)

    def load_cache(self):
        cache_path = self._get_cache_path()
        self.log.write(f'[Start loading cache from {cache_path}...]', verbose=self.verbose)     
        with open(cache_path, 'rb') as handle:
            self.cache = pickle.load(handle)
        self.log.write('[Finshed loading cache.]', verbose=self.verbose)

        self.future.cancel()
        self.executor.shutdown(wait=False)
        self.executor = None
        self.future = None
        self.running = False     

    def get_cache(self):
        if self.future is not None:
            self.future.result()  # Ensure loading is finished before accessing the cache
        return self.cache


def _load_cache_process(path, cache):
    #start = time.time()
    with open(path, 'rb') as handle:
        local_cache = pickle.load(handle)
    #print(f" ********* DONE LOADING local in {time.time() - start} seconds *********")

    #start = time.time()
    cache.update(local_cache)
    #print(f" ********* DONE COPYING shared in {time.time() - start} seconds *********")
    del local_cache

class CacheLoaderProcess(CacheLoader):
    '''
    A class for loading a cache in a separate process. It is used to load the cache in the background while the main process is running.

    Unlike CacheLoaderThread, this class is more efficient because it loads the cache in a separate process, which is not affected by the GIL.
    However, a lot of memory and time is wasted in copying the cache from the subprocess to the main process.
    '''
    def __init__(self, base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, n_cores=n_cores, log=log, verbose=verbose)
        self.manager = multiprocessing.Manager()
        self.cache = self.manager.dict()
        self.running = False
        self.process = None

    def start_loading(self):
        if self.running:
            print("Cache loading is already running. If you want to restart, please stop the current process first.")
            return
        
        cache_path = self._get_cache_path()
        
        if cache_path and os.path.exists(cache_path) is False:
            self.build_cache() # this will also load the cache
        else:
            self.running = True
            self.process = multiprocessing.Process(target=_load_cache_process, args=(cache_path, self.cache))
            self.process.start()
        
    def get_cache(self):
        if self.running:
            self.process.join()  # Wait for cache loading process to finish
            self.running = False
        return self.cache