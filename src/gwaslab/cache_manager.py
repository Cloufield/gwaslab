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


def get_cache_path(base_path):
    cache_filename = str(Path(base_path).stem) + '.pickle'
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
    cache_filename = str(Path(base_path).stem) + '.pickle'
    if os.access(os.path.dirname(base_path), os.W_OK):
        # if we have write access to the directory where the original input file is located
        return os.path.join(os.path.dirname(base_path), cache_filename)
    else:
        cache_dir = user_cache_dir(APPNAME, APPAUTHOR)
        if os.access(cache_dir, os.W_OK):
            # if we have write access to the user cache directory
            return os.path.join(cache_dir, cache_filename)
        
    raise Exception('No write access to any cache directory')


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
        cache_builder = CacheBuilder(self.base_path, ref_alt_freq=self.ref_alt_freq, n_cores=self.n_cores, log=self.log, verbose=self.verbose)
        cache_builder.start_building()
        cache_builder.save_cache(self._get_write_path())
        self._cache = cache_builder.get_cache(complete=True) # wait for all threads to finish building cache

    def load_cache(self, cache_path):
        if os.path.exists(cache_path):
            with open(self.cache_dir, 'rb') as f:
                self._cache = pickle.load(f)
        else:
            raise Exception('Cache file not found')


class CacheManager(CacheMainManager):
    def __init__(self, base_path=None, cache_loader=None, cache_process=None, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        assert sum([base_path is not None, cache_loader is not None, cache_process is not None]) == 1, 'Only one between base_path, cache_loader and cache_process should be provided'
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, n_cores=n_cores, log=log, verbose=verbose)
        self.cache_loader = cache_loader
        self.cache_process = cache_process

        if cache_loader is not None:
            assert callable(getattr(cache_loader, 'get_cache', None)), 'cache_loader must have a get_cache method'
        elif cache_process is not None:
            assert isinstance(cache_process, CacheProcess), 'cache_process must be an instance of CacheManagerProcess'
        else:
            cache_path = self._get_cache_path()
            if cache_path is not None:
                self.log(f'Start to load cache from {cache_path}...', verbose=self.verbose)
                self.load_cache(cache_path)
                self.log('Finshed loading cache.', verbose=self.verbose)
            else:
                self.log(f'Start to build cache from {base_path}...', verbose=self.verbose)
                self.build_cache()
                self.log('Finished building cache.', verbose=self.verbose)

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
                raise Exception('Cache not loaded')
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
    def __init__(self, base_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        super().__init__()
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.n_cores = n_cores
        self.log = log
        self.verbose = verbose

        self.manager = multiprocessing.Manager()
        self.input_queue = multiprocessing.Queue()  # Queue for communication between processes
        self.result_queue = multiprocessing.Queue()
        self.result_produced = multiprocessing.Value('b', True)

        if not os.path.exists(self._get_cache_path()):
            self.build_cache()

    def _get_cache_path(self):
        return get_cache_path(self.base_path)
    
    def build_cache(self):
        cache_builder = CacheBuilder(self.base_path, ref_alt_freq=self.ref_alt_freq, n_cores=self.n_cores, log=self.log, verbose=self.verbose)
        cache_builder.start_building()
        cache_builder.save_cache(self._get_write_path())

    def run(self):
        print("Start loading")
        start = time.time()
        with open(self._get_cache_path, 'rb') as handle:
            cache = pickle.load(handle)
        print(f"Done loading in {time.time() - start} seconds")
        
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


class CacheBuilder:
    def __init__(self, ref_infer, save_path, ref_alt_freq=None, n_cores=1, log=Log(), verbose=True):
        self.ref_infer = ref_infer
        self.save_path = save_path
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

        self.log(f" -Building cache on {n_cores} cores...", verbose=self.verbose)
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

    def save_cache(self):
        cache = self.get_cache(complete=True)
        self.log(f' -Saving cache to {self.save_path}', verbose=self.verbose)
        with open(self.save_path, 'wb') as f:
            pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)
        self.log(' -Cache saved', verbose=self.verbose)