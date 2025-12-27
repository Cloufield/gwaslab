from pathlib import Path
import os
import pickle
import concurrent.futures
import threading
import multiprocessing as mp
import time
import h5py

from gwaslab.info.g_Log import Log

from platformdirs import user_cache_dir
from pysam import VariantFile

APPNAME = "gwaspipe"
APPAUTHOR = "cloufield"

CACHE_EXT = '.cache'


################################################# UTILS #################################################

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

def cache_exists(path, ref_alt_freq, category='all'):
    ''' Check if the cache file exists and contains the required data '''
    found = False
    try:
        found = is_in_h5py(path, ref_alt_freq, category)
    except Exception as e:
        pass
    return found

def is_in_h5py(path, ref_alt_freq, category='all'):
    '''
    Check if the cache file exists and contains the required data.
    Raise an exception if the cache file does not exist.
    '''
    if not path or not os.path.exists(path):
        raise Exception('Cache file not found')
    
    with h5py.File(path, 'r') as f:
        if ref_alt_freq in f.keys():
            if category in f[ref_alt_freq].keys():
                if len(f[ref_alt_freq][category].keys()) > 0:
                    return True
    return False

def load_h5py_cache(path, ref_alt_freq, category='all'):
    if not path or not os.path.exists(path):
        raise Exception('Cache file not found')
    
    if not is_in_h5py(path, ref_alt_freq, category):
        raise Exception('Cache file does not contain the required data')

    _cache = {}
    with h5py.File(path, 'r') as f:
        for v in f[ref_alt_freq][category].values():
            # iterate over chromosomes
            keys = list(v['keys'].asstr()[:])
            values = list(v['values'][:])
            chrom_cache = dict(zip(keys, values)) # Combine keys and values into a dictionary
            _cache.update(chrom_cache)
    return _cache

def build_cache(base_path, ref_alt_freq=None, threads=1, return_cache=False, filter_fn=None, category='all', log=Log(), verbose=True):
    cache_builder = CacheBuilder(base_path, ref_alt_freq=ref_alt_freq, threads=threads, log=log, verbose=verbose)
    cache_builder.start_building(filter_fn=filter_fn, category=category, set_cache=return_cache) # start_building will wait for all processes to finish building cache
    if return_cache:
        return cache_builder.get_cache()

def is_palindromic(ref, alt):
    gc = (ref=="G") & (alt=="C")
    cg = (ref=="C") & (alt=="G")
    at = (ref=="A") & (alt=="T")
    ta = (ref=="T") & (alt=="A")
    palindromic = gc | cg | at | ta 
    return palindromic

def is_indel(ref, alt):
    return len(ref) != len(alt)

def filter_fn_pi(*, ref, alt):
    return is_palindromic(ref, alt) or is_indel(ref, alt)

def filter_fn_np(*, ref, alt):
    return not is_palindromic(ref, alt)

PALINDROMIC_INDEL = 'pi' # palindromic + indel
NON_PALINDROMIC = 'np' # non-palindromic

FILTER_FN = {
    PALINDROMIC_INDEL: filter_fn_pi,
    NON_PALINDROMIC: filter_fn_np
}


################################################# CACHE MANAGERs #################################################

class CacheMainManager:
    def __init__(self, base_path, ref_alt_freq=None, category='all', filter_fn=None, threads=1, log=Log(), verbose=True):
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.category = category
        self.filter_fn = filter_fn
        self.threads = threads
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
        ''' Build and load the cache'''    
        self._cache = build_cache(
            self.base_path, ref_alt_freq=self.ref_alt_freq, threads=self.threads,
            filter_fn=self.filter_fn, category=self.category,
            return_cache=True, log=self.log, verbose=self.verbose
        )

    def load_cache(self, category=None):
        if category is None:
            category = self.category
        cache_path = self._get_cache_path()
        self._cache = load_h5py_cache(cache_path, ref_alt_freq=self.ref_alt_freq, category=category)


class CacheManager(CacheMainManager):
    def __init__(self, base_path=None, cache_loader=None, cache_process=None, ref_alt_freq=None, category='all', filter_fn=None, threads=1, log=Log(), verbose=True):
        none_value = sum([cache_loader is not None, cache_process is not None])
        assert none_value in [0, 1], 'Only one between cache_loader and cache_process should be provided'
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, category=category, filter_fn=filter_fn, threads=threads, log=log, verbose=verbose)
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
                self.load_cache()
                self.log.write('Finshed loading cache.', verbose=self.verbose)
            else:
                self.log.write(f'Start building cache from {base_path}...', verbose=self.verbose)
                self.build_cache()
                self.log.write('Finished building (and loading) cache.', verbose=self.verbose)

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


class CacheProcess(mp.Process):
    '''
    A class for managing a cache in a separate process. It is used to reduce memory consumption when the cache is very large.
    This class will load the cache in a separate process and provide methods to perform operations on the cache directly on the subprocess.
    In this way, the cache is not copied to the main process, but the operations are performed on the cache in the subprocess and only the
    input and output of the operations are communicated (i.e. copied) between the main and the subprocess.
    
    This is very useful when the cache is huge (e.g. 40GB in memory) and we want to perform operations on it based on a relatively small input
    (e.g. a "small" dataframe, where small is relative to the cache size) and the output is also relatively small.
    '''
    def __init__(self, base_path, ref_alt_freq=None, category='all', filter_fn=None, threads=1, log=Log(), verbose=True):
        super().__init__()
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.filter_fn = filter_fn
        self.category = category
        self.threads = threads
        self.log = log
        self.verbose = verbose

        self.daemon = True # When parent process exits, it will attempt to terminate all of its daemonic child processes.

        self.manager = mp.Manager()
        self.input_queue = mp.Queue()  # Queue for communication between processes
        self.result_queue = mp.Queue()
        self.result_produced = mp.Value('b', True)

        cache_path = self._get_cache_path()
        if not cache_exists(cache_path, ref_alt_freq, category):
            self.build_cache()
        else:
            if threads > 1:
                self.log.warning('[CacheProcess: since the cache already exists, the parameter threads could be set to 1 without any performance loss]', verbose=self.verbose)

    def _get_cache_path(self):
        return get_cache_path(self.base_path)
    
    def build_cache(self):
        build_cache(
            self.base_path, ref_alt_freq=self.ref_alt_freq, threads=self.threads,
            filter_fn=self.filter_fn, category=self.category,
            return_cache=False, log=self.log, verbose=self.verbose
        )

    def run(self):
        cache_path = self._get_cache_path()
        self.log.write(f'[CacheProcess: Start loading cache from {cache_path}...]', verbose=self.verbose)   
        cache = load_h5py_cache(cache_path, ref_alt_freq=self.ref_alt_freq, category=self.category)
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
    
    def apply_fn(self, fn, **kwargs):
        '''
        Apply an arbitrary function to the cache. The function should take the cache as an argument,
        and all the arguments should be passed as named arguments.
        '''
        self._call_method('apply_fn', fn, **kwargs)
        return self.result_queue.get()
    
    def cache_len(self):
        self._call_method('cache_len')
        return self.result_queue.get()

    def terminate(self):
        self._call_method("terminate")


################################################# CACHE BUILDER #################################################

class CacheBuilderOld:
    def __init__(self, ref_infer, ref_alt_freq=None, threads=1, log=Log(), verbose=True):
        self.ref_infer = ref_infer
        self.ref_alt_freq = ref_alt_freq
        self.threads = threads
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
        
        threads = self.threads
        contigs = self.get_contigs()
        
        self.cancelled = False
        self.running = True

        self.log.write(f" -Building cache on {threads} cores...", verbose=self.verbose)
        self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
        self.futures = [self.executor.submit(self.build_cache, chrom) for chrom in contigs]

    def get_contigs(self):
        vcf_reader = VariantFile(self.ref_infer, drop_samples=True)
        contigs = [v.name for v in vcf_reader.header.contigs.values()]
        vcf_reader.close()
        return contigs

    def build_cache(self, chrom):
        vcf_reader = VariantFile(self.ref_infer, drop_samples=True)
        #self.log.write(f"   -Fetching contig '{chrom}'...")
        seq = vcf_reader.fetch(chrom)
        
        first = True
        for record in seq:
            if first:
                #self.log.write(f"   -Found at least one record for contig '{chrom}'...")
                first = False
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


class CacheBuilder:
    def __init__(self, ref_infer, ref_alt_freq=None, threads=1, log=Log(), verbose=True):
        self.ref_infer = ref_infer
        self.ref_alt_freq = ref_alt_freq
        self.threads = threads
        self.log = log
        self.verbose = verbose
        
        self.running = False
        self.cache = None
        
    def get_contigs(self):
        vcf_reader = VariantFile(self.ref_infer, drop_samples=True)
        contigs = [v.name for v in vcf_reader.header.contigs.values()]
        vcf_reader.close()
        return contigs
    
    def already_built(self, category):
        cache_path = get_cache_path(self.ref_infer)
        return cache_exists(cache_path, self.ref_alt_freq, category)

    def start_building(self, filter_fn=None, category='all', set_cache=True):        
        if self.running:
            print("Cache building is already running. If you want to restart, please stop the current process first.")
            return
        
        if isinstance(filter_fn, str) and filter_fn in FILTER_FN:
            filter_fn = FILTER_FN[filter_fn]
            category = filter_fn
        elif category in FILTER_FN:
            self.log.write(f" -Using the built-in filter function for category '{category}'. filter_fn will be ignored if provided.", verbose=self.verbose)
            filter_fn = FILTER_FN[category]
        
        assert filter_fn is None or category != 'all', "If filter_fn is not None, category cannot be 'all'"
        assert filter_fn is not None or category == 'all', "If category is not 'all', filter_fn must be provided"

        if self.already_built(category=category):
            # TODO: we should probably improve the checking logic, and maybe also allows to overwrite the cache
            self.log.write(f"Cache for category '{category}' and ref_alt_freq {self.ref_alt_freq} already exists. Skipping cache building", verbose=self.verbose)
            return
        
        threads = max(self.threads-1, 1) # leave one core for the watcher process
        contigs = self.get_contigs()
        
        self.running = True

        self.log.write(f" -Building cache for category '{category}' on {threads} cores...", verbose=self.verbose)

        manager = mp.Manager()
        queue = manager.Queue()
        jobs = []
        
        # Start a watcher process to handle the output of each subprocess.
        # The watcher will write the cache to the file as soon as it receives the output from the subprocess, in a safe way.
        watcher = mp.Process(target=self.handle_output, args=(queue,))
        watcher.daemon = True
        watcher.start()
        
        with mp.Pool(threads) as pool:
            for chrom in contigs:
                job = pool.apply_async(self.build_cache, args=(chrom, queue), kwds={'filter_fn': filter_fn, 'category': category})
                jobs.append(job)
            # Pool will automatically close and join when exiting the context

        queue.put('kill') # send a signal to the watcher process to stop
        watcher.join()

        if set_cache:
            self.cache = {}
            for job in jobs:
                self.cache.update(job.get()['cache'])
    
        self.running = False
    
    def build_cache(self, chrom, queue, filter_fn=None, category='all'):       
        assert filter_fn is None or category != 'all', "If filter_fn is not None, category cannot be 'all'"

        inner_cache = {}
        ref_alt_freq = self.ref_alt_freq

        vcf_reader = VariantFile(self.ref_infer, drop_samples=True)
        #self.log.write(f"   -Fetching contig '{chrom}'...", verbose=self.verbose)
        seq = vcf_reader.fetch(chrom)

        for record in seq:
            for alt in record.alts:
                if filter_fn is None or filter_fn(ref=record.ref, alt=alt):
                    key = f"{record.chrom}:{record.pos}:{record.ref}:{alt}"
                    value = record.info[ref_alt_freq][0]
                    inner_cache[key] = value

        vcf_reader.close()
        
        result = {}
        result['chrom'] = chrom
        result['ref_alt_freq'] = ref_alt_freq
        result['category'] = category
        result['cache'] = inner_cache
        queue.put(result)
        return result
    
    def handle_output(self, queue):
        ''' Function that monitors a queue and writes the cache to a file as soon as it receives the output from the subprocess.'''
        first = True
        m = queue.get() # wait for the first message, to avoid creating an empty cache file
        
        if m != 'kill':
            cache_path = get_write_path(self.ref_infer)
            with h5py.File(cache_path, mode='a') as f:
                while True:
                    if first:
                        first = False
                    else:
                        m = queue.get()
                        
                    if m == 'kill':
                        break

                    result = m
                    cache = result['cache']
                    if cache is not None and len(cache) > 0:
                        main_group = f.require_group(result['ref_alt_freq'])
                        sub_group = main_group.require_group(result['category'])
                        chrom_group = sub_group.require_group(str(result['chrom']))

                        keys_list = list(cache.keys())
                        max_len = len(max(keys_list, key=len))
                        #self.log.write(f"Writing {result['ref_alt_freq']}, {result['category']}, {str(result['chrom'])}\n")
                        keys_dataset = chrom_group.create_dataset('keys', data=keys_list, dtype=f'S{max_len}', compression="gzip", compression_opts=4)
                        values_dataset = chrom_group.create_dataset('values', data=list(cache.values()), dtype='f', compression="gzip", compression_opts=4)

    def get_cache(self):
        return self.cache


################################################# CACHE LOADERs #################################################
# Classes for loading the cache in a separate thread or process in the background while the main process is running.
# However, right now, the most efficient way to load the cache and perform operations on it is to use the CacheProcess class.

class CacheLoader:
    def __new__(cls, *args, **kwargs):
        if cls is CacheLoader:
            raise TypeError(f"You are trying to instantiate an abstract class {cls.__name__}. Please use a concrete subclass.")
        return super().__new__(cls)
    
    def __init__(self, base_path, ref_alt_freq=None, category='all', filter_fn=None, threads=1, log=Log(), verbose=True):
        self.base_path = base_path
        self.ref_alt_freq = ref_alt_freq
        self.category = category
        self.filter_fn = filter_fn
        self.threads = threads
        self.log = log
        self.verbose = verbose

    def _get_cache_path(self):
        return get_cache_path(self.base_path)
    
    def build_cache(self):       
        self.cache = build_cache(
            self.base_path, ref_alt_freq=self.ref_alt_freq, threads=self.threads,
            filter_fn=self.filter_fn, category=self.category,
            return_cache=True, log=self.log, verbose=self.verbose
        )

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
    def __init__(self, base_path, ref_alt_freq=None, category='all', filter_fn=None, threads=1, log=Log(), verbose=True):
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, category=category, filter_fn=filter_fn, threads=threads, log=log, verbose=verbose)
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
        
        if not cache_exists(cache_path, self.ref_alt_freq, self.category):
            self.log.write("Cache does not exist. Start building (and loading) cache...", verbose=self.verbose)
            self.build_cache() # this will also load the cache
        else:
            self.running = True
            self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
            self.future = self.executor.submit(self.load_cache)

    def load_cache(self):
        cache_path = self._get_cache_path()
        self.log.write(f'[Start loading cache from {cache_path}...]', verbose=self.verbose)     
        self.cache = load_h5py_cache(cache_path, ref_alt_freq=self.ref_alt_freq, category=self.category)
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


def _load_cache_process(path, ref_alt_freq, category, cache):
    #start = time.time()
    local_cache = load_h5py_cache(path, ref_alt_freq=ref_alt_freq, category=category)
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
    def __init__(self, base_path, ref_alt_freq=None, category='all', filter_fn=None, threads=1, log=Log(), verbose=True):
        super().__init__(base_path, ref_alt_freq=ref_alt_freq, category=category, filter_fn=filter_fn, threads=threads, log=log, verbose=verbose)
        self.manager = mp.Manager()
        self.cache = self.manager.dict()
        self.running = False
        self.process = None

    def start_loading(self):
        if self.running:
            print("Cache loading is already running. If you want to restart, please stop the current process first.")
            return
        
        cache_path = self._get_cache_path()
        
        if not cache_exists(cache_path, self.ref_alt_freq, self.category):
            self.log.write("Cache does not exist. Start building (and loading) cache...", verbose=self.verbose)
            self.build_cache() # this will also load the cache
        else:
            self.running = True
            self.process = mp.Process(target=_load_cache_process, args=(cache_path, self.ref_alt_freq, self.filter_fn, self.cache))
            self.process.start()
        
    def get_cache(self):
        if self.running:
            self.process.join()  # Wait for cache loading process to finish
            self.running = False
        return self.cache