import copy

def _list_func_args(func):
    return func.__code__.co_varnames

def _extract_kwargs(prefix:str, default:dict, kwargs:dict) -> dict:
    # prefix: keyword
    # default: default dict
    # kwargs: all local kwargs  + args + kwargs

    extracted = []
    extracted_single=dict()
    for key,value in kwargs.items():
        # kwargs or args
        if key=="kwargs" or key=="args":
            for key_nested,value_nested in kwargs[key].items():
                if prefix in key_nested and "arg" in key_nested:

                    if len(key_nested.split("_"))<3:
                        extracted.append(value_nested)
                    ##
                    ## prefix_arg_fontsize
                    else:
                        print(key_nested.split("_")[-1], value)
                        extracted_single[key_nested.split("_")[-1]] = value_nested
        else:
            # local kwargs
            if prefix in key and "arg" in key:
                extracted.append(value)
    if len(extracted_single.keys()) >0:
        extracted.append(extracted_single)
    merged_arg = _merge_and_sync_dic(extracted, default)
    return merged_arg

def _merge_and_sync_dic(list_of_dics:list, default:dict) -> dict:
    temp = copy.copy(default)
    for dic in list_of_dics:
        if isinstance(dic, dict):
            temp.update(dic)
    return temp

def _update_args(args=None, default_args=None):
    
    if default_args is None:
        default_args={}

    if args is None:
        # if None, return default dict
        return default_args
    else:
        # if not None, update default dict
        for key,value in args.items():
            default_args[key] = value
        return default_args

    

def _update_arg(arg=None, default_arg=None):
    if arg is None:
        # if None, return default
        return default_arg
    else:
        # if not None, return arg  
        return arg


from functools import wraps
import inspect

def resolve_overlapping_kwargs(strategy="prefer_explicit"):
    """
    Decorator that merges overlapping keyword arguments safely.
    
    Supports explicit args + one or more dicts passed via **extra dicts.

    Example:
        ax1_safe_plot(x, y, color='red', **style_dict)
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Collect and merge all dict-type kwargs (that might have been expanded from **dict)
            merged_kwargs = {}
            for key, val in list(kwargs.items()):
                # If user passes something like extra_args=<dict>, merge it directly
                if isinstance(val, dict) and key.startswith("extra_"):
                    merged_kwargs.update(val)
                    kwargs.pop(key)

            # Merge according to strategy
            if strategy == "prefer_explicit":
                combined = {**merged_kwargs, **kwargs}  # explicit wins
            elif strategy == "prefer_extra":
                combined = {**kwargs, **merged_kwargs}  # dict wins
            elif strategy == "warn":
                overlap = set(merged_kwargs) & set(kwargs)
                if overlap:
                    print(f"[Warning] Overlapping args ignored: {overlap}")
                combined = {**merged_kwargs, **kwargs}
            else:
                combined = {**merged_kwargs, **kwargs}

            return func(*args, **combined)
        return wrapper
    return decorator