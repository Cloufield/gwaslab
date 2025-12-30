from typing import Callable, Any, Dict, List, Optional, Iterable
import copy
import pandas as pd

def _list_func_args(func: Callable[..., Any]) -> tuple:
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

def _merge_and_sync_dic(list_of_dics: List[Dict[str, Any]], default: Dict[str, Any]) -> Dict[str, Any]:
    temp = copy.copy(default)
    for dic in list_of_dics:
        if isinstance(dic, dict):
            temp.update(dic)
    return temp

def _update_kwargs(
    args: Optional[Dict[str, Any]] = None,
    default_args: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    
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

    

def _update_arg(arg: Optional[Any] = None, default_arg: Optional[Any] = None) -> Any:
    if arg is None:
        # if None, return default
        return default_arg
    else:
        # if not None, return arg  
        return arg


from functools import wraps
import inspect

import inspect

import inspect
from functools import wraps

def resolve_overlapping_kwargs(strategy: str = "prefer_explicit", verbose: bool = False) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator factory for any function."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            sig = inspect.signature(func)
            bound_args = {}
            for i, (name, param) in enumerate(sig.parameters.items()):
                if i < len(args):
                    bound_args[name] = args[i]

            overlap = set(bound_args.keys()) & set(kwargs.keys())
            if overlap:
                if verbose:
                    print(f"[resolve_overlap] overlapping: {overlap}")
                if strategy == "prefer_explicit":
                    for key in overlap:
                        kwargs.pop(key, None)
                elif strategy == "prefer_kwargs":
                    for key in overlap:
                        bound_args[key] = kwargs.pop(key)

            return func(**bound_args, **kwargs)
        return wrapper
    return decorator

def remove_overlapping_kwargs(
    kwargs_dict: Dict[str, Any],
    protected_keys: Iterable[str]
) -> Dict[str, Any]:
    """
    Return a new dict with keys in `protected_keys` removed.

    Parameters
    ----------
    kwargs_dict : dict
        Original kwargs dictionary.
    protected_keys : Iterable[str]
        Keys that should NOT appear in the returned dict.

    Returns
    -------
    dict
        A cleaned dictionary with protected_keys removed.
    """
    return {k: v for k, v in kwargs_dict.items() if k not in protected_keys}


def normalize_series_inputs(keys: Optional[List[str]] = None) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            sig = inspect.signature(func)
            bound = sig.bind_partial(*args, **kwargs)
            bound.apply_defaults()
            for name, value in list(bound.arguments.items()):
                if (keys is None or name in keys) and isinstance(value, pd.Series):
                    bound.arguments[name] = value.tolist()
            # Call function with all arguments as keyword arguments
            # This preserves the original behavior while allowing Series conversion
            return func(**bound.arguments)
        return wrapper
    return decorator
