from typing import Callable, Any
import functools
import matplotlib.pyplot as plt

def add_doc(src: Callable[..., Any]) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    def deco(func):
        src_doc = src.__doc__ or ""
        func_doc = func.__doc__ or ""
        func.__doc__ = func_doc + src_doc
        return func
    return deco

def suppress_display(func: Callable[..., Any]) -> Callable[..., Any]:
    @functools.wraps(func)
    def _wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        try:
            from matplotlib.figure import Figure
            if isinstance(result, Figure):
                plt.close(result)
            elif isinstance(result, tuple) and len(result) > 0 and isinstance(result[0], Figure):
                plt.close(result[0])
            return result
        except Exception:
            return result
    return _wrapper
