"""Array coercion helpers shared by orchestration layers."""

from __future__ import annotations

from typing import Callable, TypeVar, Union

import numpy as np
import pandas as pd

T = TypeVar("T")


def as_float_array(x: Union[float, np.ndarray, pd.Series]) -> np.ndarray:
    """Coerce input to a float numpy array."""
    if isinstance(x, pd.Series):
        return np.asarray(x, dtype=np.float64)
    arr = np.asarray(x, dtype=np.float64)
    return arr


def map_array_func(
    func: Callable[..., np.ndarray],
    *args: Union[float, np.ndarray, pd.Series],
    **kwargs,
) -> Union[float, np.ndarray, pd.Series]:
    """
    Apply a numpy-only ``func`` while preserving pandas Series index when present.

    Parameters
    ----------
    func : callable
        Function accepting numpy arrays and returning numpy array or scalar.
    *args
        Positional arguments that may be Series, ndarray, or scalar.
    **kwargs
        Keyword arguments forwarded to ``func``.

    Returns
    -------
    float, numpy.ndarray, or pandas.Series
        Same container style as the first Series argument, if any.
    """
    series_arg = next((a for a in args if isinstance(a, pd.Series)), None)
    np_args = tuple(as_float_array(a) if not isinstance(a, pd.Series) else as_float_array(a) for a in args)
    result = func(*np_args, **kwargs)
    if series_arg is not None:
        index = series_arg.index
        if np.ndim(result) == 0:
            return float(result)
        return pd.Series(result, index=index, dtype=np.float64)
    if all(isinstance(a, (int, float)) for a in args) and np.ndim(result) == 0:
        return float(result)
    return result
