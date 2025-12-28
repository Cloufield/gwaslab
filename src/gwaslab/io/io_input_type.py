"""
Input type handling utilities for GWASLab functions.

Provides decorators and helpers to handle both DataFrame and Sumstats object inputs
in a unified way.
"""

import pandas as pd
from functools import wraps
import inspect


def handle_sumstats_input(param_name=None, make_copy=False):
    """
    Decorator to automatically handle DataFrame or Sumstats object input.
    
    This decorator converts a specified parameter (or first parameter) from either
    a pandas DataFrame or Sumstats object to a DataFrame, making functions work
    seamlessly with both input types.
    
    Parameters
    ----------
    param_name : str or None, optional
        Name of the parameter to convert. If None, uses the first positional argument.
        Common names: 'sumstats_or_dataframe', 'insumstats_or_dataframe', 
        'data_or_dataframe', 'common_sumstats', 'sumstats', 'insumstats', 'data'
    make_copy : bool, default=False
        If True, creates a copy of the DataFrame before passing to the function.
        Useful when the function modifies the DataFrame in-place.
    
    Returns
    -------
    decorator
        A decorator function that wraps the target function.
    
    Examples
    --------
    >>> @handle_sumstats_input()
    ... def my_function(sumstats_or_dataframe, other_param):
    ...     # sumstats_or_dataframe is now guaranteed to be a DataFrame
    ...     return sumstats_or_dataframe.head()
    
    >>> @handle_sumstats_input(param_name='data', make_copy=True)
    ... def process_data(data, log=None):
    ...     # data is a DataFrame copy
    ...     data['new_col'] = 1
    ...     return data
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get function signature
            sig = inspect.signature(func)
            param_names = list(sig.parameters.keys())
            
            if not param_names:
                # No parameters, just call the function
                return func(*args, **kwargs)
            
            # Determine which parameter to convert
            if param_name is not None:
                target_param = param_name
            else:
                # Use first parameter by default
                target_param = param_names[0]
            
            # Find the input value
            input_value = None
            
            # Check if it's in positional args
            if target_param in param_names:
                param_index = param_names.index(target_param)
                if param_index < len(args):
                    input_value = args[param_index]
                elif target_param in kwargs:
                    input_value = kwargs[target_param]
            
            # Convert Sumstats to DataFrame if needed
            if input_value is not None:
                if isinstance(input_value, pd.DataFrame):
                    converted_value = input_value.copy() if make_copy else input_value
                else:
                    # Assume it's a Sumstats object
                    try:
                        converted_value = input_value.data.copy() if make_copy else input_value.data
                    except AttributeError:
                        # Not a Sumstats object, pass through as-is
                        converted_value = input_value
                
                # Replace in args or kwargs
                if target_param in param_names:
                    param_index = param_names.index(target_param)
                    if param_index < len(args):
                        # Replace in positional args
                        new_args = list(args)
                        new_args[param_index] = converted_value
                        return func(*new_args, **kwargs)
                    elif target_param in kwargs:
                        # Replace in kwargs
                        new_kwargs = kwargs.copy()
                        new_kwargs[target_param] = converted_value
                        return func(*args, **new_kwargs)
            
            # No conversion needed or parameter not found, call function as-is
            return func(*args, **kwargs)
        
        return wrapper
    return decorator


@handle_sumstats_input()
def _get_id_column(sumstats_or_dataframe):
    """
    Internal helper function to select the appropriate ID column (SNPID or rsID).
    
    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to check for ID columns.
    
    Returns
    -------
    str
        Column name to use: "SNPID" if available, otherwise "rsID".
    """
    data = sumstats_or_dataframe
    
    if "SNPID" in data.columns:
        return "SNPID"
    else:
        return "rsID"
