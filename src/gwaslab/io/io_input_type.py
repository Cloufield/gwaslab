"""
Input type handling utilities for GWASLab functions.

Provides helpers to handle both DataFrame and Sumstats object inputs.
"""

import pandas as pd


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
    # Extract DataFrame if Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        data = sumstats_or_dataframe
    else:
        # Assume it's a Sumstats object
        try:
            data = sumstats_or_dataframe.data
        except AttributeError:
            # Not a Sumstats object, pass through as-is
            data = sumstats_or_dataframe
    
    if "SNPID" in data.columns:
        return "SNPID"
    else:
        return "rsID"
