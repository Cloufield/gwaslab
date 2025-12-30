"""
Input type handling utilities for GWASLab functions.

Provides helpers to handle both DataFrame and Sumstats object inputs.
"""

from typing import TYPE_CHECKING, Union
import pandas as pd

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _get_id_column(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame]) -> str:
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
