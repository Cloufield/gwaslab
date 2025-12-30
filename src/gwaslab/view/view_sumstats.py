from typing import TYPE_CHECKING, Union, Optional
import pandas as pd

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _view_sumstats(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], expr: Optional[str] = None) -> pd.DataFrame:
    """
    View the sumstats dataframe, optionally filtering by an expression.

    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to view.
    expr : str, optional
        A query expression string to filter the dataframe (e.g., 'P < 5e-8').
        If None, returns the original dataframe.

    Returns
    -------
    pd.DataFrame
        The filtered or original dataframe.
    """
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    if expr:
        return sumstats.query(expr)
    return sumstats

#def _head(sumstats, n=5, **kwargs):
#    """Display the first n rows of the sumstats dataframe."""
#    return sumstats.head(n)

#def _tail(sumstats, n=5, **kwargs):
#    """Display the last n rows of the sumstats dataframe."""
#    return sumstats.tail(n)

#def _random(sumstats, n=5, **kwargs):
#    """Display n random rows from the sumstats dataframe."""
#    return sumstats.sample(n=min(n, len(sumstats)), random_state=kwargs.get('random_state', None))

#def _info(sumstats, **kwargs):
#    """Display information about the sumstats dataframe."""
#    return sumstats.info()

#def _describe(sumstats, **kwargs):
#    """Display descriptive statistics of the sumstats dataframe."""
#    return sumstats.describe()