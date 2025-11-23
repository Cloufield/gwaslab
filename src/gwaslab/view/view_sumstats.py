import pandas as pd

def _head(sumstats, n=5, **kwargs):
    """Display the first n rows of the sumstats dataframe."""
    return sumstats.head(n)

def _tail(sumstats, n=5, **kwargs):
    """Display the last n rows of the sumstats dataframe."""
    return sumstats.tail(n)

#def _random(sumstats, n=5, **kwargs):
#    """Display n random rows from the sumstats dataframe."""
#    return sumstats.sample(n=min(n, len(sumstats)), random_state=kwargs.get('random_state', None))

#def _info(sumstats, **kwargs):
#    """Display information about the sumstats dataframe."""
#    return sumstats.info()

#def _describe(sumstats, **kwargs):
#    """Display descriptive statistics of the sumstats dataframe."""
#    return sumstats.describe()