import pandas as pd
from gwaslab.g_Log import Log

def _read_pipcs(data, output_prefix):
    pipcs = pd.read_csv("{}.pipcs".format(output_prefix))
    pipcs = _merge_chrpos(data,pipcs)
    return pipcs

def _merge_chrpos(data,pipcs):
    df = pd.merge(pipcs, data,on="SNPID",how="left")
    return df