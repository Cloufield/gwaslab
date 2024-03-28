import pickle
import os
import gc
from gwaslab.g_Log import Log 
import sys
from gwaslab import g_Sumstats
from gwaslab import g_Log
def dump_pickle(glsumstats,path="~/mysumstats.pickle",overwrite=False):
    glsumstats.log.write("Start to dump the Sumstats Object.")
    if overwrite==False and os.path.exists(path):
        glsumstats.log.write(" -File exists. Skipping. If you want to overwrite, please use overwrite=True.")
    else:
        with open(path, 'wb') as file:
            glsumstats.log.write(" -Dump the Sumstats Object to : ", path)
            pickle.dump(glsumstats, file)
    glsumstats.log.write("Finished dumping.")

def load_pickle(path):
    if os.path.exists(path):
        try:
            with open(path, 'rb') as file:
                glsumstats =  pickle.load(file)
                glsumstats.log.write("Loaded dumped Sumstats object created using gwaslab>=v3.4.32")    
                glsumstats.log.write("Loaded dumped Sumstats object from : ", path)
                return glsumstats
        except:
            sys.modules['gwaslab.Sumstats'] = g_Sumstats
            sys.modules['gwaslab.Log'] = g_Log
            
            with open(path, 'rb') as file:
                glsumstats =  pickle.load(file)
                glsumstats.log.write("Loaded dumped Sumstats object created using gwaslab<v3.4.32")    
                glsumstats.log.write("Loaded dumped Sumstats object from : ", path)
                return glsumstats
    else:
        Log().write("File not exists : ", path)

def load_data_from_pickle(path,usecols=None):
    data = load_pickle(path).data
    existing_cols = []
    if usecols is not None:
        for i in usecols:
            if i in data.columns:
                existing_cols.append(i)
        data = data.loc[:,existing_cols]
        gc.collect()
    return data