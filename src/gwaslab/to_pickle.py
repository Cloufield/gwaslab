import pickle
import os
import gc
from gwaslab.Log import Log 

def dump_pickle(glsumstats,path="~/mysumstats.pickle",overwrite=False):
    glsumstats.log.write("Start to dump the Sumstats Object.")
    if overwrite==False and os.path.exists(path):
        glsumstats.log.write(" -File exists. Skipping. If you want to overwrite, please use overwrite=True.")
    else:
        with open(path, 'wb') as file:
            glsumstats.log.write(" -Dump the Sumstats Object to : ", path)
            pickle.dump(glsumstats, file)
    Log().write("Finished dumping.")

def load_pickle(path):
    if os.path.exists(path):
        with open(path, 'rb') as file:
            glsumstats =  pickle.load(file)
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