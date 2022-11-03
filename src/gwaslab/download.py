from os import path
import json
import requests
import shutil
from gwaslab.Log import Log
import gzip
import os

def check_available_ref(log=Log()):
    ref_path = path.dirname(__file__) + '/data/reference.json'
    dicts = json.load(open(ref_path))
    if dicts is not None:
        for key,value in dicts.items():
            log.write(key," : ",value)
        return dicts
    else:
        log.write("No available reference files.")

def check_downloaded_ref(log=Log()):
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        with open(config_path, 'w') as f:
            f.write("")
            log.write("Config file not exists...recreated.")
    else:
        try:
            dicts = json.load(open(config_path))
            for key,value in dicts.items():
                to_remove=[]
                if not path.exists(value):
                    to_remove.append(key)
                else:
                    log.write(key," : ",value)
            for i in to_remove:
                dicts.pop(i, None)
            with open(config_path, 'w') as f:
                json.dump(dicts,f,indent=4)
            return dicts
        except:
            log.write("No records in config file.")

def get_ref_path(name):
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        log.write("Config file not exists...")
    else:
        try:
            dicts = json.load(open(config_path))
            if path.exists(dicts[name]):
                return dicts[name]
            else:
                log.write("File not exist.")
        except:
            log.write("No records in config file. Please download first.")

def download_ref(name,directory=None,local_filename=None,log=Log()):
    from_dropbox=0
    ref_path =  path.dirname(__file__) + '/data/reference.json'
    dicts = json.load(open(ref_path))
    if name in dicts.keys():
        url = dicts[name]
        log.write("Start to download ",name," ...")
        if directory is None:
            directory = path.dirname(__file__) + '/data/'
        if local_filename is None:
            local_filename = url.split('/')[-1]
        
        full_path = directory + local_filename
        log.write(" -Downloading to:",full_path)
        if local_filename.endswith("?dl=1"):
            from_dropbox=1
            local_filename = local_filename[:-5]
            if local_filename.endswith("vcf.gz"):
                tbi_url = full_path[:-5] + ".tbi?dl=1"
        full_path = download_file(url, directory=directory, local_filename=local_filename)
        update_record(name,full_path)
        if full_path.endswith("vcf.gz"):
                if from_dropbox==0:
                    tbi_url = url+".tbi"
                try:
                    download_file(tbi_url, directory=directory, local_filename=local_filename+".tbi")
                    update_record(name+"_tbi",full_path+ ".tbi")
                except:
                    pass
        if full_path.endswith("fa.gz") or full_path.endswith("fasta.gz"):
            log.write(" -Gunzipping :",full_path)
            try:
                with gzip.open(full_path, 'rb') as f_in:
                    with open(full_path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        update_record(name,full_path[:-3])
            except:
                    pass
        log.write(" -Downloaded ",name," successfully!")
    else:
        log.write(name," is not available. Please use check_available_ref() to check available reference files.")

def update_record(key,value,log=Log()):
    config_path =  path.dirname(__file__) + '/data/config.json'
    try:
        with open(config_path,'r') as f:
            dict=json.load(f)
    except:
        dict={}
    dict[key] = value
    with open(config_path, 'w') as f:
        json.dump(dict,f,indent=4)

def remove_file(name,log=Log()):
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        log.write("Config file not exists...")
    else:
        try:
            dicts = json.load(open(config_path))
            if path.exists(dicts[name]):
                os.remove(dicts[name])
                log.write("Removed :" , dicts[name])
                check_downloaded_ref()
            else:
                log.write("File not exist.")
        except:
            log.write("No records in config file. Please download first.")

def download_file(url, directory=None,local_filename=None,log=Log()):
    if directory is None:
        directory = path.dirname(__file__) + '/data/'
    if local_filename is None:
        local_filename = url.split('/')[-1]
    full_path = directory + local_filename
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        with open(full_path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)     
    return full_path