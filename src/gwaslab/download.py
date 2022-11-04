import gzip
import os
from os import path
import json
import requests
import shutil
from gwaslab.Log import Log


def check_available_ref(log=Log()):
    log.write("Start to check available reference files...")
    ref_path = path.dirname(__file__) + '/data/reference.json'
    dicts = json.load(open(ref_path))
    if dicts is not None:
        for key,value in dicts.items():
            log.write(" -",key," : ",value)
        return dicts
    else:
        log.write(" -No available reference files.")
    log.write("Finished checking available reference files...")

def check_downloaded_ref(log=Log()):
    log.write("Start to check downloaded reference files...")
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        log.write(" -Config file does not exist.")
        with open(config_path, 'w') as f:
                dict={}
                log.write(" -Recreating configuration file.")
                json.dump(dict,f,indent=4)       
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
            log.write("Finished checking downloaded reference files...")
        except:
            log.write(" -No records in config file.")
    log.write("Finished checking downloaded reference files...")

def get_path(name,log=Log()):
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
        # get url for name
        url = dicts[name]
        log.write("Start to download ",name," ...")
        # get file local path
        if directory is None:
            directory = path.dirname(__file__) + '/data/'
        if local_filename is None:
            local_filename = url.split('/')[-1]
        
        if local_filename.endswith("?dl=1"):
            # set from_dropbox indicator to 1
            from_dropbox=1
            # remove "?dl=1" suffix
            local_filename = local_filename[:-5]
        
        local_path = directory + local_filename
        log.write(" -Downloading to:",local_path)
        # download file
        download_file(url,local_path)
        # update record in config json
        update_record(name,local_path)
        
        # if vcf.gz -> check tbi
        if local_path.endswith("vcf.gz"):
                if name+"_tbi" in dicts.keys():
                    tbi_url = dicts[name+"_tbi"]
                try:
                    download_file(tbi_url, local_path+".tbi")
                    update_record(name+"_tbi",local_path+ ".tbi")
                    log.write(" -Downloading to:",local_path+".tbi")
                except:
                    pass
        # if fasta.gz or fa.gz, decompress
        if local_path.endswith("fa.gz") or local_path.endswith("fasta.gz"):
            log.write(" -gunzip :",local_path)
            try:
                with gzip.open(local_path, 'rb') as f_in:
                    with open(local_path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        update_record(name,local_path[:-3])
            except:
                    pass

        log.write("Downloaded ",name," successfully!")
    else:
        log.write(name," is not available. Please use check_available_ref() to check available reference files.")

def update_record(key,value,log=Log()):
    log.write(" -Updating record in config file...")
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
    log.write("Start to remove ",name," ...")
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

def download_file(url, file_path=None):
    if file_path is not None:
        with requests.get(url, stream=True) as r:
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)     
        return file_path