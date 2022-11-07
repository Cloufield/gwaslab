import gzip
import os
from os import path
import json
import requests
import shutil
from gwaslab.Log import Log

##################################################################################
def check_available_ref(log=Log(),verbose=True):
    if verbose : log.write("Start to check available reference files...")
    ref_path = path.dirname(__file__) + '/data/reference.json'
    try:
        dicts = json.load(open(ref_path))
        if dicts is not None:
            for key,value in dicts.items():
                if verbose :log.write(" -",key," : ",value)
            return dicts
        else:
            if verbose :log.write(" -No available reference files.")
        if verbose :log.write("Finished checking available reference files...")
    except:
        pass
    return {}

def update_available_ref(log=Log()):
    url = 'https://raw.github.com/Cloufield/gwaslab/main/src/gwaslab/data/reference.json'
    log.write("Updating available_ref list from:",url)
    r = requests.get(url, allow_redirects=True)
    data_path =  path.dirname(__file__) + '/data/reference.json'
    with open(data_path, 'wb') as file:
        file.write(r.content)
    log.write("Available_ref list has been updated!")
##################################################################################

def check_downloaded_ref(log=Log()):
    log.write("Start to check downloaded reference files...")
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        initiate_config(config_path=config_path)      
    else:
        try:
            dicts = update_config(config_path=config_path)
            return dicts
            log.write("Finished checking downloaded reference files...")
        except:
            log.write(" -No records in config file.")
    log.write("Finished checking downloaded reference files...")

##################################################################################

def get_path(name,log=Log(),verbose=True):
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        if verbose : log.write("Config file not exists...")
        if verbose : log.write("Created new config file...")
        initiate_config(config_path)      
    else:
        try:
            dicts = json.load(open(config_path))["downloaded"]
            if path.exists(dicts[name]):
                return dicts[name]
            else:
                if verbose : log.write("File not exist.")
        except:
            if verbose : log.write("No records in config file. Please download first.")
    return False

##################################################################################
def download_ref(name,
                directory=None,
                local_filename=None,
                log=Log()):
    
    from_dropbox=0
    dicts = check_available_ref(log,verbose=False)
    if name in dicts.keys():
        # get url for name
        url = dicts[name]
        log.write("Start to download ",name," ...")

        # get file local path
        if directory is None:
            directory = path.dirname(__file__) + '/data/'
        
        local_filename, from_dropbox = url_to_local_file_name(local_filename, url, from_dropbox)

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

#### helper #############################################################################################

def remove_file(name,log=Log()):
    log.write("Start to remove ",name," ...")
    config_path =  path.dirname(__file__) + '/data/config.json'
    if not path.exists(config_path):
        log.write("Config file not exists...")
    else:
        try:
            dicts = json.load(open(config_path))["downloaded"]
            if path.exists(dicts[name]):
                os.remove(dicts[name])
                log.write("Removed :" , dicts[name])
                check_downloaded_ref()
            else:
                log.write("File not exist.")
        except:
            log.write("No records in config file. Please download first.")

#### helper #############################################################################################
def update_record(  key,value,
                    log=Log(), config_path =  path.dirname(__file__) + '/data/config.json'):
    log.write(" -Updating record in config file...")
    try:
        with open(config_path,'r') as f:
            dict=json.load(f)
    except:
        dict={'default_directory':path.dirname(__file__) + '/data/','downloaded':{}}
    dict["downloaded"][key] = value
    with open(config_path, 'w') as f:
        json.dump(dict,f,indent=4)

def download_file(url, file_path=None):
    # download file from url to file_path
    if file_path is not None:
        with requests.get(url, stream=True) as r:
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)     
        return file_path

def url_to_local_file_name(local_filename, url, from_dropbox):
    if local_filename is None:
        # if local name not provided, grab it from url
        local_filename = url.split('/')[-1]
        
    if local_filename.endswith("?dl=1"):
        # if file are downloaded form dropbox
        # set from_dropbox indicator to 1
        from_dropbox=1
        # remove "?dl=1" suffix
        local_filename = local_filename[:-5]
    return local_filename, from_dropbox

##########################################################################################################
def check_and_download(name):
    # if file exsits , return file path
    if get_path(name) is not False:
        data_path = get_path(name)
        if path.exists(data_path):
            return data_path    
    
    # if file not exsits, download and return new path
    dir_path = get_default_directory()
    if not path.exists(dir_path):
        os.makedirs(dir_path)
    download_ref(name,directory = dir_path)
    data_path = get_path(name)
    return data_path

#### config ##############################################################################################
# config.json
# {
#  "default_directory":"",
#  "downloaded":{}
#}
#####################

def initiate_config(config_path=path.dirname(__file__) + '/data/config.json',log=Log()):
    log.write(" -Creating an empty config file... ")
    with open(config_path, 'w') as f:
        dict={'default_directory':path.dirname(__file__) + '/data/','downloaded':{}}
        json.dump(dict,f,indent=4) 
        log.write(" -Config file path:",path.dirname(__file__) + '/data/')
          

def update_config(config_path=path.dirname(__file__) + '/data/config.json',log=Log()):
    log.write(" -Updating config.json...")
    try:
        # if config exists
        dicts = json.load(open(config_path))
    except:
        # if config not exists
        initiate_config()
        dicts = json.load(open(config_path))
    # check if the ref file existm. If not, remove it from dicts.
    for key,value in dicts["downloaded"].items():
        to_remove=[]
        if not path.exists(value):
            to_remove.append(key)
        else:
            log.write(key," : ",value)
    for i in to_remove:
        dicts["downloaded"].pop(i, None)
    # write the dicts
    with open(config_path, 'w') as f:
        json.dump(dicts,f,indent=4)
    return dicts["downloaded"]

def set_default_directory(default_directory_path,config_path=path.dirname(__file__) + '/data/config.json'):
    try:
        # set default directory
        dicts = json.load(open(config_path))
        dicts["default_directory"] = default_directory_path
        with open(config_path, 'w') as f:
            json.dump(dicts,f,indent=4)
    except:
        # if no config file: create one and set default directory
        initiate_config()
        dicts = json.load(open(config_path))
        dicts["default_directory"] = default_directory_path
        with open(config_path, 'w') as f:
            json.dump(dicts,f,indent=4)

def get_default_directory(config_path=path.dirname(__file__) + '/data/config.json'):
    try:
        # get default directory
        dicts = json.load(open(config_path))
        return dicts["default_directory"]
    except:
        # if no config file: create one and get default directory
        initiate_config()
        dicts = json.load(open(config_path))
        return dicts["default_directory"]

##### format book ###################################################################################################
def update_formatbook(log=Log()):
    url = 'https://raw.github.com/Cloufield/formatbook/main/formatbook.json'
    log.write("Updating formatbook from:",url)
    r = requests.get(url, allow_redirects=True)
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    with open(data_path, 'wb') as file:
        file.write(r.content)
    book=json.load(open(data_path))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))
    log.write("Formatbook has been updated!")
         
def list_formats(log=Log()):
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    book=json.load(open(data_path))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))    

def check_format(fmt,log=Log()):
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    book=json.load(open(data_path))
    log.write("Available formats:",end="")
    for i in book[fmt].keys():
        log.write(i,end="")
    log.write("") 
    for i in book[fmt].values():
        log.write(i,end="")

########################################################################################################