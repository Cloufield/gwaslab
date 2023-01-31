from gwaslab.Log import Log

def _show_version(log=Log()):
    # show when loading sumstats
    log.write("GWASLab version 3.3.23 https://cloufield.github.io/gwaslab/")
    log.write("(C) 2022-2023, Yunye He, MIT License, gwaslab@gmail.com")

def gwaslab_info():
    # for output header
    dic={
   "version":"3.3.23",
   "release_date":"20230131"
    }
    return dic   