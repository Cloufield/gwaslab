from gwaslab.Log import Log

def _show_version(log=Log()):
    # show when loading sumstats
    log.write("GWASLab version 3.4.1 https://cloufield.github.io/gwaslab/")
    log.write("(C) 2022-2023, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com")

def gwaslab_info():
    # for output header
    dic={
   "version":"3.4.1",
   "release_date":"20230328"
    }
    return dic   