from gwaslab.info.g_Log import Log
import sys

def _show_version(log=Log(), verbose=True):
    # show version when loading sumstats
    log.write("GWASLab v{} https://cloufield.github.io/gwaslab/".format(gwaslab_info()["version"]),verbose=verbose)
    log.write("(C) 2022-2026, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com",verbose=verbose)
    log.write(f"Python version: {sys.version}",verbose=verbose)

def _get_version():
    # return short version string like v3.4.33
    return "v{}".format(gwaslab_info()["version"])

def gwaslab_info():
    # version meta information
    dic={
   "version":"4.0.0",
   "release_date":"20251207"
    }
    return dic   

## External tool checks moved to gwaslab.extension
