from typing import Dict, Optional
from gwaslab.info.g_Log import Log
import sys

def _show_version(log: Optional[Log] = None, verbose: bool = True) -> None:
    """Show version when loading sumstats."""
    if log is None:
        log = Log()
    log.write("GWASLab v{} https://cloufield.github.io/gwaslab/".format(gwaslab_info()["version"]),verbose=verbose)
    log.write("(C) 2022-2026, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com",verbose=verbose)
    log.write(f"Python version: {sys.version}",verbose=verbose)

def _get_version() -> str:
    """Return short version string like v3.4.33."""
    return "v{}".format(gwaslab_info()["version"])

def gwaslab_info() -> Dict[str, str]:
    """Return version meta information."""
    dic: Dict[str, str] = {
       "version":"4.0.2",
       "release_date":"20251230"
    }
    return dic   

## External tool checks moved to gwaslab.extension
