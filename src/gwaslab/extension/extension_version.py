from typing import Tuple, Dict, Any, Optional
from gwaslab.info.g_Log import Log
import shutil
import subprocess
import os
import numpy as np

def _check_tool_availability(
    tools: Tuple[str, ...] = ("tabix","bcftools"),
    log: Log = Log(),
    verbose: bool = True
) -> Dict[str, Dict[str, Any]]:
    results = {}
    for tool in tools:
        path = shutil.which(tool)
        if path is None:
            log.write(f" -{tool} not found in PATH", verbose=verbose)
            results[tool] = {"available": False, "path": None, "version": None}
            continue
        try:
            output = subprocess.check_output(f"{tool} --version", stderr=subprocess.STDOUT, shell=True, text=True)
            line = output.strip().splitlines()[0] if output else ""
            log.write(f" -{tool} version: {line}", verbose=verbose)
            results[tool] = {"available": True, "path": path, "version": line}
        except Exception:
            log.write(f" -{tool} found: {path}", verbose=verbose)
            results[tool] = {"available": True, "path": path, "version": None}
    return results

def _checking_plink_version(
    plink: Optional[str] = None,
    plink2: Optional[str] = None,
    log: Log = Log(),
    verbose: bool = True
) -> Log:
    if plink is not None:
        which_plink_script = f"{plink} --version"
    elif plink2 is not None:
        which_plink_script = f"{plink2}  --version"
    output = subprocess.check_output(which_plink_script, stderr=subprocess.STDOUT, shell=True, text=True)
    log.write(f" -PLINK version: {output.strip()}")
    return log

def _checking_r_version(r: str, log: Log = Log(), verbose: bool = True) -> Log:
    which_r_script = f"{r} --version"
    output = subprocess.check_output(which_r_script, stderr=subprocess.STDOUT, shell=True, text=True)
    log.write(f" -R version: {output.strip()}",verbose=verbose)
    return log

def _check_susie_version(r: str, log: Log = Log(), verbose: bool = True) -> Log:
    rscript = 'print(packageVersion("susieR"))'
    temp_r = f"_gwaslab_susie_temp_check_version_{np.random.randint(1, 99999999)}.R"
    with open(temp_r,"w") as file:
        file.write(rscript)
    which_susie_script = f"{r} {temp_r}"
    output = subprocess.check_output(which_susie_script, stderr=subprocess.STDOUT, shell=True, text=True)
    log.write(f" -SuSieR version: {output.strip()}",verbose=verbose)
    os.remove(temp_r)
    return log
