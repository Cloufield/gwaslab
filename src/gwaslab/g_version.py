from gwaslab.g_Log import Log
import subprocess
import os
import numpy as np 

def _show_version(log=Log()):
    # show when loading sumstats
    log.write("GWASLab v{} https://cloufield.github.io/gwaslab/".format(gwaslab_info()["version"]))
    log.write("(C) 2022-2023, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com")

def gwaslab_info():
    # for output header
    dic={
   "version":"3.4.32",
   "release_date":"20231118"
    }
    return dic   

def _checking_plink_version(v=2,log=Log()):
    if v==1:
        which_plink_script = "plink --version"
    elif v==2:
        which_plink_script = "plink2 --version" 
    output = subprocess.check_output(which_plink_script, stderr=subprocess.STDOUT, shell=True,text=True)
    log.write("   -PLINK version: {}".format(output.strip()))
    return log

def _checking_r_version(r, log):
    which_r_script = "{} --version".format(r)
    output = subprocess.check_output(which_r_script, stderr=subprocess.STDOUT, shell=True,text=True)
    log.write(" -R version: {}".format(output.strip()))
    return log

def _check_susie_version(r,log):
    rscript = 'print(packageVersion("susieR"))'
    temp_r = "_gwaslab_susie_temp_check_version_{}.R".format(np.random.randint(1, 99999999))
    with open(temp_r,"w") as file:
        file.write(rscript)
    which_susie_script = "{} {}".format(r, temp_r)
    output = subprocess.check_output(which_susie_script, stderr=subprocess.STDOUT, shell=True,text=True)
    log.write(" -SuSieR version: {}".format(output.strip()))
    os.remove(temp_r)
    return log