import matplotlib.pyplot as plt
from gwaslab.Log import Log
import time
import os.path

def save_figure(fig, save, keyword, saveargs=None, log = Log(), verbose=True):
    if saveargs is None:
        saveargs = {}
    if save:
        if verbose: log.write("Saving plot:")
        if save==True:
            default_path = get_default_path(keyword)
            fig.savefig(default_path, bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ default_path + " successfully!" )
        else:
            if os.path.exists(save):
                fig.savefig(save,bbox_inches="tight",**saveargs)
                log.write(" -Saved to "+ save + " successfully! (overwrite)" )
            else:
                fig.savefig(save,bbox_inches="tight",**saveargs)
                log.write(" -Saved to "+ save + " successfully!" )
    else:
        log.write(" -Skip saving figures!" )

def get_default_path(keyword,fmt="png"):
    path_dictionary = { 
                        "m":"manhattan",
                        "qq":"qq",
                        "mqq":"mqq",
                        "qqm":"qqm",
                        "b":"brisbane",
                        "r":"regional",
                        "trumpet_b":"trumpet_binary",
                        "trumpet_q":"trumpet_quant",
                        "power_b":"power_binary",
                        "power_q":"power_quant",
                        "power_xb":"power_x_binary",
                        "power_xq":"power_x_quant",
                        "ldscrg":"ldscrg_heatmap",
                        "miami":"miami",
                        "esc":"effect_size_comparision",
                        "afc":"allele_frequency_comparision"
                        }
    prefix = path_dictionary[keyword]
    count = 1
    for i in range(10000):
        file_path = "./gwaslab_{}_{}_{}.{}".format(prefix, time.strftime('%Y%m%d'),count,fmt) 
        if os.path.exists(file_path):
            count+=1
        else:
            break
    return file_path