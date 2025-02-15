import matplotlib.pyplot as plt
from gwaslab.g_Log import Log
import time
import os.path

def save_figure(fig, save, keyword, save_args=None, log = Log(), verbose=True):
    log.write("Start to save figure..." ,verbose=verbose)
    if save_args is None:
        save_args = {}
        
    if save:
        if save==True:
            default_path = get_default_path(keyword)
            fig.savefig(default_path, bbox_inches="tight",**save_args)
            log.write(" -Saved to "+ default_path + " successfully!" ,verbose=verbose)
        else:
            if save[-3:]=="pdf":
                if os.path.exists(save):
                    fig.savefig(save, **save_args)
                    log.write(" -Saved to "+ save + " successfully! (pdf, overwrite)" ,verbose=verbose)
                else:
                    fig.savefig(save, **save_args)
                    log.write(" -Saved to "+ save + " successfully! (pdf)" ,verbose=verbose)
            else:
                if os.path.exists(save):
                    fig.savefig(save,bbox_inches="tight",**save_args)
                    log.write(" -Saved to "+ save + " successfully! (overwrite)" ,verbose=verbose)
                else:
                    fig.savefig(save,bbox_inches="tight",**save_args)
                    log.write(" -Saved to "+ save + " successfully!" ,verbose=verbose)
    else:
        log.write(" -Skip saving figure!" ,verbose=verbose)
    log.write("Finished saving figure..." ,verbose=verbose)

def get_default_path(keyword,fmt="png"):
    path_dictionary = { 
                        "m":"manhattan",
                        "qq":"qq",
                        "mqq":"mqq",
                        "qqm":"qqm",
                        "b":"brisbane",
                        "r":"regional",
                        "stacked_r":"stacked_regional",
                        "trumpet_b":"trumpet_binary",
                        "trumpet_q":"trumpet_quant",
                        "power_b":"power_binary",
                        "power_q":"power_quant",
                        "power_xb":"power_x_binary",
                        "power_xq":"power_x_quant",
                        "ldscrg":"ldscrg_heatmap",
                        "miami":"miami",
                        "esc":"effect_size_comparision",
                        "afc":"allele_frequency_comparision",
                        "gwheatmap":"genome_wide_heatmap",
                        "scatter":"scatter",
                        "forest":"forest"
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