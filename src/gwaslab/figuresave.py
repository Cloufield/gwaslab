import matplotlib.pyplot as plt
from gwaslab.Log import Log

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
            fig.savefig(save,bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ save + " successfully!" )
    else:
        log.write(" -Skip saving figures!" )


def get_default_path(keyword):
    path_dictionary = {"trumpet":"./trumpet_plot.png",
                        "mqq":"./mqq_plot.png",
                        "regional":"./regional_plot.png",
                        "genetic_correlation":"./genetic_correlation_heatmap.png",
                        "qq":"./qq_plot.png",
                        "brisbane":"./brisbane_plot.png",
                        "miami":"./miami_plot.png",
                        "effect_size_comparision":"./effect_size_comparision_plot.png",
                        "allele_frequency_comparision":"./allele_frequency_comparision_plot.png"
                        }
    return path_dictionary[keyword]