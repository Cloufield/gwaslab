import numpy as np
from scipy.stats import norm
from gwaslab.g_Log import Log


def _get_ess(sumstats, method="metal",log=Log(),verbose=True):
    log.write("Start to estimate effective sample size (N_EFF)...", verbose=verbose)
    if type(method) is str:
        if method =="metal":
            log.write(" - Method: {} ".format(method), verbose=verbose)
            log.write(" - Referencec: {} ".format("Willer, C. J., Li, Y., & Abecasis, G. R. (2010)"), verbose=verbose)
            log.write(" - Equation: {} ".format(" N_EFF = 4 * N_CASE * N_CONTROL / (N_CASE + N_CONTROL)"), verbose=verbose)
            # Willer, C. J., Li, Y., & Abecasis, G. R. (2010). METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics, 26(17), 2190-2191.
            sumstats["N_EFF"] =  4 / (1/sumstats["N_CASE"] + 1/sumstats["N_CONTROL"])
    else:
        sumstats["N_EFF"] =  method
    log.write("Finished estimating effective sample size (N_EFF)...", verbose=verbose)
    return sumstats