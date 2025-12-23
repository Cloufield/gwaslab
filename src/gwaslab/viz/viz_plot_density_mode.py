import numpy as np
import pandas as pd

from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_plot_density import _process_density


def b_mode_setup(sumstats, mode, bwindowsizekb, chrom, pos, log=Log(), verbose=True):
    sumstats, bmean, bmedian = _process_density(
        sumstats=sumstats,
        mode=mode,
        bwindowsizekb=bwindowsizekb,
        chrom=chrom,
        pos=pos,
        verbose=verbose,
        log=log,
    )
    return sumstats, bmean, bmedian


def b_scaled_threshold(anno_sig_level):
    return anno_sig_level


def b_annotation_log(to_annotate, snpid, log=Log(), verbose=True):
    if (to_annotate.empty is not True):
        for index, row in to_annotate.iterrows():
            log.write(
                " -Annotated {} with {} at density {}".format(
                    row[snpid], row["Annotation"], row["DENSITY"]
                ),
                verbose=verbose,
            )
