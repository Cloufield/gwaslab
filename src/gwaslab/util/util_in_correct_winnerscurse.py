from typing import Union
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.algorithm.finemapping.winners_curse import winners_curse_correct


def wc_correct(
    beta: Union[float, np.ndarray],
    se: Union[float, np.ndarray],
    sig_level: float = 5e-8,
    log: Log = Log(),
    verbose: bool = True
) -> Union[float, np.ndarray]:
    """Apply winner's curse correction (delegates to ``gwaslab.algorithm``)."""
    result = winners_curse_correct(beta, se, sig_level=sig_level)
    return result


def wc_correct_test(
    beta: Union[float, np.ndarray],
    se: Union[float, np.ndarray],
    sig_level: float = 5e-8
) -> Union[float, np.ndarray]:
    """Apply winner's curse correction (delegates to ``gwaslab.algorithm``)."""
    return winners_curse_correct(beta, se, sig_level=sig_level)
