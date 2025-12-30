from typing import Union
import scipy as sp
import numpy as np
from gwaslab.info.g_Log import Log

def wc_correct(
    beta: Union[float, np.ndarray],
    se: Union[float, np.ndarray],
    sig_level: float = 5e-8,
    log: Log = Log(),
    verbose: bool = True
) -> Union[float, np.ndarray]:
    """winner's curse correction
    Args:
		beta (float): observed beta
		se (float):  observed beta se
		sig_level (float) : significance threshold
	Returns:
		float : corrected beta

    Reference:
    	Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and confidence intervals for odds ratios in genome-wide association studies. Biostatistics, 9(4), 621-634.
    """
    
    #calculate c
    c2 = sp.stats.chi2.ppf(1-sig_level,df=1)
    c = np.sqrt(c2)
    
    # def the equation to solve
    def bias(beta_T,beta_O,se):
        z = beta_T / se
        numerator = sp.stats.norm.pdf(z - c) - sp.stats.norm.pdf(- z - c)
        denominator = sp.stats.norm.cdf(z - c) + sp.stats.norm.cdf(- z - c)
        return beta_T + se * numerator / denominator - beta_O
    
    # solve the equation using brent method
    minimum = sp.optimize.brentq(lambda x : bias(x,beta,se),a=-100,b=100, maxiter=1000)
    
    return minimum

def wc_correct_test(
    beta: Union[float, np.ndarray],
    se: Union[float, np.ndarray],
    sig_level: float = 5e-8
) -> Union[float, np.ndarray]:
    """winner's curse correction
    Args:
		beta (float): observed beta
		se (float):  observed beta se
		sig_level (float) : significance threshold
	Returns:
		float : corrected beta

    Reference:
    	Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and confidence intervals for odds ratios in genome-wide association studies. Biostatistics, 9(4), 621-634.
    """

    #calculate c
    c2 = sp.stats.chi2.ppf(1-sig_level,df=1)
    c = np.sqrt(c2)
    
    # def the equation to solve
    def bias(beta_T,beta_O,se):
        z = beta_T / se
        numerator = sp.stats.norm.pdf(z - c) - sp.stats.norm.pdf(- z - c)
        denominator = sp.stats.norm.cdf(z - c) + sp.stats.norm.cdf(- z - c)
        return  beta_T + se * numerator / denominator - beta_O
    
    # solve the equation using brent method
    minimum = sp.optimize.brentq(lambda x : bias(x,beta,se),a=-100,b=100, maxiter=1000)
    
    return minimum