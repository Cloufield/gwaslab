import numpy as np
from scipy.stats import norm

def h2_obs_to_liab(h2_obs, P, K, se_obs=None):
    '''
    Adopted from ldsc
    Reference:
    Estimating Missing Heritability for Disease from Genome-wide Association Studies
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/
    
    Converts heritability on the observed scale in an ascertained sample to heritability
    on the liability scale in the population.
    Parameters
    ----------
    h2_obs : float
        Heritability on the observed scale in an ascertained sample.
    P : float in (0,1)
        Prevalence of the phenotype in the sample.
    K : float in (0,1)
        Prevalence of the phenotype in the population.
    se_obe:float
        se of h2_obs
    '''
    
    if np.isnan(P) and np.isnan(K):
        return h2_obs
    if K <= 0 or K >= 1:
        raise ValueError('K must be in the range (0,1)')
    if P <= 0 or P >= 1:
        raise ValueError('P must be in the range (0,1)')
    
    # the fraction of y that is larger than t is K
    t = norm.isf(K)
    #z the height of the normal curve at point t
    z = norm.pdf(t)
    #Equation 23
    numerator= K**2 * (1-K)**2
    denominator= P*(1-P) * z**2
    conversion_factor = numerator/denominator
    
    if se_obs:
        se_lia = conversion_factor * se_obs
    if se_obs:
        return h2_obs * conversion_factor, se_lia
    return h2_obs * conversion_factor
def h2_se_to_p(h2,se):
    z = h2 / se
    return norm.sf(abs(z))
    

def getpersnph2(sumstats,
           beta,
           af,
           verbose=True):
    if verbose: print("Calculating per-SNP heritibility by 2 * BETA * AF * (1-AF)...")
    sumstats["PerSNPh2"] = 2*sumstats[beta]**2*sumstats[af]*(1-sumstats[af])
    sumstats["PerSNPh2"]= sumstats["PerSNPh2"].round(5)
    return sumstats