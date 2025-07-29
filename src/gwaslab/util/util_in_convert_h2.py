import numpy as np
from scipy.stats import norm
from gwaslab.g_Log import Log

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
    

def _get_per_snp_r2(sumstats,
           beta="BETA",
           af="EAF",
           n = "N", 
           mode="q",
           se="SE",
           vary=1,
           ncase=None, 
           ncontrol=None, 
           prevalence=None, 
           k=1,
           log=Log(),
           adjuested=False,
           verbose=True):
    # Pierce, B. L., Ahsan, H., & VanderWeele, T. J. (2011). Power and instrument strength requirements for Mendelian randomization studies using multiple genetic variants. International journal of epidemiology, 40(3), 740-752.
    log.write("Start to calculate per-SNP heritibility...", verbose=verbose)
    if type(k) is int or type(k) is float:
       pass 
    elif k =="all":
        k = len(sumstats)

    if mode=="q":
        if beta in sumstats.columns and af in sumstats.columns:
            #Shim, H., Chasman, D. I., Smith, J. D., Mora, S., Ridker, P. M., Nickerson, D. A., ... & Stephens, M. (2015). A multivariate genome-wide association analysis of 10 LDL subfractions, and their response to statin treatment, in 1868 Caucasians. PloS one, 10(4), e0120758.
            # y = beta * x + e
            # Var(y) = beta**2*Var(x) + Var(e) 
            # Var(x) = 2*MAF*(1-MAF)
            # Var(beta * X) = beta**2 * Var(x) 
            # Var(e) = betase**2 * 2 * N * MAF * (1-MAF)
            # r2 = Var(beta * X) / Var(y)

            log.write(" -Calculating per-SNP rsq by 2 * (BETA**2) * AF * (1-AF) / Var(y)...", verbose=verbose)
            sumstats["_VAR(BETAX)"] = 2*(sumstats[beta]**2)*sumstats[af]*(1-sumstats[af])
            
            if type(vary) is int or type(vary) is float:
                log.write(" -Var(y) is provided: {}...".format(vary), verbose=verbose)
                sumstats["SNPR2"] = sumstats["_VAR(BETAX)"] / vary
            elif vary=="se":
                log.write(" -Var(y) is estimated from VAR(BETA * X), N, MAF, SE: {}...".format(vary), verbose=verbose)
                sumstats["_SIGMA2"] = sumstats[se]**2 * 2*(sumstats[n])*sumstats[af]*(1-sumstats[af])
                sumstats["SNPR2"] = sumstats["_VAR(BETAX)"] / (sumstats["_SIGMA2"] + sumstats["_VAR(BETAX)"]) 
        else:
            log.warning("Not enough information for calculation.")
    
    if mode=="b":
        if ncase not in sumstats.columns:
            sumstats["_NCASE"] = ncase     
        if ncontrol not in sumstats.columns:
            sumstats["_NCONTROL"] = ncontrol    
        if prevalence not in sumstats.columns:
            sumstats["_PREVALENCE"] = prevalence                    
        
        # equation 10 in Lee, S. H., Goddard, M. E., Wray, N. R., & Visscher, P. M. (2012). A better coefficient of determination for genetic profile analysis. Genetic epidemiology, 36(3), 214-224.
        # code : TwoSampleMR https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/add_rsq.r

        ve = (np.pi**2)/3 #3.29
        sumstats["_POPAF"] = sumstats.apply(
            lambda x: get_population_allele_frequency(af = x[af], prop=x["_NCASE"] / (x["_NCASE"]+x["_NCONTROL"]), odds_ratio=np.exp(x[beta]), prevalence=x["_PREVALENCE"]),axis=1)
        sumstats["_VG"] = (sumstats[beta]**2)* sumstats["_POPAF"]*(1- sumstats["_POPAF"])
        sumstats["SNPR2"] = sumstats["_VG"] / (sumstats["_VG"] + ve)
        
    if adjuested==True:
        sumstats["ADJUESTED_SNPR2"] = 1 - (1-sumstats["SNPR2"]) * (sumstats[n]-1) / (sumstats[n]-1 -k)   
        snpr2 = "ADJUESTED_SNPR2"
    else:
        snpr2 = "SNPR2"
    if n in sumstats.columns:
        log.write(" -Calculating F-statistic: F = [(N-k-1)/k] * (r2/1-r2)... where k = {}".format(k), verbose=verbose)
        log.write(" -For r2, {} is used.".format(snpr2), verbose=verbose)
        sumstats["F"] = sumstats[snpr2]*(sumstats[n]-1 -k)/((1-sumstats[snpr2]) * k)
        
    log.write("Finished calculating per-SNP heritability!", verbose=verbose)
    return sumstats
#
def get_population_allele_frequency(af, prop, odds_ratio, prevalence,eps=1e-15):
    #OR = (a × d)/(b × c)

    #https://stats.stackexchange.com/questions/241384/how-to-derive-2x2-cell-counts-from-contingency-table-margins-and-the-odds-ratio
    
    # ax2 + bx + c = 0 -> sovle x

    a = odds_ratio - 1 
    b = (af+prop)*(1-odds_ratio)-1
    c = odds_ratio*af*prop

    d = (b**2) - (4*a*c)
    sol1 = (-b-np.sqrt(d))/(2*a)
    sol2 = (-b+np.sqrt(d))/(2*a) 
    
    for i in [sol1, sol2]:
        #             case    control      Total
        #  allele1     a        b            af
        #  allele2     c        d            1- af
        #   total     prop      1-prop       1

        a = sol1      # AF case allele1
        c = prop-a    # AF case allele2
        b = af-a      # AF control allele1
        d = 1+a-af-prop # AF control allele2
        
        if min([a,b,c,d])>=0:
            af_case = a / (a+c)
            af_control = b / (b+d)
            pop_af = af_case*prevalence + (1-prevalence )*af_control
            return pop_af
    else:
        return np.nan