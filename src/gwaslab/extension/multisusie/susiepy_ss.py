import numpy as np
from tqdm import tqdm
import scipy.linalg as la
import sys
import numba
import copy
import random
import functools
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize


class S: 
    def __init__(
            self, pop_sizes, L, XTX_list, scaled_prior_variance, 
            residual_variance, varY, prior_weights, float_type = np.float32):
        
        #code from init_setup
        num_pop = len(pop_sizes)
        p = XTX_list[0].shape[0]
        self.alpha = np.zeros((L, p), dtype=float_type) + 1.0/p
        self.mu = np.zeros((num_pop, L, p), dtype=float_type)
        self.mu2 = np.zeros((num_pop, num_pop, L, p), dtype=float_type)
        self.Xr_list = [np.zeros(XTX.shape[0], dtype=float_type) for XTX in XTX_list]
        self.sigma2 = residual_variance.astype(float_type)
        self.pi = prior_weights.astype(float_type)
        self.n = np.array(pop_sizes)
        self.L = L
        
        #code from init_finalize
        self.V = scaled_prior_variance * varY + np.zeros((num_pop, L), dtype=float_type)
        self.V = self.V.astype(float_type)
        assert np.all(self.V >= 0)
        self.ER2 = np.zeros(num_pop, dtype=float_type) + np.nan
        self.KL = np.zeros(L, dtype=float_type) + np.nan
        self.lbf = np.zeros(L, dtype=float_type) + np.nan
        
        self.converged = False

class SER_RESULTS:
    def __init__(self, alpha, mu, mu2, lbf, lbf_model, V):
        self.alpha = alpha
        self.mu = mu
        self.mu2 = mu2
        self.lbf = lbf
        self.V = V
        self.lbf_model = lbf_model

def multisusie_rss(
    R_list,
    population_sizes,
    b_list = None,
    s_list = None,
    z_list = None,
    varY_list = None,
    rho = 0.75,
    L = 10,
    scaled_prior_variance = 0.2,
    prior_weights = None,
    standardize = False, 
    pop_spec_standardization = True,
    estimate_residual_variance = True,
    estimate_prior_variance = True,
    estimate_prior_method = 'early_EM',
    pop_spec_effect_priors = True,
    iter_before_zeroing_effects = 5,
    prior_tol = 1e-9,
    max_iter = 100,
    tol = 1e-3,
    verbose = False,
    coverage = 0.95,
    min_abs_corr = 0,
    float_type = np.float32,
    low_memory_mode = False,
    recover_R = False,
    single_population_mac_thresh = 20,
    mac_list = None,
    multi_population_maf_thresh = 0,
    maf_list = None,
    variant_ids = None
    ):
    """ Top-level function for running MultiSuSiE with summary statistics

    This function takes takes standard GWAS summary statistics, converts
    them to sufficient statistics, and runs MultiSuSiE on them.

    Parameters
    ----------
    multisusie_rss accepts two combinations of input parameters:
    1. b_list, s_list, R_list, varY_list, rho, population_sizes
        This is the preferred input format and is used in the MultiSuSiE paper.
    2. z_list, R_list, rho, population_sizes
        Here, we assume that both the genotypes and phenotypes have been 
        standardized to have variance 1 within each population.  It's extremely 
        important to censor low MAF variants in the population where they are rare. 
        This is done by setting values in z_list for these variants to 0, and 
        columns and rows corresponding to these variants to 0. multisusie_rss
        will not do this for you (but will with the first input combination) 
        because it doesn't have access to the information that would allow us to 
        approximate the minor allele frequency or minor allele count. 
        It's also kind of weird to standardize the genotypes within each 
        population beacuse variants that are rare in one population and common 
        in other other will be assigned huge effect sizes.

    R_list: length K list of PxP numpy arrays representing the LD correlation 
        matrices for each population. Each array should correspond
        to the same set of P variants, in the same order. Variants that should
        not be included in a given population due to low MAF can be assigned a
        value of np.nan.
    population_sizes: list of integers representing the number of samples in
        the GWAS for each population
    b_list: length K list of length P numpy arrays containing effect sizes,
        one for each population. Each array should correspond
        to the same set of P variants, in the same order. Variants that should
        not be included in a given population due to low MAF can be assigned a
        value of np.nan. Provide exactly one of (b_list and s_list) and z_list.
    s_list: length K list of length P numpy arrays containing effect size 
        standard errors, one for each population. Each array should correspond
        to the same set of P variants, in the same order. Variants that should
        not be included in a given population due to low MAF can be assigned a
        value of np.nan. Provide exactly one of (b_list and s_list) and z_list.
    z_list: length K list of length P numpy arrays containing Z-scores, one for
        each population. Each array should correspond to the same set of P
        variants, in the same order. Variants that should not be included in a
        given population due to low MAF can be assigned a value of np.nan.
        Provide exactly one of (b_list and s_list) and z_list.
    varY_list: length K list representing the sample variance of the outcome
        in each population
    rho: PxP numpy array representing the effect size correlation matrix. In 
        the manuscript, we show that this parameter has little impact on the
        estimated PIPs in practice.
    L: integer representing the maximum number of causal variants
    scaled_prior_variance: float representing the effect size prior variance,
        scaled by the residual variance. It's fine to set this to a number 
        larger than what you expect the squared effect size to be (like the 
        default value of 0.2) as long as estimate_prior_variance is set to True
        and estimate_prior_method is not set to None.
    prior_weights: numpy P-array of floats representing the prior probability
        of causality for each variant. Give None to use a uniform prior
    standardize: boolean, whether to adjust summmary statistics to be as if 
        genotypes were standardized to have mean 0 and variance 1
    pop_spec_standardization: boolean, if standardize is True, whether to 
        adjust summary statistics to be as if genotypes were standardized
        separately for each population, or pooled and then standardized
    estimate_residual_variance: boolean, whether to estimate the residual variance, $\\sigma^2_k$ in the manuscript
    estimate_prior_variance: boolean, whether to estimate the prior variance,
        $A^{(l)}$ in the manuscript
    estimate_prior_method: string, method to estimate the prior variance. Recommended
        values are 'early_EM' or None
    pop_spec_effect_priors: boolean, whether to estimate separate prior 
        variance parameters for each population
    iter_before_zeroing_effects: integer, number of iterations to run before
        zeroing out component-population pairs (or components if 
        pop_spec_effect_priors is False) that have a lower likelihood than a 
        null model
    prior_tol: float which places a filter on the minimum prior variance
        for a component to be included when estimating PIPs
    max_iter: integer, maximum number of iterations to run
    tol: float, after iter_before_zeroing_effects iterations, results
        are returned if the ELBO increases by less than tol in an ieration
    verbose: boolean which indicates if an progress bar should be displayed
    coverage: float representing the minimum coverage of credible sets
    min_abs_corr: float representing the minimum absolute correlation between
        any pair of variants in a credible set (purity). For each pair of variants,
        the max is taken across ancestries. In the case where min_abs_corr = 0,
        low_memory_mode = True, and recover_R = False, the purity of credible
        sets will not be calculated. 
    float_type: numpy float type used. Set to np.float32 to minimize memory
        consumption
    low_memory_mode: boolean, if True, the input R_list will be modified in place.
        Decreases memory consumption by a factor of two. BUT, THE INPUT R_list
        WILL BE OVERWRITTEN. If you need R_list, set low_memory_mode to False
    recover_R: boolean, if True, the R matrices will be recovered from XTX, 
        BUT variants with MAC/MAF estimate less than mac_filter/maf_filter will
        be censored. 
    single_population_mac_thresh: float, variants with minor allele count less
        than this value in a single population will be censored in that population.
        The variant will not be censored in the other populations. 
        If mac_list is not None, mac_list is used for the filtering.
        If mac_list is None, but maf_list is not, mac is caclulated using
        maf_list and population_sizes. Otherwise, if mac_list is None, 
        maf_list is None, and b_list and s_list are provided as input, then
        minor allele counts will be estimated under an assumption of
        Hardy-Weinberg equilibrium. 
    mac_list: length K list of floats representing the minor allele count for
        each variant in each population.
    multi_population_maf_thresh: float, variants with minor allele frequency 
        less than this value in ALL poopulations will be censored. Note that 
        it's probably faster computationally to do this before calling this 
        function. If maf_list is not None, maf_list is used for the filtering.
        If maf_list is None, but mac_list is not, maf is caclulated using
        mac_list and population_sizes. Otherwise, if maf_list is None, 
        mac_list is None, and b_list and s_list are provided as input, then
        minor allele frequencies will be estimated under an assumption of
        Hardy-Weinberg equilibrium. 
    maf_list: length K list of floats representing the minor allele frequency for
        each variant in each population. 
    variant_ids: length P list of strings representing the variant IDs. If 
        provided, sets will contain a fifth entry, containing the variant ids
        of the variants contained in each set. 
    
    Returns
    -------
    an object containing results with the following attributes:
        pip: length-P numpy array of posterior inclusion probabilities
        coef: K x P numpy array of posterior effect size estimates for 
            each variant in each population, aggregated across
            all L single effect regressions. These effects are not 
            conditional on each variant being a causal variant and thus
            account for uncertainty in the causal variant assignment.
        coef_sd: K x P numpy array of posterior effect size standard deviations
            for each variant in each population, aggregated across
            all L single effect regressions.
        sets: a list with four entries. The first contains the indices of the 
            variants contained in each set. The second contains the purity
            of each set (see min_abs_corr for calculation details). The 
            third contains the coverage of each set. The fourth contains
            whether the set has passed filtering. If variant_ids is 
            provided, a fifth entry contains the variant ids for the 
            variants contained in each set. 
        alpha: L x P numpy array of single-effect regression posterior 
            inclusion probabilities
        mu: K x L x P numpy array of single-effect regression effect size
            posterior means, conditional on each variant being the causal variant
        mu2: K x K x L x P numpy array of single-effect regression effect size
            posterior seconds moments, conditional on each variant being the 
            causal variant
        sigma2: length-K numpy array of residual variance estimates
        pi: length-P numpy array of prior inclusion probabilities
        n: length-K numpy array of sample sizes
        L: integer representing the maximum number of causal variants
        V: K x L numpy array of effect size prior variance estimates
        ER2: length-K numpy array of expected squared residuals
        KL: L x 1 numpy array of Kullback-Leibler divergences for each single
            effect regression
        lbf: L x 1 numpy array of log Bayes factors for each single effect 
            regression
        converged: boolean indicating whether the algorithm converged
        niter: integer representing the number of IBSS iterations required
            to reach convergence

            

    TODO
    ----
        - Add command line interface
        - Add more tests, probably based on tests used for SuSiER
        - Make defaults match between sufficient statistic and summary statistic
          functions
        - Implement missing functionality for individual-level multisusie
    """

    
    # provide exactly one of (b_list and s_list) and z_list
    if (z_list is not None) & ((b_list is not None) or (s_list is not None)):
        raise ValueError('provide either (b_list and s_list) or z_list, but not both. ' + \
                         'See the parameters section of help(multisusie_rss) for more information')
    if (z_list is None) & (b_list is None) & (s_list is None):
        raise ValueError('provide either (b_list and s_list) or z_list, but not both')
    if (b_list is not None) & ((s_list is None) | (varY_list is None)):
        raise ValueError('if b_list is provided, s_list and varY_list must also be provided')
    if population_sizes is None:
        raise ValueError('population_sizes must be provided')


    # check input list lengths
    K = len(R_list)
    if z_list is None:
        assert len(s_list) == K
        assert len(b_list) == K
    else:
        assert len(z_list) == K
    if varY_list is not None:
        assert len(varY_list) == K
    assert len(population_sizes) == K
    assert np.isscalar(rho) | ((rho.shape[0] == K) & (rho.shape[1] == K))

    # make copies of R if low_memory_mode=False to avoid modifying the input data 
    if low_memory_mode:
        R_list_copy = R_list
        if z_list is None:
            b_list_copy = b_list
            s_list_copy = s_list
        else:
            z_list_copy = z_list
        print('low memory mode is on. THE INPUT R MATRICES HAVE BEEN ' +
            'TRANSFORMED INTO XTX AND CENSORED BASED ON MISSINGNESS. ' +
            'THE INPUT R MATRICES HAVE BEEN CHANGED. The datatypes of ' +
            'b_list, s_list, and R_list have been mutated in place ' +
            'to match float_type')
    else:
        R_list_copy = [np.copy(R) for R in R_list]
        if z_list is None:
            b_list_copy = [np.copy(b) for b in b_list]
            s_list_copy = [np.copy(s) for s in s_list]
        else:
            z_list_copy = [np.copy(z) for z in z_list]

    # convert everything to the requested float type if necessary
    for i in range(K):
        if z_list is None:
            if b_list_copy[i].dtype != float_type:
                b_list_copy[i] = b_list_copy[i].astype(float_type, copy = False)
            if s_list_copy[i].dtype != float_type:
                s_list_copy[i] = s_list_copy[i].astype(float_type, copy = False)
        else:
            if z_list_copy[i].dtype != float_type:
                z_list_copy[i] = z_list_copy[i].astype(float_type, copy = False)
        if R_list_copy[i].dtype != float_type:
            R_list_copy[i] = R_list_copy[i].astype(float_type, copy = False)
        if (not np.isscalar(rho)) and (rho.dtype != float_type):
            rho = rho.astype(float_type, copy = True)

    population_sizes = np.array(population_sizes, dtype = float_type)
    if varY_list is None:
        varY_list = np.ones(K, dtype = float_type)
    else:
        varY_list = np.array(varY_list, dtype = float_type)
    YTY_list = []
    for i in range(K):
        YTY = (varY_list[i] * (population_sizes[i] - 1)).astype(float_type)
        YTY_list.append(YTY)


    # convert GWAS summary statistics to MultiSuSiE sufficient statistics
    XTX_list = []
    XTY_list = []
    if z_list is None:
        for i in range(K):
            # this function mutates R_list_copy[i]
            XTY = recover_XTX_and_XTY(
                b = b_list_copy[i],
                s = s_list_copy[i],
                R = R_list_copy[i],
                YTY = YTY_list[i],
                n = population_sizes[i]
            )
            XTY_list.append(XTY)
        XTX_list = R_list_copy
    else:
        for i in range(K):
            # this function mutates R_list_copy[i]
            XTX, XTY = recover_XTX_and_XTY_from_Z(
                z = z_list_copy[i],
                R = R_list_copy[i],
                n = population_sizes[i],
                float_type = float_type
            )
            XTX_list.append(XTX)
            XTY_list.append(XTY)

    # do minor allele count filtering for each population
    if single_population_mac_thresh > 0:
        for i in range(K):
            n = population_sizes[i]
            if mac_list is not None:
                mac = np.minimum(mac_list[i], 2 * n - mac_list[i])
            elif maf_list is not None:
                mac = np.minimum(2 * n * maf_list[i], 2 * n - 2 * n * maf_list[i])
            elif b_list is not None:
                var_x = np.minimum(np.diag(XTX_list[i]) / (n - 1), .5) # this can be > .5 due to HWE violations?
                maf = 1 / 2 - np.sqrt(1 - 2 * var_x) / 2
                mac = np.minimum(2 * n * maf, 2 * n - 2 * n * maf)
            else:
                print(f'Skipping MAC filtering because mac_list is not provided,' +
                      'maf_list is not provided, and z-scores are being used')
                continue
            low_mac_mask = mac < single_population_mac_thresh
            XTX_list[i][low_mac_mask, :] = 0
            XTX_list[i][:, low_mac_mask] = 0
            XTY_list[i][low_mac_mask] = 0
            print(f'Censored {np.sum(low_mac_mask)} variants in population {i}' +
                  ' due to low population-specific MAC')

    # do minor allele frequency filtering across populations. this filter
    # only censors a variant if it has low MAF in all populations
    if multi_population_maf_thresh > 0:
        low_maf_across_populations_mask = np.ones(XTX_list[0].shape[0], dtype=bool)
        for i in range(K):
            n = population_sizes[i]
            if maf_list is not None:
                maf = np.minimum(maf_list[i], 1 - maf_list[i])
            elif mac_list is not None:
                maf = np.minimum(mac_list[i] / (2 * n), 1 - mac_list[i] / (2 * n))
            elif b_list is not None:
                var_x = np.minimum(np.diag(XTX_list[i]) / (n - 1), .5)
                maf = 1 / 2 - np.sqrt(1 - 2 * var_x) / 2
            else:
                print(f'Skipping MAF filtering because mac_list is not provided,' +
                      'maf_list is not provided, and z-scores are being used')
            low_maf_across_populations_mask = low_maf_across_populations_mask & (maf < multi_population_maf_thresh)
        print(f'Censored {np.sum(low_maf_across_populations_mask)} variants in' +
            ' all populations due to low MAF across all populations')
        for i in range(K):
            XTX_list[i][low_maf_across_populations_mask, :] = 0
            XTX_list[i][:, low_maf_across_populations_mask] = 0
            XTY_list[i][low_maf_across_populations_mask] = 0



    if np.isscalar(rho):
        rho = np.ones((K, K), dtype=float_type) * rho
        rho[np.diag_indices(K)] = 1
    
    s = susie_multi_ss(
        XTX_list = XTX_list,
        XTY_list = XTY_list, 
        YTY_list = YTY_list, 
        rho = rho, 
        population_sizes = population_sizes,
        L = L, 
        scaled_prior_variance = scaled_prior_variance,
        prior_weights = prior_weights,
        standardize = standardize,
        pop_spec_standardization = pop_spec_standardization,
        estimate_residual_variance = estimate_residual_variance,
        estimate_prior_variance = estimate_prior_variance,
        estimate_prior_method = estimate_prior_method,
        prior_tol = prior_tol,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose,
        iter_before_zeroing_effects = iter_before_zeroing_effects,
        pop_spec_effect_priors = pop_spec_effect_priors,
        R_list = R_list,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        float_type = float_type,
        low_memory_mode = low_memory_mode,
        recover_R = recover_R,
        variant_ids = variant_ids
    )

    return s

susie_multi_rss = multisusie_rss

def recover_XTX_and_XTY(b, s, R, YTY, n):
    """ Recover XTX and XTY from GWAS summary statistics. INPUT R IS MUTATED 

    THIS FUNCTION MUTATES INPUT R to reduce memory consumpation. BE CAREFUL!! 
    see page 4 of the supplement of Zou et al. 2022 PLoS Genetics for a derivation

    Parameters
    ----------
    b: length-P numpy array of floats representing effect sizes
    s: length-P numpy array of floats representing standard errors
    R: PxP numpy array of floats representing the LD correlation matrix
    YTY: float representing the sample variance of the outcome,
    n: integer representing the number of samples in the GWAS

    Returns
    -------
    R: PxP numpy array representing the LD correlation matrix, with low MAF/MAC
        variants censored. Note that this is the same object as the input R
    XTY: length-P numpy array representing the X^T Y vector
    n_censored: integer representing the number of variants censored
    """

    #  see page 4 of the supplement of Zou et al. 2022 PLoS Genetics for a derivation
    sigma2 = YTY / ((b / s) ** 2 + n - 2)
    XTY = np.nan_to_num(sigma2 * b / (s ** 2))
    dR = np.nan_to_num(sigma2 / (s ** 2))
    R *= np.expand_dims(np.sqrt(dR), axis = 1)
    R *=  np.sqrt(dR)
    np.nan_to_num(R, copy = False)

    return XTY

def recover_XTX_and_XTY_from_Z(z, R, n, float_type = np.float32):
    """ Recover XTX and XTY from z scores and LD correlation matrix

    THIS FUNCTION MUTATES INPUT R to reduce memory consumpation. BE CAREFUL!! 
    This is equivalent to using standardized genotype and phenotype

    Parameters
    ----------
    z: length-P numpy array of floats representing GWAS z scores
    R: PxP numpy array of floats representing the LD correlation matrix
    n: integer representing the number of samples in the GWAS

    Returns
    -------
    XTX: PxP numpy array representing the X^T X matrix
    XTY: length-P numpy array representing the X^T Y vector
    """
    adj = (n - 1) / (z ** 2 + n - 2)
    z = np.sqrt(adj) * z
    R *= (n - 1)
    XTY = (np.sqrt(n - 1) * z).astype(float_type)
    np.nan_to_num(XTY, copy = False)
    np.nan_to_num(R, copy = False)

    return R, XTY

def susie_multi_ss(
    XTX_list, XTY_list, YTY_list,
    rho,
    population_sizes,
    L = 10,
    scaled_prior_variance=0.2,
    prior_weights=None,
    standardize = True,
    pop_spec_standardization = False,
    estimate_residual_variance=True,
    estimate_prior_variance=True,
    estimate_prior_method='early_EM',
    pop_spec_effect_priors = True,
    iter_before_zeroing_effects = 5,
    prior_tol=1e-9,
    max_iter=100,
    tol=1e-3,
    verbose=False,
    R_list = None,
    coverage = .95,
    min_abs_corr = .5,
    float_type = np.float32,
    low_memory_mode = False,
    recover_R = False,
    variant_ids = None
    ):
    """ Run MultiSuSiE on sufficient statistics

    This function runs MultiSuSiE on sufficient statistics. It is designed to 
    be run from the top-level function multisusie_rss, but can be run directly.
    It will not censor based on MAF/MAC, unlike multisusie_rss.

    Parameters
    ----------
    XTX_list: length K list of PxP numpy arrays representing the X^T X matrices
    XTY_list: length K list of length P numpy arrays representing the X^T Y vectors
    YTY_list: length K list of floats representing the sample variance of the outcome
    rho: PxP numpy array representing the effect size correlation matrix
    population_sizes: list of integers representing the number of samples in
        the GWAS for each population
    L: integer representing the maximum number of causal variants
    scaled_prior_variance: float representing the effect size prior variance,
        scaled by the residual variance. It's fine to set this to a number 
        larger than what you expect the squared effect size to be (like the 
        default value of 0.2) as long as estimate_prior_variance is set to True
        and estimate_prior_method is not set to None.
    prior_weights: numpy P-array of floats representing the prior probability
        of causality for each variant. Give None to use a uniform prior
    standardize: boolean, whether to adjust summmary statistics to be as if 
        genotypes were standardized to have mean 0 and variance 1
    pop_spec_standardization: boolean, if standardize is True, whether to 
        adjust summary statistics to be as if genotypes were standardized
        separately for each population, or pooled and then standardized
    estimate_residual_variance: boolean, whether to estimate the residual variance,
        $\\sigma^2_k$ in the manuscript
    estimate_prior_variance: boolean, whether to estimate the prior variance,
        $A^{(l)}$ in the manuscript
    estimate_prior_method: string, method to estimate the prior variance. Recommended
        values are 'early_EM' or None
    pop_spec_effect_priors: boolean, whether to estimate separate prior 
        variance parameters for each population
    iter_before_zeroing_effects: integer, number of iterations to run before
        zeroing out component-population pairs (or components if 
        pop_spec_effect_priors is False) that have a lower likelihood than a 
        null model
    prior_tol: float which places a filter on the minimum prior variance
        for a component to be included when estimating PIPs
    max_iter: integer, maximum number of iterations to run
    tol: float, after iter_before_zeroing_effects iterations, results
        are returned if the ELBO increases by less than tol in an ieration
    verbose: boolean which indicates if an progress bar should be displayed
    R_list: length K list of PxP numpy arrays representing the LD correlation.
        If set to None, and min_abs_corr > 0, the LD correlation matrices will
        be recovered from XTX_list.
    coverage: float representing the minimum coverage of credible sets
    min_abs_corr: float representing the minimum absolute correlation between
        any pair of variants in a credible set. For each pair of variants,
        the max is taken across ancestries. In the case where min_abs_corr = 0,
        low_memory_mode = True, and recover_R = False, the purity of credible
        sets will not be calculated. 
    float_type: numpy float type used. Set to np.float32 to minimize memory
        consumption
    low_memory_mode: boolean, if True, the input R_list will be modified in place.
        Decreases memory consumption by a factor of two. BUT, THE INPUT R_list
        WILL BE OVERWRITTEN. If you need R_list, set low_memory_mode to False
    recover_R: boolean, if True, the R matrices will be recovered from XTX, 
        BUT variants with MAC/MAF estimate less than mac_filter/maf_filter will
        be censored. 

    Returns
    -------
    an object containing results with the following attributes:
        alpha: L x P numpy array of single-effect regression posterior 
            inclusion probabilities
        mu: K x L x P numpy array of single-effect regression effect size
            posterior means, conditional on each variant being the causal variant
        mu2: K x K x L x P numpy array of single-effect regression effect size
            posterior seconds moments, conditional on each variant being the 
            causal variant
        sigma2: length-K numpy array of residual variance estimates
        pi: length-P numpy array of prior inclusion probabilities
        n: length-K numpy array of sample sizes
        L: integer representing the maximum number of causal variants
        V: K x L numpy array of effect size prior variance estimates
        ER2: length-K numpy array of expected squared residuals
        KL: L x 1 numpy array of Kullback-Leibler divergences for each single
            effect regression
        lbf: L x 1 numpy array of log Bayes factors for each single effect 
            regression
        converged: boolean indicating whether the algorithm converged
    """

    #check input dimensions
    assert len(XTX_list) == len(XTY_list)
    assert np.all([XTX.shape[1] == XTX_list[0].shape[1] for XTX in XTX_list])
    assert np.all([XTX.shape[0] == XTY.shape[0] for (XTX,XTY) in zip(XTX_list, XTY_list)])
    if prior_weights is not None:
        prior_weights = prior_weights.astype(float_type)

    assert not np.any([np.any(np.isnan(XTX)) for XTX in XTX_list])
        
    #compute w_pop (the relative size of each population)
    population_sizes = np.array(population_sizes)
    w_pop = (population_sizes/ population_sizes.sum()).astype(float_type)
    
    #compute rho properties
    rho = rho.astype(float_type)
    logdet_rho_sign, logdet_rho = np.linalg.slogdet(rho)
    assert logdet_rho_sign>0

    X_l2_arr = np.array([np.diag(XTX) for XTX in XTX_list], dtype = float_type)

    # calculate scaling factors
    if standardize:
        csd = np.sqrt(X_l2_arr/(np.expand_dims(population_sizes, axis = 1) - 1))
        if not pop_spec_standardization:
            csd = csd.T.dot(w_pop)
            csd = csd * np.ones((len(XTX_list), csd.shape[0]), dtype = float_type)
        for pop_i in range(len(XTX_list)):
            XTX_list[pop_i] *= (1 / csd[pop_i, :])
            XTX_list[pop_i] *= (1 / np.expand_dims(csd[pop_i, :], 1))
            XTY_list[pop_i] = XTY_list[pop_i] / csd[pop_i, :]
        X_l2_arr = np.array([np.diag(XTX) for XTX in XTX_list], dtype = float_type)
    else:
        csd = np.ones((len(XTX_list), XTX_list[0].shape[1]), dtype = float_type)
    varY = np.array([YTY/(n-1) for (YTY, n) in zip(YTY_list, population_sizes)], dtype = float_type)
    varY_pooled = (np.sum(YTY_list) / (np.sum(population_sizes) - 1)).astype(float_type)
    residual_variance = varY

    #init setup
    p = XTX_list[0].shape[1]

    if np.isscalar(scaled_prior_variance) & standardize:
        assert 0 < scaled_prior_variance <= 1
    if prior_weights is None:
        prior_weights = np.zeros(p, dtype=float_type) + 1.0/p
    else:
        prior_weights = (prior_weights / np.sum(prior_weights)).astype(float_type)
    assert prior_weights.shape[0] == p
    if p<L: L=p

    # initialize s, which tracks the current model parameter estimates
    s = S(
        population_sizes, L, XTX_list, scaled_prior_variance, residual_variance,
        varY_pooled, prior_weights, float_type = float_type
    )

    elbo = np.zeros(max_iter+1) + np.nan
    elbo[0] = -np.inf
    check_null_threshold = 0.0
    
    ### start iterations ###
    tqdm_iter = tqdm(list(range(max_iter)), disable=not verbose, file=sys.stdout)
    for i in tqdm_iter:
        tqdm_iter.set_description('iteration %d/%d'%(i+1, max_iter))

        # set the prior variance estimation method for this iteration
        if estimate_prior_method == 'early_EM':
            if i == 0: 
                current_estimate_prior_method = None
            else:
                current_estimate_prior_method = 'early_EM'
        elif (estimate_prior_method is not None) and (estimate_prior_method.split('_')[-1] == 'init'):
            if i == 0:
                current_estimate_prior_method = estimate_prior_method.split('_')[0]
            else:
                current_estimate_prior_method = 'EM'
        else:
            current_estimate_prior_method = estimate_prior_method

        if i < iter_before_zeroing_effects:
            current_check_null_threshold = -np.inf
        else:
            current_check_null_threshold = check_null_threshold

        # update all L single effect regressions. s is modified in place.
        update_each_effect(
            XTX_list, XTY_list, s, X_l2_arr, w_pop,
            rho,
            estimate_prior_variance, current_estimate_prior_method,
            check_null_threshold = current_check_null_threshold,
            pop_spec_effect_priors = pop_spec_effect_priors,
            float_type = float_type
        )

        # update the ER2 parameter of s. Used by get_objective and 
        # estimate_resdiaul_variance_func
        s.ER2 = update_ER2(XTX_list, XTY_list, YTY_list, s, X_l2_arr)

        # update the ELBO and check for convergence
        elbo[i+1] = get_objective(XTX_list, XTY_list, s, YTY_list, X_l2_arr)
        tqdm_iter.set_postfix(objective="{:.6f}".format(elbo[i]))
        if ((elbo[i+1] - elbo[i]) < tol) and (i >= (iter_before_zeroing_effects + 1)):
            s.converged = True
            tqdm_iter.close()
            break

        if estimate_residual_variance:
            s.sigma2 = estimate_residual_variance_func(
                XTX_list, XTY_list, YTY_list, s, X_l2_arr, population_sizes, 
                float_type = float_type
            )
            if np.any(s.sigma2 < 0):
                print('minimum resdiual variance less than 0. Is there mismatch between the correlation matrix and association statistics?')

    elbo = elbo[1:i+2] # Remove first (infinite) entry, and trailing NAs.
    s.elbo = elbo
    s.niter = i+1

    if not s.converged:
        print('IBSS algorithm did not converge in %d iterations'%(max_iter))
        
    s.fitted = s.Xr_list

    s.pip = susie_get_pip(s, prior_tol=prior_tol)
    s.coef = np.array([np.squeeze(np.sum(s.mu[k] * s.alpha, axis=0) / csd[k, :]) for k in range(len(XTX_list))])
    s.coef_sd = np.array([(np.squeeze(np.sqrt(np.sum(s.alpha * s.mu2[k, k] - (s.alpha*s.mu[k])**2, axis=0)) / csd[k, :])) for k in range(len(XTX_list))])

    if (low_memory_mode and (min_abs_corr > 0)) or (low_memory_mode and recover_R):
        for i in range(len(XTX_list)):
            recover_R_from_XTX(XTX_list[i], X_l2_arr[i])
        R_list = XTX_list

    s.sets = susie_get_cs(
        s = s, R_list = R_list, coverage = coverage, min_abs_corr = min_abs_corr, 
        dedup = True, n_purity = np.inf, 
        calculate_purity = (not low_memory_mode) or (min_abs_corr > 0) or recover_R
    )

    s.variant_ids = variant_ids
    if variant_ids is not None:
        cs_variant_ids = []
        for i in range(len(s.sets[0])):
            cs_variant_ids.append([variant_ids[idx] for idx in s.sets[0][i]])
        s.sets.append(cs_variant_ids)

    return s

def update_each_effect(
    XTX_list, XTY_list, s, X_l2_arr, w_pop, rho,
    estimate_prior_variance=False, estimate_prior_method='optim',
    check_null_threshold=0.0, pop_spec_effect_priors = True, float_type = np.float32):
    """ Update each single effect regression

    This calculates updated single effect regression parameter estimates

    Parameters:
    -----------
    XTX_list: length K list of PxP numpy arrays representing the X^T X matrices
    XTY_list: length K list of length P numpy arrays representing the X^T Y vectors
    s: S object representing the model parameter estimates prior to the update
    X_l2_arr: length K numpy array of floats representing the diagonal of X^T X
    w_pop: length K numpy array of floats representing the relative size of each population
    rho: PxP numpy array representing the effect size correlation matrix
    estimate_prior_variance: boolean, whether to estimate effect size prior variance
    estimate_prior_method: string, method to estimate the prior variance. Recommended
        values are 'early_EM' or None
    check_null_threshold: float representing the difference in loglikelihood that
        a component with a non-zero effect size prior variance must beat the null
        model by to be retained.
    pop_spec_effect_priors: boolean, whether to estimate separate prior effect
        size variances for each population
    float_type: numpy float type used. Set to np.float32 to minimize memory

    Returns:
    -------
    Nothing. S is modified in place.
    """

        
    if not estimate_prior_variance:
        estimate_prior_method = None
    L = s.alpha.shape[0]
    num_pop = len(XTX_list)

    for l in range(L):
        R_list = []

        # add the estimate of the current effect back to the residualized
        # XTXY vector
        for k in range(len(XTX_list)):
            s.Xr_list[k] -= XTX_list[k].dot(s.alpha[l] * s.mu[k,l])
            R_list.append(XTY_list[k] - s.Xr_list[k])
        
        # get the currrent single effect model parameter estimates after 
        # residualizing all other effects
        res = single_effect_regression(
            R_list, XTX_list, s.V[:,l], X_l2_arr, w_pop, rho,
            residual_variance=s.sigma2, prior_weights=s.pi,
            optimize_V=estimate_prior_method, 
            check_null_threshold=check_null_threshold, 
            pop_spec_effect_priors = pop_spec_effect_priors,
            alpha = s.alpha[l,:], mu2 = s.mu2[:,:,l,:],
            float_type = float_type
        )
              
        # Update the variational estimate of the posterior mean.
        s.mu[:,l,:] = res.mu
        s.alpha[l,:] = res.alpha
        s.mu2[:,:,l,:] = res.mu2
        s.V[:,l] = res.V
        s.lbf[l] = res.lbf_model
        s.KL[l] = -res.lbf_model + SER_posterior_e_loglik(X_l2_arr, R_list, s.sigma2, res.mu * res.alpha, res.mu2[range(num_pop), range(num_pop)] * res.alpha)
        for k in range(len(XTX_list)):
            s.Xr_list[k] += XTX_list[k].dot(s.alpha[l] * s.mu[k,l])

def single_effect_regression(
    XTY_list, XTX_list, V, X_l2_arr, w_pop, rho,
    residual_variance, prior_weights=None, optimize_V=None, 
    check_null_threshold=0,  pop_spec_effect_priors = True,
    alpha = None, mu2 = None, float_type = np.float32):
    """ Fit a multi-population single effect regression model

    Parameters:
    -----------
    XTY_list: length K list of length P numpy arrays representing the X^T Y vectors
    XTX_list: length K list of PxP numpy arrays representing the X^T X matrices
    V: length K numpy array of floats representing the effect size prior variance
    X_l2_arr: length K numpy array of floats representing the diagonal of X^T X
    w_pop: length K numpy array of floats representing the relative size of each population
    rho: PxP numpy array representing the effect size correlation matrix
    residual_variance: length K numpy array of floats representing the residual 
        variance in each population
    prior_weights: length P numpy array of floats representing the prior probability
        of causality for each variant
    optimize_V: string representing the method to use to optimize the effect size
    check_null_threshold: float representing the difference in loglikelihood that
        a component with a non-zero effect size prior variance must beat the null
        model by to be retained.
    pop_spec_effect_priors: boolean, whether to estimate separate prior effect
        size variances for each population
    alpha: length P numpy array of floats representing the current posterior
        inclusion probability estimates
    mu2: K x K x P numpy array of floats representing the current posterior
        second moment estimates
    float_type: numpy float type used. Set to np.float32 to minimize memory

    Returns:
    -------
    an object containing results with the following attributes:
        alpha: length P numpy array of single effect regression posterior
            inclusion probabilities
        mu: K x P numpy array of single effect regression effect sizes conditional
            on each variant being the causal variant
        mu2: K x K x P numpy array of single effect regression effect size second
            moments conditional on each variant being the causal variant
        lbf: K x P numpy array representing the log Bayes factor for each 
            variant being the causal variant
        lbf_model: float representing the log Bayes factor for the model, 
            aggregating over variants
        V: length K numpy array of floats representing the effect size prior
    """
    
    compute_lbf_params = (XTY_list, XTX_list, X_l2_arr, rho,
        residual_variance, False, float_type)

    if optimize_V not in ['EM', 'EM_corrected', None]:
        V = optimize_prior_variance(
            optimize_V, prior_weights, rho.shape[0], 
            compute_lbf_params=compute_lbf_params, alpha=alpha, post_mean2=mu2, 
            w_pop=w_pop, check_null_threshold=check_null_threshold, 
            pop_spec_effect_priors = pop_spec_effect_priors, current_V = V,
            float_type = float_type
        )
        
    # compute the log Bayes factor for each variant being the causal variant 
    # and posterior mean estimates for each variant  under the assumption 
    # that it is the causal variant
    lbf, post_mean, post_mean2 = compute_lbf(
        V = V, XTY_list = XTY_list, XTX_list = XTX_list, X_l2_arr = X_l2_arr, 
        rho = rho, residual_variance = residual_variance, return_moments=True, 
        float_type = float_type
    )

    
    # compute posterior inclusion probabilities (not combined over single effect regressions)
    maxlbf = np.max(lbf)
    w = np.exp(lbf - maxlbf)
    w_weighted = w * prior_weights
    weighted_sum_w = w_weighted.sum()
    alpha = w_weighted / weighted_sum_w
    lbf_model = maxlbf + np.log(weighted_sum_w)
        
    if optimize_V in ['EM', 'EM_corrected']:
        V = optimize_prior_variance(
            optimize_V, prior_weights, rho.shape[0], 
            compute_lbf_params=compute_lbf_params, alpha=alpha, 
            post_mean2=post_mean2, w_pop=w_pop, 
            check_null_threshold=check_null_threshold, 
            pop_spec_effect_priors = pop_spec_effect_priors, 
            float_type = float_type
        )
        
    res = SER_RESULTS(alpha=alpha, mu=post_mean, mu2=post_mean2, lbf=lbf, lbf_model=lbf_model, V=V)


    return res

def loglik(V, prior_weights, compute_lbf_params):
    lbf = compute_lbf(V, *compute_lbf_params)
    maxlbf = np.max(lbf)
    w = np.exp(lbf - maxlbf)
    w_weighted = w * prior_weights
    weighted_sum_w = w_weighted.sum()
    loglik = maxlbf + np.log(weighted_sum_w)
    return loglik

def optimize_prior_variance(
    optimize_V, prior_weights, num_pops, compute_lbf_params=None, alpha=None, post_mean2=None, w_pop=None, check_null_threshold=0, 
    pop_spec_effect_priors = False, current_V = None, float_type = np.float32, loglik_function = loglik):

    if optimize_V == 'optim':
        if pop_spec_effect_priors:
            raise Exception('estimate_prior_method="optim" with ' +
                            'pop_spec_effect_priors=True has not been implemented')
        else:
            neg_loglik_logscale = lambda lV: -loglik_function(np.array([np.exp(lV) for i in range(num_pops)]), prior_weights, compute_lbf_params)
            opt_obj = minimize_scalar(neg_loglik_logscale, bounds=(-30,15))
            lV = opt_obj.x
            V = np.exp(lV)
    elif optimize_V in ['EM', 'early_EM']:
        V = np.array([np.sum(alpha * post_mean2[i, i]) for i in range(num_pops)], dtype=float_type)
        if not pop_spec_effect_priors:
            V = (w_pop.dot(V)).astype(float_type)
    elif optimize_V == 'grid':
        if pop_spec_effect_priors:
            raise Exception('estimate_prior_method="grid" with ' +
                            'pop_spec_effect_priors=True has not been implemented')
        else:
            V_arr = np.logspace(-7, -1, 13)
            llik_arr = np.array([loglik_function(V, prior_weights, compute_lbf_params) for V in V_arr])
            V = V_arr[np.argmax(llik_arr)]
    else:
        raise ValueError('unknown optimization method')

    if not pop_spec_effect_priors:
        V = V * np.ones(post_mean2.shape[0], dtype = float_type)
        # set V exactly 0 if that beats the numerical value by check_null_threshold in loglik.
        delta_loglik = loglik_function(0, prior_weights, compute_lbf_params) + check_null_threshold - loglik_function(V, prior_weights, compute_lbf_params)
        if np.isclose(delta_loglik, 0) or delta_loglik >= 0:
            V=0
    else:
        if check_null_threshold == (-1 * np.inf):
            return V
        elif np.all(np.isclose(V, 0)):
            return 0
        # Compare our current effect prior to null models 
        else:
            V_list = [V]
            # zero out each population, one at a time
            for i in range(num_pops):
                if not np.isclose(V[i],0):
                    V_copy = V.copy()
                    V_copy[i] = 0
                    V_list.append(V_copy)
            V_list.append(np.array([np.zeros(num_pops, dtype=float_type)]))
        llik_arr = np.array([loglik_function(np.array(V), prior_weights, compute_lbf_params) for V in V_list])
        llik_arr = llik_arr + np.array([0] + [check_null_threshold for i in range(len(V_list) - 1)])
        V = V_list[np.argmax(llik_arr)]
    if isinstance(V, np.ndarray):
        V[V < 0] = 0
    elif V < 0:
        V = 0

    return V

def compute_lbf(
    V, XTY_list, XTX_list, X_l2_arr, rho,
    residual_variance, return_moments=False,  
    float_type = np.float32):
    
    num_variables = XTX_list[0].shape[0]
    num_pops = len(XTX_list)

    if np.all(np.isclose(V, 0)):
        lbf = np.zeros(num_variables, dtype=float_type)
        if return_moments:
            post_mean = np.zeros((num_pops, num_variables), dtype=float_type)
            post_mean2 = np.zeros((num_pops, num_pops, num_variables), dtype=float_type)

    elif np.any(np.isclose(V, 0)):
        nonzero_pops = np.flatnonzero(~np.isclose(V, 0))
        lbf_out = compute_lbf(
            V = V[nonzero_pops], 
            XTY_list = [XTY_list[i] for i in nonzero_pops], 
            XTX_list = [XTX_list[i] for i in nonzero_pops], 
            X_l2_arr = X_l2_arr[nonzero_pops], 
            rho = rho[nonzero_pops, :][:, nonzero_pops],
            residual_variance = residual_variance[nonzero_pops], 
            return_moments = return_moments, 
            float_type = float_type
        )
        if return_moments:
            post_mean = np.zeros((num_pops, num_variables), dtype=float_type)
            post_mean2 = np.zeros((num_pops, num_pops, num_variables), dtype=float_type)
            lbf = lbf_out[0]
            post_mean[nonzero_pops] = lbf_out[1]
            post_mean2[np.ix_(nonzero_pops, nonzero_pops)] = lbf_out[2]
        else:
            lbf = lbf_out

    else:
        XTY = np.ascontiguousarray(np.stack(XTY_list, axis = 1))
        if return_moments:
            try:
                lbf, post_mean, post_mean2 = compute_lbf_and_moments(
                    V = V, XTY = XTY, X_l2_arr = X_l2_arr, rho = rho,
                    residual_variance = residual_variance, float_type = float_type
                )
            except:
                lbf, post_mean, post_mean2 = compute_lbf_and_moments_safe(
                    V = V, XTY = XTY, X_l2_arr = X_l2_arr, rho = rho,
                    residual_variance = residual_variance, float_type = float_type
                )
        else:
            lbf = compute_lbf_no_moments(
                V = V, XTY = XTY, X_l2_arr = X_l2_arr, rho = rho,
                residual_variance = residual_variance, float_type = float_type
            )
    if return_moments:
        return lbf, post_mean, post_mean2
    else:
        return lbf

@numba.jit(nopython=True, cache=True)
def compute_lbf_no_moments(
    V, XTY, X_l2_arr, rho, residual_variance, float_type = np.float32):
    
    
    num_pops = XTY.shape[1]
    num_variables = XTY.shape[0]

    lbf = np.zeros(num_variables, dtype=float_type)
    
    YT_invD_Z = XTY / residual_variance

    #compute a (the effects covariance matrix) and its inverse and log-determinant
    A = rho * np.sqrt(np.outer(V, V))
    inv_A = np.linalg.inv(A)
    logdet_A_sign, logdetA = np.linalg.slogdet(A)

    for i in range(num_variables):
    
        #compute the diagonal of q = z.t * inv(d) * z (this is a diagonal matrix)
        # this is (n_pop)
        Q_diag = X_l2_arr[:,i]/residual_variance

        #compute log-determinent for inv(a)+q
        # this is (n_pop, n_pop)
        Ainv_plus_Q = inv_A + np.diag(Q_diag)
        logdet_Ainv_plus_Q_sign, logdet_Ainv_plus_Q = np.linalg.slogdet(Ainv_plus_Q)
        assert logdet_Ainv_plus_Q_sign>0
    
        #compute inv_ainv_plus_q_times_zt_invd_y
        inv_Ainv_plus_Q_times_ZT_invD_Y = np.linalg.solve(Ainv_plus_Q, YT_invD_Z[i,:])

        #compute log-bf for this variable
        lbf_1 = 0.5 * YT_invD_Z[i,:].dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
        lbf_2 = -0.5 * (logdetA + logdet_Ainv_plus_Q)
        lbf[i] = lbf_1 + lbf_2
        
    return lbf

@numba.jit(nopython=True, cache=True)
def compute_lbf_and_moments(
    V, XTY, X_l2_arr, rho,
    residual_variance, float_type = np.float32):

    num_pops = XTY.shape[1]
    num_variables = XTY.shape[0]

    lbf = np.zeros(num_variables, dtype=float_type)
    post_mean = np.zeros((num_pops, num_variables), dtype=float_type)
    post_mean2 = np.zeros((num_pops, num_pops, num_variables), dtype=float_type)
        

    YT_invD_Z = XTY / residual_variance

    #compute a (the effects covariance matrix) and its inverse and log-determinant
    A = rho * np.sqrt(np.outer(V, V))
    inv_A = np.linalg.inv(A)
    logdet_A_sign, logdetA = np.linalg.slogdet(A)

    for i in range(num_variables):
    
        #compute the diagonal of q = z.t * inv(d) * z (this is a diagonal matrix)
        Q_diag = X_l2_arr[:,i] / residual_variance

        #compute log-determinent for inv(a)+q
        Ainv_plus_Q = inv_A + np.diag(Q_diag)
        logdet_Ainv_plus_Q_sign, logdet_Ainv_plus_Q = np.linalg.slogdet(Ainv_plus_Q)
        assert logdet_Ainv_plus_Q_sign>0
    
        #compute inv_ainv_plus_q_times_zt_invd_y
        inv_Ainv_plus_Q_times_ZT_invD_Y = np.linalg.solve(Ainv_plus_Q, YT_invD_Z[i, :])

        #compute log-bf for this variable
        lbf_1 = 0.5 * YT_invD_Z[i, :].dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
        lbf_2 = -0.5 * (logdetA + logdet_Ainv_plus_Q)
        lbf[i] = lbf_1 + lbf_2

        #compute posterior moments for this variable
        AQ = A*Q_diag
        post_mean[:,i] = A.dot(YT_invD_Z[i, :]) - AQ.dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
        post_covar_i = A - AQ.dot(A) + AQ.dot(np.linalg.solve(Ainv_plus_Q, AQ.T))
        post_mean2[:, :, i] = np.maximum(post_covar_i + np.outer(post_mean[:,i], post_mean[:,i]), 0)

    return lbf, post_mean, post_mean2

def compute_lbf_and_moments_safe(
    V, XTY, X_l2_arr, rho,
    residual_variance, float_type = np.float32):

    num_pops = XTY.shape[1]
    num_variables = XTY.shape[0]

    lbf = np.zeros(num_variables, dtype=float_type)
    post_mean = np.zeros((num_pops, num_variables), dtype=float_type)
    post_mean2 = np.zeros((num_pops, num_pops, num_variables), dtype=float_type)
        

    YT_invD_Z = XTY / residual_variance

    #compute a (the effects covariance matrix) and its inverse and log-determinant
    A = rho * np.sqrt(np.outer(V, V))
    inv_A = np.linalg.inv(A)
    logdet_A_sign, logdetA = np.linalg.slogdet(A)

    for i in range(num_variables):

        #compute the diagonal of q = z.t * inv(d) * z (this is a diagonal matrix)
        Q_diag = X_l2_arr[:,i]/residual_variance

        #compute log-determinent for inv(a)+q
        Ainv_plus_Q = inv_A + np.diag(Q_diag)
        logdet_Ainv_plus_Q_sign, logdet_Ainv_plus_Q = np.linalg.slogdet(Ainv_plus_Q)
        assert logdet_Ainv_plus_Q_sign>0

        #compute inv_ainv_plus_q_times_zt_invd_y
        inv_Ainv_plus_Q_times_ZT_invD_Y = np.linalg.solve(Ainv_plus_Q, YT_invD_Z[i, :])

        #compute log-bf for this variable
        lbf_1 = 0.5 * YT_invD_Z[i, :].dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
        lbf_2 = -0.5 * (logdetA + logdet_Ainv_plus_Q)
        lbf[i] = lbf_1 + lbf_2

        #compute posterior moments for this variable
        AQ = A*Q_diag
        post_mean[:,i] = A.dot(YT_invD_Z[i, :]) - AQ.dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
        post_covar_i = A - AQ.dot(A) + AQ.dot(np.linalg.solve(Ainv_plus_Q, AQ.T))
        post_mean2[:, :, i] = post_covar_i + np.outer(post_mean[:,i], post_mean[:,i])

    return lbf, post_mean, post_mean2

def get_objective(XTX_list, XTY_list, s, YTY_list, X_l2_arr):
    return Eloglik(XTX_list, XTY_list, s, YTY_list, X_l2_arr) - np.sum(s.KL)

def Eloglik(XTX_list, XTY_list, s, YTY_list, X_l2_arr):
    result = -0.5 * s.n.dot(np.log(2*np.pi*s.sigma2))
    for i in range(len(XTX_list)):
        result -= 0.5/s.sigma2[i] * s.ER2[i]
    return result

def update_ER2(XTX_list, XTY_list, YTY_list, s, X_l2_arr):
    ER2 = np.zeros(len(XTX_list))
    for i in range(len(XTX_list)):
        ER2[i] = get_ER2(
            XTX_list[i], XTY_list[i], YTY_list[i], s.alpha, s.mu[i], 
            s.mu2[i, i], X_l2_arr[i]
        )
    return ER2

def get_ER2(XTX, XTY, YTY, alpha, mu, mu2, X_l2):
    B = alpha * mu # beta should be lxp
    XB2 = np.sum(np.dot(B, XTX) * B)
    betabar = np.sum(B, axis = 0)
    postb2 = alpha * mu2

    result = YTY - 2 * betabar.dot(XTY) + betabar.dot(np.dot(XTX, betabar)) - XB2 + np.sum(X_l2 * postb2) 
    return result     

def SER_posterior_e_loglik(X_l2_arr, XTY_list, s2, Eb, Eb2):
    result = 0
    for i in range(len(X_l2_arr)):
        result -= .5 / s2[i] * (-2 * XTY_list[i].dot(Eb[i]) + X_l2_arr[i].dot(Eb2[i]))
    return result

def estimate_residual_variance_func(
        XTX_list, XTY_list, YTY_list, s, X_l2_arr, population_sizes, 
        float_type = np.float32):

    sigma2_arr = np.zeros(len(XTX_list), dtype=float_type)
    for i in range(len(XTX_list)):
        sigma2_arr[i] =  s.ER2[i] / population_sizes[i]
    return sigma2_arr

def susie_get_pip(s, prior_tol=1e-9):
    include_idx = np.any(s.V > prior_tol, axis = 0)
    if not np.any(include_idx):
        return np.zeros(s.alpha.shape[1])
    res = s.alpha[include_idx, :]
    pips = 1 - np.prod(1-res, axis=0)
    return pips

def susie_get_cs(
    s, R_list, coverage = 0.95, min_abs_corr = 0.5, dedup = True, n_purity = 100,
    calculate_purity = True, X_list = None):

    if len(s.V.shape) == 2:
        include_mask = np.any(s.V > 1e-9, axis = 0)
    else:
        include_mask = s.V > 1e-9

    status = in_CS(s.alpha, coverage)
    cs = [np.argwhere(status[i, :] == 1).flatten() for i in range(status.shape[0])]
    cs = [cs[i] if include_mask[i] else [] for i in range(len(cs))]
    claimed_coverage = np.array([np.sum(s.alpha[i, cs[i]]) for i in range(len(cs))])
    include_mask = include_mask & [len(csi) > 0 for csi in cs]

    if dedup:
        cs_set = set()
        for i in range(len(cs)):
            if tuple(cs[i]) in cs_set:
                include_mask[i] = False
            elif include_mask[i]:
                cs_set.add(tuple(cs[i]))
    if not np.any(include_mask):
        return ([[] for i in range(len(include_mask))], None, None, include_mask)

    if calculate_purity:
        purity = np.array([get_purity_x(cs[i], R_list, min_abs_corr, n_purity, X_list) if include_mask[i] else np.nan for i in range(len(cs))])
        include_mask[purity < min_abs_corr] = False
    else:
        purity = np.array([np.nan for i in range(len(cs))])

    return [cs, purity, claimed_coverage, include_mask]

def in_CS(alpha, coverage):
    return np.apply_along_axis(in_CS_x, 1, alpha, coverage)

def in_CS_x(alpha, coverage):
    n = n_in_CS_x(alpha, coverage)
    o = np.argsort(alpha)[::-1]
    result = np.zeros(alpha.shape, dtype=np.int32)
    result[o[:n]] = 1
    return result

def n_in_CS(alpha, coverage):
    return np.apply_along_axis(n_in_CS_x, 1, alpha, coverage)

def n_in_CS_x(alpha, coverage):
    return np.sum(np.cumsum(np.sort(alpha)[::-1]) < coverage) + 1

def get_purity_x(cs, R_list, min_abs_cor, n_purity, X_list):
    if len(cs) > n_purity:
        cs = random.sample(cs.tolist(), n_purity)
    if R_list is None:
        R_list = [np.corrcoef(X_list[i][:,cs], rowvar = False) for i in range(len(X_list))]
    else:
        R_list = [R[cs, :][:, cs] for R in R_list]
    abs_meta_R = functools.reduce(np.maximum, [np.abs(R) for R in R_list])

    return np.min(abs_meta_R)

def recover_R_from_XTX(XTX, X_l2):
    assert((XTX[:,np.flatnonzero(X_l2 == 0)] == 0).all())
    assert((XTX[np.flatnonzero(X_l2 == 0),:] == 0).all())
    with np.errstate(divide='ignore', invalid = 'ignore'):
        XTX /= np.sqrt(np.nan_to_num(X_l2, 1))
        XTX /= np.sqrt(np.expand_dims(np.nan_to_num(X_l2, 1), 1))

