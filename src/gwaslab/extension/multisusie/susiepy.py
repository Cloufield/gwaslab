import numpy as np; np.set_printoptions(precision=3)
from tqdm import tqdm #
import scipy.stats as stats
import scipy.linalg as la
import logging
import time
import sys
from scipy.optimize import minimize_scalar
from . import susiepy_ss

#TODO:
# - implement early_EM


class S:
    def __init__(self, X_std_list, L, scaled_prior_variance, residual_variance, varY, prior_weights, float_type):
        
        #code from init_setup
        num_pop = len(X_std_list)
        p = X_std_list[0].shape[1]
        self.alpha = np.zeros((L, p), dtype=float_type) + 1.0/p
        self.mu = np.zeros((num_pop, L, p), dtype=float_type)
        #self.mu2 = np.zeros((num_pop, L, p), dtype=float_type)
        self.mu2 = np.zeros((num_pop, num_pop, L, p), dtype=float_type)
        self.Xr_list = [np.zeros(X.shape[0], dtype=float_type) for X in X_std_list]
        self.sigma2 = residual_variance.astype(float_type)
        ###assert np.isscalar(self.sigma2)
        self.pi = prior_weights.astype(float_type)
        self.n = np.array([X.shape[0] for X in X_std_list])
        
        #code from init_finalize
        self.V = scaled_prior_variance * varY + np.zeros((num_pop, L), dtype=float_type)
        self.V = self.V.astype(float_type)
        assert np.all(self.V >= 0)
        self.KL = np.zeros(L, dtype=float_type) + np.nan
        self.lbf = np.zeros(L, dtype=float_type) + np.nan
        
        self.converged = False

class SER_RESULTS:
            def __init__(self, alpha, mu, mu2, lbf, lbf_model, V, loglik):
                self.alpha = alpha
                self.mu = mu
                self.mu2 = mu2
                self.lbf = lbf
                self.lbf_model = lbf_model
                self.V = V
                self.loglik = loglik

def multisusie(
    X_list, 
    Y_list, 
    rho, 
    L,
    scaled_prior_variance = 0.2,
    prior_weights = None,
    standardize = False,
    residual_variance = None, 
    estimate_residual_variance = True,
    estimate_prior_variance = True,
    estimate_prior_method = 'early_EM',
    pop_spec_effect_priors = True,
    iter_before_zeroing_effects = 5,
    prior_tol = 1e-9,
    max_iter = 100,
    residual_variance_upperbound = np.inf,
    residual_variance_lowerbound = 0,
    tol = 1e-3,
    verbose = False,
    coverage = 0.95,
    min_abs_corr = 0,
    check_null_threshold = 0,
    float_type = np.float64,
    low_memory_mode = False,
    calculate_purity = True,
    n_purity = 100,
    variant_ids = None
    ):

    """ Top-level function for running MultiSuSiE with individual level data

    This function takes genotype and phenotype matrices and runs MultiSuSiE on them

    Parameters
    ----------
    X_list: Length K list of numpy arrays, one for each population. Rows correspond
        to samples and the columns correspond to variants. Each array should
        contain the same set of variants in the same order. Its fine if 
        some columns are constant. 
    Y_list: Length K list of one dimensional numpy arrays, one for each population. 
        Rows correspond to samples. Samples should be in the same order as in 
        X_list. 
    rho: PxP numpy array representing the effect size correlation matrix (P is the
        number of variants). In the manuscript, we show that this parameter has 
        little impact on the estimated PIPs in practice.
    L: integer representing the maximum number of causal variants
    scaled_prior_variance: float representing the effect size prior variance,
        scaled by the residual variance. It's fine to set this to a number 
        larger than what you expect the squared effect size to be (like the 
        default value of 0.2) as long as estimate_prior_variance is set to True
        and estimate_prior_method is not set to None.
    prior_weights: numpy P-array of floats representing the prior probability
        of causality for each variant. Give None to use a uniform prior
    standardize: boolnea, whether to standardize the genotypes to have
        variance of 1.
    residual_variance: Length K numpy array of floats representing the residual
        variance for each population.
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
    residual_variance_upperbound: float, upper bound on the residual variance
    residual_variance_lowerbound: float, lower bound on the residual variance
    tol: float, after iter_before_zeroing_effects iterations, results
        are returned if the ELBO increases by less than tol in an ieration
    verbose: boolean which indicates if the objective function should be printed
    coverage: float representing the minimum coverage of credible sets
    min_abs_corr: float representing the minimum absolute correlation between
        any pair of variants in a credible set (purity). For each pair of variants,
        the max is taken across ancestries. In the case where min_abs_corr = 0,
        low_memory_mode = True, and recover_R = False, the purity of credible
        sets will not be calculated. 
    variant_ids: length P list of strings representing the variant IDs. If 
        provided, sets will contain a fifth entry, containing the variant ids
        of the variants contained in each set. 

    """

    t0 = time.time()
    #check input
    assert len(X_list) == len(Y_list)
    assert np.all([X.shape[1] == X_list[0].shape[1] for X in X_list])
    assert np.all([X.shape[0] == Y.shape[0] for X,Y in zip(X_list, Y_list)])
    if prior_weights is not None:
        prior_weights = prior_weights.astype(float_type)
    
    assert not np.any([np.any(np.isnan(X)) for X in X_list])
    
    #remove missing individuals
    for i in range(len(X_list)):
        if np.any(np.isnan(Y_list[i])):
            X_list[i] = X_list[i][~np.isnan(Y_list[i])]
            Y_list[i] = Y_list[i][~np.isnan(Y_list[i])]
        
    #center Y
    mean_y_list = [Y.mean() for Y in Y_list]
    Y_list = [Y-mean_Y for Y,mean_Y in zip(Y_list, mean_y_list)]
        
    #compute w_pop (the relative size of each population)
    n_arr = np.array([X.shape[0] for X in X_list], dtype=int)
    w_pop = (n_arr / n_arr.sum()).astype(float_type)

    #compute rho properties
    rho = rho.astype(float_type)
    logdet_rho_sign, logdet_rho = np.linalg.slogdet(rho)
    assert logdet_rho_sign>0

    #compute X mean and std
    cm_arr = np.array([X.mean(axis=0) for X in X_list], dtype=float_type)
    if standardize:
        csd_arr = np.array([X.std(axis=0, ddof=1) for X in X_list], dtype=float_type)
        csd = np.sum(csd_arr * w_pop[:, np.newaxis], axis=0)
    else:
        csd = np.ones(X_list[0].shape[1], dtype=float_type)

        
    #Standardize X
    is_constant_column = np.isclose(csd, 0.0)
    csd[is_constant_column] = 1.0
    if low_memory_mode:
        X_std_list = X_list
        for pop_i in range(len(X_std_list)):
            X_std_list[pop_i] -= cm_arr[pop_i]
            X_std_list[pop_i] /= csd
    else:
        X_std_list = [(X-cm)/csd for X,cm in zip(X_list, cm_arr)]
    
    #explicitly mark that constant columns are zero
    if np.any(is_constant_column):
        for i in range(len(X_std_list)):
            X_std_list[i][:, is_constant_column] = 0.0
        
    
    #create a C-contiguous version of X
    assert np.all([X.flags['F_CONTIGUOUS'] == X_std_list[0].flags['F_CONTIGUOUS'] for X in X_std_list])
    if X_std_list[0].flags['F_CONTIGUOUS']:
        for i in range(len(X_std_list)):
            X_std_list[i] = np.ascontiguousarray(X_std_list[i])
        
    
    X_l2_arr = np.array([np.einsum('ij,ij->j', X_std, X_std) for X_std in X_std_list], dtype=float_type)
    
    
    #compute rho properties
    rho = rho.astype(float_type)
    
    #init setup
    p = X_list[0].shape[1]
    varY = np.concatenate(Y_list).var(ddof=1)
    varY_list = [Y.var(ddof=1) for Y in Y_list]
    if np.isscalar(scaled_prior_variance):
        assert 0 < scaled_prior_variance <= 1
    if residual_variance is None:
        residual_variance = np.array(varY_list, dtype=float_type)
    if prior_weights is None:
        prior_weights = np.zeros(p, dtype=float_type) + 1.0/p
    else:
        prior_weights = (prior_weights / np.sum(prior_weights)).astype(float_type)
    assert prior_weights.shape[0] == p
    if p<L: L=p
    s = S(X_std_list, L, scaled_prior_variance, residual_variance, varY, prior_weights, float_type = float_type)
    elbo = np.zeros(max_iter+1) + np.nan
    elbo[0] = -np.inf
    
    
    ### start iterations ###
    tqdm_iter = tqdm(list(range(max_iter)), disable=not verbose, file=sys.stdout)
    for i in tqdm_iter:
        tqdm_iter.set_description('iteration %d/%d'%(i+1, max_iter))

        if estimate_prior_method == 'early_EM':
            if i == 0: 
                current_estimate_prior_method = None
            else:
                current_estimate_prior_method = 'early_EM'
        else:
            current_estimate_prior_method = estimate_prior_method

        if i < iter_before_zeroing_effects:
            current_check_null_threshold = -np.inf
        else:
            current_check_null_threshold = check_null_threshold

        s = update_each_effect(
            X_std_list = X_std_list, Y_list = Y_list, s = s, X_l2_arr = X_l2_arr, 
            w_pop = w_pop, rho = rho, estimate_prior_variance = estimate_prior_variance, 
            estimate_prior_method = current_estimate_prior_method, verbose=verbose,
            float_type = float_type, pop_spec_effect_priors = pop_spec_effect_priors,
            check_null_threshold = current_check_null_threshold
        )

        #compute objective before updating residual variance
        #because part of the objective s.kl has already been computed
        #under the residual variance before the update
        elbo[i+1] = get_objective(X_std_list, Y_list, s, X_l2_arr)
        if verbose:
            logging.info('objective: %s'%(elbo[i+1]))
            print('objective: %s'%(elbo[i+1]))
        
        if (elbo[i+1] - elbo[i]) < tol:
            s.converged = True
            tqdm_iter.close()
            break
        
        if estimate_residual_variance:
            s.sigma2 = estimate_residual_variance_func(X_std_list, Y_list, s, X_l2_arr, float_type)
            s.sigma2 = np.minimum(s.sigma2, residual_variance_upperbound).astype(float_type)
            s.sigma2 = np.maximum(s.sigma2, residual_variance_lowerbound).astype(float_type)
        if verbose:
            logging.info('objective after updating sigma2: %s'%(get_objective(X_std_list, Y_list, s, X_l2_arr)))
    
    elbo = elbo[1:i+2] # Remove first (infinite) entry, and trailing NAs.
    s.elbo = elbo
    s.niter = i+1

    if verbose: 
        logging.info('done in %0.2f seconds'%(time.time() - t0))
        logging.info(f'elbo: {elbo[~np.isnan(elbo)]}')
    if not s.converged:
        logging.info('IBSS algorithm did not converge in %d iterations'%(max_iter))
    else:
        logging.info('IBSS algorithm converged in %d iterations'%(i))
        
    #zero out everything related to constant variables, just to be on the safe side
    if np.any(is_constant_column):
        s.mu[:, :, is_constant_column] = 0.0
        s.mu2[:, :, :, is_constant_column] = 0.0
        s.alpha[:, is_constant_column] = 0.0
        
    s.intercept = np.zeros(len(X_list))
    s.fitted = s.Xr_list
        
    s.pip = susie_get_pip(s, prior_tol=prior_tol)    
    s.X_column_scale_factors = csd.copy()
    s.X_column_scale_factors[is_constant_column] = 0.0
    s.coef = np.array([np.sum(s.mu[k] * s.alpha, axis=0) / csd for k in range(len(Y_list))])
    s.coef_sd = np.array([(np.sqrt(np.sum(s.alpha * s.mu2[k, k] - (s.alpha*s.mu[k])**2, axis=0)) / csd) for k in range(len(Y_list))])

    s.sets = susiepy_ss.susie_get_cs(
        s, R_list = None, coverage = coverage, min_abs_corr = min_abs_corr,
        dedup = True, n_purity = n_purity, calculate_purity = calculate_purity,
        X_list = X_std_list
    )

    s.variant_ids = variant_ids
    if variant_ids is not None:
        cs_variant_ids = []
        for i in range(len(s.sets[0])):
            cs_variant_ids.append([variant_ids[idx] for idx in s.sets[0][i]])
        s.sets.append(cs_variant_ids)

    return s


def get_ER2(X_std, Y, alpha, mu, mu2, Xr, X_l2):
    """ expected squared residuals
    """
    Xr_L = X_std.dot((alpha*mu).T)
    postb2 = alpha*mu2
    r = Y - Xr
    result = r.dot(r) - np.einsum('ij,ij->',Xr_L, Xr_L) + np.sum(X_l2.dot(postb2.T))
    return result     

def SER_posterior_e_loglik(X_std_list, Y_list, s2, Eb, Eb2, X_l2_arr, n):
    """ posterior expected loglikelihood for a SER (Eq B.6 - B.9 (after expanding the L2-norm))

    Eb the posterior mean of b (p vector) (alpha * mu)
    Eb2 the posterior second moment of b (p vector) (alpha * mu2)
    """
    result = -0.5 * n.dot(np.log(2*np.pi*s2))
    for i in range(len(Y_list)):
        Y = Y_list[i]
        X_std = X_std_list[i]
        result -= 0.5/s2[i] * (Y.dot(Y) - 2*Y.dot(X_std.dot(Eb[i])) + X_l2_arr[i].dot(Eb2[i]))
        #result -= 0.5/s2[i] * (np.sum(Y * Y) - 2*Y.dot(X_std.dot(Eb[i])) + X_l2_arr[i].dot(Eb2[i]))
        #result -= 0.5/s2[i] * (np.sum(Y * Y) - 2*Y.dot(X_std.dot(Eb[i])) + np.sum(X_l2_arr[i] * Eb2[i]))
    return result


def Eloglik(X_std_list, Y_list, s, X_l2_arr):
    """expected log-likelihood for a susie fit
    """
    result = -0.5 * s.n.dot(np.log(2*np.pi*s.sigma2))
    for i in range(len(Y_list)):
        result -= 0.5/s.sigma2[i] * get_ER2(X_std_list[i], Y_list[i], s.alpha, s.mu[i], s.mu2[i, i], s.Xr_list[i], X_l2_arr[i])
    return result

def get_objective(X_std_list, Y_list, s, X_l2_arr):
    return Eloglik(X_std_list, Y_list, s, X_l2_arr) - np.sum(s.KL)

def estimate_residual_variance_func(X_std_list, Y_list, s, X_l2_arr, float_type):
    sigma2_arr = np.zeros(len(Y_list), dtype=float_type)
    for i in range(len(Y_list)):
        sigma2_arr[i] =  get_ER2(X_std_list[i], Y_list[i], s.alpha, s.mu[i], s.mu2[i, i], s.Xr_list[i], X_l2_arr[i]) / s.n[i]
    return sigma2_arr

def loglik_f(V, prior_weights, compute_lbf_params):
    lbf = compute_lbf(V, *compute_lbf_params)
    maxlbf = np.max(lbf)
    w = np.exp(lbf - maxlbf)
    w_weighted = w * prior_weights
    weighted_sum_w = w_weighted.sum()
    loglik = maxlbf + np.log(weighted_sum_w)
    return loglik

#def optimize_prior_variance(optimize_V, prior_weights, compute_lbf_params=None, alpha=None, post_mean2=None, w_pop=None, check_null_threshold=0, float_type = np.float64):
#    
#    if optimize_V == 'optim':
#        if pop_spec_effect_priors:
#            raise Exception('estimate_prior_method="optim" with ' +
#        neg_loglik_logscale = lambda lV: -loglik_f(np.exp(lV), prior_weights, compute_lbf_params)
#        opt_obj = minimize_scalar(neg_loglik_logscale, bounds=(-30,15))
#        lV = opt_obj.x
#        V = np.exp(lV)
#    elif optimize_V == 'EM':
#        V_arr = np.array([np.sum(alpha * post_mean2[i]) for i in range(post_mean2.shape[0])], dtype=float_type)
#        if not pop_spec_effect_priors:
#            V = (w_pop.dot(V_arr)).astype(float_type)
#    else:
#        raise ValueError('unknown optimization method')
#        
#    # set V exactly 0 if that beats the numerical value by check_null_threshold in loglik.
#    #if check_null_threshold>0:
#    if not pop_spec_effect_priors:
#
#    if loglik_f(0, prior_weights, compute_lbf_params) + check_null_threshold >= loglik_f(V, prior_weights, compute_lbf_params):
#        V=0
#    return V
    
def compute_lbf_1pop(V, X_std, Y, X_l2,
        residual_variance,
        return_moments=False,
        verbose=False,
        float_type = np.float64
        ):
    
    Xty = Y.dot(X_std)
    betahat = np.zeros(Xty.shape[0], dtype=float_type)
    shat2 = np.zeros(Xty.shape[0], dtype=float_type)
    lbf = np.zeros(Xty.shape[0], dtype=float_type)
    nz = (X_l2!=0)
    betahat[nz] = Xty[nz] / X_l2[nz]
    shat2[nz] = residual_variance / X_l2[nz]
    lbf[nz] = stats.norm(0, np.sqrt(V+shat2[nz])).logpdf(betahat[nz]) - stats.norm(0, np.sqrt(shat2[nz])).logpdf(betahat[nz]) #equation A.3
    if not return_moments: return lbf
    
    post_var = np.zeros((1, Xty.shape[0]), dtype=float_type)
    if not np.isclose(V, 0):
        post_var[0,nz] = 1.0 / (1.0/V + 1.0/shat2[nz]) # posterior variance
    post_mean = (1.0/residual_variance) * post_var * Xty
    post_mean2 = post_var + post_mean**2 # second moment

    return lbf, post_mean, post_mean2
    
def compute_lbf(V, Y_list, X_std_list, X_l2_arr,
        rho,
        residual_variance,
        return_moments=False,
        verbose=False,
        float_type=np.float64
        ):

    num_pops = len(Y_list)
    
    
    num_variables = X_std_list[0].shape[1]
    lbf = np.zeros(num_variables, dtype=float_type)
    
    if return_moments:
        post_mean = np.zeros((num_pops, num_variables), dtype=float_type)
        post_mean2 = np.zeros((num_pops, num_pops, num_variables), dtype=float_type)

    if np.all(np.isclose(V, 0)):
        pass
    elif np.any(np.isclose(V, 0)):
        nonzero_pops = np.flatnonzero(~np.isclose(V, 0))
        lbf_out = compute_lbf(
            V = V[nonzero_pops], 
            Y_list = [Y_list[i] for i in nonzero_pops], 
            X_std_list = [X_std_list[i] for i in nonzero_pops], 
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
        #compute YT_invD_Z = Y.T * inv(D) * Z - i.e., the inner product of Y on all X in each population, scaled by 1/sigma2 for that population
        YT_invD_Z = np.array([Y.dot(X)/sigma2 for Y,X,sigma2 in zip(Y_list, X_std_list, residual_variance)], dtype=float_type)
    
        #compute A (the effects covariance matrix) and its inverse and log-determinant
        A = rho * np.sqrt(np.outer(V, V))
        inv_A = np.linalg.inv(A)
        logdet_A_sign, logdetA = np.linalg.slogdet(A)
    
        for i in range(num_variables):
        
            #compute the diagonal of Q = Z.T * inv(D) * Z (this is a diagonal matrix)
            Q_diag = X_l2_arr[:,i] / residual_variance

            #compute log-determinent for inv(A)+Q
            Ainv_plus_Q = inv_A + np.diag(Q_diag)
            logdet_Ainv_plus_Q_sign, logdet_Ainv_plus_Q = np.linalg.slogdet(Ainv_plus_Q)
            assert logdet_Ainv_plus_Q_sign>0
        
            #compute inv_Ainv_plus_Q_times_ZT_invD_Y
            inv_Ainv_plus_Q_times_ZT_invD_Y = np.linalg.solve(Ainv_plus_Q, YT_invD_Z[:,i])

            #compute log-BF for this variable
            lbf_1 = 0.5 * YT_invD_Z[:,i].dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
            lbf_2 = -0.5 * (logdetA + logdet_Ainv_plus_Q)
            lbf[i] = lbf_1 + lbf_2
        
            #compute posterior moments for this variable
            if return_moments:
                AQ = A*Q_diag
                post_mean[:,i] = A.dot(YT_invD_Z[:,i]) - AQ.dot(inv_Ainv_plus_Q_times_ZT_invD_Y)
                post_covar_i = A - AQ.dot(A) + AQ.dot(np.linalg.solve(Ainv_plus_Q, AQ.T)) # AQ is symmetric...
                post_mean2[:,:,i] = np.maximum(post_covar_i + np.outer(post_mean[:,i], post_mean[:,i]), 0)
        
    if return_moments:
        return lbf, post_mean, post_mean2
    else:
        return lbf
        
def single_effect_regression(Y_list, X_std_list, V, X_l2_arr, w_pop,
        rho,
        residual_variance, 
        prior_weights=None,
        optimize_V=None,
        check_null_threshold=0,
        pop_spec_effect_priors = True,
        alpha = None,
        mu2 = None,
        verbose=False,
        float_type = np.float64
        ):

    
    #optimize V if needed (V is sigma_0^2 in the paper)
    compute_lbf_params = (Y_list, X_std_list, X_l2_arr, rho, residual_variance, False, verbose, float_type)
    if optimize_V not in ['EM', None]:
        V = susiepy_ss.optimize_prior_variance(
            optimize_V = optimize_V, prior_weights = prior_weights, num_pops = rho.shape[0], 
            compute_lbf_params=compute_lbf_params, alpha=alpha, post_mean2=mu2, w_pop=w_pop, 
            check_null_threshold=check_null_threshold, pop_spec_effect_priors = pop_spec_effect_priors, 
            float_type = float_type, loglik_function = loglik_f)
        
    #compute lbf (log Bayes-factors)
    lbf, post_mean, post_mean2 = compute_lbf(V, Y_list, X_std_list, X_l2_arr, rho, residual_variance, return_moments=True, verbose=verbose)
    
    #compute alpha as defined in Appendix A.2
    maxlbf = np.max(lbf)
    w = np.exp(lbf - maxlbf)
    w_weighted = w * prior_weights
    weighted_sum_w = w_weighted.sum()
    alpha = w_weighted / weighted_sum_w
    
    #compute log-likelihood (equation A.5)
    lbf_model = maxlbf + np.log(weighted_sum_w)
    loglik = lbf_model + np.sum([np.sum(stats.norm(0, np.sqrt(sigma2)).logpdf(Y)) for sigma2,Y in zip(residual_variance, Y_list)])
        
    if optimize_V == 'EM':
        V = susiepy_ss.optimize_prior_variance(
            optimize_V = optimize_V, prior_weights = prior_weights, num_pops = rho.shape[0],
            compute_lbf_params=compute_lbf_params, alpha=alpha, post_mean2=post_mean2, w_pop=w_pop, 
            check_null_threshold=check_null_threshold, pop_spec_effect_priors = pop_spec_effect_priors,
            float_type=float_type, loglik_function = loglik_f)
        
    res = SER_RESULTS(alpha=alpha, mu=post_mean, mu2=post_mean2, lbf=lbf, lbf_model=lbf_model, V=V, loglik=loglik)
    return res

def update_each_effect(X_std_list, Y_list, s, X_l2_arr, w_pop,
        rho,
        estimate_prior_variance=False, 
        estimate_prior_method='optim',
        check_null_threshold=0.0,
        verbose=False,
        pop_spec_effect_priors = True,
        float_type = np.float64
        ):

        
    if not estimate_prior_variance:
        estimate_prior_method = None
    L = s.alpha.shape[0]
    num_pop = rho.shape[0]
    
    tqdm_L = tqdm(range(L), disable=not verbose, file=sys.stdout)
    for l in tqdm_L:
        tqdm_L.set_description('csnp %d/%d'%(l+1, L))
        R_list = []
        for k in range(len(X_std_list)):
            s.Xr_list[k] -= X_std_list[k].dot(s.alpha[l] * s.mu[k,l])
            R_list.append(Y_list[k] - s.Xr_list[k])
        res = single_effect_regression(
            R_list, X_std_list, s.V[:,l], X_l2_arr, w_pop,
            rho, residual_variance=s.sigma2, prior_weights=s.pi,
            optimize_V=estimate_prior_method, check_null_threshold=check_null_threshold,
            pop_spec_effect_priors = pop_spec_effect_priors, verbose=verbose, 
            alpha = s.alpha[l,:], mu2 = s.mu2[:,:,l,:],
            float_type = float_type
        )
              
        # Update the variational estimate of the posterior mean.
        s.mu[:,l,:] = res.mu
        s.alpha[l,:] = res.alpha
        s.mu2[:,:,l,:] = res.mu2
        s.V[:,l] = res.V
        s.lbf[l] = res.lbf_model
        s.KL[l] = -res.loglik + SER_posterior_e_loglik(X_std_list, R_list, s.sigma2, res.mu * res.alpha, res.mu2[range(num_pop), range(num_pop)] * res.alpha, X_l2_arr, s.n)
        for k in range(len(X_std_list)):
            s.Xr_list[k] += X_std_list[k].dot(s.alpha[l] * s.mu[k,l])

        
    return s
        
def susie_get_pip(s, prior_tol=1e-9):
    include_idx = np.any(s.V > prior_tol, axis = 0)
    if not np.any(include_idx):
        return np.zeros(s.alpha.shape[1])
    res = s.alpha[include_idx, :]
    pips = 1 - np.prod(1-res, axis=0)
    return pips
        
susie_multi = multisusie
