#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.

"""
import numpy as np
from scipy import linalg 
from numpy import random
from gwaslab.extension.prscs.prscs_gigrnd2 import gigrnd
from gwaslab.info.g_Log import Log
import time
def mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, write_psi, write_pst, seed, log, debug=False):
    log.write('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)
    
    # Debug: initial state
    if debug:
        log.write('[DEBUG] Initialization: p={}, n_blk={}, n_iter={}, n_burnin={}'.format(
            len(sst_dict['SNP']), len(ld_blk), n_iter, n_burnin))

    # derived stats - use efficient array creation with explicit dtype
    # Original: np.array(..., ndmin=2).T
    # Optimized: np.asarray(...).reshape(-1, 1) is more efficient
    beta_mrg = np.asarray(sst_dict['BETA'], dtype=np.float64).reshape(-1, 1)
    maf = np.asarray(sst_dict['MAF'], dtype=np.float64).reshape(-1, 1)
    n_pst = int((n_iter - n_burnin) / thin)
    p = len(sst_dict['SNP'])
    n_blk = len(ld_blk)

    # Pre-compute maf scaling factor if needed (for per-allele conversion)
    # This avoids recomputing it later in the conversion step
    if beta_std == 'FALSE':
        maf_scale = np.sqrt(2.0 * maf * (1.0 - maf), dtype=np.float64)
    else:
        maf_scale = None

    # initialization - explicit dtype for numerical stability
    beta = np.zeros((p, 1), dtype=np.float64)
    psi = np.ones((p, 1), dtype=np.float64)
    sigma = np.float64(1.0)
    
    if phi == None:
        phi = np.float64(1.0)
        phi_updt = True
    else:
        phi = np.float64(phi)
        phi_updt = False

    if write_pst == 'TRUE':
        beta_pst = np.zeros((p, n_pst), dtype=np.float64)

    beta_est = np.zeros((p, 1), dtype=np.float64)
    psi_est = np.zeros((p, 1), dtype=np.float64)
    sigma_est = np.float64(0.0)
    phi_est = np.float64(0.0)

    # MCMC
    pp = 0
    start_time = time.time() 
    for itr in range(1,n_iter+1):
        if itr ==2:
            loop_time = time.time() - start_time 
            log.write(" -Estimated time: {} mins".format((loop_time*n_iter)/60))
        
        if itr % 100 == 0:
            log.write('--- iter-' + str(itr) + ' ---')
        elif itr % 100 > 2:
            log.write('-', end="", show_time=False)
        elif itr % 100 ==2:
            log.write('-', end="")

        mm = 0; quad = 0.0
        
        # Pre-compute sqrt(sigma/n) once per iteration for efficiency
        sqrt_sigma_n = np.sqrt(sigma / n)

        for kk in range(n_blk):
            if blk_size[kk] == 0:
                continue
            
            # Use slice for efficient numpy indexing instead of range
            blk_start = mm
            blk_end = mm + blk_size[kk]
            blk_len = blk_size[kk]
            idx_blk = slice(blk_start, blk_end)
            
            # Optimize diagonal matrix creation: add to diagonal directly
            # Original: dinvt = ld_blk[kk] + np.diag(1.0/psi[idx_blk].T[0])
            dinvt = ld_blk[kk].copy()  # Copy to avoid modifying original
            psi_blk = psi[idx_blk, 0]  # Extract 1D array efficiently
            # Protect against division by zero
            psi_inv = np.reciprocal(np.maximum(psi_blk, 1e-12))
            # Add to diagonal using diag_indices (more efficient than creating full diagonal matrix)
            diag_idx = np.diag_indices(blk_len)
            dinvt[diag_idx] += psi_inv
            
            # Cholesky decomposition
            dinvt_chol = linalg.cholesky(dinvt, check_finite=False)
            
            # Sample beta for this block
            beta_mrg_blk = beta_mrg[idx_blk]
            # Solve: dinvt_chol^T * x = beta_mrg_blk
            beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg_blk, trans='T')
            # Add noise (vectorized)
            beta_tmp += sqrt_sigma_n * random.randn(blk_len, 1)
            # Solve: dinvt_chol * beta = beta_tmp
            beta[idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
            
            # Compute quadratic form: beta^T * dinvt * beta
            # Original: quad += np.dot(np.dot(beta[idx_blk].T, dinvt), beta[idx_blk])
            # Optimized: use einsum for better performance (vectorized)
            beta_blk = beta[idx_blk, 0]  # Extract 1D array
            quad += np.einsum('i,ij,j', beta_blk, dinvt, beta_blk)
            
            mm += blk_len

        # Update sigma - use vectorized sum operations
        # Original: err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))
        n_half = n / 2.0
        beta_beta_mrg = np.sum(beta * beta_mrg)  # Vectorized sum
        beta_sq_psi = np.sum(beta ** 2 / psi)    # Vectorized sum
        err = max(n_half * (1.0 - 2.0 * beta_beta_mrg + quad), n_half * beta_sq_psi)

        # Update sigma - pre-compute (n+p)/2.0
        n_p_half = (n + p) / 2.0
        sigma = 1.0 / random.gamma(n_p_half, 1.0 / err)
        
        # Update delta - vectorized operation
        # Original: delta = random.gamma(a+b, 1.0/(psi+phi))
        a_plus_b = a + b
        psi_plus_phi = psi + phi
        delta = random.gamma(a_plus_b, 1.0 / psi_plus_phi)
        
        # Update psi - pre-compute values for efficiency (loop necessary due to gigrnd)
        # Original: for jj in range(p): psi[jj] = gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)
        a_minus_half = a - 0.5
        beta_sq_scaled = (beta[:, 0] ** 2) * n / sigma  # Pre-compute all values (vectorized)
        delta_flat = delta[:, 0]  # Extract 1D array for indexing
        for jj in range(p):
            psi[jj, 0] = gigrnd(a_minus_half, 2.0 * delta_flat[jj], beta_sq_scaled[jj])
        
        # Cap psi at 1.0 - use vectorized clip operation
        # Original: psi[psi>1] = 1.0
        np.clip(psi, None, 1.0, out=psi)
        
        if phi_updt == True:
            w = random.gamma(1.0, 1.0 / (phi + 1.0))
            delta_sum = np.sum(delta)  # Vectorized sum
            phi = random.gamma(p * b + 0.5, 1.0 / (delta_sum + w))
        
        # Debug: log intermediate values at key iterations
        if debug and (itr % 50 == 0 or itr <= 5):
            beta_sum = np.sum(beta)
            beta_mean = np.mean(beta)
            beta_std_val = np.std(beta)
            psi_mean = np.mean(psi)
            psi_min = np.min(psi)
            psi_max = np.max(psi)
            log.write('[DEBUG-ORIG] iter={}: sigma={:.6e}, phi={:.6e}, beta_sum={:.6e}, beta_mean={:.6e}, beta_std={:.6e}, psi_mean={:.6e}, psi_range=[{:.6e}, {:.6e}]'.format(
                itr, sigma, phi, beta_sum, beta_mean, beta_std_val, psi_mean, psi_min, psi_max))

        # posterior
        if (itr>n_burnin) and (itr % thin == 0):
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst

            if write_pst == 'TRUE':
                beta_pst[:,[pp]] = beta
                pp += 1

    # convert standardized beta to per-allele beta
    if beta_std == 'FALSE':
        # Use pre-computed maf_scale instead of recomputing
        beta_est /= maf_scale

        if write_pst == 'TRUE':
            beta_pst /= maf_scale


    # write posterior effect sizes
    if phi_updt == True:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
    else:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)

    with open(eff_file, 'w') as ff:
        if write_pst == 'TRUE':
            for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_pst):
                ff.write(('%d\t%s\t%d\t%s\t%s' + '\t%.6e'*n_pst + '\n') % (chrom, snp, bp, a1, a2, *beta))
        else:
            for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_est):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # write posterior estimates of psi
    if write_psi == 'TRUE':
        if phi_updt == True:
            psi_file = out_dir + '_pst_psi_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
        else:
            psi_file = out_dir + '_pst_psi_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)

        with open(psi_file, 'w') as ff:
            for snp, psi in zip(sst_dict['SNP'], psi_est):
                ff.write('%s\t%.6e\n' % (snp, psi))

    # print estimated phi
    if phi_updt == True:
        log.write('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )
    
    # Debug: final state
    if debug:
        beta_est_sum = np.sum(beta_est)
        beta_est_mean = np.mean(beta_est)
        beta_est_std = np.std(beta_est)
        psi_est_mean = np.mean(psi_est)
        log.write('[DEBUG-ORIG] Final: beta_est_sum={:.6e}, beta_est_mean={:.6e}, beta_est_std={:.6e}, psi_est_mean={:.6e}, sigma_est={:.6e}, phi_est={:.6e}'.format(
            beta_est_sum, beta_est_mean, beta_est_std, psi_est_mean, sigma_est, phi_est))

    log.write('... Done ...')
