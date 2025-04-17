#!/usr/bin/env python

"""
PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel.

Reference: T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.
           Nature Communications, 10:1776, 2019.


Usage:
python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR
                [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR
                 --chrom=CHROM --write_psi=WRITE_PSI --write_pst=WRITE_POSTERIOR_SAMPLES --seed=SEED]

"""


import os
import sys
import getopt

import gwaslab.prscs_parse_genet as parse_genet
import gwaslab.prscs_mcmc_gtb as mcmc_gtb
import gwaslab.prscs_gigrnd as gigrnd


def _run_prscs(
         ref_dir=None,
         bim_prefix=None,
         sst_file=None,
            a= 1, 
            b= 0.5, 
            phi= None, 
            n_gwas= None,
            n_iter= 1000, 
            n_burnin= 500, 
            thin= 5, 
            out_dir= "./", 
            chrom= range(1,23),
            beta_std= 'FALSE', 
            write_psi= 'FALSE', 
            write_pst= 'FALSE', 
            seed= None,
            log=None,
         **kwargs):
    ## T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.Nature Communications, 10:1776, 2019.
    sst_file = sst_file.rename(columns={"rsID":"SNP","POS":"BP"})
    log.write("Start to runnig PRScs...")
    param_dict = {'ref_dir': ref_dir, 
                    'bim_prefix': bim_prefix, 
                    'a': a, 
                    'b': b, 
                    'phi': phi, 
                    'n_gwas': n_gwas,
                    'n_iter': n_iter, 
                    'n_burnin': n_burnin, 
                    'thin':thin, 
                    'out_dir': out_dir, 
                    'chrom': chrom,
                    'beta_std': beta_std, 
                    'write_psi': write_psi, 
                    'write_pst': write_pst, 
                    'seed': seed}
    
    for chrom in param_dict['chrom']:
        log.write('##### process chromosome %d #####' % int(chrom))

        if '1kg' in os.path.basename(param_dict['ref_dir']):
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_1kg_hm3', int(chrom), log)
        elif 'ukbb' in os.path.basename(param_dict['ref_dir']):
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_ukbb_hm3', int(chrom), log)

        vld_dict = parse_genet.parse_bim(param_dict['bim_prefix'], int(chrom))

        sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, sst_file, param_dict['n_gwas'], log)

        ld_blk, blk_size = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict, int(chrom), log)

        mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, param_dict['n_gwas'], ld_blk, blk_size,
            param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], int(chrom), param_dict['out_dir'], param_dict['beta_std'],
	    param_dict['write_psi'], param_dict['write_pst'], param_dict['seed'], log)

    log.write("Finished!")

