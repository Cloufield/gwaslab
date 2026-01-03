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

from typing import TYPE_CHECKING, Optional, Union, Any
import os
import sys
import getopt
import pandas as pd
from multiprocessing import Pool, cpu_count

import gwaslab.extension.prscs.prscs_parse_genet as parse_genet
import gwaslab.extension.prscs.prscs_mcmc_gtb as mcmc_gtb
import gwaslab.extension.prscs.precs_mcmc_gtb2 as mcmc_gtb2
import gwaslab.extension.prscs.prscs_gigrnd as gigrnd

if TYPE_CHECKING:
    from gwaslab.info.g_Log import Log


def _process_chromosome(args):
    """
    Process a single chromosome. This function is designed to be called in parallel.
    
    Parameters:
    -----------
    args : tuple
        Tuple containing all arguments needed to process a chromosome:
        (chrom, ref_dir, bim_prefix, sst_file, param_dict, mcmc2, debug, seed_offset)
    
    Returns:
    --------
    chrom : int
        Chromosome number (for result tracking)
    """
    chrom, ref_dir, bim_prefix, sst_file, param_dict, mcmc2, debug, seed_offset = args
    
    # Create a simple logger for this process (log object may not be pickleable)
    # We'll write to a temporary log or suppress logging in parallel mode
    try:
        from gwaslab.info.g_Log import Log
        log = Log(verbose=False)  # Create a minimal log for this process
    except:
        class DummyLog:
            def write(self, *args, **kwargs):
                pass
        log = DummyLog()
    
    log.write('##### process chromosome %d #####' % int(chrom))

    # Parse reference panel
    if '1kg' in os.path.basename(ref_dir):
        ref_dict = parse_genet.parse_ref(ref_dir + '/snpinfo_1kg_hm3', int(chrom), log)
    elif 'ukbb' in os.path.basename(ref_dir):
        ref_dict = parse_genet.parse_ref(ref_dir + '/snpinfo_ukbb_hm3', int(chrom), log)
    else:
        log.write('Warning: Unknown reference panel type')
        return chrom

    vld_dict = parse_genet.parse_bim(bim_prefix, int(chrom))

    sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, sst_file, param_dict['n_gwas'], log)

    ld_blk, blk_size = parse_genet.parse_ldblk(ref_dir, sst_dict, int(chrom), log)

    # Adjust seed for this chromosome to ensure reproducibility
    # Each chromosome gets a different seed offset
    chrom_seed = None
    if param_dict['seed'] is not None:
        chrom_seed = param_dict['seed'] + seed_offset

    # Use optimized mcmc2 if requested, otherwise use original mcmc
    if mcmc2:
        mcmc_gtb2.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, 
                      param_dict['n_gwas'], ld_blk, blk_size,
                      param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], 
                      int(chrom), param_dict['out_dir'], param_dict['beta_std'],
                      param_dict['write_psi'], param_dict['write_pst'], chrom_seed, log, debug)
    else:
        mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, 
                     param_dict['n_gwas'], ld_blk, blk_size,
                     param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], 
                     int(chrom), param_dict['out_dir'], param_dict['beta_std'],
                     param_dict['write_psi'], param_dict['write_pst'], chrom_seed, log)
    
    return chrom


def _run_prscs(
         ref_dir: Optional[str] = None,
         bim_prefix: Optional[str] = None,
         sst_file: Optional[pd.DataFrame] = None,
            a: float = 1, 
            b: float = 0.5, 
            phi: Optional[float] = None, 
            n_gwas: Optional[int] = None,
            n_iter: int = 1000, 
            n_burnin: int = 500, 
            thin: int = 5, 
            out_dir: str = "./", 
            chrom: Union[range, list] = range(1,23),
            beta_std: str = 'FALSE', 
            write_psi: str = 'FALSE', 
            write_pst: str = 'FALSE', 
            seed: Optional[int] = None,
            log: Optional['Log'] = None,
            mcmc2=False,
            threads: int = 1,
         **kwargs: Any) -> None:
    ## T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.Nature Communications, 10:1776, 2019.
    sst_file = sst_file.rename(columns={"rsID":"SNP","POS":"BP"})
    log.write("Start to runnig PRScs...")
    
    # Prepare parameters
    param_dict = {'ref_dir': ref_dir, 
                    'bim_prefix': bim_prefix, 
                    'a': a, 
                    'b': b, 
                    'phi': phi, 
                    'n_gwas': n_gwas,
                    'n_iter': n_iter, 
                    'n_burnin': n_burnin, 
                    'thin': thin, 
                    'out_dir': out_dir, 
                    'beta_std': beta_std, 
                    'write_psi': write_psi, 
                    'write_pst': write_pst, 
                    'seed': seed}
    
    # Get debug flag
    debug = kwargs.get('debug', False)
    
    # Determine number of threads
    chrom_list = list(chrom) if isinstance(chrom, range) else chrom
    n_chrom = len(chrom_list)
    
    # Limit threads to number of chromosomes and available CPUs
    max_threads = min(threads, n_chrom, cpu_count())
    
    if max_threads > 1:
        log.write("Running PRS-CS in parallel mode with {} threads for {} chromosomes".format(max_threads, n_chrom))
        
        # Prepare arguments for each chromosome
        # Each chromosome gets a unique seed offset to ensure reproducibility
        args_list = []
        for idx, chrom in enumerate(chrom_list):
            args_list.append((chrom, ref_dir, bim_prefix, sst_file, param_dict, 
                            mcmc2, debug, idx * 1000))  # Offset seed by 1000 per chromosome
        
        # Process chromosomes in parallel
        with Pool(processes=max_threads) as pool:
            results = pool.map(_process_chromosome, args_list)
        
        log.write("Completed processing {} chromosomes in parallel".format(len(results)))
    else:
        # Sequential processing (original behavior)
        log.write("Running PRS-CS in sequential mode")
        for idx, chrom in enumerate(chrom_list):
            log.write('##### process chromosome %d #####' % int(chrom))

            if '1kg' in os.path.basename(ref_dir):
                ref_dict = parse_genet.parse_ref(ref_dir + '/snpinfo_1kg_hm3', int(chrom), log)
            elif 'ukbb' in os.path.basename(ref_dir):
                ref_dict = parse_genet.parse_ref(ref_dir + '/snpinfo_ukbb_hm3', int(chrom), log)

            vld_dict = parse_genet.parse_bim(bim_prefix, int(chrom))

            sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, sst_file, param_dict['n_gwas'], log)

            ld_blk, blk_size = parse_genet.parse_ldblk(ref_dir, sst_dict, int(chrom), log)

            # Adjust seed for this chromosome
            chrom_seed = None
            if seed is not None:
                chrom_seed = seed + idx * 1000

            # Use optimized mcmc2 if requested, otherwise use original mcmc
            if mcmc2:
                mcmc_gtb2.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, 
                              param_dict['n_gwas'], ld_blk, blk_size,
                              param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], 
                              int(chrom), param_dict['out_dir'], param_dict['beta_std'],
                              param_dict['write_psi'], param_dict['write_pst'], chrom_seed, log, debug)
            else:
                mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, 
                             param_dict['n_gwas'], ld_blk, blk_size,
                             param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], 
                             int(chrom), param_dict['out_dir'], param_dict['beta_std'],
                             param_dict['write_psi'], param_dict['write_pst'], chrom_seed, log)

    log.write("Finished!")

