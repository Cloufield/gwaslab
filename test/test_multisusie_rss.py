#!/usr/bin/env python3
"""
Test script for MultiSuSiE RSS integration with GWASLab.

This script tests the run_multisusie_rss method using sample data from
the MultiSuSiE repository.
"""

import os
import sys
import unittest
import numpy as np
import pandas as pd
import scipy.stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.g_Sumstats import Sumstats
from gwaslab.g_SumstatsMulti import SumstatsMulti

# Path to MultiSuSiE example data - can be set via environment variable
MULTISUSIE_DATA_DIR = os.environ.get(
    'MULTISUSIE_DATA_DIR',
    '/home/yunye/work/github/MultiSuSiE/example_data'
)


class TestMultiSuSiERSS(unittest.TestCase):
    """Test MultiSuSiE RSS integration with GWASLab."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data once for all tests."""
        # Check if data directory exists
        if not os.path.exists(MULTISUSIE_DATA_DIR):
            raise unittest.SkipTest(
                f"MultiSuSiE data directory not found: {MULTISUSIE_DATA_DIR}. "
                "Set MULTISUSIE_DATA_DIR environment variable to point to the data directory."
            )
        
        # Load data
        cls.geno_list = cls.load_genotype_data()
        cls.beta_hat_list, cls.se_list, cls.n_list, cls.varY_list = cls.create_summary_statistics(cls.geno_list)
        cls.R_list = cls.create_ld_matrices(cls.geno_list)
        cls.snp_meta = cls.load_snp_metadata()
        cls.sumstats_objects = cls.create_sumstats_objects(
            cls.beta_hat_list, cls.se_list, cls.n_list, cls.snp_meta
        )
        cls.sumstats_multi = SumstatsMulti(
            sumstatsObjects=cls.sumstats_objects,
            group_name="TestGroup",
            verbose=False
        )
    
    @staticmethod
    def load_genotype_data():
        """Load genotype data from MultiSuSiE example files."""
        geno_YRI = np.loadtxt(os.path.join(MULTISUSIE_DATA_DIR, 'geno_YRI.txt'))
        geno_CEU = np.loadtxt(os.path.join(MULTISUSIE_DATA_DIR, 'geno_CEU.txt'))
        geno_JPT = np.loadtxt(os.path.join(MULTISUSIE_DATA_DIR, 'geno_JPT.txt'))
        return [geno_YRI, geno_CEU, geno_JPT]
    
    @staticmethod
    def create_summary_statistics(geno_list):
        """Create summary statistics from genotype data."""
        # Create simulated effect sizes
        beta_YRI = np.zeros(40)
        beta_CEU = np.zeros(40)
        beta_JPT = np.zeros(40)
        beta_YRI[10] = 0.75
        beta_CEU[10] = 1.0
        beta_JPT[10] = 0.5
        beta_YRI[3] = 0.5
        beta_CEU[3] = 0.5
        beta_JPT[3] = 0.5
        beta_YRI[38] = 1.0
        beta_CEU[38] = 0.0
        beta_JPT[38] = 0.0
        beta_list = [beta_YRI, beta_CEU, beta_JPT]
        
        # Simulate phenotypes
        rng = np.random.default_rng(2)
        y_list = [geno.dot(beta) + rng.standard_normal(geno.shape[0]) 
                  for (geno, beta) in zip(geno_list, beta_list)]
        y_list = [y - np.mean(y) for y in y_list]
        
        # Calculate summary statistics
        XTY_list = [geno.T.dot(y) for (geno, y) in zip(geno_list, y_list)]
        XTX_diagonal_list = [np.diagonal(geno.T.dot(geno)) for geno in geno_list]
        
        with np.errstate(divide='ignore', invalid='ignore'):
            beta_hat_list = [XTY / XTX_diag 
                            for (XTY, XTX_diag) in zip(XTY_list, XTX_diagonal_list)]
        
        N_list = [geno.shape[0] for geno in geno_list]
        # Calculate residuals using beta_hat (not true beta)
        residuals_list = [np.expand_dims(y, 1) - (geno * beta_hat) 
                          for (y, geno, beta_hat) in zip(y_list, geno_list, beta_hat_list)]
        sum_of_squared_residuals_list = [np.sum(resid ** 2, axis=0) 
                                         for resid in residuals_list]
        se_list = [np.sqrt(ssr / ((N - 2) * XTX)) 
                   for (ssr, N, XTX) in zip(sum_of_squared_residuals_list, N_list, XTX_diagonal_list)]
        
        varY_list = [np.var(y, ddof=1) for y in y_list]
        
        return beta_hat_list, se_list, N_list, varY_list
    
    @staticmethod
    def create_ld_matrices(geno_list):
        """Create LD correlation matrices from genotype data."""
        with np.errstate(divide='ignore', invalid='ignore'):
            R_list = [np.corrcoef(geno, rowvar=False) for geno in geno_list]
        return R_list
    
    @staticmethod
    def load_snp_metadata():
        """Load SNP metadata."""
        # Format: rsid,CHR,0.0,POS,EA,NEA
        snp_meta = pd.read_csv(
            os.path.join(MULTISUSIE_DATA_DIR, 'snp_meta.csv'),
            header=None,
            names=['SNPID', 'CHR', 'DUMMY', 'POS', 'EA', 'NEA']
        )
        # Drop the dummy column (0.0)
        snp_meta = snp_meta.drop(columns=['DUMMY'])
        return snp_meta
    
    @staticmethod
    def create_sumstats_objects(beta_hat_list, se_list, n_list, snp_meta):
        """Create Sumstats objects for each population."""
        sumstats_objects = []
        study_names = ['YRI', 'CEU', 'JPT']
        
        for i, (beta_hat, se, n, study_name) in enumerate(zip(beta_hat_list, se_list, n_list, study_names)):
            # Calculate p-values from z-scores
            z_scores = beta_hat / se
            p_values = 2 * (1 - np.abs(scipy.stats.norm.cdf(np.abs(z_scores))))
            
            df = pd.DataFrame({
                'SNPID': snp_meta['SNPID'].values,
                'CHR': snp_meta['CHR'].values,
                'POS': snp_meta['POS'].values,
                'EA': snp_meta['EA'].values,
                'NEA': snp_meta['NEA'].values,
                'BETA': beta_hat,
                'SE': se,
                'N': n,
                'P': p_values
            })
            
            # Create Sumstats object
            sumstats = Sumstats(
                sumstats=df,
                study=study_name,
                verbose=False
            )
            
            sumstats_objects.append(sumstats)
        
        return sumstats_objects
    
    def test_sumstats_multi_creation(self):
        """Test that SumstatsMulti object is created correctly."""
        self.assertIsNotNone(self.sumstats_multi)
        self.assertEqual(self.sumstats_multi.meta['gwaslab']['number_of_studies'], 3)
        self.assertGreater(len(self.sumstats_multi.data), 0)
        
        # Check that expected columns exist
        expected_cols = ['SNPID', 'BETA_1', 'BETA_2', 'BETA_3', 'SE_1', 'SE_2', 'SE_3']
        for col in expected_cols:
            self.assertIn(col, self.sumstats_multi.data.columns, f"Missing column: {col}")
    
    def test_run_multisusie_rss(self):
        """Test running MultiSuSiE RSS finemapping."""
        # Run MultiSuSiE RSS
        results = self.sumstats_multi.run_multisusie_rss(
            R_list=self.R_list,
            varY_list=self.varY_list,
            rho=0.75,
            L=10,
            scaled_prior_variance=0.2,
            max_iter=100,
            verbose=False,
            float_type=np.float64,
            single_population_mac_thresh=10,
            low_memory_mode=False,
            recover_R=False
        )
        
        # Assert results are not None
        self.assertIsNotNone(results, "Results should not be None")
        
        # Assert results is a DataFrame
        self.assertIsInstance(results, pd.DataFrame, "Results should be a DataFrame")
        
        # Assert results has expected columns
        expected_cols = ['SNPID', 'PIP', 'CREDIBLE_SET_INDEX']
        for col in expected_cols:
            self.assertIn(col, results.columns, f"Results missing column: {col}")
        
        # Assert results has the expected number of variants
        self.assertEqual(len(results), len(self.snp_meta), 
                        f"Results should have {len(self.snp_meta)} variants")
        
        # Assert PIP values are in valid range [0, 1]
        self.assertTrue(
            results['PIP'].between(0, 1, inclusive='both').all(),
            "All PIP values should be between 0 and 1"
        )
        
        # Assert CREDIBLE_SET_INDEX is non-negative
        self.assertTrue(
            (results['CREDIBLE_SET_INDEX'] >= 0).all(),
            "All CREDIBLE_SET_INDEX values should be non-negative"
        )
        
        # Check PIPs for causal variants
        # Expected from example.ipynb Cell 15 (with mac_thresh=10): 
        # 0.4655837427342536, 0.9999781922032149, 0.969634716536287
        causal_snpids = ['rs28475450', 'rs2765023', 'rs9439458']  # variants 3, 10, 38
        expected_pips = [0.465584, 0.999978, 0.969635]
        tolerance = 0.01  # Allow 1% tolerance for numerical differences
        
        for snpid, expected_pip in zip(causal_snpids, expected_pips):
            variant_result = results[results['SNPID'] == snpid]
            self.assertGreater(
                len(variant_result), 0,
                f"Causal variant {snpid} should be in results"
            )
            pip = variant_result.iloc[0]['PIP']
            self.assertAlmostEqual(
                pip, expected_pip, places=3,
                msg=f"PIP for {snpid} should be approximately {expected_pip}, got {pip}"
            )
        
        # Assert that at least one credible set was found
        num_credible_sets = len(results[results['CREDIBLE_SET_INDEX'] > 0].groupby('CREDIBLE_SET_INDEX'))
        self.assertGreater(
            num_credible_sets, 0,
            "At least one credible set should be found"
        )
        
        # Assert that the top variant by PIP has high PIP (> 0.9)
        max_pip = results['PIP'].max()
        self.assertGreater(
            max_pip, 0.9,
            f"Top variant should have PIP > 0.9, got {max_pip}"
        )


if __name__ == "__main__":
    unittest.main()
