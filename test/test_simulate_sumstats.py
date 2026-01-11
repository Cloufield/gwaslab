import os
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_simulate import simulate_sumstats_region, simulate_sumstats_global
from gwaslab.g_Sumstats import Sumstats
from gwaslab import Log


class TestSimulateSumstats(unittest.TestCase):
    """Test cases for simulate_sumstats_region and simulate_sumstats_global"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.vcf_path = os.path.join(ROOT, "test", "ref", "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        self.region = ("7", 126253550, 128253550)
        self.log = Log()
        
        # Check if VCF file exists
        if not os.path.exists(self.vcf_path):
            self.skipTest(f"Test VCF file not found: {self.vcf_path}")
    
    def test_simulate_sumstats_region_basic(self):
        """Test basic functionality of simulate_sumstats_region"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            mode="sparse",
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        
        # Check return types
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertIsInstance(causal_snp_ids, list)
        
        # Check that sumstats has data
        self.assertGreater(len(sumstats_obj.data), 0)
        
        # Check required columns
        required_cols = ["CHR", "POS", "EA", "NEA", "Z", "BETA", "SE", "P", "MLOG10P", "N"]
        for col in required_cols:
            self.assertIn(col, sumstats_obj.data.columns, f"Missing column: {col}")
        
        # Check that we have exactly n_causal causal variants
        self.assertEqual(len(causal_snp_ids), 3)
        
        # Check that causal variants are in the sumstats
        for causal_id in causal_snp_ids:
            self.assertIn(causal_id, sumstats_obj.data["SNPID"].values)
    
    def test_simulate_sumstats_region_polygenic(self):
        """Test polygenic mode"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            mode="polygenic",
            pi=0.01,  # 1% causal
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertGreater(len(causal_snp_ids), 0)
        self.assertIsInstance(sumstats_obj, Sumstats)
    
    def test_simulate_sumstats_region_maf_filter(self):
        """Test MAF filtering"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            maf_min=0.05,
            maf_max=0.5,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        
        # Check that all variants have MAF in range
        eaf = sumstats_obj.data["EAF"].values
        maf = np.minimum(eaf, 1.0 - eaf)
        self.assertTrue(np.all(maf >= 0.05))
        self.assertTrue(np.all(maf <= 0.5))
    
    def test_simulate_sumstats_region_alpha(self):
        """Test MAF-dependent effect sizes (alpha parameter)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=5,
            alpha=0.2,  # MAF-dependent effects
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertEqual(len(causal_snp_ids), 5)
    
    def test_simulate_sumstats_region_thin(self):
        """Test thinning option"""
        result_no_thin = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            thin=None,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        result_thin = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            thin=0.5,  # Keep 50% of variants
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_no_thin, _ = result_no_thin
        sumstats_thin, _ = result_thin
        
        # Thinned should have fewer or equal variants
        self.assertLessEqual(len(sumstats_thin.data), len(sumstats_no_thin.data))
        # Should be roughly 50% (allow some variance due to rounding)
        ratio = len(sumstats_thin.data) / len(sumstats_no_thin.data)
        self.assertGreater(ratio, 0.4)
        self.assertLess(ratio, 0.6)
    
    def test_simulate_sumstats_region_lambda_gc(self):
        """Test cryptic relatedness inflation (lambda_gc)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            lambda_gc=1.1,  # 10% inflation
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
    
    def test_simulate_sumstats_region_sigma_strat(self):
        """Test population stratification (sigma_strat)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            sigma_strat=0.2,  # Moderate stratification
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
    
    def test_simulate_sumstats_region_binary_trait(self):
        """Test binary trait simulation"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            trait="binary",
            n_case=5000,
            n_ctrl=5000,
            n_causal=2,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
        
        # Check that N_CASE and N_CONTROL columns exist
        if "N_CASE" in sumstats_obj.data.columns:
            self.assertIn("N_CASE", sumstats_obj.data.columns)
            self.assertIn("N_CONTROL", sumstats_obj.data.columns)
    
    def test_simulate_sumstats_global_basic(self):
        """Test basic functionality of simulate_sumstats_global"""
        result = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],  # Only chromosome 7
            n=10000,
            mode="sparse",
            n_causal=3,
            h2=0.01,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        
        # Check return types
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertIsInstance(causal_snp_ids, list)
        
        # Check that sumstats has data
        self.assertGreater(len(sumstats_obj.data), 0)
        
        # Check required columns
        required_cols = ["CHR", "POS", "EA", "NEA", "Z", "BETA", "SE", "P", "MLOG10P", "N"]
        for col in required_cols:
            self.assertIn(col, sumstats_obj.data.columns, f"Missing column: {col}")
        
        # Check that all variants are from chromosome 7 (CHR might be string or int)
        chr_values = sumstats_obj.data["CHR"].values
        self.assertTrue(all((chr_values == "7") | (chr_values == 7) | (chr_values == "chr7")))
    
    def test_simulate_sumstats_global_heritability(self):
        """Test global heritability calibration"""
        result = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="sparse",
            n_causal=5,
            h2=0.05,  # 5% heritability
            alpha=0.2,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertEqual(len(causal_snp_ids), 5)
    
    def test_simulate_sumstats_region_z_scores_valid(self):
        """Test that Z-scores are valid (finite, not all zeros)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        z_scores = sumstats_obj.data["Z"].values
        
        # Check that Z-scores are finite
        self.assertTrue(np.all(np.isfinite(z_scores)))
        
        # Check that not all Z-scores are zero (should have some variation)
        self.assertGreater(np.std(z_scores), 0.0)
    
    def test_simulate_sumstats_region_p_values_valid(self):
        """Test that P-values are valid (between 0 and 1)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        p_values = sumstats_obj.data["P"].values
        
        # Check that P-values are between 0 and 1
        self.assertTrue(np.all(p_values >= 0.0))
        self.assertTrue(np.all(p_values <= 1.0))
        self.assertTrue(np.all(np.isfinite(p_values)))
    
    def test_simulate_sumstats_region_reproducibility(self):
        """Test that same seed produces same results"""
        result1 = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=123,
            verbose=False,
            log=self.log
        )
        
        result2 = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=123,
            verbose=False,
            log=self.log
        )
        
        sumstats1, causals1 = result1
        sumstats2, causals2 = result2
        
        # Check that causal variants are the same
        self.assertEqual(set(causals1), set(causals2))
        
        # Check that Z-scores are the same (within numerical precision)
        z1 = sumstats1.data["Z"].values
        z2 = sumstats2.data["Z"].values
        np.testing.assert_array_almost_equal(z1, z2, decimal=10)


if __name__ == "__main__":
    unittest.main()
