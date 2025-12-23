import os
import sys
import unittest

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.g_Sumstats import Sumstats


class TestLDSC(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.test_data_path = os.path.join(ROOT, "test", "ref", "bbj_t2d_hm3_chr7_variants.txt.gz")
        # LDSC expects path with @ as chromosome placeholder (e.g., /path/to/@ becomes /path/to/1, /path/to/2, etc.)
        self.ldscore_path = os.path.join(ROOT, "test", "ref", "eas_ldscores", "@")
        
        # Load test sumstats
        self.sumstats = Sumstats(
            sumstats=self.test_data_path,
            fmt="gwaslab",
            verbose=False
        )
        
        # Set build to match LD scores (assuming hg19/GRCh37)
        self.sumstats.set_build("19", verbose=False)

    def test_estimate_h2_by_ldsc_basic(self):
        """Test basic heritability estimation using LDSC"""
        # Run LDSC heritability estimation
        result = self.sumstats.estimate_h2_by_ldsc(
            ref_ld_chr=self.ldscore_path,
            w_ld_chr=self.ldscore_path,
            verbose=False
        )
        
        # Check that results are stored
        import pandas as pd
        self.assertIsNotNone(self.sumstats.ldsc_h2, "ldsc_h2 should be set")
        self.assertIsInstance(self.sumstats.ldsc_h2, pd.DataFrame, "ldsc_h2 should be a DataFrame")
        
        # Check that result is a DataFrame (or None if coefficients not requested)
        if result is not None:
            self.assertIsInstance(result, pd.DataFrame, "Result should be a DataFrame")
            if len(result) > 0:
                self.assertGreater(len(result), 0, "Result DataFrame should have rows")
        
        # Check that ldsc_h2 has content
        self.assertGreater(len(self.sumstats.ldsc_h2), 0, "ldsc_h2 DataFrame should have rows")

    def test_estimate_h2_by_ldsc_with_coefficients(self):
        """Test heritability estimation with coefficient results"""
        # Run LDSC with print_coefficients enabled
        result = self.sumstats.estimate_h2_by_ldsc(
            ref_ld_chr=self.ldscore_path,
            w_ld_chr=self.ldscore_path,
            print_coefficients="ldsc",
            verbose=False
        )
        
        # Check that results are stored
        self.assertIsNotNone(self.sumstats.ldsc_h2, "ldsc_h2 should be set")
        
        # Check coefficient results if available
        if self.sumstats.ldsc_h2_results is not None:
            import pandas as pd
            self.assertIsInstance(self.sumstats.ldsc_h2_results, pd.DataFrame, 
                                "ldsc_h2_results should be a DataFrame")

    def test_estimate_h2_by_ldsc_with_munge(self):
        """Test heritability estimation with munging"""
        # Create a fresh sumstats object for this test
        sumstats = Sumstats(
            sumstats=self.test_data_path,
            fmt="gwaslab",
            verbose=False
        )
        sumstats.set_build("19", verbose=False)
        
        # Run LDSC with munging
        result = sumstats.estimate_h2_by_ldsc(
            ref_ld_chr=self.ldscore_path,
            w_ld_chr=self.ldscore_path,
            munge=True,
            munge_kwargs={"info": 0.9, "maf": 0.01},
            verbose=False
        )
        
        # Check that results are stored
        self.assertIsNotNone(sumstats.ldsc_h2, "ldsc_h2 should be set after munging")

    def test_estimate_h2_by_ldsc_build_parameter(self):
        """Test that build parameter works correctly"""
        sumstats = Sumstats(
            sumstats=self.test_data_path,
            fmt="gwaslab",
            verbose=False
        )
        
        # Test with explicit build parameter
        result = sumstats.estimate_h2_by_ldsc(
            ref_ld_chr=self.ldscore_path,
            w_ld_chr=self.ldscore_path,
            build="19",
            verbose=False
        )
        
        # Check that results are stored
        self.assertIsNotNone(sumstats.ldsc_h2, "ldsc_h2 should be set")

    def test_estimate_h2_by_ldsc_verbose_false(self):
        """Test that verbose=False works without errors"""
        sumstats = Sumstats(
            sumstats=self.test_data_path,
            fmt="gwaslab",
            verbose=False
        )
        sumstats.set_build("19", verbose=False)
        
        # Should not raise any errors
        try:
            result = sumstats.estimate_h2_by_ldsc(
                ref_ld_chr=self.ldscore_path,
                w_ld_chr=self.ldscore_path,
                verbose=False
            )
            success = True
        except Exception as e:
            success = False
            self.fail(f"estimate_h2_by_ldsc raised an exception: {e}")
        
        self.assertTrue(success, "estimate_h2_by_ldsc should complete without errors")


if __name__ == "__main__":
    unittest.main()

