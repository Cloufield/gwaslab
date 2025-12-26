"""
Test suite for get_proxy function to find LD proxy variants.

Tests cover:
- Basic proxy finding with VCF reference
- Finding proxies with include_all option
- Handling SNPs not found in VCF
- Edge cases with empty results
- Multiple SNPs in snplist
"""

import os
import sys
import unittest
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.g_Sumstats import Sumstats

# Test data paths
TEST_DIR = os.path.join(os.path.dirname(__file__), "ref")
SUMMSTATS_PATH = os.path.join(TEST_DIR, "bbj_t2d_hm3_chr7_variants.txt.gz")
VCF_PATH = os.path.join(TEST_DIR, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")


class TestGetProxyBasic(unittest.TestCase):
    """Test basic get_proxy functionality."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests."""
        # Check if test data files exist
        if not os.path.exists(SUMMSTATS_PATH):
            raise unittest.SkipTest(f"Test data file not found: {SUMMSTATS_PATH}")
        if not os.path.exists(VCF_PATH):
            raise unittest.SkipTest(f"Test data file not found: {VCF_PATH}")
        
        # Load sumstats
        cls.sumstats = Sumstats(
            sumstats=SUMMSTATS_PATH,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            build="19",
            verbose=False
        )
        cls.sumstats.basic_check(verbose=False)
        
        # VCF region: chr7:126253550-128253550
        # Filter sumstats to SNPs in VCF region for more reliable tests
        cls.vcf_start = 126253550
        cls.vcf_end = 128253550
        cls.sumstats_in_region = cls.sumstats.data[
            (cls.sumstats.data["CHR"] == 7) &
            (cls.sumstats.data["POS"] >= cls.vcf_start) &
            (cls.sumstats.data["POS"] <= cls.vcf_end)
        ].copy()
    
    def test_get_proxy_basic(self):
        """Test basic get_proxy functionality with a single SNP."""
        # Get a SNP from the sumstats in VCF region
        if len(self.sumstats_in_region) == 0:
            self.skipTest("No SNPs available in VCF region")
        
        test_snp = self.sumstats_in_region["SNPID"].iloc[0]
        result = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.8,
            verbose=False
        )
        
        # Check result structure
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        
        # If proxies found, check columns
        if len(result) > 0:
            expected_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "RSQ", "LD_REF_VARIANT"]
            for col in expected_cols:
                self.assertIn(col, result.columns, f"Column {col} missing from result")
            
            # Check that RSQ values are above threshold
            if "RSQ" in result.columns:
                self.assertTrue(
                    (result["RSQ"] > 0.8).all(),
                    "All RSQ values should be above threshold"
                )
    
    def test_get_proxy_multiple_snps(self):
        """Test get_proxy with multiple SNPs in snplist."""
        # Get multiple SNPs from the sumstats in VCF region
        if len(self.sumstats_in_region) < 2:
            self.skipTest("Not enough SNPs available in VCF region")
        
        test_snps = self.sumstats_in_region["SNPID"].head(3).tolist()
        
        result = self.sumstats.get_proxy(
            snplist=test_snps,
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.8,
            verbose=False
        )
        
        # Check result structure
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        
        # If proxies found, check that LD_REF_VARIANT references the input SNPs
        if len(result) > 0 and "LD_REF_VARIANT" in result.columns:
            ref_variants = result["LD_REF_VARIANT"].unique()
            for ref_var in ref_variants:
                self.assertIn(ref_var, test_snps, f"LD_REF_VARIANT {ref_var} not in input snplist")
    
    def test_get_proxy_include_all(self):
        """Test get_proxy with include_all=True option."""
        # Get a SNP from the sumstats in VCF region
        if len(self.sumstats_in_region) == 0:
            self.skipTest("No SNPs available in VCF region")
        
        test_snp = self.sumstats_in_region["SNPID"].iloc[0]
        
        # Test with include_all=False (default)
        result_default = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.8,
            include_all=False,
            verbose=False
        )
        
        # Test with include_all=True
        result_all = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.8,
            include_all=True,
            verbose=False
        )
        
        # Check that include_all may return more or equal proxies
        self.assertIsNotNone(result_default)
        self.assertIsNotNone(result_all)
        self.assertIsInstance(result_default, pd.DataFrame)
        self.assertIsInstance(result_all, pd.DataFrame)
        
        # include_all should return at least as many proxies as default
        # (it includes variants not in sumstats)
        if len(result_all) > 0:
            self.assertGreaterEqual(
                len(result_all),
                len(result_default),
                "include_all=True should return at least as many proxies as include_all=False"
            )
    
    def test_get_proxy_empty_snplist(self):
        """Test get_proxy with empty snplist."""
        result = self.sumstats.get_proxy(
            snplist=[],
            vcf_path=VCF_PATH,
            verbose=False
        )
        
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 0)
    
    def test_get_proxy_snp_not_in_sumstats(self):
        """Test get_proxy with SNP not in sumstats."""
        # Use a SNP ID that's unlikely to be in the sumstats
        fake_snp = "7:999999999_A_G"
        
        result = self.sumstats.get_proxy(
            snplist=[fake_snp],
            vcf_path=VCF_PATH,
            verbose=False
        )
        
        # Should return empty DataFrame when SNP not found
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 0)
    
    def test_get_proxy_different_ld_threshold(self):
        """Test get_proxy with different LD thresholds."""
        if len(self.sumstats_in_region) == 0:
            self.skipTest("No SNPs available in VCF region")
        
        test_snp = self.sumstats_in_region["SNPID"].iloc[0]
        
        # Test with higher threshold (should return fewer proxies)
        result_high = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.9,
            verbose=False
        )
        
        # Test with lower threshold (should return more proxies)
        result_low = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.7,
            verbose=False
        )
        
        self.assertIsNotNone(result_high)
        self.assertIsNotNone(result_low)
        
        # Lower threshold should return at least as many proxies
        if len(result_low) > 0:
            self.assertGreaterEqual(
                len(result_low),
                len(result_high),
                "Lower LD threshold should return at least as many proxies"
            )
    
    def test_get_proxy_different_window_size(self):
        """Test get_proxy with different window sizes."""
        if len(self.sumstats_in_region) == 0:
            self.skipTest("No SNPs available in VCF region")
        
        test_snp = self.sumstats_in_region["SNPID"].iloc[0]
        
        # Test with smaller window
        result_small = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=100,
            ld_threshold=0.8,
            verbose=False
        )
        
        # Test with larger window
        result_large = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=1000,
            ld_threshold=0.8,
            verbose=False
        )
        
        self.assertIsNotNone(result_small)
        self.assertIsNotNone(result_large)
        
        # Larger window should return at least as many proxies
        if len(result_large) > 0:
            self.assertGreaterEqual(
                len(result_large),
                len(result_small),
                "Larger window should return at least as many proxies"
            )
    
    def test_get_proxy_result_sorted_by_rsq(self):
        """Test that get_proxy results are sorted by RSQ in descending order."""
        if len(self.sumstats_in_region) == 0:
            self.skipTest("No SNPs available in VCF region")
        
        test_snp = self.sumstats_in_region["SNPID"].iloc[0]
        result = self.sumstats.get_proxy(
            snplist=[test_snp],
            vcf_path=VCF_PATH,
            windowsizekb=500,
            ld_threshold=0.8,
            verbose=False
        )
        
        if len(result) > 1 and "RSQ" in result.columns:
            # Check that RSQ values are sorted in descending order
            rsq_values = result["RSQ"].values
            is_sorted = all(rsq_values[i] >= rsq_values[i+1] for i in range(len(rsq_values)-1))
            self.assertTrue(
                is_sorted,
                "RSQ values should be sorted in descending order"
            )


if __name__ == "__main__":
    # Run tests
    unittest.main(verbosity=2)

