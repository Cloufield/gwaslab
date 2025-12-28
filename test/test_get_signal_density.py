"""
Test suite for _get_signal_density2 function with comprehensive test cases.

Tests cover:
- Basic density calculation (without sig_sumstats)
- Conditional density calculation (with sig_sumstats)
- Edge cases (empty data, single variant, multiple chromosomes)
- Different window sizes
- Missing values handling
- Boundary conditions
"""

import os
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_get_density import _get_signal_density2
from gwaslab.info.g_Log import Log


class TestGetSignalDensity2(unittest.TestCase):
    """Test cases for _get_signal_density2 function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.log.verbose = False
    
    def test_basic_density_single_chromosome(self):
        """Test basic density calculation on a single chromosome."""
        # Create variants on chromosome 1, spaced 50kb apart
        df = pd.DataFrame({
            "CHR": [1] * 5,
            "POS": [100000, 150000, 200000, 250000, 300000],
            "SNPID": [f"1:{pos}_A_G" for pos in [100000, 150000, 200000, 250000, 300000]],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Check that DENSITY column was added
        self.assertIn("DENSITY", result.columns)
        
        # With 100kb window (symmetric), each variant sees variants in [pos-wsize, pos+wsize]
        # Variant at 100k: sees variants in [0, 200k] = 100k, 150k, 200k = 3 total - 1 self = 2
        # Variant at 150k: sees variants in [50k, 250k] = 100k, 150k, 200k, 250k = 4 total - 1 self = 3
        # Variant at 200k: sees variants in [100k, 300k] = 100k, 150k, 200k, 250k, 300k = 5 total - 1 self = 4
        # Variant at 250k: sees variants in [150k, 350k] = 150k, 200k, 250k, 300k = 4 total - 1 self = 3
        # Variant at 300k: sees variants in [200k, 400k] = 200k, 250k, 300k = 3 total - 1 self = 2
        
        densities = result["DENSITY"].values
        self.assertEqual(densities[0], 2)  # 100k
        self.assertEqual(densities[1], 3)  # 150k
        self.assertEqual(densities[2], 4)  # 200k
        self.assertEqual(densities[3], 3)  # 250k
        self.assertEqual(densities[4], 2)  # 300k
    
    def test_basic_density_multiple_chromosomes(self):
        """Test basic density calculation across multiple chromosomes."""
        df = pd.DataFrame({
            "CHR": [1, 1, 2, 2, 2],
            "POS": [100000, 150000, 100000, 150000, 200000],
            "SNPID": ["1:100000_A_G", "1:150000_A_G", "2:100000_A_G", "2:150000_A_G", "2:200000_A_G"],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Check that DENSITY column was added
        self.assertIn("DENSITY", result.columns)
        
        # Chromosome 1: variants at 100k and 150k (50kb apart, within 100kb window)
        # Variant at 100k: sees variants in [0, 200k] = 100k, 150k = 2 total - 1 self = 1
        # Variant at 150k: sees variants in [50k, 250k] = 100k, 150k = 2 total - 1 self = 1
        chr1_mask = result["CHR"] == 1
        chr1_densities = result.loc[chr1_mask, "DENSITY"].values
        self.assertEqual(chr1_densities[0], 1)  # 100k sees 150k
        self.assertEqual(chr1_densities[1], 1)  # 150k sees 100k
        
        # Chromosome 2: variants at 100k, 150k, 200k
        # Variant at 100k: sees variants in [0, 200k] = 100k, 150k, 200k = 3 total - 1 self = 2
        # Variant at 150k: sees variants in [50k, 250k] = 100k, 150k, 200k = 3 total - 1 self = 2
        # Variant at 200k: sees variants in [100k, 300k] = 100k, 150k, 200k = 3 total - 1 self = 2
        chr2_mask = result["CHR"] == 2
        chr2_densities = result.loc[chr2_mask, "DENSITY"].values
        self.assertEqual(chr2_densities[0], 2)  # 100k sees 150k and 200k
        self.assertEqual(chr2_densities[1], 2)  # 150k sees 100k and 200k
        self.assertEqual(chr2_densities[2], 2)  # 200k sees 100k and 150k
    
    def test_basic_density_large_window(self):
        """Test density calculation with a large window size."""
        df = pd.DataFrame({
            "CHR": [1] * 5,
            "POS": [100000, 200000, 300000, 400000, 500000],
            "SNPID": [f"1:{pos}_A_G" for pos in [100000, 200000, 300000, 400000, 500000]],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=500,  # 500kb window
            log=self.log,
            verbose=False
        )
        
        # With 500kb window, all variants should see all others (5 variants - 1 self = 4)
        densities = result["DENSITY"].values
        self.assertTrue(all(d == 4 for d in densities))
    
    def test_basic_density_small_window(self):
        """Test density calculation with a small window size."""
        df = pd.DataFrame({
            "CHR": [1] * 5,
            "POS": [100000, 200000, 300000, 400000, 500000],
            "SNPID": [f"1:{pos}_A_G" for pos in [100000, 200000, 300000, 400000, 500000]],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=50,  # 50kb window
            log=self.log,
            verbose=False
        )
        
        # With 50kb window, variants 100kb apart should not see each other
        # All variants should have density 0 (no neighbors within 50kb)
        densities = result["DENSITY"].values
        self.assertTrue(all(d == 0 for d in densities))
    
    def test_conditional_density_with_sig_sumstats(self):
        """Test conditional density calculation with significant variants."""
        # All variants
        df = pd.DataFrame({
            "CHR": [1] * 5,
            "POS": [100000, 150000, 200000, 250000, 300000],
            "SNPID": [f"1:{pos}_A_G" for pos in [100000, 150000, 200000, 250000, 300000]],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        
        # Significant variants (only at 100k and 200k)
        sig_df = pd.DataFrame({
            "CHR": [1, 1],
            "POS": [100000, 200000],
            "SNPID": ["1:100000_A_G", "1:200000_A_G"],
            "EA": ["A", "A"],
            "NEA": ["G", "G"]
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            sig_sumstats=sig_df,
            log=self.log,
            verbose=False
        )
        
        # With 100kb window (symmetric), each sig variant adds +1 to all variants in [sig_pos-wsize, sig_pos+wsize]
        # Sig at 100k: affects variants in [0, 200k] = 100k, 150k, 200k
        # Sig at 200k: affects variants in [100k, 300k] = 100k, 150k, 200k, 250k, 300k
        # - Variant at 100k: within window of sig at 100k = +1, within window of sig at 200k = +1, total = 2
        # - Variant at 150k: within window of sig at 100k = +1, within window of sig at 200k = +1, total = 2
        # - Variant at 200k: within window of sig at 100k = +1, within window of sig at 200k = +1, total = 2
        # - Variant at 250k: within window of sig at 200k = +1, total = 1
        # - Variant at 300k: within window of sig at 200k = +1, total = 1
        
        densities = result["DENSITY"].values
        self.assertEqual(densities[0], 2)  # 100k
        self.assertEqual(densities[1], 2)  # 150k
        self.assertEqual(densities[2], 2)  # 200k
        self.assertEqual(densities[3], 1)  # 250k
        self.assertEqual(densities[4], 1)  # 300k
    
    def test_conditional_density_multiple_chromosomes(self):
        """Test conditional density with significant variants across chromosomes."""
        df = pd.DataFrame({
            "CHR": [1, 1, 2, 2],
            "POS": [100000, 150000, 100000, 150000],
            "SNPID": ["1:100000_A_G", "1:150000_A_G", "2:100000_A_G", "2:150000_A_G"],
            "EA": ["A"] * 4,
            "NEA": ["G"] * 4
        })
        
        # Significant variant only on chromosome 1
        sig_df = pd.DataFrame({
            "CHR": [1],
            "POS": [100000],
            "SNPID": ["1:100000_A_G"],
            "EA": ["A"],
            "NEA": ["G"]
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            sig_sumstats=sig_df,
            log=self.log,
            verbose=False
        )
        
        # Chromosome 1 variants should have density > 0
        chr1_mask = result["CHR"] == 1
        chr1_densities = result.loc[chr1_mask, "DENSITY"].values
        self.assertTrue(all(d > 0 for d in chr1_densities))
        
        # Chromosome 2 variants should have density 0 (no sig variants on chr2)
        chr2_mask = result["CHR"] == 2
        chr2_densities = result.loc[chr2_mask, "DENSITY"].values
        self.assertTrue(all(d == 0 for d in chr2_densities))
    
    def test_single_variant(self):
        """Test density calculation with a single variant."""
        df = pd.DataFrame({
            "CHR": [1],
            "POS": [100000],
            "SNPID": ["1:100000_A_G"],
            "EA": ["A"],
            "NEA": ["G"]
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Single variant should have density 0 (no neighbors)
        self.assertEqual(result["DENSITY"].iloc[0], 0)
    
    def test_empty_dataframe(self):
        """Test density calculation with empty DataFrame."""
        df = pd.DataFrame(columns=["CHR", "POS", "SNPID", "EA", "NEA"])
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Should return empty DataFrame with DENSITY column
        self.assertIn("DENSITY", result.columns)
        self.assertEqual(len(result), 0)
    
    def test_missing_values(self):
        """Test density calculation with missing CHR or POS values."""
        df = pd.DataFrame({
            "CHR": [1, np.nan, 1, 1],
            "POS": [100000, 150000, np.nan, 200000],
            "SNPID": ["1:100000_A_G", "1:150000_A_G", "1:200000_A_G", "1:300000_A_G"],
            "EA": ["A"] * 4,
            "NEA": ["G"] * 4
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Variants with missing CHR or POS should have NA density
        # Note: After sort_values with ignore_index=True, indices are reset
        # So we need to check by position in the result
        result_sorted = result.sort_values(["CHR", "POS"], na_position="last")
        
        # Find rows with missing values
        missing_chr_idx = result_sorted[result_sorted["CHR"].isna()].index
        missing_pos_idx = result_sorted[result_sorted["POS"].isna()].index
        
        if len(missing_chr_idx) > 0:
            self.assertTrue(pd.isna(result_sorted.loc[missing_chr_idx[0], "DENSITY"]))
        if len(missing_pos_idx) > 0:
            self.assertTrue(pd.isna(result_sorted.loc[missing_pos_idx[0], "DENSITY"]))
        
        # Valid variants should have numeric density
        valid_mask = result_sorted["CHR"].notna() & result_sorted["POS"].notna()
        valid_densities = result_sorted.loc[valid_mask, "DENSITY"]
        self.assertTrue(all(not pd.isna(d) for d in valid_densities))
    
    def test_boundary_conditions_exact_window(self):
        """Test density calculation at exact window boundaries."""
        df = pd.DataFrame({
            "CHR": [1] * 3,
            "POS": [100000, 200000, 300000],  # Exactly 100kb apart
            "SNPID": ["1:100000_A_G", "1:200000_A_G", "1:300000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # With 100kb window, variants exactly 100kb apart should see each other
        # Variant at 100k: should see 200k (exactly 100kb away) = 1
        # Variant at 200k: should see 100k and 300k (both exactly 100kb away) = 2
        # Variant at 300k: should see 200k (exactly 100kb away) = 1
        densities = result["DENSITY"].values
        self.assertEqual(densities[0], 1)
        self.assertEqual(densities[1], 2)
        self.assertEqual(densities[2], 1)
    
    def test_conditional_density_no_overlap(self):
        """Test conditional density when significant variants are far from all variants."""
        df = pd.DataFrame({
            "CHR": [1] * 3,
            "POS": [100000, 200000, 300000],
            "SNPID": ["1:100000_A_G", "1:200000_A_G", "1:300000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        
        # Significant variant far away (1Mb away)
        sig_df = pd.DataFrame({
            "CHR": [1],
            "POS": [1100000],  # 1Mb away
            "SNPID": ["1:1100000_A_G"],
            "EA": ["A"],
            "NEA": ["G"]
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,  # 100kb window
            sig_sumstats=sig_df,
            log=self.log,
            verbose=False
        )
        
        # All variants should have density 0 (sig variant is too far)
        densities = result["DENSITY"].values
        self.assertTrue(all(d == 0 for d in densities))
    
    def test_conditional_density_all_variants_sig(self):
        """Test conditional density when all variants are significant."""
        df = pd.DataFrame({
            "CHR": [1] * 3,
            "POS": [100000, 150000, 200000],
            "SNPID": ["1:100000_A_G", "1:150000_A_G", "1:200000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        
        # All variants are significant
        sig_df = df.copy()
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            sig_sumstats=sig_df,
            log=self.log,
            verbose=False
        )
        
        # Each sig variant adds +1 to all variants within its window
        # Sig at 100k: affects [0, 200k] = 100k, 150k, 200k
        # Sig at 150k: affects [50k, 250k] = 100k, 150k, 200k
        # Sig at 200k: affects [100k, 300k] = 100k, 150k, 200k
        # - Variant at 100k: within window of all 3 sig variants = 3
        # - Variant at 150k: within window of all 3 sig variants = 3
        # - Variant at 200k: within window of all 3 sig variants = 3
        densities = result["DENSITY"].values
        self.assertEqual(densities[0], 3)
        self.assertEqual(densities[1], 3)
        self.assertEqual(densities[2], 3)
    
    def test_different_snpid_column_name(self):
        """Test density calculation with different SNPID column name."""
        df = pd.DataFrame({
            "CHR": [1] * 3,
            "POS": [100000, 150000, 200000],
            "rsID": ["rs1", "rs2", "rs3"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        
        result = _get_signal_density2(
            df,
            snpid="rsID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Should work with different column name
        self.assertIn("DENSITY", result.columns)
        self.assertEqual(len(result), 3)
    
    def test_return_type(self):
        """Test that function returns DataFrame with DENSITY column."""
        df = pd.DataFrame({
            "CHR": [1] * 3,
            "POS": [100000, 150000, 200000],
            "SNPID": ["1:100000_A_G", "1:150000_A_G", "1:200000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        
        result = _get_signal_density2(
            df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            bwindowsizekb=100,
            log=self.log,
            verbose=False
        )
        
        # Should return DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        
        # Should have DENSITY column
        self.assertIn("DENSITY", result.columns)
        
        # Should preserve original columns
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("SNPID", result.columns)
        
        # DENSITY should be Int32 type
        self.assertEqual(result["DENSITY"].dtype, "Int32")


if __name__ == "__main__":
    unittest.main()

