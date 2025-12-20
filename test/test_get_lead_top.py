"""
Test suite for get_lead and get_top functions with comprehensive corner cases.

Tests cover:
- All variants are significant
- Single row
- Duplicated positions
- Empty DataFrame
- No significant variants
- Identical P values
- NaN values
- Multiple chromosomes
- Edge cases for get_top
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

from gwaslab.g_Sumstats import Sumstats


class TestGetLeadCornerCases(unittest.TestCase):
    """Test corner cases for get_lead function."""
    
    def test_single_row(self):
        """Test get_lead with a single row."""
        df = pd.DataFrame({
            "CHR": [1],
            "POS": [1000000],
            "P": [1e-10],
            "MLOG10P": [10.0],
            "SNPID": ["1:1000000_A_G"],
            "EA": ["A"],
            "NEA": ["G"]
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P", 
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["SNPID"], "1:1000000_A_G")
    
    def test_all_significant(self):
        """Test get_lead when all variants are significant."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 2, 2],
            "POS": [1000000, 1001000, 2000000, 500000, 1500000],
            "P": [1e-10, 1e-9, 1e-8, 1e-11, 1e-9],
            "MLOG10P": [10.0, 9.0, 8.0, 11.0, 9.0],
            "SNPID": [f"1:{pos}_A_G" if chr == 1 else f"2:{pos}_A_G" 
                     for chr, pos in zip([1, 1, 1, 2, 2], [1000000, 1001000, 2000000, 500000, 1500000])],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=1e-7, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should cluster variants within window and select leads
        self.assertGreater(len(result), 0)
        self.assertLessEqual(len(result), len(df))
    
    def test_duplicated_positions(self):
        """Test get_lead with duplicated positions."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "POS": [1000000, 1000000, 1000000, 2000000],
            "P": [1e-10, 1e-9, 1e-8, 1e-11],
            "MLOG10P": [10.0, 9.0, 8.0, 11.0],
            "SNPID": ["1:1000000_A_G", "1:1000000_T_C", "1:1000000_G_A", "1:2000000_A_G"],
            "EA": ["A", "T", "G", "A"],
            "NEA": ["G", "C", "A", "G"]
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should handle duplicates and select the most significant one at each position
        self.assertGreater(len(result), 0)
        # The most significant at 1000000 should be included
        pos_1m_results = result[result["POS"] == 1000000]
        if len(pos_1m_results) > 0:
            # Should have the most significant P value
            self.assertEqual(pos_1m_results.iloc[0]["P"], 1e-10)
    
    def test_no_significant_variants(self):
        """Test get_lead when no variants are significant."""
        df = pd.DataFrame({
            "CHR": [1, 1, 2],
            "POS": [1000000, 2000000, 500000],
            "P": [0.1, 0.05, 0.01],
            "MLOG10P": [1.0, 1.3, 2.0],
            "SNPID": ["1:1000000_A_G", "1:2000000_A_G", "2:500000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        # Should return empty DataFrame or None
        if result is not None:
            self.assertEqual(len(result), 0)
    
    def test_identical_p_values(self):
        """Test get_lead with identical P values."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 1001000, 2000000],
            "P": [1e-10, 1e-10, 1e-10],
            "MLOG10P": [10.0, 10.0, 10.0],
            "SNPID": ["1:1000000_A_G", "1:1001000_A_G", "1:2000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should still select leads even with identical P values
        self.assertGreater(len(result), 0)
    
    def test_nan_values(self):
        """Test get_lead with NaN values in P and MLOG10P."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "POS": [1000000, 2000000, 3000000, 4000000],
            "P": [1e-10, np.nan, 1e-9, 1e-8],
            "MLOG10P": [10.0, np.nan, 9.0, 8.0],
            "SNPID": [f"1:{pos}_A_G" for pos in [1000000, 2000000, 3000000, 4000000]],
            "EA": ["A"] * 4,
            "NEA": ["G"] * 4
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should only include variants with valid P values
        self.assertGreater(len(result), 0)
        # Should not include NaN values
        self.assertFalse(result["P"].isna().any())
    
    def test_multiple_chromosomes(self):
        """Test get_lead across multiple chromosomes."""
        df = pd.DataFrame({
            "CHR": [1, 1, 2, 2, 3],
            "POS": [1000000, 2000000, 500000, 1500000, 800000],
            "P": [1e-10, 1e-9, 1e-11, 1e-8, 1e-9],
            "MLOG10P": [10.0, 9.0, 11.0, 8.0, 9.0],
            "SNPID": [f"{chr}:{pos}_A_G" for chr, pos in 
                      zip([1, 1, 2, 2, 3], [1000000, 2000000, 500000, 1500000, 800000])],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should handle multiple chromosomes correctly
        self.assertGreater(len(result), 0)
        unique_chrs = result["CHR"].unique()
        self.assertGreater(len(unique_chrs), 0)
    
    def test_empty_dataframe(self):
        """Test get_lead with empty DataFrame."""
        df = pd.DataFrame(columns=["CHR", "POS", "P", "MLOG10P", "SNPID", "EA", "NEA"])
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        # Should return empty DataFrame or None
        if result is not None:
            self.assertEqual(len(result), 0)
    
    def test_use_p_flag(self):
        """Test get_lead with use_p=True flag."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 1001000, 2000000],
            "P": [1e-10, 1e-9, 1e-8],
            "MLOG10P": [10.0, 9.0, 8.0],
            "SNPID": ["1:1000000_A_G", "1:1001000_A_G", "1:2000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=500, use_p=True, verbose=False)
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)
    
    def test_small_window_size(self):
        """Test get_lead with very small window size."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 1000100, 2000000],  # First two within 1kb
            "P": [1e-10, 1e-9, 1e-8],
            "MLOG10P": [10.0, 9.0, 8.0],
            "SNPID": ["1:1000000_A_G", "1:1000100_A_G", "1:2000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
                     snpid="SNPID", ea="EA", nea="NEA", verbose=False)
        result = gl.get_lead(sig_level=5e-8, windowsizekb=1, verbose=False)
        self.assertIsNotNone(result)
        # With 1kb window, first two should be in same cluster
        self.assertGreater(len(result), 0)


class TestGetTopCornerCases(unittest.TestCase):
    """Test corner cases for get_top function."""
    
    def test_single_row(self):
        """Test get_top with a single row."""
        df = pd.DataFrame({
            "CHR": [1],
            "POS": [1000000],
            "DENSITY": [10.5],
            "SNPID": ["1:1000000_A_G"],
            "EA": ["A"],
            "NEA": ["G"]
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID", 
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["SNPID"], "1:1000000_A_G")
    
    def test_all_same_density(self):
        """Test get_top when all variants have the same density."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 1001000, 2000000],
            "DENSITY": [5.0, 5.0, 5.0],
            "SNPID": ["1:1000000_A_G", "1:1001000_A_G", "1:2000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should still select leads even with identical density
        self.assertGreater(len(result), 0)
    
    def test_duplicated_positions(self):
        """Test get_top with duplicated positions."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "POS": [1000000, 1000000, 1000000, 2000000],
            "DENSITY": [10.0, 8.0, 9.0, 11.0],
            "SNPID": ["1:1000000_A_G", "1:1000000_T_C", "1:1000000_G_A", "1:2000000_A_G"],
            "EA": ["A", "T", "G", "A"],
            "NEA": ["G", "C", "A", "G"]
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should handle duplicates and select the one with highest density
        self.assertGreater(len(result), 0)
        pos_1m_results = result[result["POS"] == 1000000]
        if len(pos_1m_results) > 0:
            # Should have the highest density
            self.assertEqual(pos_1m_results.iloc[0]["DENSITY"], 10.0)
    
    def test_no_variants_above_threshold(self):
        """Test get_top when no variants pass the threshold."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 2000000, 3000000],
            "DENSITY": [1.0, 2.0, 3.0],
            "SNPID": ["1:1000000_A_G", "1:2000000_A_G", "1:3000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", threshold=100.0, windowsizekb=500, verbose=False)
        # Should return None or empty DataFrame when threshold too high
        if result is not None:
            self.assertEqual(len(result), 0)
    
    def test_nan_values(self):
        """Test get_top with NaN values in the metric column."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "POS": [1000000, 2000000, 3000000, 4000000],
            "DENSITY": [10.0, np.nan, 9.0, 8.0],
            "SNPID": [f"1:{pos}_A_G" for pos in [1000000, 2000000, 3000000, 4000000]],
            "EA": ["A"] * 4,
            "NEA": ["G"] * 4
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should only include variants with valid density values
        self.assertGreater(len(result), 0)
        # Should not include NaN values
        self.assertFalse(result["DENSITY"].isna().any())
    
    def test_multiple_chromosomes(self):
        """Test get_top across multiple chromosomes."""
        df = pd.DataFrame({
            "CHR": [1, 1, 2, 2, 3],
            "POS": [1000000, 2000000, 500000, 1500000, 800000],
            "DENSITY": [10.0, 9.0, 11.0, 8.0, 9.5],
            "SNPID": [f"{chr}:{pos}_A_G" for chr, pos in 
                      zip([1, 1, 2, 2, 3], [1000000, 2000000, 500000, 1500000, 800000])],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should handle multiple chromosomes correctly
        self.assertGreater(len(result), 0)
        unique_chrs = result["CHR"].unique()
        self.assertGreater(len(unique_chrs), 0)
    
    def test_empty_dataframe(self):
        """Test get_top with empty DataFrame."""
        df = pd.DataFrame(columns=["CHR", "POS", "DENSITY", "SNPID", "EA", "NEA"])
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        # Should return None or empty DataFrame
        if result is not None:
            self.assertEqual(len(result), 0)
    
    def test_custom_metric_column(self):
        """Test get_top with a custom metric column."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 1001000, 2000000],
            "CUSTOM_METRIC": [100.0, 90.0, 80.0],
            "SNPID": ["1:1000000_A_G", "1:1001000_A_G", "1:2000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="CUSTOM_METRIC", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)
        self.assertIn("CUSTOM_METRIC", result.columns)
    
    def test_small_window_size(self):
        """Test get_top with very small window size."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000000, 1000100, 2000000],  # First two within 1kb
            "DENSITY": [10.0, 9.0, 8.0],
            "SNPID": ["1:1000000_A_G", "1:1000100_A_G", "1:2000000_A_G"],
            "EA": ["A"] * 3,
            "NEA": ["G"] * 3
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=1, verbose=False)
        self.assertIsNotNone(result)
        # With 1kb window, first two should be in same cluster
        self.assertGreater(len(result), 0)
    
    def test_all_variants_high_density(self):
        """Test get_top when all variants have high density."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 2, 2],
            "POS": [1000000, 1001000, 2000000, 500000, 1500000],
            "DENSITY": [50.0, 45.0, 40.0, 55.0, 48.0],
            "SNPID": [f"{chr}:{pos}_A_G" if chr == 1 else f"{chr}:{pos}_A_G"
                     for chr, pos in zip([1, 1, 1, 2, 2], [1000000, 1001000, 2000000, 500000, 1500000])],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5
        })
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", snpid="SNPID",
                     ea="EA", nea="NEA", verbose=False)
        result = gl.get_top(by="DENSITY", windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        # Should cluster variants within window and select top ones
        self.assertGreater(len(result), 0)
        self.assertLessEqual(len(result), len(df))


if __name__ == "__main__":
    unittest.main()

