import os
import sys
import unittest
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl


class TestSumstatsMultiKeepAllVariants(unittest.TestCase):
    """Test that SumstatsMulti uses the same merging logic as SumstatsPair"""
    
    def setUp(self):
        """Set up test data with overlapping and non-overlapping variants"""
        # Study 1: variants at positions 100, 200, 300
        self.sumstats1 = gl.Sumstats(
            pd.DataFrame({
                "SNPID": ["1:100:A:G", "1:200:G:C", "1:300:T:C"],
                "CHR": [1, 1, 1],
                "POS": [100, 200, 300],
                "EA": ["A", "G", "T"],
                "NEA": ["G", "C", "C"],
                "BETA": [0.1, -0.2, 0.3],
                "SE": [0.05, 0.1, 0.15],
                "P": [0.05, 0.05, 0.05],
                "N": [1000, 1000, 1000],
                "EAF": [0.3, 0.4, 0.5],
            }),
            verbose=False
        )
        
        # Study 2: variants at positions 100, 200, 400 (overlaps at 100, 200; unique at 400)
        self.sumstats2 = gl.Sumstats(
            pd.DataFrame({
                "SNPID": ["1:100:A:G", "1:200:G:C", "1:400:C:T"],
                "CHR": [1, 1, 1],
                "POS": [100, 200, 400],
                "EA": ["A", "G", "C"],
                "NEA": ["G", "C", "T"],
                "BETA": [0.12, -0.18, 0.05],
                "SE": [0.06, 0.11, 0.02],
                "P": [0.04, 0.10, 0.01],
                "N": [2000, 2000, 2000],
                "EAF": [0.32, 0.38, 0.2],
            }),
            verbose=False
        )
        
        # Study 3: variant at position 500 (unique)
        self.sumstats3 = gl.Sumstats(
            pd.DataFrame({
                "SNPID": ["1:500:A:T"],
                "CHR": [1],
                "POS": [500],
                "EA": ["A"],
                "NEA": ["T"],
                "BETA": [0.15],
                "SE": [0.08],
                "P": [0.06],
                "N": [1500],
                "EAF": [0.35],
            }),
            verbose=False
        )

    def test_keep_all_variants_false(self):
        """Test that keep_all_variants=False only keeps overlapping variants"""
        multi = gl.SumstatsMulti(
            [self.sumstats1, self.sumstats2],
            keep_all_variants=False,
            verbose=False
        )
        
        # Should only have variants present in both studies (100, 200)
        result_positions = set(multi.data["POS"].values)
        self.assertEqual(result_positions, {100, 200}, 
                        "Should only have overlapping variants when keep_all_variants=False")
        self.assertEqual(len(multi.data), 2, 
                        "Should have 2 variants (both present in both studies)")

    def test_keep_all_variants_true(self):
        """Test that keep_all_variants=True keeps all variants"""
        multi = gl.SumstatsMulti(
            [self.sumstats1, self.sumstats2],
            keep_all_variants=True,
            verbose=False
        )
        
        # Should have all variants from both studies (100, 200, 300, 400)
        result_positions = set(multi.data["POS"].values)
        self.assertEqual(result_positions, {100, 200, 300, 400}, 
                        "Should have all variants when keep_all_variants=True")
        self.assertEqual(len(multi.data), 4, 
                        "Should have 4 variants (2 overlapping + 1 from study1 + 1 from study2)")

    def test_default_keep_all_variants(self):
        """Test that default behavior (keep_all_variants=True) keeps all variants"""
        multi = gl.SumstatsMulti(
            [self.sumstats1, self.sumstats2],
            verbose=False
        )
        
        # Should have all variants from both studies (100, 200, 300, 400) by default
        result_positions = set(multi.data["POS"].values)
        self.assertEqual(result_positions, {100, 200, 300, 400}, 
                        "Should have all variants by default (keep_all_variants=True)")
        self.assertEqual(len(multi.data), 4, 
                        "Should have 4 variants (2 overlapping + 1 from study1 + 1 from study2)")

    def test_three_studies_keep_all_variants(self):
        """Test keep_all_variants with three studies"""
        multi = gl.SumstatsMulti(
            [self.sumstats1, self.sumstats2, self.sumstats3],
            keep_all_variants=True,
            verbose=False
        )
        
        # Should have all unique variants: 100, 200, 300, 400, 500
        result_positions = set(multi.data["POS"].values)
        self.assertEqual(result_positions, {100, 200, 300, 400, 500}, 
                        "Should have all variants from all three studies")
        self.assertEqual(len(multi.data), 5, 
                        "Should have 5 variants total")

    def test_consistency_with_pair(self):
        """Test that SumstatsMulti and SumstatsPair produce consistent results"""
        # Create pair
        pair = gl.SumstatsPair(
            self.sumstats1, 
            self.sumstats2, 
            keep_all_variants=False,
            verbose=False
        )
        
        # Create multi
        multi = gl.SumstatsMulti(
            [self.sumstats1, self.sumstats2],
            keep_all_variants=False,
            verbose=False
        )
        
        # Both should have same number of variants
        self.assertEqual(len(pair.data), len(multi.data),
                        "SumstatsPair and SumstatsMulti should have same number of variants")
        
        # Both should have same positions
        pair_positions = set(pair.data["POS"].values)
        multi_positions = set(multi.data["POS"].values)
        self.assertEqual(pair_positions, multi_positions,
                        "SumstatsPair and SumstatsMulti should have same variant positions")


if __name__ == "__main__":
    unittest.main()
