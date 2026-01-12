import os
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_meta import meta_analyze_multi
from gwaslab.info.g_Log import Log


class TestMetaAnalyzeMulti(unittest.TestCase):
    def setUp(self):
        """Set up test data"""
        # Create test data with 2 studies
        self.sumstats_multi = pd.DataFrame({
            "SNPID": ["1:100:A:G", "1:200:G:C", "2:300:T:C", "2:400:C:T", "1:500:A:T"],
            "CHR": [1, 1, 2, 2, 1],
            "POS": [100, 200, 300, 400, 500],
            "EA": ["A", "G", "T", "C", "A"],
            "NEA": ["G", "C", "C", "T", "T"],
            "BETA_1": [0.1, -0.2, 0.3, 0.05, 0.15],
            "SE_1": [0.05, 0.1, 0.15, 0.02, 0.08],
            "N_1": [1000, 1000, 1000, 1000, 1000],
            "EAF_1": [0.3, 0.4, 0.5, 0.2, 0.35],
            "BETA_2": [0.12, -0.18, 0.28, np.nan, np.nan],
            "SE_2": [0.06, 0.11, 0.16, np.nan, np.nan],
            "N_2": [2000, 2000, 2000, np.nan, np.nan],
            "EAF_2": [0.32, 0.38, 0.52, np.nan, np.nan],
        })

    def test_basic_meta_analysis(self):
        """Test basic fixed-effect meta-analysis"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=False
        )
        
        # Check result is Sumstats object
        self.assertIsNotNone(result)
        self.assertIsNotNone(result.data)
        
        # Check expected columns
        expected_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF", "Z", "Q", "DOF", "DIRECTION"]
        for col in expected_cols:
            self.assertIn(col, result.data.columns, f"Column {col} not found")
        
        # Check that variants present in both studies have DOF=1 (2 studies - 1 = 1)
        # DOF is initialized to -1, then +1 for each study with data
        # For 2 studies: DOF = -1 + 1 + 1 = 1
        variants_both = result.data[
            (result.data["SNPID"] == "1:100:A:G") | 
            (result.data["SNPID"] == "1:200:G:C") | 
            (result.data["SNPID"] == "2:300:T:C")
        ]
        if len(variants_both) > 0:
            self.assertTrue(all(variants_both["DOF"] == 1), 
                         f"Variants in both studies should have DOF=1 (2 studies - 1), got {variants_both['DOF'].values}")
        
        # Check that variants only in study 1 have DOF=0 (1 study - 1 = 0)
        variants_one = result.data[result.data["SNPID"] == "1:500:A:T"]
        if len(variants_one) > 0:
            self.assertEqual(variants_one["DOF"].iloc[0], 0, 
                            f"Variant only in study 1 should have DOF=0 (1 study - 1), got {variants_one['DOF'].iloc[0]}")

    def test_na_handling(self):
        """Test that NA values are handled correctly"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=False
        )
        
        # Variants with NA in study 2 should only use study 1 data
        variant_na = result.data[result.data["SNPID"] == "2:400:C:T"]
        if len(variant_na) > 0:
            # DOF should be 0 (1 study - 1 = 0, since study 2 has NA so not counted)
            dof_value = variant_na["DOF"].iloc[0]
            self.assertEqual(dof_value, 0, 
                            f"Variant with NA in study 2 should have DOF=0 (1 study - 1), got {dof_value}")
            self.assertEqual(variant_na["N"].iloc[0], 1000, "N should only come from study 1")
            # BETA and SE should match study 1 (within floating point precision)
            self.assertAlmostEqual(variant_na["BETA"].iloc[0], 0.05, places=5)
            self.assertAlmostEqual(variant_na["SE"].iloc[0], 0.02, places=5)

    def test_weighted_mean_calculation(self):
        """Test that weighted mean is calculated correctly"""
        # For variant 1:100:A:G
        # Study 1: BETA=0.1, SE=0.05, weight = 1/(0.05^2) = 400
        # Study 2: BETA=0.12, SE=0.06, weight = 1/(0.06^2) = 277.78
        # Weighted mean = (0.1*400 + 0.12*277.78) / (400 + 277.78) ≈ 0.108
        
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=False
        )
        
        variant = result.data[result.data["SNPID"] == "1:100:A:G"]
        if len(variant) > 0:
            # Check that BETA is reasonable (weighted average)
            self.assertGreater(variant["BETA"].iloc[0], 0.1)
            self.assertLess(variant["BETA"].iloc[0], 0.12)
            # Check that SE is smaller than individual SEs (meta-analysis should be more precise)
            self.assertLess(variant["SE"].iloc[0], 0.05)
            self.assertLess(variant["SE"].iloc[0], 0.06)

    def test_duplicate_handling(self):
        """Test that duplicate variants are handled correctly"""
        # Add duplicate variant
        sumstats_with_dup = pd.concat([
            self.sumstats_multi,
            pd.DataFrame({
                "SNPID": ["1:100:A:G"],  # Duplicate
                "CHR": [1],
                "POS": [100],
                "EA": ["A"],
                "NEA": ["G"],
                "BETA_1": [0.11],  # Slightly different
                "SE_1": [0.051],
                "N_1": [1001],
                "EAF_1": [0.31],
                "BETA_2": [0.13],
                "SE_2": [0.061],
                "N_2": [2001],
                "EAF_2": [0.33],
            })
        ], ignore_index=True)
        
        result = meta_analyze_multi(
            sumstats_with_dup,
            nstudy=2,
            random_effects=False
        )
        
        # Should only have one row for 1:100:A:G
        variant_count = (result.data["SNPID"] == "1:100:A:G").sum()
        self.assertEqual(variant_count, 1, "Duplicate variants should be removed")

    def test_random_effects(self):
        """Test random-effects meta-analysis"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=True
        )
        
        # Check random effects columns exist
        self.assertIn("BETA_RANDOM", result.data.columns)
        self.assertIn("SE_RANDOM", result.data.columns)
        self.assertIn("Z_RANDOM", result.data.columns)
        self.assertIn("P_RANDOM", result.data.columns)
        
        # Random effects SE should be >= fixed effects SE
        variants_both = result.data[
            (result.data["SNPID"] == "1:100:A:G") | 
            (result.data["SNPID"] == "1:200:G:C") | 
            (result.data["SNPID"] == "2:300:T:C")
        ]
        if len(variants_both) > 0:
            # Random effects typically have larger SE due to between-study variance
            for idx, row in variants_both.iterrows():
                if pd.notna(row["SE_RANDOM"]):
                    self.assertGreaterEqual(row["SE_RANDOM"], row["SE"], 
                                          "Random effects SE should be >= fixed effects SE")

    def test_single_study(self):
        """Test meta-analysis with single study"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=1,
            random_effects=False
        )
        
        # All variants should have DOF=0 (1 study - 1 = 0)
        # DOF is initialized to -1, then +1 for study 1 = 0
        dof_values = result.data["DOF"].values
        self.assertTrue(all(dof_values == 0), 
                        f"Single study should have DOF=0 (1 study - 1), got {dof_values}")
        
        # BETA and SE should match study 1 values
        variant = result.data[result.data["SNPID"] == "1:100:A:G"]
        if len(variant) > 0:
            self.assertAlmostEqual(variant["BETA"].iloc[0], 0.1, places=5)
            self.assertAlmostEqual(variant["SE"].iloc[0], 0.05, places=5)

    def test_all_na_in_study(self):
        """Test handling when all values are NA for a study"""
        sumstats_all_na = pd.DataFrame({
            "SNPID": ["1:100:A:G"],
            "CHR": [1],
            "POS": [100],
            "EA": ["A"],
            "NEA": ["G"],
            "BETA_1": [0.1],
            "SE_1": [0.05],
            "N_1": [1000],
            "EAF_1": [0.3],
            "BETA_2": [np.nan],
            "SE_2": [np.nan],
            "N_2": [np.nan],
            "EAF_2": [np.nan],
        })
        
        result = meta_analyze_multi(
            sumstats_all_na,
            nstudy=2,
            random_effects=False
        )
        
        variant = result.data[result.data["SNPID"] == "1:100:A:G"]
        if len(variant) > 0:
            # DOF should be 0 (1 study - 1 = 0, since study 2 has all NA so not counted)
            dof_value = variant["DOF"].iloc[0]
            self.assertEqual(dof_value, 0, f"Should have DOF=0 (1 study - 1), got {dof_value}")
            self.assertAlmostEqual(variant["BETA"].iloc[0], 0.1, places=5)

    def test_heterogeneity_statistics(self):
        """Test that heterogeneity statistics are calculated"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=False
        )
        
        # Check heterogeneity columns
        self.assertIn("Q", result.data.columns)
        self.assertIn("P_HET", result.data.columns)
        self.assertIn("I2", result.data.columns)
        
        # Variants in both studies should have Q, P_HET, I2
        variants_both = result.data[
            (result.data["SNPID"] == "1:100:A:G") | 
            (result.data["SNPID"] == "1:200:G:C") | 
            (result.data["SNPID"] == "2:300:T:C")
        ]
        if len(variants_both) > 0:
            self.assertTrue(all(pd.notna(variants_both["Q"])), "Q should be calculated")
            self.assertTrue(all(pd.notna(variants_both["P_HET"])), "P_HET should be calculated")
            self.assertTrue(all(pd.notna(variants_both["I2"])), "I2 should be calculated")
            # I2 should be between 0 and 1 (or 0 and 100 if percentage)
            self.assertTrue(all((variants_both["I2"] >= 0) & (variants_both["I2"] <= 1)), 
                          "I2 should be between 0 and 1")

    def test_direction_column(self):
        """Test DIRECTION column is correctly populated"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=False
        )
        
        self.assertIn("DIRECTION", result.data.columns)
        
        # Check that direction is set for variants
        variant_pos = result.data[result.data["SNPID"] == "1:100:A:G"]
        variant_neg = result.data[result.data["SNPID"] == "1:200:G:C"]
        
        if len(variant_pos) > 0:
            self.assertIn("+", variant_pos["DIRECTION"].iloc[0], "Positive BETA should have +")
        if len(variant_neg) > 0:
            self.assertIn("-", variant_neg["DIRECTION"].iloc[0], "Negative BETA should have -")

    def test_eaf_calculation(self):
        """Test that EAF is calculated as weighted average"""
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=False
        )
        
        # For variant 1:100:A:G
        # Study 1: N=1000, EAF=0.3 -> EA_N = 300, NEA_N = 700
        # Study 2: N=2000, EAF=0.32 -> EA_N = 640, NEA_N = 1360
        # Total: EA_N = 940, NEA_N = 2060, EAF = 940/(940+2060) = 0.3133
        
        variant = result.data[result.data["SNPID"] == "1:100:A:G"]
        if len(variant) > 0:
            eaf = variant["EAF"].iloc[0]
            self.assertGreater(eaf, 0.3)
            self.assertLess(eaf, 0.32)
            # Should be closer to study 2 since it has larger N

    def test_w2_sum_calculation(self):
        """Test that _W2_SUM is calculated correctly (sum of squared weights)"""
        # This is indirectly tested through random effects, but we can verify
        # that random effects calculation doesn't fail
        result = meta_analyze_multi(
            self.sumstats_multi.copy(),
            nstudy=2,
            random_effects=True
        )
        
        # If _W2_SUM was wrong, tau^2 calculation would be wrong
        # Check that random effects results are reasonable
        variants_both = result.data[
            (result.data["SNPID"] == "1:100:A:G") | 
            (result.data["SNPID"] == "1:200:G:C")
        ]
        if len(variants_both) > 0:
            # Random effects BETA should be similar to fixed effects (not wildly different)
            for idx, row in variants_both.iterrows():
                if pd.notna(row["BETA_RANDOM"]) and pd.notna(row["BETA"]):
                    diff = abs(row["BETA_RANDOM"] - row["BETA"])
                    self.assertLess(diff, 0.1, "Random and fixed effects BETA should be similar")


if __name__ == "__main__":
    unittest.main()
