import os
import sys
import unittest
import numpy as np
import pandas as pd
import polars as pl

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_meta import meta_analyze_multi
from gwaslab.util.util_in_meta_polars import meta_analyze_polars
from gwaslab.info.g_Log import Log


class TestMetaAnalyzeConsistency(unittest.TestCase):
    """Test consistency between pandas and polars meta-analysis implementations"""
    
    def setUp(self):
        """Set up test data"""
        # Create test data with 2 studies, some overlapping variants, some with NA
        self.sumstats_multi_pd = pd.DataFrame({
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
        
        # Convert to polars for polars version
        self.sumstats_multi_pl = pl.from_pandas(self.sumstats_multi_pd)

    def compare_results(self, result_pd, result_pl, rtol=1e-5, atol=1e-8):
        """Compare results from pandas and polars implementations"""
        # Convert polars to pandas for comparison
        result_pl_pd = result_pl.to_pandas()
        
        # Sort both by CHR, POS for consistent comparison
        result_pd = result_pd.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result_pl_pd = result_pl_pd.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Check same number of variants
        self.assertEqual(len(result_pd), len(result_pl_pd), 
                         "Should have same number of variants")
        
        # Compare key columns
        numeric_cols = ["BETA", "SE", "P", "Z", "N", "EAF", "DOF", "Q", "I2"]
        for col in numeric_cols:
            if col in result_pd.columns and col in result_pl_pd.columns:
                pd_vals = result_pd[col].values
                pl_vals = result_pl_pd[col].values
                
                # Handle NaN values
                pd_nan = pd.isna(pd_vals)
                pl_nan = pd.isna(pl_vals)
                
                # Check NaN positions match
                self.assertTrue(np.array_equal(pd_nan, pl_nan),
                              f"NaN positions should match for {col}")
                
                # Compare non-NaN values
                if not np.all(pd_nan):
                    pd_vals_clean = pd_vals[~pd_nan]
                    pl_vals_clean = pl_vals[~pl_nan]
                    np.testing.assert_allclose(pd_vals_clean, pl_vals_clean, 
                                              rtol=rtol, atol=atol,
                                              err_msg=f"Values don't match for {col}")
        
        # Compare string columns
        string_cols = ["SNPID", "EA", "NEA", "DIRECTION"]
        for col in string_cols:
            if col in result_pd.columns and col in result_pl_pd.columns:
                pd_vals = result_pd[col].values
                pl_vals = result_pl_pd[col].values
                # Convert to string and handle NaN
                pd_vals_str = [str(v) if pd.notna(v) else "nan" for v in pd_vals]
                pl_vals_str = [str(v) if pd.notna(v) else "nan" for v in pl_vals]
                self.assertEqual(pd_vals_str, pl_vals_str,
                               f"Values don't match for {col}")

    def test_fixed_effects_consistency(self):
        """Test that fixed-effects meta-analysis produces consistent results"""
        # Run pandas version
        result_pd = meta_analyze_multi(
            self.sumstats_multi_pd.copy(),
            nstudy=2,
            random_effects=False
        )
        
        # Run polars version
        result_pl = meta_analyze_polars(
            self.sumstats_multi_pl.clone(),
            nstudy=2,
            random_effects=False
        )
        
        # Compare results
        self.compare_results(result_pd.data, result_pl)

    def test_random_effects_consistency(self):
        """Test that random-effects meta-analysis produces consistent results"""
        # Run pandas version
        result_pd = meta_analyze_multi(
            self.sumstats_multi_pd.copy(),
            nstudy=2,
            random_effects=True
        )
        
        # Run polars version
        result_pl = meta_analyze_polars(
            self.sumstats_multi_pl.clone(),
            nstudy=2,
            random_effects=True
        )
        
        # Compare fixed-effects results
        self.compare_results(result_pd.data, result_pl)
        
        # Compare random-effects specific columns
        result_pl_pd = result_pl.to_pandas()
        result_pd_sorted = result_pd.data.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result_pl_pd_sorted = result_pl_pd.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        random_cols = ["BETA_RANDOM", "SE_RANDOM", "Z_RANDOM", "P_RANDOM"]
        for col in random_cols:
            if col in result_pd_sorted.columns and col in result_pl_pd_sorted.columns:
                pd_vals = result_pd_sorted[col].values
                pl_vals = result_pl_pd_sorted[col].values
                
                # Handle NaN values
                # Note: pandas version only updates random effects for variants with _R2 > 0,
                # while polars updates all variants. For variants with tau^2 = 0, both should
                # give the same result (random effects = fixed effects)
                pd_nan = pd.isna(pd_vals)
                pl_nan = pd.isna(pl_vals)
                
                # For variants where both have values, they should match
                # Exclude inf values which can occur for single-study variants
                both_not_nan = ~pd_nan & ~pl_nan
                if np.any(both_not_nan):
                    pd_vals_clean = pd_vals[both_not_nan]
                    pl_vals_clean = pl_vals[both_not_nan]
                    
                    # Filter out inf values (can occur for single-study variants with DOF=0)
                    finite_mask = np.isfinite(pd_vals_clean) & np.isfinite(pl_vals_clean)
                    if np.any(finite_mask):
                        pd_finite = pd_vals_clean[finite_mask]
                        pl_finite = pl_vals_clean[finite_mask]
                        np.testing.assert_allclose(pd_finite, pl_finite, 
                                                  rtol=1e-5, atol=1e-8,
                                                  err_msg=f"Values don't match for {col} where both are finite")
                
                # For variants where pandas has NaN but polars doesn't, check if polars value
                # matches fixed effects (this happens when tau^2 = 0)
                pd_nan_pl_not = pd_nan & ~pl_nan
                if np.any(pd_nan_pl_not):
                    # Check if polars value matches fixed effects
                    for idx in np.where(pd_nan_pl_not)[0]:
                        if col == "BETA_RANDOM" and "BETA" in result_pd_sorted.columns:
                            # Random effects should equal fixed effects when tau^2 = 0
                            fixed_val = result_pd_sorted["BETA"].iloc[idx]
                            random_val = pl_vals[idx]
                            np.testing.assert_allclose(fixed_val, random_val, 
                                                      rtol=1e-5, atol=1e-8,
                                                      err_msg=f"BETA_RANDOM should equal BETA when tau^2=0 for variant {idx}")
                        elif col == "SE_RANDOM" and "SE" in result_pd_sorted.columns:
                            fixed_val = result_pd_sorted["SE"].iloc[idx]
                            random_val = pl_vals[idx]
                            np.testing.assert_allclose(fixed_val, random_val, 
                                                      rtol=1e-5, atol=1e-8,
                                                      err_msg=f"SE_RANDOM should equal SE when tau^2=0 for variant {idx}")

    def test_single_study_consistency(self):
        """Test consistency with single study"""
        # Run pandas version
        result_pd = meta_analyze_multi(
            self.sumstats_multi_pd.copy(),
            nstudy=1,
            random_effects=False
        )
        
        # Run polars version
        result_pl = meta_analyze_polars(
            self.sumstats_multi_pl.clone(),
            nstudy=1,
            random_effects=False
        )
        
        # Compare results
        self.compare_results(result_pd.data, result_pl)

    def test_na_handling_consistency(self):
        """Test that NA handling is consistent between implementations"""
        # Create data with more NA patterns
        sumstats_na_pd = pd.DataFrame({
            "SNPID": ["1:100:A:G", "1:200:G:C", "1:300:T:C"],
            "CHR": [1, 1, 1],
            "POS": [100, 200, 300],
            "EA": ["A", "G", "T"],
            "NEA": ["G", "C", "C"],
            "BETA_1": [0.1, -0.2, 0.3],
            "SE_1": [0.05, 0.1, 0.15],
            "N_1": [1000, 1000, 1000],
            "EAF_1": [0.3, 0.4, 0.5],
            "BETA_2": [0.12, np.nan, np.nan],  # Second variant has NA in study 2
            "SE_2": [0.06, np.nan, np.nan],
            "N_2": [2000, np.nan, np.nan],
            "EAF_2": [0.32, np.nan, np.nan],
        })
        
        sumstats_na_pl = pl.from_pandas(sumstats_na_pd)
        
        # Run both versions
        result_pd = meta_analyze_multi(
            sumstats_na_pd.copy(),
            nstudy=2,
            random_effects=False
        )
        
        result_pl = meta_analyze_polars(
            sumstats_na_pl.clone(),
            nstudy=2,
            random_effects=False
        )
        
        # Compare results
        self.compare_results(result_pd.data, result_pl)

    def test_three_studies_consistency(self):
        """Test consistency with three studies"""
        # Create data with 3 studies
        sumstats_3_pd = pd.DataFrame({
            "SNPID": ["1:100:A:G", "1:200:G:C"],
            "CHR": [1, 1],
            "POS": [100, 200],
            "EA": ["A", "G"],
            "NEA": ["G", "C"],
            "BETA_1": [0.1, -0.2],
            "SE_1": [0.05, 0.1],
            "N_1": [1000, 1000],
            "EAF_1": [0.3, 0.4],
            "BETA_2": [0.12, -0.18],
            "SE_2": [0.06, 0.11],
            "N_2": [2000, 2000],
            "EAF_2": [0.32, 0.38],
            "BETA_3": [0.11, -0.19],
            "SE_3": [0.055, 0.105],
            "N_3": [1500, 1500],
            "EAF_3": [0.31, 0.39],
        })
        
        sumstats_3_pl = pl.from_pandas(sumstats_3_pd)
        
        # Run both versions
        result_pd = meta_analyze_multi(
            sumstats_3_pd.copy(),
            nstudy=3,
            random_effects=False
        )
        
        result_pl = meta_analyze_polars(
            sumstats_3_pl.clone(),
            nstudy=3,
            random_effects=False
        )
        
        # Compare results
        self.compare_results(result_pd.data, result_pl)

    def test_w2_sum_consistency(self):
        """Test that _W2_SUM calculation is consistent (tests the bug fix)"""
        # This test ensures both implementations correctly calculate _W2_SUM
        # We can't directly access _W2_SUM in the results, but we can verify
        # through random effects results which depend on it
        
        result_pd = meta_analyze_multi(
            self.sumstats_multi_pd.copy(),
            nstudy=2,
            random_effects=True
        )
        
        result_pl = meta_analyze_polars(
            self.sumstats_multi_pl.clone(),
            nstudy=2,
            random_effects=True
        )
        
        # Compare random effects results (which depend on _W2_SUM)
        result_pl_pd = result_pl.to_pandas()
        result_pd_sorted = result_pd.data.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result_pl_pd_sorted = result_pl_pd.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Check that variants with heterogeneity have consistent random effects
        # Variants in both studies should have tau^2 > 0 (if there's heterogeneity)
        variants_both = result_pd_sorted[
            (result_pd_sorted["SNPID"] == "1:100:A:G") | 
            (result_pd_sorted["SNPID"] == "1:200:G:C") | 
            (result_pd_sorted["SNPID"] == "2:300:T:C")
        ]
        
        if len(variants_both) > 0:
            for idx, row in variants_both.iterrows():
                snpid = row["SNPID"]
                pl_row = result_pl_pd_sorted[result_pl_pd_sorted["SNPID"] == snpid]
                
                if len(pl_row) > 0:
                    # Compare BETA_RANDOM and SE_RANDOM
                    if pd.notna(row["BETA_RANDOM"]) and pd.notna(pl_row["BETA_RANDOM"].iloc[0]):
                        np.testing.assert_allclose(
                            row["BETA_RANDOM"], 
                            pl_row["BETA_RANDOM"].iloc[0],
                            rtol=1e-5, atol=1e-8,
                            err_msg=f"BETA_RANDOM doesn't match for {snpid}"
                        )
                    
                    if pd.notna(row["SE_RANDOM"]) and pd.notna(pl_row["SE_RANDOM"].iloc[0]):
                        np.testing.assert_allclose(
                            row["SE_RANDOM"], 
                            pl_row["SE_RANDOM"].iloc[0],
                            rtol=1e-5, atol=1e-8,
                            err_msg=f"SE_RANDOM doesn't match for {snpid}"
                        )

    def test_heterogeneity_statistics_consistency(self):
        """Test that heterogeneity statistics (Q, I2, P_HET) are consistent"""
        result_pd = meta_analyze_multi(
            self.sumstats_multi_pd.copy(),
            nstudy=2,
            random_effects=False
        )
        
        result_pl = meta_analyze_polars(
            self.sumstats_multi_pl.clone(),
            nstudy=2,
            random_effects=False
        )
        
        result_pl_pd = result_pl.to_pandas()
        result_pd_sorted = result_pd.data.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result_pl_pd_sorted = result_pl_pd.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare Q, I2, P_HET for variants in both studies
        het_cols = ["Q", "I2", "P_HET"]
        for col in het_cols:
            if col in result_pd_sorted.columns and col in result_pl_pd_sorted.columns:
                pd_vals = result_pd_sorted[col].values
                pl_vals = result_pl_pd_sorted[col].values
                
                # Handle NaN values
                pd_nan = pd.isna(pd_vals)
                pl_nan = pd.isna(pl_vals)
                self.assertTrue(np.array_equal(pd_nan, pl_nan),
                              f"NaN positions should match for {col}")
                
                # Compare non-NaN values
                if not np.all(pd_nan):
                    pd_vals_clean = pd_vals[~pd_nan]
                    pl_vals_clean = pl_vals[~pl_nan]
                    np.testing.assert_allclose(pd_vals_clean, pl_vals_clean, 
                                              rtol=1e-5, atol=1e-8,
                                              err_msg=f"Values don't match for {col}")


if __name__ == "__main__":
    unittest.main()
