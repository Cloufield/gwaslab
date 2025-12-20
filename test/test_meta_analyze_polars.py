import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
import polars as pl

from gwaslab.g_Sumstats_polars import Sumstatsp
from gwaslab.g_SumstatsMulti import SumstatsMulti
from gwaslab.util.util_in_meta_polars import meta_analyze_polars


def make_sumstatsp(df, study_name):
    """Helper function to create a Sumstatsp object"""
    return Sumstatsp(
        sumstats=df,
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        snpid="SNPID",
        beta="BETA",
        se="SE",
        p="P",
        n="N",
        eaf="EAF",
        verbose=False
    )


class TestMetaAnalyzePolars(unittest.TestCase):
    def setUp(self):
        """Set up test data with two studies"""
        # Study 1 data
        df1 = pd.DataFrame({
            "CHR": [1, 1, 2, 2],
            "POS": [100, 200, 300, 400],
            "EA": ["A", "G", "T", "C"],
            "NEA": ["G", "C", "C", "T"],
            "SNPID": ["1:100:A:G", "1:200:G:C", "2:300:T:C", "2:400:C:T"],
            "BETA": [0.1, -0.2, 0.3, 0.05],
            "SE": [0.05, 0.1, 0.15, 0.02],
            "P": [0.05, 0.05, 0.05, 0.01],
            "N": [1000, 1000, 1000, 1000],
            "EAF": [0.3, 0.4, 0.5, 0.2],
        })
        
        # Study 2 data - overlapping variants
        df2 = pd.DataFrame({
            "CHR": [1, 1, 2],
            "POS": [100, 200, 300],
            "EA": ["A", "G", "T"],
            "NEA": ["G", "C", "C"],
            "SNPID": ["1:100:A:G", "1:200:G:C", "2:300:T:C"],
            "BETA": [0.12, -0.18, 0.28],
            "SE": [0.06, 0.11, 0.16],
            "P": [0.04, 0.10, 0.08],
            "N": [2000, 2000, 2000],
            "EAF": [0.32, 0.38, 0.52],
        })
        
        self.s1 = make_sumstatsp(df1, study_name="Study1")
        self.s2 = make_sumstatsp(df2, study_name="Study2")

    def test_meta_analyze_polars_with_sumstatsmulti(self):
        """Test meta_analyze_polars through SumstatsMulti.run_meta_analysis()"""
        # Create SumstatsMulti with polars engine
        sm = SumstatsMulti(
            sumstatsObjects=[self.s1, self.s2],
            engine="polars",
            verbose=False
        )
        
        # Run meta-analysis
        result = sm.run_meta_analysis(random_effects=False)
        
        # Check that result is a polars DataFrame
        self.assertIsInstance(result, pl.DataFrame)
        
        # Check that it has data
        self.assertGreater(result.height, 0)
        
        # Check that expected columns are present
        expected_cols = ["CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF", "Z", "Q", "DOF", "DIRECTION"]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Column {col} not found in result")
        
        # Check that meta-analysis columns are computed
        self.assertIn("BETA", result.columns)
        self.assertIn("SE", result.columns)
        self.assertIn("P", result.columns)
        self.assertIn("Z", result.columns)
        
        # Check that study-specific columns are removed (should not have _1, _2 suffixes)
        for col in result.columns:
            self.assertFalse(col.endswith("_1") or col.endswith("_2"), 
                           f"Study-specific column {col} should be removed")

    def test_meta_analyze_polars_with_random_effects(self):
        """Test meta_analyze_polars with random effects model"""
        sm = SumstatsMulti(
            sumstatsObjects=[self.s1, self.s2],
            engine="polars",
            verbose=False
        )
        
        # Run meta-analysis with random effects
        result = sm.run_meta_analysis(random_effects=True)
        
        # Check that result is a polars DataFrame
        self.assertIsInstance(result, pl.DataFrame)
        
        # Check that random effects columns are present
        expected_cols = ["BETA_RANDOM", "SE_RANDOM", "Z_RANDOM", "P_RANDOM"]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Random effects column {col} not found")
        
        # Check that fixed effects columns are also present
        self.assertIn("BETA", result.columns)
        self.assertIn("SE", result.columns)

    def test_meta_analyze_polars_direct_call(self):
        """Test meta_analyze_polars function directly with polars DataFrame"""
        # Create SumstatsMulti to get the merged data structure
        sm = SumstatsMulti(
            sumstatsObjects=[self.s1, self.s2],
            engine="polars",
            verbose=False
        )
        
        # Call meta_analyze_polars directly (note: function doesn't accept verbose parameter)
        result = meta_analyze_polars(
            sumstats_multi=sm.data,
            random_effects=False,
            nstudy=2,
            log=sm.log
        )
        
        # Check that result is a polars DataFrame
        self.assertIsInstance(result, pl.DataFrame)
        
        # Check that it has data
        self.assertGreater(result.height, 0)
        
        # Check that expected columns are present
        self.assertIn("BETA", result.columns)
        self.assertIn("SE", result.columns)
        self.assertIn("P", result.columns)
        self.assertIn("Z", result.columns)
        self.assertIn("N", result.columns)
        self.assertIn("EAF", result.columns)

    def test_meta_analyze_polars_with_three_studies(self):
        """Test meta_analyze_polars with three studies"""
        # Create third study
        df3 = pd.DataFrame({
            "CHR": [1, 2],
            "POS": [100, 300],
            "EA": ["A", "T"],
            "NEA": ["G", "C"],
            "SNPID": ["1:100:A:G", "2:300:T:C"],
            "BETA": [0.11, 0.29],
            "SE": [0.07, 0.17],
            "P": [0.12, 0.09],
            "N": [1500, 1500],
            "EAF": [0.31, 0.51],
        })
        
        s3 = make_sumstatsp(df3, study_name="Study3")
        
        sm = SumstatsMulti(
            sumstatsObjects=[self.s1, self.s2, s3],
            engine="polars",
            verbose=False
        )
        
        # Run meta-analysis
        result = sm.run_meta_analysis(random_effects=False)
        
        # Check that result is a polars DataFrame
        self.assertIsInstance(result, pl.DataFrame)
        
        # Check that it has data
        self.assertGreater(result.height, 0)
        
        # Check that meta-analysis columns are computed
        self.assertIn("BETA", result.columns)
        self.assertIn("SE", result.columns)
        self.assertIn("P", result.columns)
        self.assertIn("DOF", result.columns)
        
        # Check that DOF is correct (should be number of studies with non-null BETA)
        # For overlapping variants, DOF should be 3
        if result.height > 0:
            dof_values = result["DOF"].to_list()
            # At least some variants should have DOF >= 2
            self.assertTrue(any(d >= 2 for d in dof_values if d is not None))

    def test_meta_analyze_polars_direction_column(self):
        """Test that DIRECTION column is correctly computed"""
        sm = SumstatsMulti(
            sumstatsObjects=[self.s1, self.s2],
            engine="polars",
            verbose=False
        )
        
        result = sm.run_meta_analysis(random_effects=False)
        
        # Check that DIRECTION column exists
        self.assertIn("DIRECTION", result.columns)
        
        # Check that DIRECTION values are strings
        if result.height > 0:
            direction_values = result["DIRECTION"].to_list()
            self.assertTrue(all(isinstance(d, str) for d in direction_values if d is not None))

    def test_meta_analyze_polars_statistics_range(self):
        """Test that computed statistics are in reasonable ranges"""
        sm = SumstatsMulti(
            sumstatsObjects=[self.s1, self.s2],
            engine="polars",
            verbose=False
        )
        
        result = sm.run_meta_analysis(random_effects=False)
        
        if result.height > 0:
            # Check that P-values are between 0 and 1
            p_values = result["P"].to_list()
            self.assertTrue(all(0 <= p <= 1 for p in p_values if p is not None))
            
            # Check that SE is positive
            se_values = result["SE"].to_list()
            self.assertTrue(all(se > 0 for se in se_values if se is not None))
            
            # Check that EAF is between 0 and 1
            eaf_values = result["EAF"].to_list()
            self.assertTrue(all(0 <= eaf <= 1 for eaf in eaf_values if eaf is not None))


if __name__ == "__main__":
    unittest.main()

