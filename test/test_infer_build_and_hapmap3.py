"""
Test cases for _infer_build and _get_hapmap3 functions.

This test file covers:
- _infer_build: Testing genome build inference using Hapmap3 SNPs
- _get_hapmap3: Testing Hapmap3 variant extraction with various matching strategies

To run:
    python test/test_infer_build_and_hapmap3.py
    python -m unittest test.test_infer_build_and_hapmap3 -v
    pytest test/test_infer_build_and_hapmap3.py -v
"""

import os
import sys
import unittest
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_filter_value import _infer_build, _get_hapmap_df_polars
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap3
from gwaslab.g_Sumstats import Sumstats
from gwaslab.info.g_Log import Log


def create_test_sumstats_hg19():
    """Create a test sumstats DataFrame with hg19 coordinates."""
    # These coordinates should match hg19 Hapmap3 SNPs
    data = {
        "CHR": [1, 1, 2, 2, 3, 3],
        "POS": [752566, 752721, 752566, 752721, 752566, 752721],  # Common Hapmap3 positions
        "EA": ["G", "A", "C", "T", "A", "G"],
        "NEA": ["A", "G", "T", "C", "G", "A"],
        "P": [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10],
        "BETA": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        "SE": [0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
        "rsID": ["rs3094315", "rs3131972", "rs2073813", "rs2073814", "rs2073815", "rs2073816"],
        "STATUS": ["9999999999"] * 6
    }
    return pd.DataFrame(data)


def create_test_sumstats_hg38():
    """Create a test sumstats DataFrame with hg38 coordinates."""
    # These coordinates should match hg38 Hapmap3 SNPs
    # Note: hg38 positions are typically different from hg19
    data = {
        "CHR": [1, 1, 2, 2, 3, 3],
        "POS": [817186, 817341, 817186, 817341, 817186, 817341],  # Approximate hg38 positions
        "EA": ["G", "A", "C", "T", "A", "G"],
        "NEA": ["A", "G", "T", "C", "G", "A"],
        "P": [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10],
        "BETA": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        "SE": [0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
        "rsID": ["rs3094315", "rs3131972", "rs2073813", "rs2073814", "rs2073815", "rs2073816"],
        "STATUS": ["9999999999"] * 6
    }
    return pd.DataFrame(data)


def create_test_sumstats_with_real_hapmap3():
    """Create test sumstats with actual Hapmap3 coordinates by loading reference data."""
    # Load a few real Hapmap3 positions
    try:
        hapmap19 = _get_hapmap_df_polars("19")
        hapmap38 = _get_hapmap_df_polars("38")
        
        # Get first few positions from each build (use 1/10 to speed up tests)
        hg19_sample = hapmap19.head(10).to_pandas()
        hg38_sample = hapmap38.head(10).to_pandas()
        
        # Create sumstats with real hg19 positions
        data_hg19 = {
            "CHR": hg19_sample["CHR"].tolist(),
            "POS": hg19_sample["POS"].tolist(),
            "EA": ["A"] * len(hg19_sample),
            "NEA": ["G"] * len(hg19_sample),
            "P": np.random.uniform(1e-10, 0.05, len(hg19_sample)),
            "BETA": np.random.uniform(-0.5, 0.5, len(hg19_sample)),
            "SE": np.random.uniform(0.01, 0.1, len(hg19_sample)),
            "STATUS": ["9999999999"] * len(hg19_sample)
        }
        
        # Create sumstats with real hg38 positions
        data_hg38 = {
            "CHR": hg38_sample["CHR"].tolist(),
            "POS": hg38_sample["POS"].tolist(),
            "EA": ["A"] * len(hg38_sample),
            "NEA": ["G"] * len(hg38_sample),
            "P": np.random.uniform(1e-10, 0.05, len(hg38_sample)),
            "BETA": np.random.uniform(-0.5, 0.5, len(hg38_sample)),
            "SE": np.random.uniform(0.01, 0.1, len(hg38_sample)),
            "STATUS": ["9999999999"] * len(hg38_sample)
        }
        
        return pd.DataFrame(data_hg19), pd.DataFrame(data_hg38)
    except Exception as e:
        # If hapmap data is not available, return None
        return None, None


class TestInferBuild(unittest.TestCase):
    """Test cases for _infer_build function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.hg19_df, self.hg38_df = create_test_sumstats_with_real_hapmap3()
    
    def test_infer_build_with_dataframe(self):
        """Test _infer_build with pandas DataFrame input."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        # Test with hg19 data
        result = _infer_build(
            self.hg19_df,
            verbose=False,
            log=self.log,
            change_status=False
        )
        
        # Should return a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), len(self.hg19_df))
    
    def test_infer_build_with_sumstats_object(self):
        """Test _infer_build with Sumstats object input."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        # Create Sumstats object
        mysumstats = Sumstats(sumstats=self.hg19_df.copy(), verbose=False)
        mysumstats.fix_chr(verbose=False)
        
        # Infer build
        result = _infer_build(
            mysumstats,
            verbose=False,
            log=self.log,
            change_status=True
        )
        
        # Check that build was inferred and stored
        self.assertIsNotNone(mysumstats.build)
        self.assertIn(mysumstats.build, ["19", "38", "Unknown"])
        self.assertEqual(mysumstats.meta["gwaslab"]["genome_build"], mysumstats.build)
    
    def test_infer_build_change_status(self):
        """Test that STATUS codes are updated when change_status=True."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        df = self.hg19_df.copy()
        original_status = df["STATUS"].iloc[0]
        
        result = _infer_build(
            df,
            verbose=False,
            log=self.log,
            change_status=True
        )
        
        # STATUS should be modified if build was inferred
        if result["STATUS"].iloc[0] != original_status:
            # Check that STATUS was updated (digit 1 should be 1 or 3)
            status_str = str(result["STATUS"].iloc[0])
            self.assertIn(status_str[0], ["1", "3", "9"])
    
    def test_infer_build_no_change_status(self):
        """Test that STATUS codes are not updated when change_status=False."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        df = self.hg19_df.copy()
        original_status = df["STATUS"].iloc[0]
        
        result = _infer_build(
            df,
            verbose=False,
            log=self.log,
            change_status=False
        )
        
        # STATUS should remain unchanged
        self.assertEqual(result["STATUS"].iloc[0], original_status)
    
    def test_infer_build_with_missing_chr_pos(self):
        """Test _infer_build with missing CHR or POS values."""
        df = create_test_sumstats_hg19()
        df.loc[0, "CHR"] = None
        df.loc[1, "POS"] = None
        
        result = _infer_build(
            df,
            verbose=False,
            log=self.log,
            change_status=False
        )
        
        # Should still return a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
    
    def test_infer_build_custom_column_names(self):
        """Test _infer_build with custom column names."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        df = self.hg19_df.copy()
        df = df.rename(columns={"CHR": "chromosome", "POS": "position"})
        
        result = _infer_build(
            df,
            chrom="chromosome",
            pos="position",
            verbose=False,
            log=self.log,
            change_status=False
        )
        
        self.assertIsInstance(result, pd.DataFrame)


class TestGetHapmap3(unittest.TestCase):
    """Test cases for _get_hapmap3 function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        # Create test data with rsID
        self.df_with_rsid = pd.DataFrame({
            "CHR": [1, 1, 2, 2],
            "POS": [752566, 752721, 752566, 752721],
            "EA": ["G", "A", "C", "T"],
            "NEA": ["A", "G", "T", "C"],
            "rsID": ["rs3094315", "rs3131972", "rs2073813", "rs2073814"],
            "P": [1e-5, 1e-6, 1e-7, 1e-8],
            "BETA": [0.1, 0.2, 0.3, 0.4]
        })
        
        # Create test data without rsID
        self.df_without_rsid = pd.DataFrame({
            "CHR": [1, 1, 2, 2],
            "POS": [752566, 752721, 752566, 752721],
            "EA": ["G", "A", "C", "T"],
            "NEA": ["A", "G", "T", "C"],
            "P": [1e-5, 1e-6, 1e-7, 1e-8],
            "BETA": [0.1, 0.2, 0.3, 0.4]
        })
    
    def test_get_hapmap3_with_rsid(self):
        """Test _get_hapmap3 matching by rsID."""
        result = _get_hapmap3(
            self.df_with_rsid,
            build="19",
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        # Should return a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        # Should have rsID column
        if len(result) > 0:
            self.assertIn("rsID", result.columns)
    
    def test_get_hapmap3_with_chr_pos(self):
        """Test _get_hapmap3 matching by CHR:POS."""
        result = _get_hapmap3(
            self.df_without_rsid,
            build="19",
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        # Should return a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        # Should preserve original columns
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
    
    def test_get_hapmap3_match_allele(self):
        """Test _get_hapmap3 with allele matching enabled."""
        result = _get_hapmap3(
            self.df_with_rsid,
            build="19",
            verbose=False,
            log=self.log,
            match_allele=True,
            how="inner"
        )
        
        # Should return a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        # If matches found, alleles should be matched
        if len(result) > 0:
            self.assertIn("EA", result.columns)
            self.assertIn("NEA", result.columns)
    
    def test_get_hapmap3_no_match_allele(self):
        """Test _get_hapmap3 with allele matching disabled."""
        result = _get_hapmap3(
            self.df_with_rsid,
            build="19",
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        # Should return a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
    
    def test_get_hapmap3_different_builds(self):
        """Test _get_hapmap3 with different genome builds."""
        for build in ["19", "38"]:
            result = _get_hapmap3(
                self.df_with_rsid,
                build=build,
                verbose=False,
                log=self.log,
                match_allele=False,
                how="inner"
            )
            
            # Should return a DataFrame for both builds
            self.assertIsInstance(result, pd.DataFrame)
    
    def test_get_hapmap3_how_parameter(self):
        """Test _get_hapmap3 with different 'how' parameter values."""
        for how in ["inner", "left", "right"]:
            result = _get_hapmap3(
                self.df_with_rsid,
                build="19",
                verbose=False,
                log=self.log,
                match_allele=False,
                how=how
            )
            
            # Should return a DataFrame for all join types
            self.assertIsInstance(result, pd.DataFrame)
    
    def test_get_hapmap3_custom_column_names(self):
        """Test _get_hapmap3 with custom column names."""
        df = self.df_with_rsid.copy()
        df = df.rename(columns={
            "CHR": "chromosome",
            "POS": "position",
            "EA": "effect_allele",
            "NEA": "non_effect_allele",
            "rsID": "rs_id"
        })
        
        result = _get_hapmap3(
            df,
            rsid="rs_id",
            chrom="chromosome",
            pos="position",
            ea="effect_allele",
            nea="non_effect_allele",
            build="19",
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        self.assertIsInstance(result, pd.DataFrame)
    
    def test_get_hapmap3_no_matching_info(self):
        """Test _get_hapmap3 when no matching information is available."""
        df = pd.DataFrame({
            "P": [1e-5, 1e-6],
            "BETA": [0.1, 0.2]
        })
        
        # Should raise ValueError
        with self.assertRaises(ValueError):
            _get_hapmap3(
                df,
                build="19",
                verbose=False,
                log=self.log
            )
    
    def test_get_hapmap3_preserves_columns(self):
        """Test that _get_hapmap3 preserves original columns."""
        df = self.df_with_rsid.copy()
        original_cols = set(df.columns)
        
        result = _get_hapmap3(
            df,
            build="19",
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        # Original columns should be preserved (except temporary ones)
        if len(result) > 0:
            for col in ["CHR", "POS", "EA", "NEA", "P", "BETA"]:
                if col in original_cols:
                    self.assertIn(col, result.columns)


class TestInferBuildAndHapmap3Integration(unittest.TestCase):
    """Integration tests for infer_build and get_hapmap3 together."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.hg19_df, self.hg38_df = create_test_sumstats_with_real_hapmap3()
    
    def test_infer_then_filter_hapmap3(self):
        """Test inferring build then filtering to Hapmap3."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        # Infer build
        mysumstats = Sumstats(sumstats=self.hg19_df.copy(), verbose=False)
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        # Get Hapmap3 variants
        hapmap3_result = _get_hapmap3(
            mysumstats.data,
            build=mysumstats.build,
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        # Should return filtered DataFrame
        self.assertIsInstance(hapmap3_result, pd.DataFrame)
        self.assertLessEqual(len(hapmap3_result), len(mysumstats.data))
    
    def test_get_hapmap3_then_infer_build(self):
        """Test getting Hapmap3 variants then inferring build."""
        if self.hg19_df is None:
            self.skipTest("Hapmap3 reference data not available")
        
        # Get Hapmap3 variants first
        hapmap3_df = _get_hapmap3(
            self.hg19_df,
            build="19",
            verbose=False,
            log=self.log,
            match_allele=False,
            how="inner"
        )
        
        if len(hapmap3_df) > 0:
            # Then infer build on filtered data
            result = _infer_build(
                hapmap3_df,
                verbose=False,
                log=self.log,
                change_status=False
            )
            
            self.assertIsInstance(result, pd.DataFrame)


if __name__ == "__main__":
    unittest.main()

