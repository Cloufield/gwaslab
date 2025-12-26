"""
Test suite for GSF (GWASLab Standard Format) I/O functions.
Tests load_gsf() and to_gsf() methods.
"""
import os
import sys
import unittest
import tempfile
import shutil
from pathlib import Path

import pandas as pd
import numpy as np

# Add parent directory to path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl
from gwaslab.g_Sumstats import Sumstats


def make_test_sumstats(n_variants=100):
    """Create test sumstats DataFrame."""
    np.random.seed(42)
    rows = []
    for i in range(n_variants):
        chr_ = np.random.choice([1, 2, 3, 7, 22])
        pos_ = 1000000 + i * 1000
        snpid = f"{chr_}:{pos_}_A_G"
        rows.append({
            "CHR": chr_,
            "POS": pos_,
            "EA": "A",
            "NEA": "G",
            "SNPID": snpid,
            "rsID": f"rs{i+1}",
            "EAF": np.random.uniform(0.1, 0.9),
            "BETA": np.random.normal(0, 0.1),
            "SE": np.random.uniform(0.01, 0.1),
            "P": np.random.uniform(1e-8, 1.0),
            "N": np.random.randint(1000, 100000),
        })
    return pd.DataFrame(rows)


class TestGSFIO(unittest.TestCase):
    """Test GSF I/O functionality."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests."""
        cls.test_dir = Path(__file__).parent
        cls.temp_dir = cls.test_dir / "tmp" / f"gwaslab_gsf_test_{os.getpid()}"
        cls.temp_dir.mkdir(parents=True, exist_ok=True)
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        if cls.temp_dir.exists():
            shutil.rmtree(cls.temp_dir)
    
    def setUp(self):
        """Set up for each test."""
        # Create test sumstats
        df = make_test_sumstats(n_variants=100)
        self.sumstats = Sumstats(sumstats=df, verbose=False)
        self.sumstats.set_build("19", verbose=False)
        # Add some metadata
        self.sumstats.meta["test_meta"] = {"key": "value", "number": 42}
    
    def tearDown(self):
        """Clean up after each test."""
        # Remove any test files
        for file in self.temp_dir.glob("*"):
            if file.is_file():
                file.unlink()
            elif file.is_dir():
                shutil.rmtree(file)
    
    def test_to_gsf_basic(self):
        """Test basic to_gsf() functionality."""
        gsf_path = self.temp_dir / "test_basic.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Check file exists
        self.assertTrue(gsf_path.exists(), "GSF file should be created")
        self.assertGreater(gsf_path.stat().st_size, 0, "GSF file should not be empty")
    
    def test_load_gsf_basic(self):
        """Test basic load_gsf() functionality."""
        gsf_path = self.temp_dir / "test_load_basic.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load back
        loaded = gl.load_gsf(str(gsf_path), verbose=False)
        
        # Check data integrity
        self.assertEqual(len(loaded.data), len(self.sumstats.data), "Number of variants should match")
        self.assertListEqual(list(loaded.data.columns), list(self.sumstats.data.columns), "Columns should match")
        
        # Compare values (dtypes may differ: categories become strings in GSF)
        for col in self.sumstats.data.columns:
            if isinstance(self.sumstats.data[col].dtype, pd.CategoricalDtype):
                # Convert category to string for comparison
                original_values = self.sumstats.data[col].astype(str).values
            else:
                original_values = self.sumstats.data[col].values
            
            if pd.api.types.is_string_dtype(loaded.data[col]):
                loaded_values = loaded.data[col].astype(str).values
            else:
                loaded_values = loaded.data[col].values
            
            np.testing.assert_array_equal(loaded_values, original_values, 
                                        err_msg=f"Column {col} values should match")
    
    def test_metadata_preservation(self):
        """Test that metadata is preserved in GSF files."""
        gsf_path = self.temp_dir / "test_metadata.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load back
        loaded = gl.load_gsf(str(gsf_path), verbose=False)
        
        # Check metadata
        self.assertIn("test_meta", loaded.meta, "Custom metadata should be preserved")
        self.assertEqual(loaded.meta["test_meta"]["key"], "value", "Metadata values should match")
        self.assertEqual(loaded.meta["test_meta"]["number"], 42, "Metadata numbers should match")
        self.assertEqual(loaded._build, "19", "Genome build should be preserved")
    
    def test_load_gsf_column_selection(self):
        """Test loading specific columns."""
        gsf_path = self.temp_dir / "test_columns.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load only specific columns
        columns = ["CHR", "POS", "P", "BETA"]
        loaded = gl.load_gsf(str(gsf_path), columns=columns, verbose=False)
        
        # Check columns
        self.assertEqual(set(loaded.data.columns), set(columns), "Only selected columns should be loaded")
        self.assertEqual(len(loaded.data), len(self.sumstats.data), "Number of variants should match")
    
    def test_load_gsf_filter_equality(self):
        """Test filtering with equality operator."""
        gsf_path = self.temp_dir / "test_filter_eq.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Filter for CHR == 7
        loaded = gl.load_gsf(str(gsf_path), filters="CHR == 7", verbose=False)
        
        # Check filter worked
        self.assertTrue(all(loaded.data["CHR"] == 7), "All variants should have CHR == 7")
        expected_count = len(self.sumstats.data[self.sumstats.data["CHR"] == 7])
        self.assertEqual(len(loaded.data), expected_count, "Filtered count should match")
    
    def test_load_gsf_filter_inequality(self):
        """Test filtering with inequality operators."""
        gsf_path = self.temp_dir / "test_filter_ineq.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Filter for P < 0.05
        loaded = gl.load_gsf(str(gsf_path), filters="P < 0.05", verbose=False)
        
        # Check filter worked
        self.assertTrue(all(loaded.data["P"] < 0.05), "All variants should have P < 0.05")
        expected_count = len(self.sumstats.data[self.sumstats.data["P"] < 0.05])
        self.assertEqual(len(loaded.data), expected_count, "Filtered count should match")
    
    def test_load_gsf_filter_and(self):
        """Test filtering with AND operator."""
        gsf_path = self.temp_dir / "test_filter_and.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Filter for CHR == 7 & P < 0.05
        loaded = gl.load_gsf(str(gsf_path), filters="CHR == 7 & P < 0.05", verbose=False)
        
        # Check filter worked
        self.assertTrue(all(loaded.data["CHR"] == 7), "All variants should have CHR == 7")
        self.assertTrue(all(loaded.data["P"] < 0.05), "All variants should have P < 0.05")
        expected = self.sumstats.data[(self.sumstats.data["CHR"] == 7) & (self.sumstats.data["P"] < 0.05)]
        self.assertEqual(len(loaded.data), len(expected), "Filtered count should match")
    
    def test_load_gsf_filter_in(self):
        """Test filtering with 'in' operator."""
        gsf_path = self.temp_dir / "test_filter_in.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Filter for CHR in [1, 2, 3]
        loaded = gl.load_gsf(str(gsf_path), filters="CHR in [1, 2, 3]", verbose=False)
        
        # Check filter worked
        self.assertTrue(all(loaded.data["CHR"].isin([1, 2, 3])), "All variants should have CHR in [1, 2, 3]")
        expected = self.sumstats.data[self.sumstats.data["CHR"].isin([1, 2, 3])]
        self.assertEqual(len(loaded.data), len(expected), "Filtered count should match")
    
    def test_load_gsf_filter_combined(self):
        """Test filtering with column selection and filters."""
        gsf_path = self.temp_dir / "test_filter_combined.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load with both column selection and filter
        columns = ["CHR", "POS", "P", "BETA"]
        loaded = gl.load_gsf(str(gsf_path), columns=columns, filters="P < 0.01", verbose=False)
        
        # Check both worked
        self.assertEqual(set(loaded.data.columns), set(columns), "Only selected columns should be loaded")
        self.assertTrue(all(loaded.data["P"] < 0.01), "All variants should have P < 0.01")
    
    def test_to_gsf_partitioned(self):
        """Test saving partitioned GSF files."""
        gsf_dir = self.temp_dir / "test_partitioned"
        gsf_dir.mkdir(exist_ok=True)
        
        # Save partitioned by CHR
        self.sumstats.to_gsf(str(gsf_dir), partition_cols=["CHR"], verbose=False)
        
        # Check partition directories exist
        self.assertTrue(gsf_dir.exists(), "Partition directory should exist")
        # Check that some partition files exist
        parquet_files = list(gsf_dir.rglob("*.parquet"))
        self.assertGreater(len(parquet_files), 0, "Partition files should be created")
    
    def test_load_gsf_partitioned(self):
        """Test loading partitioned GSF files."""
        gsf_dir = self.temp_dir / "test_partitioned_load"
        gsf_dir.mkdir(exist_ok=True)
        
        # Save partitioned by CHR
        self.sumstats.to_gsf(str(gsf_dir), partition_cols=["CHR"], verbose=False)
        
        # Load back
        loaded = gl.load_gsf(str(gsf_dir), verbose=False)
        
        # Check data integrity
        self.assertEqual(len(loaded.data), len(self.sumstats.data), "Number of variants should match")
        
        # Sort both dataframes by SNPID for comparison (order may differ after partitioning)
        # Use SNPID as a unique identifier
        if "SNPID" in loaded.data.columns and "SNPID" in self.sumstats.data.columns:
            loaded_sorted = loaded.data.sort_values(by="SNPID").reset_index(drop=True)
            original_sorted = self.sumstats.data.sort_values(by="SNPID").reset_index(drop=True)
            
            # Compare column by column (handling dtype differences)
            for col in original_sorted.columns:
                if col in loaded_sorted.columns:
                    if isinstance(original_sorted[col].dtype, pd.CategoricalDtype):
                        original_values = original_sorted[col].astype(str).values
                    else:
                        original_values = original_sorted[col].values
                    
                    if pd.api.types.is_string_dtype(loaded_sorted[col]):
                        loaded_values = loaded_sorted[col].astype(str).values
                    else:
                        loaded_values = loaded_sorted[col].values
                    
                    np.testing.assert_array_equal(loaded_values, original_values,
                                                err_msg=f"Column {col} should match in partitioned data")
        else:
            # Fallback: just check that we got the right number of rows
            self.assertEqual(len(loaded.data), len(self.sumstats.data), 
                          "Partitioned data should have same number of rows")
    
    def test_to_gsf_compression(self):
        """Test different compression options."""
        compressions = ["zstd", "snappy", "gzip"]
        
        for comp in compressions:
            gsf_path = self.temp_dir / f"test_compression_{comp}.gsf"
            self.sumstats.to_gsf(str(gsf_path), compression=comp, verbose=False)
            
            # Check file exists
            self.assertTrue(gsf_path.exists(), f"GSF file with {comp} compression should be created")
            
            # Load back and verify
            loaded = gl.load_gsf(str(gsf_path), verbose=False)
            self.assertEqual(len(loaded.data), len(self.sumstats.data), 
                           f"Data integrity should be preserved with {comp} compression")
    
    def test_load_gsf_filter_type_casting(self):
        """Test that filters work even when columns are stored as strings."""
        gsf_path = self.temp_dir / "test_type_cast.gsf"
        
        # Create sumstats with CHR as category (will be converted to string in GSF)
        df = self.sumstats.data.copy()
        df["CHR"] = df["CHR"].astype("category")
        sumstats_cat = Sumstats(sumstats=df, verbose=False)
        sumstats_cat.to_gsf(str(gsf_path), verbose=False)
        
        # Try to filter with numeric comparison (should handle type casting)
        loaded = gl.load_gsf(str(gsf_path), filters="CHR == 7", verbose=False)
        
        # Check filter worked (even though CHR was stored as string)
        if len(loaded.data) > 0:
            self.assertTrue(all(loaded.data["CHR"] == 7), "Filter should work with type casting")
    
    def test_load_gsf_empty_filter(self):
        """Test loading without filters returns all data."""
        gsf_path = self.temp_dir / "test_no_filter.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load without filter
        loaded = gl.load_gsf(str(gsf_path), verbose=False)
        
        # Should have all data
        self.assertEqual(len(loaded.data), len(self.sumstats.data), "Should load all data without filter")
    
    def test_load_gsf_nonexistent_file(self):
        """Test that loading non-existent file raises error."""
        gsf_path = self.temp_dir / "nonexistent.gsf"
        
        with self.assertRaises(FileNotFoundError):
            gl.load_gsf(str(gsf_path), verbose=False)
    
    def test_load_gsf_invalid_filter(self):
        """Test that invalid filter syntax raises error."""
        gsf_path = self.temp_dir / "test_invalid_filter.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Invalid filter syntax
        with self.assertRaises((ValueError, IOError)):
            gl.load_gsf(str(gsf_path), filters="invalid filter syntax", verbose=False)
    
    def test_to_gsf_large_dataset(self):
        """Test saving and loading larger dataset."""
        # Create larger dataset
        df = make_test_sumstats(n_variants=1000)
        large_sumstats = Sumstats(sumstats=df, verbose=False)
        
        gsf_path = self.temp_dir / "test_large.gsf"
        large_sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load back
        loaded = gl.load_gsf(str(gsf_path), verbose=False)
        
        # Check integrity
        self.assertEqual(len(loaded.data), len(large_sumstats.data), "Large dataset should load correctly")
    
    def test_to_gsf_preserves_dtypes(self):
        """Test that data types are preserved appropriately."""
        gsf_path = self.temp_dir / "test_dtypes.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Load back
        loaded = gl.load_gsf(str(gsf_path), verbose=False)
        
        # Check that numeric columns are numeric
        self.assertTrue(pd.api.types.is_numeric_dtype(loaded.data["CHR"]), "CHR should be numeric")
        self.assertTrue(pd.api.types.is_numeric_dtype(loaded.data["POS"]), "POS should be numeric")
        self.assertTrue(pd.api.types.is_numeric_dtype(loaded.data["P"]), "P should be numeric")
        self.assertTrue(pd.api.types.is_numeric_dtype(loaded.data["BETA"]), "BETA should be numeric")
    
    def test_load_gsf_filter_multiple_conditions(self):
        """Test filtering with multiple conditions."""
        gsf_path = self.temp_dir / "test_multi_filter.gsf"
        self.sumstats.to_gsf(str(gsf_path), verbose=False)
        
        # Filter with multiple conditions
        loaded = gl.load_gsf(str(gsf_path), filters="CHR == 7 & P < 0.01 & BETA > 0", verbose=False)
        
        # Check all conditions met
        if len(loaded.data) > 0:
            self.assertTrue(all(loaded.data["CHR"] == 7), "CHR condition should be met")
            self.assertTrue(all(loaded.data["P"] < 0.01), "P condition should be met")
            self.assertTrue(all(loaded.data["BETA"] > 0), "BETA condition should be met")


if __name__ == "__main__":
    unittest.main()

