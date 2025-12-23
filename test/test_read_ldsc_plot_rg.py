import os
import sys
import unittest
import pandas as pd
import numpy as np
import tempfile
import shutil

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.io.io_read_ldsc import read_ldsc
from gwaslab.viz.viz_plot_rg_heatmap import plot_rg


class TestReadLDSCAndPlotRG(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.test_ref_dir = os.path.join(ROOT, "test", "ref", "ldsc_logs")
        
        # Paths to simulated LDSC log files
        self.h2_log_1 = os.path.join(self.test_ref_dir, "simulated_h2_1.log")
        self.h2_log_2 = os.path.join(self.test_ref_dir, "simulated_h2_2.log")
        self.rg_log = os.path.join(self.test_ref_dir, "simulated_rg.log")
        self.rg_log_multiple = os.path.join(self.test_ref_dir, "simulated_rg_multiple.log")

    def test_read_ldsc_h2_mode_single_file(self):
        """Test reading single LDSC h2 log file"""
        result = read_ldsc([self.h2_log_1], mode="h2")
        
        # Check result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame, "Result should be a DataFrame")
        
        # Check expected columns
        expected_columns = ['Filename', 'h2_obs', 'h2_se', 'Lambda_gc', 'Mean_chi2', 
                           'Intercept', 'Intercept_se', 'Ratio', 'Ratio_se']
        self.assertEqual(list(result.columns), expected_columns, 
                        "Result should have expected columns")
        
        # Check single row
        self.assertEqual(len(result), 1, "Should have one row")
        
        # Check values
        self.assertEqual(result.iloc[0]['Filename'], 'simulated_h2_1.log')
        self.assertEqual(result.iloc[0]['h2_obs'], '0.1583')
        self.assertEqual(result.iloc[0]['h2_se'], '0.0281')
        self.assertEqual(result.iloc[0]['Lambda_gc'], '1.1523')
        self.assertEqual(result.iloc[0]['Mean_chi2'], '1.2843')
        self.assertEqual(result.iloc[0]['Intercept'], '1.0563')
        self.assertEqual(result.iloc[0]['Intercept_se'], '0.0114')
        self.assertEqual(result.iloc[0]['Ratio'], '0.1981')
        self.assertEqual(result.iloc[0]['Ratio_se'], '0.0402')

    def test_read_ldsc_h2_mode_multiple_files(self):
        """Test reading multiple LDSC h2 log files"""
        result = read_ldsc([self.h2_log_1, self.h2_log_2], mode="h2")
        
        # Check result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame, "Result should be a DataFrame")
        
        # Check two rows
        self.assertEqual(len(result), 2, "Should have two rows")
        
        # Check filenames
        filenames = result['Filename'].tolist()
        self.assertIn('simulated_h2_1.log', filenames)
        self.assertIn('simulated_h2_2.log', filenames)
        
        # Check values from second file
        row2 = result[result['Filename'] == 'simulated_h2_2.log'].iloc[0]
        self.assertEqual(row2['h2_obs'], '0.0743')
        self.assertEqual(row2['h2_se'], '0.0123')

    def test_read_ldsc_rg_mode_single_pair(self):
        """Test reading LDSC rg log file with single trait pair"""
        result = read_ldsc([self.rg_log], mode="rg")
        
        # Check result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame, "Result should be a DataFrame")
        
        # Check expected columns
        expected_columns = ['p1', 'p2', 'rg', 'se', 'z', 'p', 'h2_obs', 'h2_obs_se',
                           'h2_int', 'h2_int_se', 'gcov_int', 'gcov_int_se']
        self.assertEqual(list(result.columns), expected_columns,
                        "Result should have expected columns")
        
        # Check single row
        self.assertEqual(len(result), 1, "Should have one row")
        
        # Check values
        self.assertEqual(result.iloc[0]['p1'], 'trait1.sumstats.gz')
        self.assertEqual(result.iloc[0]['p2'], 'trait2.sumstats.gz')
        self.assertAlmostEqual(float(result.iloc[0]['rg']), 0.1601, places=4)
        self.assertAlmostEqual(float(result.iloc[0]['se']), 0.1821, places=4)
        self.assertAlmostEqual(float(result.iloc[0]['p']), 0.3792, places=4)

    def test_read_ldsc_rg_mode_multiple_pairs(self):
        """Test reading LDSC rg log file with multiple trait pairs"""
        result = read_ldsc([self.rg_log_multiple], mode="rg")
        
        # Check result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame, "Result should be a DataFrame")
        
        # Check three rows
        self.assertEqual(len(result), 3, "Should have three rows")
        
        # Check that numeric columns are properly converted
        self.assertTrue(pd.api.types.is_numeric_dtype(result['rg']))
        self.assertTrue(pd.api.types.is_numeric_dtype(result['p']))
        
        # Check values from first pair
        row1 = result[(result['p1'] == 'trait1.sumstats.gz') & 
                     (result['p2'] == 'trait2.sumstats.gz')].iloc[0]
        self.assertAlmostEqual(float(row1['rg']), 0.1601, places=4)
        
        # Check values from second pair
        row2 = result[(result['p1'] == 'trait1.sumstats.gz') & 
                     (result['p2'] == 'trait3.sumstats.gz')].iloc[0]
        self.assertAlmostEqual(float(row2['rg']), 0.4523, places=4)
        self.assertAlmostEqual(float(row2['p']), 0.0002, places=4)

    def test_plot_rg_basic(self):
        """Test basic plot_rg functionality"""
        # Read rg data
        rg_data = read_ldsc([self.rg_log_multiple], mode="rg")
        
        # Plot should not raise errors
        try:
            fig, ax, log, df = plot_rg(rg_data, verbose=False)
            success = True
        except Exception as e:
            success = False
            self.fail(f"plot_rg raised an exception: {e}")
        
        self.assertTrue(success, "plot_rg should complete without errors")
        self.assertIsNotNone(fig, "plot_rg should return a figure")
        self.assertIsNotNone(ax, "plot_rg should return an axes")
        self.assertIsNotNone(df, "plot_rg should return a DataFrame")

    def test_plot_rg_with_custom_parameters(self):
        """Test plot_rg with custom parameters"""
        # Read rg data
        rg_data = read_ldsc([self.rg_log_multiple], mode="rg")
        
        # Plot with custom parameters
        try:
            fig, ax, log, df = plot_rg(
                rg_data,
                sig_levels=[0.05, 0.01],
                panno=True,
                corrections=["non", "fdr", "bon"],
                verbose=False
            )
            success = True
        except Exception as e:
            success = False
            self.fail(f"plot_rg with custom parameters raised an exception: {e}")
        
        self.assertTrue(success, "plot_rg with custom parameters should complete")
        self.assertIsNotNone(fig, "plot_rg should return a figure")

    def test_plot_rg_with_single_pair(self):
        """Test plot_rg with single trait pair"""
        # Read rg data
        rg_data = read_ldsc([self.rg_log], mode="rg")
        
        # Plot should work even with single pair
        try:
            fig, ax, log, df = plot_rg(rg_data, verbose=False)
            success = True
        except Exception as e:
            success = False
            self.fail(f"plot_rg with single pair raised an exception: {e}")
        
        self.assertTrue(success, "plot_rg with single pair should complete")

    def test_plot_rg_dataframe_validation(self):
        """Test that plot_rg validates required columns"""
        # Create DataFrame with missing required columns
        invalid_data = pd.DataFrame({
            'p1': ['trait1'],
            'p2': ['trait2']
            # Missing 'rg' and 'p' columns
        })
        
        # Should raise an error or handle gracefully
        with self.assertRaises((KeyError, ValueError)):
            plot_rg(invalid_data, verbose=False)

    def test_plot_rg_with_na_values(self):
        """Test that plot_rg handles NA values correctly"""
        # Read rg data
        rg_data = read_ldsc([self.rg_log_multiple], mode="rg")
        
        # Add some NA values
        rg_data_with_na = rg_data.copy()
        rg_data_with_na.loc[0, 'p'] = np.nan
        
        # Should handle NA values gracefully
        try:
            fig, ax, log, df = plot_rg(rg_data_with_na, verbose=False)
            success = True
        except Exception as e:
            success = False
            self.fail(f"plot_rg with NA values raised an exception: {e}")
        
        self.assertTrue(success, "plot_rg should handle NA values")

    def test_read_ldsc_empty_filelist(self):
        """Test read_ldsc with empty file list"""
        result = read_ldsc([], mode="h2")
        
        # Should return empty DataFrame with correct columns
        self.assertIsInstance(result, pd.DataFrame)
        expected_columns = ['Filename', 'h2_obs', 'h2_se', 'Lambda_gc', 'Mean_chi2', 
                           'Intercept', 'Intercept_se', 'Ratio', 'Ratio_se']
        self.assertEqual(list(result.columns), expected_columns)
        self.assertEqual(len(result), 0)

    def test_plot_rg_annotation_options(self):
        """Test plot_rg with different annotation options"""
        # Read rg data
        rg_data = read_ldsc([self.rg_log_multiple], mode="rg")
        
        # Test with rganno option
        try:
            fig, ax, log, df = plot_rg(
                rg_data,
                rganno="rg",  # Show rg values
                panno=True,
                verbose=False
            )
            success = True
        except Exception as e:
            success = False
            self.fail(f"plot_rg with rganno raised an exception: {e}")
        
        self.assertTrue(success, "plot_rg with rganno should complete")


if __name__ == "__main__":
    unittest.main()

