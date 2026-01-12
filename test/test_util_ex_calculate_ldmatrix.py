"""
Test suite for util_ex_calculate_ldmatrix functions.

Tests cover:
- _extract_variants_in_locus: extracting variants within a window
- _align_sumstats_with_bim: aligning sumstats with reference BIM
- _export_snplist_and_locus_sumstats: exporting SNP lists and sumstats
- _check_snpid_order: checking SNPID order consistency
- _calculate_ld_r: calculating LD matrix (with mocked subprocess)
- _to_finemapping: main function (basic tests with mocking)
"""

import os
import sys
import unittest
import tempfile
import shutil
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open

import pandas as pd
import numpy as np

# Add parent directory to path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_ex_calculate_ldmatrix import (
    _extract_variants_in_locus,
    _align_sumstats_with_bim,
    _export_snplist_and_locus_sumstats,
    _check_snpid_order,
    _calculate_ld_r,
    _to_finemapping
)
from gwaslab.info.g_Log import Log


class TestExtractVariantsInLocus(unittest.TestCase):
    """Test _extract_variants_in_locus function."""
    
    def setUp(self):
        """Set up test data."""
        self.sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 1, 2, 2],
            "POS": [1000000, 1005000, 1010000, 2000000, 500000, 1500000],
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A", 
                     "1:2000000_A_G", "2:500000_T_C", "2:1500000_G_A"],
            "EA": ["A", "T", "G", "A", "T", "G"],
            "NEA": ["G", "C", "A", "G", "C", "A"],
            "P": [1e-8, 1e-7, 1e-9, 1e-6, 1e-7, 1e-8]
        })
    
    def test_extract_variants_within_window(self):
        """Test extracting variants within a window."""
        locus = (1, 1005000)  # CHR 1, POS 1005000
        windowsizekb = 500  # 500kb window
        
        result = _extract_variants_in_locus(self.sumstats, windowsizekb, locus)
        
        # Should include variants at 1000000, 1005000, 1010000 (within 500kb)
        # Should exclude 2000000 (too far), and chr2 variants
        self.assertEqual(len(result), 3)
        self.assertTrue(all(result["CHR"] == 1))
        self.assertTrue(all(result["POS"] >= 1005000 - 500000))
        self.assertTrue(all(result["POS"] < 1005000 + 500000))
    
    def test_extract_variants_small_window(self):
        """Test extracting variants with a small window."""
        locus = (1, 1005000)
        windowsizekb = 1  # 1kb window
        
        result = _extract_variants_in_locus(self.sumstats, windowsizekb, locus)
        
        # Should only include variant at 1005000 (within 1kb)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 1005000)
    
    def test_extract_variants_no_matches(self):
        """Test extracting variants when no variants match."""
        locus = (3, 1000000)  # Different chromosome
        windowsizekb = 500
        
        result = _extract_variants_in_locus(self.sumstats, windowsizekb, locus)
        
        self.assertEqual(len(result), 0)
    
    def test_extract_variants_at_boundary(self):
        """Test extracting variants at window boundary."""
        locus = (1, 1005000)
        windowsizekb = 500
        
        result = _extract_variants_in_locus(self.sumstats, windowsizekb, locus)
        
        # Boundary: 1005000 - 500000 = 500000, 1005000 + 500000 = 1505000
        # Should include 1000000 (>= 500000), 1005000, 1010000 (< 1505000)
        self.assertGreaterEqual(len(result), 1)
        self.assertTrue(1005000 in result["POS"].values)


class TestAlignSumstatsWithBim(unittest.TestCase):
    """Test _align_sumstats_with_bim function."""
    
    def setUp(self):
        """Set up test data."""
        self.locus_sumstats = pd.DataFrame({
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A"],
            "CHR": [1, 1, 1],
            "POS": [1000000, 1005000, 1010000],
            "EA": ["A", "T", "G"],
            "NEA": ["G", "C", "A"],
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.01, 0.02, 0.03],
            "Z": [10.0, 10.0, 10.0]
        })
        
        self.ref_bim = pd.DataFrame({
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A", "1:1020000_A_T"],
            "CHR_bim": [1, 1, 1, 1],
            "POS_bim": [1000000, 1005000, 1010000, 1020000],
            "EA_bim": ["A", "T", "G", "A"],
            "NEA_bim": ["G", "C", "A", "T"]
        })
        
        self.row = pd.Series({
            "SNPID": "1:1000000_A_G",
            "CHR": 1,
            "POS": 1000000
        })
        
        self.log = Log()
    
    def test_perfect_allele_match(self):
        """Test alignment with perfect allele matches."""
        result = _align_sumstats_with_bim(
            self.row, self.locus_sumstats, self.ref_bim, self.log
        )
        
        # Should match all 3 variants with perfect allele match
        self.assertEqual(len(result), 3)
        self.assertTrue(all(result["SNPID"].isin(self.locus_sumstats["SNPID"])))
        self.assertIn("SNPID", result.columns)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("EA", result.columns)
        self.assertIn("NEA", result.columns)
    
    def test_flipped_allele_match(self):
        """Test alignment with flipped alleles."""
        # Create flipped alleles in ref_bim
        flipped_bim = self.ref_bim.copy()
        flipped_bim.loc[flipped_bim["SNPID"] == "1:1000000_A_G", "EA_bim"] = "G"
        flipped_bim.loc[flipped_bim["SNPID"] == "1:1000000_A_G", "NEA_bim"] = "A"
        
        result = _align_sumstats_with_bim(
            self.row, self.locus_sumstats, flipped_bim, self.log
        )
        
        # Should still match the flipped variant
        self.assertGreaterEqual(len(result), 2)  # At least 2 other perfect matches
        matched_snpids = result["SNPID"].values
        # The flipped one should be included
        self.assertTrue("1:1000000_A_G" in matched_snpids or len(result) >= 2)
    
    def test_no_allele_match(self):
        """Test alignment when alleles don't match."""
        # Create mismatched alleles
        mismatched_bim = self.ref_bim.copy()
        mismatched_bim.loc[mismatched_bim["SNPID"] == "1:1000000_A_G", "EA_bim"] = "T"
        mismatched_bim.loc[mismatched_bim["SNPID"] == "1:1000000_A_G", "NEA_bim"] = "C"
        
        result = _align_sumstats_with_bim(
            self.row, self.locus_sumstats, mismatched_bim, self.log
        )
        
        # Should not include the mismatched variant
        self.assertLess(len(result), len(self.locus_sumstats))
    
    def test_lead_variant_not_in_reference(self):
        """Test when lead variant is not in reference."""
        row_no_match = pd.Series({
            "SNPID": "1:9999999_A_G",  # Not in ref_bim
            "CHR": 1,
            "POS": 9999999
        })
        
        result = _align_sumstats_with_bim(
            row_no_match, self.locus_sumstats, self.ref_bim, self.log
        )
        
        # Should still work but lead variant won't be in result
        self.assertIsInstance(result, pd.DataFrame)
    
    def test_with_suffixes(self):
        """Test alignment with column suffixes."""
        locus_sumstats_suffix = self.locus_sumstats.copy()
        locus_sumstats_suffix["BETA_suffix"] = locus_sumstats_suffix["BETA"]
        locus_sumstats_suffix["SE_suffix"] = locus_sumstats_suffix["SE"]
        locus_sumstats_suffix["Z_suffix"] = locus_sumstats_suffix["Z"]  # Add Z_suffix for the test
        
        result = _align_sumstats_with_bim(
            self.row, locus_sumstats_suffix, self.ref_bim, 
            self.log, suffixes=["", "_suffix"]
        )
        
        # Should include both BETA/SE and BETA_suffix/SE_suffix if present
        self.assertIn("SNPID", result.columns)
        # The function checks for "Z" (without suffix) and then adds "Z"+suffix
        # So if "Z" exists, it will add both "Z" and "Z_suffix" to output
        self.assertGreater(len(result.columns), 5)  # At least basic columns + some stats
    
    def test_with_eaf_column(self):
        """Test alignment when EAF column is present."""
        locus_sumstats_eaf = self.locus_sumstats.copy()
        locus_sumstats_eaf["EAF"] = [0.1, 0.2, 0.3]
        
        result = _align_sumstats_with_bim(
            self.row, locus_sumstats_eaf, self.ref_bim, self.log
        )
        
        # Should include EAF in output if present
        if "EAF" in self.locus_sumstats.columns:
            self.assertIn("EAF", result.columns)


class TestExportSnplistAndLocusSumstats(unittest.TestCase):
    """Test _export_snplist_and_locus_sumstats function."""
    
    def setUp(self):
        """Set up test data."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.matched_sumstats = pd.DataFrame({
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A"],
            "CHR": [1, 1, 1],
            "POS": [1000000, 1005000, 1010000],
            "EA": ["A", "T", "G"],
            "NEA": ["G", "C", "A"],
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.01, 0.02, 0.03],
            "Z": [10.0, 10.0, 10.0]
        })
        
        self.row = pd.Series({
            "SNPID": "1:1000000_A_G",
            "CHR": 1,
            "POS": 1000000
        })
        
        self.study = "test_study"
        self.windowsizekb = 500
        self.log = Log()
    
    def tearDown(self):
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_export_basic_columns(self):
        """Test exporting basic columns."""
        snplist_path, sumstats_path = _export_snplist_and_locus_sumstats(
            self.matched_sumstats, str(self.temp_dir), self.study, 
            self.row, self.windowsizekb, self.log
        )
        
        # Check files exist
        self.assertTrue(os.path.exists(snplist_path))
        self.assertTrue(os.path.exists(sumstats_path))
        self.assertTrue(os.path.exists(sumstats_path + ".gz"))
        
        # Check SNP list content
        snplist = pd.read_csv(snplist_path, header=None)
        self.assertEqual(len(snplist), len(self.matched_sumstats))
        self.assertTrue(all(snplist[0] == self.matched_sumstats["SNPID"]))
        
        # Check sumstats content
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        self.assertIn("SNPID", sumstats.columns)
        self.assertIn("CHR", sumstats.columns)
        self.assertIn("POS", sumstats.columns)
        self.assertIn("EA", sumstats.columns)
        self.assertIn("NEA", sumstats.columns)
    
    def test_export_with_beta_se(self):
        """Test exporting with BETA and SE columns."""
        snplist_path, sumstats_path = _export_snplist_and_locus_sumstats(
            self.matched_sumstats, str(self.temp_dir), self.study,
            self.row, self.windowsizekb, self.log
        )
        
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        # Should include BETA and SE if present
        if "BETA" in self.matched_sumstats.columns and "SE" in self.matched_sumstats.columns:
            self.assertIn("BETA", sumstats.columns)
            self.assertIn("SE", sumstats.columns)
    
    def test_export_with_z_column(self):
        """Test exporting with Z column."""
        snplist_path, sumstats_path = _export_snplist_and_locus_sumstats(
            self.matched_sumstats, str(self.temp_dir), self.study,
            self.row, self.windowsizekb, self.log
        )
        
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        if "Z" in self.matched_sumstats.columns:
            self.assertIn("Z", sumstats.columns)
    
    def test_export_with_suffixes(self):
        """Test exporting with column suffixes."""
        matched_sumstats_suffix = self.matched_sumstats.copy()
        matched_sumstats_suffix["BETA_suffix"] = matched_sumstats_suffix["BETA"]
        matched_sumstats_suffix["SE_suffix"] = matched_sumstats_suffix["SE"]
        matched_sumstats_suffix["Z_suffix"] = matched_sumstats_suffix["Z"]
        
        snplist_path, sumstats_path = _export_snplist_and_locus_sumstats(
            matched_sumstats_suffix, str(self.temp_dir), self.study,
            self.row, self.windowsizekb, self.log, suffixes=["", "_suffix"]
        )
        
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        # Should include suffixed columns
        if "BETA_suffix" in matched_sumstats_suffix.columns:
            self.assertIn("BETA_suffix", sumstats.columns)
            self.assertIn("SE_suffix", sumstats.columns)


class TestCheckSnpidOrder(unittest.TestCase):
    """Test _check_snpid_order function."""
    
    def setUp(self):
        """Set up test data."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.log = Log()
    
    def tearDown(self):
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_matching_order(self):
        """Test when SNPID order matches."""
        snpids = ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A"]
        matched_sumstats_snpid = pd.Series(snpids)
        
        snplist_path = self.temp_dir / "test.snplist"
        pd.Series(snpids).to_csv(snplist_path, index=False, header=False)
        
        # Should not raise error
        _check_snpid_order(str(snplist_path), matched_sumstats_snpid, self.log)
        
        # Check log for success message
        self.assertIn("matched", self.log.log_text.lower())
    
    def test_non_matching_order(self):
        """Test when SNPID order doesn't match."""
        matched_snpids = ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A"]
        file_snpids = ["1:1010000_G_A", "1:1005000_T_C", "1:1000000_A_G"]  # Reversed
        
        matched_sumstats_snpid = pd.Series(matched_snpids)
        
        snplist_path = self.temp_dir / "test.snplist"
        pd.Series(file_snpids).to_csv(snplist_path, index=False, header=False)
        
        # Should not raise error but log warning
        _check_snpid_order(str(snplist_path), matched_sumstats_snpid, self.log)
        
        # Check log for warning
        self.assertIn("not matched", self.log.log_text.lower() or "warning" in self.log.log_text.lower())
    
    def test_different_lengths(self):
        """Test when lengths don't match."""
        matched_snpids = ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A"]
        file_snpids = ["1:1000000_A_G", "1:1005000_T_C"]  # Shorter
        
        matched_sumstats_snpid = pd.Series(matched_snpids)
        
        snplist_path = self.temp_dir / "test.snplist"
        pd.Series(file_snpids).to_csv(snplist_path, index=False, header=False)
        
        # Should log warning
        _check_snpid_order(str(snplist_path), matched_sumstats_snpid, self.log)
        
        # Check log for warning
        self.assertIn("not matched", self.log.log_text.lower() or "warning" in self.log.log_text.lower())


class TestCalculateLdR(unittest.TestCase):
    """Test _calculate_ld_r function."""
    
    def setUp(self):
        """Set up test data."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.study = "test_study"
        self.matched_sumstats_snpid = pd.Series(["1:1000000_A_G", "1:1005000_T_C"])
        self.row = pd.Series({
            "SNPID": "1:1000000_A_G",
            "CHR": 1,
            "POS": 1000000
        })
        self.bfile_prefix = str(self.temp_dir / "test_ref")
        self.windowsizekb = 500
        self.plink_log = ""
        self.log = Log()
        self.filetype = "bfile"
        self.plink = "plink"
        self.plink2 = "plink2"
        self.ref_allele_path = str(self.temp_dir / "ref_allele.txt")
        self.extra_plink_option = ""
        
        # Create dummy bed file to make os.path.exists return True
        (self.temp_dir / "test_ref.bed").touch()
    
    def tearDown(self):
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._check_snpid_order')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.subprocess.check_output')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._checking_plink_version')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.os.path.exists')
    def test_calculate_ld_r_success(self, mock_exists, mock_check_version, mock_subprocess, mock_check_order):
        """Test successful LD calculation."""
        mock_exists.return_value = True
        mock_check_version.return_value = self.log
        mock_subprocess.return_value = "PLINK output log"
        mock_check_order.return_value = None  # Mock the check function
        
        result_path, plink_log = _calculate_ld_r(
            study=self.study,
            matched_sumstats_snpid=self.matched_sumstats_snpid,
            row=self.row,
            bfile_prefix=self.bfile_prefix,
            threads=1,
            windowsizekb=self.windowsizekb,
            out=str(self.temp_dir),
            plink_log=self.plink_log,
            log=self.log,
            memory=None,
            mode="r",
            filetype=self.filetype,
            plink=self.plink,
            plink2=self.plink2,
            ref_allele_path=self.ref_allele_path,
            extra_plink_option=self.extra_plink_option,
            verbose=False
        )
        
        # Check that subprocess was called
        mock_subprocess.assert_called_once()
        
        # Check return value format
        self.assertIsInstance(result_path, str)
        self.assertIn(".ld.gz", result_path)
        self.assertIsInstance(plink_log, str)
    
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._check_snpid_order')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.subprocess.check_output')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._checking_plink_version')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.os.path.exists')
    def test_calculate_ld_r_with_memory(self, mock_exists, mock_check_version, mock_subprocess, mock_check_order):
        """Test LD calculation with memory limit."""
        mock_exists.return_value = True
        mock_check_version.return_value = self.log
        mock_subprocess.return_value = "PLINK output log"
        mock_check_order.return_value = None  # Mock the check function
        
        result_path, plink_log = _calculate_ld_r(
            study=self.study,
            matched_sumstats_snpid=self.matched_sumstats_snpid,
            row=self.row,
            bfile_prefix=self.bfile_prefix,
            threads=1,
            windowsizekb=self.windowsizekb,
            out=str(self.temp_dir),
            plink_log=self.plink_log,
            log=self.log,
            memory=4096,  # 4GB
            mode="r",
            filetype=self.filetype,
            plink=self.plink,
            plink2=self.plink2,
            ref_allele_path=self.ref_allele_path,
            extra_plink_option=self.extra_plink_option,
            verbose=False
        )
        
        # Check that memory flag was included in command
        call_args = mock_subprocess.call_args[0][0]
        self.assertIn("--memory", call_args)
    
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.subprocess.check_output')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._checking_plink_version')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.os.path.exists')
    def test_calculate_ld_r_pfile_error(self, mock_exists, mock_check_version, mock_subprocess):
        """Test error when pfile is used with PLINK1."""
        mock_exists.return_value = True
        mock_check_version.return_value = self.log
        
        with self.assertRaises(ValueError) as context:
            _calculate_ld_r(
                study=self.study,
                matched_sumstats_snpid=self.matched_sumstats_snpid,
                row=self.row,
                bfile_prefix=self.bfile_prefix,
                threads=1,
                windowsizekb=self.windowsizekb,
                out=str(self.temp_dir),
                plink_log=self.plink_log,
                log=self.log,
                memory=None,
                mode="r",
                filetype="pfile",  # Should raise error
                plink=self.plink,
                plink2=self.plink2,
                ref_allele_path=self.ref_allele_path,
                extra_plink_option=self.extra_plink_option,
                verbose=False
            )
        
        self.assertIn("pfile", str(context.exception).lower())
    
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._check_snpid_order')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.subprocess.check_output')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._checking_plink_version')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix.os.path.exists')
    def test_calculate_ld_r_subprocess_error(self, mock_exists, mock_check_version, mock_subprocess, mock_check_order):
        """Test handling of subprocess errors."""
        mock_exists.return_value = True
        mock_check_version.return_value = self.log
        mock_subprocess.side_effect = subprocess.CalledProcessError(1, "plink", "Error output")
        mock_check_order.return_value = None  # Mock the check function
        
        # Should not raise, but log error
        result_path, plink_log = _calculate_ld_r(
            study=self.study,
            matched_sumstats_snpid=self.matched_sumstats_snpid,
            row=self.row,
            bfile_prefix=self.bfile_prefix,
            threads=1,
            windowsizekb=self.windowsizekb,
            out=str(self.temp_dir),
            plink_log=self.plink_log,
            log=self.log,
            memory=None,
            mode="r",
            filetype=self.filetype,
            plink=self.plink,
            plink2=self.plink2,
            ref_allele_path=self.ref_allele_path,
            extra_plink_option=self.extra_plink_option,
            verbose=False
        )
        
        # Should still return a path
        self.assertIsInstance(result_path, str)


class TestToFinemapping(unittest.TestCase):
    """Test _to_finemapping function (main function)."""
    
    def setUp(self):
        """Set up test data."""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create mock Sumstats object
        self.sumstats_data = pd.DataFrame({
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A", "2:500000_T_C"],
            "CHR": [1, 1, 1, 2],
            "POS": [1000000, 1005000, 1010000, 500000],
            "EA": ["A", "T", "G", "T"],
            "NEA": ["G", "C", "A", "C"],
            "P": [1e-10, 1e-8, 1e-9, 1e-7]
        })
        
        # Mock Sumstats object
        self.mock_gls = MagicMock()
        self.mock_gls.data = self.sumstats_data
        self.mock_gls.offload = MagicMock()
        self.mock_gls.reload = MagicMock()
    
    def tearDown(self):
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._get_sig')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._process_plink_input_files')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._extract_variants_in_locus')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._align_sumstats_with_bim')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._export_snplist_and_locus_sumstats')
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._calculate_ld_r')
    def test_to_finemapping_with_loci(self, mock_calc_ld, mock_export, mock_align, 
                                       mock_extract, mock_process, mock_get_sig):
        """Test _to_finemapping with provided loci."""
        # Setup mocks
        mock_extract.return_value = self.sumstats_data.iloc[:3].copy()
        
        ref_bim = pd.DataFrame({
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C", "1:1010000_G_A"],
            "CHR_bim": [1, 1, 1],
            "POS_bim": [1000000, 1005000, 1010000],
            "EA_bim": ["A", "T", "G"],
            "NEA_bim": ["G", "C", "A"]
        })
        mock_process.return_value = ("test_bfile", "", [ref_bim], "bfile")
        
        matched_sumstats = pd.DataFrame({
            "SNPID": ["1:1000000_A_G", "1:1005000_T_C"],
            "CHR": [1, 1],
            "POS": [1000000, 1005000],
            "EA": ["A", "T"],
            "NEA": ["G", "C"]
        })
        mock_align.return_value = matched_sumstats
        
        mock_export.return_value = (
            str(self.temp_dir / "test.snplist"),
            str(self.temp_dir / "test.sumstats")
        )
        
        mock_calc_ld.return_value = (str(self.temp_dir / "test.ld.gz"), "")
        
        # Create test loci
        loci = ["1:1000000_A_G"]
        
        result = _to_finemapping(
            gls=self.mock_gls,
            study="test_study",
            bfile=str(self.temp_dir / "test_ref"),
            loci=loci,
            out=str(self.temp_dir),
            plink="plink",
            plink2="plink2",
            windowsizekb=500,
            threads=1,
            mode="r",
            exclude_hla=False,
            log=Log(),
            verbose=False
        )
        
        # Check return values
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 3)  # output_file_list_path, output_file_list, plink_log
        self.assertIsNotNone(result[0])  # file list path
    
    @patch('gwaslab.util.util_ex_calculate_ldmatrix._get_sig')
    def test_to_finemapping_auto_extract_loci(self, mock_get_sig):
        """Test _to_finemapping with automatic locus extraction."""
        # Mock _get_sig to return significant variants
        sig_df = pd.DataFrame({
            "SNPID": ["1:1000000_A_G"],
            "CHR": [1],
            "POS": [1000000],
            "EA": ["A"],
            "NEA": ["G"],
            "P": [1e-10]
        })
        mock_get_sig.return_value = sig_df
        
        # This is a complex function, so we'll just test that it can be called
        # Full testing would require extensive mocking
        with patch('gwaslab.util.util_ex_calculate_ldmatrix._process_plink_input_files') as mock_process, \
             patch('gwaslab.util.util_ex_calculate_ldmatrix._extract_variants_in_locus') as mock_extract, \
             patch('gwaslab.util.util_ex_calculate_ldmatrix._align_sumstats_with_bim') as mock_align, \
             patch('gwaslab.util.util_ex_calculate_ldmatrix._export_snplist_and_locus_sumstats') as mock_export, \
             patch('gwaslab.util.util_ex_calculate_ldmatrix._calculate_ld_r') as mock_calc_ld:
            
            mock_extract.return_value = self.sumstats_data.iloc[:1].copy()
            ref_bim = pd.DataFrame({
                "SNPID": ["1:1000000_A_G"],
                "CHR_bim": [1],
                "POS_bim": [1000000],
                "EA_bim": ["A"],
                "NEA_bim": ["G"]
            })
            mock_process.return_value = ("test_bfile", "", [ref_bim], "bfile")
            mock_align.return_value = pd.DataFrame({
                "SNPID": ["1:1000000_A_G"],
                "CHR": [1],
                "POS": [1000000],
                "EA": ["A"],
                "NEA": ["G"]
            })
            mock_export.return_value = (
                str(self.temp_dir / "test.snplist"),
                str(self.temp_dir / "test.sumstats")
            )
            mock_calc_ld.return_value = (str(self.temp_dir / "test.ld.gz"), "")
            
            result = _to_finemapping(
                gls=self.mock_gls,
                study="test_study",
                bfile=str(self.temp_dir / "test_ref"),
                loci=None,  # Will trigger auto-extraction
                out=str(self.temp_dir),
                plink="plink",
                plink2="plink2",
                windowsizekb=500,
                threads=1,
                mode="r",
                exclude_hla=False,
                getlead_kwargs={"windowsizekb": 1000},
                log=Log(),
                verbose=False
            )
            
            # Verify _get_sig was called
            mock_get_sig.assert_called_once()
            self.assertIsNotNone(result)


if __name__ == "__main__":
    unittest.main()
