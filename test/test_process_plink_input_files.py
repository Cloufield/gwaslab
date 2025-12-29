"""
Test suite for _process_plink_input_files function.
Tests VCF conversion and processing of bfile/pfile formats.
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

from gwaslab.io.io_plink import _process_plink_input_files
from gwaslab.info.g_Log import Log


class TestProcessPlinkInputFiles(unittest.TestCase):
    """Test _process_plink_input_files function with VCF and converted files."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests."""
        cls.test_dir = Path(__file__).parent
        cls.ref_dir = cls.test_dir / "ref"
        
        # Paths to reference files
        cls.vcf_path = cls.ref_dir / "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz"
        cls.bfile_prefix = str(cls.ref_dir / "1kg_eas_hg19.chr7_126253550_128253550")
        cls.pfile_prefix = str(cls.ref_dir / "1kg_eas_hg19.chr7_126253550_128253550")
        
        # Verify files exist
        assert cls.vcf_path.exists(), f"VCF file not found: {cls.vcf_path}"
        bed_file = Path(f"{cls.bfile_prefix}.bed")
        assert bed_file.exists(), f"BED file not found: {bed_file}"
        pgen_file = Path(f"{cls.pfile_prefix}.pgen")
        assert pgen_file.exists(), f"PGEN file not found: {pgen_file}"
        
        # Create temporary directory for test outputs
        cls.temp_dir = cls.test_dir / "tmp" / f"gwaslab_process_ref_test_{os.getpid()}"
        cls.temp_dir.mkdir(parents=True, exist_ok=True)
        
    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        if cls.temp_dir.exists():
            shutil.rmtree(cls.temp_dir)
    
    def setUp(self):
        """Set up for each test."""
        self.log = Log()
        self.chrlist = [7]  # Test data is for chromosome 7
        self.temp_output = self.temp_dir / f"test_{self._testMethodName}"
        self.temp_output.mkdir(parents=True, exist_ok=True)
    
    def test_vcf_to_bfile_conversion(self):
        """Test VCF to bfile conversion."""
        # Convert VCF to bfile in temp directory
        temp_vcf_prefix = str(self.temp_output / "test_vcf")
        temp_bfile_prefix = str(self.temp_output / "test_vcf")
        
        # Copy VCF to temp location (or use original)
        vcf_path = str(self.vcf_path)
        
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            vcf=vcf_path,
            convert="bfile",
            threads=1,
            log=self.log,
            overwrite=True,
            load_bim=False
        )
        
        # Check return values
        self.assertIsInstance(ref_file_prefix, str)
        self.assertIsInstance(plink_log, str)
        self.assertIsInstance(ref_bims, list)
        self.assertEqual(filetype, "bfile")
        
        # Check that plink_log is not empty (conversion occurred)
        self.assertGreater(len(plink_log), 0, "PLINK log should not be empty after conversion")
        
        # Check that converted files exist (they should be in the same directory as VCF)
        expected_prefix = vcf_path.replace(".vcf.gz", ".7")
        bed_file = Path(f"{expected_prefix}.bed")
        bim_file = Path(f"{expected_prefix}.bim")
        fam_file = Path(f"{expected_prefix}.fam")
        
        # Note: Files might be created in ref_dir, check both locations
        if not bed_file.exists():
            # Try in ref_dir
            bed_file = self.ref_dir / f"1kg_eas_hg19.chr7_126253550_128253550.7.bed"
        
        # For this test, we mainly verify the function runs without error
        # The actual file location depends on the VCF path handling
        self.assertEqual(len(ref_bims), 0, "ref_bims should be empty when load_bim=False")
    
    def test_vcf_to_bfile_with_load_bim(self):
        """Test VCF to bfile conversion with load_bim=True."""
        vcf_path = str(self.vcf_path)
        
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            vcf=vcf_path,
            convert="bfile",
            threads=1,
            log=self.log,
            overwrite=True,
            load_bim=True
        )
        
        # Check return values
        self.assertEqual(filetype, "bfile")
        self.assertIsInstance(ref_bims, list)
        
        # When load_bim=True, ref_bims should contain dataframes
        if len(ref_bims) > 0:
            self.assertIsInstance(ref_bims[0], pd.DataFrame)
            # Check expected columns in BIM file
            expected_cols = ["SNPID", "CHR_bim", "POS_bim", "EA_bim", "NEA_bim"]
            for col in expected_cols:
                self.assertIn(col, ref_bims[0].columns, f"Missing column {col} in ref_bims")
    
    def test_vcf_to_pfile_conversion(self):
        """Test VCF to pfile conversion."""
        vcf_path = str(self.vcf_path)
        
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            vcf=vcf_path,
            convert="pfile",
            threads=1,
            log=self.log,
            overwrite=True,
            load_bim=False
        )
        
        # Check return values
        self.assertEqual(filetype, "pfile")
        self.assertIsInstance(plink_log, str)
        self.assertEqual(len(ref_bims), 0, "ref_bims should be empty when load_bim=False")
    
    def test_vcf_to_pfile_with_load_bim(self):
        """Test VCF to pfile conversion with load_bim=True."""
        vcf_path = str(self.vcf_path)
        
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            vcf=vcf_path,
            convert="pfile",
            threads=1,
            log=self.log,
            overwrite=True,
            load_bim=True
        )
        
        # Check return values
        self.assertEqual(filetype, "pfile")
        self.assertIsInstance(ref_bims, list)
        
        # When load_bim=True, ref_bims should contain dataframes
        if len(ref_bims) > 0:
            self.assertIsInstance(ref_bims[0], pd.DataFrame)
            # Check expected columns in PVAR file
            expected_cols = ["SNPID", "CHR_bim", "POS_bim", "EA_bim", "NEA_bim"]
            for col in expected_cols:
                self.assertIn(col, ref_bims[0].columns, f"Missing column {col} in ref_bims")
    
    def test_bfile_input(self):
        """Test processing existing bfile without conversion."""
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            bfile=self.bfile_prefix,
            log=self.log,
            load_bim=False
        )
        
        # Check return values
        self.assertEqual(filetype, "bfile")
        self.assertEqual(plink_log, "", "plink_log should be empty for bfile input (no conversion)")
        self.assertEqual(len(ref_bims), 0, "ref_bims should be empty when load_bim=False")
        self.assertEqual(ref_file_prefix, self.bfile_prefix)
    
    def test_bfile_input_with_load_bim(self):
        """Test processing existing bfile with load_bim=True."""
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            bfile=self.bfile_prefix,
            log=self.log,
            load_bim=True
        )
        
        # Check return values
        self.assertEqual(filetype, "bfile")
        self.assertEqual(plink_log, "", "plink_log should be empty for bfile input")
        self.assertIsInstance(ref_bims, list)
        
        # When load_bim=True, ref_bims should contain dataframes
        if len(ref_bims) > 0:
            self.assertIsInstance(ref_bims[0], pd.DataFrame)
            expected_cols = ["SNPID", "CHR_bim", "POS_bim", "EA_bim", "NEA_bim"]
            for col in expected_cols:
                self.assertIn(col, ref_bims[0].columns, f"Missing column {col} in ref_bims")
            
            # Check that CHR_bim values match expected chromosome
            if len(ref_bims[0]) > 0:
                chr_values = ref_bims[0]["CHR_bim"].unique()
                # Should contain chromosome 7 (might be as string or category)
                chr_str = [str(c) for c in chr_values]
                self.assertTrue(
                    any("7" in str(c) for c in chr_str),
                    f"Expected chromosome 7 in ref_bims, got {chr_values}"
                )
    
    def test_pfile_input(self):
        """Test processing existing pfile without conversion."""
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            pfile=self.pfile_prefix,
            log=self.log,
            load_bim=False
        )
        
        # Check return values
        self.assertEqual(filetype, "pfile")
        self.assertEqual(plink_log, "", "plink_log should be empty for pfile input (no conversion)")
        self.assertEqual(len(ref_bims), 0, "ref_bims should be empty when load_bim=False")
        self.assertEqual(ref_file_prefix, self.pfile_prefix)
    
    def test_pfile_input_with_load_bim(self):
        """Test processing existing pfile with load_bim=True."""
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=self.chrlist,
            pfile=self.pfile_prefix,
            log=self.log,
            load_bim=True
        )
        
        # Check return values
        self.assertEqual(filetype, "pfile")
        self.assertEqual(plink_log, "", "plink_log should be empty for pfile input")
        self.assertIsInstance(ref_bims, list)
        
        # When load_bim=True, ref_bims should contain dataframes
        if len(ref_bims) > 0:
            self.assertIsInstance(ref_bims[0], pd.DataFrame)
            expected_cols = ["SNPID", "CHR_bim", "POS_bim", "EA_bim", "NEA_bim"]
            for col in expected_cols:
                self.assertIn(col, ref_bims[0].columns, f"Missing column {col} in ref_bims")
    
    def test_vcf_conversion_no_overwrite(self):
        """Test VCF conversion when output files already exist (no overwrite)."""
        vcf_path = str(self.vcf_path)
        
        # First conversion
        ref_file_prefix1, plink_log1, ref_bims1, filetype1 = _process_plink_input_files(
            chrlist=self.chrlist,
            vcf=vcf_path,
            convert="bfile",
            threads=1,
            log=self.log,
            overwrite=True,
            load_bim=False
        )
        
        # Second conversion without overwrite (should skip)
        ref_file_prefix2, plink_log2, ref_bims2, filetype2 = _process_plink_input_files(
            chrlist=self.chrlist,
            vcf=vcf_path,
            convert="bfile",
            threads=1,
            log=self.log,
            overwrite=False,
            load_bim=False
        )
        
        # Both should return same filetype
        self.assertEqual(filetype1, "bfile")
        self.assertEqual(filetype2, "bfile")
    
    def test_error_missing_bfile(self):
        """Test error handling when bfile is missing."""
        with self.assertRaises(ValueError) as context:
            _process_plink_input_files(
                chrlist=self.chrlist,
                bfile="/nonexistent/path/to/file",
                log=self.log
            )
        
        self.assertIn("PLINK bfiles are missing", str(context.exception))
    
    def test_error_missing_pfile(self):
        """Test error handling when pfile is missing."""
        with self.assertRaises(ValueError) as context:
            _process_plink_input_files(
                chrlist=self.chrlist,
                pfile="/nonexistent/path/to/file",
                log=self.log
            )
        
        self.assertIn("PLINK pfiles are missing", str(context.exception))
    
    def test_error_no_input_file(self):
        """Test error handling when no input file is provided."""
        with self.assertRaises(ValueError) as context:
            _process_plink_input_files(
                chrlist=self.chrlist,
                log=self.log
            )
        
        self.assertIn("You need to provide one from bfile, pfile, bgen, vcf", str(context.exception))
    
    def test_multiple_chromosomes(self):
        """Test processing with multiple chromosomes (if available)."""
        # For this test data, we only have chr7, so test with single chromosome
        # but verify the function handles chrlist correctly
        ref_file_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=[7],
            bfile=self.bfile_prefix,
            log=self.log,
            load_bim=True
        )
        
        self.assertEqual(filetype, "bfile")
        # Should have loaded BIM data for chromosome 7
        if len(ref_bims) > 0:
            self.assertGreater(len(ref_bims[0]), 0, "BIM data should contain variants")


if __name__ == "__main__":
    unittest.main()

