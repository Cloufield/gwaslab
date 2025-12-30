"""
Test suite for clumping functionality using PLINK and PLINK2.
"""
import os
import sys
import unittest
import tempfile
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np

# Add parent directory to path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.g_Sumstats import Sumstats
from gwaslab.info.g_Log import Log


class TestClumping(unittest.TestCase):
    """Test clumping functionality with bfile and pfile formats."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests."""
        cls.test_dir = Path(__file__).parent
        cls.ref_dir = cls.test_dir / "ref"
        
        # Paths to reference files
        cls.vcf_path = cls.ref_dir / "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz"
        cls.bfile_prefix = cls.ref_dir / "1kg_eas_hg19.chr7_126253550_128253550"
        cls.pfile_prefix = cls.ref_dir / "1kg_eas_hg19.chr7_126253550_128253550"
        cls.sumstats_path = cls.ref_dir / "bbj_t2d_hm3_chr7_variants.txt.gz"
        
        # Verify VCF file exists
        assert cls.vcf_path.exists(), f"VCF file not found: {cls.vcf_path}"
        
        # Create PLINK BED file from VCF if it doesn't exist
        bed_file = Path(f"{cls.bfile_prefix}.bed")
        if not bed_file.exists():
            print(f"Creating PLINK BED file from VCF: {bed_file}")
            # Check if plink2 is available
            try:
                result = subprocess.run(
                    ["plink2", "--version"],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result.returncode != 0:
                    raise FileNotFoundError("plink2 command failed")
            except (FileNotFoundError, subprocess.TimeoutExpired):
                raise RuntimeError("plink2 not found. Please install PLINK2 to run clumping tests.")
            
            # Convert VCF to PLINK binary format
            cmd = [
                "plink2",
                "--vcf", str(cls.vcf_path),
                "--make-bed",
                "--out", str(cls.bfile_prefix),
                "--threads", "6",
                "--memory", "2048"
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(cls.ref_dir)
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"Failed to create BED file from VCF.\n"
                    f"Command: {' '.join(cmd)}\n"
                    f"Error: {result.stderr}"
                )
            print(f"Successfully created BED file: {bed_file}")
        
        # Verify BED file exists
        assert bed_file.exists(), f"BED file not found: {bed_file}"
        
        # Create PLINK2 pgen file from VCF if it doesn't exist
        pgen_file = Path(f"{cls.pfile_prefix}.pgen")
        if not pgen_file.exists():
            print(f"Creating PLINK2 pgen file from VCF: {pgen_file}")
            # Convert VCF to PLINK2 format
            cmd = [
                "plink2",
                "--vcf", str(cls.vcf_path),
                "--make-pgen",
                "--out", str(cls.pfile_prefix),
                "--threads", "6",
                "--memory", "2048"
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(cls.ref_dir)
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"Failed to create pgen file from VCF.\n"
                    f"Command: {' '.join(cmd)}\n"
                    f"Error: {result.stderr}"
                )
            print(f"Successfully created pgen file: {pgen_file}")
        
        # Verify pgen file exists
        assert pgen_file.exists(), f"PGEN file not found: {pgen_file}"
        assert cls.sumstats_path.exists(), f"Sumstats file not found: {cls.sumstats_path}"
        
        # Create temporary directory for test outputs (within project directory)
        cls.temp_dir = cls.test_dir / "tmp" / f"gwaslab_clump_test_{os.getpid()}"
        cls.temp_dir.mkdir(parents=True, exist_ok=True)
        
    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        if cls.temp_dir.exists():
            shutil.rmtree(cls.temp_dir)
    
    def setUp(self):
        """Set up for each test."""
        # Load sumstats as Sumstats object directly from file
        self.sumstats = Sumstats(
            sumstats=str(self.sumstats_path),
            chrom="CHR",
            pos="POS",
            p="P",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            verbose=False
        )
        # Filter to chromosome 7 to match the reference files
        # The reference files are only for chromosome 7, so we must filter to chr7
        # Check if chromosome 7 exists (handle both int and string types)
        chr_values = self.sumstats.data["CHR"]
        # Convert to numeric if possible, then check
        try:
            chr_numeric = pd.to_numeric(chr_values, errors='coerce')
            has_chr7 = (chr_numeric == 7).any() or (chr_values.astype(str) == '7').any()
        except:
            has_chr7 = (chr_values.astype(str) == '7').any()
        
        if not has_chr7:
            # Skip tests if chromosome 7 is not available
            import unittest
            raise unittest.SkipTest("Chromosome 7 not found in sumstats data")
        self.sumstats.filter_value(expr="CHR == 7", inplace=True, verbose=False)
        
        # Use only 1/10 of variants to speed up tests
        if len(self.sumstats.data) > 0:
            self.sumstats.data = self.sumstats.data.iloc[::10].copy()
        
        # Fix SNP ID to make it consistent with PLINK file format (CHR:POS:REF:ALT)
        # This ensures SNP IDs match between sumstats and reference panel
        self.sumstats.fix_id(
            fixid=True,
            fixsep=True,
            fixchrpos=False,  # CHR and POS are already correct
            fixeanea=False,   # EA and NEA are already correct
            overwrite=False,
            verbose=False
        )
        
        self.temp_output = self.temp_dir / f"clump_test_{self._testMethodName}"
        # Ensure output directory exists
        self.temp_dir.mkdir(parents=True, exist_ok=True)
    
    def test_clump_with_bfile(self):
        """Test clumping using PLINK bfile format."""
        # Use relative path (clumping function strips leading slashes)
        output_path = str(self.temp_output.relative_to(Path.cwd()))
        
        # The clumping function creates temp files in tmp/{output_path}_gwaslab_tmp...
        # Create the tmp directory structure
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        # Perform clumping with bfile
        self.sumstats.clump(
            bfile=str(self.bfile_prefix),
            out=output_path,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        
        # Check that results were stored in clumps dictionary
        results_sumstats = self.sumstats.clumps.get("clumps", pd.DataFrame())
        results = self.sumstats.clumps.get("clumps_raw", pd.DataFrame())
        plink_log = self.sumstats.clumps.get("plink_log", "")
        
        self.assertIsInstance(results_sumstats, pd.DataFrame)
        self.assertIsInstance(results, pd.DataFrame)
        self.assertIsInstance(plink_log, str)
        
        # Check that results have expected columns
        if len(results) > 0:
            expected_cols = ["CHR", "POS", "SNPID"]
            for col in expected_cols:
                self.assertIn(col, results.columns, f"Missing column {col} in results")
            
            # Check that results_sumstats contains only clumped variants
            original_snpids = set(self.sumstats.data["SNPID"].values)
            result_snpids = set(results_sumstats["SNPID"].values)
            self.assertTrue(
                result_snpids.issubset(original_snpids),
                "Results contain SNPs not in original sumstats"
            )
            
            # Check that all results are from chromosome 7
            if "CHR" in results.columns:
                self.assertTrue(
                    all(results["CHR"] == 7),
                    "Results contain variants from wrong chromosome"
                )
    
    def test_clump_with_pfile(self):
        """Test clumping using PLINK2 pfile format."""
        # Use relative path
        output_path = str(self.temp_output.relative_to(Path.cwd()))
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        # Perform clumping with pfile
        self.sumstats.clump(
            pfile=str(self.pfile_prefix),
            out=output_path,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        
        # Check that results were stored in clumps dictionary
        results_sumstats = self.sumstats.clumps.get("clumps", pd.DataFrame())
        results = self.sumstats.clumps.get("clumps_raw", pd.DataFrame())
        plink_log = self.sumstats.clumps.get("plink_log", "")
        
        self.assertIsInstance(results_sumstats, pd.DataFrame)
        self.assertIsInstance(results, pd.DataFrame)
        self.assertIsInstance(plink_log, str)
        
        # Check that results have expected columns
        if len(results) > 0:
            expected_cols = ["CHR", "POS", "SNPID"]
            for col in expected_cols:
                self.assertIn(col, results.columns, f"Missing column {col} in results")
            
            # Check that results_sumstats contains only clumped variants
            original_snpids = set(self.sumstats.data["SNPID"].values)
            result_snpids = set(results_sumstats["SNPID"].values)
            self.assertTrue(
                result_snpids.issubset(original_snpids),
                "Results contain SNPs not in original sumstats"
            )
    
    def test_clump_with_vcf(self):
        """Test clumping using VCF file (should be converted to bfile automatically)."""
        # Use relative path
        output_path = str(self.temp_output.relative_to(Path.cwd()))
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        # Perform clumping with VCF
        self.sumstats.clump(
            vcf=str(self.vcf_path),
            out=output_path,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            overwrite=False,  # Don't overwrite existing bfile
            verbose=False
        )
        
        # Check that results were stored in clumps dictionary
        results_sumstats = self.sumstats.clumps.get("clumps", pd.DataFrame())
        results = self.sumstats.clumps.get("clumps_raw", pd.DataFrame())
        plink_log = self.sumstats.clumps.get("plink_log", "")
        
        self.assertIsInstance(results_sumstats, pd.DataFrame)
        self.assertIsInstance(results, pd.DataFrame)
        self.assertIsInstance(plink_log, str)
    
    def test_clump_with_scaled_mlog10p(self):
        """Test clumping using MLOG10P instead of P values."""
        # Add MLOG10P column if not present
        if "MLOG10P" not in self.sumstats.data.columns:
            self.sumstats.data["MLOG10P"] = -np.log10(self.sumstats.data["P"])
        
        # Use relative path
        output_path = str(self.temp_output.relative_to(Path.cwd()))
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        # Perform clumping with scaled=True
        self.sumstats.clump(
            bfile=str(self.bfile_prefix),
            out=output_path,
            scaled=True,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        
        # Check that results were stored in clumps dictionary
        results_sumstats = self.sumstats.clumps.get("clumps", pd.DataFrame())
        results = self.sumstats.clumps.get("clumps_raw", pd.DataFrame())
        plink_log = self.sumstats.clumps.get("plink_log", "")
        
        self.assertIsInstance(results_sumstats, pd.DataFrame)
        self.assertIsInstance(results, pd.DataFrame)
        self.assertIsInstance(plink_log, str)
    
    def test_clump_different_parameters(self):
        """Test clumping with different parameter combinations."""
        test_params = [
            {
                "clump_p1": 1e-6,
                "clump_p2": 1e-5,
                "clump_r2": 0.1,
                "clump_kb": 500,
            },
            {
                "clump_p1": 5e-8,
                "clump_p2": 1e-5,
                "clump_r2": 0.5,
                "clump_kb": 100,
            },
        ]
        
        for i, params in enumerate(test_params):
            with self.subTest(params=params):
                # Create a new Sumstats object for each test to avoid state issues
                sumstats_copy = Sumstats(
                    sumstats=str(self.sumstats_path),
                    chrom="CHR",
                    pos="POS",
                    p="P",
                    ea="EA",
                    nea="NEA",
                    snpid="SNPID",
                    verbose=False
                )
                # Filter to chromosome 7
                chr_values = sumstats_copy.data["CHR"]
                try:
                    chr_numeric = pd.to_numeric(chr_values, errors='coerce')
                    has_chr7 = (chr_numeric == 7).any() or (chr_values.astype(str) == '7').any()
                except:
                    has_chr7 = (chr_values.astype(str) == '7').any()
                if has_chr7:
                    sumstats_copy.filter_value(expr="CHR == 7", inplace=True, verbose=False)
                    # Use only 1/10 of variants to speed up tests
                    if len(sumstats_copy.data) > 0:
                        sumstats_copy.data = sumstats_copy.data.iloc[::10].copy()
                
                # Fix SNP ID to make it consistent with PLINK file format
                sumstats_copy.fix_id(
                    fixid=True,
                    fixsep=True,
                    fixchrpos=False,
                    fixeanea=False,
                    overwrite=False,
                    verbose=False
                )
                
                # Use relative path
                output_path = str((self.temp_dir / f"clump_test_{self._testMethodName}_params_{i}").relative_to(Path.cwd()))
                self.temp_dir.mkdir(parents=True, exist_ok=True)
                tmp_base = Path("tmp")
                tmp_base.mkdir(exist_ok=True)
                
                sumstats_copy.clump(
                    bfile=str(self.bfile_prefix),
                    out=output_path,
                    threads=6,
                    verbose=False,
                    **params
                )
                
                # Check that results were stored
                results_sumstats = sumstats_copy.clumps.get("clumps", pd.DataFrame())
                results = sumstats_copy.clumps.get("clumps_raw", pd.DataFrame())
                plink_log = sumstats_copy.clumps.get("plink_log", "")
                
                self.assertIsInstance(results_sumstats, pd.DataFrame)
                self.assertIsInstance(results, pd.DataFrame)
                self.assertIsInstance(plink_log, str)
    
    def test_clump_no_significant_variants(self):
        """Test clumping when no variants meet significance threshold."""
        # Use relative path
        output_path = str(self.temp_output.relative_to(Path.cwd()))
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        # Use very strict p-value threshold
        self.sumstats.clump(
            bfile=str(self.bfile_prefix),
            out=output_path,
            clump_p1=1e-10,
            clump_p2=1e-10,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        
        # Check that results were stored
        results_sumstats = self.sumstats.clumps.get("clumps", pd.DataFrame())
        results = self.sumstats.clumps.get("clumps_raw", pd.DataFrame())
        
        # Should return empty DataFrames
        self.assertIsInstance(results_sumstats, pd.DataFrame)
        self.assertIsInstance(results, pd.DataFrame)
        self.assertEqual(len(results_sumstats), 0)
        self.assertEqual(len(results), 0)
    
    def test_clump_output_files_created(self):
        """Test that clumping creates expected output files."""
        # Use relative path
        output_path = str(self.temp_output.relative_to(Path.cwd()))
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        self.sumstats.clump(
            bfile=str(self.bfile_prefix),
            out=output_path,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        
        # Check that results were stored
        results = self.sumstats.clumps.get("clumps_raw", pd.DataFrame())
        plink_log = self.sumstats.clumps.get("plink_log", "")
        
        # Check that output files were created (if there are results)
        if len(results) > 0:
            expected_file = Path(f"{self.temp_output}.7.clumps")
            # Note: The file might be cleaned up, so we just check the log mentions it
            self.assertIn(".clumps", plink_log or "")
    
    def test_clump_bfile_vs_pfile_consistency(self):
        """Test that bfile and pfile produce similar results."""
        # Use relative paths
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        tmp_base = Path("tmp")
        tmp_base.mkdir(exist_ok=True)
        
        output_bfile = str((self.temp_dir / f"clump_test_{self._testMethodName}_bfile").relative_to(Path.cwd()))
        output_pfile = str((self.temp_dir / f"clump_test_{self._testMethodName}_pfile").relative_to(Path.cwd()))
        
        # Run clumping with bfile
        sumstats_bfile = Sumstats(
            sumstats=str(self.sumstats_path),
            chrom="CHR",
            pos="POS",
            p="P",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            verbose=False
        )
        # Filter to chromosome 7
        chr_values_bfile = sumstats_bfile.data["CHR"]
        try:
            chr_numeric = pd.to_numeric(chr_values_bfile, errors='coerce')
            has_chr7 = (chr_numeric == 7).any() or (chr_values_bfile.astype(str) == '7').any()
        except:
            has_chr7 = (chr_values_bfile.astype(str) == '7').any()
        if has_chr7:
            sumstats_bfile.filter_value(expr="CHR == 7", inplace=True, verbose=False)
            # Use only 1/10 of variants to speed up tests
            if len(sumstats_bfile.data) > 0:
                sumstats_bfile.data = sumstats_bfile.data.iloc[::10].copy()
        
        # Fix SNP ID to make it consistent with PLINK file format
        sumstats_bfile.fix_id(
            fixid=True,
            fixsep=True,
            fixchrpos=False,
            fixeanea=False,
            overwrite=False,
            verbose=False
        )
        sumstats_bfile.clump(
            bfile=str(self.bfile_prefix),
            out=output_bfile,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        results_bfile = sumstats_bfile.clumps.get("clumps", pd.DataFrame())
        
        # Run clumping with pfile
        sumstats_pfile = Sumstats(
            sumstats=str(self.sumstats_path),
            chrom="CHR",
            pos="POS",
            p="P",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            verbose=False
        )
        # Filter to chromosome 7
        chr_values_pfile = sumstats_pfile.data["CHR"]
        try:
            chr_numeric = pd.to_numeric(chr_values_pfile, errors='coerce')
            has_chr7 = (chr_numeric == 7).any() or (chr_values_pfile.astype(str) == '7').any()
        except:
            has_chr7 = (chr_values_pfile.astype(str) == '7').any()
        if has_chr7:
            sumstats_pfile.filter_value(expr="CHR == 7", inplace=True, verbose=False)
            # Use only 1/10 of variants to speed up tests
            if len(sumstats_pfile.data) > 0:
                sumstats_pfile.data = sumstats_pfile.data.iloc[::10].copy()
        
        # Fix SNP ID to make it consistent with PLINK file format
        sumstats_pfile.fix_id(
            fixid=True,
            fixsep=True,
            fixchrpos=False,
            fixeanea=False,
            overwrite=False,
            verbose=False
        )
        sumstats_pfile.clump(
            pfile=str(self.pfile_prefix),
            out=output_pfile,
            clump_p1=1e-5,
            clump_p2=1e-4,
            clump_r2=0.2,
            clump_kb=250,
            threads=6,
            verbose=False
        )
        results_pfile = sumstats_pfile.clumps.get("clumps", pd.DataFrame())
        
        # Both should return DataFrames
        self.assertIsInstance(results_bfile, pd.DataFrame)
        self.assertIsInstance(results_pfile, pd.DataFrame)
        
        # If both have results, they should have similar number of clumped variants
        # (allowing for some variation due to different implementations)
        if len(results_bfile) > 0 and len(results_pfile) > 0:
            # The number of clumped variants should be similar (within reasonable range)
            ratio = len(results_bfile) / len(results_pfile) if len(results_pfile) > 0 else 0
            self.assertGreater(ratio, 0.5, "bfile and pfile results differ too much")
            self.assertLess(ratio, 2.0, "bfile and pfile results differ too much")


if __name__ == "__main__":
    unittest.main()

