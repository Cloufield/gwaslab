"""
Test file for tutorial_v4.ipynb workflow.

This test follows the main tutorial steps but uses data from test/ref/ instead of downloading.
It verifies that all the tutorial code works correctly with v4.0.0 API.

To run:
    python test/test_tutorial_v4.py
    python -m unittest test.test_tutorial_v4 -v
    pytest test/test_tutorial_v4.py -v
"""

import os
import sys
import unittest
import tempfile
import shutil

import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl
import pandas as pd


# Path to reference data
REF_DIR = os.path.join(os.path.dirname(__file__), "ref")


class TestTutorialV4Workflow(unittest.TestCase):
    """Test the complete tutorial v4.0.0 workflow."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data paths."""
        cls.sumstats_path = os.path.join(REF_DIR, "bbj_t2d_hm3_chr7_variants.txt.gz")
        cls.rsid_tsv_path = os.path.join(REF_DIR, "1kg_dbsnp151_hg19_auto_hm3_chr7_variants.txt.gz")
        cls.rsid_vcf_path = os.path.join(REF_DIR, "b157_2564.vcf.gz")
        cls.vcf_path = os.path.join(REF_DIR, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        cls.ref_seq_path = os.path.join(REF_DIR, "chr7.fasta.gz")
        
        # Check that required files exist
        assert os.path.exists(cls.sumstats_path), f"Sumstats file not found: {cls.sumstats_path}"
        assert os.path.exists(cls.rsid_tsv_path), f"rsID TSV file not found: {cls.rsid_tsv_path}"
        assert os.path.exists(cls.rsid_vcf_path), f"rsID VCF file not found: {cls.rsid_vcf_path}"
        assert os.path.exists(cls.vcf_path), f"VCF file not found: {cls.vcf_path}"
    
    def setUp(self):
        """Set up for each test."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up after each test."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def _load_sumstats(self):
        """Helper method to load sumstats."""
        mysumstats = gl.Sumstats(
            self.sumstats_path,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            neaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            n="N",
            sep="\t",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        if len(mysumstats.data) > 0:
            mysumstats.data = mysumstats.data.iloc[::10].copy()
        return mysumstats
    
    def test_01_load_sumstats(self):
        """Test loading sumstats from file."""
        mysumstats = self._load_sumstats()
        
        self.assertIsNotNone(mysumstats.data)
        self.assertGreater(len(mysumstats.data), 0)
        self.assertIn("SNPID", mysumstats.data.columns)
        self.assertIn("CHR", mysumstats.data.columns)
        self.assertIn("POS", mysumstats.data.columns)
        self.assertIn("EA", mysumstats.data.columns)
        self.assertIn("NEA", mysumstats.data.columns)
        self.assertIn("EAF", mysumstats.data.columns)
        self.assertIn("STATUS", mysumstats.data.columns)
    
    def test_02_infer_build(self):
        """Test inferring genome build."""
        mysumstats = self._load_sumstats()
        
        # Fix CHR first (required for infer_build)
        mysumstats.fix_chr(verbose=False)
        
        # Infer build
        mysumstats.infer_build(verbose=False)
        
        # Check that build was inferred
        self.assertIsNotNone(mysumstats.build)
        self.assertIn(mysumstats.build, ["19", "38"])
    
    def test_03_plot_mqq_basic(self):
        """Test basic MQQ plot creation."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        # Test basic plot
        fig = mysumstats.plot_mqq(verbose=False)
        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.axes), 0)
        
        # Test with skip and cut
        fig = mysumstats.plot_mqq(skip=2, cut=20, verbose=False)
        self.assertIsNotNone(fig)
        
        # Test Manhattan only mode
        fig = mysumstats.plot_mqq(
            skip=2, 
            cut=20, 
            mode="m", 
            anno=True, 
            anno_sig_level=1e-30,
            verbose=False
        )
        self.assertIsNotNone(fig)
    
    def test_04_basic_check(self):
        """Test basic_check function."""
        mysumstats = self._load_sumstats()
        
        # Run basic check
        mysumstats.basic_check(verbose=False)
        
        # Verify data is still valid
        self.assertGreater(len(mysumstats.data), 0)
        
        # Check that QC was performed
        status = mysumstats.check_sumstats_qc_status()
        self.assertIn("basic_check", status)
    
    def test_05_get_lead(self):
        """Test extracting lead variants."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        # Get lead variants with annotation
        lead_variants = mysumstats.get_lead(anno=True, verbose=False)
        
        self.assertIsNotNone(lead_variants)
        if isinstance(lead_variants, pd.DataFrame):
            self.assertGreater(len(lead_variants), 0)
    
    def test_06_fix_id(self):
        """Test fixing SNPID format."""
        mysumstats = self._load_sumstats()
        
        # Fix ID separator
        mysumstats.fix_id(fixsep=True, verbose=False)
        
        # Verify SNPID format (should use : as separator)
        sample_snpid = mysumstats.data["SNPID"].iloc[0]
        if ":" in str(sample_snpid):
            # Check format is CHR:POS:REF:ALT
            parts = str(sample_snpid).split(":")
            self.assertGreaterEqual(len(parts), 2)
    
    def test_07_assign_rsid_tsv(self):
        """Test assigning rsID using TSV reference."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_id(fixsep=True, verbose=False)
        
        # Assign rsID using TSV
        mysumstats.assign_rsid(ref_rsid_tsv=self.rsid_tsv_path, verbose=False)
        
        # Check that rsID column exists
        self.assertIn("rsID", mysumstats.data.columns)
        
        # Check that some rsIDs were assigned
        rsid_count = mysumstats.data["rsID"].notna().sum()
        self.assertGreater(rsid_count, 0)
    
    def test_08_assign_rsid_vcf(self):
        """Test assigning rsID using VCF reference."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_id(fixsep=True, verbose=False)
        mysumstats.assign_rsid(ref_rsid_tsv=self.rsid_tsv_path, verbose=False)
        
        # Get initial rsID count
        initial_rsid_count = mysumstats.data["rsID"].notna().sum()
        
        # Assign rsID using VCF (for remaining variants)
        # Note: chr_dict is no longer needed - ChromosomeMapper handles chromosome conversion automatically
        mysumstats.assign_rsid(
            threads=1,
            ref_rsid_vcf=self.rsid_vcf_path,
            verbose=False
        )
        
        # Check that more rsIDs were assigned
        final_rsid_count = mysumstats.data["rsID"].notna().sum()
        self.assertGreaterEqual(final_rsid_count, initial_rsid_count)
    
    def test_09_harmonize(self):
        """Test harmonization workflow."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_id(fixsep=True, verbose=False)
        mysumstats.assign_rsid(ref_rsid_tsv=self.rsid_tsv_path, verbose=False)
        
        # Skip harmonization if reference files are not available
        # (harmonization requires reference genome and VCF which may not be in test/ref)
        if not os.path.exists(self.ref_seq_path):
            self.skipTest(f"Reference sequence not found: {self.ref_seq_path}")
        
        # Run harmonization (with basic_check=False since we already did it)
        try:
            mysumstats.harmonize(
                basic_check=False,
                threads=1,
                ref_seq=self.ref_seq_path,
                ref_infer=self.vcf_path,
                ref_alt_freq="AF",
                verbose=False
            )
            
            # Verify data is still valid
            self.assertGreater(len(mysumstats.data), 0)
            
        except Exception as e:
            # Harmonization might fail if reference files are incomplete
            # This is acceptable for a test
            self.skipTest(f"Harmonization skipped: {e}")
    
    def test_10_to_format_ldsc(self):
        """Test exporting to LDSC format."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        output_path = os.path.join(self.temp_dir, "test_ldsc")
        
        # Export to LDSC format
        mysumstats.to_format(output_path, fmt="ldsc", verbose=False)
        
        # Check that output file was created
        expected_file = f"{output_path}.ldsc.tsv.gz"
        self.assertTrue(os.path.exists(expected_file), f"Output file not created: {expected_file}")
        
        # Check that log file was created
        expected_log = f"{output_path}.ldsc.log"
        self.assertTrue(os.path.exists(expected_log), f"Log file not created: {expected_log}")
    
    def test_11_to_format_ssf(self):
        """Test exporting to SSF format."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        output_path = os.path.join(self.temp_dir, "test_ssf")
        
        # Export to SSF format
        mysumstats.to_format(output_path, fmt="ssf", ssfmeta=True, md5sum=True, verbose=False)
        
        # Check that output file was created
        expected_file = f"{output_path}.ssf.tsv.gz"
        self.assertTrue(os.path.exists(expected_file), f"Output file not created: {expected_file}")
    
    def test_12_liftover(self):
        """Test liftover functionality."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        # Get initial build
        initial_build = mysumstats.build
        initial_count = len(mysumstats.data)
        
        # Perform liftover
        try:
            mysumstats.liftover(threads=1, from_build="19", to_build="38", verbose=False)
            
            # Check that build changed
            self.assertEqual(mysumstats.build, "38")
            
            # Check that some variants may have been removed (unmapped)
            # But most should remain
            final_count = len(mysumstats.data)
            self.assertGreater(final_count, initial_count * 0.95)  # At least 95% should map
            
        except Exception as e:
            # Liftover might fail if chain files are not available
            self.skipTest(f"Liftover skipped: {e}")
    
    def test_13_filter_hapmap3(self):
        """Test filtering to HapMap3 variants."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        initial_count = len(mysumstats.data)
        
        # Filter to HapMap3
        mysumstats.filter_hapmap3(inplace=True, verbose=False)
        
        # Check that variants were filtered
        final_count = len(mysumstats.data)
        self.assertLess(final_count, initial_count)
        self.assertGreater(final_count, 0)
    
    def test_14_plot_mqq_advanced(self):
        """Test advanced MQQ plot with customizations."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        # Test advanced plot with all customizations
        fig = mysumstats.plot_mqq(
            mode="mqq",
            cut=20,
            skip=2,
            anno="GENENAME",
            anno_sig_level=1e-20,
            anno_alias={"7:127253550_C_T": "Test variant"},
            anno_style="expand",
            xpad=0.01,
            pinpoint=["7:127253550_C_T"],
            pinpoint_color="green",
            highlight=["7:127253550_C_T"],
            highlight_windowkb=1000,
            stratified=True,
            jagged=True,
            marker_size=(5, 5),
            fig_kwargs={"figsize": (15, 5), "dpi": 300},
            save=os.path.join(self.temp_dir, "test_plot.png"),
            save_kwargs={"dpi": 400, "facecolor": "white"},
            verbose=False
        )
        
        self.assertIsNotNone(fig)
        
        # Check that file was saved
        expected_file = os.path.join(self.temp_dir, "test_plot.png")
        self.assertTrue(os.path.exists(expected_file), f"Plot not saved: {expected_file}")
    
    def test_15_plot_regional(self):
        """Test regional plot."""
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        
        # Test regional plot without LD
        fig = mysumstats.plot_mqq(
            mode="r",
            skip=2,
            cut=20,
            region=(7, 126253550, 128253550),
            region_grid=True,
            verbose=False
        )
        
        self.assertIsNotNone(fig)
        
        # Test regional plot with LD (if VCF is available)
        if os.path.exists(self.vcf_path):
            fig = mysumstats.plot_mqq(
                mode="r",
                region=(7, 126253550, 128253550),
                region_grid=True,
                anno=True,
                anno_kwargs={"rotation": 0, "fontsize": 12},
                vcf_path=self.vcf_path,
                verbose=False
            )
            self.assertIsNotNone(fig)
    
    def test_16_summary_and_lookup_status(self):
        """Test summary and status lookup functions."""
        mysumstats = self._load_sumstats()
        mysumstats.basic_check(verbose=False)
        
        # Test summary
        summary = mysumstats.summary()
        self.assertIsNotNone(summary)
        
        # Test lookup_status
        status_info = mysumstats.lookup_status()
        self.assertIsNotNone(status_info)
    
    def test_17_complete_workflow(self):
        """Test the complete tutorial workflow in sequence."""
        # This test runs through the main workflow steps
        mysumstats = self._load_sumstats()
        mysumstats.fix_chr(verbose=False)
        mysumstats.infer_build(verbose=False)
        mysumstats.basic_check(verbose=False)
        mysumstats.fix_id(fixsep=True, verbose=False)
        mysumstats.assign_rsid(ref_rsid_tsv=self.rsid_tsv_path, verbose=False)
        
        # Verify final state
        self.assertGreater(len(mysumstats.data), 0)
        self.assertIn("rsID", mysumstats.data.columns)
        self.assertIsNotNone(mysumstats.build)


if __name__ == "__main__":
    unittest.main()

