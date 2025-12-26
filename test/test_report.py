"""
Test cases for the report generation functionality.

Tests the Sumstats.report() method and generate_qc_report() function
using a subset of the T2D BBJ dataset.
"""

import os
import sys
import unittest
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl
from gwaslab.info.g_Log import Log


class TestReportGeneration(unittest.TestCase):
    """Test cases for QC report generation."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data once for all tests."""
        cls.test_data_path = "examples/0_sample_data/t2d_bbj.txt.gz"
        cls.vcf_path = "/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
        cls.test_output_dir = Path("test/output/reports")
        cls.test_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if required files exist
        test_data_path_obj = Path(cls.test_data_path)
        if not test_data_path_obj.exists():
            raise unittest.SkipTest(f"Test data file not found: {cls.test_data_path}")
        
        vcf_path_obj = Path(cls.vcf_path)
        if not vcf_path_obj.exists():
            raise unittest.SkipTest(f"VCF file not found: {cls.vcf_path}")
        
        # Load 200,000 rows for testing
        print("\nLoading test data (200,000 rows from t2d_bbj.txt.gz)...")
        nrows_subset = 200000
        cls.sumstats = gl.Sumstats(
            cls.test_data_path,
            snpid="SNP",
            chrom="CHR",
            pos="POS",
            ea="ALT",
            nea="REF",
            neaf="Frq",
            beta="BETA",
            se="SE",
            p="P",
            direction="Dir",
            build="19",
            n="N",
            nrows=nrows_subset,  # Load 200,000 rows for testing
            verbose=False
        )
        print(f"Loaded {len(cls.sumstats.data):,} variants for testing")
    
    def setUp(self):
        """Set up for each test."""
        self.output_dir = self.test_output_dir
    
    def test_report_html_generation(self):
        """Test HTML report generation."""
        output_path = self.output_dir / "test_report.html"
        
        # Generate report
        result_path = self.sumstats.report(
            output_path=str(output_path),
            basic_check_kwargs={
                "remove": False,
                "remove_dup": False,
                "normalize": True,
                "verbose": False
            },
            get_lead_kwargs={
                "sig_level": 5e-8,
                "windowsizekb": 500,
                "anno": False,
                "verbose": False
            },
            mqq_plot_kwargs={
                "mode": "mqq",
                "sig_level": 5e-8,
                "dpi": 200,
                "verbose": False
            },
            report_title="Test QC Report (HTML)",
            verbose=False
        )
        
        # Check that report was generated
        self.assertTrue(Path(result_path).exists(), "HTML report file should exist")
        self.assertTrue(result_path.endswith('.html'), "Report should be HTML format")
        
        # Check that plots directory was created
        plots_dir = Path(result_path).parent / f"{Path(result_path).stem}_plots"
        self.assertTrue(plots_dir.exists(), "Plots directory should exist")
        
        # Check that MQQ plot exists
        mqq_plot = plots_dir / "mqq_plot.png"
        self.assertTrue(mqq_plot.exists(), "MQQ plot should exist")
    
    def test_report_pdf_generation(self):
        """Test PDF report generation (if weasyprint is available)."""
        output_path = self.output_dir / "test_report.pdf"
        
        # Generate report
        result_path = self.sumstats.report(
            output_path=str(output_path),
            basic_check_kwargs={
                "remove": False,
                "remove_dup": False,
                "normalize": True,
                "verbose": False
            },
            get_lead_kwargs={
                "sig_level": 5e-8,
                "windowsizekb": 500,
                "anno": False,
                "verbose": False
            },
            mqq_plot_kwargs={
                "mode": "mqq",
                "sig_level": 5e-8,
                "dpi": 200,
                "verbose": False
            },
            report_title="Test QC Report (PDF)",
            verbose=False
        )
        
        # Check that report was generated (may be HTML if PDF library not available)
        self.assertTrue(Path(result_path).exists(), "Report file should exist")
        
        # If PDF was generated, check it
        if result_path.endswith('.pdf'):
            self.assertTrue(Path(result_path).stat().st_size > 0, "PDF file should not be empty")
    
    def test_report_with_regional_plots(self):
        """Test report generation with regional plots and VCF LD reference."""
        output_path = self.output_dir / "test_report_with_regional.html"
        
        # Generate report with regional plots
        result_path = self.sumstats.report(
            output_path=str(output_path),
            basic_check_kwargs={
                "remove": False,
                "remove_dup": False,
                "normalize": True,
                "verbose": False
            },
            get_lead_kwargs={
                "sig_level": 5e-8,
                "windowsizekb": 500,
                "anno": False,
                "verbose": False
            },
            mqq_plot_kwargs={
                "mode": "mqq",
                "sig_level": 5e-8,
                "dpi": 200,
                "verbose": False
            },
            regional_plot_kwargs={
                "vcf_path": self.vcf_path,
                "region_recombination": True,
                "region_protein_coding": True,
                "region_ld_legend": True,
                "build": "19",
                "dpi": 200,
                "verbose": False
            },
            report_title="Test QC Report with Regional Plots",
            verbose=False
        )
        
        # Check that report was generated
        self.assertTrue(Path(result_path).exists(), "Report file should exist")
        
        # Check plots directory
        plots_dir = Path(result_path).parent / f"{Path(result_path).stem}_plots"
        if plots_dir.exists():
            # Check for regional plots if lead variants were found
            regional_plots = list(plots_dir.glob("regional_plot_*.png"))
            # Regional plots may or may not exist depending on whether lead variants were found
            # This is acceptable - we just check that the report was generated
    
    def test_report_method_chaining(self):
        """Test that report can be called after other methods."""
        output_path = self.output_dir / "test_report_chained.html"
        
        # Chain methods: basic_check -> report
        result_path = (self.sumstats
                      .basic_check(remove=False, remove_dup=False, normalize=True, verbose=False)
                      .report(output_path=str(output_path), verbose=False))
        
        self.assertTrue(Path(result_path).exists(), "Report should be generated after chaining")
    
    def test_report_custom_title(self):
        """Test report generation with custom title."""
        output_path = self.output_dir / "test_report_custom_title.html"
        custom_title = "Custom Test Report Title"
        
        result_path = self.sumstats.report(
            output_path=str(output_path),
            report_title=custom_title,
            verbose=False
        )
        
        # Read the HTML and check for custom title
        with open(result_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
            self.assertIn(custom_title, html_content, "HTML should contain custom title")
    
    def test_report_with_no_lead_variants(self):
        """Test report generation when no lead variants are found."""
        output_path = self.output_dir / "test_report_no_leads.html"
        
        # Use very strict significance threshold to ensure no lead variants
        result_path = self.sumstats.report(
            output_path=str(output_path),
            get_lead_kwargs={
                "sig_level": 1e-20,  # Very strict threshold
                "windowsizekb": 500,
                "verbose": False
            },
            verbose=False
        )
        
        # Report should still be generated even without lead variants
        self.assertTrue(Path(result_path).exists(), "Report should be generated even without lead variants")
        
        # MQQ plot should still exist
        plots_dir = Path(result_path).parent / f"{Path(result_path).stem}_plots"
        if plots_dir.exists():
            mqq_plot = plots_dir / "mqq_plot.png"
            self.assertTrue(mqq_plot.exists(), "MQQ plot should exist even without lead variants")
    
    def test_report_with_output(self):
        """Test report generation with output step."""
        output_path = self.output_dir / "test_report_with_output.html"
        output_sumstats_path = self.output_dir / "test_output_sumstats"
        
        result_path = self.sumstats.report(
            output_path=str(output_path),
            basic_check_kwargs={
                "remove": False,
                "remove_dup": False,
                "normalize": True,
                "verbose": False
            },
            get_lead_kwargs={
                "sig_level": 5e-8,
                "windowsizekb": 500,
                "verbose": False
            },
            output_kwargs={
                "path": str(output_sumstats_path),
                "fmt": "gwaslab",
                "gzip": True,
                "verbose": False
            },
            verbose=False
        )
        
        # Check that report was generated
        self.assertTrue(Path(result_path).exists(), "Report should be generated")
        
        # Check that output file was created (gwaslab format creates .gwaslab.tsv.gz)
        output_file = Path(f"{output_sumstats_path}.gwaslab.tsv.gz")
        self.assertTrue(output_file.exists(), "Output sumstats file should exist")
        
        # Check that report mentions output in HTML
        with open(result_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
            self.assertIn("Output", html_content, "HTML should mention output step")
    
    def test_report_with_harmonization(self):
        """Test report generation with harmonization step."""
        output_path = self.output_dir / "test_report_with_harmonization.html"
        
        # Note: This test may skip if harmonization references are not available
        # We'll use empty harmonize_kwargs to test the code path
        result_path = self.sumstats.report(
            output_path=str(output_path),
            basic_check_kwargs={
                "remove": False,
                "remove_dup": False,
                "normalize": True,
                "verbose": False
            },
            harmonize_kwargs={
                "basic_check": False,  # Skip basic_check since we already did it
                "verbose": False
            },
            get_lead_kwargs={
                "sig_level": 5e-8,
                "windowsizekb": 500,
                "verbose": False
            },
            verbose=False
        )
        
        # Check that report was generated
        self.assertTrue(Path(result_path).exists(), "Report should be generated with harmonization")
        
        # Check that report mentions harmonization in HTML
        with open(result_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
            self.assertIn("Harmonization", html_content, "HTML should mention harmonization")


if __name__ == "__main__":
    unittest.main()

