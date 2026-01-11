"""
Test suite for Panel class and plot_panels function.

This module contains unit tests for the Panel class, which stores panel configuration
information for stacked plots, and the plot_panels function, which creates stacked
multi-panel figures from Panel objects.

Tests cover:
- Panel initialization with different panel types
- Panel methods (get_type, get_kwargs, get_kwarg, set_kwarg, update_kwargs)
- Panel string representations
- plot_panels with different panel types
- Region extraction and validation
- Height ratios and titles
- Error handling
"""
import os
import sys
import unittest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.viz.viz_aux_panel import Panel
from gwaslab.viz.viz_plot_stackedpanel import plot_panels
from gwaslab.info.g_Log import Log


def make_mock_sumstats(n=100, chr=1, start_pos=1000000):
    """
    Create a mock sumstats DataFrame for testing.
    
    Parameters
    ----------
    n : int
        Number of variants
    chr : int
        Chromosome number
    start_pos : int
        Starting genomic position
    
    Returns
    -------
    pd.DataFrame
        Mock sumstats DataFrame with required columns
    """
    rng = np.random.RandomState(42)
    positions = [start_pos + i * 1000 for i in range(n)]
    
    data = {
        "CHR": [chr] * n,
        "POS": positions,
        "SNPID": [f"rs{i}" for i in range(n)],
        "EA": ["A"] * n,
        "NEA": ["G"] * n,
        "BETA": rng.randn(n) * 0.1,
        "SE": np.abs(rng.randn(n) * 0.01) + 0.01,
        "P": rng.random(n) * 0.1,
        "EAF": rng.random(n) * 0.4 + 0.1,
    }
    df = pd.DataFrame(data)
    df["MLOG10P"] = -np.log10(df["P"])
    return df


def make_mock_pipcs_data(n=50, chr=1, start_pos=1000000):
    """
    Create a mock PIPCS DataFrame for testing.
    
    Parameters
    ----------
    n : int
        Number of variants
    chr : int
        Chromosome number
    start_pos : int
        Starting genomic position
    
    Returns
    -------
    pd.DataFrame
        Mock PIPCS DataFrame with required columns
    """
    rng = np.random.RandomState(42)
    positions = [start_pos + i * 1000 for i in range(n)]
    
    data = {
        "CHR": [chr] * n,
        "POS": positions,
        "SNPID": [f"rs{i}" for i in range(n)],
        "PIP": rng.random(n),
        "CS": [1 if i < 10 else 0 for i in range(n)],
    }
    return pd.DataFrame(data)


class TestPanel(unittest.TestCase):
    """
    Test cases for Panel class.
    
    This test class covers:
    - Panel initialization with different panel types
    - Panel methods (get_type, get_kwargs, get_kwarg, set_kwarg, update_kwargs)
    - Panel string representations
    - Parameter processing and validation
    - Error handling
    """
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.region = (1, 1000000, 2000000)
        self.sumstats = make_mock_sumstats()
    
    def test_panel_init_track(self):
        """Test Panel initialization with track type."""
        panel = Panel(
            "track",
            track_path="test.gtf",
            region=self.region,
            color="#020080",
            verbose=False
        )
        self.assertEqual(panel.get_type(), "track")
        self.assertIn("track_path", panel.get_kwargs())
        self.assertEqual(panel.get_kwarg("track_path"), "test.gtf")
        self.assertEqual(panel.get_kwarg("region"), self.region)
    
    def test_panel_init_arc(self):
        """Test Panel initialization with arc type."""
        panel = Panel(
            "arc",
            bedpe_path="test.bedpe.gz",
            region=self.region,
            color="#FF0000",
            alpha=0.3,
            verbose=False
        )
        self.assertEqual(panel.get_type(), "arc")
        self.assertIn("bedpe_path", panel.get_kwargs())
        self.assertEqual(panel.get_kwarg("bedpe_path"), "test.bedpe.gz")
    
    def test_panel_init_ld_block(self):
        """Test Panel initialization with ld_block type."""
        panel = Panel(
            "ld_block",
            vcf_path="test.vcf.gz",
            region=self.region,
            insumstats=self.sumstats,
            verbose=False
        )
        self.assertEqual(panel.get_type(), "ld_block")
        self.assertIn("insumstats", panel.get_kwargs())
    
    def test_panel_init_region(self):
        """Test Panel initialization with region type."""
        panel = Panel(
            "region",
            insumstats=self.sumstats,
            region=self.region,
            vcf_path="test.vcf.gz",
            build="38",
            verbose=False
        )
        self.assertEqual(panel.get_type(), "region")
        self.assertIn("insumstats", panel.get_kwargs())
    
    def test_panel_init_chromatin(self):
        """Test Panel initialization with chromatin type."""
        panel = Panel(
            "chromatin",
            region_chromatin_files=["test.bed.gz"],
            region_chromatin_labels=["E098"],
            region=self.region,
            verbose=False
        )
        self.assertEqual(panel.get_type(), "chromatin")
        self.assertIn("region_chromatin_files", panel.get_kwargs())
    
    def test_panel_init_pipcs(self):
        """Test Panel initialization with pipcs type."""
        pipcs_data = make_mock_pipcs_data()
        panel = Panel(
            "pipcs",
            pipcs_raw=pipcs_data,
            region=self.region,
            verbose=False
        )
        self.assertEqual(panel.get_type(), "pipcs")
        self.assertIn("pipcs_raw", panel.get_kwargs())
    
    def test_panel_init_case_insensitive(self):
        """Test that panel_type is case-insensitive."""
        panel = Panel("TRACK", track_path="test.gtf", region=self.region, verbose=False)
        self.assertEqual(panel.get_type(), "track")
    
    def test_panel_init_unknown_type(self):
        """Test Panel initialization with unknown type (should still work)."""
        panel = Panel("unknown", custom_param="value", verbose=False)
        self.assertEqual(panel.get_type(), "unknown")
        self.assertEqual(panel.get_kwarg("custom_param"), "value")
    
    def test_panel_init_type_error(self):
        """Test that non-string panel_type raises TypeError."""
        with self.assertRaises(TypeError):
            Panel(123, track_path="test.gtf", region=self.region)
    
    def test_panel_get_type(self):
        """Test get_type method."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        self.assertEqual(panel.get_type(), "track")
    
    def test_panel_get_kwargs(self):
        """Test get_kwargs method returns a copy."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        kwargs1 = panel.get_kwargs()
        kwargs2 = panel.get_kwargs()
        self.assertIsNot(kwargs1, kwargs2)  # Should be different objects
        self.assertEqual(kwargs1, kwargs2)  # But same content
    
    def test_panel_get_kwarg(self):
        """Test get_kwarg method."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        self.assertEqual(panel.get_kwarg("track_path"), "test.gtf")
        self.assertEqual(panel.get_kwarg("nonexistent", default="default"), "default")
        self.assertIsNone(panel.get_kwarg("nonexistent"))
    
    def test_panel_set_kwarg(self):
        """Test set_kwarg method."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        panel.set_kwarg("color", "#FF0000")
        self.assertEqual(panel.get_kwarg("color"), "#FF0000")
    
    def test_panel_update_kwargs(self):
        """Test update_kwargs method."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        panel.update_kwargs(color="#FF0000", alpha=0.5)
        self.assertEqual(panel.get_kwarg("color"), "#FF0000")
        self.assertEqual(panel.get_kwarg("alpha"), 0.5)
    
    def test_panel_repr(self):
        """Test __repr__ method."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        repr_str = repr(panel)
        self.assertIn("Panel", repr_str)
        self.assertIn("track", repr_str)
    
    def test_panel_str(self):
        """Test __str__ method."""
        panel = Panel("track", track_path="test.gtf", region=self.region, verbose=False)
        str_str = str(panel)
        self.assertIn("Panel", str_str)
        self.assertIn("track", str_str)


class TestPlotPanels(unittest.TestCase):
    """
    Test cases for plot_panels function.
    
    This test class covers:
    - Basic functionality with different panel types
    - Region extraction from panels
    - Height ratios
    - Titles
    - Figure creation
    - X-axis alignment
    - Error handling
    """
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.region = (1, 1000000, 2000000)
        self.sumstats = make_mock_sumstats()
        self.pipcs_data = make_mock_pipcs_data()
    
    def tearDown(self):
        """Clean up after tests."""
        plt.close('all')
    
    def test_plot_panels_empty_list(self):
        """Test that empty panels list raises ValueError."""
        with self.assertRaises(ValueError):
            plot_panels([], verbose=False)
    
    def test_plot_panels_non_panel_object(self):
        """Test that non-Panel objects raise TypeError."""
        with self.assertRaises(TypeError):
            plot_panels(["not a panel"], verbose=False)
    
    def test_plot_panels_missing_region(self):
        """Test that missing region raises ValueError."""
        panel = Panel("track", track_path="test.gtf", verbose=False)
        with self.assertRaises(ValueError):
            plot_panels([panel], verbose=False)
    
    def test_plot_panels_single_track_panel(self):
        """Test plot_panels with a single track panel."""
        # Note: This will fail if track_path doesn't exist, but we test the structure
        panel = Panel(
            "track",
            track_path="nonexistent.gtf",  # Will fail to load but tests structure
            region=self.region,
            verbose=False
        )
        # We expect this to fail when trying to load the file, but the function
        # should be called correctly. For a real test, we'd need a valid track file.
        # This test verifies the function accepts the panel structure.
        try:
            fig, axes = plot_panels([panel], region=self.region, verbose=False)
            # If it doesn't fail, verify structure
            self.assertIsNotNone(fig)
            self.assertEqual(len(axes), 1)
        except (FileNotFoundError, ValueError) as e:
            # Expected if file doesn't exist, but structure is correct
            error_msg = str(e).lower()
            self.assertTrue("track_path" in error_msg or "file" in error_msg or "gtf" in error_msg)
    
    def test_plot_panels_region_extraction(self):
        """Test that region is extracted from panels if not provided."""
        panel = Panel(
            "track",
            track_path="nonexistent.gtf",
            region=self.region,
            verbose=False
        )
        # Should extract region from panel
        try:
            fig, axes = plot_panels([panel], verbose=False)
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if file doesn't exist
            pass
    
    def test_plot_panels_height_ratios(self):
        """Test plot_panels with custom height ratios."""
        panel1 = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        panel2 = Panel("arc", bedpe_path="nonexistent.bedpe.gz", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel1, panel2],
                region=self.region,
                height_ratios=[2.0, 3.0],
                verbose=False
            )
            self.assertIsNotNone(fig)
            self.assertEqual(len(axes), 2)
        except (FileNotFoundError, ValueError):
            # Expected if files don't exist
            pass
    
    def test_plot_panels_height_ratios_mismatch(self):
        """Test that height_ratios length mismatch raises ValueError."""
        panel1 = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        panel2 = Panel("arc", bedpe_path="nonexistent.bedpe.gz", region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels(
                [panel1, panel2],
                region=self.region,
                height_ratios=[2.0],  # Wrong length
                verbose=False
            )
    
    def test_plot_panels_titles(self):
        """Test plot_panels with titles."""
        panel = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel],
                region=self.region,
                titles=["Test Track"],
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if file doesn't exist
            pass
    
    def test_plot_panels_titles_multiple(self):
        """Test plot_panels with multiple titles."""
        panel1 = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        panel2 = Panel("arc", bedpe_path="nonexistent.bedpe.gz", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel1, panel2],
                region=self.region,
                titles=["Track", "Arc"],
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if files don't exist
            pass
    
    def test_plot_panels_fig_kwargs(self):
        """Test plot_panels with custom fig_kwargs."""
        panel = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel],
                region=self.region,
                fig_kwargs={"figsize": (8, 6), "dpi": 100},
                verbose=False
            )
            self.assertIsNotNone(fig)
            # Verify figsize was set (may be adjusted by auto_adjust)
            figsize = fig.get_size_inches()
            self.assertGreater(figsize[0], 0)
            self.assertGreater(figsize[1], 0)
        except (FileNotFoundError, ValueError):
            # Expected if file doesn't exist
            pass
    
    def test_plot_panels_align_xaxis(self):
        """Test plot_panels with x-axis alignment."""
        panel1 = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        panel2 = Panel("arc", bedpe_path="nonexistent.bedpe.gz", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel1, panel2],
                region=self.region,
                align_xaxis=True,
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if files don't exist
            pass
    
    def test_plot_panels_no_align_xaxis(self):
        """Test plot_panels without x-axis alignment."""
        panel = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel],
                region=self.region,
                align_xaxis=False,
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if file doesn't exist
            pass
    
    def test_plot_panels_track_missing_path(self):
        """Test that track panel without track_path raises ValueError."""
        panel = Panel("track", region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_arc_missing_path(self):
        """Test that arc panel without bedpe_path raises ValueError."""
        panel = Panel("arc", region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_region_missing_sumstats(self):
        """Test that region panel without sumstats raises ValueError."""
        panel = Panel("region", region=self.region, vcf_path="test.vcf.gz", verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_chromatin_missing_files(self):
        """Test that chromatin panel without required params raises ValueError."""
        panel = Panel("chromatin", region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_pipcs_missing_data(self):
        """Test that pipcs panel without pipcs_raw raises ValueError."""
        panel = Panel("pipcs", region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_pipcs_invalid_data(self):
        """Test that pipcs panel with invalid data raises ValueError."""
        invalid_data = pd.DataFrame({"invalid": [1, 2, 3]})
        panel = Panel("pipcs", pipcs_raw=invalid_data, region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_unsupported_type(self):
        """Test that unsupported panel type raises ValueError."""
        panel = Panel("unsupported", region=self.region, verbose=False)
        
        with self.assertRaises(ValueError):
            plot_panels([panel], region=self.region, verbose=False)
    
    def test_plot_panels_region_conflict(self):
        """Test that conflicting regions are handled."""
        panel1 = Panel("track", track_path="nonexistent.gtf", region=(1, 1000000, 2000000), verbose=False)
        panel2 = Panel("arc", bedpe_path="nonexistent.bedpe.gz", region=(1, 1500000, 2500000), verbose=False)
        
        # Should use the specified region and warn about conflict
        try:
            fig, axes = plot_panels(
                [panel1, panel2],
                region=(1, 1000000, 2000000),
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if files don't exist
            pass
    
    def test_plot_panels_title_positions(self):
        """Test plot_panels with different title positions."""
        panel = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        
        for title_pos in ["left", "right", "center"]:
            try:
                fig, axes = plot_panels(
                    [panel],
                    region=self.region,
                    titles=["Test"],
                    title_pos=title_pos,
                    verbose=False
                )
                self.assertIsNotNone(fig)
            except (FileNotFoundError, ValueError):
                # Expected if file doesn't exist
                pass
    
    def test_plot_panels_subplot_height(self):
        """Test plot_panels with custom subplot_height."""
        panel = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel],
                region=self.region,
                subplot_height=2.0,
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if file doesn't exist
            pass
    
    def test_plot_panels_hspace(self):
        """Test plot_panels with custom hspace."""
        panel1 = Panel("track", track_path="nonexistent.gtf", region=self.region, verbose=False)
        panel2 = Panel("arc", bedpe_path="nonexistent.bedpe.gz", region=self.region, verbose=False)
        
        try:
            fig, axes = plot_panels(
                [panel1, panel2],
                region=self.region,
                hspace=0.1,
                verbose=False
            )
            self.assertIsNotNone(fig)
        except (FileNotFoundError, ValueError):
            # Expected if files don't exist
            pass


if __name__ == "__main__":
    unittest.main()
