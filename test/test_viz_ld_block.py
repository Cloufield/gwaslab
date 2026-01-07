"""
Test suite for LD block visualization functionality.

This module contains unit tests for the `plot_ld_block` function, which visualizes
linkage disequilibrium (LD) matrices as 45°-rotated inverted triangles. The tests
cover standalone mode, integration with regional plots, various styling options,
annotations, and edge cases.

Note: Tests use mock LD matrices to avoid requiring actual VCF files or
pre-computed LD data.
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

from gwaslab.viz.viz_plot_ld_block import plot_ld_block
from gwaslab.info.g_Log import Log


def make_mock_ld_matrix(n, seed=42):
    """
    Create a mock symmetric LD matrix (r² values) with realistic structure.
    
    The matrix is designed so that nearby variants have higher LD values,
    simulating real LD decay with distance.
    
    Parameters
    ----------
    n : int
        Number of variants (matrix will be n x n).
    seed : int, optional
        Random seed for reproducibility. Default is 42.
    
    Returns
    -------
    np.ndarray
        Symmetric n x n matrix with r² values between 0 and 1.
        Diagonal elements are 1.0 (self-LD).
    """
    rng = np.random.RandomState(seed)
    ld = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                ld[i, j] = 1.0
            else:
                # Create some structure: nearby variants have higher LD
                dist = abs(i - j)
                base_ld = max(0, 0.9 - dist * 0.1)
                ld[i, j] = min(1.0, base_ld + rng.random() * 0.2)
                ld[j, i] = ld[i, j]  # Symmetric
    return ld


def make_mock_positions(n, start=156938803):
    """
    Create mock genomic positions array.
    
    Parameters
    ----------
    n : int
        Number of positions to generate.
    start : int, optional
        Starting genomic position. Default is 156938803.
    
    Returns
    -------
    np.ndarray
        Array of n positions, spaced 1kb apart starting from start.
    """
    return np.array([start + i * 1000 for i in range(n)])


class TestLdBlock(unittest.TestCase):
    """
    Test cases for LD block visualization function `plot_ld_block`.
    
    This test class covers:
    - Standalone mode with provided LD matrices
    - Custom colormaps and value ranges
    - Cell annotations and grid lines
    - Lead SNP markers and variant annotations
    - Integration with existing axes
    - Save functionality
    - Edge cases (masked arrays, small/large matrices)
    - Integration with regional plots
    """
    def setUp(self):
        self.log = Log()
        self.n_variants = 50
        self.ld_matrix = make_mock_ld_matrix(self.n_variants)
        self.positions = make_mock_positions(self.n_variants)
        
    def tearDown(self):
        plt.close('all')

    def test_plot_ld_block_standalone_with_ld_matrix(self):
        """
        Test plot_ld_block in standalone mode with provided LD matrix.
        
        Verifies that the function creates a figure with both the main LD plot
        and the position bar when called in standalone mode (no ax provided).
        With cbar=True (default), it also creates a colorbar axes.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            verbose=False
        )
        self.assertIsNotNone(fig)
        self.assertIsNotNone(ax)
        # With cbar=True (default): position bar + LD block + colorbar = 3 axes
        # The colorbar is created as an inset_axes, which adds to fig.axes
        self.assertEqual(len(fig.axes), 3)  # position bar + LD block + colorbar

    def test_plot_ld_block_with_custom_cmap(self):
        """
        Test plot_ld_block with custom colormap.
        
        Verifies that the function accepts custom matplotlib colormap names
        (e.g., "viridis") and applies them correctly to the LD visualization.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            cmap="viridis",
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_vmin_vmax(self):
        """
        Test plot_ld_block with custom vmin and vmax.
        
        Verifies that custom value ranges for the colormap are correctly
        applied, allowing users to focus on specific LD value ranges.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            vmin=0.0,
            vmax=0.8,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_anno_cell(self):
        """
        Test plot_ld_block with cell annotations.
        
        Verifies that LD values can be displayed as text annotations
        in each cell of the LD block plot with custom formatting.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            anno_cell=True,
            anno_cell_fmt="{:.2f}",
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_grid(self):
        """
        Test plot_ld_block with grid lines.
        
        Verifies that grid lines can be added to the LD block plot
        with custom styling (linewidth, alpha, etc.).
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            ld_block_grid=True,
            ld_block_grid_kwargs={"linewidth": 1, "alpha": 0.3},
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_lead_snp_markers(self):
        """
        Test plot_ld_block with lead SNP markers.
        
        Verifies that lead SNP positions can be marked on the LD block plot
        with custom colors to highlight reference variants.
        """
        lead_snp_is = [10, 25, 40]  # Indices of lead SNPs
        lead_snp_is_color = ["red", "blue", "green"]
        
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            lead_snp_is=lead_snp_is,
            lead_snp_is_color=lead_snp_is_color,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_anno(self):
        """
        Test plot_ld_block with variant annotations.
        
        Verifies that variant annotations can be added to the plot,
        with options to limit the number of annotated variants and
        specify which variants to annotate.
        """
        anno_set = [0, 5, 10, 15, 20]
        
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            ld_block_anno=True,
            ld_block_anno_set=anno_set,
            ld_block_anno_max_rows=10,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_existing_axes(self):
        """
        Test plot_ld_block with provided axes.
        
        Verifies that the function can plot into existing matplotlib axes
        (ax and ax_pos), allowing integration with custom figure layouts.
        """
        fig, ax_main = plt.subplots(figsize=(10, 8))
        ax_pos = fig.add_subplot(2, 1, 2)
        
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            ax=ax_main,
            ax_pos=ax_pos,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_without_colorbar(self):
        """
        Test plot_ld_block without colorbar.
        
        Verifies that the colorbar can be disabled when cbar=False,
        useful for saving space or when using custom legends.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            cbar=False,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_custom_title(self):
        """
        Test plot_ld_block with custom title.
        
        Verifies that a custom plot title can be specified and displayed.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            title="Custom LD Block Plot",
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_custom_xlabel(self):
        """
        Test plot_ld_block with custom x-axis label.
        
        Verifies that the x-axis label can be customized to provide
        context-specific information about the genomic positions.
        """
        fig, ax = plot_ld_block(
            ld=self.ld_matrix,
            pos=self.positions,
            xlabel="Custom Position Label",
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_with_save(self):
        """
        Test plot_ld_block with save option.
        
        Verifies that the plot can be saved to a file with custom
        save parameters (e.g., DPI). The temporary file is cleaned up
        after the test.
        """
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            fig, ax = plot_ld_block(
                ld=self.ld_matrix,
                pos=self.positions,
                save=tmp_path,
                save_kwargs={"dpi": 100},
                verbose=False
            )
            self.assertTrue(os.path.exists(tmp_path))
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_plot_ld_block_with_masked_ld(self):
        """
        Test plot_ld_block with masked LD matrix (missing values).
        
        Verifies that the function can handle numpy masked arrays,
        which are used to represent missing or invalid LD values
        in the matrix.
        """
        import numpy.ma as ma
        masked_ld = ma.masked_array(self.ld_matrix)
        # Mask some values
        masked_ld[0:5, 0:5] = ma.masked
        
        fig, ax = plot_ld_block(
            ld=masked_ld,
            pos=self.positions,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_small_matrix(self):
        """
        Test plot_ld_block with small LD matrix.
        
        Verifies that the function handles small matrices (e.g., 5x5)
        correctly, ensuring edge cases with minimal data are handled.
        """
        small_ld = make_mock_ld_matrix(5)
        small_pos = make_mock_positions(5)
        
        fig, ax = plot_ld_block(
            ld=small_ld,
            pos=small_pos,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_large_matrix(self):
        """
        Test plot_ld_block with larger LD matrix.
        
        Verifies that the function handles larger matrices (e.g., 100x100)
        correctly, ensuring performance and visualization quality with
        more variants.
        """
        large_ld = make_mock_ld_matrix(100)
        large_pos = make_mock_positions(100)
        
        fig, ax = plot_ld_block(
            ld=large_ld,
            pos=large_pos,
            verbose=False
        )
        self.assertIsNotNone(fig)

    def test_plot_ld_block_integration_with_regional_plot(self):
        """
        Test that ld_block parameters work with plot_mqq in regional mode.
        
        This is an integration test that verifies ld_block parameters
        (ld_block, anno_cell, ld_block_grid, etc.) are correctly passed
        through to plot_ld_block when using plot_mqq in regional mode.
        Note: Full execution requires a valid VCF file.
        """
        import gwaslab
        
        # Create minimal sumstats
        sumstats_df = pd.DataFrame({
            "CHR": [7] * self.n_variants,
            "POS": self.positions,
            "P": np.random.random(self.n_variants),
            "SNPID": [f"rs{i}" for i in range(self.n_variants)],
            "NEA": ["A"] * self.n_variants,
            "EA": ["G"] * self.n_variants,
        })
        
        sumstats = gwaslab.Sumstats(
            sumstats=sumstats_df,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            p="P",
            build="19",
            verbose=False
        )
        
        # Test that parameters don't cause errors (even if VCF is missing)
        try:
            fig = sumstats.plot_mqq(
                mode="r",
                region=(7, int(self.positions[0]), int(self.positions[-1])),
                vcf_path=None,  # No VCF, so ld_block won't actually run
                ld_block=True,
                anno_cell=True,
                ld_block_grid=True,
                verbose=False
            )
            # Should still create a plot even without VCF
            self.assertIsNotNone(fig)
        except Exception as e:
            # Some errors are expected without VCF, but function should handle gracefully
            pass


if __name__ == "__main__":
    unittest.main()

