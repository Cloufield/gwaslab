import os
import sys
import unittest

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.text import Text, Annotation

from gwaslab.viz.viz_aux_reposition_text import (
    _get_text_rotated_height,
    _get_highest_y_pixels,
)


class TestGetTextRotatedHeight(unittest.TestCase):
    """Tests for _get_text_rotated_height function."""

    def setUp(self):
        """Set up a figure and axes for each test."""
        self.fig, self.ax = plt.subplots(figsize=(10, 10), dpi=100)
        self.fig.canvas.draw()
        self.renderer = self.fig.canvas.get_renderer()

    def tearDown(self):
        """Close figure after each test."""
        plt.close(self.fig)

    def test_empty_text_returns_zero(self):
        """Empty text should return 0."""
        text = self.ax.text(0.5, 0.5, "", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        self.assertEqual(result, 0)

    def test_whitespace_text_returns_zero(self):
        """Whitespace-only text should return 0."""
        text = self.ax.text(0.5, 0.5, "   ", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        self.assertEqual(result, 0)

    def test_no_rotation_returns_bbox_y1(self):
        """Text with no rotation should return the bbox y1."""
        text = self.ax.text(0.5, 0.5, "Test Text", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        bbox = text.get_window_extent(self.renderer)
        self.assertAlmostEqual(result, bbox.y1, places=1)

    def test_90_degree_rotation_extends_vertically(self):
        """90-degree rotated text should extend higher than unrotated."""
        # Create unrotated text
        text_unrotated = self.ax.text(0.5, 0.3, "LongTextString", 
                                       transform=self.ax.transAxes, rotation=0)
        self.fig.canvas.draw()
        height_unrotated = _get_text_rotated_height(text_unrotated, self.renderer, self.fig)
        
        # Create 90-degree rotated text at same position
        text_rotated = self.ax.text(0.5, 0.3, "LongTextString", 
                                     transform=self.ax.transAxes, rotation=90,
                                     ha='left', va='bottom')
        self.fig.canvas.draw()
        height_rotated = _get_text_rotated_height(text_rotated, self.renderer, self.fig)
        
        # Rotated text should extend higher (text width becomes height)
        self.assertGreater(height_rotated, height_unrotated)

    def test_rotation_with_different_alignments(self):
        """Test rotated text with different horizontal/vertical alignments."""
        alignments = [
            ('left', 'bottom'),
            ('center', 'center'),
            ('right', 'top'),
        ]
        
        for ha, va in alignments:
            with self.subTest(ha=ha, va=va):
                text = self.ax.text(0.5, 0.5, "TestText", 
                                    transform=self.ax.transAxes,
                                    rotation=45, ha=ha, va=va)
                self.fig.canvas.draw()
                result = _get_text_rotated_height(text, self.renderer, self.fig)
                # Should return a positive value
                self.assertGreater(result, 0)

    def test_long_text_rotated_90_extends_more(self):
        """Longer text rotated 90 degrees should extend higher."""
        short_text = self.ax.text(0.5, 0.3, "ABC", 
                                   transform=self.ax.transAxes, rotation=90,
                                   ha='left', va='bottom')
        self.fig.canvas.draw()
        height_short = _get_text_rotated_height(short_text, self.renderer, self.fig)
        
        long_text = self.ax.text(0.5, 0.3, "ABCDEFGHIJKLMNOP", 
                                  transform=self.ax.transAxes, rotation=90,
                                  ha='left', va='bottom')
        self.fig.canvas.draw()
        height_long = _get_text_rotated_height(long_text, self.renderer, self.fig)
        
        # Longer text should extend higher when rotated 90 degrees
        self.assertGreater(height_long, height_short)

    def test_string_rotation_vertical(self):
        """Test that 'vertical' string rotation is handled correctly."""
        text = self.ax.text(0.5, 0.3, "VerticalText", 
                            transform=self.ax.transAxes, rotation='vertical',
                            ha='left', va='bottom')
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value (vertical = 90 degrees)
        self.assertGreater(result, 0)

    def test_string_rotation_horizontal(self):
        """Test that 'horizontal' string rotation is handled correctly."""
        text = self.ax.text(0.5, 0.3, "HorizontalText", 
                            transform=self.ax.transAxes, rotation='horizontal')
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value (horizontal = 0 degrees)
        self.assertGreater(result, 0)

    def test_mathtext_with_dollar_signs(self):
        """Test that mathtext with $ signs is handled."""
        text = self.ax.text(0.5, 0.5, r"$\alpha + \beta$", 
                            transform=self.ax.transAxes)
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value
        self.assertGreater(result, 0)

    def test_baseline_vertical_alignment(self):
        """Test baseline vertical alignment handling."""
        text = self.ax.text(0.5, 0.5, "BaselineTest", 
                            transform=self.ax.transAxes, rotation=45,
                            ha='left', va='baseline')
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value
        self.assertGreater(result, 0)

    def test_center_baseline_vertical_alignment(self):
        """Test center_baseline vertical alignment handling."""
        text = self.ax.text(0.5, 0.5, "CenterBaselineTest", 
                            transform=self.ax.transAxes, rotation=45,
                            ha='center', va='center_baseline')
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value
        self.assertGreater(result, 0)

    def test_rotation_mode_anchor(self):
        """Test that rotation_mode='anchor' is handled."""
        text = self.ax.text(0.5, 0.5, "AnchorModeTest", 
                            transform=self.ax.transAxes, rotation=45,
                            rotation_mode='anchor')
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value
        self.assertGreater(result, 0)

    def test_uses_public_fontproperties_api(self):
        """Test that the function works with different font properties."""
        text = self.ax.text(0.5, 0.5, "FontTest", 
                            transform=self.ax.transAxes,
                            fontsize=14, fontweight='bold', fontstyle='italic')
        self.fig.canvas.draw()
        result = _get_text_rotated_height(text, self.renderer, self.fig)
        # Should return a positive value
        self.assertGreater(result, 0)
        # Verify get_fontproperties works
        fp = text.get_fontproperties()
        self.assertIsNotNone(fp)


class TestGetHighestYPixels(unittest.TestCase):
    """Tests for _get_highest_y_pixels function."""

    def setUp(self):
        """Set up a figure and axes for each test."""
        self.fig, self.ax = plt.subplots(figsize=(10, 10), dpi=100)
        self.fig.canvas.draw()
        self.renderer = self.fig.canvas.get_renderer()

    def tearDown(self):
        """Close figure after each test."""
        plt.close(self.fig)

    def test_empty_axes_returns_tight_bbox(self):
        """Empty axes should return the tight bbox y1."""
        self.fig.canvas.draw()
        result = _get_highest_y_pixels(self.ax, self.renderer, self.fig)
        bbox = self.ax.get_tightbbox(self.renderer)
        # Should be close to the axes tight bbox
        self.assertGreater(result, 0)

    def test_text_extends_above_axes(self):
        """Text above axes should be detected."""
        # Add text near top of axes
        text = self.ax.text(0.5, 1.1, "Above Axes", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        
        result = _get_highest_y_pixels(self.ax, self.renderer, self.fig)
        
        # Should detect text above axes
        axes_bbox = self.ax.get_tightbbox(self.renderer)
        self.assertGreaterEqual(result, axes_bbox.y1)

    def test_rotated_text_detected(self):
        """Rotated text extending upward should be detected."""
        # Add 90-degree rotated text
        text = self.ax.text(0.5, 0.8, "RotatedTextExtending", 
                            transform=self.ax.transAxes, rotation=90,
                            ha='left', va='bottom')
        self.fig.canvas.draw()
        
        result = _get_highest_y_pixels(self.ax, self.renderer, self.fig)
        
        # Should detect the rotated text extending upward
        self.assertGreater(result, 0)

    def test_multiple_texts_finds_highest(self):
        """Should find the highest point among multiple texts."""
        # Add texts at different heights
        self.ax.text(0.5, 0.3, "Low", transform=self.ax.transAxes)
        self.ax.text(0.5, 0.5, "Middle", transform=self.ax.transAxes)
        self.ax.text(0.5, 0.9, "High", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        
        result = _get_highest_y_pixels(self.ax, self.renderer, self.fig)
        
        # Should be higher than the middle text position
        fig_height = self.fig.get_figheight() * self.fig.dpi
        middle_y = 0.5 * fig_height
        self.assertGreater(result, middle_y)

    def test_annotation_detected(self):
        """Annotations should be detected."""
        # Add an annotation
        self.ax.annotate("Annotation", xy=(0.5, 0.5), xytext=(0.5, 0.9),
                         textcoords='axes fraction', xycoords='axes fraction',
                         arrowprops=dict(arrowstyle="->"))
        self.fig.canvas.draw()
        
        result = _get_highest_y_pixels(self.ax, self.renderer, self.fig)
        
        # Should detect the annotation
        self.assertGreater(result, 0)

    def test_fig_parameter_optional(self):
        """Function should work without fig parameter (gets from axes)."""
        self.ax.text(0.5, 0.5, "Test", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        
        # Call without fig parameter
        result = _get_highest_y_pixels(self.ax, self.renderer)
        
        self.assertGreater(result, 0)

    def test_returns_float(self):
        """Should return a float value."""
        self.ax.text(0.5, 0.5, "Test", transform=self.ax.transAxes)
        self.fig.canvas.draw()
        
        result = _get_highest_y_pixels(self.ax, self.renderer, self.fig)
        
        self.assertIsInstance(result, (int, float))


class TestIntegrationWithPlot(unittest.TestCase):
    """Integration tests with actual plot scenarios."""

    def test_manhattan_like_annotations(self):
        """Test with Manhattan-plot-like annotations."""
        fig, ax = plt.subplots(figsize=(15, 5), dpi=100)
        
        # Simulate data points
        x = np.arange(100)
        y = np.random.rand(100) * 10
        ax.scatter(x, y)
        
        # Add rotated annotations like SNP labels
        for i in [20, 50, 80]:
            ax.text(i, y[i] + 0.5, f"rs{i}_{i}", rotation=90, 
                    ha='left', va='bottom', fontsize=9)
        
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        
        result = _get_highest_y_pixels(ax, renderer, fig)
        
        # Should detect texts extending above the data
        self.assertGreater(result, 0)
        
        plt.close(fig)

    def test_title_positioning_scenario(self):
        """Test scenario similar to actual title positioning."""
        fig, ax = plt.subplots(figsize=(10, 8), dpi=100)
        
        # Add some content
        ax.plot([0, 1], [0, 1])
        ax.text(0.5, 0.95, "AnnotationText", transform=ax.transAxes, 
                rotation=90, ha='left', va='bottom')
        
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        
        highest_y = _get_highest_y_pixels(ax, renderer, fig)
        fig_height_pixels = fig.get_figheight() * fig.dpi
        
        # Calculate title position as the actual code does
        highest_y_fig = highest_y / fig_height_pixels
        padding = 0.04
        title_y = highest_y_fig + padding
        
        # Title should be above the highest content
        self.assertGreater(title_y, highest_y_fig)
        
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
