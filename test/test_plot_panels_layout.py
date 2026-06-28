"""Layout behavior tests for plot_panels height_ratios."""

import os
import sys
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
FIXTURES_DIR = os.path.join(ROOT, "test", "fixtures")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
if FIXTURES_DIR not in sys.path:
    sys.path.insert(0, FIXTURES_DIR)

import gwaslab as gl
from gwaslab.viz.viz_plot_stackedpanel import (
    _add_panel_title,
    _apply_shared_font_defaults,
    _left_label_extent_fig_x,
)
from panel_demo_fixtures import setup_demo


def _axis_heights(fig):
    return [ax.get_position().height for ax in fig.axes[:8]]


class TestPlotPanelsLayout(unittest.TestCase):
    def setUp(self):
        self.demo = setup_demo()
        self.region = self.demo["region"]
        self.mysumstats = gl.Sumstats(
            self.demo["sumstats"],
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            build="38",
            verbose=False,
        )

    def tearDown(self):
        plt.close("all")

    def _genomic_panels(self):
        return [
            gl.Panel("track", track_path=self.demo["gtf"], region=self.region, feature_type="gene", verbose=False),
            gl.Panel("arc", bedpe_path=self.demo["bedpe"], region=self.region, verbose=False),
            gl.Panel(
                "chromatin",
                region_chromatin_files=[self.demo["chromatin_bed"]],
                region_chromatin_labels=["Simulated"],
                region=self.region,
                verbose=False,
            ),
            gl.Panel("pipcs", pipcs_raw=self.demo["pipcs"], region=self.region, verbose=False),
            gl.Panel(
                "ld_block",
                vcf_path=self.demo["vcf"],
                insumstats=self.demo["ld_sumstats"],
                region=self.region,
                tabix=None,
                verbose=False,
            ),
            self.mysumstats.Panel(
                "region",
                region=self.region,
                vcf_path=self.demo["vcf"],
                gtf_path=self.demo["gtf"],
                tabix=None,
                build="38",
                verbose=False,
            ),
        ]

    def test_single_panel_without_height_ratios(self):
        panels = self._genomic_panels()[:1]
        fig, axes = gl.plot_panels(
            panels,
            region=self.region,
            fig_kwargs={"figsize": (8, 3), "dpi": 80},
            verbose=False,
        )
        self.assertIsNotNone(fig)
        self.assertEqual(len(axes), 1)

    def test_user_figsize_width_not_capped(self):
        panels = self._genomic_panels()[:1]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            fig_kwargs={"figsize": (80, 10), "dpi": 80},
            verbose=False,
        )
        w, h = fig.get_size_inches()
        self.assertAlmostEqual(w, 80.0, places=3)
        self.assertAlmostEqual(h, 10.0, places=3)

    def test_custom_height_ratios_preserve_region_gene_ratio(self):
        panels = self._genomic_panels()
        ratios = [1.5, 1.5, 1.5, 1.5, 5.0, 10.0]
        fig, axes = gl.plot_panels(
            panels,
            region=self.region,
            height_ratios=ratios,
            hspace=0.12,
            fig_kwargs={"figsize": (12, 20), "dpi": 80},
            verbose=False,
        )
        self.assertEqual(len(axes), 8)
        heights = _axis_heights(fig)
        # region scatter (index 6) vs gene track (index 7): expanded 10 vs 5
        self.assertAlmostEqual(heights[6] / heights[7], 2.0, places=2)
        # fig size unchanged after ld_block + alignment
        w, h = fig.get_size_inches()
        self.assertAlmostEqual(w, 12.0, places=3)
        self.assertAlmostEqual(h, 20.0, places=3)

    def test_shared_font_defaults_respect_user_overrides(self):
        kw = {"fontsize": 14, "font_family": "Times"}
        _apply_shared_font_defaults(kw, "chromatin", 9, "Arial")
        self.assertEqual(kw["fontsize"], 14)
        self.assertEqual(kw["font_family"], "Times")

    def test_shared_font_defaults_inject_region_sizes(self):
        kw = {}
        _apply_shared_font_defaults(kw, "region", 9, "Arial")
        self.assertEqual(kw["anno_fontsize"], 9)
        self.assertEqual(kw["cbar_fontsize"], 9)
        self.assertEqual(kw["title_fontsize"], 9)
        self.assertEqual(kw["track_font_family"], "Arial")

    def test_left_title_in_margin_left_of_ylabel(self):
        fig, ax = plt.subplots()
        ax.set_ylabel("Signal", rotation=90)
        from gwaslab.viz.viz_plot_stackedpanel import _apply_panel_titles
        _apply_panel_titles(fig, [(ax, "Panel A")], "left", {"fontsize": 10, "family": "Arial"})
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        ylab_x = _left_label_extent_fig_x(ax, fig, renderer)
        titles = [t for t in fig.texts if t.get_text() == "Panel A"]
        self.assertEqual(len(titles), 1)
        self.assertEqual(titles[0].get_ha(), "left")
        t_bb = titles[0].get_window_extent(renderer).transformed(fig.transFigure.inverted())
        self.assertLessEqual(t_bb.x1, ylab_x)

    def test_plot_panels_margin_titles(self):
        panels = self._genomic_panels()[:1]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            titles=["Genes"],
            fig_kwargs={"figsize": (8, 3), "dpi": 80},
            verbose=False,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        titles = [t for t in fig.texts if t.get_text() == "Genes"]
        self.assertEqual(len(titles), 1)
        ax = fig.axes[0]
        ylab_x = _left_label_extent_fig_x(ax, fig, renderer)
        t_bb = titles[0].get_window_extent(renderer).transformed(fig.transFigure.inverted())
        self.assertLessEqual(t_bb.x1, ylab_x)

    def test_arc_panel_has_full_spines(self):
        panels = self._genomic_panels()[1:2]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            fig_kwargs={"figsize": (8, 3), "dpi": 80},
            verbose=False,
        )
        ax = fig.axes[0]
        for side in ("left", "right", "top", "bottom"):
            self.assertTrue(ax.spines[side].get_visible())

    def _assert_colorbar_inside_panel(self, fig, panel_ax, pad_frac=0.03):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        cbar_axes = [
            ax for ax in fig.axes
            if ax is not panel_ax and getattr(panel_ax, "_gwaslab_arc_cbar_axes", None) is ax
        ]
        self.assertEqual(len(cbar_axes), 1)
        cax = cbar_axes[0]
        panel_bb = panel_ax.get_position()
        pad_x = pad_frac * panel_bb.width
        pad_y = pad_frac * panel_bb.height
        tight = cax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())
        eps = 1e-6
        self.assertLessEqual(tight.y1, panel_bb.y1 - pad_y + eps)
        self.assertGreaterEqual(tight.y0, panel_bb.y0 + pad_y - eps)
        self.assertLessEqual(tight.x1, panel_bb.x1 - pad_x + eps)
        self.assertGreaterEqual(tight.x0, panel_bb.x0 + pad_x - eps)

    def test_arc_colorbar_stays_in_panel_after_alignment(self):
        panels = self._genomic_panels()[:2]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            fig_kwargs={"figsize": (8, 6), "dpi": 80},
            verbose=False,
        )
        self._assert_colorbar_inside_panel(fig, fig.axes[1])

    def test_arc_colorbar_stays_in_single_panel(self):
        panels = self._genomic_panels()[1:2]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            fig_kwargs={"figsize": (8, 3), "dpi": 80},
            verbose=False,
        )
        self._assert_colorbar_inside_panel(fig, fig.axes[0])

    def test_plot_panels_titles_left_aligned_column(self):
        panels = self._genomic_panels()[:2]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            titles=["Genes", "Arc"],
            fig_kwargs={"figsize": (8, 5), "dpi": 80},
            verbose=False,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        titles = [t for t in fig.texts if t.get_text() in ("Genes", "Arc")]
        self.assertEqual(len(titles), 2)
        x0s = [
            t.get_window_extent(renderer).transformed(fig.transFigure.inverted()).x0
            for t in titles
        ]
        self.assertAlmostEqual(x0s[0], x0s[1], places=3)
        for t, ax in zip(titles, fig.axes[:2]):
            ylab_x = _left_label_extent_fig_x(ax, fig, renderer)
            t_bb = t.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            self.assertLessEqual(t_bb.x1, ylab_x)


if __name__ == "__main__":
    unittest.main()
