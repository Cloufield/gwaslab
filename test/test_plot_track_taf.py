"""Tests for track panel taf text_y_offset (taf[4]) behavior."""

import os
import sys
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
FIXTURES_DIR = os.path.join(ROOT, "test", "fixtures")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
if FIXTURES_DIR not in sys.path:
    sys.path.insert(0, FIXTURES_DIR)

import gwaslab as gl
from gwaslab.viz.viz_plot_track import plot_track
from panel_demo_fixtures import setup_demo


class TestPlotTrackTaf(unittest.TestCase):
    def setUp(self):
        self.demo = setup_demo()
        self.region = self.demo["region"]

    def tearDown(self):
        plt.close("all")

    def test_text_y_offset_preserved_for_middle_labels(self):
        fig, ax = plt.subplots(figsize=(8, 3), dpi=80)
        plot_track(
            self.demo["gtf"],
            self.region,
            ax=ax,
            fig=fig,
            feature_type="gene",
            taf=[4, 0, 0.95, 1, 6],
            verbose=False,
        )
        ys = sorted(t.get_position()[1] for t in ax.texts if t.get_text())
        self.assertTrue(ys)
        for y in ys:
            self.assertAlmostEqual(y, 6.0, places=5)

    def test_ylim_includes_gene_labels(self):
        fig, ax = plt.subplots(figsize=(8, 3), dpi=80)
        plot_track(
            self.demo["gtf"],
            self.region,
            ax=ax,
            fig=fig,
            feature_type="gene",
            taf=[4, 0, 0.95, 1, 6],
            verbose=False,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = ax.transData.inverted()
        ymin, ymax = ax.get_ylim()
        for text in ax.texts:
            if not text.get_text().strip():
                continue
            bb = text.get_window_extent(renderer).transformed(inv)
            self.assertGreaterEqual(bb.y0, ymin)
            self.assertLessEqual(bb.y1, ymax)

    def test_text_y_offset_via_panel(self):
        panel = gl.Panel(
            "track",
            track_path=self.demo["gtf"],
            region=self.region,
            feature_type="gene",
            taf=[4, 0, 0.95, 1, 3],
            verbose=False,
        )
        fig, _ = gl.plot_panels(
            [panel],
            region=self.region,
            fig_kwargs={"figsize": (8, 3), "dpi": 80},
            verbose=False,
        )
        ys = sorted(t.get_position()[1] for t in fig.axes[0].texts if t.get_text())
        self.assertTrue(ys)
        for y in ys:
            self.assertAlmostEqual(y, 3.0, places=5)


if __name__ == "__main__":
    unittest.main()
