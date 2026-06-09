"""Tests for AlphaGenome track bundles and ag_* panels."""

import os
import sys
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.viz.viz_aux_panel import Panel
from gwaslab.viz.viz_aux_track_bundle import (
    OverlayBundle,
    TrackBundle,
)
from gwaslab.viz.viz_plot_stackedpanel import plot_panels


def _mock_track_bundle(n_bins=100, n_tracks=2, chr=1, start=1_000_000):
    end = start + n_bins * 10
    values = np.random.RandomState(0).random((n_bins, n_tracks))
    meta = pd.DataFrame({
        "name": [f"track_{i}" for i in range(n_tracks)],
        "strand": ["+"] * n_tracks,
        "biosample_name": ["lung"] * n_tracks,
    })
    return TrackBundle(
        values=values,
        metadata=meta,
        resolution=10,
        region=(chr, start, end),
    )


def _mock_sumstats(n=50, chr=1, start=1_000_000):
    rng = np.random.RandomState(42)
    pos = [start + i * 1000 for i in range(n)]
    df = pd.DataFrame({
        "CHR": [chr] * n,
        "POS": pos,
        "SNPID": [f"rs{i}" for i in range(n)],
        "rsID": [f"rs{i}" for i in range(n)],
        "EA": ["A"] * n,
        "NEA": ["G"] * n,
        "P": rng.random(n) * 0.01,
        "MLOG10P": -np.log10(rng.random(n) * 0.01),
    })
    return df


class TestTrackBundle(unittest.TestCase):
    def test_num_axes(self):
        b = _mock_track_bundle(n_tracks=3)
        self.assertEqual(b.num_axes, 3)


class TestPanelSubplotCount(unittest.TestCase):
    def test_ag_tracks_subplot_count(self):
        b = _mock_track_bundle(n_tracks=2)
        p = Panel("ag_tracks", bundle=b, region=b.region, verbose=False)
        self.assertEqual(p.get_subplot_count(), 2)

    def test_deferred_defaults_one(self):
        p = Panel("ag_tracks", ag_spec={"output": "RNA_SEQ"}, region=(1, 1, 2), verbose=False)
        self.assertTrue(p.is_deferred())
        self.assertEqual(p.get_subplot_count(), 1)


class TestPlotAgPanels(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_plot_ag_tracks_panel(self):
        region = (1, 1_000_000, 1_001_000)
        b = _mock_track_bundle(n_tracks=2, start=region[1])
        panel = Panel("ag_tracks", bundle=b, region=region, verbose=False)
        fig, axes = plot_panels(
            [panel],
            region=region,
            verbose=False,
            fig_kwargs={"figsize": (6, 4), "dpi": 80},
        )
        self.assertIsNotNone(fig)
        self.assertEqual(len(axes), 2)

    def test_variant_positions_line(self):
        region = (1, 1_000_000, 1_001_000)
        b = _mock_track_bundle(n_tracks=1, start=region[1])
        panel = Panel("ag_tracks", bundle=b, region=region, verbose=False)
        fig, axes = plot_panels(
            [panel],
            region=region,
            variant_positions=[1_000_500],
            verbose=False,
            fig_kwargs={"figsize": (6, 3), "dpi": 80},
        )
        self.assertIsNotNone(fig)


if __name__ == "__main__":
    unittest.main()
