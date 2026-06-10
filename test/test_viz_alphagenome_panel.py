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
from gwaslab.viz.viz_aux_track_bundle import OverlayBundle, TrackBundle, bundle_num_axes
from gwaslab.viz.viz_plot_stackedpanel import plot_panels
from gwaslab.viz.viz_plot_alphagenome import (
    DEFAULT_AG_YLABEL_TEMPLATE,
    _default_ag_ylabel,
    _modality_from_row,
    _style_ag_panel_spines,
    _ylabel_from_metadata,
    finalize_ag_panel_spines,
    plot_ag_overlay,
)


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


def _mock_overlay_bundle(n_bins=100, n_tracks=1, start=1_000_000, alt_scale=1.2):
    end = start + n_bins * 10
    meta = pd.DataFrame({
        "name": ["RNA-seq"] * n_tracks,
        "strand": ["+"] * n_tracks,
        "biosample_name": ["lung"] * n_tracks,
    })
    base = np.linspace(0.2, 0.8, n_bins)
    ref_tb = TrackBundle(
        values=base.reshape(-1, n_tracks),
        metadata=meta,
        resolution=10,
        region=(1, start, end),
    )
    alt_tb = TrackBundle(
        values=(base * alt_scale).reshape(-1, n_tracks),
        metadata=meta.copy(),
        resolution=10,
        region=(1, start, end),
    )
    return OverlayBundle(
        tracks={"REF": ref_tb, "ALT": alt_tb},
        region=(1, start, end),
        variant_pos=start + 500,
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


class TestAgSpines(unittest.TestCase):
    def test_style_ag_panel_spines_full_box(self):
        _, ax = plt.subplots()
        _style_ag_panel_spines(ax)
        for side in ("top", "right", "left", "bottom"):
            self.assertTrue(ax.spines[side].get_visible())
        plt.close("all")

    def test_finalize_ag_panel_spines_all_axes(self):
        fig, axes = plt.subplots(2, 1)
        finalize_ag_panel_spines(axes)
        for ax in axes:
            for side in ("top", "right", "left", "bottom"):
                self.assertTrue(ax.spines[side].get_visible())
        plt.close(fig)


class TestAgYlabel(unittest.TestCase):
    def test_default_three_line_ylabel(self):
        row = {
            "biosample_name": "lung",
            "ontology_curie": "UBERON:0002048",
            "name": "polyA plus RNA-seq",
            "strand": "+",
        }
        self.assertEqual(
            _default_ag_ylabel(row),
            "lung\nUBERON:0002048\npolyA plus RNA-seq",
        )

    def test_default_ylabel_skips_missing_fields(self):
        row = {"biosample_name": "lung", "name": "ATAC", "strand": "+"}
        self.assertEqual(_default_ag_ylabel(row), "lung\nATAC")

    def test_modality_prefers_assay_title_over_name(self):
        row = {
            "name": "UBERON:0002107 total RNA-seq",
            "Assay title": "total RNA-seq",
            "ontology_curie": "UBERON:0002107",
            "biosample_name": "liver",
        }
        self.assertEqual(_modality_from_row(row), "total RNA-seq")
        self.assertEqual(
            _default_ag_ylabel(row),
            "liver\nUBERON:0002107\ntotal RNA-seq",
        )

    def test_modality_strips_ontology_from_name_fallback(self):
        row = {
            "name": "UBERON:0002048 polyA plus RNA-seq",
            "ontology_curie": "UBERON:0002048",
            "biosample_name": "lung",
        }
        self.assertEqual(_modality_from_row(row), "polyA plus RNA-seq")
        self.assertEqual(
            _default_ag_ylabel(row),
            "lung\nUBERON:0002048\npolyA plus RNA-seq",
        )

    def test_ylabel_from_metadata_uses_template(self):
        meta = pd.DataFrame({
            "biosample_name": ["lung"],
            "ontology_curie": ["UBERON:0002048"],
            "name": ["RNA-seq"],
            "strand": ["+"],
        })
        ylab = _ylabel_from_metadata(meta, 0, DEFAULT_AG_YLABEL_TEMPLATE)
        self.assertEqual(ylab, "lung\nUBERON:0002048\nRNA-seq")


class TestPanelSubplotCount(unittest.TestCase):
    def test_ag_tracks_subplot_count(self):
        b = _mock_track_bundle(n_tracks=2)
        p = Panel("ag_tracks", bundle=b, region=b.region, verbose=False)
        self.assertEqual(p.get_subplot_count(), 2)

    def test_ag_overlay_subplot_count_splits_ref_alt(self):
        bundle = _mock_overlay_bundle(n_tracks=2)
        self.assertEqual(bundle_num_axes(bundle), 4)
        p = Panel("ag_overlay", bundle=bundle, region=bundle.region, verbose=False)
        self.assertEqual(p.get_subplot_count(), 4)

    def test_deferred_defaults_one(self):
        p = Panel("ag_tracks", ag_spec={"output": "RNA_SEQ"}, region=(1, 1, 2), verbose=False)
        self.assertTrue(p.is_deferred())
        self.assertEqual(p.get_subplot_count(), 1)

    def test_deferred_ag_overlay_defaults_two(self):
        p = Panel("ag_overlay", ag_spec={"output": "RNA_SEQ", "mode": "overlay"}, region=(1, 1, 2), verbose=False)
        self.assertTrue(p.is_deferred())
        self.assertEqual(p.get_subplot_count(), 2)


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

    def test_plot_ag_tracks_ylabel_fontsize(self):
        region = (1, 1_000_000, 1_001_000)
        meta = pd.DataFrame({
            "name": ["RNA-seq"],
            "strand": ["+"],
            "biosample_name": ["lung"],
            "ontology_curie": ["UBERON:0002048"],
        })
        values = np.random.RandomState(0).random((100, 1))
        b = TrackBundle(
            values=values,
            metadata=meta,
            resolution=10,
            region=(1, region[1], region[1] + 1000),
        )
        panel = Panel("ag_tracks", bundle=b, region=region, verbose=False)
        fig, axes = plot_panels(
            [panel],
            region=region,
            verbose=False,
            fig_kwargs={"figsize": (6, 3), "dpi": 80},
        )
        self.assertEqual(axes[0].get_ylabel(), "lung\nUBERON:0002048\nRNA-seq")
        self.assertEqual(axes[0].yaxis.label.get_size(), 9.0)

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

    def test_plot_ag_overlay_interleaves_ref_alt_per_modality(self):
        bundle = _mock_overlay_bundle(n_tracks=2)
        fig, axes = plt.subplots(4, 1, figsize=(6, 6))
        plot_ag_overlay(bundle, axes, verbose=False)
        alleles = [ax.get_ylabel().split("\n", 1)[0] for ax in axes]
        self.assertEqual(alleles, ["REF", "ALT", "REF", "ALT"])

    def test_plot_ag_overlay_draws_ref_and_alt_on_separate_axes(self):
        bundle = _mock_overlay_bundle()
        fig, axes = plt.subplots(2, 1, figsize=(6, 4))
        plot_ag_overlay(bundle, axes, verbose=False)
        self.assertEqual(len(axes[0].get_lines()), 2)
        self.assertEqual(len(axes[1].get_lines()), 2)
        self.assertTrue(axes[0].get_ylabel().startswith("REF\n"))
        self.assertTrue(axes[1].get_ylabel().startswith("ALT\n"))

    def test_plot_ag_overlay_panel(self):
        region = (1, 1_000_000, 1_001_000)
        bundle = _mock_overlay_bundle(start=region[1])
        panel = Panel("ag_overlay", bundle=bundle, region=region, verbose=False)
        fig, axes = plot_panels(
            [panel],
            region=region,
            verbose=False,
            fig_kwargs={"figsize": (6, 4), "dpi": 80},
        )
        self.assertEqual(len(axes), 2)
        self.assertEqual(len(axes[0].get_lines()), 2)
        self.assertEqual(len(axes[1].get_lines()), 2)


if __name__ == "__main__":
    unittest.main()
