import os
import sys
import tempfile
import unittest

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
SCRIPTS_DIR = os.path.join(ROOT, "scripts")
FIXTURES_DIR = os.path.join(ROOT, "test", "fixtures")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)
if FIXTURES_DIR not in sys.path:
    sys.path.insert(0, FIXTURES_DIR)

import gwaslab as gl
from gwaslab.g_Sumstats import Sumstats
from gwaslab.viz.viz_aux_sankey import (
    MAF_CATEGORY_ORDER,
    MAF_PALETTE,
    _categorize_maf,
    _categorize_p,
    assign_sankey_colors,
    maf_bin_counts,
    prepare_sankey_data,
)
from gwaslab.viz.viz_plot_sankey import plot_sankey


def _make_fixture(n: int = 100) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    eaf = rng.uniform(0.0005, 0.5, size=n)
    p = 10 ** (-rng.uniform(4, 12, size=n))
    beta = rng.normal(0, 0.2, size=n)
    return pd.DataFrame({"EAF": eaf, "P": p, "BETA": beta})


class TestSankeyAux(unittest.TestCase):
    def test_preset_maf_matches_summarize_bins(self):
        df = _make_fixture(200)
        expected = maf_bin_counts(df["EAF"])
        categorized = _categorize_maf(df["EAF"]).dropna()
        for label in MAF_CATEGORY_ORDER:
            self.assertEqual(int((categorized == label).sum()), expected[label])

    def test_build_links_maf_to_p(self):
        df = _make_fixture(100)
        maf = _categorize_maf(df["EAF"])
        p_cat = _categorize_p(df["P"])
        manual = (
            pd.DataFrame({"MAF": maf, "P": p_cat})
            .dropna()
            .groupby(["MAF", "P"])
            .size()
            .reset_index(name="value")
        )

        _, links_df, _, _, _, _ = prepare_sankey_data(df, columns=["MAF", "P"])
        built = links_df.copy()
        built["MAF"] = built["source"].str.split("|", n=1).str[1]
        built["P"] = built["target"].str.split("|", n=1).str[1]
        built = (
            built.groupby(["MAF", "P"], as_index=False)["value"]
            .sum()
            .sort_values(["MAF", "P"])
            .reset_index(drop=True)
        )
        manual = manual.sort_values(["MAF", "P"]).reset_index(drop=True)
        pd.testing.assert_frame_equal(built, manual.astype({"value": float}), check_dtype=False)

    def test_fine_links_have_color_key(self):
        df = pd.DataFrame(
            {
                "EAF": [0.2, 0.2, 0.003, 0.003],
                "P": [1e-10, 1e-4, 1e-10, 1e-4],
                "BETA": [0.2, 0.2, 0.2, 0.2],
            }
        )
        _, links_df, _, _, _, _ = prepare_sankey_data(df, columns=["MAF", "P", "BETA"])
        self.assertIn("color_key", links_df.columns)
        hop1 = links_df[links_df["stage_from"] == 1]
        dup_targets = hop1.groupby(["source", "target"]).size()
        self.assertTrue((dup_targets > 1).any())

    def test_maf_palette_order(self):
        df = _make_fixture(50)
        _, _, _, stage_names, work, _ = prepare_sankey_data(df, columns=["MAF", "P"])
        flow_colors, _ = assign_sankey_colors(stage_names, work, palette="auto")
        present = [c for c in MAF_CATEGORY_ORDER if c in flow_colors]
        if len(present) >= 2:
            idx0 = MAF_CATEGORY_ORDER.index(present[0])
            idx1 = MAF_CATEGORY_ORDER.index(present[1])
            pal0 = MAF_PALETTE[min(idx0, len(MAF_PALETTE) - 1)]
            pal1 = MAF_PALETTE[min(idx1, len(MAF_PALETTE) - 1)]
            self.assertEqual(flow_colors[present[0]], pal0)
            self.assertEqual(flow_colors[present[1]], pal1)

    def test_stacked_node_bands(self):
        df = _make_fixture(120)
        nodes_df, _, node_bands_df, _, _, _ = prepare_sankey_data(
            df, columns=["MAF", "P", "BETA"]
        )
        for node_id, node in nodes_df.set_index("node_id").iterrows():
            bands = node_bands_df[node_bands_df["node_id"] == node_id]
            if bands.empty:
                continue
            band_height = (bands["y1"] - bands["y0"]).sum()
            node_height = node["y1"] - node["y0"]
            self.assertAlmostEqual(band_height, node_height, places=5)

    def test_unknown_stage_raises(self):
        df = _make_fixture(10)
        with self.assertRaisesRegex(ValueError, "Unknown stage"):
            prepare_sankey_data(df, columns=["MAF", "not_a_real_stage"])


class TestPlotSankey(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_plot_sankey_smoke(self):
        df = _make_fixture(80)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "sankey.png")
            fig, ax, tables = plot_sankey(
                df,
                columns=["MAF", "P", "BETA"],
                title="test sankey",
                save=out,
                verbose=False,
            )
            self.assertIsNotNone(fig)
            self.assertIsNotNone(ax)
            self.assertGreater(os.path.getsize(out), 0)
            self.assertFalse(tables["nodes"].empty)
            self.assertFalse(tables["links"].empty)

    def test_sumstats_method(self):
        df = _make_fixture(60)
        ss = Sumstats(
            sumstats=df,
            chrom="CHR",
            pos="POS",
            p="P",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            verbose=False,
        )
        ss.data["CHR"] = 1
        ss.data["POS"] = np.arange(1, len(ss.data) + 1)
        ss.data["EA"] = "A"
        ss.data["NEA"] = "G"
        ss.data["SNPID"] = [f"rs{i}" for i in range(len(ss.data))]

        fig, ax, tables = ss.plot_sankey(columns=["MAF", "P"], verbose=False)
        self.assertIsNotNone(fig)
        self.assertIsNotNone(ax)
        self.assertEqual(tables["stages"], ["MAF", "P"])

    def test_gl_plot_sankey_wrapper(self):
        df = _make_fixture(40)
        fig, ax, tables = gl.plot_sankey(df, columns=["MAF", "P"], verbose=False)
        self.assertIsNotNone(fig)
        self.assertIsNotNone(ax)
        self.assertIn("nodes", tables)


class TestDiseaseSubtypeFixture(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_disease_subtype_sankey_example(self):
        from sankey_demo_data import (  # noqa: WPS433
            simulate_disease_subtype_sumstats,
            overall_signal_colors,
        )

        df = simulate_disease_subtype_sumstats(n_variants=500, seed=0)
        self.assertIn("overall_signal", df.columns)
        self.assertIn("subtype_signal", df.columns)
        self.assertGreater((df["overall_signal"] == "Overall GW significant").sum(), 0)

        fig, ax, tables = plot_sankey(
            df,
            columns=["overall_signal", "subtype_signal", "MAF"],
            colors=overall_signal_colors(),
            verbose=False,
        )
        self.assertIsNotNone(fig)
        self.assertFalse(tables["links"].empty)
        self.assertIn("Overall GW significant", tables["flow_colors"])


if __name__ == "__main__":
    unittest.main()
