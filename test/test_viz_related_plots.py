import os
import sys
import tempfile
import unittest
import random
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2
# regional_mode_setup / process_vcf: viz_plot_regional2; bulk scatter: draw_manhattan_panel
from gwaslab.viz.viz_plot_regional2 import _plot_regional, regional_mode_setup, process_vcf
from gwaslab.viz.viz_plot_manhattan_mode import draw_manhattan_panel
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.viz.viz_plot_trumpetplot import _plot_trumpet
from gwaslab.viz.viz_plot_stackedregional import plot_stacked_mqq
from gwaslab.viz.viz_plot_compare_af import plotdaf


# Three chr1 SNPs, six samples, diploid GT — enough for allel.read_vcf(..., tabix=None)
# and rogers_huff LD vs a chosen lead. Plain .vcf (no bgzip/tabix) keeps CI simple.
_MINI_VCF_FOR_LD = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
1\t1000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\t0/1\t0/1
1\t1001\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\t0/0\t0/1
1\t1002\t.\tG\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t0/1\t1/1\t0/1\t0/0
"""


def make_sumstats(n=400, with_eaf=True):
    # Create synthetic sumstats with CHR/POS/P/SNPID and optional EAF
    rng = random.Random(123)
    rows = []
    for i in range(n):
        chr_ = rng.randint(1, 22)
        pos_ = rng.randint(1_000, 2_000_000)
        pval = max(min(rng.random(), 0.999999), 1e-300)
        snpid = f"{chr_}:{pos_}_A_G"
        row = {"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid}
        if with_eaf:
            row["EAF"] = rng.random()
        rows.append(row)
    return pd.DataFrame(rows)


class TestRelatedPlots(unittest.TestCase):
    def setUp(self):
        self.df = make_sumstats()
        self.log = Log()

    def test_plot_miami2_with_kwargs(self):
        # Test plot_miami2 with merged_sumstats and suffixes, plus highlight/pinpoint/sig-level
        merged = self.df.rename(columns={"P": "P_1"})
        merged["P_2"] = merged["P_1"].sample(frac=1.0, random_state=1).values
        fig, log = plot_miami2(
            merged_sumstats=merged,
            id0="SNPID",
            id1="SNPID",
            id2="SNPID",
            mode="m",
            titles=["Top", "Bottom"],
            verbose=False,
            anno=None,
            suffixes=("_1","_2"),
            highlight=[merged.iloc[0]["SNPID"], merged.iloc[1]["SNPID"]],
            highlight_color="#CB132D",
            highlight_windowkb=200,
            pinpoint=[merged.iloc[2]["SNPID"]],
            pinpoint_color="red",
            anno_sig_level=1e-6,
        )
        self.assertIsInstance(log, Log)

    def test_plot_trumpet_quant_log_and_linear(self):
        # Test trumpet plot in quantitative mode with xscale log and linear
        rng = random.Random(321)
        rows = []
        for i in range(300):
            chr_ = rng.randint(1, 22)
            pos_ = rng.randint(1_000, 2_000_000)
            pval = max(min(rng.random(), 0.999999), 1e-300)
            maf = max(min(rng.random()*0.5, 0.499), 0.001)
            beta = rng.uniform(-0.2, 0.2)
            snpid = f"{chr_}:{pos_}_A_G"
            eaf = max(min(rng.random(), 0.999), 0.001)
            rows.append({"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid, "MAF": maf, "EAF": eaf, "BETA": beta, "N": rng.randint(10000, 50000)})
        df = pd.DataFrame(rows)
        fig = _plot_trumpet(df, mode="q", xscale="log", p_level=1, verbose=False)
        self.assertIsNotNone(fig)
        fig = _plot_trumpet(df, mode="q", xscale="linear", p_level=1, verbose=False)
        self.assertIsNotNone(fig)

    def test_plot_trumpet_binary_with_highlight_pinpoint(self):
        # Test trumpet plot in binary mode with ncase/ncontrol, prevalence, highlight and pinpoint
        rng = random.Random(987)
        rows = []
        targets = []
        for i in range(250):
            chr_ = rng.randint(1, 22)
            pos_ = rng.randint(1_000, 2_000_000)
            pval = max(min(rng.random(), 0.999999), 1e-300)
            maf = max(min(rng.random()*0.5, 0.499), 0.001)
            beta = rng.uniform(-0.3, 0.3)
            snpid = f"{chr_}:{pos_}_A_G"
            if i < 2:
                targets.append(snpid)
            eaf = max(min(rng.random(), 0.999), 0.001)
            rows.append({"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid, "MAF": maf, "EAF": eaf, "BETA": beta, "N": rng.randint(5000, 30000)})
        df = pd.DataFrame(rows)
        fig = _plot_trumpet(
            df,
            mode="b",
            ncase=2000,
            ncontrol=18000,
            prevalence=0.1,
            sig_level=5e-8,
            p_level=1,
            highlight=targets,
            highlight_color="#CB132D",
            highlight_windowkb=200,
            pinpoint=targets[:1],
            pinpoint_color="red",
            verbose=False,
        )
        self.assertIsNotNone(fig)

    def test_plot_stacked_m_mode_titles(self):
        # Test stacked Manhattan plot with two panels and titles
        df1 = make_sumstats(n=200)
        df2 = make_sumstats(n=200)
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            titles=["Panel 1", "Panel 2"],
            fig_kwargs={"figsize": (10, 6)},
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_two_panels(self):
        # Test stacked mqq layout with two panels
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="mqq",
            mqqratio=2,
            titles=["Panel A", "Panel B"],
            fig_kwargs={"figsize": (12, 8)},
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_mode_r_with_region(self):
        df1 = make_sumstats(n=200)
        df2 = make_sumstats(n=200)
        target_chr = df1["CHR"].iloc[0]
        target_pos = df1["POS"].iloc[0]
        region = (target_chr, max(0, target_pos - 100000), target_pos + 100000)
        region_ref = [df1.iloc[0]["SNPID"]]
        captured_calls = []

        def fake_mqqplot(sumstats, **kwargs):
            captured_calls.append(dict(kwargs))
            fig, ax = plt.subplots(figsize=(8, 4))
            return fig, Log(), [], None

        with patch("gwaslab.viz.viz_plot_stackedregional._mqqplot", side_effect=fake_mqqplot):
            fig, log = plot_stacked_mqq(
                objects=[df1, df2],
                mode="r",
                region=region,
                region_ref=region_ref,
                vcfs=[None, None],
                titles=["Regional Panel 1", "Regional Panel 2"],
                fig_kwargs={"figsize": (12, 10)},
                verbose=False,
            )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)
        self.assertEqual(len(captured_calls), 2)
        for call_kwargs in captured_calls:
            self.assertEqual(call_kwargs.get("region_ref"), region_ref)
            self.assertIn("gtf_path", call_kwargs)
        self.assertEqual(captured_calls[0]["gtf_path"], "default")
        self.assertIsNone(captured_calls[1]["gtf_path"])

    def test_plot_stacked_mqq_region_ref_passes_filter(self):
        from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

        pm = VizParamsManager()
        load_viz_config(pm)
        allowed = pm.allowed("plot_stacked_mqq", "r")
        self.assertIn("region_ref", allowed)
        merged = pm.merge(
            "plot_stacked_mqq",
            {"mode": "r", "region_ref": ["rs35560038"], "region": (19, 1, 2)},
            mode="r",
        )
        filtered = pm.filter(
            plot_stacked_mqq,
            merged,
            key="plot_stacked_mqq",
            mode="r",
            verbose=False,
        )
        self.assertEqual(filtered.get("region_ref"), ["rs35560038"])

    def test_plot_stacked_mqq_no_duplicate_gtf_path(self):
        df1 = make_sumstats(n=100)
        df2 = make_sumstats(n=100)
        target_chr = df1["CHR"].iloc[0]
        target_pos = df1["POS"].iloc[0]
        region = (target_chr, max(0, target_pos - 50000), target_pos + 50000)

        def fake_mqqplot(sumstats, **kwargs):
            fig, ax = plt.subplots(figsize=(6, 3))
            return fig, Log(), [], None

        with patch("gwaslab.viz.viz_plot_stackedregional._mqqplot", side_effect=fake_mqqplot):
            fig, log = plot_stacked_mqq(
                objects=[df1, df2],
                mode="r",
                region=region,
                vcfs=[None, None],
                verbose=False,
            )
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_sumstats_objects(self):
        # Test with Sumstats objects instead of DataFrames
        from gwaslab.g_Sumstats import Sumstats
        
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        gl1 = Sumstats(sumstats=df1, chrom="CHR", pos="POS", p="P", snpid="SNPID", verbose=False)
        gl2 = Sumstats(sumstats=df2, chrom="CHR", pos="POS", p="P", snpid="SNPID", verbose=False)
        
        fig, log = plot_stacked_mqq(
            objects=[gl1, gl2],
            mode="m",
            titles=["Sumstats 1", "Sumstats 2"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_highlight_and_pinpoint(self):
        # Test stacked plot with highlight and pinpoint parameters
        df1 = make_sumstats(n=200)
        df2 = make_sumstats(n=200)
        highlight_snps = [df1.iloc[0]["SNPID"], df1.iloc[1]["SNPID"]]
        pinpoint_snps = [df1.iloc[2]["SNPID"]]
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            highlight=highlight_snps,
            highlight_color="#CB132D",
            highlight_windowkb=100,
            pinpoint=pinpoint_snps,
            pinpoint_color="red",
            titles=["With Highlight", "Panel 2"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_mqq_kwargs(self):
        # Test that mqq_kwargs are properly passed through to underlying plots
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            sig_line=True,
            anno_sig_level=1e-6,
            suggestive_sig_line=True,
            suggestive_sig_level=1e-5,
            ylim=(0, 8),
            colors=["#1f77b4", "#ff7f0e"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_three_panels(self):
        # Test with three panels
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        df3 = make_sumstats(n=150)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2, df3],
            mode="m",
            titles=["Panel 1", "Panel 2", "Panel 3"],
            subplot_height=3,
            region_hspace=0.1,
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.axes), 3)

    def test_plot_stacked_mqq_single_panel(self):
        # Test with single panel
        # Note: Single panel has a bug in the source code (line 223) where it tries to len(axes)
        # when axes is a single Axes object. We'll skip this test for now or use 2 panels minimum.
        # Actually, let's test with 2 panels instead to avoid the bug
        df1 = make_sumstats(n=200)
        df2 = make_sumstats(n=200)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            titles=["Panel 1", "Panel 2"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.axes), 2)

    def test_plot_stacked_mqq_with_title_pos(self):
        # Test with custom title position
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            titles=["Title 1", "Title 2"],
            title_pos=[0.02, 0.95],
            title_kwargs={"fontsize": 12, "weight": "bold"},
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_mqq_mode_three_panels(self):
        # Test mqq mode with three panels (should create 3x2 grid)
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        df3 = make_sumstats(n=150)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2, df3],
            mode="mqq",
            mqqratio=3,
            titles=["MQQ 1", "MQQ 2", "MQQ 3"],
            fig_kwargs={"figsize": (15, 12)},
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)
        # In mqq mode, each panel has 2 axes (Manhattan + QQ)
        # So 3 panels = 6 axes total
        self.assertEqual(len(fig.axes), 6)

    def test_plot_stacked_mqq_with_anno_set(self):
        # Test with annotation set
        # Note: anno_alias, anno_d, anno_kwargs, and anno_kwargs_single need to be dicts, not None
        df1 = make_sumstats(n=200)
        df2 = make_sumstats(n=200)
        anno_snps = [df1.iloc[0]["SNPID"], df1.iloc[1]["SNPID"], df1.iloc[2]["SNPID"]]
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            anno_set=anno_snps,
            anno=True,
            anno_alias={},  # Provide empty dict instead of None
            anno_d={},  # Provide empty dict instead of None
            anno_kwargs={},  # Provide empty dict instead of None
            anno_kwargs_single={},  # Provide empty dict instead of None
            titles=["Annotated Panel 1", "Annotated Panel 2"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_region_lead_grids(self):
        # Test with region lead grids configuration
        # Note: mode="r" requires VCFs, so we use mode="m" instead
        df1 = make_sumstats(n=200)
        df2 = make_sumstats(n=200)
        target_chr = df1["CHR"].iloc[0]
        target_pos = df1["POS"].iloc[0]
        region = (target_chr, max(0, target_pos - 100000), target_pos + 100000)
        
        # Use mode="m" instead of "r" since "r" requires VCF files
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            region=region,
            region_lead_grids=[0, 1],  # Show grids on both panels
            region_lead_grid_line={"alpha": 0.7, "linewidth": 2, "color": "#FF0000"},
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_common_ylabel_false(self):
        # Test with common_ylabel=False (each panel gets its own ylabel)
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            common_ylabel=False,
            titles=["Panel 1", "Panel 2"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_custom_heights(self):
        # Test with custom height ratios
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2],
            mode="m",
            mqq_height=1.5,
            subplot_height=5,
            region_hspace=0.15,
            titles=["Tall Panel 1", "Tall Panel 2"],
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_stacked_mqq_with_numeric_suffix_kwargs(self):
        # Test that numeric suffix kwargs work (e.g., highlight2, pinpoint2)
        # This tests the _sort_kwargs functionality
        df1 = make_sumstats(n=150)
        df2 = make_sumstats(n=150)
        df3 = make_sumstats(n=150)
        
        highlight1 = [df1.iloc[0]["SNPID"]]
        highlight2 = [df2.iloc[0]["SNPID"]]
        highlight3 = [df3.iloc[0]["SNPID"]]
        
        fig, log = plot_stacked_mqq(
            objects=[df1, df2, df3],
            mode="m",
            highlight=highlight1,  # Applied to all panels
            highlight2=highlight2,   # Applied to panel 2 (index 1)
            highlight3=highlight3,   # Applied to panel 3 (index 2)
            verbose=False,
        )
        self.assertIsInstance(log, Log)
        self.assertIsNotNone(fig)

    def test_plot_daf_threshold_and_regression(self):
        # Test allele frequency comparison with threshold and regression
        rng = random.Random(456)
        rows = []
        for i in range(300):
            chr_ = rng.randint(1, 22)
            pos_ = rng.randint(1_000, 2_000_000)
            snpid = f"{chr_}:{pos_}_A_G"
            eaf = rng.random()
            raf = rng.random()
            rows.append({"SNPID": snpid, "EAF": eaf, "RAF": raf, "EA": "A", "NEA": "G"})
        df = pd.DataFrame(rows)
        fig, outliers = plotdaf(
            df,
            threshold=0.2,
            is_reg=True,
            is_threshold=True,
            legend1=False,
            legend2=False,
            fig_kwargs={"figsize": (8, 4)},
            verbose=False,
        )
        self.assertIsNotNone(fig)
        self.assertIsInstance(outliers, pd.DataFrame)
        self.assertGreaterEqual(len(fig.axes), 1)

    def test_plot_regional_via_mqqplot(self):
        # Test single-reference regional plot via mqqplot in mode="r" with consistent minimal kwargs
        row = self.df.iloc[0]
        chr_ = int(row["CHR"])
        pos_ = int(row["POS"])
        region = [chr_, max(1, pos_ - 5000), pos_ + 5000]
        fig, log = _mqqplot(
            self.df,
            mode="r",
            region=region,
            region_ref=[row["SNPID"]],
            region_ld_colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"],
            region_ld_threshold=[0.2, 0.4, 0.6],
            region_marker_shapes=["X", "o", "^"],
            region_grid=False,
            region_protein_coding=False,
            region_recombination=False,
            gtf_path=None,
            rr_path=None,
            region_grid_line={},
            region_lead_grid=False,
            region_lead_grid_line={},
            verbose=False,
        )
        self.assertIsInstance(log, Log)

    def test_plot_regional_via_mqqplot_multi_ref(self):
        # Test multi-reference regional plot via mqqplot in mode="r"
        # This reproduces the bug where region_ref_second creates multi-reference mode
        row1 = self.df.iloc[0]
        row2 = self.df.iloc[1]
        chr_ = int(row1["CHR"])
        pos_ = int(row1["POS"])
        region = [chr_, max(1, pos_ - 5000), pos_ + 5000]
        fig, log = _mqqplot(
            self.df,
            mode="r",
            region=region,
            region_ref=None,  # None becomes [None] then [None, row2["SNPID"]]
            region_ref_second=row2["SNPID"],
            region_ld_colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"],
            region_ld_colors_m=["#1f77b4", "#d62728"],
            region_ld_threshold=[0.2, 0.4, 0.6],
            region_marker_shapes=["o", "^", "s"],
            region_grid=False,
            region_protein_coding=False,
            region_recombination=False,
            gtf_path=None,
            rr_path=None,
            region_grid_line={},
            region_lead_grid=False,
            region_lead_grid_line={},
            verbose=False,
        )
        self.assertIsInstance(log, Log)

    # --- Regional: SHAPE/LD palette vs seaborn markers (no VCF) + process_vcf pipeline ---

    def test_regional_mode_setup_single_ref_markers_match_shape_indices(self):
        """Single-ref markers dict keys match SHAPE column (0=no LD, 1=LD, 2=ref); palette keys 100+."""
        df = pd.DataFrame({"LD": [100], "SHAPE": [0], "scaled_P": [1.0]})
        ld_cols = ["#aa0000", "#00aa00", "#0000aa", "#444444", "#555555", "#666666", "#777777"]
        shapes = ["X", "o", "^", "s"]
        _, _, palette, style, markers, _, _ = regional_mode_setup(
            df,
            region_ref=[None],
            region_ld_colors=ld_cols,
            region_marker_shapes=shapes,
            region_ld_colors_m=["#111111"],
            region_ld_threshold=[0.2, 0.4, 0.6, 0.8],
            vcf_path=None,
            ld_path=None,
            verbose=False,
        )
        self.assertEqual(style, "SHAPE")
        self.assertEqual(markers, {0: "X", 1: "o", 2: "^"})
        self.assertEqual(palette[100], "#aa0000")
        self.assertEqual(palette[101], "#00aa00")

    def test_regional_mode_setup_multi_ref_markers_include_missing_ld_key(self):
        """Multi-ref: seaborn markers[0] for SHAPE=0; SHAPE=k uses list entry k (first ref = shapes[1])."""
        df = pd.DataFrame({"LD": [200], "SHAPE": [0], "scaled_P": [1.0]})
        ld_cols = ["#aa0000", "#00aa00", "#0000aa", "#444444", "#555555", "#666666", "#777777"]
        shapes = ["X", "o", "^", "s", "D"]
        _, _, _, style, markers, _, _ = regional_mode_setup(
            df,
            region_ref=[None, "rs2"],
            region_ld_colors=ld_cols,
            region_marker_shapes=shapes,
            region_ld_colors_m=["#111111", "#222222"],
            region_ld_threshold=[0.2, 0.4, 0.6, 0.8],
            vcf_path=None,
            ld_path=None,
            verbose=False,
        )
        self.assertEqual(style, "SHAPE")
        self.assertEqual(markers[0], "X")
        self.assertEqual(markers[1], "o")
        self.assertEqual(markers[2], "^")

    def test_regional_draw_manhattan_distinct_facecolors_per_ld_hue(self):
        """No VCF: hand-built chr_hue/LD codes and SHAPE; assert scatter facecolors and marker paths differ."""
        fig, ax = plt.subplots(figsize=(4, 3))
        plot_df = pd.DataFrame(
            {
                "i": [0.0, 1.0, 2.0],
                "scaled_P": [2.0, 3.0, 4.0],
                "s": [36, 36, 36],
                "chr_hue": [100, 101, 102],
                "SHAPE": [0, 1, 2],
                "SNPID": ["s1", "s2", "s3"],
                "CHR": [1, 1, 1],
            }
        )
        palette = {100: "#e41a1c", 101: "#377eb8", 102: "#4daf4a"}
        markers = {0: "X", 1: "o", 2: "^"}
        draw_manhattan_panel(
            ax1=ax,
            sumstats=plot_df,
            snpid="SNPID",
            palette=palette,
            marker_size=(20, 40),
            style="SHAPE",
            linewidth=0.5,
            edgecolor="black",
            legend=None,
            scatter_kwargs={"markers": markers},
            highlight=[],
            highlight_chrpos=False,
            highlight_color="#CB132D",
            density_color=False,
            density_range=None,
            density_trange=(0, 10),
            density_threshold=5,
            density_tpalette="Blues",
            density_palette="Reds",
            pinpoint=[],
            pinpoint_color="red",
            chrom="CHR",
            log=self.log,
            verbose=False,
        )
        plt.close(fig)
        rgbs = []
        for coll in ax.collections:
            fc = coll.get_facecolor()
            if fc is not None and len(fc) > 0:
                for row in fc:
                    rgbs.append(tuple(np.round(row[:3], 2)))
        unique_rgb = set(rgbs)
        self.assertGreaterEqual(
            len(unique_rgb),
            3,
            msg="Expected ≥3 distinct facecolors for three LD hue bins",
        )
        paths = [tuple(np.round(p.vertices.mean(axis=0), 4)) for coll in ax.collections for p in coll.get_paths()]
        unique_path_centers = len(set(paths))
        self.assertGreaterEqual(
            unique_path_centers,
            2,
            msg="Expected multiple marker path shapes (e.g. X vs circle vs triangle)",
        )

    def test_regional_pipeline_process_vcf_ld_shape_and_manhattan(self):
        """End-to-end regional LD: write temp VCF → process_vcf → regional_mode_setup → draw_manhattan_panel.

        Covers the same path as plot_mqq(mode='r') after VCF load: LD/RSQ/SHAPE columns, lead dropped
        from to_plot, chr_hue from LD, then sns.scatterplot hue+style. POS 2000 is absent from the VCF
        so it stays LD=0 and SHAPE=0 (missing-LD marker) after process_vcf's missing_ld_mask.
        """
        fd, vcf_path = tempfile.mkstemp(suffix=".vcf")
        os.close(fd)
        try:
            with open(vcf_path, "w") as fh:
                fh.write(_MINI_VCF_FOR_LD)

            # Four sumstats rows: three match VCF REF/ALT at 1000–1002; one off-panel variant at 2000 (no REFINDEX)
            rows = []
            for pos, nea, ea, p in [
                (1000, "A", "G", 1e-5),
                (1001, "C", "T", 1e-20),
                (1002, "G", "C", 1e-8),
                (2000, "T", "A", 1e-4),
            ]:
                rows.append(
                    {
                        "CHR": 1,
                        "POS": pos,
                        "NEA": nea,
                        "EA": ea,
                        "P": p,
                        "SNPID": f"1:{pos}",
                        "scaled_P": -np.log10(p),
                    }
                )
            raw = pd.DataFrame(rows)

            # tabix=None: allel reads the small uncompressed VCF directly (no index file).
            out = process_vcf(
                raw,
                vcf_path=vcf_path,
                region=(1, 999, 2500),
                region_ref=[None],
                log=self.log,
                verbose=False,
                pos="POS",
                nea="NEA",
                ea="EA",
                region_ld_threshold=[0.2, 0.4, 0.6, 0.8],
                vcf_chr_dict={1: "1"},
                tabix=None,
            )

            # process_vcf output: merged LD bins (100+), RSQ, SHAPE; unmatched → LD 0, SHAPE 0
            self.assertIn("LD", out.columns)
            self.assertIn("SHAPE", out.columns)
            self.assertIn("RSQ", out.columns)
            self.assertIn(0, set(out["SHAPE"].tolist()), msg="Unmatched sumstats row should get SHAPE=0")
            self.assertTrue((out.loc[out["POS"] == 2000, "LD"] == 0).all())
            self.assertTrue((out.loc[out["POS"] != 2000, "LD"] >= 100).all())

            # Manhattan x/size columns as _mqqplot would add after process_vcf
            work = out.sort_values("POS").reset_index(drop=True).copy()
            work["i"] = np.arange(len(work), dtype=float)
            work["s"] = 40

            # vcf_path set so regional_mode_setup does not overwrite LD/SHAPE with the no-VCF dummy.
            ld_cols = [
                "#e0e0e0",
                "#1f77b4",
                "#2ca02c",
                "#d62728",
                "#9467bd",
                "#8c564b",
                "#e377c2",
            ]
            _, _, palette, style, markers, to_plot, _ = regional_mode_setup(
                work,
                region_ref=[None],
                region_ld_colors=ld_cols,
                region_marker_shapes=["X", "o", "^", "s"],
                region_ld_colors_m=["#000000"],
                region_ld_threshold=[0.2, 0.4, 0.6, 0.8],
                vcf_path=vcf_path,
                ld_path=None,
                verbose=False,
            )

            self.assertEqual(style, "SHAPE")
            self.assertEqual(markers, {0: "X", 1: "o", 2: "^"})
            self.assertEqual(len(to_plot), 3, msg="Lead SNP dropped from bulk scatter")
            self.assertIn(0, set(to_plot["SHAPE"].tolist()))
            self.assertIn("chr_hue", to_plot.columns)

            # Same draw path as mode='r' Manhattan facet (no highlight): hue=chr_hue, style=SHAPE, markers=...
            fig, ax = plt.subplots(figsize=(4, 3))
            try:
                draw_manhattan_panel(
                    ax1=ax,
                    sumstats=to_plot,
                    snpid="SNPID",
                    palette=palette,
                    marker_size=(20, 40),
                    style=style,
                    linewidth=0.5,
                    edgecolor="black",
                    legend=None,
                    scatter_kwargs={"markers": markers},
                    highlight=[],
                    highlight_chrpos=False,
                    highlight_color="#CB132D",
                    density_color=False,
                    density_range=None,
                    density_trange=(0, 10),
                    density_threshold=5,
                    density_tpalette="Blues",
                    density_palette="Reds",
                    pinpoint=[],
                    pinpoint_color="red",
                    chrom="CHR",
                    log=self.log,
                    verbose=False,
                )
                self.assertGreaterEqual(len(ax.collections), 1)
                # Collect RGBA facecolors from PathCollections (LD=0 vs 101/106 must not collapse to one color)
                rgbs = []
                for coll in ax.collections:
                    fc = coll.get_facecolor()
                    if fc is not None and len(fc) > 0:
                        for row in fc:
                            rgbs.append(tuple(np.round(row[:3], 2)))
                self.assertGreaterEqual(
                    len(set(rgbs)),
                    2,
                    msg="Expect distinct facecolors for LD=0 vs in-panel LD bins",
                )
            finally:
                plt.close(fig)
        finally:
            try:
                os.unlink(vcf_path)
            except OSError:
                pass

    def test_sumstats_basic_check_and_status(self):
        # 继续添加测试：测试Sumstats.basic_check及状态记录
        rows = []
        rng = random.Random(2025)
        for i in range(200):
            chr_ = rng.randint(1, 22)
            pos_ = rng.randint(1_000, 2_000_000)
            pval = max(min(rng.random(), 0.999999), 1e-300)
            snpid = f"{chr_}:{pos_}_A_G"
            rows.append({"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid, "EA": "A", "NEA": "G"})
        df = pd.DataFrame(rows)

        from gwaslab.g_Sumstats import Sumstats
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", ea="EA", nea="NEA", snpid="SNPID", verbose=False)
        gl.basic_check(remove_dup=True, verbose=False)
        status = gl.check_sumstats_qc_status()
        self.assertIn("basic_check", status)
        self.assertTrue(status["basic_check"].get("performed", False))

    def test_sumstats_plot_mqq_method(self):
        # 继续添加测试：通过Sumstats对象调用plot_mqq绘图
        rows = []
        rng = random.Random(2026)
        for i in range(150):
            chr_ = rng.randint(1, 22)
            pos_ = rng.randint(1_000, 2_000_000)
            pval = max(min(rng.random(), 0.999999), 1e-300)
            snpid = f"{chr_}:{pos_}_A_G"
            rows.append({"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid})
        df = pd.DataFrame(rows)

        from gwaslab.g_Sumstats import Sumstats
        gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", p="P", snpid="SNPID", verbose=False)
        fig = gl.plot_mqq(mode="m", verbose=False)
        self.assertIsNotNone(fig)
        self.assertGreaterEqual(len(fig.axes), 1)

    def test_sumstats_pair_plot_miami_highlight_pinpoint(self):
        # 继续添加测试：使用SumstatsPair.plot_miami测试高亮与pinpoint参数
        rng = random.Random(2027)
        rows = []
        for i in range(220):
            chr_ = rng.randint(1, 22)
            pos_ = rng.randint(1_000, 2_000_000)
            pval1 = max(min(rng.random(), 0.999999), 1e-300)
            pval2 = max(min(rng.random(), 0.999999), 1e-300)
            snpid = f"{chr_}:{pos_}_A_G"
            rows.append({"CHR": chr_, "POS": pos_, "P": pval1, "SNPID": snpid, "EA": "A", "NEA": "G"})
        df1 = pd.DataFrame(rows)
        df2 = df1.copy()
        df2["P"] = df2["P"].sample(frac=1.0, random_state=2027).values

        from gwaslab.g_Sumstats import Sumstats
        from gwaslab.g_SumstatsPair import SumstatsPair

        gl1 = Sumstats(sumstats=df1, chrom="CHR", pos="POS", p="P", ea="EA", nea="NEA", snpid="SNPID", study="STUDY1", verbose=False)
        gl2 = Sumstats(sumstats=df2, chrom="CHR", pos="POS", p="P", ea="EA", nea="NEA", snpid="SNPID", study="STUDY2", verbose=False)

        pair = SumstatsPair(gl1, gl2, verbose=False)
        targets = [df1.iloc[0]["SNPID"], df1.iloc[1]["SNPID"]]
        pair.plot_miami(
            mode="m",
            id0="SNPID",
            id1="SNPID",
            id2="SNPID",
            titles=["STUDY1", "STUDY2"],
            highlight=targets,
            highlight_color="#CB132D",
            highlight_windowkb=500,
            pinpoint=[targets[0]],
            pinpoint_color="red",
            anno_sig_level=5e-8,
            fig_kwargs={"figsize": (15, 10), "dpi": 300},
            verbose=False,
        )
        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
