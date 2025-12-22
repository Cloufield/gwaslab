import os
import sys
import unittest
import random

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
import numpy as np

import gwaslab
from gwaslab.info.g_Log import Log


def make_sumstats(n=500, with_eaf=True, with_density=False):
    # Create a synthetic sumstats DataFrame with CHR/POS/P/SNPID and optional EAF/DENSITY
    rng = random.Random(42)
    rows = []
    for i in range(n):
        chr_ = rng.randint(1, 22)
        pos_ = rng.randint(1_000, 2_000_000)
        pval = max(min(rng.random(), 0.999999), 1e-300)
        snpid = f"{chr_}:{pos_}_A_G"
        row = {"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid}
        if with_eaf:
            row["EAF"] = rng.random()
        if with_density:
            row["DENSITY"] = rng.random() * 10
        rows.append(row)
    return pd.DataFrame(rows)


class TestMQQPlotOptions(unittest.TestCase):
    def setUp(self):
        self.df = make_sumstats(n=300, with_eaf=True, with_density=True)
        self.sumstats = gwaslab.Sumstats(sumstats=self.df, snpid="SNPID", chrom="CHR", pos="POS", p="P", eaf="EAF", build="19", verbose=False)
        self.log = Log()

    def test_m_mode_basic(self):
        # Test mode="m" (Manhattan) basic plotting
        fig = self.sumstats.plot_manhattan(mode="m", verbose=False)
        self.assertEqual(len(fig.axes), 1)

    def test_mqq_mode_has_two_axes(self):
        # Test mode="mqq" combined Manhattan+QQ layout (two axes)
        fig = self.sumstats.plot_mqq(mode="mqq", verbose=False)
        self.assertEqual(len(fig.axes), 2)

    def test_qq_mode_basic(self):
        # Test mode="qq" basic QQ plotting
        fig = self.sumstats.plot_qq(mode="qq", verbose=False)
        self.assertEqual(len(fig.axes), 1)

    def test_qq_xlim_and_labels(self):
        # Test QQ plot axis limits and custom x tick labels (numeric)
        fig = self.sumstats.plot_qq(
            mode="qq",
            qq_xlim=(0, 4),
            qq_xlabels=[0, 1, 2, 3, 4],
            verbose=False,
        )
        ax = fig.axes[0]
        x0, x1 = ax.get_xlim()
        self.assertAlmostEqual(x0, 0, places=2)
        self.assertAlmostEqual(x1, 4, places=2)

    def test_b_mode_basic(self):
        # Test mode="b" (density/Brisbane) basic plotting
        fig = self.sumstats.plot_snp_density(mode="b", verbose=False)
        self.assertGreaterEqual(len(fig.axes), 1)

    def test_highlight_and_pinpoint(self):
        # Test Manhattan highlighting by SNPID and pinpoint single variant
        target = self.df.iloc[0]["SNPID"]
        fig = self.sumstats.plot_manhattan(
            mode="m",
            highlight=[target],
            highlight_windowkb=100,
            pinpoint=[target],
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_stratified_qq(self):
        # Test stratified QQ by MAF bins with custom colors
        fig = self.sumstats.plot_qq(
            mode="qq",
            stratified=True,
            maf_bins=[[0, 0.01], [0.01, 0.05], [0.05, 0.25], [0.25, 0.5]],
            maf_bin_colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"],
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_ylim_option_applied(self):
        # Test y-axis limits via ylim
        fig = self.sumstats.plot_manhattan(mode="m", ylim=(0, 5), verbose=False)
        y0, y1 = fig.axes[0].get_ylim()
        self.assertAlmostEqual(y0, 0, places=3)
        self.assertAlmostEqual(y1, 5, places=3)

    def test_sig_line_and_additional_lines(self):
        # Test significance line, suggestive line, and additional horizontal lines
        fig = self.sumstats.plot_manhattan(
            mode="m",
            sig_line=True,
            sig_level_plot=1e-6,
            suggestive_sig_line=True,
            suggestive_sig_level=1e-5,
            additional_line=[1e-4, 1e-3],
            additional_line_color=["#888", "#777"],
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_scaled_mlog10p(self):
        # Test scaled=True using precomputed MLOG10P values
        df = self.df.copy()
        df["MLOG10P"] = -np.log10(df["P"].clip(lower=1e-300))
        sumstats = gwaslab.Sumstats(sumstats=df, snpid="SNPID", chrom="CHR", pos="POS", p="P", eaf="EAF", build="19", verbose=False)
        fig = sumstats.plot_manhattan(mode="m", scaled=True, verbose=False)
        self.assertEqual(len(fig.axes), 1)

    def test_x_padding_options(self):
        # Test x-axis padding options: xpad, xpadl, xpadr
        fig = self.sumstats.plot_manhattan(mode="m", xpad=0.02, verbose=False)
        self.assertEqual(len(fig.axes), 1)
        fig = self.sumstats.plot_manhattan(mode="m", xpadl=0.02, xpadr=None, verbose=False)
        self.assertEqual(len(fig.axes), 1)
        fig = self.sumstats.plot_manhattan(mode="m", xpad=None, xpadr=0.02, verbose=False)
        self.assertEqual(len(fig.axes), 1)

    def test_highlight_chrpos_mode(self):
        # Test highlighting by CHR:POS using highlight_chrpos=True
        # Use CHR:POS highlighting
        row = self.df.iloc[1]
        highlight_chrpos = [[int(row["CHR"]), int(row["POS"] )]]
        fig = self.sumstats.plot_manhattan(
            mode="m",
            highlight=highlight_chrpos,
            highlight_chrpos=True,
            highlight_windowkb=50,
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_pinpoint_multiple(self):
        # Test pinpointing multiple variants by SNPID
        targets = self.df["SNPID"].head(2).tolist()
        fig = self.sumstats.plot_manhattan(
            mode="m",
            pinpoint=targets,
            pinpoint_color="red",
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_density_b_mode(self):
        # Test density options in mode="b": density_color/range/threshold
        fig = self.sumstats.plot_snp_density(
            mode="b",
            density_color=True,
            density_range=(0, 10),
            density_threshold=3,
            verbose=False,
        )
        self.assertGreaterEqual(len(fig.axes), 1)

    def test_save_png(self):
        # Test saving PNG via save path
        fig = self.sumstats.plot_manhattan(mode="m", save="./_tmp_m.png", verbose=False)
        self.assertTrue(os.path.exists("./_tmp_m.png"))
        os.remove("./_tmp_m.png")

    def test_save_pdf(self):
        # Test saving PDF via save path
        fig = self.sumstats.plot_manhattan(mode="m", save="./_tmp_m.pdf", verbose=False)
        self.assertTrue(os.path.exists("./_tmp_m.pdf"))
        os.remove("./_tmp_m.pdf")

    def test_anno_true_chrpos(self):
        # Test annotation with anno=True (annotate by chr:pos)
        fig = self.sumstats.plot_manhattan(
            mode="m",
            anno=True,
            anno_set=[self.df.iloc[0]["SNPID"], self.df.iloc[1]["SNPID"]],
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_anno_custom_column(self):
        # Test annotation using a custom column name via anno="Annotation"
        df = self.df.copy()
        df["Annotation"] = [f"ANN{i}" for i in range(len(df))]
        sumstats = gwaslab.Sumstats(sumstats=df, snpid="SNPID", chrom="CHR", pos="POS", p="P", eaf="EAF", build="19", verbose=False)
        fig = sumstats.plot_manhattan(
            mode="m",
            anno="Annotation",
            anno_set=[df.iloc[0]["SNPID"], df.iloc[1]["SNPID"], df.iloc[2]["SNPID"]],
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_anno_style_variants(self):
        # Test different annotation styles: right/tight/expand
        target = self.df.iloc[0]["SNPID"]
        for style in ["right", "tight", "expand"]:
            fig = self.sumstats.plot_manhattan(
                mode="m",
                anno=True,
                anno_set=[target],
                anno_style=style,
                verbose=False,
            )
            self.assertEqual(len(fig.axes), 1)

    def test_anno_kwargs_and_single(self):
        # Test anno_kwargs for global styling and anno_kwargs_single for per-variant overrides
        t0 = self.df.iloc[0]["SNPID"]
        t1 = self.df.iloc[1]["SNPID"]
        fig = self.sumstats.plot_manhattan(
            mode="m",
            anno=True,
            anno_set=[t0, t1],
            anno_kwargs={"color": "blue", "fontsize": 8},
            anno_kwargs_single={t1: {"color": "green"}},
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_anno_arm_controls(self):
        # Test arm controls: anno_fixed_arm_length, arm_scale, anno_height, anno_d direction map
        targets = [self.df.iloc[i]["SNPID"] for i in range(3)]
        fig = self.sumstats.plot_manhattan(
            mode="m",
            anno=True,
            anno_set=targets,
            anno_fixed_arm_length=0.02,
            arm_scale=1.2,
            anno_height=1.1,
            anno_d={"0": "l", "1": "r"},
            repel_force=0.04,
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_titles_for_mqq(self):
        # Test titles for mqq layout: mtitle and qtitle
        fig = self.sumstats.plot_mqq(
            mode="mqq",
            mtitle="Manhattan Title",
            qtitle="QQ Title",
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 2)

    def test_cut_and_skip(self):
        # Test cut and skip parameters; just ensure no exceptions
        fig = self.sumstats.plot_manhattan(
            mode="m",
            cut=2,
            skip=1,
            cut_log=False,
            verbose=False,
        )
        self.assertEqual(len(fig.axes), 1)

    def test_vector_save_svg(self):
        # Test saving as SVG (vector), which adjusts dpi and rasterization settings
        fig = self.sumstats.plot_manhattan(mode="m", save="./_tmp_m.svg", verbose=False)
        self.assertTrue(os.path.exists("./_tmp_m.svg"))
        os.remove("./_tmp_m.svg")

    def test_anno_kwargs_none_fix(self):
        # Test that anno=True works even if anno_kwargs is implicitly None (bug fix verification)
        # This was causing AttributeError: 'NoneType' object has no attribute 'items'
        fig = self.sumstats.plot_mqq(
            mode="m",
            anno=True,
            anno_kwargs=None, # Explicitly testing None
            verbose=False
        )
        self.assertEqual(len(fig.axes), 1)

    def test_regional_plot_basic(self):
        # Test regional plot mode="r"
        # We need to define a region (chr, start, end)
        # Using a region that likely contains some data from make_sumstats
        # make_sumstats(n=300) creates random data. Let's find a region with points.
        target_chr = self.df["CHR"].iloc[0]
        target_pos = self.df["POS"].iloc[0]
        region = (target_chr, max(0, target_pos - 100000), target_pos + 100000)
        
        fig = self.sumstats.plot_mqq(
            mode="r",
            region=region,
            verbose=False,
        )
        # Check if fig is a Figure
        self.assertIsInstance(fig, matplotlib.figure.Figure)

    def test_regional_plot_with_highlight(self):
        # Test regional plot with highlight
        target_chr = self.df["CHR"].iloc[0]
        target_pos = self.df["POS"].iloc[0]
        region = (target_chr, max(0, target_pos - 100000), target_pos + 100000)
        highlight_snp = self.df["SNPID"].iloc[0]
        
        fig = self.sumstats.plot_mqq(
            mode="r",
            region=region,
            highlight=[highlight_snp],
            verbose=False
        )
        self.assertIsInstance(fig, matplotlib.figure.Figure)

    def test_mqq_with_anno(self):
         # Test mqq mode with annotation
        fig = self.sumstats.plot_mqq(
            mode="mqq",
            anno=True,
            anno_set=[self.df["SNPID"].iloc[0]],
            verbose=False
        )
        self.assertEqual(len(fig.axes), 2)


if __name__ == "__main__":
    unittest.main()
