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
import matplotlib.pyplot as plt

from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2
from gwaslab.viz.viz_plot_regional2 import _plot_regional
from gwaslab.g_Log import Log
from gwaslab.viz.viz_plot_mqqplot import mqqplot
from gwaslab.viz.viz_plot_trumpetplot import plottrumpet
from gwaslab.viz.viz_plot_stackedregional import plot_stacked_mqq
from gwaslab.viz.viz_plot_compare_af import plotdaf


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
            sig_level_plot=1e-6,
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
        fig = plottrumpet(df, mode="q", xscale="log", p_level=1, verbose=False)
        self.assertIsNotNone(fig)
        fig = plottrumpet(df, mode="q", xscale="linear", p_level=1, verbose=False)
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
        fig = plottrumpet(
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
        # Test regional plot via mqqplot in mode="r" with consistent minimal kwargs
        row = self.df.iloc[0]
        chr_ = int(row["CHR"])
        pos_ = int(row["POS"])
        region = [chr_, max(1, pos_ - 5000), pos_ + 5000]
        fig, log = mqqplot(
            self.df,
            mode="r",
            region=region,
            region_ref=[row["SNPID"], None],
            region_ld_colors=["#1f77b4", "#ff7f0e", "#2ca02c"],
            region_ld_colors_m=["#1f77b4", "#ff7f0e"],
            region_ld_threshold=[0.2, 0.4, 0.6],
            region_marker_shapes=["o", "s"],
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
            sig_level_plot=5e-8,
            fig_kwargs={"figsize": (15, 10), "dpi": 300},
            verbose=False,
        )
        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
