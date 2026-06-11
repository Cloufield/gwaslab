import os
import sys
import tempfile
import unittest

import matplotlib

matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab import Sumstats, plot_lead_overlap
from gwaslab.info.g_Log import Log


def make_sumstats(rows, study_name):
    df = pd.DataFrame(rows)
    ss = Sumstats(
        sumstats=df,
        snpid="SNPID",
        chrom="CHR",
        pos="POS",
        p="P",
        verbose=False,
    )
    ss.meta["gwaslab"]["study_name"] = study_name
    return ss


class TestLeadOverlapPlot(unittest.TestCase):
    def test_two_study_venn_overlap_counts(self):
        ss1 = make_sumstats(
            [
                {"SNPID": "s1_shared", "CHR": 1, "POS": 1_000_000, "P": 1e-10, "GENE": "GENEA"},
                {"SNPID": "s1_only", "CHR": 1, "POS": 10_000_000, "P": 1e-9, "GENE": "GENEB"},
                {"SNPID": "s1_noise", "CHR": 2, "POS": 50_000, "P": 0.5, "GENE": "NOISE"},
            ],
            "Study1",
        )
        ss2 = make_sumstats(
            [
                {"SNPID": "s2_shared", "CHR": 1, "POS": 1_100_000, "P": 1e-12, "GENE": "GENEA"},
                {"SNPID": "s2_only", "CHR": 2, "POS": 2_000_000, "P": 1e-9, "GENE": "GENEC"},
            ],
            "Study2",
        )

        overlap_df, fig, log = plot_lead_overlap(
            objects=[ss1, ss2],
            titles=["A", "B"],
            anno=False,
            windowsizekb_for_overlap=500,
            verbose=False,
        )

        counts = overlap_df["MEMBERSHIP_KEY"].value_counts().to_dict()
        self.assertEqual(counts["1|1"], 1)
        self.assertEqual(counts["1|0"], 1)
        self.assertEqual(counts["0|1"], 1)
        self.assertIsNotNone(fig)
        self.assertIsInstance(log, Log)

    def test_three_study_venn_overlap_counts(self):
        objects = [
            make_sumstats([{"SNPID": "s1", "CHR": 1, "POS": 1_000_000, "P": 1e-10}], "S1"),
            make_sumstats([{"SNPID": "s2", "CHR": 1, "POS": 1_050_000, "P": 1e-10}], "S2"),
            make_sumstats([{"SNPID": "s3", "CHR": 1, "POS": 1_100_000, "P": 1e-10}], "S3"),
        ]

        overlap_df, fig, _ = plot_lead_overlap(
            objects=objects,
            titles=["S1", "S2", "S3"],
            anno=False,
            windowsizekb_for_overlap=500,
            verbose=False,
        )

        self.assertEqual(len(overlap_df), 1)
        self.assertEqual(overlap_df.iloc[0]["MEMBERSHIP_KEY"], "1|1|1")
        self.assertIsNotNone(fig)

    def test_auto_uses_upset_for_four_studies(self):
        objects = [
            make_sumstats([{"SNPID": f"s{i}", "CHR": 1, "POS": 1_000_000 + i * 10_000, "P": 1e-10}], f"S{i}")
            for i in range(4)
        ]

        overlap_df, fig, _ = plot_lead_overlap(
            objects=objects,
            titles=["S0", "S1", "S2", "S3"],
            anno=False,
            mode="auto",
            windowsizekb_for_overlap=500,
            verbose=False,
        )

        self.assertEqual(len(overlap_df), 1)
        self.assertEqual(overlap_df.iloc[0]["N_STUDIES"], 4)
        self.assertEqual(overlap_df.iloc[0]["SET_ID"], "Set1")
        set_list = overlap_df.attrs["set_list"]
        self.assertEqual(set_list.iloc[0]["SET_ID"], "Set1")
        self.assertEqual(set_list.iloc[0]["MEMBERSHIP_KEY"], "1|1|1|1")
        self.assertEqual(list(fig.axes[1].get_xticks()), [])
        connector_lines = fig.axes[1].get_lines()
        self.assertEqual(len(connector_lines), 1)
        self.assertEqual(connector_lines[0].get_solid_capstyle(), "butt")
        self.assertEqual(connector_lines[0].get_color(), "#597FBD")
        self.assertGreater(connector_lines[0].get_ydata()[0], 0)
        self.assertLess(connector_lines[0].get_ydata()[1], 3)
        self.assertIsNotNone(fig)

    def test_save_to_explicit_path(self):
        ss1 = make_sumstats([{"SNPID": "s1", "CHR": 1, "POS": 1_000_000, "P": 1e-10}], "S1")
        ss2 = make_sumstats([{"SNPID": "s2", "CHR": 2, "POS": 2_000_000, "P": 1e-10}], "S2")

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "lead_overlap.png")
            plot_lead_overlap(
                objects=[ss1, ss2],
                titles=["S1", "S2"],
                anno=False,
                save=path,
                verbose=False,
            )
            self.assertTrue(os.path.exists(path))

    def test_bad_inputs(self):
        ss1 = make_sumstats([{"SNPID": "s1", "CHR": 1, "POS": 1_000_000, "P": 1e-10}], "S1")
        with self.assertRaises(ValueError):
            plot_lead_overlap(objects=[ss1], anno=False, verbose=False)
        with self.assertRaises(ValueError):
            plot_lead_overlap(objects=[ss1, ss1, ss1, ss1], mode="venn", anno=False, verbose=False)


if __name__ == "__main__":
    unittest.main()
