import os
import sys
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.g_Sumstats import Sumstats
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_build import _process_build
from gwaslab.viz.viz_plot_regional2 import process_gtf

_EMPTY_GTF = pd.DataFrame(
    columns=["seqname", "start", "end", "strand", "feature", "gene_biotype", "gene_id", "gene_name"]
)
_REGION = (2, 54513738, 56513738)


def _empty_gtf(*args, **kwargs):
    gtf = _EMPTY_GTF.copy()
    if kwargs.get("if_return_path"):
        return gtf, None
    return gtf


def _region_sumstats():
    df = pd.DataFrame(
        {
            "CHR": [2, 2, 2],
            "POS": [55000000, 55500000, 56000000],
            "P": [1e-5, 1e-6, 1e-4],
            "SNPID": ["2:55000000_A_G", "2:55500000_A_G", "2:56000000_A_G"],
        }
    )
    ss = Sumstats(df, fmt="gwaslab", verbose=False)
    ss._build = "99"
    return ss


class TestFillBuild(unittest.TestCase):
    def setUp(self):
        self.ss = _region_sumstats()

    def test_process_build_maps_37_to_19(self):
        self.assertEqual(_process_build("37", verbose=False), "19")

    def test_fill_build_prefers_call_time_over_object(self):
        params = self.ss._fill_build({}, build="37")
        self.assertEqual(params["build"], "37")
        self.assertEqual(self.ss._build, "99")

    def test_fill_build_falls_back_to_object_when_unset(self):
        params = self.ss._fill_build({}, build=None)
        self.assertEqual(params["build"], "99")

    def test_fill_build_keeps_existing_params_when_call_time_is_none(self):
        params = self.ss._fill_build({"build": "38"}, build=None)
        self.assertEqual(params["build"], "38")


class TestPlotMqqHonorsBuild(unittest.TestCase):
    def setUp(self):
        self.ss = _region_sumstats()
        self.plot_kwargs = dict(
            region=_REGION,
            gtf_path="default",
            rr_path=None,
            region_recombination=False,
            region_grid=False,
            region_lead_grid=False,
            verbose=False,
        )

    def _assert_gtf_build(self, mock_gtf, expected):
        mock_gtf.assert_called()
        self.assertEqual(mock_gtf.call_args.kwargs.get("build"), expected)

    def test_plot_mqq_build_37_reaches_get_gtf_as_19(self):
        with patch("gwaslab.viz.viz_plot_regional2.get_gtf", side_effect=_empty_gtf) as mock_gtf:
            self.ss.plot_mqq(mode="r", build="37", **self.plot_kwargs)
        self._assert_gtf_build(mock_gtf, "19")
        self.assertEqual(self.ss._build, "99")

    def test_plot_mqq_build_38_reaches_get_gtf(self):
        with patch("gwaslab.viz.viz_plot_regional2.get_gtf", side_effect=_empty_gtf) as mock_gtf:
            self.ss.plot_mqq(mode="r", build="38", **self.plot_kwargs)
        self._assert_gtf_build(mock_gtf, "38")

    def test_plot_region_build_37_reaches_get_gtf_as_19(self):
        with patch("gwaslab.viz.viz_plot_regional2.get_gtf", side_effect=_empty_gtf) as mock_gtf:
            self.ss.plot_region(build="37", **self.plot_kwargs)
        self._assert_gtf_build(mock_gtf, "19")

    def test_plot_mqq_without_build_keeps_unknown(self):
        with patch("gwaslab.viz.viz_plot_regional2.get_gtf", side_effect=_empty_gtf) as mock_gtf:
            self.ss.plot_mqq(mode="r", **self.plot_kwargs)
        self._assert_gtf_build(mock_gtf, "99")


class TestProcessGtfUnknownBuild(unittest.TestCase):
    def test_unknown_build_warns_and_returns_empty(self):
        log = Log()
        with patch("gwaslab.viz.viz_plot_regional2.get_gtf", side_effect=_empty_gtf):
            genes, exons = process_gtf(
                gtf_path="default",
                region=_REGION,
                region_flank_factor=0.2,
                build="99",
                region_protein_coding=True,
                gtf_chr_dict=get_number_to_chr(),
                gtf_gene_name=None,
                verbose=False,
                log=log,
            )
        self.assertIn("build=99", log.log_text)
        self.assertIn("infer_build", log.log_text)
        self.assertEqual(len(genes), 0)
        self.assertEqual(len(exons), 0)

    def test_known_build_does_not_warn(self):
        log = Log()
        with patch("gwaslab.viz.viz_plot_regional2.get_gtf", side_effect=_empty_gtf):
            process_gtf(
                gtf_path="default",
                region=_REGION,
                region_flank_factor=0.2,
                build="19",
                region_protein_coding=True,
                gtf_chr_dict=get_number_to_chr(),
                gtf_gene_name=None,
                verbose=False,
                log=log,
            )
        self.assertNotIn("build=99", log.log_text)


class TestGetLeadHonorsBuild(unittest.TestCase):
    def test_get_lead_passes_call_time_build(self):
        ss = _region_sumstats()
        with patch("gwaslab.g_Sumstats._get_sig", return_value=ss.data) as mock_sig:
            ss.get_lead(build="37")
        self.assertEqual(mock_sig.call_args.kwargs.get("build"), "37")

    def test_get_top_passes_call_time_build(self):
        ss = _region_sumstats()
        with patch("gwaslab.g_Sumstats._get_top", return_value=ss.data) as mock_top:
            ss.get_top(build="38")
        self.assertEqual(mock_top.call_args.kwargs.get("build"), "38")


if __name__ == "__main__":
    unittest.main()
