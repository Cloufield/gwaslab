"""overwrite + TSV/VCF cascade for assign_rsid.

overwrite applies once to existing rsIDs; TSV then VCF fill remaining.
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm.hm_harmonize_sumstats import (
    _parallelize_assign_rsid,
    _prepare_rsid_overwrite,
)


STATUS_OK = 9980000
STATUS_BAD = 9981000  # digit 4 != 0


def _fake_assign_rsid_single(
    sumstats,
    path,
    rsid="rsID",
    chr="CHR",
    pos="POS",
    ref="NEA",
    alt="EA",
    mapper=None,
):
    return pd.Series(
        [f"rs_vcf_{int(c)}_{int(p)}" for c, p in zip(sumstats[chr], sumstats[pos])],
        index=sumstats.index,
        dtype="string",
    )


def _fixture_df():
    # A: user rsID, in TSV and VCF
    # B: user rsID, VCF only
    # C: empty, in TSV
    # D: empty, VCF only
    # E: valid rsID, bad STATUS, SNPID in TSV
    return pd.DataFrame(
        {
            "SNPID": ["1:100:A:G", "1:200:A:G", "1:300:A:G", "1:400:A:G", "1:500:A:G"],
            "CHR": [1, 1, 1, 1, 1],
            "POS": [100, 200, 300, 400, 500],
            "NEA": ["A", "A", "A", "A", "A"],
            "EA": ["G", "G", "G", "G", "G"],
            "STATUS": [STATUS_OK, STATUS_OK, STATUS_OK, STATUS_OK, STATUS_BAD],
            "rsID": pd.Series(
                ["rs_user_A", "rs_user_B", pd.NA, pd.NA, "rs_keep_E"],
                dtype="string",
            ),
        }
    )


def _write_tsv(path):
    pd.DataFrame(
        {
            "SNPID": ["1:100:A:G", "1:300:A:G", "1:500:A:G"],
            "rsID": ["rs_tsv_A", "rs_tsv_C", "rs_tsv_E"],
        }
    ).to_csv(path, sep="\t", index=False)


class TestPrepareRsidOverwrite(unittest.TestCase):
    def test_all_clears_status_ok_only(self):
        df = _fixture_df()
        out = _prepare_rsid_overwrite(df.copy(), overwrite="all", verbose=False)
        self.assertTrue(pd.isna(out.loc[0, "rsID"]))
        self.assertTrue(pd.isna(out.loc[1, "rsID"]))
        self.assertEqual(out.loc[4, "rsID"], "rs_keep_E")

    def test_empty_does_not_clear(self):
        df = _fixture_df()
        out = _prepare_rsid_overwrite(df.copy(), overwrite="empty", verbose=False)
        self.assertEqual(out.loc[0, "rsID"], "rs_user_A")
        self.assertTrue(pd.isna(out.loc[2, "rsID"]))

    def test_invalid_string_raises(self):
        with self.assertRaises(ValueError):
            _prepare_rsid_overwrite(_fixture_df(), overwrite="nope", verbose=False)


class TestTsvAssignOverwrite(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
        self.tmp.close()
        _write_tsv(self.tmp.name)

    def tearDown(self):
        try:
            os.remove(self.tmp.name)
        except OSError:
            pass

    def test_empty_does_not_clobber_existing(self):
        out = _parallelize_assign_rsid(
            _fixture_df(),
            path=self.tmp.name,
            ref_mode="tsv",
            overwrite="empty",
            verbose=False,
            log_run_plan=False,
        )
        self.assertEqual(out.loc[0, "rsID"], "rs_user_A")
        self.assertEqual(out.loc[2, "rsID"], "rs_tsv_C")
        self.assertEqual(out.loc[4, "rsID"], "rs_keep_E")

    def test_all_overwrites_status_ok_only(self):
        out = _parallelize_assign_rsid(
            _fixture_df(),
            path=self.tmp.name,
            ref_mode="tsv",
            overwrite="all",
            verbose=False,
            log_run_plan=False,
        )
        self.assertEqual(out.loc[0, "rsID"], "rs_tsv_A")
        self.assertEqual(out.loc[2, "rsID"], "rs_tsv_C")
        self.assertEqual(out.loc[4, "rsID"], "rs_keep_E")


class TestDualRefAssignRsid(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
        self.tmp.close()
        _write_tsv(self.tmp.name)

    def tearDown(self):
        try:
            os.remove(self.tmp.name)
        except OSError:
            pass

    def _sumstats(self):
        raw = _fixture_df()
        ss = Sumstats(raw, fmt="gwaslab", verbose=False)
        ss.data["STATUS"] = raw["STATUS"].to_numpy()
        return ss

    def _run(self, overwrite):
        ss = self._sumstats()
        patches = [
            patch(
                "gwaslab.hm.hm_harmonize_sumstats.assign_rsid_single",
                _fake_assign_rsid_single,
            ),
            patch(
                "gwaslab.bd.bd_chromosome_mapper.ChromosomeMapper.detect_reference_format",
                return_value=("numeric", ""),
            ),
        ]
        with patches[0], patches[1]:
            ss.assign_rsid(
                ref_rsid_tsv=self.tmp.name,
                ref_rsid_vcf="dummy.vcf.gz",
                overwrite=overwrite,
                verbose=False,
                log_run_plan=False,
            )
        return ss.data

    def test_empty_cascade(self):
        out = self._run("empty")
        self.assertEqual(out.loc[0, "rsID"], "rs_user_A")
        self.assertEqual(out.loc[1, "rsID"], "rs_user_B")
        self.assertEqual(out.loc[2, "rsID"], "rs_tsv_C")
        self.assertEqual(out.loc[3, "rsID"], "rs_vcf_1_400")
        self.assertEqual(out.loc[4, "rsID"], "rs_keep_E")

    def test_all_clears_then_tsv_wins_overlap(self):
        out = self._run("all")
        self.assertEqual(out.loc[0, "rsID"], "rs_tsv_A")
        self.assertEqual(out.loc[1, "rsID"], "rs_vcf_1_200")
        self.assertEqual(out.loc[2, "rsID"], "rs_tsv_C")
        self.assertEqual(out.loc[3, "rsID"], "rs_vcf_1_400")
        self.assertEqual(out.loc[4, "rsID"], "rs_keep_E")


class TestHarmonizeSweepOverwriteKwargs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
        self.tmp.close()
        _write_tsv(self.tmp.name)

    def tearDown(self):
        try:
            os.remove(self.tmp.name)
        except OSError:
            pass

    def test_sweep_receives_empty_after_prepare(self):
        ss = Sumstats(_fixture_df(), fmt="gwaslab", verbose=False)
        ss.data["STATUS"] = _fixture_df()["STATUS"].to_numpy()
        captured = {}

        def _fake_assign_rsid(sumstats_obj, **kwargs):
            captured["overwrite"] = kwargs.get("overwrite")
            if hasattr(sumstats_obj, "data"):
                return sumstats_obj.data
            return sumstats_obj

        with patch("gwaslab.g_Sumstats._assign_rsid", side_effect=_fake_assign_rsid):
            ss.harmonize(
                basic_check=False,
                ref_rsid_tsv=self.tmp.name,
                ref_rsid_vcf="dummy.vcf.gz",
                sweep_mode=True,
                assign_rsid_kwargs={"overwrite": "all", "log_run_plan": False},
                verbose=False,
            )
        self.assertEqual(captured.get("overwrite"), "empty")
        self.assertEqual(ss.data.loc[0, "rsID"], "rs_tsv_A")
        self.assertEqual(ss.data.loc[4, "rsID"], "rs_keep_E")


if __name__ == "__main__":
    unittest.main()
