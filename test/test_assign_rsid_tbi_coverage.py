"""Tests for incomplete dbSNP TBI coverage checks in rsID assignment.

Covers chromosome-interval -r strings, lookup-stream POS gaps, and
per-chromosome assignment-rate outliers. Does not require a real NCBI TBI.
"""
import io
import os
import shlex
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

from gwaslab.hm.hm_assign_rsid import (
    _bcf_chr_region,
    _chr_assignment_rate_outliers,
    _chr_lookup_coverage_issue,
    _first_lookup_pos,
    _handle_chr_lookup_worker_result,
    _warn_chr_assignment_outliers,
    _worker_bcf_lookup,
)
from gwaslab.info.g_Log import Log


def _redirect_path_from_cmd(cmd):
    """Parse the path after ``>`` in a ``cmd > 'file'`` shell string."""
    return shlex.split(cmd.rsplit(">", 1)[-1].strip())[0]


_WORKER_FMT = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n"


class DummyMapper:
    def reference_to_sumstats(self, chr_val, reference_file=None):
        return chr_val


class TestBcfChrRegion(unittest.TestCase):
    def test_min_max_interval(self):
        self.assertEqual(_bcf_chr_region(14, [100, 500, 200]), "14:100-500")

    def test_refseq_contig(self):
        self.assertEqual(
            _bcf_chr_region("NC_000014.8", pd.Series([30, 10])),
            "NC_000014.8:10-30",
        )

    def test_single_position(self):
        self.assertEqual(_bcf_chr_region("7", [126253551]), "7:126253551-126253551")

    def test_worker_uses_interval_region(self):
        captured = {}

        def fake_check_call(cmd, shell=True):
            captured["cmd"] = cmd

        df = pd.DataFrame({"CHR": ["14", "14"], "POS": [100, 500]})
        args = (14, df, "/tmp/fake.vcf.gz", _WORKER_FMT, True, DummyMapper(), 1)
        with patch("gwaslab.hm.hm_assign_rsid.subprocess.check_call", fake_check_call):
            _worker_bcf_lookup(args)
        self.assertIn("cmd", captured)
        self.assertIn("-r 14:100-500", captured["cmd"])
        self.assertNotIn("-r 14 -T", captured["cmd"])

    def test_worker_empty_lookup_tuple(self):
        df = pd.DataFrame({"CHR": ["14", "14"], "POS": [100, 500]})
        args = (14, df, "/tmp/fake.vcf.gz", _WORKER_FMT, True, DummyMapper(), 1)
        with patch("gwaslab.hm.hm_assign_rsid.subprocess.check_call"):
            result = _worker_bcf_lookup(args)
        lookup_path, orig_chr, n_targets, min_target_pos, min_lookup_pos, n_before = result
        self.assertIsNone(lookup_path)
        self.assertEqual(str(orig_chr), "14")
        self.assertEqual(n_targets, 2)
        self.assertEqual(min_target_pos, 100)
        self.assertIsNone(min_lookup_pos)
        self.assertEqual(n_before, 2)

    def test_worker_truncated_stream_tuple(self):
        def fake_check_call(cmd, shell=True):
            path = _redirect_path_from_cmd(cmd)
            with open(path, "w") as out_f:
                out_f.write("14\t50000000\tA\tG\trs1\n")

        df = pd.DataFrame(
            {"CHR": ["14", "14", "14"], "POS": [1_000_000, 2_000_000, 50_000_000]}
        )
        args = (14, df, "/tmp/fake.vcf.gz", _WORKER_FMT, True, DummyMapper(), 1)
        lookup_path = None
        try:
            with patch(
                "gwaslab.hm.hm_assign_rsid.subprocess.check_call", fake_check_call
            ):
                result = _worker_bcf_lookup(args)
            lookup_path, orig_chr, n_targets, min_target_pos, min_lookup_pos, n_before = (
                result
            )
            self.assertEqual(str(orig_chr), "14")
            self.assertEqual(n_targets, 3)
            self.assertEqual(min_target_pos, 1_000_000)
            self.assertEqual(min_lookup_pos, 50_000_000)
            self.assertEqual(n_before, 2)
            self.assertIsNotNone(lookup_path)
            self.assertTrue(os.path.exists(lookup_path))
        finally:
            if lookup_path and os.path.exists(lookup_path):
                os.remove(lookup_path)


class TestLookupCoverageIssue(unittest.TestCase):
    def test_pos_gap_warns(self):
        # 40% of targets sit before first returned POS (50 Mb vs 1 Mb)
        msg = _chr_lookup_coverage_issue(
            chr_label="14",
            n_targets=1000,
            min_target_pos=1_000_000,
            min_lookup_pos=50_000_000,
            n_targets_before=400,
            vcf_path="/path/to/dbsnp.vcf.gz",
        )
        self.assertIsNotNone(msg)
        self.assertIn("Chromosome 14", msg)
        self.assertIn("50,000,000", msg)
        self.assertIn("bcftools index -t -f /path/to/dbsnp.vcf.gz", msg)

    def test_matching_mins_no_warning(self):
        msg = _chr_lookup_coverage_issue(
            chr_label="14",
            n_targets=1000,
            min_target_pos=1_000_000,
            min_lookup_pos=1_000_000,
            n_targets_before=0,
        )
        self.assertIsNone(msg)

    def test_small_gap_below_threshold_no_warning(self):
        msg = _chr_lookup_coverage_issue(
            chr_label="14",
            n_targets=1000,
            min_target_pos=1_000_000,
            min_lookup_pos=2_000_000,
            n_targets_before=40,
        )
        self.assertIsNone(msg)

    def test_empty_lookup_warns(self):
        msg = _chr_lookup_coverage_issue(
            chr_label="14",
            n_targets=200,
            min_target_pos=1_000_000,
            min_lookup_pos=None,
            n_targets_before=200,
            vcf_path="dbsnp.vcf.gz",
        )
        self.assertIsNotNone(msg)
        self.assertIn("no records", msg)
        self.assertIn("bcftools index -t -f dbsnp.vcf.gz", msg)

    def test_first_lookup_pos(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tmp:
            tmp.write("14\t50331648\tA\tG\trs1\n")
            tmp.write("14\t50332000\tC\tT\trs2\n")
            path = tmp.name
        try:
            self.assertEqual(_first_lookup_pos(path), 50331648)
        finally:
            Path(path).unlink(missing_ok=True)

    def test_first_lookup_pos_empty(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tmp:
            path = tmp.name
        try:
            self.assertIsNone(_first_lookup_pos(path))
        finally:
            Path(path).unlink(missing_ok=True)


class TestHandleChrLookupWorkerResult(unittest.TestCase):
    def _warning_text(self, log):
        return "".join(
            e["message"] for e in log.log_entries if e.get("type") == "warning"
        )

    def test_pos_gap_logs_warning(self):
        with tempfile.NamedTemporaryFile("w", suffix=".lookup.tsv", delete=False) as tmp:
            tmp.write("14\t50000000\tA\tG\trs1\n")
            path = tmp.name
        log = Log()
        out_f = io.StringIO()
        result = (path, "14", 1000, 1_000_000, 50_000_000, 400)
        n_rows, chr_label = _handle_chr_lookup_worker_result(
            result,
            out_f,
            DummyMapper(),
            "/tmp/dbsnp.vcf.gz",
            log,
            verbose=False,
        )
        self.assertEqual(chr_label, "14")
        self.assertEqual(n_rows, 1)
        self.assertFalse(os.path.exists(path))
        warning_text = self._warning_text(log)
        self.assertIn("Chromosome 14", warning_text)
        self.assertIn("bcftools index -t -f /tmp/dbsnp.vcf.gz", warning_text)
        self.assertIn("50,000,000", warning_text)

    def test_empty_lookup_logs_no_records(self):
        log = Log()
        out_f = io.StringIO()
        result = (None, "14", 200, 1_000_000, None, 200)
        n_rows, chr_label = _handle_chr_lookup_worker_result(
            result,
            out_f,
            DummyMapper(),
            "/tmp/dbsnp.vcf.gz",
            log,
            verbose=False,
        )
        self.assertEqual(chr_label, "14")
        self.assertEqual(n_rows, 0)
        warning_text = self._warning_text(log)
        self.assertIn("no records", warning_text)
        self.assertIn("bcftools index -t -f /tmp/dbsnp.vcf.gz", warning_text)


def _rates_from_fill(n_per_chr, fill_by_chr):
    n_eligible = pd.Series(n_per_chr)
    n_filled = pd.Series(
        {chr_val: int(round(n_per_chr[chr_val] * fill_by_chr[chr_val])) for chr_val in n_per_chr}
    )
    return n_eligible, n_filled


class TestChrAssignmentRateOutliers(unittest.TestCase):
    def test_chr14_outlier_warns(self):
        n_per_chr = {i: 1000 for i in range(1, 23)}
        fill = {i: 0.996 for i in range(1, 23)}
        fill[14] = 0.65
        n_eligible, n_filled = _rates_from_fill(n_per_chr, fill)
        outliers = _chr_assignment_rate_outliers(n_eligible, n_filled)
        chrs = [row["chr"] for row in outliers]
        self.assertEqual(chrs, [14])
        self.assertAlmostEqual(outliers[0]["rate"], 0.65, places=3)

    def test_similar_rates_no_warning(self):
        n_per_chr = {i: 1000 for i in range(1, 23)}
        fill = {i: 0.996 for i in range(1, 23)}
        n_eligible, n_filled = _rates_from_fill(n_per_chr, fill)
        self.assertEqual(_chr_assignment_rate_outliers(n_eligible, n_filled), [])

    def test_small_n_ignored(self):
        n_per_chr = {i: 499 for i in range(1, 23)}
        fill = {i: 0.996 for i in range(1, 23)}
        fill[14] = 0.10
        n_eligible, n_filled = _rates_from_fill(n_per_chr, fill)
        self.assertEqual(_chr_assignment_rate_outliers(n_eligible, n_filled), [])

    def test_warn_helper_logs(self):
        rows = []
        for chrom in range(1, 23):
            n = 600
            n_filled = 390 if chrom == 14 else 598
            chr_df = pd.DataFrame(
                {
                    "CHR": [chrom] * n,
                    "rsID": ["rs1"] * n_filled + [pd.NA] * (n - n_filled),
                }
            )
            rows.append(chr_df)
        sumstats = pd.concat(rows, ignore_index=True)
        sumstats["rsID"] = sumstats["rsID"].astype("string")
        log = Log()
        outliers = _warn_chr_assignment_outliers(
            sumstats,
            to_assign_mask=pd.Series(True, index=sumstats.index),
            log=log,
            verbose=False,
            vcf_path="/tmp/dbsnp.vcf.gz",
        )
        self.assertEqual([row["chr"] for row in outliers], [14])
        warning_text = "".join(
            e["message"] for e in log.log_entries if e.get("type") == "warning"
        )
        self.assertIn("Chromosome 14", warning_text)
        self.assertIn("bcftools index -t -f /tmp/dbsnp.vcf.gz", warning_text)


if __name__ == "__main__":
    unittest.main()
