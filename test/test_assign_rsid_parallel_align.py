"""Regression tests for issue #225: VCF assign_rsid row alignment.

In 4.2.0–4.2.1, parallel VCF workers used imap_unordered and wrote rsIDs
back with .values, which permuted rsIDs across variants.
"""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gwaslab.qc.qc_fix_sumstats import _df_split
from gwaslab.hm.hm_harmonize_sumstats import (
    _apply_assigned_rsid,
    _parallelize_assign_rsid,
)


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
    """Encode source CHR:POS in the rsID so permutation is detectable."""
    return pd.Series(
        [f"rs{int(c)}_{int(p)}" for c, p in zip(sumstats[chr], sumstats[pos])],
        index=sumstats.index,
        dtype="string",
    )


class ReverseImapPool:
    """In-process Pool stand-in that yields chunks last-first (imap_unordered)."""

    def __init__(self, _n):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def imap(self, func, iterable):
        results = [func(item) for item in iterable]
        return reversed(results)


def _issue225_sumstats(n=12000):
    return pd.DataFrame(
        {
            "CHR": (np.arange(n) % 22) + 1,
            "POS": np.arange(n) + 1000,
            "NEA": "A",
            "EA": "G",
            "STATUS": 9980000,
            "rsID": pd.Series(pd.NA, index=range(n), dtype="string"),
        }
    )


def _rsid_chr(series):
    return series.str.extract(r"^rs(\d+)_", expand=False).astype("Int64")


def _values_write_like_421(sumstats, assigned_rsid, rsid="rsID"):
    to_assign = sumstats[rsid].isna()
    sumstats.loc[to_assign, rsid] = assigned_rsid.values
    return sumstats


class TestApplyAssignedRsidIndexAlign(unittest.TestCase):
    def _sumstats_and_unordered_workers(self):
        sumstats = pd.DataFrame(
            {
                "CHR": [8, 8, 13, 2],
                "POS": [100, 200, 300, 400],
                "rsID": pd.Series([pd.NA, pd.NA, pd.NA, pd.NA], dtype="string"),
            }
        )
        later_chunk = pd.Series(
            ["rs11148602", "rs3131956"],
            index=[2, 3],
            name="rsID",
        )
        earlier_chunk = pd.Series(
            ["rsChr8a", "rsChr8b"],
            index=[0, 1],
            name="rsID",
        )
        unordered = pd.concat([later_chunk, earlier_chunk])
        expected = pd.Series(
            ["rsChr8a", "rsChr8b", "rs11148602", "rs3131956"],
            index=sumstats.index,
        )
        return sumstats, unordered, expected

    def test_index_align_survives_completion_order_concat(self):
        sumstats, unordered, expected = self._sumstats_and_unordered_workers()
        out = _apply_assigned_rsid(sumstats.copy(), unordered, rsid="rsID")
        pd.testing.assert_series_equal(
            out["rsID"].astype("string"),
            expected.astype("string"),
            check_names=False,
        )

    def test_positional_values_assignment_is_the_225_bug(self):
        sumstats, unordered, expected = self._sumstats_and_unordered_workers()
        buggy = sumstats.copy()
        to_assign = buggy["rsID"].isna()
        buggy.loc[to_assign, "rsID"] = unordered.values
        self.assertTrue((buggy["rsID"] != expected.astype("string")).any())
        self.assertEqual(buggy.loc[0, "rsID"], "rs11148602")
        self.assertEqual(buggy.loc[0, "CHR"], 8)

    def test_empty_series_is_a_no_op(self):
        sumstats = pd.DataFrame({"rsID": pd.Series([pd.NA], dtype="string")})
        out = _apply_assigned_rsid(sumstats.copy(), pd.Series(dtype="string"), rsid="rsID")
        self.assertTrue(out["rsID"].isna().all())


class TestParallelAssignRsidIssue225(unittest.TestCase):
    """Drive `_parallelize_assign_rsid` with reversed worker completion.

    4.2.0–4.2.1 used imap_unordered + ``.values``. That path must scramble
    CHR vs rsID; the index-aligned write must not.
    """

    def _run(self, apply_fn=None):
        df = _issue225_sumstats(12000)
        patches = [
            patch("gwaslab.hm.hm_harmonize_sumstats.Pool", ReverseImapPool),
            patch(
                "gwaslab.hm.hm_harmonize_sumstats.assign_rsid_single",
                _fake_assign_rsid_single,
            ),
            patch(
                "gwaslab.bd.bd_chromosome_mapper.ChromosomeMapper.detect_reference_format",
                return_value=("numeric", ""),
            ),
        ]
        if apply_fn is not None:
            patches.append(
                patch(
                    "gwaslab.hm.hm_harmonize_sumstats._apply_assigned_rsid",
                    apply_fn,
                )
            )
        with patches[0], patches[1], patches[2]:
            if apply_fn is not None:
                with patches[3]:
                    out = _parallelize_assign_rsid(
                        df,
                        path="dummy.vcf.gz",
                        ref_mode="vcf",
                        threads=4,
                        verbose=False,
                        log_run_plan=False,
                    )
            else:
                out = _parallelize_assign_rsid(
                    df,
                    path="dummy.vcf.gz",
                    ref_mode="vcf",
                    threads=4,
                    verbose=False,
                    log_run_plan=False,
                )
        mismatched = int((_rsid_chr(out["rsID"]) != out["CHR"]).sum())
        return out, mismatched

    def test_index_align_survives_reversed_workers(self):
        out, mismatched = self._run()
        self.assertEqual(mismatched, 0)
        self.assertEqual(
            out.loc[0, "rsID"],
            f"rs{int(out.loc[0, 'CHR'])}_{int(out.loc[0, 'POS'])}",
        )

    def test_values_write_with_reversed_workers_is_issue_225(self):
        out, mismatched = self._run(apply_fn=_values_write_like_421)
        self.assertGreater(mismatched, 0)
        got_chr = int(_rsid_chr(out["rsID"]).iloc[0])
        row_chr = int(out.loc[0, "CHR"])
        self.assertNotEqual(got_chr, row_chr)


def _reversed_worker_concat(df, n_chunks=4):
    """Same split as production, concat last chunk first (imap_unordered)."""
    return pd.concat(_df_split(df, n_chunks)[::-1])


class TestLatentValuesWriteIfUnordered(unittest.TestCase):
    """Document the same .values hazard still used by infer_strand / normalize.

    Production still uses ordered pool.map, so these paths are OK today.
    If concat order is reversed, loc[mask] = series.values scrambles rows.
    """

    def test_infer_strand_status_values_scrambles_when_reversed(self):
        n = 8
        sumstats = pd.DataFrame(
            {
                "CHR": [8, 8, 13, 13, 2, 2, 1, 1],
                "POS": np.arange(n) + 100,
                "STATUS": np.arange(n, dtype=np.int64),
            }
        )
        mask = pd.Series([True] * n)
        inferred = pd.Series(
            10_000 + np.arange(n),
            index=sumstats.index,
            name="STATUS",
        )
        unordered = _reversed_worker_concat(inferred.to_frame("STATUS"))["STATUS"]
        buggy = sumstats.copy()
        buggy.loc[mask, "STATUS"] = unordered.values
        self.assertFalse((buggy["STATUS"].to_numpy() == inferred.to_numpy()).all())
        self.assertEqual(int(buggy.loc[0, "STATUS"]), int(inferred.iloc[n - n // 4]))

        safe = sumstats.copy()
        safe.loc[unordered.index, "STATUS"] = unordered
        self.assertTrue((safe["STATUS"].to_numpy() == inferred.to_numpy()).all())

    def test_infer_strand_indel_status_values_scrambles_when_reversed(self):
        n = 8
        sumstats = pd.DataFrame({"STATUS": np.zeros(n, dtype=np.int64)})
        inferred = pd.Series(np.arange(n) + 300, index=sumstats.index)
        unordered = _reversed_worker_concat(inferred.to_frame("STATUS"))["STATUS"]
        buggy = sumstats.copy()
        buggy.loc[sumstats.index, "STATUS"] = unordered.values
        self.assertNotEqual(int(buggy.loc[0, "STATUS"]), int(inferred.iloc[0]))
        self.assertEqual(int(buggy.loc[0, "STATUS"]), int(inferred.iloc[6]))

    def test_normalize_allele_values_scrambles_when_reversed(self):
        n = 8
        cols = ["POS", "NEA", "EA", "STATUS"]
        sumstats = pd.DataFrame(
            {
                "POS": np.arange(n) + 1000,
                "NEA": list("AAAAAAAA"),
                "EA": list("CCCCCCCC"),
                "STATUS": np.arange(n),
            }
        )
        normalized = sumstats[cols].copy()
        normalized["POS"] = normalized["POS"] + 50
        normalized["NEA"] = list("TTTTTTTT")
        unordered = _reversed_worker_concat(normalized)
        buggy = sumstats.copy()
        buggy.loc[sumstats.index, cols] = unordered.values
        self.assertNotEqual(int(buggy.loc[0, "POS"]), int(normalized.loc[0, "POS"]))
        self.assertEqual(int(buggy.loc[0, "POS"]), int(normalized.loc[6, "POS"]))
        self.assertEqual(buggy.loc[0, "NEA"], normalized.loc[6, "NEA"])

        safe = sumstats.copy()
        safe.loc[unordered.index, cols] = unordered
        pd.testing.assert_frame_equal(
            safe[cols].reset_index(drop=True),
            normalized[cols].reset_index(drop=True),
        )


if __name__ == "__main__":
    unittest.main()
