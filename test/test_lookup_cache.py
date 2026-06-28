"""Tests for persistent lookup cache (hm_lookup_cache)."""
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "test/ref/1kg_eas_hg19.chr7_126253550_128253550.vcf.gz"


class TestLookupCachePaths(unittest.TestCase):
    def test_bundle_path_under_gwaslab(self):
        from gwaslab.bd.bd_download import get_default_directory
        from gwaslab.hm.hm_lookup_cache import lookup_cache_bundle_path

        gd = Path(get_default_directory())
        bundle = lookup_cache_bundle_path(str(gd / "EAS.test.vcf.gz"), ("AF",))
        self.assertEqual(bundle.parent, gd / "lookup")
        self.assertIn("lookup.AF", bundle.name)


class TestTargetDelta(unittest.TestCase):
    def test_delta_subset_empty(self):
        from gwaslab.hm.hm_lookup_cache import (
            compute_target_delta_by_chr,
            init_bundle_from_extract_build,
            is_lookup_bundle,
        )

        with tempfile.TemporaryDirectory() as td:
            bundle = Path(td) / "bundle"
            build = Path(td) / "build.tsv"
            build.write_text(
                "CHR\tPOS\tREF\tALT\trsID\n"
                "7\t126253551\tG\tT\trs1\n"
                "7\t126253567\tC\tT\trs2\n"
            )
            init_bundle_from_extract_build(
                str(build), bundle, str(REF), ("rsID",)
            )
            self.assertTrue(is_lookup_bundle(bundle))
            targets = pd.DataFrame({"CHR": [7, 7], "POS": [126253551, 126253567]})
            delta = compute_target_delta_by_chr(targets, bundle)
            self.assertTrue(delta.empty)

    def test_delta_new_positions(self):
        from gwaslab.hm.hm_lookup_cache import (
            compute_target_delta_by_chr,
            init_bundle_from_extract_build,
        )

        with tempfile.TemporaryDirectory() as td:
            bundle = Path(td) / "bundle"
            build = Path(td) / "build.tsv"
            build.write_text(
                "CHR\tPOS\tREF\tALT\trsID\n7\t100\tG\tT\trs1\n"
            )
            init_bundle_from_extract_build(
                str(build), bundle, str(REF), ("rsID",)
            )
            targets = pd.DataFrame({"CHR": [7, 7], "POS": [100, 200]})
            delta = compute_target_delta_by_chr(targets, bundle)
            self.assertEqual(len(delta), 1)
            self.assertEqual(int(delta.iloc[0]["POS"]), 200)


class TestResolveLookupCacheHit(unittest.TestCase):
    @unittest.skipUnless(REF.exists() and shutil.which("bcftools"), "needs ref+bcftools")
    def test_second_resolve_skips_extract(self):
        from gwaslab.bd.bd_download import get_default_directory
        from gwaslab.hm.hm_lookup_cache import (
            lookup_cache_bundle_path,
            resolve_lookup_for_vcf,
        )
        from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
        from gwaslab.info.g_Log import Log

        targets = pd.DataFrame(
            {"CHR": [7, 7, 7], "POS": [126253551, 126253567, 126253570]}
        )
        mapper = ChromosomeMapper(log=Log(), verbose=False)
        mapper.detect_reference_format(str(REF))
        log = Log()
        bundle = lookup_cache_bundle_path(str(REF), ("rsID",))
        if bundle.is_dir():
            shutil.rmtree(bundle)

        with mock.patch(
            "gwaslab.hm.hm_assign_rsid._extract_lookup_table_from_vcf_bcf"
        ) as mock_extract:
            def _fake_extract(**kwargs):
                out = kwargs["out_lookup"]
                Path(out).write_text(
                    "CHR\tPOS\tREF\tALT\trsID\n"
                    "7\t126253551\tG\tT\trs1\n"
                    "7\t126253567\tC\tT\trs2\n"
                    "7\t126253570\tG\tA\trs3\n"
                )
                return out, False

            mock_extract.side_effect = _fake_extract
            path1, _ = resolve_lookup_for_vcf(
                vcf_path=str(REF),
                targets=targets,
                assign_cols=("rsID",),
                mapper=mapper,
                reuse_lookup=True,
                extract_threads=1,
                verbose=False,
                log_run_plan=False,
                log=log,
            )
            self.assertEqual(mock_extract.call_count, 1)

            path2, _ = resolve_lookup_for_vcf(
                vcf_path=str(REF),
                targets=targets,
                assign_cols=("rsID",),
                mapper=mapper,
                reuse_lookup=True,
                extract_threads=1,
                verbose=False,
                log_run_plan=False,
                log=log,
            )
            self.assertEqual(mock_extract.call_count, 1)
            self.assertTrue(str(path2).endswith("lookup.rsID") or "lookup.rsID" in str(path2))

        # cleanup test bundle under ~/.gwaslab/lookup
        if bundle.is_dir():
            shutil.rmtree(bundle)


@unittest.skipUnless(
    REF.exists() and shutil.which("bcftools"),
    "needs ref+bcftools",
)
class TestBundleAssignParity(unittest.TestCase):
    def test_bundle_assign_matches_parquet(self):
        from gwaslab.hm.hm_assign_rsid import (
            _assign_from_lookup,
            _extract_lookup_table_from_vcf_bcf,
        )
        from gwaslab.hm.hm_lookup_cache import (
            init_bundle_from_extract_build,
            is_lookup_bundle,
        )
        from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
        from gwaslab.info.g_Log import Log

        targets = pd.DataFrame(
            {"CHR": [7, 7, 7], "POS": [126253551, 126253567, 126253570]}
        )
        sumstats = pd.DataFrame(
            {
                "CHR": [7, 7, 7],
                "POS": [126253551, 126253567, 126253570],
                "EA": ["T", "T", "A"],
                "NEA": ["G", "C", "G"],
            }
        )
        mapper = ChromosomeMapper(log=Log(), verbose=False)

        parquet_path = tempfile.NamedTemporaryFile(
            suffix=".lookup.parquet", delete=False
        ).name
        with tempfile.TemporaryDirectory() as td:
            bundle = Path(td) / "bundle"
            build_path = Path(td) / "build.tsv"
            try:
                _extract_lookup_table_from_vcf_bcf(
                    vcf_path=str(REF),
                    sumstats=targets,
                    assign_cols=["rsID"],
                    out_lookup=parquet_path,
                    extract_threads=1,
                    verbose=False,
                    log_run_plan=False,
                )
                # rebuild build tsv from parquet for bundle init test
                df = pd.read_parquet(parquet_path)
                df.to_csv(build_path, sep="\t", index=False)
                init_bundle_from_extract_build(
                    str(build_path), bundle, str(REF), ("rsID",)
                )
                self.assertTrue(is_lookup_bundle(bundle))

                out_pq = sumstats.copy()
                out_bd = sumstats.copy()
                _assign_from_lookup(
                    out_pq,
                    lookup_table=parquet_path,
                    assign_cols=("rsID",),
                    verbose=False,
                    log_run_plan=False,
                    mapper=mapper,
                )
                _assign_from_lookup(
                    out_bd,
                    lookup_table=str(bundle),
                    assign_cols=("rsID",),
                    verbose=False,
                    log_run_plan=False,
                    mapper=mapper,
                    reference_file=str(REF),
                )
                self.assertIn("rsID", out_pq.columns)
                self.assertTrue(out_pq["rsID"].equals(out_bd["rsID"]))
            finally:
                if os.path.exists(parquet_path):
                    os.remove(parquet_path)
                build_parquet = f"{parquet_path}.build"
                if os.path.exists(build_parquet):
                    os.remove(build_parquet)


if __name__ == "__main__":
    unittest.main()
