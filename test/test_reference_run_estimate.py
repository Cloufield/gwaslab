"""Tests for reference run estimation and progress utilities."""
import json
import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ["GWASLAB_CALIBRATION_ROOT"] = str(
    ROOT / "benchmarks" / "reference_run_calibration" / "runs"
)

from gwaslab.util.util_in_reference_run import (
    ComputeProfile,
    CpuTier,
    format_bytes,
    format_duration,
    estimate_run_plan,
    load_calibration_runs,
    parse_harmonize_log_steps,
    plan_harmonize_reference_steps,
    RunProgressTracker,
    StorageProfile,
    _hybrid_estimate,
    _parse_duration_token,
)
from gwaslab.info.g_Log import Log


class TestReferenceRunFormat(unittest.TestCase):
    def test_format_duration(self):
        self.assertEqual(format_duration(45), "45s")
        self.assertIn("m", format_duration(125))

    def test_format_bytes(self):
        self.assertIn("MB", format_bytes(1024 ** 2))


class TestReferenceRunEstimate(unittest.TestCase):
    def test_hdd_slower_than_ssd(self):
        compute = ComputeProfile(logical_cores=8, tier=CpuTier.BASELINE)
        ssd = estimate_run_plan(
            "bcf_lookup_extract",
            storage_profile=StorageProfile.SSD,
            compute_profile=compute,
            threads=6,
            task_count=22,
            target_variants=100000,
            ref_info=type("I", (), {"size_bytes": 10 ** 9, "storage_profile": StorageProfile.SSD})(),
        )
        hdd = estimate_run_plan(
            "bcf_lookup_extract",
            storage_profile=StorageProfile.HDD,
            compute_profile=compute,
            threads=6,
            task_count=22,
            target_variants=100000,
            ref_info=type("I", (), {"size_bytes": 10 ** 9, "storage_profile": StorageProfile.HDD})(),
        )
        self.assertGreater(hdd.total_estimate_low, ssd.total_estimate_low)

    def test_more_threads_faster_parallel(self):
        compute = ComputeProfile(logical_cores=12, tier=CpuTier.BASELINE)
        t1 = estimate_run_plan(
            "tabix_vcf_query",
            storage_profile=StorageProfile.SSD,
            compute_profile=compute,
            threads=1,
            task_count=4,
            target_variants=50000,
        )
        t6 = estimate_run_plan(
            "tabix_vcf_query",
            storage_profile=StorageProfile.SSD,
            compute_profile=compute,
            threads=6,
            task_count=4,
            target_variants=50000,
        )
        self.assertLess(t6.total_estimate_high, t1.total_estimate_high)

    def test_effective_workers_cap(self):
        compute = ComputeProfile(logical_cores=4, tier=CpuTier.BASELINE)
        self.assertEqual(compute.effective_workers(16, 22), 4)

    def test_parse_duration_token(self):
        self.assertEqual(_parse_duration_token("9m 16s"), 556.0)
        self.assertEqual(_parse_duration_token("60.197s"), 60.197)

    def test_ci_calibration_not_used_for_full_scale_bcf_lookup(self):
        compute = ComputeProfile(logical_cores=12, tier=CpuTier.BASELINE)
        ref_info = type(
            "I",
            (),
            {
                "size_bytes": 12 * 1024 ** 3,
                "storage_profile": StorageProfile.HDD,
                "indexed": True,
            },
        )()
        plan = estimate_run_plan(
            "bcf_lookup_extract",
            storage_profile=StorageProfile.HDD,
            compute_profile=compute,
            threads=4,
            task_count=23,
            target_variants=12_557_761,
            ref_info=ref_info,
        )
        self.assertLess(plan.total_estimate_high, 3600.0)

    def test_full_scale_bcf_lookup_heuristic_order_of_minutes(self):
        compute = ComputeProfile(logical_cores=12, tier=CpuTier.BASELINE)
        ref_info = type(
            "I",
            (),
            {
                "size_bytes": 2_955_303_556,
                "storage_profile": StorageProfile.HDD,
                "indexed": True,
            },
        )()
        plan = estimate_run_plan(
            "bcf_lookup_extract",
            storage_profile=StorageProfile.HDD,
            compute_profile=compute,
            threads=4,
            task_count=23,
            target_variants=12_557_761,
            ref_info=ref_info,
        )
        self.assertLess(plan.total_estimate_high, 1800.0)
        self.assertGreater(plan.total_estimate_low, 120.0)

    def test_fasta_check_estimate_minutes_not_hours(self):
        compute = ComputeProfile(logical_cores=12, tier=CpuTier.BASELINE)
        ref_info = type(
            "I",
            (),
            {
                "size_bytes": 905 * 1024 ** 2,
                "storage_profile": StorageProfile.HDD,
                "indexed": True,
            },
        )()
        plan = estimate_run_plan(
            "fasta_load_check",
            storage_profile=StorageProfile.HDD,
            compute_profile=compute,
            chroms_loaded=23,
            variants_to_check=12_557_761,
            ref_info=ref_info,
        )
        self.assertLess(plan.total_estimate_high, 600.0)

    def test_harmonize_plan_total_under_one_hour(self):
        plan = plan_harmonize_reference_steps(
            ref_seq="/home/yunye/.gwaslab/hg19.fa.gz",
            ref_infer="/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz",
            ref_rsid_vcf="/home/yunye/CommonData/Reference/dbsnp/ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.25.gz.strip.bcf",
            sweep_mode=True,
            threads=4,
            extract_threads=4,
            sumstats_rows=12_557_761,
            target_infer=12_557_761,
            target_rsid=12_557_761,
            variants_check_ref=12_557_761,
            chroms_in_sumstats=23,
            storage_profile="hdd",
        )
        self.assertLess(plan.total_estimate_high, 3600.0)


class TestHarmonizeLogParse(unittest.TestCase):
    def test_parse_full_bbj_log(self):
        log_path = ROOT / "benchmarks/reference_run_calibration/runs/full/full_bbj_harmonize_sweep.log"
        if not log_path.is_file():
            self.skipTest("full benchmark log missing")
        steps = parse_harmonize_log_steps(log_path.read_text())
        names = [s["name"] for s in steps]
        self.assertIn("check_ref", names)
        self.assertIn("infer_extract", names)
        self.assertIn("rsid_extract", names)
        infer = next(s for s in steps if s["name"] == "infer_extract")
        self.assertGreater(infer["wall_seconds"], 400.0)
        self.assertLess(infer["wall_seconds"], 800.0)


class TestRunProgressTracker(unittest.TestCase):
    def test_eta_after_updates(self):
        tracker = RunProgressTracker(total=10, label="chromosomes")
        tracker.update("chr1")
        tracker.update("chr2")
        eta = tracker.eta_seconds()
        self.assertIsNotNone(eta)
        self.assertGreater(eta, 0)

    def test_log_progress(self):
        log = Log()
        tracker = RunProgressTracker(total=5, label="workers")
        tracker.update("w1")
        tracker.log_progress(log, verbose=False)
        self.assertIn("Progress", log.log_text)


class TestCalibrationJson(unittest.TestCase):
    def test_load_ci_runs_if_present(self):
        runs = load_calibration_runs(force=True)
        ci_dir = ROOT / "benchmarks" / "reference_run_calibration" / "runs" / "ci"
        if ci_dir.is_dir():
            for p in ci_dir.glob("*.json"):
                with open(p) as f:
                    data = json.load(f)
                self.assertIn("operation", data)
                self.assertIn("timing", data)

    def test_bundled_ci_in_package(self):
        from gwaslab.util.util_in_reference_run import _bundled_calibration_runs

        bundled = _bundled_calibration_runs()
        self.assertIsNotNone(bundled)
        self.assertTrue(any(bundled.glob("ci/*.json")))

    def test_load_merges_without_env(self):
        old = os.environ.pop("GWASLAB_CALIBRATION_ROOT", None)
        try:
            runs = load_calibration_runs(force=True)
            self.assertGreaterEqual(len(runs), 5)
        finally:
            if old is not None:
                os.environ["GWASLAB_CALIBRATION_ROOT"] = old
            load_calibration_runs(force=True)


if __name__ == "__main__":
    unittest.main()
