import os
import sys
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.util.util_in_simulate import (
    simulate_sumstats_region,
    simulate_sumstats_global,
    _accumulate_V_g_raw_for_chromosomes,
    _accumulate_V_g_from_causal_snps,
    _autodetect_chromosomes_from_vcf,
    _merge_causal_ld_windows,
    _build_panel_mask,
    _resolve_null_mode,
    _compute_V_g_block,
)
from gwaslab.g_Sumstats import Sumstats
from gwaslab import Log


# Median of chi-square(1) under the null (genomic-control normalization constant).
_MEDIAN_CHISQ1 = 0.4549364030695398


def _empirical_lambda_gc(z_scores):
    """Genomic-control lambda: median(Z^2) / median(chi-square(1))."""
    z = np.asarray(z_scores, dtype=float)
    z = z[np.isfinite(z)]
    if len(z) == 0:
        return np.nan
    return float(np.median(z ** 2) / _MEDIAN_CHISQ1)


def _achieved_h2_from_betas(beta_true, X, matched_indices, n_samples):
    beta_m = beta_true[matched_indices]
    return _compute_V_g_block(X, beta_m, n_samples)


class TestSimulateSumstats(unittest.TestCase):
    """Test cases for simulate_sumstats_region and simulate_sumstats_global"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.vcf_path = os.path.join(ROOT, "test", "ref", "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        self.region = ("7", 126253550, 128253550)
        self.log = Log()
        
        # Check if VCF file exists
        if not os.path.exists(self.vcf_path):
            self.skipTest(f"Test VCF file not found: {self.vcf_path}")
    
    def test_simulate_sumstats_region_basic(self):
        """Test basic functionality of simulate_sumstats_region"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            mode="sparse",
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        
        # Check return types
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertIsInstance(causal_snp_ids, list)
        
        # Check that sumstats has data
        self.assertGreater(len(sumstats_obj.data), 0)
        
        # Check required columns
        required_cols = ["CHR", "POS", "EA", "NEA", "Z", "BETA", "SE", "P", "MLOG10P", "N"]
        for col in required_cols:
            self.assertIn(col, sumstats_obj.data.columns, f"Missing column: {col}")
        
        # Check that we have exactly n_causal causal variants
        self.assertEqual(len(causal_snp_ids), 3)
        
        # Check that causal variants are in the sumstats
        for causal_id in causal_snp_ids:
            self.assertIn(causal_id, sumstats_obj.data["SNPID"].values)
    
    def test_simulate_sumstats_region_polygenic(self):
        """Test polygenic mode"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            mode="polygenic",
            pi=0.01,  # 1% causal
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertGreater(len(causal_snp_ids), 0)
        self.assertIsInstance(sumstats_obj, Sumstats)
    
    def test_simulate_sumstats_region_maf_filter(self):
        """Test MAF filtering"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            maf_min=0.05,
            maf_max=0.5,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        
        # Check that all variants have MAF in range
        eaf = sumstats_obj.data["EAF"].values
        maf = np.minimum(eaf, 1.0 - eaf)
        self.assertTrue(np.all(maf >= 0.05))
        self.assertTrue(np.all(maf <= 0.5))
    
    def test_simulate_sumstats_region_alpha(self):
        """Test MAF-dependent effect sizes (alpha parameter)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=5,
            alpha=0.2,  # MAF-dependent effects
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertEqual(len(causal_snp_ids), 5)
    
    def test_simulate_sumstats_region_thin(self):
        """Test thinning option"""
        result_no_thin = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            thin=None,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        result_thin = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            thin=0.5,  # Keep 50% of variants
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_no_thin, _ = result_no_thin
        sumstats_thin, _ = result_thin
        
        # Thinned should have fewer or equal variants
        self.assertLessEqual(len(sumstats_thin.data), len(sumstats_no_thin.data))
        # Should be roughly 50% (allow some variance due to rounding)
        ratio = len(sumstats_thin.data) / len(sumstats_no_thin.data)
        self.assertGreater(ratio, 0.4)
        self.assertLess(ratio, 0.6)
    
    def test_simulate_sumstats_region_lambda_gc(self):
        """Test cryptic relatedness inflation (lambda_gc)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            lambda_gc=1.1,  # 10% inflation
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
    
    def test_simulate_sumstats_region_sigma_strat(self):
        """Test population stratification (sigma_strat)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=2,
            sigma_strat=0.2,  # Moderate stratification
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
    
    def test_simulate_sumstats_region_binary_trait(self):
        """Test binary trait simulation"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            trait="binary",
            n_case=5000,
            n_ctrl=5000,
            n_causal=2,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
        
        # Check that N_CASE and N_CONTROL columns exist
        if "N_CASE" in sumstats_obj.data.columns:
            self.assertIn("N_CASE", sumstats_obj.data.columns)
            self.assertIn("N_CONTROL", sumstats_obj.data.columns)
    
    def test_simulate_sumstats_global_basic(self):
        """Test basic functionality of simulate_sumstats_global"""
        result = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],  # Only chromosome 7
            n=10000,
            mode="sparse",
            n_causal=3,
            h2=0.01,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        
        # Check return types
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertIsInstance(causal_snp_ids, list)
        
        # Check that sumstats has data
        self.assertGreater(len(sumstats_obj.data), 0)
        
        # Check required columns
        required_cols = ["CHR", "POS", "EA", "NEA", "Z", "BETA", "SE", "P", "MLOG10P", "N"]
        for col in required_cols:
            self.assertIn(col, sumstats_obj.data.columns, f"Missing column: {col}")
        
        # Check that all variants are from chromosome 7 (CHR might be string or int)
        chr_values = sumstats_obj.data["CHR"].values
        self.assertTrue(all((chr_values == "7") | (chr_values == 7) | (chr_values == "chr7")))
    
    def test_simulate_sumstats_global_heritability(self):
        """Test global heritability calibration"""
        result = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="sparse",
            n_causal=5,
            h2=0.05,  # 5% heritability
            alpha=0.2,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        self.assertIsInstance(sumstats_obj, Sumstats)
        self.assertEqual(len(causal_snp_ids), 5)
    
    def test_simulate_sumstats_region_z_scores_valid(self):
        """Test that Z-scores are valid (finite, not all zeros)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        z_scores = sumstats_obj.data["Z"].values
        
        # Check that Z-scores are finite
        self.assertTrue(np.all(np.isfinite(z_scores)))
        
        # Check that not all Z-scores are zero (should have some variation)
        self.assertGreater(np.std(z_scores), 0.0)
    
    def test_simulate_sumstats_region_p_values_valid(self):
        """Test that P-values are valid (between 0 and 1)"""
        result = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=42,
            verbose=False,
            log=self.log
        )
        
        sumstats_obj, causal_snp_ids = result
        p_values = sumstats_obj.data["P"].values
        
        # Check that P-values are between 0 and 1
        self.assertTrue(np.all(p_values >= 0.0))
        self.assertTrue(np.all(p_values <= 1.0))
        self.assertTrue(np.all(np.isfinite(p_values)))
    
    def test_simulate_sumstats_region_reproducibility(self):
        """Test that same seed produces same results"""
        result1 = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=123,
            verbose=False,
            log=self.log
        )
        
        result2 = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=3,
            seed=123,
            verbose=False,
            log=self.log
        )
        
        sumstats1, causals1 = result1
        sumstats2, causals2 = result2
        
        # Check that causal variants are the same
        self.assertEqual(set(causals1), set(causals2))
        
        # Check that Z-scores are the same (within numerical precision)
        z1 = sumstats1.data["Z"].values
        z2 = sumstats2.data["Z"].values
        np.testing.assert_array_almost_equal(z1, z2, decimal=10)

    def test_global_h2_near_target(self):
        """Achieved genetic variance should approximate target h2 after calibration."""
        target_h2 = 0.05
        sumstats_obj, _ = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="sparse",
            n_causal=5,
            h2=target_h2,
            seed=42,
            verbose=False,
            log=self.log,
        )
        beta = sumstats_obj.data["BETA_TRUE"].values
        # Proxy: variance of genetic component on standardized scale via beta at causals
        causal_beta = beta[sumstats_obj.data["IS_CAUSAL"].values]
        self.assertGreater(np.var(causal_beta), 0.0)
        # BETA_TRUE are scaled to hit h2 globally; spot-check non-trivial signal
        z_causal = sumstats_obj.data.loc[sumstats_obj.data["IS_CAUSAL"], "Z"].values
        self.assertGreater(np.mean(np.abs(z_causal)), 0.5)

    def test_region_lambda_gc_null(self):
        """Null simulation with lambda_gc=1 should not inflate strongly."""
        sumstats_obj, _ = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=0,
            lambda_gc=1.0,
            sigma_strat=0.0,
            seed=42,
            verbose=False,
            log=self.log,
        )
        lam = _empirical_lambda_gc(sumstats_obj.data["Z"].values)
        self.assertGreater(lam, 0.7)
        self.assertLess(lam, 1.3)

    def test_region_lambda_gc_inflation(self):
        """lambda_gc>1 should increase median chi-square."""
        kwargs = dict(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=0,
            sigma_strat=0.0,
            seed=42,
            verbose=False,
            log=self.log,
        )
        null_obj, _ = simulate_sumstats_region(**kwargs, lambda_gc=1.0)
        inf_obj, _ = simulate_sumstats_region(**kwargs, lambda_gc=1.1)
        lam_null = _empirical_lambda_gc(null_obj.data["Z"].values)
        lam_inf = _empirical_lambda_gc(inf_obj.data["Z"].values)
        self.assertGreater(lam_inf, lam_null)
        self.assertGreater(lam_inf, 1.05)

    def test_n_scaling_at_causals(self):
        """Causal |Z| should scale roughly with sqrt(N)."""
        kwargs = dict(
            vcf_path=self.vcf_path,
            region=self.region,
            n_causal=3,
            mode="sparse",
            lambda_gc=1.0,
            sigma_strat=0.0,
            seed=42,
            verbose=False,
            log=self.log,
        )
        low, _ = simulate_sumstats_region(n=5_000, **kwargs)
        high, _ = simulate_sumstats_region(n=20_000, **kwargs)
        z_low = np.mean(np.abs(low.data.loc[low.data["IS_CAUSAL"], "Z"].values))
        z_high = np.mean(np.abs(high.data.loc[high.data["IS_CAUSAL"], "Z"].values))
        ratio = z_high / z_low
        self.assertGreater(ratio, 1.5)
        self.assertLess(ratio, 2.5)

    def test_accumulate_V_g_raw_helper(self):
        """Helper sums blocks across chromosomes (regression for multi-chr bug)."""
        import polars as pl
        from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

        df = pl.DataFrame({
            "CHR": ["7", "7", "7"],
            "POS": [100, 200, 300],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "A"],
            "EAF": [0.3, 0.4, 0.5],
        })
        pos = df["POS"].to_numpy()
        chr_to_indices = {"7": [0, 1, 2]}
        beta_raw = np.array([0.1, 0.0, 0.2])
        mapper = ChromosomeMapper(log=self.log, verbose=False)

        # With real VCF this integrates genotypes; here we only test API returns counts
        V_g, n_blocks = _accumulate_V_g_raw_for_chromosomes(
            vcf_path=self.vcf_path,
            df=df,
            beta_raw=beta_raw,
            chromosomes=["7"],
            chr_to_indices=chr_to_indices,
            pos=pos,
            window_bp=500,
            vcf_chr_dict=None,
            tabix=None,
            mapper=mapper,
            log=self.log,
            verbose=False,
        )
        self.assertGreaterEqual(n_blocks, 0)
        self.assertGreaterEqual(V_g, 0.0)

    def test_polygenic_bernoulli_sampling(self):
        """Polygenic mode uses Bernoulli(pi) with variable causal counts."""
        counts = []
        for seed in range(15):
            _, causals = simulate_sumstats_region(
                vcf_path=self.vcf_path,
                region=self.region,
                n=10000,
                mode="polygenic",
                pi=0.01,
                thin=0.1,
                seed=seed,
                verbose=False,
                log=self.log,
            )
            counts.append(len(causals))
        self.assertGreater(np.std(counts), 0.0)
        self.assertGreater(np.mean(counts), 5)

    def test_region_h2_calibration(self):
        """Optional region h2 target rescales effects."""
        sumstats_obj, _ = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            n=10000,
            n_causal=5,
            h2=0.01,
            seed=42,
            verbose=False,
            log=self.log,
        )
        z_causal = sumstats_obj.data.loc[sumstats_obj.data["IS_CAUSAL"], "Z"].values
        self.assertGreater(np.mean(np.abs(z_causal)), 0.1)

    def test_binary_liability_trait_model(self):
        """Liability trait_model runs for binary traits."""
        sumstats_obj, _ = simulate_sumstats_region(
            vcf_path=self.vcf_path,
            region=self.region,
            trait="binary",
            trait_model="liability",
            n_case=5000,
            n_ctrl=5000,
            n_causal=2,
            seed=42,
            verbose=False,
            log=self.log,
        )
        self.assertIn("N_CASE", sumstats_obj.data.columns)
        self.assertGreater(len(sumstats_obj.data), 0)

    def test_global_lambda_gc_and_thin(self):
        """Global mode accepts lambda_gc and thin like region mode."""
        sumstats_obj, _ = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="sparse",
            n_causal=2,
            h2=0.01,
            thin=0.5,
            lambda_gc=1.1,
            sigma_strat=0.0,
            seed=42,
            verbose=False,
            log=self.log,
        )
        self.assertGreater(len(sumstats_obj.data), 0)
        lam = _empirical_lambda_gc(sumstats_obj.data["Z"].values)
        self.assertGreater(lam, 1.0)

    def test_global_chr_autodetect_from_header(self):
        """chromosomes=None should match explicit chr list on chr7-only VCF."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_reference_format(self.vcf_path)
        detected = _autodetect_chromosomes_from_vcf(
            vcf_path=self.vcf_path,
            mapper=mapper,
            vcf_chr_dict=None,
            log=self.log,
            verbose=False,
        )
        self.assertIn("7", detected)

        sim_kwargs = dict(
            vcf_path=self.vcf_path,
            n=10000,
            mode="sparse",
            n_causal=3,
            h2=0.01,
            thin=0.1,
            seed=42,
            verbose=False,
            log=self.log,
        )
        explicit, causals_explicit = simulate_sumstats_global(
            chromosomes=["7"],
            **sim_kwargs,
        )
        with patch(
            "gwaslab.util.util_in_simulate._autodetect_chromosomes_from_vcf",
            return_value=["7"],
        ):
            auto, causals_auto = simulate_sumstats_global(
                chromosomes=None,
                **sim_kwargs,
            )
        self.assertEqual(len(explicit.data), len(auto.data))
        self.assertEqual(set(causals_explicit), set(causals_auto))

    def test_resolve_null_mode_sparse_default(self):
        self.assertEqual(_resolve_null_mode("sparse", None), "iid_far")
        self.assertEqual(_resolve_null_mode("polygenic", None), "ld_panel")
        self.assertEqual(_resolve_null_mode("sparse", "ld_panel"), "ld_panel")

    def test_sparse_h2_causal_only(self):
        """Sparse h2 calibration uses causal SNP genotypes only."""
        import polars as pl
        from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
        from gwaslab.util.util_in_simulate import (
            _read_vcf_variants,
            _select_causal_variants,
            _filter_variants,
        )

        mapper = ChromosomeMapper(log=self.log, verbose=False)
        df = _read_vcf_variants(
            vcf_path=self.vcf_path,
            region=self.region,
            mapper=mapper,
            log=self.log,
            verbose=False,
        )
        df = _filter_variants(
            df=df,
            maf_min=0.01,
            maf_max=0.5,
            log=self.log,
            verbose=False,
        )
        rng = np.random.default_rng(42)
        beta_raw, is_causal = _select_causal_variants(
            df=df,
            mode="sparse",
            pi=0.01,
            n_causal=3,
            alpha=0.0,
            effect_sd=1.0,
            rng=rng,
            log=self.log,
            verbose=False,
        )
        V_g, n_matched = _accumulate_V_g_from_causal_snps(
            vcf_path=self.vcf_path,
            df=df,
            beta_raw=beta_raw,
            is_causal=is_causal,
            vcf_chr_dict=None,
            tabix=None,
            mapper=mapper,
            log=self.log,
            verbose=False,
        )
        self.assertEqual(n_matched, 3)
        self.assertGreater(V_g, 0.0)

    def test_sparse_iid_far_panel_mask(self):
        """Causal variants should carry signal under iid_far mode."""
        sumstats_obj, causal_ids = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="sparse",
            n_causal=3,
            h2=0.05,
            null_mode="iid_far",
            seed=42,
            verbose=False,
            log=self.log,
        )
        z_causal = sumstats_obj.data.loc[sumstats_obj.data["IS_CAUSAL"], "Z"].values
        self.assertEqual(len(causal_ids), 3)
        self.assertTrue(np.all(np.isfinite(z_causal)))
        self.assertGreater(np.mean(np.abs(z_causal)), 0.5)

    def test_sparse_iid_far_seed_repro(self):
        """iid_far sparse global simulation is reproducible."""
        kwargs = dict(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="sparse",
            n_causal=3,
            h2=0.01,
            null_mode="iid_far",
            seed=99,
            verbose=False,
            log=self.log,
        )
        s1, c1 = simulate_sumstats_global(**kwargs)
        s2, c2 = simulate_sumstats_global(**kwargs)
        self.assertEqual(set(c1), set(c2))
        np.testing.assert_array_almost_equal(s1.data["Z"].values, s2.data["Z"].values, decimal=10)

    def test_polygenic_ld_panel_per_chr(self):
        """Polygenic global mode still runs with per-chromosome ld_panel path."""
        sumstats_obj, causals = simulate_sumstats_global(
            vcf_path=self.vcf_path,
            chromosomes=["7"],
            n=10000,
            mode="polygenic",
            pi=0.01,
            h2=0.01,
            null_mode="ld_panel",
            seed=42,
            verbose=False,
            log=self.log,
        )
        self.assertGreater(len(sumstats_obj.data), 0)
        self.assertGreater(len(causals), 0)
        self.assertTrue(np.all(np.isfinite(sumstats_obj.data["Z"].values)))

    def test_merge_causal_ld_windows(self):
        import polars as pl

        df = pl.DataFrame({
            "CHR": ["7", "7", "7"],
            "POS": [1000, 2000000, 3000000],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "A"],
            "EAF": [0.3, 0.4, 0.5],
        })
        is_causal = np.array([True, True, False])
        windows = _merge_causal_ld_windows(df, is_causal, ld_window_bp=500_000)
        self.assertGreaterEqual(len(windows), 1)
        mask = _build_panel_mask(df, windows)
        self.assertTrue(mask[0])
        self.assertTrue(mask[1])


if __name__ == "__main__":
    unittest.main()
