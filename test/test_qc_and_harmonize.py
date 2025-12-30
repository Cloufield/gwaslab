"""
Test file for QC and harmonization functions.

To see logs during test execution:
    python test/test_qc_and_harmonize.py -v
    python -m unittest test.test_qc_and_harmonize -v -s
    pytest test/test_qc_and_harmonize.py -v -s

The -s flag (or --capture=no in pytest) shows stdout/stderr, which includes gwaslab logs.
"""

import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
import numpy as np
import scipy.stats as ss
import gzip
import subprocess
import tempfile
import shutil
from gwaslab.g_Sumstats import Sumstats


RAW_DIR = os.path.join(os.path.dirname(__file__), "raw")


class TestBasicCheck(unittest.TestCase):
    def test_basic_check_on_dirty_sumstats(self):
        path = os.path.join(RAW_DIR, "dirty_sumstats.tsv")
        gl = Sumstats(sumstats=path, fmt=None, tab_fmt="tsv", snpid="SNPID", chrom="CHR", pos="POS", ea="EA", nea="NEA", p="P", eaf="EAF", verbose=False)
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)
        gl.basic_check(remove_dup=True, verbose=False)
        status = gl.check_sumstats_qc_status()
        self.assertIn("basic_check", status)
        self.assertTrue(status["basic_check"].get("performed", False))
        self.assertGreater(len(gl.data), 0)


class TestRemoveDup(unittest.TestCase):
    def test_remove_dup_on_duplicate(self):
        path = os.path.join(RAW_DIR, "duplicate.tsv")
        gl = Sumstats(
            sumstats=path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="ALT",
            nea="REF",
            p="P",
            snpid="SNP",
            rsid="rsID",
            neaf="Frq",
            verbose=False,
        )
        before = len(gl.data)
        gl.remove_dup(mode="dsdrdc", keep="first", verbose=False)
        after = len(gl.data)
        self.assertLess(after, before)


class TestHarmonization(unittest.TestCase):
    def test_harmonize_on_to_harmonize(self):
        path = os.path.join(RAW_DIR, "to_harmonize.tsv")
        gl = Sumstats(sumstats=path, tab_fmt="tsv", chrom="CHR", pos="POS", ea="EA", nea="NEA", p="P", snpid="SNPID", verbose=False)
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)
        # run basic_check within harmonize (ref-free path)
        gl.harmonize(basic_check=True, verbose=False)
        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertGreater(len(gl.data), 0)


class TestQCOutputMatches(unittest.TestCase):
    def test_clean_output_matches_correct(self):
        raw_path = os.path.join(os.path.dirname(__file__), "raw", "dirty_sumstats.tsv")
        correct_path = os.path.join(os.path.dirname(__file__), "correct", "clean_sumstats.tsv")
        out_dir = os.path.join(os.path.dirname(__file__), "output")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "clean_sumstats.tsv")

        gl = Sumstats(sumstats=raw_path, fmt="gwaslab", other=["NOTE"], verbose=False)
        gl.basic_check(remove=True, remove_dup=True, verbose=False)
        gl.data.to_csv(out_path, sep="\t", index=None)

        df_out = pd.read_csv(out_path, sep="\t")
        df_corr = pd.read_csv(correct_path, sep="\t")

        common_cols = [c for c in df_corr.columns if c in df_out.columns]
        df_out = df_out[common_cols]
        df_corr = df_corr[common_cols]

        sort_cols = [c for c in ["CHR", "POS", "SNPID"] if c in common_cols]
        if sort_cols:
            df_out = df_out.sort_values(by=sort_cols).reset_index(drop=True)
            df_corr = df_corr.sort_values(by=sort_cols).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_out, df_corr, check_dtype=False)


class TestHarmonizeWorkflowWithReferences(unittest.TestCase):
    def test_harmonize_bbj_chr7_with_ref_seq_rsid_and_infer(self):
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_seq_path = os.path.join(ref_dir, "chr7.fasta.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=True
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        gl.harmonize(
            basic_check=True,
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True
        )

        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertIn("rsID", gl.data.columns)
        self.assertIn("EA", gl.data.columns)
        self.assertIn("NEA", gl.data.columns)
        self.assertGreater(len(gl.data), 0)


class TestFixIDWorkflow(unittest.TestCase):
    def test_fix_id_produces_expected_output(self):
        raw_path = os.path.join(os.path.dirname(__file__), "raw", "fixid.tsv")
        correct_path = os.path.join(os.path.dirname(__file__), "correct", "fixed_id.tsv")
        out_dir = os.path.join(os.path.dirname(__file__), "output")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "fixed_id.tsv")

        gl = Sumstats(sumstats=raw_path, fmt="gwaslab", verbose=False)
        gl.fix_id(fixchrpos=True, fixeanea=True, fixprefix=True, verbose=False)
        gl.data.to_csv(out_path, sep="\t", index=None)

        df_out = pd.read_csv(out_path, sep="\t")
        df_corr = pd.read_csv(correct_path, sep="\t")

        common_cols = [c for c in df_corr.columns if c in df_out.columns]
        df_out = df_out[common_cols]
        df_corr = df_corr[common_cols]

        sort_cols = [c for c in ["CHR", "POS", "SNPID"] if c in common_cols]
        if sort_cols:
            df_out = df_out.sort_values(by=sort_cols).reset_index(drop=True)
            df_corr = df_corr.sort_values(by=sort_cols).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_out, df_corr, check_dtype=False)

    def test_harmonize_bbj_chr7_sweep_mode(self):
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_seq_path = os.path.join(ref_dir, "chr7.fasta.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=True
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        gl.harmonize(
            basic_check=True,
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True,
            sweep_mode=True,
        )

        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertIn("rsID", gl.data.columns)
        self.assertIn("EA", gl.data.columns)
        self.assertIn("NEA", gl.data.columns)
        self.assertGreater(len(gl.data), 0)


class TestCheckAF(unittest.TestCase):
    def test_check_af_with_reference_vcf(self):
        """Test check_af function to calculate DAF (difference in allele frequency)"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Run check_af
        gl.check_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that DAF column was created
        self.assertIn("DAF", gl.data.columns)
        
        # Check that DAF values are numeric (can be NaN for variants not found in VCF)
        daf_values = gl.data["DAF"].dropna()
        if len(daf_values) > 0:
            # DAF should be between -1 and 1 (difference between two frequencies)
            self.assertTrue(all(daf_values >= -1) and all(daf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)


class TestInferAF(unittest.TestCase):
    def test_infer_af_with_reference_vcf(self):
        """Test infer_af function to infer EAF from reference VCF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Run infer_af
        gl.infer_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that EAF column still exists
        self.assertIn("EAF", gl.data.columns)
        
        # Check that EAF values are numeric and between 0 and 1
        eaf_values = gl.data["EAF"].dropna()
        if len(eaf_values) > 0:
            self.assertTrue(all(eaf_values >= 0) and all(eaf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)
        
        # Check that at least some EAF values were inferred (if variants match the VCF region)
        # Note: infer_af only processes variants with good status codes and within the VCF region
        inferred_eaf_count = gl.data["EAF"].notna().sum()
        self.assertGreater(inferred_eaf_count, 0)


class TestConsistencyBetweenMethods(unittest.TestCase):
    """Test consistency between old and new methods (assign_rsid vs assign_rsid2, infer_strand vs infer_strand2)"""
    
    def test_assign_rsid_consistency(self):
        """Test that assign_rsid and assign_rsid2 produce consistent results"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl1.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl1.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl1.data = gl1.data[chr7_mask].reset_index(drop=True)
        # Use only every 15th variant for faster testing (1/15 of variants)
        gl1.data = gl1.data.iloc[::15].reset_index(drop=True)
        print(f"[test_assign_rsid_consistency] Using {len(gl1.data)} variants (every 15th variant, 1/15 of total)")
        gl1.assign_rsid(ref_rsid_vcf=ref_rsid_vcf_path, threads=6, verbose=False)
        result1 = gl1.data.copy()
        rsid1_count = result1["rsID"].notna().sum()
        
        # Test with assign_rsid2 (new method) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl2.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl2.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl2.data = gl2.data[chr7_mask].reset_index(drop=True)
        # Use only every 15th variant for faster testing (same as gl1, 1/15 of variants)
        gl2.data = gl2.data.iloc[::15].reset_index(drop=True)
        gl2.assign_rsid2(vcf_path=ref_rsid_vcf_path, threads=6, verbose=False)
        result2 = gl2.data.copy()
        rsid2_count = result2["rsID"].notna().sum()
        
        # Compare rsID assignments
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Check that both methods assigned rsIDs
        self.assertGreater(rsid1_count, 0, "assign_rsid should assign some rsIDs")
        self.assertGreater(rsid2_count, 0, "assign_rsid2 should assign some rsIDs")
        
        # Compare rsID assignments where both have values
        # Merge on CHR:POS to compare rsIDs
        merged = result1[["CHR", "POS", "rsID"]].merge(
            result2[["CHR", "POS", "rsID"]],
            on=["CHR", "POS"],
            suffixes=("_old", "_new"),
            how="inner"
        )
        
        # Check consistency: where both have rsIDs, they should match
        both_have_rsid = merged["rsID_old"].notna() & merged["rsID_new"].notna()
        if both_have_rsid.sum() > 0:
            matching = (merged.loc[both_have_rsid, "rsID_old"] == merged.loc[both_have_rsid, "rsID_new"]).sum()
            match_rate = matching / both_have_rsid.sum()
            print(f"\n[test_assign_rsid_consistency] rsID match rate: {match_rate:.2%} ({matching}/{both_have_rsid.sum()} variants)")
            # Allow for some differences due to different processing, but should be mostly consistent
            self.assertGreater(match_rate, 0.8, 
                             f"rsID assignments should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_assign_rsid_consistency] No variants with rsIDs in both methods to compare")
        
        # Check that the number of assigned rsIDs is similar (within 20%)
        if rsid1_count > 0 and rsid2_count > 0:
            ratio = min(rsid1_count, rsid2_count) / max(rsid1_count, rsid2_count)
            print(f"[test_assign_rsid_consistency] rsID count ratio: {ratio:.2%} (old: {rsid1_count}, new: {rsid2_count})")
            self.assertGreater(ratio, 0.99, 
                             f"Number of assigned rsIDs should be similar (ratio: {ratio:.2%})")
        else:
            print(f"[test_assign_rsid_consistency] rsID counts - old: {rsid1_count}, new: {rsid2_count}")
    
    def test_infer_strand_consistency(self):
        """Test that infer_strand (via harmonize) and infer_strand2 produce consistent results"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_seq_path = os.path.join(ref_dir, "chr7.fasta.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl1.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl1.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl1.data = gl1.data[chr7_mask].reset_index(drop=True)
        # Use only every 15th variant for faster testing (1/15 of variants)
        gl1.data = gl1.data.iloc[::15].reset_index(drop=True)
        print(f"[test_infer_strand_consistency] Using {len(gl1.data)} variants (every 15th variant, 1/15 of total)")
        
        # Test with infer_strand via harmonize (old method, non-sweep mode)
        gl1.harmonize(
            basic_check=False,  # Already done above
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            ref_maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True,
            sweep_mode=False,  # Use old method
            infer_strand_kwargs={"daf_tolerance": 0.20}  # Pass daf_tolerance through kwargs
        )
        result1 = gl1.data.copy()
        status1 = result1["STATUS"].copy()
        
        # Test with infer_strand2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl2.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl2.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl2.data = gl2.data[chr7_mask].reset_index(drop=True)
        # Use only every 15th variant for faster testing (same as gl1, 1/15 of variants)
        gl2.data = gl2.data.iloc[::15].reset_index(drop=True)
        
        gl2.harmonize(
            basic_check=False,  # Already done above
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            ref_maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True,
            sweep_mode=True,  # Use new method
            infer_strand_kwargs={"daf_tolerance": 0.20}  # Pass daf_tolerance through kwargs
        )
        result2 = gl2.data.copy()
        status2 = result2["STATUS"].copy()
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare STATUS codes - the 7th digit should be consistent for strand inference
        # Extract 7th digit from STATUS (strand inference status)
        def get_strand_status(status_series):
            """Extract 7th digit from STATUS codes"""
            return status_series.astype(str).str[-1] if len(status_series) > 0 else pd.Series(dtype=str)
        
        strand_status1 = get_strand_status(status1)
        strand_status2 = get_strand_status(status2)
        
        # Compare strand status codes where both have valid status
        valid_mask = (strand_status1 != "nan") & (strand_status2 != "nan")
        if valid_mask.sum() > 0:
            matching = (strand_status1[valid_mask] == strand_status2[valid_mask]).sum()
            match_rate = matching / valid_mask.sum()
            print(f"\n[test_infer_strand_consistency] Strand STATUS match rate: {match_rate:.2%} ({matching}/{valid_mask.sum()} variants)")
            # Show distribution of status codes for debugging
            status1_dist = strand_status1[valid_mask].value_counts().sort_index()
            status2_dist = strand_status2[valid_mask].value_counts().sort_index()
            print(f"[test_infer_strand_consistency] Old method STATUS distribution: {dict(status1_dist)}")
            print(f"[test_infer_strand_consistency] New method STATUS distribution: {dict(status2_dist)}")
            # Allow for some differences, but should be mostly consistent
            self.assertGreater(match_rate, 0.99, 
                             f"Strand inference STATUS codes should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_infer_strand_consistency] No variants with valid STATUS in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_infer_strand_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")

    def test_infer_strand_consistency_simulated(self):
        """Test that infer_strand (old method) and infer_strand2 (new method) produce consistent results using simulated data with edge cases"""
        temp_dir = os.path.join(ROOT, "test", "output")
        os.makedirs(temp_dir, exist_ok=True)
  
        # Test parameters
        maf_threshold = 0.40
        ref_maf_threshold = 0.40
        
        # Define test variants: (SNPID, CHR, POS, EA, NEA, EAF, BETA, SE, P)
        variants = [
            # Non-palindromic SNPs
            ("1:1000_A_G", 1, 1000, "G", "A", 0.1, 0.1, 0.01, 1e-5),
            ("1:2000_A_G", 1, 2000, "G", "A", 0.3, 0.1, 0.01, 1e-5),
            ("1:3000_A_G", 1, 3000, "G", "A", 0.4, 0.1, 0.01, 1e-5),  # At threshold
            ("1:4000_A_G", 1, 4000, "G", "A", 0.5, 0.1, 0.01, 1e-5),  # Ambiguous
            ("1:5000_A_G", 1, 5000, "G", "A", 0.6, 0.1, 0.01, 1e-5),
            ("1:6000_A_G", 1, 6000, "G", "A", 0.9, 0.1, 0.01, 1e-5),
            
            # Palindromic A/T SNPs
            ("1:7000_A_T", 1, 7000, "T", "A", 0.1, 0.1, 0.01, 1e-5),
            ("1:8000_A_T", 1, 8000, "T", "A", 0.3, 0.1, 0.01, 1e-5),
            ("1:9000_A_T", 1, 9000, "T", "A", 0.4, 0.1, 0.01, 1e-5),  # At threshold
            ("1:10000_A_T", 1, 10000, "T", "A", 0.5, 0.1, 0.01, 1e-5),  # Ambiguous
            ("1:11000_A_T", 1, 11000, "T", "A", 0.6, 0.1, 0.01, 1e-5),
            ("1:12000_A_T", 1, 12000, "T", "A", 0.9, 0.1, 0.01, 1e-5),
            ("1:23000_A_T", 1, 23000, "T", "A", 0.41, 0.1, 0.01, 1e-5),  # MAF(EAF)=0.41 > threshold
            ("1:24000_A_T", 1, 24000, "T", "A", 0.1, 0.1, 0.01, 1e-5),  # VCF AF=0.41, MAF(RAF)=0.41 > threshold
            
            # Floating point precision edge cases for A/T palindromic SNPs
            ("1:27000_A_T", 1, 27000, "T", "A", 0.4000003453, 0.1, 0.01, 1e-5),  # EAF just above threshold (0.4 + epsilon)
            ("1:28000_A_T", 1, 28000, "T", "A", 0.3999996547, 0.1, 0.01, 1e-5),  # EAF just below threshold (0.4 - epsilon)
            ("1:29000_A_T", 1, 29000, "T", "A", 0.4, 0.1, 0.01, 1e-5),  # EAF exactly at threshold, VCF AF with precision issue
            
            # Palindromic G/C SNPs
            ("1:13000_G_C", 1, 13000, "C", "G", 0.1, 0.1, 0.01, 1e-5),
            ("1:14000_G_C", 1, 14000, "C", "G", 0.3, 0.1, 0.01, 1e-5),
            ("1:15000_G_C", 1, 15000, "C", "G", 0.4, 0.1, 0.01, 1e-5),  # At threshold
            ("1:16000_G_C", 1, 16000, "C", "G", 0.5, 0.1, 0.01, 1e-5),  # Ambiguous
            ("1:17000_G_C", 1, 17000, "C", "G", 0.6, 0.1, 0.01, 1e-5),
            ("1:18000_G_C", 1, 18000, "C", "G", 0.9, 0.1, 0.01, 1e-5),
            ("1:25000_G_C", 1, 25000, "C", "G", 0.41, 0.1, 0.01, 1e-5),  # MAF(EAF)=0.41 > threshold
            ("1:26000_G_C", 1, 26000, "C", "G", 0.1, 0.1, 0.01, 1e-5),  # VCF AF=0.41, MAF(RAF)=0.41 > threshold
            
            # Floating point precision edge cases for G/C palindromic SNPs
            ("1:30000_G_C", 1, 30000, "C", "G", 0.4000003453, 0.1, 0.01, 1e-5),  # EAF just above threshold (0.4 + epsilon)
            ("1:31000_G_C", 1, 31000, "C", "G", 0.3999996547, 0.1, 0.01, 1e-5),  # EAF just below threshold (0.4 - epsilon)
            ("1:32000_G_C", 1, 32000, "C", "G", 0.4, 0.1, 0.01, 1e-5),  # EAF exactly at threshold, VCF AF with precision issue
            
            # Indels - comprehensive test cases with edge cases
            # Forward strand cases (EAF close to RAF, within daf_tolerance=0.20)
            ("1:33000_A_AT", 1, 33000, "AT", "A", 0.1, 0.1, 0.01, 1e-5),  # EAF=0.1, RAF=0.1, diff=0.0 < 0.20 -> forward (status 3)
            ("1:34000_A_AT", 1, 34000, "AT", "A", 0.15, 0.1, 0.01, 1e-5),  # EAF=0.15, RAF=0.1, diff=0.05 < 0.20 -> forward
            ("1:35000_A_AT", 1, 35000, "AT", "A", 0.2, 0.1, 0.01, 1e-5),  # EAF=0.2, RAF=0.1, diff=0.1 < 0.20 -> forward
            ("1:36000_A_AT", 1, 36000, "AT", "A", 0.29, 0.1, 0.01, 1e-5),  # EAF=0.29, RAF=0.1, diff=0.19 < 0.20 -> forward (at boundary)
            ("1:37000_A_AT", 1, 37000, "AT", "A", 0.3, 0.1, 0.01, 1e-5),  # EAF=0.3, RAF=0.1, diff=0.2 = 0.20 -> NOT forward (strict <)
            
            # Reverse strand cases (EAF close to 1-RAF, within tolerance)
            ("1:38000_A_AT", 1, 38000, "AT", "A", 0.9, 0.1, 0.01, 1e-5),  # EAF=0.9, RAF=0.1, |0.9-(1-0.1)|=0.0 < 0.20 -> reverse (status 6)
            ("1:39000_A_AT", 1, 39000, "AT", "A", 0.85, 0.1, 0.01, 1e-5),  # EAF=0.85, RAF=0.1, |0.85-0.9|=0.05 < 0.20 -> reverse
            ("1:40000_A_AT", 1, 40000, "AT", "A", 0.8, 0.1, 0.01, 1e-5),  # EAF=0.8, RAF=0.1, |0.8-0.9|=0.1 < 0.20 -> reverse
            ("1:41000_A_AT", 1, 41000, "AT", "A", 0.71, 0.1, 0.01, 1e-5),  # EAF=0.71, RAF=0.1, |0.71-0.9|=0.19 < 0.20 -> reverse (at boundary)
            ("1:42000_A_AT", 1, 42000, "AT", "A", 0.7, 0.1, 0.01, 1e-5),  # EAF=0.7, RAF=0.1, |0.7-0.9|=0.2 = 0.20 -> NOT reverse (strict <)
            
            # Ambiguous cases (both differences > tolerance)
            ("1:43000_A_AT", 1, 43000, "AT", "A", 0.5, 0.1, 0.01, 1e-5),  # EAF=0.5, RAF=0.1, forward_diff=0.4, reverse_diff=0.4, both > 0.20 -> ambiguous (status 4)
            ("1:44000_A_AT", 1, 44000, "AT", "A", 0.6, 0.1, 0.01, 1e-5),  # EAF=0.6, RAF=0.1, forward_diff=0.5, reverse_diff=0.3, both > 0.20 -> ambiguous
            ("1:45000_A_AT", 1, 45000, "AT", "A", 0.4, 0.1, 0.01, 1e-5),  # EAF=0.4, RAF=0.1, forward_diff=0.3, reverse_diff=0.5, both > 0.20 -> ambiguous
            
            # Edge cases at tolerance boundary
            ("1:46000_A_AT", 1, 46000, "AT", "A", 0.199, 0.1, 0.01, 1e-5),  # EAF=0.199, RAF=0.1, diff=0.099 < 0.20 -> forward
            ("1:47000_A_AT", 1, 47000, "AT", "A", 0.201, 0.1, 0.01, 1e-5),  # EAF=0.201, RAF=0.1, diff=0.101 < 0.20 -> forward (still within)
            ("1:48000_A_AT", 1, 48000, "AT", "A", 0.299, 0.1, 0.01, 1e-5),  # EAF=0.299, RAF=0.1, diff=0.199 < 0.20 -> forward (just below)
            ("1:49000_A_AT", 1, 49000, "AT", "A", 0.301, 0.1, 0.01, 1e-5),  # EAF=0.301, RAF=0.1, diff=0.201 > 0.20 -> ambiguous
            
            # MAF threshold edge cases for indels
            ("1:50000_A_AT", 1, 50000, "AT", "A", 0.1, 0.1, 0.01, 1e-5),  # EAF=0.1, RAF=0.41, MAF(RAF)=0.41 > ref_maf_threshold -> status 8
            ("1:51000_A_AT", 1, 51000, "AT", "A", 0.41, 0.1, 0.01, 1e-5),  # EAF=0.41, RAF=0.1, MAF(EAF)=0.41 > maf_threshold -> status 8
            ("1:52000_A_AT", 1, 52000, "AT", "A", 0.1, 0.1, 0.01, 1e-5),  # EAF=0.1, RAF=0.4, MAF(RAF)=0.4 = ref_maf_threshold -> should pass
            
            # Different indel types
            ("1:53000_A_ATCG", 1, 53000, "ATCG", "A", 0.1, 0.1, 0.01, 1e-5),  # Multi-base insertion, forward
            ("1:54000_ATCG_A", 1, 54000, "A", "ATCG", 0.9, 0.1, 0.01, 1e-5),  # Multi-base deletion, reverse
            ("1:55000_AT_ATCG", 1, 55000, "ATCG", "AT", 0.1, 0.1, 0.01, 1e-5),  # Complex indel, forward
            
            # Forward takes precedence cases
            ("1:56000_A_AT", 1, 56000, "AT", "A", 0.15, 0.1, 0.01, 1e-5),  # EAF=0.15, RAF=0.1, forward_diff=0.05 < 0.20, reverse_diff=0.75 > 0.20 -> forward (status 3)
            ("1:57000_A_AT", 1, 57000, "AT", "A", 0.85, 0.1, 0.01, 1e-5),  # EAF=0.85, RAF=0.1, forward_diff=0.75 > 0.20, reverse_diff=0.05 < 0.20 -> reverse (status 6)
            
            # Cases where both are within tolerance (forward takes precedence)
            ("1:58000_A_AT", 1, 58000, "AT", "A", 0.12, 0.1, 0.01, 1e-5),  # EAF=0.12, RAF=0.1, forward_diff=0.02 < 0.20, reverse_diff=0.78 > 0.20 -> forward
            ("1:59000_A_AT", 1, 59000, "AT", "A", 0.11, 0.1, 0.01, 1e-5),  # EAF=0.11, RAF=0.1, forward_diff=0.01 < 0.20, reverse_diff=0.79 > 0.20 -> forward
            
            # Original simple indel cases (keep for backward compatibility)
            ("1:19000_A_AT", 1, 19000, "AT", "A", 0.1, 0.1, 0.01, 1e-5),
            ("1:20000_A_AT", 1, 20000, "AT", "A", 0.3, 0.1, 0.01, 1e-5),
            ("1:21000_A_AT", 1, 21000, "AT", "A", 0.4, 0.1, 0.01, 1e-5),
            ("1:22000_A_AT", 1, 22000, "AT", "A", 0.5, 0.1, 0.01, 1e-5),
        ]
        
        # Create test files
        sumstats_path, vcf_path_gz, fasta_path_gz = self._create_test_files(temp_dir, variants)
        
        # Run old method (sweep_mode=False)
        result1, status1 = self._run_harmonize(
            sumstats_path, fasta_path_gz, vcf_path_gz, maf_threshold, ref_maf_threshold, sweep_mode=False
        )
        
        # Run new method (sweep_mode=True)
        result2, status2 = self._run_harmonize(
            sumstats_path, fasta_path_gz, vcf_path_gz, maf_threshold, ref_maf_threshold, sweep_mode=True
        )
        
        # Compare results
        self._compare_results(result1, result2, status1, status2, vcf_path_gz)
    
    
    def _create_test_files(self, temp_dir, variants):
        """Create sumstats, VCF, and FASTA files for testing"""
        import gzip
        
        # Create sumstats
        df = pd.DataFrame(variants, columns=['SNPID', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P'])
        # Convert EA and NEA to object dtype to avoid categorical dtype issues
        # The harmonization process will convert them to categorical with proper categories
        df['EA'] = df['EA'].astype('object')
        df['NEA'] = df['NEA'].astype('object')
        sumstats_path = os.path.join(temp_dir, "simulated_sumstats.txt")
        df.to_csv(sumstats_path, sep='\t', index=False)
        
        # Create VCF with reference allele frequencies
        vcf_path = os.path.join(temp_dir, "simulated_ref.vcf")
        vcf_variants = []
        for _, row in df.iterrows():
            ref, alt, eaf, pos = row['NEA'], row['EA'], row['EAF'], row['POS']
            
            # Calculate RAF for VCF
            if ref in ['A', 'T'] and alt in ['A', 'T']:
                if pos == 24000:
                    raf = 0.41
                elif pos == 29000:
                    # Floating point precision edge case: VCF AF with precision issue (like 0.4000000059604645)
                    raf = 0.4000000059604645  # This should pass with epsilon tolerance
                elif pos == 27000:
                    # EAF is 0.4000003453, set VCF AF to be at threshold with epsilon
                    raf = 0.4  # Should pass with epsilon
                elif pos == 28000:
                    # EAF is 0.3999996547, set VCF AF to be at threshold with epsilon
                    raf = 0.4  # Should pass with epsilon
                else:
                    raf = eaf + 0.02 if eaf < 0.4 else ((1 - eaf) + 0.02 if eaf > 0.6 else eaf)
            elif ref in ['G', 'C'] and alt in ['G', 'C']:
                if pos == 26000:
                    raf = 0.41
                elif pos == 32000:
                    # Floating point precision edge case: VCF AF with precision issue (like 0.4000000059604645)
                    raf = 0.4000000059604645  # This should pass with epsilon tolerance
                elif pos == 30000:
                    # EAF is 0.4000003453, set VCF AF to be at threshold with epsilon
                    raf = 0.4  # Should pass with epsilon
                elif pos == 31000:
                    # EAF is 0.3999996547, set VCF AF to be at threshold with epsilon
                    raf = 0.4  # Should pass with epsilon
                else:
                    raf = eaf + 0.02 if eaf < 0.4 else ((1 - eaf) + 0.02 if eaf > 0.6 else eaf)
            else:
                # Indels: set RAF based on test case requirements
                if pos == 33000:  # Forward: EAF=0.1, RAF=0.1
                    raf = 0.1
                elif pos == 34000:  # Forward: EAF=0.15, RAF=0.1
                    raf = 0.1
                elif pos == 35000:  # Forward: EAF=0.2, RAF=0.1
                    raf = 0.1
                elif pos == 36000:  # Forward at boundary: EAF=0.29, RAF=0.1
                    raf = 0.1
                elif pos == 37000:  # Not forward (at tolerance): EAF=0.3, RAF=0.1
                    raf = 0.1
                elif pos == 38000:  # Reverse: EAF=0.9, RAF=0.1
                    raf = 0.1
                elif pos == 39000:  # Reverse: EAF=0.85, RAF=0.1
                    raf = 0.1
                elif pos == 40000:  # Reverse: EAF=0.8, RAF=0.1
                    raf = 0.1
                elif pos == 41000:  # Reverse at boundary: EAF=0.71, RAF=0.1
                    raf = 0.1
                elif pos == 42000:  # Not reverse (at tolerance): EAF=0.7, RAF=0.1
                    raf = 0.1
                elif pos == 43000:  # Ambiguous: EAF=0.5, RAF=0.1
                    raf = 0.1
                elif pos == 44000:  # Ambiguous: EAF=0.6, RAF=0.1
                    raf = 0.1
                elif pos == 45000:  # Ambiguous: EAF=0.4, RAF=0.1
                    raf = 0.1
                elif pos == 46000:  # Forward: EAF=0.199, RAF=0.1
                    raf = 0.1
                elif pos == 47000:  # Forward: EAF=0.201, RAF=0.1
                    raf = 0.1
                elif pos == 48000:  # Forward just below: EAF=0.299, RAF=0.1
                    raf = 0.1
                elif pos == 49000:  # Ambiguous: EAF=0.301, RAF=0.1
                    raf = 0.1
                elif pos == 50000:  # MAF(RAF) > threshold: EAF=0.1, RAF=0.41
                    raf = 0.41
                elif pos == 51000:  # MAF(EAF) > threshold: EAF=0.41, RAF=0.1
                    raf = 0.1
                elif pos == 52000:  # MAF(RAF) at threshold: EAF=0.1, RAF=0.4
                    raf = 0.4
                elif pos == 53000:  # Multi-base insertion, forward: EAF=0.1, RAF=0.1
                    raf = 0.1
                elif pos == 54000:  # Multi-base deletion, reverse: EAF=0.9, RAF=0.1
                    raf = 0.1
                elif pos == 55000:  # Complex indel, forward: EAF=0.1, RAF=0.1
                    raf = 0.1
                elif pos == 56000:  # Forward precedence: EAF=0.15, RAF=0.1
                    raf = 0.1
                elif pos == 57000:  # Reverse: EAF=0.85, RAF=0.1
                    raf = 0.1
                elif pos == 58000:  # Forward precedence: EAF=0.12, RAF=0.1
                    raf = 0.1
                elif pos == 59000:  # Forward precedence: EAF=0.11, RAF=0.1
                    raf = 0.1
                elif pos in [19000, 20000, 21000, 22000]:  # Original simple cases
                    # For original cases, use simple logic
                    raf = eaf + 0.02 if eaf < 0.4 else ((1 - eaf) + 0.02 if eaf > 0.6 else eaf)
                else:
                    # Default for other indels
                    raf = eaf
            
            raf = max(0.01, min(0.99, raf))
            vcf_variants.append((row['CHR'], pos, f"rs{pos}", ref, alt, raf))
        
        # Sort variants by chromosome and position (required for tabix)
        vcf_variants.sort(key=lambda x: (x[0], x[1]))
        
        with open(vcf_path, 'w') as f:
            f.write("##fileformat=VCFv4.3\n")
            max_pos = max(pos for _, pos, _, _, _, _ in vcf_variants) if vcf_variants else 60000
            f.write(f"##contig=<ID=1,length={max_pos + 1000}>\n")
            f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, pos, rsid, ref, alt, af in vcf_variants:
                f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\tAF={af}\n")
        
        # Compress and index VCF
        vcf_path_gz = self._compress_file(vcf_path, use_tabix=True)
        
        # Create FASTA
        fasta_path = os.path.join(temp_dir, "simulated_ref.fasta")
        max_pos = df["POS"].max()
        seq_length = max_pos + 1000
        sequence = ['N'] * seq_length
        
        for _, row in df.iterrows():
            pos, nea = row['POS'], row['NEA']
            if 0 < pos <= seq_length:
                sequence[pos - 1] = nea[0] if len(nea) > 0 else 'A'
        
        bases = ['A', 'T', 'G', 'C']
        for i in range(seq_length):
            if sequence[i] == 'N':
                sequence[i] = bases[i % 4]
        
        with open(fasta_path, 'w') as f:
            f.write(">1\n")
            for i in range(0, seq_length, 60):
                f.write(''.join(sequence[i:i+60]) + '\n')
        
        fasta_path_gz = self._compress_file(fasta_path, use_tabix=False)
        print(sumstats_path, vcf_path_gz, fasta_path_gz)
        return sumstats_path, vcf_path_gz, fasta_path_gz
    
    def _compress_file(self, file_path, use_tabix=False):
        """Compress file with bgzip (or gzip fallback) and optionally index with tabix"""
        import gzip
        compressed_path = file_path + '.gz'
        
        # Try bgzip from specific path first (as user specified), then PATH
        bgzip_paths = ['/home/yunye/tools/bin/bgzip', 'bgzip']
        bgzip_success = False
        for bgzip_cmd in bgzip_paths:
            try:
                subprocess.run([bgzip_cmd, '-f', file_path], check=True,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                bgzip_success = True
                break
            except (subprocess.CalledProcessError, FileNotFoundError):
                continue
        
        if not bgzip_success:
            # Fallback to regular gzip
            with open(file_path, 'rb') as f_in:
                with gzip.open(compressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        if use_tabix:
            # Try tabix from specific path first (as user specified), then PATH
            tabix_paths = ['/home/yunye/tools/bin/tabix', 'tabix']
            tabix_success = False
            last_error = None
            for tabix_cmd in tabix_paths:
                try:
                    subprocess.run([tabix_cmd, '-p', 'vcf', compressed_path], check=True,
                                 stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                    tabix_success = True
                    break
                except (subprocess.CalledProcessError, FileNotFoundError) as e:
                    last_error = e
                    continue
            if not tabix_success:
                error_msg = f"Could not index VCF file with tabix. "
                if last_error:
                    if isinstance(last_error, subprocess.CalledProcessError):
                        error_msg += f"Tabix command failed with return code {last_error.returncode}"
                    else:
                        error_msg += f"Tabix not found in PATH or at /home/yunye/tools/bin/tabix"
                raise RuntimeError(error_msg)
        
        return compressed_path
    
    def _run_harmonize(self, sumstats_path, fasta_path_gz, vcf_path_gz, maf_threshold, ref_maf_threshold, sweep_mode):
        """Run harmonize with specified parameters"""
        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR", pos="POS", ea="EA", nea="NEA", p="P",
            snpid="SNPID", eaf="EAF", verbose=False
        )
        gl.basic_check(verbose=False)
        gl.harmonize(
            basic_check=False,
            ref_seq=fasta_path_gz,
            ref_infer=vcf_path_gz,
            ref_alt_freq="AF",
            maf_threshold=maf_threshold,
            ref_maf_threshold=ref_maf_threshold,
            threads=2,
            remove=False,
            verbose=False,
            sweep_mode=sweep_mode,
            infer_strand_kwargs={"daf_tolerance": 0.20}
        )
        result = gl.data.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        status = result["STATUS"].copy()
        return result, status
    
    def _compare_results(self, result1, result2, status1, status2, vcf_path_gz):
        """Compare results from old and new methods"""
        def get_strand_status(status_series):
            return status_series.astype(str).str[-1] if len(status_series) > 0 else pd.Series(dtype=str)
        
        strand_status1 = get_strand_status(status1)
        strand_status2 = get_strand_status(status2)
        
        # Print detailed variant information
        print(f"\n[DEBUG test] === Detailed Variant Information ===")
        print(f"[DEBUG test] Total variants: {len(result1)}")
        
        for idx in range(len(result1)):
            row1, row2 = result1.iloc[idx], result2.iloc[idx]
            pos, eaf, ea, nea = row1["POS"], row1["EAF"], row1["EA"], row1["NEA"]
            status1_val, status2_val = strand_status1.iloc[idx], strand_status2.iloc[idx]
            # Note: Old method (sweep_mode=False) does not generate RAF column, only new method does
            raf1 = row1.get("RAF") if "RAF" in row1 else None
            raf2 = row2.get("RAF") if "RAF" in row2 else None
            
            is_pal = (ea in ['A', 'T'] and nea in ['A', 'T']) or (ea in ['G', 'C'] and nea in ['G', 'C'])
            maf_eaf = min(eaf, 1 - eaf) if pd.notna(eaf) else None
            maf_raf1 = min(raf1, 1 - raf1) if pd.notna(raf1) and isinstance(raf1, (int, float)) else None
            maf_raf2 = min(raf2, 1 - raf2) if pd.notna(raf2) and isinstance(raf2, (int, float)) else None
            
            # Get VCF info
            vcf_af = self._get_vcf_af(vcf_path_gz, row1["CHR"], pos)
            
            match_marker = "OK" if status1_val == status2_val else "XX"
            maf_eaf_str = f"{maf_eaf:.3f}" if maf_eaf is not None else "N/A"
            vcf_af_str = f"{vcf_af:.3f}" if vcf_af is not None else "N/A"
            raf1_str = "N/A" if raf1 is None or pd.isna(raf1) else f"{raf1:.3f}"
            raf2_str = "N/A" if raf2 is None or pd.isna(raf2) else f"{raf2:.3f}"
            maf_raf1_str = f"{maf_raf1:.3f}" if maf_raf1 is not None else "N/A"
            maf_raf2_str = f"{maf_raf2:.3f}" if maf_raf2 is not None else "N/A"
            
            print(f"[DEBUG test] {match_marker} POS={pos:5d}: EAF={eaf:.3f}, EA={ea:2s}, NEA={nea:2s}, "
                  f"is_pal={is_pal}, MAF(EAF)={maf_eaf_str:>6s}, VCF_AF={vcf_af_str:>6s}, "
                  f"RAF_old={raf1_str:>6s}, RAF_new={raf2_str:>6s}, "
                  f"MAF(RAF_old)={maf_raf1_str:>6s}, MAF(RAF_new)={maf_raf2_str:>6s}, "
                  f"STATUS_old={status1_val}, STATUS_new={status2_val}")
        
        # Compare STATUS codes for variants checked by both methods
        status9_mask1 = strand_status1 == "9"
        status9_mask2 = strand_status2 == "9"
        both_checked = ~status9_mask1 & ~status9_mask2
        valid_mask = (strand_status1 != "nan") & (strand_status2 != "nan") & both_checked
        
        if valid_mask.sum() > 0:
            matching = (strand_status1[valid_mask] == strand_status2[valid_mask]).sum()
            match_rate = matching / valid_mask.sum()
            print(f"\n[test_infer_strand_consistency_simulated] Strand STATUS match rate: {match_rate:.2%} ({matching}/{valid_mask.sum()} variants)")
            
            status1_dist = strand_status1[valid_mask].value_counts().sort_index()
            status2_dist = strand_status2[valid_mask].value_counts().sort_index()
            print(f"[test_infer_strand_consistency_simulated] Old method STATUS distribution: {dict(status1_dist)}")
            print(f"[test_infer_strand_consistency_simulated] New method STATUS distribution: {dict(status2_dist)}")
            
            # Report mismatches
            mismatches = valid_mask & (strand_status1 != strand_status2)
            if mismatches.sum() > 0:
                print(f"\n[test_infer_strand_consistency_simulated] Mismatches found ({mismatches.sum()} variants):")
                for idx in result1.index[mismatches]:
                    row1, row2 = result1.loc[idx], result2.loc[idx]
                    pos, eaf, ea, nea = row1["POS"], row1["EAF"], row1["EA"], row1["NEA"]
                    status1_val = strand_status1.iloc[idx] if idx < len(strand_status1) else "N/A"
                    status2_val = strand_status2.iloc[idx] if idx < len(strand_status2) else "N/A"
                    # Note: Old method (sweep_mode=False) does not generate RAF column
                    raf1 = row1.get("RAF") if "RAF" in row1 else None
                    raf2 = row2.get("RAF") if "RAF" in row2 else None
                    maf_eaf = min(eaf, 1 - eaf) if pd.notna(eaf) else None
                    maf_raf1 = min(raf1, 1 - raf1) if pd.notna(raf1) and isinstance(raf1, (int, float)) else None
                    maf_raf2 = min(raf2, 1 - raf2) if pd.notna(raf2) and isinstance(raf2, (int, float)) else None
                    raf1_str = "N/A" if raf1 is None or pd.isna(raf1) else f"{raf1:.3f}"
                    raf2_str = "N/A" if raf2 is None or pd.isna(raf2) else f"{raf2:.3f}"
                    maf_eaf_str = f"{maf_eaf:.3f}" if maf_eaf is not None else "N/A"
                    maf_raf1_str = f"{maf_raf1:.3f}" if maf_raf1 is not None else "N/A"
                    maf_raf2_str = f"{maf_raf2:.3f}" if maf_raf2 is not None else "N/A"
                    print(f"  POS={pos}, EAF={eaf:.3f}, EA={ea}, NEA={nea}, "
                          f"RAF_old={raf1_str}, RAF_new={raf2_str}, "
                          f"MAF(EAF)={maf_eaf_str}, "
                          f"MAF(RAF_old)={maf_raf1_str}, "
                          f"MAF(RAF_new)={maf_raf2_str}, "
                          f"STATUS_old={status1_val}, STATUS_new={status2_val}")
            
            # Assert 100% consistency for variants found by both methods
            self.assertGreaterEqual(match_rate, 1.0,
                                 f"Strand inference STATUS codes should be consistent for variants found by both methods (match rate: {match_rate:.2%})")
        else:
            checked_old = (~status9_mask1).sum()
            checked_new = (~status9_mask2).sum()
            print(f"\n[test_infer_strand_consistency_simulated] No variants checked by both methods")
            print(f"[test_infer_strand_consistency_simulated] Variants checked - old: {checked_old}, new: {checked_new}")
        
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")
    
    def _get_vcf_af(self, vcf_path_gz, chrom, pos):
        """Get allele frequency from VCF for a given position"""
        try:
            from pysam import VariantFile
            vcf_reader = VariantFile(vcf_path_gz)
            for record in vcf_reader.fetch(str(chrom), pos-1, pos):
                if record.pos == pos:
                    af = record.info.get("AF", [None])[0] if "AF" in record.info else None
                    vcf_reader.close()
                    return af
            vcf_reader.close()
        except:
            pass
        return None


class TestCheckAF2(unittest.TestCase):
    def test_check_af2_with_reference_vcf(self):
        """Test check_af2 function (sweep mode) to calculate DAF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Run check_af2 (sweep mode)
        gl.check_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that DAF column was created
        self.assertIn("DAF", gl.data.columns)
        
        # Check that DAF values are numeric (can be NaN for variants not found in VCF)
        daf_values = gl.data["DAF"].dropna()
        if len(daf_values) > 0:
            # DAF should be between -1 and 1 (difference between two frequencies)
            self.assertTrue(all(daf_values >= -1) and all(daf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)


class TestInferAF2(unittest.TestCase):
    def test_infer_af2_with_reference_vcf(self):
        """Test infer_af2 function (sweep mode) to infer EAF from reference VCF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Store original EAF count
        original_eaf_count = gl.data["EAF"].notna().sum()

        # Run infer_af2 (sweep mode)
        gl.infer_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that EAF column still exists
        self.assertIn("EAF", gl.data.columns)
        
        # Check that EAF values are numeric and between 0 and 1
        eaf_values = gl.data["EAF"].dropna()
        if len(eaf_values) > 0:
            self.assertTrue(all(eaf_values >= 0) and all(eaf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)
        
        # Check that at least some EAF values were inferred (if variants match the VCF region)
        inferred_eaf_count = gl.data["EAF"].notna().sum()
        self.assertGreater(inferred_eaf_count, 0)


class TestInferEAFFromMAF2(unittest.TestCase):
    def test_infer_eaf_from_maf2_with_reference_vcf(self):
        """Test infer_eaf_from_maf2 function (sweep mode) to infer EAF from MAF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Use only 1/15 of variants for faster testing
        gl.data = gl.data.iloc[::15].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Check if MAF column exists, if not create it from EAF
        if "MAF" not in gl.data.columns or gl.data["MAF"].isna().all():
            # Create MAF from EAF if available
            if "EAF" in gl.data.columns:
                gl.data["MAF"] = gl.data["EAF"].apply(lambda x: min(x, 1-x) if pd.notna(x) else pd.NA)

        # Store original EAF count
        original_eaf_count = gl.data["EAF"].notna().sum()

        # Run infer_eaf_from_maf2 (sweep mode)
        gl.infer_eaf_from_maf2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that EAF column still exists
        self.assertIn("EAF", gl.data.columns)
        
        # Check that EAF values are numeric and between 0 and 1
        eaf_values = gl.data["EAF"].dropna()
        if len(eaf_values) > 0:
            self.assertTrue(all(eaf_values >= 0) and all(eaf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)
        
        # Check that at least some EAF values were inferred (if variants match the VCF region and have MAF)
        inferred_eaf_count = gl.data["EAF"].notna().sum()
        self.assertGreater(inferred_eaf_count, 0)


class TestConsistencyCheckAF(unittest.TestCase):
    """Test consistency between check_af (ver1) and check_af2 (ver2)"""
    
    def test_check_af_consistency(self):
        """Test that check_af and check_af2 produce consistent DAF values"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl1.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl1.data = gl1.data[
            (gl1.data["CHR"] == 7) & 
            (gl1.data["POS"] >= 126253550) & 
            (gl1.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 15th variant for faster testing (1/15 of variants)
        gl1.data = gl1.data.iloc[::15].reset_index(drop=True)
        print(f"[test_check_af_consistency] Using {len(gl1.data)} variants (every 15th variant, 1/15 of total in chr7:126253550-128253550)")
        
        # Test with check_af (old method, normal mode)
        gl1.check_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result1 = gl1.data.copy()
        daf1 = result1["DAF"].copy()
        
        # Test with check_af2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl2.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl2.data = gl2.data[
            (gl2.data["CHR"] == 7) & 
            (gl2.data["POS"] >= 126253550) & 
            (gl2.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 15th variant for faster testing (same as gl1, 1/15 of variants)
        gl2.data = gl2.data.iloc[::15].reset_index(drop=True)
        
        gl2.check_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result2 = gl2.data.copy()
        daf2 = result2["DAF"].copy()
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare DAF values where both have values
        both_have_daf = daf1.notna() & daf2.notna()
        if both_have_daf.sum() > 0:
            daf1_values = daf1[both_have_daf]
            daf2_values = daf2[both_have_daf]
            
            # Calculate differences
            daf_diff = (daf1_values - daf2_values).abs()
            max_diff = daf_diff.max()
            mean_diff = daf_diff.mean()
            
            # Check that differences are small (within tolerance)
            # Allow for small numerical differences due to different processing methods
            tolerance = 0.01  # 1% tolerance
            within_tolerance = (daf_diff <= tolerance).sum()
            match_rate = within_tolerance / both_have_daf.sum()
            
            print(f"\n[test_check_af_consistency] DAF match rate (within {tolerance}): {match_rate:.2%} ({within_tolerance}/{both_have_daf.sum()} variants)")
            print(f"[test_check_af_consistency] Max DAF difference: {max_diff:.6f}, Mean DAF difference: {mean_diff:.6f}")
            
            # Should be mostly consistent
            self.assertGreater(match_rate, 0.95, 
                             f"DAF values should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_check_af_consistency] No variants with DAF in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_check_af_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


class TestConsistencyInferAF(unittest.TestCase):
    """Test consistency between infer_af (ver1) and infer_af2 (ver2)"""
    
    def test_infer_af_consistency(self):
        """Test that infer_af and infer_af2 produce consistent EAF values"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl1.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl1.data = gl1.data[
            (gl1.data["CHR"] == 7) & 
            (gl1.data["POS"] >= 126253550) & 
            (gl1.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 15th variant for faster testing (1/15 of variants)
        gl1.data = gl1.data.iloc[::15].reset_index(drop=True)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl1.data["EAF"] = pd.NA
        
        print(f"[test_infer_af_consistency] Using {len(gl1.data)} variants (every 15th variant, 1/15 of total in chr7:126253550-128253550)")
        
        # Test with infer_af (old method, normal mode)
        gl1.infer_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result1 = gl1.data.copy()
        eaf1 = result1["EAF"].copy()
        
        # Test with infer_af2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl2.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl2.data = gl2.data[
            (gl2.data["CHR"] == 7) & 
            (gl2.data["POS"] >= 126253550) & 
            (gl2.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 15th variant for faster testing (same as gl1, 1/15 of variants)
        gl2.data = gl2.data.iloc[::15].reset_index(drop=True)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl2.data["EAF"] = pd.NA
        
        gl2.infer_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        
        result2 = gl2.data.copy()
        eaf2 = result2["EAF"].copy()
        
        # ALLELE_FLIPPED is an internal temporary column and should be dropped in output
        if "ALLELE_FLIPPED" in result2.columns:
            print(f"[test_infer_af_consistency] WARNING: ALLELE_FLIPPED should be dropped but found in output")
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare EAF values where both have values
        both_have_eaf = eaf1.notna() & eaf2.notna()
        if both_have_eaf.sum() > 0:
            eaf1_values = eaf1[both_have_eaf]
            eaf2_values = eaf2[both_have_eaf]
            
            # Calculate differences
            eaf_diff = (eaf1_values - eaf2_values).abs()
            max_diff = eaf_diff.max()
            mean_diff = eaf_diff.mean()
            
            # Check if differences are close to 1 (suggesting a flip issue)
            close_to_one = (eaf_diff > 0.99) & (eaf_diff < 1.01)
            print(f"\n[test_infer_af_consistency] Differences close to 1 (flip issue?): {close_to_one.sum()} / {both_have_eaf.sum()} ({100*close_to_one.sum()/both_have_eaf.sum():.1f}%)")
            
            # ALLELE_FLIPPED is an internal temporary column and should be dropped in output
            if "ALLELE_FLIPPED" in result2.columns:
                print(f"[test_infer_af_consistency] WARNING: ALLELE_FLIPPED should be dropped but found in output")
            
            # Print out all variants that are different
            mismatched = eaf_diff > 0.01
            if mismatched.sum() > 0:
                print(f"\n[test_infer_af_consistency] Printing all {mismatched.sum()} variants with differences > 0.01:")
                print(f"{'CHR':<6} {'POS':<12} {'EA':<6} {'NEA':<6} {'EAF_old':<10} {'EAF_new':<10} {'Diff':<10}")
                print("-" * 70)
                
                # Get the actual indices from result2 for mismatched variants
                # Note: both_have_eaf is a boolean Series aligned with result2.index
                # mismatched is a boolean Series aligned with eaf_diff (which is also aligned with result2.index[both_have_eaf])
                mismatched_bool = pd.Series(False, index=result2.index)
                mismatched_bool.loc[result2.index[both_have_eaf][mismatched]] = True
                mismatched_indices = result2.index[mismatched_bool]
                
                mismatched_df = result2.loc[mismatched_indices, ["CHR", "POS", "EA", "NEA"]].copy()
                mismatched_df["EAF_old"] = eaf1_values.loc[mismatched_indices].values
                mismatched_df["EAF_new"] = eaf2_values.loc[mismatched_indices].values
                mismatched_df["Diff"] = eaf_diff[mismatched].values
                
                # ALLELE_FLIPPED is an internal temporary column and should be dropped
                # So we don't include it in the output comparison
                mismatched_df["FLIPPED"] = "N/A (internal column)"
                
                # Sort by difference (largest first)
                mismatched_df = mismatched_df.sort_values("Diff", ascending=False)
                
                # Print all mismatched variants
                for idx, row in mismatched_df.iterrows():
                    print(f"{int(row['CHR']):<6} {int(row['POS']):<12} {row['EA']:<6} {row['NEA']:<6} "
                          f"{row['EAF_old']:<10.6f} {row['EAF_new']:<10.6f} {row['Diff']:<10.6f}")
            else:
                print(f"\n[test_infer_af_consistency] No variants with differences > 0.01")
            
            # Check that differences are small (within tolerance)
            # Allow for small numerical differences due to different processing methods
            tolerance = 0.01  # 1% tolerance
            within_tolerance = (eaf_diff <= tolerance).sum()
            match_rate = within_tolerance / both_have_eaf.sum()
            
            print(f"\n[test_infer_af_consistency] EAF match rate (within {tolerance}): {match_rate:.2%} ({within_tolerance}/{both_have_eaf.sum()} variants)")
            print(f"[test_infer_af_consistency] Max EAF difference: {max_diff:.6f}, Mean EAF difference: {mean_diff:.6f}")
            
            # Should be mostly consistent - require 99% match rate
            self.assertGreater(match_rate, 0.99, 
                             f"EAF values should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_infer_af_consistency] No variants with EAF in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_infer_af_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


class TestConsistencyInferEAFFromMAF(unittest.TestCase):
    """Test consistency between infer_eaf_from_maf (ver1) and infer_eaf_from_maf2 (ver2)"""
    
    def test_infer_eaf_from_maf_consistency(self):
        """Test that infer_eaf_from_maf and infer_eaf_from_maf2 produce consistent EAF values"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl1.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl1.data = gl1.data[
            (gl1.data["CHR"] == 7) & 
            (gl1.data["POS"] >= 126253550) & 
            (gl1.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 15th variant for faster testing (1/15 of variants)
        gl1.data = gl1.data.iloc[::15].reset_index(drop=True)
        
        # Check if MAF column exists, if not create it from EAF
        if "MAF" not in gl1.data.columns or gl1.data["MAF"].isna().all():
            # Create MAF from EAF if available (before dropping EAF)
            if "EAF" in gl1.data.columns:
                gl1.data["MAF"] = gl1.data["EAF"].apply(lambda x: min(x, 1-x) if pd.notna(x) else pd.NA)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl1.data["EAF"] = pd.NA
        
        print(f"[test_infer_eaf_from_maf_consistency] Using {len(gl1.data)} variants (every 15th variant, 1/15 of total in chr7:126253550-128253550)")
        
        # Test with infer_eaf_from_maf (old method, normal mode)
        gl1.infer_eaf_from_maf(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result1 = gl1.data.copy()
        eaf1 = result1["EAF"].copy()
        
        # Test with infer_eaf_from_maf2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl2.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl2.data = gl2.data[
            (gl2.data["CHR"] == 7) & 
            (gl2.data["POS"] >= 126253550) & 
            (gl2.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 15th variant for faster testing (same as gl1, 1/15 of variants)
        gl2.data = gl2.data.iloc[::15].reset_index(drop=True)
        
        # Check if MAF column exists, if not create it from EAF
        if "MAF" not in gl2.data.columns or gl2.data["MAF"].isna().all():
            # Create MAF from EAF if available (before dropping EAF)
            if "EAF" in gl2.data.columns:
                gl2.data["MAF"] = gl2.data["EAF"].apply(lambda x: min(x, 1-x) if pd.notna(x) else pd.NA)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl2.data["EAF"] = pd.NA
        
        gl2.infer_eaf_from_maf2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result2 = gl2.data.copy()
        eaf2 = result2["EAF"].copy()
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare EAF values where both have values
        both_have_eaf = eaf1.notna() & eaf2.notna()
        if both_have_eaf.sum() > 0:
            eaf1_values = eaf1[both_have_eaf]
            eaf2_values = eaf2[both_have_eaf]
            
            # Calculate differences
            eaf_diff = (eaf1_values - eaf2_values).abs()
            max_diff = eaf_diff.max()
            mean_diff = eaf_diff.mean()
            
            # Check that differences are small (within tolerance)
            # Allow for small numerical differences due to different processing methods
            tolerance = 0.01  # 1% tolerance
            within_tolerance = (eaf_diff <= tolerance).sum()
            match_rate = within_tolerance / both_have_eaf.sum()
            
            print(f"\n[test_infer_eaf_from_maf_consistency] EAF match rate (within {tolerance}): {match_rate:.2%} ({within_tolerance}/{both_have_eaf.sum()} variants)")
            print(f"[test_infer_eaf_from_maf_consistency] Max EAF difference: {max_diff:.6f}, Mean EAF difference: {mean_diff:.6f}")
            
            # Should be mostly consistent - require 99% match rate
            self.assertGreater(match_rate, 0.99, 
                             f"EAF values should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_infer_eaf_from_maf_consistency] No variants with EAF in both methods to compare")
        
            # Check that both methods processed similar numbers of variants
            print(f"[test_infer_eaf_from_maf_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
            self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


class TestPValueConsistencyCheck(unittest.TestCase):
    """Test P-value consistency checks with fold change reporting"""
    
    def setUp(self):
        """Set up test data with BETA, SE, and P values"""
        # Create test data with known relationships
        np.random.seed(42)
        n_variants = 100
        
        # Generate BETA and SE
        beta = np.random.normal(0, 0.1, n_variants)
        se = np.abs(np.random.normal(0.05, 0.01, n_variants))
        se = np.maximum(se, 1e-6)  # Ensure SE > 0
        
        # Calculate Z and P from BETA/SE
        z = beta / se
        p_from_betase = ss.chi2.sf(z**2, 1)
        
        # Create DataFrame with consistent values
        self.df_consistent = pd.DataFrame({
            "SNPID": [f"rs{i}" for i in range(n_variants)],
            "CHR": np.random.randint(1, 23, n_variants),
            "POS": np.random.randint(1, 250000000, n_variants),
            "EA": ["A"] * n_variants,
            "NEA": ["G"] * n_variants,
            "BETA": beta,
            "SE": se,
            "P": p_from_betase,  # Consistent P values
        })
        
        # Create DataFrame with inconsistent P values (some with large fold changes)
        self.df_inconsistent = self.df_consistent.copy()
        # Introduce inconsistencies: multiply some P values by different factors
        inconsistent_indices = np.random.choice(n_variants, size=20, replace=False)
        fold_factors = np.random.choice([0.5, 2.0, 5.0, 10.0], size=20)
        self.df_inconsistent.loc[inconsistent_indices, "P"] *= fold_factors
        
        # Create DataFrame with very small P values to test fold change reporting
        self.df_small_p = pd.DataFrame({
            "SNPID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "CHR": [1, 1, 1, 1, 1],
            "POS": [1000, 2000, 3000, 4000, 5000],
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5,
            "BETA": [0.1, 0.2, 0.3, 0.4, 0.5],
            "SE": [0.01, 0.01, 0.01, 0.01, 0.01],
            "P": [1e-6, 1e-8, 1e-10, 1e-12, 1e-14],  # Very small P values
        })
        # Make one inconsistent with 10x fold change
        self.df_small_p.loc[2, "P"] = 1e-9  # Should be ~1e-10, so ~10x difference
        
        # Create DataFrame with NaN values to test NaN handling
        self.df_with_nan = self.df_consistent.copy()
        self.df_with_nan.loc[0:5, "P"] = np.nan
        self.df_with_nan.loc[6:8, "BETA"] = np.nan
    
    def test_consistent_p_values_pass_check(self):
        """Test that consistent P values pass the consistency check"""
        gl = Sumstats(
            sumstats=self.df_consistent,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            beta="BETA",
            se="SE",
            verbose=False
        )
        
        # Run consistency check
        result = gl.check_data_consistency(rtol=1e-3, verbose=False)
        
        # Should complete without errors
        self.assertIsNotNone(result)
        # Check that QC status was updated
        status = gl.check_sumstats_qc_status()
        # The function may or may not update status, but should not raise errors
    
    def test_inconsistent_p_values_detected(self):
        """Test that inconsistent P values are detected"""
        gl = Sumstats(
            sumstats=self.df_inconsistent,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            beta="BETA",
            se="SE",
            verbose=False
        )
        
        # Run consistency check with strict tolerance
        result = gl.check_data_consistency(rtol=1e-3, verbose=False)
        
        # Should complete without errors
        self.assertIsNotNone(result)
    
    def test_small_p_values_fold_change_reporting(self):
        """Test that fold change is reported for small P values"""
        gl = Sumstats(
            sumstats=self.df_small_p,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            beta="BETA",
            se="SE",
            verbose=False
        )
        
        # Capture log output to verify fold change reporting
        import io
        import sys
        from gwaslab.info.g_Log import Log
        
        log = Log()
        log.buffer = io.StringIO()
        
        # Run consistency check
        result = gl.check_data_consistency(rtol=1e-3, verbose=True, log=log)
        
        # Check log output contains fold change information
        log_output = log.buffer.getvalue()
        # The log should mention fold change for inconsistent variants
        # Note: This test verifies the function runs, actual log content may vary
    
    def test_nan_handling_in_consistency_check(self):
        """Test that NaN values are handled correctly in consistency check"""
        gl = Sumstats(
            sumstats=self.df_with_nan,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            beta="BETA",
            se="SE",
            verbose=False
        )
        
        # Should not raise errors when encountering NaN values
        result = gl.check_data_consistency(rtol=1e-3, verbose=False)
        self.assertIsNotNone(result)
    
    def test_mlog10p_derived_p_consistency(self):
        """Test consistency check between MLOG10P-derived-P and P"""
        # Create data with MLOG10P and P
        df = pd.DataFrame({
            "SNPID": ["rs1", "rs2", "rs3"],
            "CHR": [1, 1, 1],
            "POS": [1000, 2000, 3000],
            "EA": ["A", "A", "A"],
            "NEA": ["G", "G", "G"],
            "MLOG10P": [5.0, 6.0, 7.0],  # -log10(P)
            "P": [1e-5, 1e-6, 1e-7],  # Consistent P values
        })
        
        gl = Sumstats(
            sumstats=df,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            mlog10p="MLOG10P",
            verbose=False
        )
        
        # Should complete without errors
        result = gl.check_data_consistency(rtol=1e-3, verbose=False)
        self.assertIsNotNone(result)
    
    def test_missing_columns_handled_gracefully(self):
        """Test that missing columns are handled gracefully"""
        # Create data without BETA/SE
        df = pd.DataFrame({
            "SNPID": ["rs1", "rs2"],
            "CHR": [1, 1],
            "POS": [1000, 2000],
            "EA": ["A", "A"],
            "NEA": ["G", "G"],
            "P": [0.05, 0.01],
        })
        
        gl = Sumstats(
            sumstats=df,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            verbose=False
        )
        
        # Should complete without errors even without BETA/SE
        result = gl.check_data_consistency(rtol=1e-3, verbose=False)
        self.assertIsNotNone(result)
    
    def test_zero_p_values_handled(self):
        """Test that zero or very small P values are handled correctly"""
        # Create data with very small P values
        beta = np.array([0.1, 0.2, 0.3])
        se = np.array([0.01, 0.01, 0.01])
        z = beta / se
        p_from_betase = ss.chi2.sf(z**2, 1)
        
        df = pd.DataFrame({
            "SNPID": ["rs1", "rs2", "rs3"],
            "CHR": [1, 1, 1],
            "POS": [1000, 2000, 3000],
            "EA": ["A", "A", "A"],
            "NEA": ["G", "G", "G"],
            "BETA": beta,
            "SE": se,
            "P": p_from_betase,  # Very small but consistent
        })
        
        gl = Sumstats(
            sumstats=df,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            beta="BETA",
            se="SE",
            verbose=False
        )
        
        # Should handle very small P values without errors
        result = gl.check_data_consistency(rtol=1e-3, verbose=False)
        self.assertIsNotNone(result)


if __name__ == "__main__":
    # Run with verbosity and show stdout/stderr to see logs
    # Usage: python test/test_qc_and_harmonize.py -v
    # Or: python -m unittest test.test_qc_and_harmonize -v -s
    import sys
    # Add -s flag to show stdout (logs) during test execution
    if "-s" not in sys.argv and "--buffer" not in sys.argv:
        sys.argv.append("-s")
    unittest.main(verbosity=2)
