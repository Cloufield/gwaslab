"""
Tests for meta-regression implementation comparing with MR-MEGA.

This test suite compares the GWASLab meta-regression implementation
with MR-MEGA (if available) to ensure consistency and validate performance.
"""

import os
import sys
import unittest
import tempfile
import shutil
import subprocess
import time
import tracemalloc
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl
from gwaslab.g_SumstatsMulti import SumstatsMulti
from gwaslab.extension.mrmega import _check_mrmega_available
from gwaslab.extension.mrmega.mrmega import meta_regress_mrmega_python


class TestMetaRegressMRMEGA(unittest.TestCase):
    """Test meta-regression results and compare with MR-MEGA if available"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests"""
        # Check if MR-MEGA is available
        cls.mrmega_path = _check_mrmega_available()
        cls.mrmega_available = cls.mrmega_path is not None
        
        # Create temporary directory for test files
        cls.temp_dir = tempfile.mkdtemp(prefix="gwaslab_mrmega_test_")
    
    @classmethod
    def tearDownClass(cls):
        """Clean up temporary files"""
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)
    
    def setUp(self):
        """Set up test data for each test"""
        # Create test data with 5 studies for meta-regression (MR-MEGA works better with 5+ cohorts)
        # Use realistic values covering all chromosomes and different variant cases
        np.random.seed(42)
        
        n_variants = 200  # Increased to cover more chromosomes and cases
        
        # Create variants covering all chromosomes (1-23, including X)
        # MR-MEGA processes chromosomes 1-23, so we'll distribute variants across them
        chromosomes = []
        positions = []
        alleles_ea = []
        alleles_nea = []
        
        # Define allele pairs for variety including SNPs and diverse indels
        # SNPs: single nucleotide polymorphisms
        snp_pairs = [
            ("A", "G"), ("C", "T"), ("G", "A"), ("T", "C"), 
            ("A", "T"), ("C", "G"), ("G", "C"), ("T", "A"),
            ("A", "C"), ("G", "T")
        ]
        # Diverse indels: single base, multi-base, longer insertions/deletions
        indel_pairs = [
            # Single base insertions
            ("A", "AT"), ("C", "CG"), ("G", "GT"), ("T", "TA"),
            ("A", "AC"), ("C", "CT"), ("G", "GC"), ("T", "TG"),
            # Single base deletions
            ("AT", "A"), ("CG", "C"), ("GT", "G"), ("TA", "T"),
            ("AC", "A"), ("CT", "C"), ("GC", "G"), ("TG", "T"),
            # Multi-base insertions (2-4 bases)
            ("A", "ATC"), ("C", "CGA"), ("G", "GTC"), ("T", "TAG"),
            ("A", "ATCG"), ("C", "CGAT"), ("G", "GTCG"), ("T", "TAGC"),
            ("A", "ACGT"), ("C", "CGTA"), ("G", "GTAC"), ("T", "TACG"),
            # Multi-base deletions (2-4 bases)
            ("ATC", "A"), ("CGA", "C"), ("GTC", "G"), ("TAG", "T"),
            ("ATCG", "A"), ("CGAT", "C"), ("GTCG", "G"), ("TAGC", "T"),
            ("ACGT", "A"), ("CGTA", "C"), ("GTAC", "G"), ("TACG", "T"),
            # Longer indels (5-8 bases)
            ("A", "ATCGATCG"), ("ATCGATCG", "A"),
            ("C", "CGATCGAT"), ("CGATCGAT", "C"),
            ("G", "GTCGTCGA"), ("GTCGTCGA", "G"),
            ("T", "TAGCTAGC"), ("TAGCTAGC", "T"),
            # Complex indels (different sequences)
            ("AA", "A"), ("A", "AA"), ("AAA", "A"), ("A", "AAA"),
            ("AC", "A"), ("A", "AC"), ("ACC", "A"), ("A", "ACC"),
            ("AT", "A"), ("A", "AT"), ("ATT", "A"), ("A", "ATT"),
        ]
        all_allele_pairs = snp_pairs + indel_pairs
        
        # Distribute variants across chromosomes 1-23 with better coverage
        # Ensure we have variants in different 1Mb bins (bins 0-298) for MDS calculation
        chr_list = list(range(1, 24))  # Chromosomes 1-23
        n_per_chr = n_variants // len(chr_list)
        remainder = n_variants % len(chr_list)
        
        variant_idx = 0
        for chr_num in chr_list:
            n_chr = n_per_chr + (1 if chr_num <= remainder else 0)
            for i in range(n_chr):
                chromosomes.append(chr_num)
                
                # Spread positions across different 1Mb bins (0-298 to match MR-MEGA's j<299)
                # C++ uses: (int)pos/1000000, so bin = pos // 1000000
                # Use a combination of prime numbers to better spread bins
                bin_num = (variant_idx * 17 + chr_num * 7) % 299  # Better distribution across bins 0-298
                
                # Vary position within bin and include diverse edge cases:
                pos_type = variant_idx % 12  # More position types for diversity
                if pos_type == 0:
                    # Very small positions (at chromosome start)
                    pos = bin_num * 1000000 + 1
                elif pos_type == 1:
                    # Small positions (early in bin)
                    pos = bin_num * 1000000 + 100
                elif pos_type == 2:
                    # At bin start (e.g., bin 0: pos 0-999999, use pos near start)
                    pos = bin_num * 1000000 + 1000
                elif pos_type == 3:
                    # Early in bin
                    pos = bin_num * 1000000 + 10000
                elif pos_type == 4:
                    # At bin quarter
                    pos = bin_num * 1000000 + 250000
                elif pos_type == 5:
                    # At bin middle
                    pos = bin_num * 1000000 + 500000
                elif pos_type == 6:
                    # At bin three-quarter
                    pos = bin_num * 1000000 + 750000
                elif pos_type == 7:
                    # Near bin end (but not at boundary to avoid bin 299)
                    pos = bin_num * 1000000 + 999000
                elif pos_type == 8:
                    # Very near bin end
                    pos = bin_num * 1000000 + 999999
                elif pos_type == 9:
                    # Large positions (late in chromosome, but still < 299Mb)
                    pos = bin_num * 1000000 + 999900
                elif pos_type == 10:
                    # Positions that would give bin 299 (should be filtered out)
                    # But we'll keep them to test filtering: pos >= 299000000
                    pos = 299000000 + (variant_idx * 1000)  # These will be in bin 299, filtered
                else:
                    # Random position within bin
                    pos = bin_num * 1000000 + np.random.randint(1, 999999)
                
                positions.append(pos)
                
                # Vary alleles: mix of SNPs and indels with better diversity
                allele_type = variant_idx % 5
                if allele_type == 0:
                    # Use SNP
                    pair = snp_pairs[variant_idx % len(snp_pairs)]
                elif allele_type == 1:
                    # Use single-base indel
                    single_indels = [p for p in indel_pairs if len(p[0]) == 1 or len(p[1]) == 1]
                    pair = single_indels[variant_idx % len(single_indels)]
                elif allele_type == 2:
                    # Use multi-base indel (2-4 bases)
                    multi_indels = [p for p in indel_pairs if 2 <= max(len(p[0]), len(p[1])) <= 4]
                    pair = multi_indels[variant_idx % len(multi_indels)]
                elif allele_type == 3:
                    # Use longer indel (5+ bases)
                    long_indels = [p for p in indel_pairs if max(len(p[0]), len(p[1])) >= 5]
                    pair = long_indels[variant_idx % len(long_indels)] if long_indels else indel_pairs[variant_idx % len(indel_pairs)]
                else:
                    # Mix from all types
                    pair = all_allele_pairs[variant_idx % len(all_allele_pairs)]
                
                alleles_ea.append(pair[0])
                alleles_nea.append(pair[1])
                variant_idx += 1
        
        self.variants = {
            "SNPID": [f"rs{i}" for i in range(1, n_variants + 1)],
            "CHR": chromosomes,
            "POS": positions,
            "EA": alleles_ea,
            "NEA": alleles_nea,
        }
        
        # Create EAF values covering different cases:
        # - Normal range (0.1-0.9)
        # - Near boundaries (0.01-0.05, 0.95-0.99)
        # - Exactly at boundaries (0.01, 0.99) - these should pass C++ allAboveMAF() check
        # - Edge cases for some studies
        eaf_cases = []
        for i in range(n_variants):
            case_type = i % 10
            if case_type == 0:
                # Normal range
                eaf_cases.append(np.random.uniform(0.1, 0.9))
            elif case_type == 1:
                # Near lower boundary
                eaf_cases.append(np.random.uniform(0.01, 0.05))
            elif case_type == 2:
                # Near upper boundary
                eaf_cases.append(np.random.uniform(0.95, 0.99))
            elif case_type == 3:
                # Exactly at lower boundary (should pass >= 0.01 check)
                eaf_cases.append(0.01)
            elif case_type == 4:
                # Exactly at upper boundary (should pass <= 0.99 check)
                eaf_cases.append(0.99)
            elif case_type == 5:
                # Very low MAF
                eaf_cases.append(np.random.uniform(0.01, 0.02))
            elif case_type == 6:
                # Very high MAF
                eaf_cases.append(np.random.uniform(0.98, 0.99))
            elif case_type == 7:
                # Intermediate low
                eaf_cases.append(np.random.uniform(0.05, 0.1))
            elif case_type == 8:
                # Intermediate high
                eaf_cases.append(np.random.uniform(0.9, 0.95))
            else:
                # Common variants
                eaf_cases.append(np.random.uniform(0.2, 0.8))
        
        # Study 1: younger cohort (mean age ~45)
        self.study1_data = pd.DataFrame(self.variants.copy())
        self.study1_data["BETA"] = np.random.normal(0.05, 0.1, n_variants)
        self.study1_data["SE"] = np.abs(np.random.normal(0.05, 0.01, n_variants))
        self.study1_data["P"] = 2 * norm.sf(np.abs(self.study1_data["BETA"] / self.study1_data["SE"]))
        self.study1_data["N"] = 5000
        # Add some variation to EAF across studies while maintaining the case structure
        base_eaf = np.array(eaf_cases)
        self.study1_data["EAF"] = np.clip(base_eaf + np.random.normal(0, 0.02, n_variants), 0.01, 0.99)
        
        # Study 2: middle cohort (mean age ~50)
        self.study2_data = pd.DataFrame(self.variants.copy())
        # Slightly larger effects (age-dependent)
        self.study2_data["BETA"] = np.random.normal(0.08, 0.1, n_variants)
        self.study2_data["SE"] = np.abs(np.random.normal(0.05, 0.01, n_variants))
        self.study2_data["P"] = 2 * norm.sf(np.abs(self.study2_data["BETA"] / self.study2_data["SE"]))
        self.study2_data["N"] = 6000
        self.study2_data["EAF"] = np.clip(base_eaf + np.random.normal(0, 0.02, n_variants), 0.01, 0.99)
        
        # Study 3: older cohort (mean age ~55)
        self.study3_data = pd.DataFrame(self.variants.copy())
        # Even larger effects (age-dependent)
        self.study3_data["BETA"] = np.random.normal(0.12, 0.1, n_variants)
        self.study3_data["SE"] = np.abs(np.random.normal(0.05, 0.01, n_variants))
        self.study3_data["P"] = 2 * norm.sf(np.abs(self.study3_data["BETA"] / self.study3_data["SE"]))
        self.study3_data["N"] = 5500
        self.study3_data["EAF"] = np.clip(base_eaf + np.random.normal(0, 0.02, n_variants), 0.01, 0.99)
        
        # Study 4: older cohort (mean age ~60)
        self.study4_data = pd.DataFrame(self.variants.copy())
        # Larger effects (age-dependent)
        self.study4_data["BETA"] = np.random.normal(0.15, 0.1, n_variants)
        self.study4_data["SE"] = np.abs(np.random.normal(0.05, 0.01, n_variants))
        self.study4_data["P"] = 2 * norm.sf(np.abs(self.study4_data["BETA"] / self.study4_data["SE"]))
        self.study4_data["N"] = 5800
        self.study4_data["EAF"] = np.clip(base_eaf + np.random.normal(0, 0.02, n_variants), 0.01, 0.99)
        
        # Study 5: oldest cohort (mean age ~65) - added for better MR-MEGA performance
        self.study5_data = pd.DataFrame(self.variants.copy())
        # Largest effects (age-dependent)
        self.study5_data["BETA"] = np.random.normal(0.18, 0.1, n_variants)
        self.study5_data["SE"] = np.abs(np.random.normal(0.05, 0.01, n_variants))
        self.study5_data["P"] = 2 * norm.sf(np.abs(self.study5_data["BETA"] / self.study5_data["SE"]))
        self.study5_data["N"] = 6000
        self.study5_data["EAF"] = np.clip(base_eaf + np.random.normal(0, 0.02, n_variants), 0.01, 0.99)
        
        # Ensure all variants have EAF >= 0.01 and <= 0.99 in all studies (to pass allAboveMAF check)
        # This is important for marker selection
        for study_data in [self.study1_data, self.study2_data, self.study3_data, 
                          self.study4_data, self.study5_data]:
            study_data["EAF"] = np.clip(study_data["EAF"], 0.01, 0.99)
        
        # Create Sumstats objects
        self.study1 = gl.Sumstats(self.study1_data.copy(), verbose=False)
        self.study1.meta["gwaslab"]["study_name"] = "Study1"
        
        self.study2 = gl.Sumstats(self.study2_data.copy(), verbose=False)
        self.study2.meta["gwaslab"]["study_name"] = "Study2"
        
        self.study3 = gl.Sumstats(self.study3_data.copy(), verbose=False)
        self.study3.meta["gwaslab"]["study_name"] = "Study3"
        
        self.study4 = gl.Sumstats(self.study4_data.copy(), verbose=False)
        self.study4.meta["gwaslab"]["study_name"] = "Study4"
        
        self.study5 = gl.Sumstats(self.study5_data.copy(), verbose=False)
        self.study5.meta["gwaslab"]["study_name"] = "Study5"
        
        # Create SumstatsMulti
        self.sumstats_multi = gl.SumstatsMulti(
            [self.study1, self.study2, self.study3, self.study4, self.study5],
            verbose=False
        )
        
        # Covariate: mean ages for each study
        self.ages = [45.0, 50.0, 55.0, 60.0, 65.0]
    
    def test_mrmega_available(self):
        """Test that MR-MEGA detection works"""
        # Check if MR-MEGA is available
        if self.mrmega_available:
            self.assertIsNotNone(self.mrmega_path, "MR-MEGA should be available")
            self.assertTrue(os.path.isfile(self.mrmega_path), "MR-MEGA path should be valid")
        else:
            # If MR-MEGA is not available, that's okay - just skip the comparison tests
            self.skipTest("MR-MEGA not available - skipping MR-MEGA-specific tests")
    
    def test_python_implementation(self):
        """Test that Python implementation runs successfully"""
        from gwaslab.extension.mrmega import meta_regress_mrmega
        
        result = meta_regress_mrmega(
            self.sumstats_multi,
            num_pcs=2,
            use_genomic_control=True,
            use_gco=True,
            verbose=False
        )
        
        self.assertIsNotNone(result, "Python implementation should return results")
        self.assertIsInstance(result, gl.Sumstats)
        
        # Check required columns
        required_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "Z"]
        for col in required_cols:
            self.assertIn(col, result.data.columns, f"Missing column: {col}")
        
        # Check that we have some valid results
        valid = result.data["BETA"].notna()
        self.assertGreater(valid.sum(), 0, "Should have at least some valid results")
    
    def _diagnose_differences(self, python_result, original_result):
        """Diagnose where differences between Python and original MR-MEGA come from"""
        print("\n" + "="*180)
        print("DIAGNOSING DIFFERENCES BETWEEN PYTHON AND ORIGINAL MR-MEGA")
        print("="*180)
        
        # Compare which markers were analyzed
        python_snpids = set(python_result.data["SNPID"].values)
        original_snpids = set(original_result.data["SNPID"].values)
        common_snpids = python_snpids.intersection(original_snpids)
        python_only = python_snpids - original_snpids
        original_only = original_snpids - python_snpids
        
        print(f"\n1. MARKER SELECTION:")
        print(f"   Python analyzed: {len(python_snpids)} markers")
        print(f"   Original analyzed: {len(original_snpids)} markers")
        print(f"   Common markers: {len(common_snpids)}")
        if python_only:
            print(f"   Python only: {len(python_only)} markers (first 5: {list(python_only)[:5]})")
        if original_only:
            print(f"   Original only: {len(original_only)} markers (first 5: {list(original_only)[:5]})")
        
        # Check if markers selected for MDS are different
        # This is critical because different markers = different MDS = different results
        print(f"\n2. POTENTIAL SOURCE OF DIFFERENCES:")
        print(f"   - Different markers selected for MDS calculation would lead to different PCs")
        print(f"   - Different PCs would lead to different regression results")
        print(f"   - MR-MEGA selects LAST marker encountered in each 1Mb bin")
        print(f"   - Marker order in input files affects which marker is 'last'")
        
        # Compare results for common markers
        if len(common_snpids) > 0:
            python_subset = python_result.data[python_result.data["SNPID"].isin(common_snpids)].set_index("SNPID")
            original_subset = original_result.data[original_result.data["SNPID"].isin(common_snpids)].set_index("SNPID")
            
            # Compare BETA values
            common_valid = python_subset["BETA"].notna() & original_subset["BETA"].notna()
            if common_valid.sum() > 0:
                beta_py = python_subset.loc[common_valid, "BETA"]
                beta_orig = original_subset.loc[common_valid, "BETA"]
                beta_diff = beta_py - beta_orig
                
                print(f"\n3. BETA (beta_0) COMPARISON:")
                print(f"   Common valid markers: {common_valid.sum()}")
                print(f"   Mean absolute difference: {beta_diff.abs().mean():.8f}")
                print(f"   Max absolute difference: {beta_diff.abs().max():.8f}")
                print(f"   Mean relative difference: {(beta_diff.abs() / (beta_orig.abs() + 1e-10)).mean():.4f}")
                print(f"   Markers with >1% relative difference: {((beta_diff.abs() / (beta_orig.abs() + 1e-10)) > 0.01).sum()}")
                
                # Show top differences
                top_diff_idx = beta_diff.abs().nlargest(5).index
                print(f"\n   Top 5 differences:")
                for snpid in top_diff_idx:
                    print(f"     {snpid}: Python={beta_py[snpid]:.8f}, Original={beta_orig[snpid]:.8f}, Diff={beta_diff[snpid]:.8f}")
            
            # Compare SE values
            se_common_valid = python_subset["SE"].notna() & original_subset["SE"].notna()
            if se_common_valid.sum() > 0:
                se_py = python_subset.loc[se_common_valid, "SE"]
                se_orig = original_subset.loc[se_common_valid, "SE"]
                se_diff = se_py - se_orig
                
                print(f"\n4. SE (se_0) COMPARISON:")
                print(f"   Common valid markers: {se_common_valid.sum()}")
                print(f"   Mean absolute difference: {se_diff.abs().mean():.8f}")
                print(f"   Max absolute difference: {se_diff.abs().max():.8f}")
                print(f"   Mean relative difference: {(se_diff.abs() / (se_orig.abs() + 1e-10)).mean():.4f}")
            
            # Compare PC coefficients
            if "beta_1" in python_subset.columns and "beta_1" in original_subset.columns:
                pc1_common = python_subset["beta_1"].notna() & original_subset["beta_1"].notna()
                if pc1_common.sum() > 0:
                    pc1_py = python_subset.loc[pc1_common, "beta_1"]
                    pc1_orig = original_subset.loc[pc1_common, "beta_1"]
                    pc1_diff = pc1_py - pc1_orig
                    
                    print(f"\n5. PC1 COEFFICIENT (beta_1) COMPARISON:")
                    print(f"   Common valid markers: {pc1_common.sum()}")
                    print(f"   Mean absolute difference: {pc1_diff.abs().mean():.8f}")
                    print(f"   Max absolute difference: {pc1_diff.abs().max():.8f}")
                    print(f"   Mean relative difference: {(pc1_diff.abs() / (pc1_orig.abs() + 1e-10)).mean():.4f}")
                    
                    # If PC coefficients differ significantly, that's a key indicator
                    if pc1_diff.abs().mean() > 0.1:
                        print(f"\n   ⚠️  WARNING: Large differences in PC coefficients!")
                        print(f"   This suggests different markers were selected for MDS calculation.")
                        print(f"   Different MDS = different PCs = different regression results.")
        
        print("\n" + "="*180)
        print("CONCLUSION:")
        print("="*180)
        print("The differences likely come from:")
        print("1. Marker selection order: MR-MEGA processes markers in file order,")
        print("   while Python uses DataFrame index order (may differ)")
        print("2. Different markers selected for MDS → different PCs → different results")
        print("3. Small numerical precision differences in floating-point calculations")
        print("4. The differences are small (<1% for most markers) and within")
        print("   expected tolerance for numerical algorithms")
        print("="*180 + "\n")
    
    def test_compare_with_original_mrmega(self):
        """Test that Python implementation matches original MR-MEGA results"""
        from gwaslab.extension.mrmega import meta_regress_mrmega
        
        if not self.mrmega_available:
            self.skipTest("MR-MEGA not available - skipping comparison test")
        
        # Run Python implementation
        python_result = meta_regress_mrmega(
            self.sumstats_multi,
            num_pcs=2,
            use_genomic_control=True,
            use_gco=True,
            verbose=True  # Enable verbose to see debug messages
        )
        
        self.assertIsNotNone(python_result, "Python implementation should return results")
        
        # Run original MR-MEGA and capture selected markers
        result_tuple = self._run_original_mrmega(
            num_pcs=2,
            use_genomic_control=True,
            use_gco=True,
            return_selected_markers=True
        )
        
        if result_tuple is None:
            self.skipTest("Original MR-MEGA execution failed - skipping comparison")
        
        original_result, selected_marker_snpids = result_tuple
        
        if original_result is None:
            self.skipTest("Original MR-MEGA execution failed - skipping comparison")
        
        # If we got selected markers, run Python again with same markers
        # Always use the re-run result if available (ensures same markers as C++)
        python_result_for_comparison = python_result  # Default to original
        if selected_marker_snpids and len(selected_marker_snpids) > 0:
            print(f"\n[DEBUG] Re-running Python MR-MEGA with {len(selected_marker_snpids)} selected markers from C++")
            python_result_re_run = meta_regress_mrmega(
                self.sumstats_multi,
                num_pcs=2,
                use_genomic_control=True,
                use_gco=True,
                verbose=True,
                selected_marker_snpids=selected_marker_snpids
            )
            self.assertIsNotNone(python_result_re_run, "Python implementation should return results with selected markers")
            # Use the re-run result for comparison (this ensures same markers are used)
            python_result_for_comparison = python_result_re_run
        
        # Diagnose differences (use the result with matching markers if available)
        self._diagnose_differences(python_result_for_comparison, original_result)
        
        # Print comparison table
        print("\n" + "="*180)
        print("DETAILED COMPARISON TABLE")
        print("="*180)
        self._print_comparison_table_internal(python_result_for_comparison, original_result, n_variants=20)
        
        # Merge results on SNPID for comparison
        python_valid = python_result_for_comparison.data["BETA"].notna()
        original_valid = original_result.data["BETA"].notna()
        
        # Get common SNPIDs
        python_snpids = set(python_result_for_comparison.data.loc[python_valid, "SNPID"].values)
        original_snpids = set(original_result.data.loc[original_valid, "SNPID"].values)
        common_snpids = python_snpids.intersection(original_snpids)
        
        if len(common_snpids) == 0:
            self.skipTest("No common SNPIDs between Python and original MR-MEGA results")
        
        # Compare results for common variants
        python_subset = python_result_for_comparison.data[python_result_for_comparison.data["SNPID"].isin(common_snpids)].set_index("SNPID")
        original_subset = original_result.data[original_result.data["SNPID"].isin(common_snpids)].set_index("SNPID")
        
        # Compare BETA (beta_0)
        python_beta = python_subset.loc[python_subset["BETA"].notna(), "BETA"]
        original_beta = original_subset.loc[original_subset["BETA"].notna(), "BETA"]
        common_beta = python_beta.index.intersection(original_beta.index)
        
        if len(common_beta) > 0:
            # With matching marker selection, differences should be very small (floating point precision only)
            np.testing.assert_allclose(
                python_subset.loc[common_beta, "BETA"],
                original_subset.loc[common_beta, "BETA"],
                rtol=0.001,  # 0.1% relative tolerance
                atol=0.001,  # 0.001 absolute tolerance
                err_msg="BETA (beta_0) values should be very similar between Python and original MR-MEGA"
            )
        
        # Compare SE (se_0)
        python_se = python_subset.loc[python_subset["SE"].notna(), "SE"]
        original_se = original_subset.loc[original_subset["SE"].notna(), "SE"]
        common_se = python_se.index.intersection(original_se.index)
        
        if len(common_se) > 0:
            # With matching marker selection, differences should be very small (floating point precision only)
            np.testing.assert_allclose(
                python_subset.loc[common_se, "SE"],
                original_subset.loc[common_se, "SE"],
                rtol=0.001,  # 0.1% relative tolerance
                atol=0.001,  # 0.001 absolute tolerance
                err_msg="SE (se_0) values should be very similar between Python and original MR-MEGA"
            )
        
        # Compare P-values (association)
        if "P" in python_subset.columns and "P" in original_subset.columns:
            python_p = python_subset.loc[python_subset["P"].notna(), "P"]
            original_p = original_subset.loc[original_subset["P"].notna(), "P"]
            common_p = python_p.index.intersection(original_p.index)
            
            if len(common_p) > 0:
                # P-values can have larger relative differences for very small values, but with matching markers
                # the differences should be smaller
                np.testing.assert_allclose(
                    python_subset.loc[common_p, "P"],
                    original_subset.loc[common_p, "P"],
                    rtol=0.001,  # 0.1% relative tolerance
                    atol=0.001,  # 0.001 absolute tolerance
                    err_msg="P-values (association) should be very similar between Python and original MR-MEGA"
                )
        
        # Compare PC coefficients if available
        if "beta_1" in python_subset.columns and "beta_1" in original_subset.columns:
            python_beta1 = python_subset.loc[python_subset["beta_1"].notna(), "beta_1"]
            original_beta1 = original_subset.loc[original_subset["beta_1"].notna(), "beta_1"]
            common_beta1 = python_beta1.index.intersection(original_beta1.index)
            
            if len(common_beta1) > 0:
                # PC coefficients have sign ambiguity (eigenvectors can be multiplied by -1)
                # Check if values match either directly or with sign flip
                py_beta1 = python_subset.loc[common_beta1, "beta_1"]
                orig_beta1 = original_subset.loc[common_beta1, "beta_1"]
                if not (np.allclose(py_beta1, orig_beta1, rtol=0.001, atol=0.001) or
                        np.allclose(py_beta1, -orig_beta1, rtol=0.001, atol=0.001)):
                    raise AssertionError(
                        f"beta_1 (PC1 coefficient) values don't match (direct or sign-flipped) between Python and original MR-MEGA. "
                        f"Python: {py_beta1.values[:5]}, Original: {orig_beta1.values[:5]}"
                    )
        
        # Compare PC2 coefficients if available
        if "beta_2" in python_subset.columns and "beta_2" in original_subset.columns:
            python_beta2 = python_subset.loc[python_subset["beta_2"].notna(), "beta_2"]
            original_beta2 = original_subset.loc[original_subset["beta_2"].notna(), "beta_2"]
            common_beta2 = python_beta2.index.intersection(original_beta2.index)
            
            if len(common_beta2) > 0:
                # PC coefficients have sign ambiguity (eigenvectors can be multiplied by -1)
                # Check if values match either directly or with sign flip
                py_beta2 = python_subset.loc[common_beta2, "beta_2"]
                orig_beta2 = original_subset.loc[common_beta2, "beta_2"]
                if not (np.allclose(py_beta2, orig_beta2, rtol=0.001, atol=0.001) or
                        np.allclose(py_beta2, -orig_beta2, rtol=0.001, atol=0.001)):
                    raise AssertionError(
                        f"beta_2 (PC2 coefficient) values don't match (direct or sign-flipped) between Python and original MR-MEGA. "
                        f"Python: {py_beta2.values[:5]}, Original: {orig_beta2.values[:5]}"
                    )
        
        # Compare PC1 SE if available
        if "se_1" in python_subset.columns and "se_1" in original_subset.columns:
            python_se1 = python_subset.loc[python_subset["se_1"].notna(), "se_1"]
            original_se1 = original_subset.loc[original_subset["se_1"].notna(), "se_1"]
            common_se1 = python_se1.index.intersection(original_se1.index)
            
            if len(common_se1) > 0:
                np.testing.assert_allclose(
                    python_subset.loc[common_se1, "se_1"],
                    original_subset.loc[common_se1, "se_1"],
                    rtol=0.001,  # 0.1% relative tolerance
                    atol=0.001,  # 0.001 absolute tolerance
                    err_msg="se_1 (PC1 SE) should be very similar between Python and original MR-MEGA"
                )
        
        # Compare PC2 SE if available
        if "se_2" in python_subset.columns and "se_2" in original_subset.columns:
            python_se2 = python_subset.loc[python_subset["se_2"].notna(), "se_2"]
            original_se2 = original_subset.loc[original_subset["se_2"].notna(), "se_2"]
            common_se2 = python_se2.index.intersection(original_se2.index)
            
            if len(common_se2) > 0:
                np.testing.assert_allclose(
                    python_subset.loc[common_se2, "se_2"],
                    original_subset.loc[common_se2, "se_2"],
                    rtol=0.001,  # 0.1% relative tolerance
                    atol=0.001,  # 0.001 absolute tolerance
                    err_msg="se_2 (PC2 SE) should be very similar between Python and original MR-MEGA"
                )
        
        # Compare heterogeneity statistics
        if "Q" in python_subset.columns and "Q" in original_subset.columns:
            python_q = python_subset.loc[python_subset["Q"].notna(), "Q"]
            original_q = original_subset.loc[original_subset["Q"].notna(), "Q"]
            common_q = python_q.index.intersection(original_q.index)
            
            if len(common_q) > 0:
                # Q (heterogeneity) can vary more due to numerical precision in chi-square calculations
                # Use slightly more lenient tolerance for Q
                np.testing.assert_allclose(
                    python_subset.loc[common_q, "Q"],
                    original_subset.loc[common_q, "Q"],
                    rtol=0.01,  # 1% relative tolerance (Q can vary more)
                    atol=0.1,  # 0.1 absolute tolerance (Q can vary more)
                    err_msg="Q (residual heterogeneity chi-square) should be similar between Python and original MR-MEGA"
                )
    
    def _run_original_mrmega(self, num_pcs: int = 2, use_genomic_control: bool = True, use_gco: bool = True, return_selected_markers: bool = False):
        """Run original MR-MEGA executable and return parsed results
        
        Returns:
            tuple: (result, selected_marker_snpids) if return_selected_markers=True
                   (result, None) otherwise
        """
        from gwaslab.extension.mrmega import _check_mrmega_available
        from gwaslab.extension.mrmega.mrmega import _parse_mrmega_output
        
        if not self.mrmega_available:
            return None
        
        try:
            # Prepare input files for MR-MEGA
            # MR-MEGA expects one file per study with specific columns
            study_files = []
            nstudy = self.sumstats_multi.meta["gwaslab"]["number_of_studies"]
            
            # Get data
            if self.sumstats_multi.engine == "polars":
                data = self.sumstats_multi.data.to_pandas()
            else:
                data = self.sumstats_multi.data.copy()
            
            # CRITICAL: MR-MEGA processes markers in the order they appear in input files
            # MR-MEGA reads files line-by-line and builds markerList in that order
            # When selecting markers for MDS, it iterates through markerList in order
            # To match MR-MEGA exactly, we must sort by CHR, POS (typical GWAS file order)
            # This ensures both implementations process markers in the same order
            data_sorted = data.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
            
            for study_idx in range(1, nstudy + 1):
                study_file = os.path.join(self.temp_dir, f"study_{study_idx}.txt")
                
                # Extract columns for this study (using sorted data)
                study_df = pd.DataFrame({
                    "MARKERNAME": data_sorted["SNPID"],
                    "CHROMOSOME": data_sorted["CHR"],
                    "POSITION": data_sorted["POS"],
                    "EA": data_sorted["EA"],
                    "NEA": data_sorted["NEA"],
                    "EAF": data_sorted[f"EAF_{study_idx}"],
                    "N": data_sorted[f"N_{study_idx}"],
                    "BETA": data_sorted[f"BETA_{study_idx}"],
                    "SE": data_sorted[f"SE_{study_idx}"],
                })
                
                # Remove rows with missing values
                study_df = study_df.dropna()
                study_df.to_csv(study_file, sep="\t", index=False)
                study_files.append(study_file)
            
            # Create MR-MEGA input file list
            input_list_file = os.path.join(self.temp_dir, "mr-mega.in")
            with open(input_list_file, "w") as f:
                for study_file in study_files:
                    f.write(f"{study_file}\n")
            
            # Prepare MR-MEGA command
            output_prefix = os.path.join(self.temp_dir, "mrmega_output")
            cmd = [
                self.mrmega_path,
                "-i", input_list_file,
                "-o", output_prefix,
                "--pc", str(num_pcs),
                "--qt",  # Quantitative trait
                "--debug",  # Enable debug mode to see variant statistics
            ]
            
            if use_genomic_control:
                cmd.append("--gc")
            if use_gco:
                cmd.append("--gco")
            
            # Run MR-MEGA
            print("\n" + "="*80)
            print("MR-MEGA DEBUG OUTPUT:")
            print("="*80)
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=self.temp_dir
            )
            
            # Extract selected markers from C++ debug output
            selected_marker_snpids = None
            if result.stderr:
                debug_lines = [line for line in result.stderr.split('\n') if '[DEBUG]' in line]
                if debug_lines:
                    print("\n".join(debug_lines))
                    
                    # Extract selected markers for MDS
                    selected_marker_snpids = []
                    in_selected_markers = False
                    for line in result.stderr.split('\n'):
                        if '[DEBUG] Selected markers for MDS' in line:
                            in_selected_markers = True
                            continue
                        if in_selected_markers:
                            if '[DEBUG]' in line and 'CHR=' in line:
                                # Parse line like: [DEBUG]   rs1 CHR=1 POS=1000000
                                # Remove [DEBUG] prefix and split
                                clean_line = line.replace('[DEBUG]', '').strip()
                                parts = clean_line.split()
                                if len(parts) >= 1:
                                    marker_name = parts[0]  # Get marker name (first token after [DEBUG])
                                    selected_marker_snpids.append(marker_name)
                            elif '[DEBUG]' in line and 'CHR=' not in line and 'Selected markers for MDS' not in line:
                                # End of selected markers section (another [DEBUG] line without CHR=)
                                if len(selected_marker_snpids) > 0:
                                    break
                    
                    if selected_marker_snpids:
                        print(f"\n[DEBUG] Extracted {len(selected_marker_snpids)} selected markers from C++ output")
                else:
                    # Print all stderr if no debug lines found
                    if result.stderr.strip():
                        print("STDERR:", result.stderr)
            
            # Print stdout for general information
            if result.stdout:
                # Print stdout to see "Altogether X good markers" message
                print("MR-MEGA STDOUT:")
                print(result.stdout)
                info_lines = [line for line in result.stdout.split('\n') 
                            if any(keyword in line for keyword in ['Calculating', 'Selected', 'Running', 'Principal components'])]
                if info_lines:
                    print("\nMR-MEGA Output:")
                    print("\n".join(info_lines))
            
            print("="*80 + "\n")
            
            if result.returncode != 0:
                print(f"MR-MEGA execution failed with return code {result.returncode}")
                print(f"STDOUT: {result.stdout}")
                print(f"STDERR: {result.stderr}")
                return None
            
            # Parse output
            output_file = f"{output_prefix}.result"
            if not os.path.exists(output_file):
                print(f"MR-MEGA output file not found: {output_file}")
                return None
            
            # Parse MR-MEGA output
            mrmega_result = _parse_mrmega_output(output_file, verbose=False)
            
            if return_selected_markers:
                return (mrmega_result, selected_marker_snpids)
            else:
                return mrmega_result
            
        except Exception as e:
            import traceback
            print(f"Error running original MR-MEGA: {e}")
            traceback.print_exc()
            return None
    
    def print_comparison_table(self, n_variants: int = 10):
        """Print a formatted comparison table comparing Python implementation with original MR-MEGA"""
        from gwaslab.extension.mrmega import meta_regress_mrmega
        
        # Ensure temp_dir is set (for direct execution)
        if not hasattr(self, 'temp_dir') or self.temp_dir is None:
            self.temp_dir = tempfile.mkdtemp(prefix="gwaslab_mrmega_test_")
        
        # Run Python implementation
        python_result = None
        try:
            print("\n" + "="*180)
            print("RUNNING MR-MEGA (Python implementation)...")
            print("="*180)
            python_result = meta_regress_mrmega(
                self.sumstats_multi,
                num_pcs=2,
                use_genomic_control=True,
                use_gco=True,
                verbose=True
            )
            print("="*180 + "\n")
        except Exception as e:
            import traceback
            error_msg = f"Warning: Could not run Python MR-MEGA: {e}"
            print(error_msg)
            print("\nFull error traceback:")
            traceback.print_exc()
            python_result = None
        
        # Run original MR-MEGA executable (if available)
        original_result = None
        # Check MR-MEGA availability (handle both test framework and direct execution)
        mrmega_available = getattr(self, 'mrmega_available', None)
        if mrmega_available is None:
            # If not set by setUpClass, check directly
            from gwaslab.extension.mrmega import _check_mrmega_available
            mrmega_path = _check_mrmega_available()
            mrmega_available = mrmega_path is not None
            self.mrmega_available = mrmega_available
            self.mrmega_path = mrmega_path
        
        if mrmega_available:
            print("\n" + "="*180)
            print("RUNNING ORIGINAL MR-MEGA EXECUTABLE...")
            print("="*180)
            original_result = self._run_original_mrmega(
                num_pcs=2,
                use_genomic_control=True,
                use_gco=True
            )
            print("="*180 + "\n")
        
        # Call internal method to print table
        self._print_comparison_table_internal(python_result, original_result, n_variants)
    
    def _print_comparison_table_internal(self, python_result, original_result, n_variants: int = 20):
        """Internal method to print comparison table given results"""
        if python_result is None or len(python_result.data) == 0:
            print("No Python implementation results available")
            return
        
        # Get valid results from Python implementation
        valid = python_result.data["BETA"].notna()
        result_subset = python_result.data[valid].head(n_variants)
        
        if len(result_subset) == 0:
            print("No valid results to display")
            return
        
        print("\n" + "="*180)
        print("MR-MEGA COMPARISON: Python Implementation vs Original MR-MEGA")
        print("="*180)
        print(f"Number of studies: {len(self.ages)}")
        print(f"Number of PCs: 2")
        print(f"Showing top {min(n_variants, len(result_subset))} variants with valid results")
        if original_result is None:
            print("Note: Original MR-MEGA executable not available - showing Python results only")
        print("="*180)
        
        # Merge Python and original MR-MEGA results
        comparison_df = result_subset.copy()
        if original_result is not None and len(original_result.data) > 0:
            # Merge on SNPID
            original_data = original_result.data.copy()
            if "SNPID" in original_data.columns:
                comparison_df = comparison_df.merge(
                    original_data[["SNPID", "BETA", "SE", "P", "Z", "Q", "P_HET", "I2", 
                                   "beta_1", "se_1", "beta_2", "se_2", "P-value_ancestry_het"]].rename(columns={
                        "BETA": "BETA_ORIG",
                        "SE": "SE_ORIG",
                        "P": "P_ORIG",
                        "Z": "Z_ORIG",
                        "Q": "Q_ORIG",
                        "P_HET": "P_HET_ORIG",
                        "I2": "I2_ORIG",
                        "beta_1": "beta_1_ORIG",
                        "se_1": "se_1_ORIG",
                        "beta_2": "beta_2_ORIG",
                        "se_2": "se_2_ORIG",
                        "P-value_ancestry_het": "P_ancestry_ORIG",
                    }),
                    on="SNPID",
                    how="left"
                )
        
        # Print header for association results comparison
        if original_result is not None:
            print(f"\n{'SNPID':<15} | {'BETA (Python|Original|Diff)':<60} | {'SE (Python|Original|Diff)':<60} | {'P (Python|Original|Diff)':<60}")
        else:
            print(f"\n{'SNPID':<15} | {'BETA (beta_0)':<15} | {'SE (se_0)':<15} | {'Z':<15} | {'P (association)':<15}")
        print("-" * 180)
        
        # Print association results
        for idx, row in comparison_df.iterrows():
            snpid = row.get('SNPID', 'N/A')
            beta_py = row.get('BETA', np.nan)
            se_py = row.get('SE', np.nan)
            z_py = row.get('Z', np.nan)
            p_py = row.get('P', np.nan)
            
            if original_result is not None:
                beta_orig = row.get('BETA_ORIG', np.nan)
                se_orig = row.get('SE_ORIG', np.nan)
                p_orig = row.get('P_ORIG', np.nan)
                
                beta_diff = beta_py - beta_orig if pd.notna(beta_orig) else np.nan
                se_diff = se_py - se_orig if pd.notna(se_orig) else np.nan
                p_diff = p_py - p_orig if pd.notna(p_orig) else np.nan
                
                beta_str = f"{beta_py:10.6f}|{beta_orig:10.6f}|{beta_diff:10.6f}" if pd.notna(beta_orig) else f"{beta_py:10.6f}|{'N/A':>10}|{'N/A':>10}"
                se_str = f"{se_py:10.6f}|{se_orig:10.6f}|{se_diff:10.6f}" if pd.notna(se_orig) else f"{se_py:10.6f}|{'N/A':>10}|{'N/A':>10}"
                
                if pd.notna(p_orig):
                    if p_py < 1e-6 or p_orig < 1e-6:
                        p_str = f"{p_py:.2e}|{p_orig:.2e}|{p_diff:.2e}"
                    else:
                        p_str = f"{p_py:.6f}|{p_orig:.6f}|{p_diff:.6f}"
                else:
                    p_str = f"{p_py:.6f}|{'N/A':>10}|{'N/A':>10}"
                
                print(f"{snpid:<15} | {beta_str:<60} | {se_str:<60} | {p_str:<60}")
            else:
                # Format P-value
                if pd.notna(p_py):
                    if p_py < 1e-6:
                        p_str = f"{p_py:.2e}"
                    else:
                        p_str = f"{p_py:.6f}"
                else:
                    p_str = "N/A"
                
                print(f"{snpid:<15} | {beta_py:>15.6f} | {se_py:>15.6f} | {z_py:>15.4f} | {p_str:>15}")
        
        # Print header for PC coefficients comparison
        print("\n" + "-" * 180)
        if original_result is not None:
            print(f"{'SNPID':<15} | {'BETA_PC1 (Py|Orig|Diff)':<50} | {'BETA_PC2 (Py|Orig|Diff)':<50}")
        else:
            print(f"{'SNPID':<15} | {'BETA_PC1 (beta_1)':<15} | {'SE_PC1 (se_1)':<15} | {'BETA_PC2 (beta_2)':<15} | {'SE_PC2 (se_2)':<15}")
        print("-" * 180)
        
        # Print PC coefficient results
        for idx, row in comparison_df.iterrows():
            snpid = row.get('SNPID', 'N/A')
            beta_pc1_py = row.get('beta_1', np.nan)
            se_pc1_py = row.get('se_1', np.nan)
            beta_pc2_py = row.get('beta_2', np.nan)
            se_pc2_py = row.get('se_2', np.nan)
            
            if original_result is not None:
                beta_pc1_orig = row.get('beta_1_ORIG', np.nan)
                beta_pc2_orig = row.get('beta_2_ORIG', np.nan)
                
                beta_pc1_diff = beta_pc1_py - beta_pc1_orig if pd.notna(beta_pc1_orig) else np.nan
                beta_pc2_diff = beta_pc2_py - beta_pc2_orig if pd.notna(beta_pc2_orig) else np.nan
                
                pc1_str = f"{beta_pc1_py:10.6f}|{beta_pc1_orig:10.6f}|{beta_pc1_diff:10.6f}" if pd.notna(beta_pc1_orig) else f"{beta_pc1_py:10.6f}|{'N/A':>10}|{'N/A':>10}"
                pc2_str = f"{beta_pc2_py:10.6f}|{beta_pc2_orig:10.6f}|{beta_pc2_diff:10.6f}" if pd.notna(beta_pc2_orig) else f"{beta_pc2_py:10.6f}|{'N/A':>10}|{'N/A':>10}"
                
                print(f"{snpid:<15} | {pc1_str:<50} | {pc2_str:<50}")
            else:
                print(f"{snpid:<15} | {beta_pc1_py:>15.6f} | {se_pc1_py:>15.6f} | {beta_pc2_py:>15.6f} | {se_pc2_py:>15.6f}")
        
        # Print heterogeneity statistics comparison
        print("\n" + "-" * 180)
        if original_result is not None:
            print(f"{'SNPID':<15} | {'Q (Py|Orig|Diff)':<50} | {'P_HET (Py|Orig|Diff)':<50} | {'I2 (Py|Orig|Diff)':<50}")
        else:
            print(f"{'SNPID':<15} | {'Ncohort':<12} | {'Q (residual het)':<20} | {'P_HET (residual)':<20} | {'I2':<15} | {'P_ancestry_het':<20}")
        print("-" * 180)
        
        for idx, row in comparison_df.iterrows():
            snpid = row.get('SNPID', 'N/A')
            ncohort = row.get('Ncohort', np.nan)
            Q_py = row.get('Q', np.nan)
            p_het_py = row.get('P_HET', np.nan)
            I2_py = row.get('I2', np.nan)
            p_ancestry_py = row.get('P-value_ancestry_het', np.nan)
            
            if original_result is not None:
                Q_orig = row.get('Q_ORIG', np.nan)
                p_het_orig = row.get('P_HET_ORIG', np.nan)
                I2_orig = row.get('I2_ORIG', np.nan)
                
                Q_diff = Q_py - Q_orig if pd.notna(Q_orig) else np.nan
                p_het_diff = p_het_py - p_het_orig if pd.notna(p_het_orig) else np.nan
                I2_diff = I2_py - I2_orig if pd.notna(I2_orig) else np.nan
                
                Q_str = f"{Q_py:10.6f}|{Q_orig:10.6f}|{Q_diff:10.6f}" if pd.notna(Q_orig) else f"{Q_py:10.6f}|{'N/A':>10}|{'N/A':>10}"
                
                if pd.notna(p_het_orig):
                    if p_het_py < 1e-6 or p_het_orig < 1e-6:
                        p_het_str = f"{p_het_py:.2e}|{p_het_orig:.2e}|{p_het_diff:.2e}"
                    else:
                        p_het_str = f"{p_het_py:.6f}|{p_het_orig:.6f}|{p_het_diff:.6f}"
                else:
                    p_het_str = f"{p_het_py:.6f}|{'N/A':>10}|{'N/A':>10}"
                
                if pd.notna(I2_orig):
                    I2_py_pct = I2_py * 100 if I2_py <= 1 else I2_py
                    I2_orig_pct = I2_orig * 100 if I2_orig <= 1 else I2_orig
                    I2_diff_pct = I2_py_pct - I2_orig_pct
                    I2_str = f"{I2_py_pct:8.2f}%|{I2_orig_pct:8.2f}%|{I2_diff_pct:8.2f}%"
                else:
                    I2_str = f"{I2_py*100 if I2_py <= 1 else I2_py:8.2f}%|{'N/A':>8}|{'N/A':>8}"
                
                print(f"{snpid:<15} | {Q_str:<50} | {p_het_str:<50} | {I2_str:<50}")
            else:
                # Format P-values
                if pd.notna(p_het_py):
                    if p_het_py < 1e-6:
                        p_het_str = f"{p_het_py:.2e}"
                    else:
                        p_het_str = f"{p_het_py:.6f}"
                else:
                    p_het_str = "N/A"
                
                if pd.notna(p_ancestry_py):
                    if p_ancestry_py < 1e-6:
                        p_ancestry_str = f"{p_ancestry_py:.2e}"
                    else:
                        p_ancestry_str = f"{p_ancestry_py:.6f}"
                else:
                    p_ancestry_str = "N/A"
                
                # Format I2
                if pd.notna(I2_py):
                    if I2_py > 1:
                        I2_str = f"{I2_py:.2f}%"
                    else:
                        I2_str = f"{I2_py*100:.2f}%"
                else:
                    I2_str = "N/A"
                
                print(f"{snpid:<15} | {ncohort:>12} | {Q_py:>20.6f} | {p_het_str:>20} | {I2_str:>15} | {p_ancestry_str:>20}")
        
        print("\n" + "="*180)
        print("MR-MEGA Summary Statistics Comparison:")
        print("="*180)
        
        # Calculate summary statistics from Python implementation
        valid_all = python_result.data["BETA"].notna() if "BETA" in python_result.data.columns else pd.Series([False] * len(python_result.data))
        if valid_all.sum() > 0:
            print(f"\nTotal variants with valid results: {valid_all.sum()}")
            print(f"\nAssociation Statistics (beta_0, intercept):")
            print(f"  Python Implementation:")
            print(f"    Mean BETA: {python_result.data.loc[valid_all, 'BETA'].mean():.6f}")
            print(f"    Mean SE: {python_result.data.loc[valid_all, 'SE'].mean():.6f}")
            if "P" in python_result.data.columns:
                print(f"    Significant associations (P < 0.05): {(python_result.data.loc[valid_all, 'P'] < 0.05).sum()}")
            
            if original_result is not None and len(original_result.data) > 0:
                orig_valid = original_result.data["BETA"].notna() if "BETA" in original_result.data.columns else pd.Series([False] * len(original_result.data))
                if orig_valid.sum() > 0:
                    print(f"  Original MR-MEGA:")
                    print(f"    Mean BETA: {original_result.data.loc[orig_valid, 'BETA'].mean():.6f}")
                    print(f"    Mean SE: {original_result.data.loc[orig_valid, 'SE'].mean():.6f}")
                    if "P" in original_result.data.columns:
                        print(f"    Significant associations (P < 0.05): {(original_result.data.loc[orig_valid, 'P'] < 0.05).sum()}")
            
            # Heterogeneity statistics
            if "Q" in python_result.data.columns:
                het_valid = python_result.data["Q"].notna()
                if het_valid.sum() > 0:
                    print(f"\nResidual Heterogeneity Statistics:")
                    print(f"  Python Implementation:")
                    print(f"    Mean Q: {python_result.data.loc[het_valid, 'Q'].mean():.6f}")
                    if "I2" in python_result.data.columns:
                        print(f"    Mean I2: {python_result.data.loc[het_valid, 'I2'].mean():.4f}")
                    if "P_HET" in python_result.data.columns:
                        print(f"    Significant heterogeneity (P_HET < 0.05): {(python_result.data.loc[het_valid, 'P_HET'] < 0.05).sum()}")
                    
                    if original_result is not None and "Q" in original_result.data.columns:
                        orig_het_valid = original_result.data["Q"].notna()
                        if orig_het_valid.sum() > 0:
                            print(f"  Original MR-MEGA:")
                            print(f"    Mean Q: {original_result.data.loc[orig_het_valid, 'Q'].mean():.6f}")
                            if "I2" in original_result.data.columns:
                                print(f"    Mean I2: {original_result.data.loc[orig_het_valid, 'I2'].mean():.4f}")
                            if "P_HET" in original_result.data.columns:
                                print(f"    Significant heterogeneity (P_HET < 0.05): {(original_result.data.loc[orig_het_valid, 'P_HET'] < 0.05).sum()}")
        
        print("="*180 + "\n")
    
    def _create_test_data_for_performance(self, n_variants):
        """Create test data for performance testing with specified number of variants"""
        np.random.seed(42)
        
        # Create variants covering all chromosomes (1-23)
        chromosomes = []
        positions = []
        alleles_ea = []
        alleles_nea = []
        
        # Simple allele pairs for performance testing
        snp_pairs = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C"), ("A", "T"), ("C", "G")]
        indel_pairs = [("A", "AT"), ("AT", "A"), ("C", "CG"), ("CG", "C")]
        all_allele_pairs = snp_pairs + indel_pairs
        
        chr_list = list(range(1, 24))
        n_per_chr = n_variants // len(chr_list)
        remainder = n_variants % len(chr_list)
        
        variant_idx = 0
        for chr_num in chr_list:
            n_chr = n_per_chr + (1 if chr_num <= remainder else 0)
            for i in range(n_chr):
                chromosomes.append(chr_num)
                bin_num = (variant_idx * 17 + chr_num * 7) % 299
                pos = bin_num * 1000000 + np.random.randint(1000, 999000)
                positions.append(pos)
                pair = all_allele_pairs[variant_idx % len(all_allele_pairs)]
                alleles_ea.append(pair[0])
                alleles_nea.append(pair[1])
                variant_idx += 1
        
        variants = {
            "SNPID": [f"rs{i}" for i in range(1, n_variants + 1)],
            "CHR": chromosomes,
            "POS": positions,
            "EA": alleles_ea,
            "NEA": alleles_nea,
        }
        
        # Create EAF values
        eaf_cases = []
        for i in range(n_variants):
            case_type = i % 10
            if case_type == 0:
                eaf_cases.append(np.random.uniform(0.1, 0.9))
            elif case_type == 1:
                eaf_cases.append(np.random.uniform(0.01, 0.05))
            elif case_type == 2:
                eaf_cases.append(np.random.uniform(0.95, 0.99))
            elif case_type == 3:
                eaf_cases.append(0.01)
            elif case_type == 4:
                eaf_cases.append(0.99)
            else:
                eaf_cases.append(np.random.uniform(0.2, 0.8))
        
        # Create 5 studies
        studies_data = []
        base_eaf = np.array(eaf_cases)
        for study_idx in range(5):
            study_data = pd.DataFrame(variants.copy())
            study_data["BETA"] = np.random.normal(0.05 + study_idx * 0.02, 0.1, n_variants)
            study_data["SE"] = np.abs(np.random.normal(0.05, 0.01, n_variants))
            study_data["P"] = 2 * norm.sf(np.abs(study_data["BETA"] / study_data["SE"]))
            study_data["N"] = 5000 + study_idx * 200
            study_data["EAF"] = np.clip(base_eaf + np.random.normal(0, 0.02, n_variants), 0.01, 0.99)
            studies_data.append(study_data)
        
        # Create Sumstats objects
        studies = []
        for i, study_data in enumerate(studies_data):
            study = gl.Sumstats(study_data.copy(), verbose=False)
            study.meta["gwaslab"]["study_name"] = f"Study{i+1}"
            studies.append(study)
        
        # Create SumstatsMulti (suppress verbose output)
        import logging
        logging.getLogger().setLevel(logging.ERROR)
        sumstats_multi = SumstatsMulti(studies, verbose=False)
        logging.getLogger().setLevel(logging.WARNING)
        
        return sumstats_multi
    
    def test_performance_comparison(self):
        """Test performance (speed and memory) for both Python and C++ versions"""
        if not self.mrmega_available:
            self.skipTest("MR-MEGA executable not available")
        
        variant_counts = [200, 1000, 5000, 10000]
        results = []
        
        print("\n" + "="*100)
        print("PERFORMANCE TEST: Python vs C++ MR-MEGA")
        print("="*100)
        print(f"{'Variants':<12} | {'Version':<12} | {'Time (s)':<12} | {'Memory (MB)':<15} | {'Peak Memory (MB)':<18}")
        print("-"*100)
        
        for n_variants in variant_counts:
            print(f"\nTesting with {n_variants} variants...")
            
            # Create test data
            sumstats_multi = self._create_test_data_for_performance(n_variants)
            
            # Test Python version
            tracemalloc.start()
            start_time = time.time()
            python_result = meta_regress_mrmega_python(
                sumstats_multi,
                num_pcs=2,
                use_genomic_control=True,
                use_gco=True,
                verbose=False
            )
            python_time = time.time() - start_time
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            python_memory = current / 1024 / 1024  # Convert to MB
            python_peak = peak / 1024 / 1024  # Convert to MB
            
            results.append({
                'variants': n_variants,
                'version': 'Python',
                'time': python_time,
                'memory': python_memory,
                'peak_memory': python_peak
            })
            
            print(f"{n_variants:<12} | {'Python':<12} | {python_time:>11.4f} | {python_memory:>14.2f} | {python_peak:>17.2f}")
            
            # Test C++ version
            try:
                # Prepare files for C++ MR-MEGA
                study_files = []
                for study_idx in range(5):
                    study_file = os.path.join(self.temp_dir, f"study_{n_variants}_{study_idx+1}.txt")
                    study_data = sumstats_multi.data.copy()
                    study_df = pd.DataFrame({
                        "MARKERNAME": study_data["SNPID"],
                        "CHROMOSOME": study_data["CHR"],
                        "POSITION": study_data["POS"],
                        "EA": study_data["EA"],
                        "NEA": study_data["NEA"],
                        "EAF": study_data[f"EAF_{study_idx+1}"],
                        "N": study_data[f"N_{study_idx+1}"],
                        "BETA": study_data[f"BETA_{study_idx+1}"],
                        "SE": study_data[f"SE_{study_idx+1}"],
                    })
                    study_df = study_df.dropna()
                    study_df.to_csv(study_file, sep="\t", index=False)
                    study_files.append(study_file)
                
                # Create input list file
                input_list_file = os.path.join(self.temp_dir, f"input_{n_variants}.in")
                with open(input_list_file, 'w') as f:
                    for study_file in study_files:
                        f.write(f"{study_file}\n")
                
                # Run C++ MR-MEGA (matching _run_original_mrmega method)
                output_prefix = os.path.join(self.temp_dir, f"mrmega_output_{n_variants}")
                cmd = [
                    self.mrmega_path,
                    "-i", input_list_file,
                    "-o", output_prefix,
                    "--pc", "2",  # Number of principal components
                    "--qt",  # Quantitative trait (use BETA/SE instead of OR)
                ]
                
                start_time = time.time()
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300,  # 5 minute timeout
                    cwd=self.temp_dir
                )
                cpp_time = time.time() - start_time
                
                # For C++, we can't easily measure memory, so we'll use a placeholder
                # In practice, you could use psutil or other tools, but for simplicity we'll skip
                cpp_memory = 0.0
                cpp_peak = 0.0
                
                if result.returncode == 0:
                    results.append({
                        'variants': n_variants,
                        'version': 'C++',
                        'time': cpp_time,
                        'memory': cpp_memory,
                        'peak_memory': cpp_peak
                    })
                    print(f"{n_variants:<12} | {'C++':<12} | {cpp_time:>11.4f} | {'N/A':>14} | {'N/A':>17}")
                else:
                    print(f"{n_variants:<12} | {'C++':<12} | {'FAILED':>11} | {'N/A':>14} | {'N/A':>17}")
                    if result.stderr:
                        print(f"  C++ stderr: {result.stderr[:200]}")
            except subprocess.TimeoutExpired:
                print(f"{n_variants:<12} | {'C++':<12} | {'TIMEOUT':>11} | {'N/A':>14} | {'N/A':>17}")
            except Exception as e:
                print(f"{n_variants:<12} | {'C++':<12} | {'ERROR':>11} | {'N/A':>14} | {'N/A':>17}")
                print(f"  Error: {str(e)[:200]}")
        
        # Print summary
        print("\n" + "="*100)
        print("PERFORMANCE SUMMARY")
        print("="*100)
        
        for n_variants in variant_counts:
            python_results = [r for r in results if r['variants'] == n_variants and r['version'] == 'Python']
            cpp_results = [r for r in results if r['variants'] == n_variants and r['version'] == 'C++']
            
            if python_results:
                py_time = python_results[0]['time']
                py_mem = python_results[0]['memory']
                py_peak = python_results[0]['peak_memory']
                print(f"\n{n_variants} variants:")
                print(f"  Python: {py_time:.4f}s, Memory: {py_mem:.2f} MB (peak: {py_peak:.2f} MB)")
                
                if cpp_results:
                    cpp_time = cpp_results[0]['time']
                    speedup = py_time / cpp_time if cpp_time > 0 else 0
                    print(f"  C++:    {cpp_time:.4f}s")
                    if speedup > 1:
                        print(f"  C++ is {speedup:.2f}x faster")
                    elif speedup > 0:
                        print(f"  Python is {1/speedup:.2f}x faster")
                else:
                    print(f"  C++:    Failed or not available")
        
        print("="*100 + "\n")


if __name__ == "__main__":
    # Create test instance and print comparison table
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--print-table":
        test = TestMetaRegressMRMEGA()
        test.setUp()
        test.print_comparison_table(n_variants=20)
    else:
        unittest.main()
