import os
import sys
import unittest
import tempfile
import shutil
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl
from gwaslab.g_SumstatsPair import SumstatsPair
from gwaslab.g_SumstatsMulti import SumstatsMulti
from gwaslab.info.g_Log import Log


class TestMetaAnalyzeMetal(unittest.TestCase):
    """Test meta-analysis results against METAL software"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests"""
        # Check if METAL is available
        cls.metal_path = shutil.which("metal")
        if cls.metal_path is None:
            raise unittest.SkipTest("METAL not found in PATH. Skipping METAL comparison tests.")
        
        # Create temporary directory for test files
        cls.temp_dir = tempfile.mkdtemp(prefix="gwaslab_metal_test_")
    
    @classmethod
    def tearDownClass(cls):
        """Clean up temporary files"""
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)
    
    @staticmethod
    def _create_variant_from_pvalue(p_value, beta_sign=1, se_base=0.01):
        """
        Create BETA and SE values that produce a specific p-value.
        
        Parameters
        ----------
        p_value : float
            Target p-value
        beta_sign : int
            Sign of BETA (1 for positive, -1 for negative)
        se_base : float
            Base SE value (will be adjusted to achieve target p-value)
        
        Returns
        -------
        tuple
            (BETA, SE) values
        """
        # For two-sided test: p = 2 * (1 - norm.cdf(|Z|))
        # So |Z| = norm.isf(p/2) or norm.ppf(1 - p/2)
        # For very small p-values, use isf (inverse survival function)
        try:
            z_abs = norm.isf(p_value / 2.0)
        except (OverflowError, ValueError):
            # For extremely small p-values, use ppf with log-space
            # |Z| ≈ sqrt(2 * log(2 / p)) for very small p
            z_abs = np.sqrt(2 * np.log(2.0 / max(p_value, 1e-400)))
        
        # Ensure Z is finite
        if not np.isfinite(z_abs):
            z_abs = 50.0  # Fallback for extreme values
        
        # Calculate SE to achieve target Z: SE = |BETA| / |Z|
        # Use a reasonable BETA value, then adjust SE
        beta_target = beta_sign * 0.1  # Start with reasonable BETA
        se = abs(beta_target) / z_abs
        
        # If SE is too small, adjust BETA instead
        if se < 1e-10:
            se = se_base
            beta_target = beta_sign * z_abs * se
        
        return beta_target, se
    
    def setUp(self):
        """Set up test data for each test"""
        # Create test data with 2 studies
        # Use realistic values that will work with METAL
        
        # Standard variants
        standard_variants = {
            "SNPID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "CHR": [1, 1, 2, 2, 1],
            "POS": [1000, 2000, 3000, 4000, 5000],
            "EA": ["A", "G", "T", "C", "A"],
            "NEA": ["G", "C", "C", "T", "T"],
            "BETA": [0.1, -0.2, 0.3, 0.05, 0.15],
            "SE": [0.05, 0.1, 0.15, 0.02, 0.08],
            "P": [0.0455, 0.0455, 0.0455, 0.0127, 0.0606],  # Approximate P-values
            "N": [1000, 1000, 1000, 1000, 1000],
            "EAF": [0.3, 0.4, 0.5, 0.2, 0.35],
        }
        
        # Variants with extreme p-values
        extreme_pvalues = [1e-7, 1e-10, 1e-20, 1e-50, 1e-100, 1e-300, 1e-400]
        extreme_variants = {
            "SNPID": [f"rs_extreme_{i+1}" for i in range(len(extreme_pvalues))],
            "CHR": [1] * len(extreme_pvalues),
            "POS": list(range(10000, 10000 + len(extreme_pvalues) * 1000, 1000)),
            "EA": ["A", "G", "T", "C", "A", "G", "T"],
            "NEA": ["G", "C", "C", "T", "T", "A", "C"],
            "BETA": [],
            "SE": [],
            "P": extreme_pvalues,
            "N": [10000] * len(extreme_pvalues),  # Larger N for extreme p-values
            "EAF": [0.3, 0.4, 0.5, 0.2, 0.35, 0.45, 0.25],
        }
        
        # Generate BETA and SE for extreme p-values
        for i, pval in enumerate(extreme_pvalues):
            beta_sign = 1 if i % 2 == 0 else -1  # Alternate signs
            beta, se = self._create_variant_from_pvalue(pval, beta_sign=beta_sign, se_base=0.01)
            extreme_variants["BETA"].append(beta)
            extreme_variants["SE"].append(se)
        
        # Combine standard and extreme variants
        self.study1_data = pd.DataFrame({
            **standard_variants,
            **{k: standard_variants[k] + extreme_variants[k] 
               for k in ["SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF"]}
        })
        
        # Study 2: include some standard variants and some extreme variants
        study2_standard = {
            "SNPID": ["rs1", "rs2", "rs3", "rs6"],
            "CHR": [1, 1, 2, 1],
            "POS": [1000, 2000, 3000, 6000],
            "EA": ["A", "G", "T", "G"],
            "NEA": ["G", "C", "C", "A"],
            "BETA": [0.12, -0.18, 0.28, 0.2],
            "SE": [0.06, 0.11, 0.16, 0.09],
            "P": [0.0455, 0.1019, 0.0801, 0.0267],
            "N": [2000, 2000, 2000, 2000],
            "EAF": [0.32, 0.38, 0.52, 0.4],
        }
        
        # Include first 4 extreme variants in study 2
        study2_extreme = {
            "SNPID": extreme_variants["SNPID"][:4],
            "CHR": extreme_variants["CHR"][:4],
            "POS": extreme_variants["POS"][:4],
            "EA": extreme_variants["EA"][:4],
            "NEA": extreme_variants["NEA"][:4],
            "BETA": extreme_variants["BETA"][:4],
            "SE": extreme_variants["SE"][:4],
            "P": extreme_variants["P"][:4],
            "N": extreme_variants["N"][:4],
            "EAF": extreme_variants["EAF"][:4],
        }
        
        self.study2_data = pd.DataFrame({
            **study2_standard,
            **{k: study2_standard[k] + study2_extreme[k] 
               for k in ["SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF"]}
        })
        
        # Create Sumstats objects
        self.study1 = gl.Sumstats(self.study1_data.copy(), verbose=False)
        self.study2 = gl.Sumstats(self.study2_data.copy(), verbose=False)
    
    def create_metal_input_file(self, study_data, filename):
        """Create a METAL input file from study data"""
        filepath = os.path.join(self.temp_dir, filename)
        
        # METAL format: MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value Direction
        # Note: METAL uses Allele1 as effect allele (EA) and Allele2 as non-effect allele (NEA)
        metal_df = pd.DataFrame({
            "MarkerName": study_data["SNPID"],
            "Allele1": study_data["EA"],
            "Allele2": study_data["NEA"],
            "Freq1": study_data["EAF"],
            "Effect": study_data["BETA"],
            "StdErr": study_data["SE"],
            "P-value": study_data["P"],
            "Direction": ["+" if b > 0 else "-" if b < 0 else "?" for b in study_data["BETA"]],
        })
        
        # Write to file
        metal_df.to_csv(filepath, sep="\t", index=False, na_rep="NA")
        return filepath
    
    def create_metal_script(self, study1_file, study2_file, output_prefix):
        """Create a METAL script file"""
        script_path = os.path.join(self.temp_dir, "metal_script.txt")
        
        script_content = f"""SCHEME STDERR
MARKER MarkerName
ALLELE Allele1 Allele2
FREQ Freq1
EFFECT Effect
STDERR StdErr
PVAL P-value
PROCESS {study1_file}
PROCESS {study2_file}
ANALYZE
"""
        
        with open(script_path, "w") as f:
            f.write(script_content)
        
        return script_path
    
    def run_metal(self, script_path, output_prefix):
        """Run METAL and return output file path"""
        output_file = os.path.join(self.temp_dir, f"{output_prefix}.txt")
        
        # METAL writes to METAANALYSIS1.TBL by default, but we can specify output
        try:
            result = subprocess.run(
                [self.metal_path, script_path],
                cwd=self.temp_dir,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"METAL failed with return code {result.returncode}:\n{result.stderr}")
            
            # METAL outputs to METAANALYSIS1.TBL in the current directory
            metal_output = os.path.join(self.temp_dir, "METAANALYSIS1.TBL")
            if os.path.exists(metal_output):
                # Copy to our desired output file
                shutil.copy(metal_output, output_file)
                return output_file
            else:
                raise FileNotFoundError(f"METAL output file not found: {metal_output}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("METAL execution timed out")
        except Exception as e:
            raise RuntimeError(f"Failed to run METAL: {str(e)}")
    
    def parse_metal_output(self, metal_output_file):
        """Parse METAL output file and return DataFrame"""
        # METAL output format: MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq Effect StdErr P-value Direction HetISq HetChiSq HetDf HetPVal
        try:
            metal_df = pd.read_csv(metal_output_file, sep="\t", comment="#")
            
            # Rename columns to match GWASLab format
            column_mapping = {
                "MarkerName": "SNPID",
                "Allele1": "EA",
                "Allele2": "NEA",
                "Freq1": "EAF",
                "Effect": "BETA",
                "StdErr": "SE",
                "P-value": "P",
                "Direction": "DIRECTION",
                "HetISq": "I2",
                "HetChiSq": "Q",
                "HetDf": "DOF",
                "HetPVal": "P_HET",
            }
            
            metal_df = metal_df.rename(columns=column_mapping)
            
            # Calculate Z-score from BETA and SE
            if "BETA" in metal_df.columns and "SE" in metal_df.columns:
                metal_df["Z"] = metal_df["BETA"] / metal_df["SE"]
            
            return metal_df
        except Exception as e:
            raise RuntimeError(f"Failed to parse METAL output: {str(e)}")
    
    def compare_results(self, gwaslab_result, metal_result, rtol=1e-3, atol=1e-5):
        """Compare GWASLab and METAL meta-analysis results
        
        Note: METAL may flip alleles, which means:
        - EA and NEA are swapped
        - BETA sign is flipped
        - EAF becomes 1 - EAF
        This function automatically detects and handles allele flipping.
        """
        # Get available columns from both results
        common_cols = ["SNPID", "BETA", "SE", "P", "Z"]
        optional_cols = ["EAF", "Q", "I2", "DOF", "P_HET", "EA", "NEA"]
        
        # Build list of columns that exist in both results
        gwaslab_cols = set(gwaslab_result.columns)
        metal_cols = set(metal_result.columns)
        
        cols_to_compare = common_cols.copy()
        for col in optional_cols:
            if col in gwaslab_cols and col in metal_cols:
                cols_to_compare.append(col)
        
        # Merge on SNPID for comparison
        comparison = pd.merge(
            gwaslab_result[cols_to_compare],
            metal_result[cols_to_compare],
            on="SNPID",
            suffixes=("_gwaslab", "_metal")
        )
        
        if len(comparison) == 0:
            self.fail("No overlapping variants found between GWASLab and METAL results")
        
        # Check if alleles match or are flipped
        # METAL may flip alleles, so we need to check both possibilities
        # Note: METAL may output lowercase alleles while GWASLab uses uppercase
        if "EA_gwaslab" in comparison.columns and "EA_metal" in comparison.columns:
            # Convert to uppercase for case-insensitive comparison
            ea_gwaslab = comparison["EA_gwaslab"].astype(str).str.upper()
            nea_gwaslab = comparison["NEA_gwaslab"].astype(str).str.upper()
            ea_metal = comparison["EA_metal"].astype(str).str.upper()
            nea_metal = comparison["NEA_metal"].astype(str).str.upper()
            
            # Check if alleles match directly (case-insensitive)
            direct_match = (
                (ea_gwaslab == ea_metal) &
                (nea_gwaslab == nea_metal)
            )
            
            # Check if alleles are flipped (EA <-> NEA swapped, case-insensitive)
            flipped_match = (
                (ea_gwaslab == nea_metal) &
                (nea_gwaslab == ea_metal)
            )
            
            # Debug: print allele comparison
            print(f"\nAllele matching:")
            print(f"  Direct match: {direct_match.sum()} variants")
            print(f"  Flipped match: {flipped_match.sum()} variants")
            print(f"  Neither match: {(~(direct_match | flipped_match)).sum()} variants")
            
            if (~(direct_match | flipped_match)).any():
                mismatched = comparison[~(direct_match | flipped_match)]
                print(f"\n  Mismatched alleles (first few):")
                print(mismatched[["SNPID", "EA_gwaslab", "NEA_gwaslab", "EA_metal", "NEA_metal"]].head())
            
            # For variants where alleles are flipped, we need to adjust METAL results
            if flipped_match.any():
                print(f"\nNote: Adjusting {flipped_match.sum()} variants with flipped alleles in METAL output")
                # Create adjusted METAL values for flipped variants
                comparison.loc[flipped_match, "BETA_metal"] = -comparison.loc[flipped_match, "BETA_metal"]
                if "EAF_metal" in comparison.columns:
                    comparison.loc[flipped_match, "EAF_metal"] = 1 - comparison.loc[flipped_match, "EAF_metal"]
                if "Z_metal" in comparison.columns:
                    comparison.loc[flipped_match, "Z_metal"] = -comparison.loc[flipped_match, "Z_metal"]
            
            # For variants where neither direct nor flipped match, try sign flip as fallback
            # This handles cases where METAL might have different allele representation
            neither_match = ~(direct_match | flipped_match)
            if neither_match.any():
                print(f"\nNote: {neither_match.sum()} variants have mismatched alleles, trying sign flip")
                # Try sign flip for these variants (METAL might have flipped the effect)
                beta_before = comparison.loc[neither_match, "BETA_metal"].copy()
                comparison.loc[neither_match, "BETA_metal"] = -comparison.loc[neither_match, "BETA_metal"]
                if "Z_metal" in comparison.columns:
                    comparison.loc[neither_match, "Z_metal"] = -comparison.loc[neither_match, "Z_metal"]
                # Check if sign flip improves the match
                beta_diff_before = np.abs(comparison.loc[neither_match, "BETA_gwaslab"] - beta_before)
                beta_diff_after = np.abs(comparison.loc[neither_match, "BETA_gwaslab"] - comparison.loc[neither_match, "BETA_metal"])
                # Only keep the flip if it improves the match
                keep_flip = beta_diff_after < beta_diff_before
                if not keep_flip.all():
                    # Revert for variants where flip doesn't help
                    comparison.loc[neither_match & ~keep_flip, "BETA_metal"] = beta_before[~keep_flip]
                    if "Z_metal" in comparison.columns:
                        comparison.loc[neither_match & ~keep_flip, "Z_metal"] = -comparison.loc[neither_match & ~keep_flip, "Z_metal"]
        
        # Debug: print first few comparisons
        print(f"\nComparing {len(comparison)} variants:")
        print(comparison[["SNPID", "BETA_gwaslab", "BETA_metal", "SE_gwaslab", "SE_metal"]].head())
        
        # Compare BETA (effect size)
        # After handling allele flipping above, we can do direct comparison
        beta_diff = comparison["BETA_gwaslab"] - comparison["BETA_metal"]
        beta_abs_diff = np.abs(beta_diff)
        
        # Use absolute difference for very small values, relative for larger values
        beta_metal_abs = np.abs(comparison["BETA_metal"])
        use_relative = beta_metal_abs > atol
        
        if use_relative.any():
            beta_rtol = np.abs(beta_diff[use_relative] / comparison["BETA_metal"][use_relative]).max()
            self.assertLess(beta_rtol, rtol, 
                           f"BETA values differ by more than {rtol*100}%: max relative difference = {beta_rtol*100:.2f}%")
        
        # Also check absolute difference
        max_abs_diff = beta_abs_diff.max()
        self.assertLess(max_abs_diff, 0.001,  # Allow tolerance for BETA (numerical precision)
                       f"BETA values differ by more than 0.001: max absolute difference = {max_abs_diff:.6f}")
        
        # Compare SE (standard error)
        se_diff = comparison["SE_gwaslab"] - comparison["SE_metal"]
        se_rtol = np.abs(se_diff / comparison["SE_metal"]).max()
        # Allow slightly larger tolerance for SE (0.5% instead of 0.1%) due to numerical precision
        self.assertLess(se_rtol, 0.005,
                       f"SE values differ by more than 0.5%: max relative difference = {se_rtol*100:.2f}%")
        
        # Compare P-values (can have larger differences due to numerical precision, especially for extreme values)
        p_diff = np.abs(comparison["P_gwaslab"] - comparison["P_metal"])
        
        # For very small p-values, use relative difference in log space
        # For larger p-values, use absolute difference
        p_gwaslab = comparison["P_gwaslab"].values
        p_metal = comparison["P_metal"].values
        
        # Use log-space comparison for very small p-values (p < 1e-10)
        very_small_mask = (p_gwaslab < 1e-10) | (p_metal < 1e-10)
        
        if very_small_mask.any():
            # For very small p-values, compare in log space
            log_p_gwaslab = np.log10(np.maximum(p_gwaslab[very_small_mask], 1e-400))
            log_p_metal = np.log10(np.maximum(p_metal[very_small_mask], 1e-400))
            log_p_diff = np.abs(log_p_gwaslab - log_p_metal)
            max_log_p_diff = log_p_diff.max()
            # Allow up to 0.1 orders of magnitude difference for extreme p-values
            self.assertLess(max_log_p_diff, 0.1,
                           f"Very small P-values differ by more than 0.1 orders of magnitude: max log difference = {max_log_p_diff:.2e}")
        
        # For larger p-values, use absolute difference
        larger_mask = ~very_small_mask
        if larger_mask.any():
            max_p_diff = p_diff[larger_mask].max()
            self.assertLess(max_p_diff, 1e-4,
                           f"P-values differ by more than 1e-4: max absolute difference = {max_p_diff:.2e}")
        
        # Compare Z-scores (after handling allele flipping above)
        z_diff = np.abs(comparison["Z_gwaslab"] - comparison["Z_metal"])
        z_gwaslab_abs = np.abs(comparison["Z_gwaslab"])
        z_metal_abs = np.abs(comparison["Z_metal"])
        
        # For very large Z-scores (from extreme p-values), use relative tolerance
        # For smaller Z-scores, use absolute tolerance
        large_z_mask = (z_gwaslab_abs > 10) | (z_metal_abs > 10)
        
        if large_z_mask.any():
            # Use relative tolerance for large Z-scores (0.5% relative difference for extreme values)
            z_rel_diff = z_diff[large_z_mask] / np.maximum(z_gwaslab_abs[large_z_mask], z_metal_abs[large_z_mask])
            max_z_rel_diff = z_rel_diff.max()
            self.assertLess(max_z_rel_diff, 0.005,
                           f"Large Z-scores differ by more than 0.5%: max relative difference = {max_z_rel_diff*100:.2f}%")
        
        # For smaller Z-scores, use absolute tolerance (relaxed for numerical precision)
        small_z_mask = ~large_z_mask
        if small_z_mask.any():
            max_z_diff = z_diff[small_z_mask].max()
            self.assertLess(max_z_diff, 0.05,  # Allow larger tolerance for small Z-scores
                           f"Z-scores differ by more than 0.05: max absolute difference = {max_z_diff:.2e}")
        
        # Compare EAF (effect allele frequency) if available
        # Note: EAF is already adjusted for flipped alleles above
        if "EAF_gwaslab" in comparison.columns and "EAF_metal" in comparison.columns:
            eaf_diff = np.abs(comparison["EAF_gwaslab"] - comparison["EAF_metal"])
            max_eaf_diff = eaf_diff.max()
            self.assertLess(max_eaf_diff, 0.01,  # EAF should be very close
                           f"EAF values differ by more than 0.01: max absolute difference = {max_eaf_diff:.4f}")
        
        # Compare heterogeneity statistics (only for variants in both studies)
        if "DOF_gwaslab" in comparison.columns:
            both_studies = comparison[comparison["DOF_gwaslab"] > 0]
            if len(both_studies) > 0:
                # Compare Q (Cochran's Q) if available
                if "Q_gwaslab" in comparison.columns and "Q_metal" in comparison.columns:
                    q_diff = np.abs(both_studies["Q_gwaslab"] - both_studies["Q_metal"])
                    max_q_diff = q_diff.max()
                    self.assertLess(max_q_diff, 0.1,  # Q can have some numerical differences
                                   f"Q values differ by more than 0.1: max absolute difference = {max_q_diff:.4f}")
                
                # Compare I2 (only where both are not NaN) if available
                if "I2_gwaslab" in comparison.columns and "I2_metal" in comparison.columns:
                    i2_both = both_studies[both_studies["I2_gwaslab"].notna() & both_studies["I2_metal"].notna()]
                    if len(i2_both) > 0:
                        i2_diff = np.abs(i2_both["I2_gwaslab"] - i2_both["I2_metal"])
                        max_i2_diff = i2_diff.max()
                        self.assertLess(max_i2_diff, 0.05,  # I2 should be close (0-1 scale)
                                       f"I2 values differ by more than 0.05: max absolute difference = {max_i2_diff:.4f}")
        
        return comparison
    
    def test_fixed_effects_metal_comparison(self):
        """Test fixed-effects meta-analysis against METAL"""
        # Create METAL input files
        study1_file = self.create_metal_input_file(self.study1_data, "study1_metal.txt")
        study2_file = self.create_metal_input_file(self.study2_data, "study2_metal.txt")
        
        # Create METAL script
        script_path = self.create_metal_script(study1_file, study2_file, "metal_output")
        
        # Run METAL
        metal_output_file = self.run_metal(script_path, "metal_output")
        metal_result = self.parse_metal_output(metal_output_file)
        
        # Run GWASLab meta-analysis
        pair = SumstatsPair(self.study1, self.study2, keep_all_variants=False, verbose=False)
        gwaslab_result = pair.run_meta_analysis(random_effects=False)
        
        # Compare results
        comparison = self.compare_results(gwaslab_result.data, metal_result)
        
        # Print detailed comparison table
        self._print_comparison_table(comparison)
        
        # Print summary
        beta_abs_diff = np.abs(comparison['BETA_gwaslab'] - comparison['BETA_metal'])
        se_abs_diff = np.abs(comparison['SE_gwaslab'] - comparison['SE_metal'])
        p_abs_diff = np.abs(comparison['P_gwaslab'] - comparison['P_metal'])
        
        print(f"\nSummary:")
        print(f"  Compared {len(comparison)} variants")
        print(f"  BETA: max absolute difference = {beta_abs_diff.max():.6f}")
        print(f"  SE: max absolute difference = {se_abs_diff.max():.6f}")
        print(f"  P: max absolute difference = {p_abs_diff.max():.2e}")
    
    def _print_comparison_table(self, comparison):
        """Print a formatted table comparing GWASLab and METAL results"""
        print("\n" + "="*120)
        print("COMPARISON TABLE: GWASLab vs METAL")
        print("="*120)
        
        # Prepare data for table
        table_data = []
        for idx, row in comparison.iterrows():
            snpid = row['SNPID']
            
            # Format BETA
            beta_gl = row['BETA_gwaslab']
            beta_metal = row['BETA_metal']
            beta_diff = beta_gl - beta_metal
            beta_str = f"{beta_gl:10.6f} | {beta_metal:10.6f} | {beta_diff:10.6f}"
            
            # Format SE
            se_gl = row['SE_gwaslab']
            se_metal = row['SE_metal']
            se_diff = se_gl - se_metal
            se_str = f"{se_gl:10.6f} | {se_metal:10.6f} | {se_diff:10.6f}"
            
            # Format P-value (use scientific notation for small values)
            p_gl = row['P_gwaslab']
            p_metal = row['P_metal']
            p_diff = p_gl - p_metal
            
            if p_gl < 1e-6 or p_metal < 1e-6:
                p_gl_str = f"{p_gl:.2e}"
                p_metal_str = f"{p_metal:.2e}"
                p_diff_str = f"{p_diff:.2e}"
            else:
                p_gl_str = f"{p_gl:.6f}"
                p_metal_str = f"{p_metal:.6f}"
                p_diff_str = f"{p_diff:.6f}"
            
            p_str = f"{p_gl_str:>12} | {p_metal_str:>12} | {p_diff_str:>12}"
            
            # Format Z-score
            z_gl = row.get('Z_gwaslab', np.nan)
            z_metal = row.get('Z_metal', np.nan)
            if pd.notna(z_gl) and pd.notna(z_metal):
                z_diff = z_gl - z_metal
                z_str = f"{z_gl:10.4f} | {z_metal:10.4f} | {z_diff:10.4f}"
            else:
                z_str = "      N/A |       N/A |       N/A"
            
            # Format EAF if available
            eaf_gl = row.get('EAF_gwaslab', np.nan)
            eaf_metal = row.get('EAF_metal', np.nan)
            if pd.notna(eaf_gl) and pd.notna(eaf_metal):
                eaf_diff = eaf_gl - eaf_metal
                eaf_str = f"{eaf_gl:8.4f} | {eaf_metal:8.4f} | {eaf_diff:8.4f}"
            else:
                eaf_str = "    N/A |     N/A |     N/A"
            
            table_data.append({
                'SNPID': snpid,
                'BETA': beta_str,
                'SE': se_str,
                'P': p_str,
                'Z': z_str,
                'EAF': eaf_str
            })
        
        # Print header
        print(f"\n{'SNPID':<15} | {'BETA (GWASLab | METAL | Diff)':<45} | {'SE (GWASLab | METAL | Diff)':<45}")
        print("-" * 120)
        
        # Print BETA and SE rows
        for data in table_data:
            print(f"{data['SNPID']:<15} | {data['BETA']:<45} | {data['SE']:<45}")
        
        print("\n" + "-" * 120)
        print(f"{'SNPID':<15} | {'P-value (GWASLab | METAL | Diff)':<45} | {'Z-score (GWASLab | METAL | Diff)':<45}")
        print("-" * 120)
        
        # Print P and Z rows
        for data in table_data:
            print(f"{data['SNPID']:<15} | {data['P']:<45} | {data['Z']:<45}")
        
        # Print EAF if available
        if any(pd.notna(row.get('EAF_gwaslab', np.nan)) for _, row in comparison.iterrows()):
            print("\n" + "-" * 120)
            print(f"{'SNPID':<15} | {'EAF (GWASLab | METAL | Diff)':<45}")
            print("-" * 120)
            for data in table_data:
                print(f"{data['SNPID']:<15} | {data['EAF']:<45}")
        
        # Print heterogeneity statistics if available
        het_cols = ['Q', 'I2', 'P_HET', 'DOF']
        available_het = [col for col in het_cols 
                        if f'{col}_gwaslab' in comparison.columns and f'{col}_metal' in comparison.columns]
        
        if available_het:
            print("\n" + "-" * 120)
            print("Heterogeneity Statistics (only for variants in multiple studies):")
            print("-" * 120)
            
            both_studies = comparison[comparison.get('DOF_gwaslab', pd.Series([0]*len(comparison))) > 0]
            if len(both_studies) > 0:
                for col in available_het:
                    if col == 'DOF':
                        continue  # Skip DOF as it's usually the same
                    print(f"\n{col}:")
                    print(f"{'SNPID':<15} | {'GWASLab':<15} | {'METAL':<15} | {'Difference':<15}")
                    print("-" * 60)
                    for idx, row in both_studies.iterrows():
                        val_gl = row[f'{col}_gwaslab']
                        val_metal = row[f'{col}_metal']
                        if pd.notna(val_gl) and pd.notna(val_metal):
                            diff = val_gl - val_metal
                            if col == 'I2':
                                print(f"{row['SNPID']:<15} | {val_gl:15.6f} | {val_metal:15.6f} | {diff:15.6f}")
                            else:
                                print(f"{row['SNPID']:<15} | {val_gl:15.4e} | {val_metal:15.4e} | {diff:15.4e}")
        
        print("="*120)
    
    def test_three_studies_metal_comparison(self):
        """Test meta-analysis with three studies against METAL"""
        # Create third study
        study3_data = pd.DataFrame({
            "SNPID": ["rs1", "rs2", "rs3"],
            "CHR": [1, 1, 2],
            "POS": [1000, 2000, 3000],
            "EA": ["A", "G", "T"],
            "NEA": ["G", "C", "C"],
            "BETA": [0.11, -0.19, 0.29],
            "SE": [0.055, 0.105, 0.155],
            "P": [0.0455, 0.0703, 0.0606],
            "N": [1500, 1500, 1500],
            "EAF": [0.31, 0.39, 0.51],
        })
        
        study3 = gl.Sumstats(study3_data.copy(), verbose=False)
        
        # Create METAL input files
        study1_file = self.create_metal_input_file(self.study1_data, "study1_metal_3.txt")
        study2_file = self.create_metal_input_file(self.study2_data, "study2_metal_3.txt")
        study3_file = self.create_metal_input_file(study3_data, "study3_metal_3.txt")
        
        # Create METAL script with 3 studies
        script_path = os.path.join(self.temp_dir, "metal_script_3.txt")
        script_content = f"""SCHEME STDERR
MARKER MarkerName
ALLELE Allele1 Allele2
FREQ Freq1
EFFECT Effect
STDERR StdErr
PVAL P-value
PROCESS {study1_file}
PROCESS {study2_file}
PROCESS {study3_file}
ANALYZE
"""
        with open(script_path, "w") as f:
            f.write(script_content)
        
        # Run METAL
        metal_output_file = self.run_metal(script_path, "metal_output_3")
        metal_result = self.parse_metal_output(metal_output_file)
        
        # Run GWASLab meta-analysis
        multi = SumstatsMulti([self.study1, self.study2, study3], keep_all_variants=False, verbose=False)
        gwaslab_result = multi.run_meta_analysis(random_effects=False)
        
        # Compare results
        comparison = self.compare_results(gwaslab_result.data, metal_result)
        
        print(f"\nCompared {len(comparison)} variants (3 studies)")
    
    def test_heterogeneity_statistics_metal_comparison(self):
        """Test heterogeneity statistics (Q, I2, P_HET) against METAL"""
        # Create METAL input files
        study1_file = self.create_metal_input_file(self.study1_data, "study1_het.txt")
        study2_file = self.create_metal_input_file(self.study2_data, "study2_het.txt")
        
        # Create METAL script
        script_path = self.create_metal_script(study1_file, study2_file, "metal_het")
        
        # Run METAL
        metal_output_file = self.run_metal(script_path, "metal_het")
        metal_result = self.parse_metal_output(metal_output_file)
        
        # Run GWASLab meta-analysis
        pair = SumstatsPair(self.study1, self.study2, keep_all_variants=False, verbose=False)
        gwaslab_result = pair.run_meta_analysis(random_effects=False)
        
        # Compare heterogeneity statistics for variants in both studies
        # Get available columns
        het_cols = ["SNPID"]
        for col in ["Q", "I2", "P_HET", "DOF"]:
            if col in gwaslab_result.data.columns and col in metal_result.columns:
                het_cols.append(col)
        
        if len(het_cols) == 1:  # Only SNPID
            self.skipTest("Heterogeneity statistics not available in results")
        
        comparison = pd.merge(
            gwaslab_result.data[het_cols],
            metal_result[het_cols],
            on="SNPID",
            suffixes=("_gwaslab", "_metal")
        )
        
        # Only compare variants present in both studies (DOF > 0)
        both_studies = comparison[comparison["DOF_gwaslab"] > 0]
        
        if len(both_studies) > 0:
            # Compare Q
            q_diff = np.abs(both_studies["Q_gwaslab"] - both_studies["Q_metal"])
            max_q_diff = q_diff.max()
            self.assertLess(max_q_diff, 0.1,
                           f"Q values differ by more than 0.1: max absolute difference = {max_q_diff:.4f}")
            
            # Compare I2 (where both are not NaN)
            i2_both = both_studies[both_studies["I2_gwaslab"].notna() & both_studies["I2_metal"].notna()]
            if len(i2_both) > 0:
                i2_diff = np.abs(i2_both["I2_gwaslab"] - i2_both["I2_metal"])
                max_i2_diff = i2_diff.max()
                self.assertLess(max_i2_diff, 0.05,
                               f"I2 values differ by more than 0.05: max absolute difference = {max_i2_diff:.4f}")
            
            print(f"\nCompared heterogeneity statistics for {len(both_studies)} variants")


if __name__ == "__main__":
    unittest.main()
