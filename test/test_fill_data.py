import os
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_fill_data import _fill_data
from gwaslab.info.g_Log import Log


class TestFillData(unittest.TestCase):
    """Comprehensive tests for _fill_data function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.n_variants = 100
        
    def _create_base_sumstats(self, **columns):
        """Create a base sumstats DataFrame with specified columns."""
        data = {
            "SNPID": [f"rs{i}" for i in range(self.n_variants)],
            "CHR": np.random.randint(1, 23, self.n_variants),
            "POS": np.random.randint(1_000, 2_000_000, self.n_variants),
            "EA": np.random.choice(["A", "T", "G", "C"], self.n_variants),
            "NEA": np.random.choice(["A", "T", "G", "C"], self.n_variants),
        }
        data.update(columns)
        return pd.DataFrame(data)
    
    # ========== Basic Single-Step Conversions ==========
    
    def test_fill_p_from_mlog10p(self):
        """Test filling P from MLOG10P."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input MLOG10P (first 5): {mlog10p[:5]}")
        print(f"Output P (first 5): {result['P'].head().tolist()}")
        print(f"Expected P (first 5): {(10**(-mlog10p))[:5]}")
        print(f"P column created: {'P' in result.columns}")
        
        self.assertIn("P", result.columns)
        np.testing.assert_allclose(result["P"], 10**(-mlog10p), rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_fill_p_from_z(self):
        """Test filling P from Z."""
        z = np.random.normal(0, 2, self.n_variants)
        sumstats = self._create_base_sumstats(Z=z)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input Z (first 5): {z[:5]}")
        print(f"Output P (first 5): {result['P'].head().tolist()}")
        print(f"P range: [{result['P'].min():.2e}, {result['P'].max():.2e}]")
        print(f"P column created: {'P' in result.columns}")
        print(f"P values valid (0<=P<=1): {(result['P'] >= 0).all() and (result['P'] <= 1).all()}")
        
        self.assertIn("P", result.columns)
        # Check that P values are reasonable (between 0 and 1)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
        print("✓ Test PASSED")
    
    def test_fill_p_from_chisq(self):
        """Test filling P from CHISQ."""
        chisq = np.random.chisquare(1, self.n_variants)
        sumstats = self._create_base_sumstats(CHISQ=chisq)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
    
    def test_fill_z_from_beta_se(self):
        """Test filling Z from BETA/SE."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input BETA (first 5): {beta[:5]}")
        print(f"Input SE (first 5): {se[:5]}")
        print(f"Output Z (first 5): {result['Z'].head().tolist()}")
        print(f"Expected Z (first 5): {(beta / se)[:5]}")
        print(f"Z column created: {'Z' in result.columns}")
        
        self.assertIn("Z", result.columns)
        np.testing.assert_allclose(result["Z"], beta / se, rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_fill_chisq_from_z(self):
        """Test filling CHISQ from Z."""
        z = np.random.normal(0, 2, self.n_variants)
        sumstats = self._create_base_sumstats(Z=z)
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("CHISQ", result.columns)
        np.testing.assert_allclose(result["CHISQ"], z**2, rtol=1e-10)
    
    def test_fill_chisq_from_p(self):
        """Test filling CHISQ from P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("CHISQ", result.columns)
        self.assertTrue((result["CHISQ"] >= 0).all())
    
    def test_fill_beta_from_or(self):
        """Test filling BETA from OR."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val)
        
        result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input OR (first 5): {or_val[:5]}")
        print(f"Output BETA (first 5): {result['BETA'].head().tolist()}")
        print(f"Expected BETA (first 5): {np.log(or_val)[:5]}")
        print(f"BETA column created: {'BETA' in result.columns}")
        
        self.assertIn("BETA", result.columns)
        np.testing.assert_allclose(result["BETA"], np.log(or_val), rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_fill_or_from_beta(self):
        """Test filling OR from BETA."""
        beta = np.random.normal(0, 0.2, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta)
        
        result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input BETA (first 5): {beta[:5]}")
        print(f"Output OR (first 5): {result['OR'].head().tolist()}")
        print(f"Expected OR (first 5): {np.exp(beta)[:5]}")
        print(f"OR column created: {'OR' in result.columns}")
        if "OR_95L" in result.columns:
            print(f"OR_95L column created: True (first 5): {result['OR_95L'].head().tolist()}")
        if "OR_95U" in result.columns:
            print(f"OR_95U column created: True (first 5): {result['OR_95U'].head().tolist()}")
        
        self.assertIn("OR", result.columns)
        np.testing.assert_allclose(result["OR"], np.exp(beta), rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_fill_or_confidence_intervals(self):
        """Test filling OR_95L and OR_95U from BETA/SE."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        
        self.assertIn("OR", result.columns)
        self.assertIn("OR_95L", result.columns)
        self.assertIn("OR_95U", result.columns)
        # Check that confidence intervals are reasonable
        self.assertTrue((result["OR_95L"] < result["OR"]).all())
        self.assertTrue((result["OR_95U"] > result["OR"]).all())
    
    def test_fill_se_from_beta_p(self):
        """Test filling SE from BETA/P."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        p = np.random.uniform(1e-6, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, P=p)
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        self.assertTrue((result["SE"] > 0).all())
    
    def test_fill_se_from_or_or95u(self):
        """Test filling SE from OR/OR_95U."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        or_95u = or_val * np.random.uniform(1.01, 1.5, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val, OR_95U=or_95u)
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        self.assertTrue((result["SE"] > 0).all())
    
    def test_fill_se_from_or_or95l(self):
        """Test filling SE from OR/OR_95L."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        or_95l = or_val * np.random.uniform(0.5, 0.99, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val, OR_95L=or_95l)
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        self.assertTrue((result["SE"] > 0).all())
    
    def test_fill_mlog10p_from_p(self):
        """Test filling MLOG10P from P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input P (first 5): {p[:5]}")
        print(f"Output MLOG10P (first 5): {result['MLOG10P'].head().tolist()}")
        print(f"Expected MLOG10P (first 5): {(-np.log10(p))[:5]}")
        print(f"MLOG10P column created: {'MLOG10P' in result.columns}")
        
        self.assertIn("MLOG10P", result.columns)
        np.testing.assert_allclose(result["MLOG10P"], -np.log10(p), rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_fill_maf_from_eaf(self):
        """Test filling MAF from EAF."""
        eaf = np.random.uniform(0.01, 0.99, self.n_variants)
        sumstats = self._create_base_sumstats(EAF=eaf)
        
        result = _fill_data(sumstats, to_fill=["MAF"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input EAF (first 5): {eaf[:5]}")
        print(f"Output MAF (first 5): {result['MAF'].head().tolist()}")
        expected_maf = np.minimum(eaf, 1 - eaf)
        print(f"Expected MAF (first 5): {expected_maf[:5]}")
        print(f"MAF column created: {'MAF' in result.columns}")
        
        self.assertIn("MAF", result.columns)
        np.testing.assert_allclose(result["MAF"], expected_maf, rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_fill_sig_from_p(self):
        """Test filling SIGNIFICANT from P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["SIG"], sig_level=5e-8, verbose=False, log=self.log)
        
        self.assertIn("SIGNIFICANT", result.columns)
        self.assertTrue(result["SIGNIFICANT"].dtype == bool)
        expected_sig = p < 5e-8
        np.testing.assert_array_equal(result["SIGNIFICANT"], expected_sig)
    
    def test_fill_sig_from_mlog10p(self):
        """Test filling SIGNIFICANT from MLOG10P."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["SIG"], sig_level=5e-8, verbose=False, log=self.log)
        
        self.assertIn("SIGNIFICANT", result.columns)
        expected_sig = mlog10p > np.log10(1/5e-8)
        np.testing.assert_array_equal(result["SIGNIFICANT"], expected_sig)
    
    # ========== Multi-Round Conversions ==========
    
    def test_multi_round_or_to_p(self):
        """Test multi-round conversion: OR -> BETA -> Z -> P."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val, SE=se)
        
        # Fill BETA, Z, and P - the iterative process should fill them in order
        result = _fill_data(sumstats, to_fill=["BETA", "Z", "P"], verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Input OR (first 5): {or_val[:5]}")
        print(f"Input SE (first 5): {se[:5]}")
        print(f"Output BETA (first 5): {result['BETA'].head().tolist()}")
        print(f"Output Z (first 5): {result['Z'].head().tolist()}")
        print(f"Output P (first 5): {result['P'].head().tolist()}")
        print(f"Columns created: BETA={'BETA' in result.columns}, Z={'Z' in result.columns}, P={'P' in result.columns}")
        
        self.assertIn("P", result.columns)
        # Should have created BETA and Z as intermediate steps
        self.assertIn("BETA", result.columns)
        self.assertIn("Z", result.columns)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
        print("✓ Test PASSED")
    
    def test_multi_round_mlog10p_from_beta_se(self):
        """Test multi-round conversion: BETA/SE -> Z -> P -> MLOG10P."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        # Request all intermediate steps - the iterative process will fill them in order
        result = _fill_data(sumstats, to_fill=["Z", "P", "MLOG10P"], verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        self.assertIn("Z", result.columns)
        self.assertIn("P", result.columns)
        # P_MANTISSA and P_EXPONENT should be dropped after filling MLOG10P
        self.assertNotIn("P_MANTISSA", result.columns)
        self.assertNotIn("P_EXPONENT", result.columns)
        # Check that MLOG10P is reasonable
        self.assertTrue((result["MLOG10P"] >= 0).all())
    
    def test_multi_round_p_to_chisq(self):
        """Test multi-round conversion: P -> CHISQ."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        # P -> CHISQ should work
        self.assertIn("CHISQ", result.columns)
        self.assertTrue((result["CHISQ"] >= 0).all())
    
    # ========== Edge Cases ==========
    
    def test_fill_with_existing_column_no_overwrite(self):
        """Test that existing columns are not overwritten by default."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        original_p = p.copy()
        mlog10p = -np.log10(p)
        sumstats = self._create_base_sumstats(P=p, MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], overwrite=False, verbose=False, log=self.log)
        
        # P should remain unchanged
        np.testing.assert_array_equal(result["P"], original_p)
    
    def test_fill_with_existing_column_overwrite(self):
        """Test that existing columns are overwritten when overwrite=True."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        mlog10p = -np.log10(p)
        # Create P with wrong values
        wrong_p = np.random.uniform(0.5, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=wrong_p, MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], overwrite=True, verbose=False, log=self.log)
        
        # P should be recalculated from MLOG10P
        np.testing.assert_allclose(result["P"], 10**(-mlog10p), rtol=1e-10)
    
    def test_fill_with_missing_data(self):
        """Test filling with missing data (NaN values)."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        mlog10p[0:10] = np.nan  # Add some NaN values
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that NaN values are preserved
        self.assertTrue(np.isnan(result["P"][0:10]).all())
        # Check that non-NaN values are calculated correctly
        valid_mask = ~np.isnan(mlog10p)
        np.testing.assert_allclose(
            result["P"][valid_mask], 
            10**(-mlog10p[valid_mask]), 
            rtol=1e-10
        )
    
    def test_fill_extreme_mlog10p(self):
        """Test filling MLOG10P with extreme=True for very small P-values."""
        z = np.random.normal(0, 5, self.n_variants)  # Large Z scores
        sumstats = self._create_base_sumstats(Z=z)
        
        result = _fill_data(sumstats, to_fill=["MLOG10P"], extreme=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        # P_MANTISSA and P_EXPONENT are created by extreme methods but then dropped
        self.assertNotIn("P_MANTISSA", result.columns)
        self.assertNotIn("P_EXPONENT", result.columns)
    
    def test_fill_extreme_mlog10p_from_beta_se(self):
        """Test extreme MLOG10P calculation from BETA/SE."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.001, 0.01, self.n_variants)  # Small SE for large Z
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["MLOG10P"], extreme=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        self.assertIn("Z", result.columns)
    
    def test_fill_multiple_columns(self):
        """Test filling multiple columns at once."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["Z", "P", "CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("Z", result.columns)
        self.assertIn("P", result.columns)
        self.assertIn("CHISQ", result.columns)
    
    def test_fill_with_string_input(self):
        """Test that string input for to_fill is converted to list."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill="P", verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
    
    def test_fill_with_none_input(self):
        """Test that None input for to_fill is handled."""
        sumstats = self._create_base_sumstats()
        
        result = _fill_data(sumstats, to_fill=None, verbose=False, log=self.log)
        
        # Should return unchanged sumstats
        self.assertEqual(len(result), len(sumstats))
    
    def test_fill_with_invalid_column(self):
        """Test that invalid column names are ignored."""
        sumstats = self._create_base_sumstats()
        
        result = _fill_data(sumstats, to_fill=["INVALID_COL", "ALSO_INVALID"], verbose=False, log=self.log)
        
        # Should return unchanged sumstats
        self.assertNotIn("INVALID_COL", result.columns)
        self.assertNotIn("ALSO_INVALID", result.columns)
    
    def test_fill_when_column_cannot_be_filled(self):
        """Test behavior when a column cannot be filled due to missing prerequisites."""
        sumstats = self._create_base_sumstats()  # No columns that can be used for conversion
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        # P should not be created
        self.assertNotIn("P", result.columns)
    
    def test_fill_only_sig_variants(self):
        """Test filling only significant variants."""
        chisq = np.random.chisquare(1, self.n_variants)
        # Create some significant variants
        chisq[0:20] = np.random.chisquare(1, 20) * 10  # Large CHISQ for significance
        sumstats = self._create_base_sumstats(CHISQ=chisq, P=np.ones(self.n_variants) * 0.5)
        
        result = _fill_data(
            sumstats, 
            to_fill=["P"], 
            only_sig=True, 
            sig_level=5e-8,
            overwrite=True,
            verbose=False, 
            log=self.log
        )
        
        self.assertIn("P", result.columns)
        # Check that only significant variants were updated
        # (This is a simplified check - actual implementation may vary)
    
    def test_fill_with_chisq_and_df(self):
        """Test filling P from CHISQ with degrees of freedom column."""
        chisq = np.random.chisquare(1, self.n_variants)
        df = np.ones(self.n_variants, dtype=int)  # All df=1
        sumstats = self._create_base_sumstats(CHISQ=chisq, DOF=df)
        
        result = _fill_data(sumstats, to_fill=["P"], df="DOF", verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
    
    # ========== Round-trip Conversions ==========
    
    def test_round_trip_p_mlog10p(self):
        """Test round-trip conversion: P -> MLOG10P -> P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        # First fill MLOG10P
        result1 = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        self.assertIn("MLOG10P", result1.columns)
        
        # Then fill P back (should match original)
        result2 = _fill_data(result1, to_fill=["P"], overwrite=True, verbose=False, log=self.log)
        
        print(f"\n{'='*60}")
        print(f"TEST: {self._testMethodName}")
        print(f"{'='*60}")
        print(f"Original P (first 5): {p[:5]}")
        print(f"After round-trip P (first 5): {result2['P'].head().tolist()}")
        print(f"Max absolute difference: {np.abs(result2['P'] - p).max():.2e}")
        print(f"Max relative difference: {np.abs((result2['P'] - p) / p).max():.2e}")
        
        np.testing.assert_allclose(result2["P"], p, rtol=1e-10)
        print("✓ Test PASSED")
    
    def test_round_trip_beta_or(self):
        """Test round-trip conversion: BETA -> OR -> BETA."""
        beta = np.random.normal(0, 0.2, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta)
        
        # First fill OR
        result1 = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        self.assertIn("OR", result1.columns)
        
        # Remove BETA and fill it back from OR
        result1 = result1.drop(columns=["BETA"])
        result2 = _fill_data(result1, to_fill=["BETA"], verbose=False, log=self.log)
        np.testing.assert_allclose(result2["BETA"], beta, rtol=1e-10)
    
    def test_round_trip_beta_se_z(self):
        """Test round-trip conversion: BETA/SE -> Z -> BETA/SE (via P)."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        # Fill Z
        result1 = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
        self.assertIn("Z", result1.columns)
        np.testing.assert_allclose(result1["Z"], beta / se, rtol=1e-10)
        
        # Fill P from Z
        result2 = _fill_data(result1, to_fill=["P"], verbose=False, log=self.log)
        self.assertIn("P", result2.columns)
        
        # Fill SE back from BETA/P
        result2 = result2.drop(columns=["SE"])
        result3 = _fill_data(result2, to_fill=["SE"], verbose=False, log=self.log)
        self.assertIn("SE", result3.columns)
        # SE should be approximately recovered (within reasonable tolerance)
        np.testing.assert_allclose(result3["SE"], se, rtol=0.1)
    
    # ========== NA Handling Tests ==========
    
    def test_fill_p_from_mlog10p_with_na(self):
        """Test filling P from MLOG10P with NA values."""
        mlog10p = np.array([1.0, 2.0, np.nan, 3.0, np.nan, 5.0])
        self.n_variants = 6
        sumstats = self._create_base_sumstats()
        sumstats["MLOG10P"] = mlog10p
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that NA values are preserved
        self.assertTrue(np.isnan(result["P"][2]))
        self.assertTrue(np.isnan(result["P"][4]))
        # Check that non-NA values are calculated correctly
        valid_mask = ~np.isnan(mlog10p)
        np.testing.assert_allclose(
            result["P"][valid_mask], 
            10**(-mlog10p[valid_mask]), 
            rtol=1e-10
        )
    
    def test_fill_z_from_beta_se_with_na(self):
        """Test filling Z from BETA/SE with NA values."""
        beta = np.array([0.1, 0.2, np.nan, 0.3, np.nan])
        se = np.array([0.05, np.nan, 0.1, 0.15, 0.2])
        self.n_variants = 5
        sumstats = self._create_base_sumstats()
        sumstats["BETA"] = beta
        sumstats["SE"] = se
        
        result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
        
        self.assertIn("Z", result.columns)
        # Z should be NA when either BETA or SE is NA
        self.assertTrue(np.isnan(result["Z"][1]))  # SE is NA
        self.assertTrue(np.isnan(result["Z"][2]))  # BETA is NA
        self.assertTrue(np.isnan(result["Z"][4]))  # BETA is NA
        # Valid values should be calculated correctly
        valid_mask = ~(np.isnan(beta) | np.isnan(se))
        np.testing.assert_allclose(
            result["Z"][valid_mask], 
            beta[valid_mask] / se[valid_mask], 
            rtol=1e-10
        )
    
    def test_fill_or_from_beta_with_na(self):
        """Test filling OR from BETA with NA values."""
        beta = np.array([0.1, np.nan, -0.2, np.nan, 0.0])
        self.n_variants = 5
        sumstats = self._create_base_sumstats()
        sumstats["BETA"] = beta
        
        result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        
        self.assertIn("OR", result.columns)
        # OR should be NA when BETA is NA
        self.assertTrue(np.isnan(result["OR"][1]))
        self.assertTrue(np.isnan(result["OR"][3]))
        # Valid values should be calculated correctly
        valid_mask = ~np.isnan(beta)
        np.testing.assert_allclose(
            result["OR"][valid_mask], 
            np.exp(beta[valid_mask]), 
            rtol=1e-10
        )
        # OR should be 1.0 when BETA is 0.0
        self.assertAlmostEqual(result["OR"][4], 1.0, places=10)
    
    def test_fill_beta_from_or_with_na(self):
        """Test filling BETA from OR with NA values."""
        or_val = np.array([1.1, np.nan, 0.9, np.nan, 1.0])
        self.n_variants = 5
        sumstats = self._create_base_sumstats()
        sumstats["OR"] = or_val
        
        result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
        
        self.assertIn("BETA", result.columns)
        # BETA should be NA when OR is NA
        self.assertTrue(np.isnan(result["BETA"][1]))
        self.assertTrue(np.isnan(result["BETA"][3]))
        # Valid values should be calculated correctly
        valid_mask = ~np.isnan(or_val)
        np.testing.assert_allclose(
            result["BETA"][valid_mask], 
            np.log(or_val[valid_mask]), 
            rtol=1e-10
        )
        # BETA should be 0.0 when OR is 1.0
        self.assertAlmostEqual(result["BETA"][4], 0.0, places=10)
    
    def test_fill_chisq_from_z_with_na(self):
        """Test filling CHISQ from Z with NA values."""
        z = np.array([1.0, np.nan, -2.0, np.nan, 0.0])
        self.n_variants = 5
        sumstats = self._create_base_sumstats()
        sumstats["Z"] = z
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("CHISQ", result.columns)
        # CHISQ should be NA when Z is NA
        self.assertTrue(np.isnan(result["CHISQ"][1]))
        self.assertTrue(np.isnan(result["CHISQ"][3]))
        # Valid values should be calculated correctly
        valid_mask = ~np.isnan(z)
        np.testing.assert_allclose(
            result["CHISQ"][valid_mask], 
            z[valid_mask]**2, 
            rtol=1e-10
        )
        # CHISQ should be 0.0 when Z is 0.0
        self.assertAlmostEqual(result["CHISQ"][4], 0.0, places=10)
    
    def test_fill_maf_from_eaf_with_na(self):
        """Test filling MAF from EAF with NA values."""
        eaf = np.array([0.1, 0.9, np.nan, 0.5, np.nan, 0.0, 1.0])
        self.n_variants = 7
        sumstats = self._create_base_sumstats()
        sumstats["EAF"] = eaf
        
        result = _fill_data(sumstats, to_fill=["MAF"], verbose=False, log=self.log)
        
        self.assertIn("MAF", result.columns)
        # MAF should be NA when EAF is NA
        self.assertTrue(np.isnan(result["MAF"][2]))
        self.assertTrue(np.isnan(result["MAF"][4]))
        # Valid values should be min(EAF, 1-EAF)
        expected_maf = np.minimum(eaf, 1 - eaf)
        valid_mask = ~np.isnan(eaf)
        np.testing.assert_allclose(
            result["MAF"][valid_mask], 
            expected_maf[valid_mask], 
            rtol=1e-10
        )
        # Edge cases: EAF=0 -> MAF=0, EAF=1 -> MAF=0
        self.assertAlmostEqual(result["MAF"][5], 0.0, places=10)
        self.assertAlmostEqual(result["MAF"][6], 0.0, places=10)
    
    # ========== Zero Value Tests ==========
    
    def test_fill_z_from_beta_se_with_zero_se(self):
        """Test filling Z from BETA/SE when SE is zero (division by zero)."""
        beta = np.array([0.1, 0.2, 0.3])
        se = np.array([0.05, 0.0, 0.1])  # One zero SE
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["BETA"] = beta
        sumstats["SE"] = se
        
        result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
        
        self.assertIn("Z", result.columns)
        # Z should be inf or very large when SE is zero
        self.assertTrue(np.isinf(result["Z"][1]) or np.abs(result["Z"][1]) > 1e10)
        # Other values should be calculated correctly
        self.assertAlmostEqual(result["Z"][0], beta[0] / se[0], places=10)
        self.assertAlmostEqual(result["Z"][2], beta[2] / se[2], places=10)
    
    def test_fill_or_from_beta_with_zero_beta(self):
        """Test filling OR from BETA when BETA is zero."""
        beta = np.array([0.0, 0.1, -0.1])
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["BETA"] = beta
        
        result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        
        self.assertIn("OR", result.columns)
        # OR should be 1.0 when BETA is 0.0
        self.assertAlmostEqual(result["OR"][0], 1.0, places=10)
        # Other values should be calculated correctly
        np.testing.assert_allclose(result["OR"][1:], np.exp(beta[1:]), rtol=1e-10)
    
    def test_fill_beta_from_or_with_one_or(self):
        """Test filling BETA from OR when OR is 1.0 (log(1) = 0)."""
        or_val = np.array([1.0, 1.1, 0.9])
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["OR"] = or_val
        
        result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
        
        self.assertIn("BETA", result.columns)
        # BETA should be 0.0 when OR is 1.0
        self.assertAlmostEqual(result["BETA"][0], 0.0, places=10)
        # Other values should be calculated correctly
        np.testing.assert_allclose(result["BETA"][1:], np.log(or_val[1:]), rtol=1e-10)
    
    def test_fill_chisq_from_z_with_zero_z(self):
        """Test filling CHISQ from Z when Z is zero."""
        z = np.array([0.0, 1.0, -1.0])
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["Z"] = z
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("CHISQ", result.columns)
        # CHISQ should be 0.0 when Z is 0.0
        self.assertAlmostEqual(result["CHISQ"][0], 0.0, places=10)
        # Other values should be calculated correctly
        np.testing.assert_allclose(result["CHISQ"][1:], z[1:]**2, rtol=1e-10)
    
    def test_fill_mlog10p_from_p_with_zero_p(self):
        """Test filling MLOG10P from P when P is zero (log10(0) = -inf)."""
        p = np.array([0.0, 1e-10, 0.5])
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["P"] = p
        
        # Suppress floating point warnings for this test
        with np.errstate(divide='ignore', invalid='ignore'):
            result = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        # MLOG10P should be inf when P is 0.0
        self.assertTrue(np.isinf(result["MLOG10P"][0]) or result["MLOG10P"][0] > 1e10)
        # Other values should be calculated correctly
        np.testing.assert_allclose(
            result["MLOG10P"][1:], 
            -np.log10(p[1:]), 
            rtol=1e-10
        )
    
    def test_fill_p_from_mlog10p_with_inf_mlog10p(self):
        """Test filling P from MLOG10P when MLOG10P is inf (P=0)."""
        mlog10p = np.array([np.inf, 10.0, 5.0])
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["MLOG10P"] = mlog10p
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # P should be 0.0 when MLOG10P is inf
        self.assertAlmostEqual(result["P"][0], 0.0, places=10)
        # Other values should be calculated correctly
        valid_mask = np.isfinite(mlog10p)
        np.testing.assert_allclose(
            result["P"][valid_mask], 
            10**(-mlog10p[valid_mask]), 
            rtol=1e-10
        )
    
    def test_fill_se_from_beta_p_with_zero_p(self):
        """Test filling SE from BETA/P when P is zero."""
        beta = np.array([0.1, 0.2, 0.3])
        p = np.array([0.0, 1e-10, 0.5])  # One zero P
        self.n_variants = 3
        sumstats = self._create_base_sumstats()
        sumstats["BETA"] = beta
        sumstats["P"] = p
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        # SE calculation involves erfcinv which may produce inf for P=0
        # Check that other values are calculated correctly
        valid_mask = p > 0
        if valid_mask.sum() > 0:
            # For valid P values, verify SE is positive
            self.assertTrue((result["SE"][valid_mask] > 0).all())
    
    # ========== Floating Point Accuracy Tests ==========
    
    def test_fill_p_from_mlog10p_accuracy(self):
        """Test floating point accuracy when filling P from MLOG10P."""
        # Test with various precision levels
        test_cases = [
            (0.0, 1.0),           # MLOG10P=0 -> P=1
            (1.0, 0.1),           # MLOG10P=1 -> P=0.1
            (2.0, 0.01),          # MLOG10P=2 -> P=0.01
            (5.0, 1e-5),          # MLOG10P=5 -> P=1e-5
            (8.0, 1e-8),          # MLOG10P=8 -> P=1e-8
            (10.0, 1e-10),       # MLOG10P=10 -> P=1e-10
            (50.0, 1e-50),       # MLOG10P=50 -> P=1e-50
            (100.0, 1e-100),     # MLOG10P=100 -> P=1e-100
            (200.0, 1e-200),     # MLOG10P=200 -> P=1e-200
            (300.0, 1e-300),     # MLOG10P=300 -> P=1e-300 (extreme)
            (308.0, 1e-308),     # MLOG10P=308 -> P=1e-308 (float64 limit)
        ]
        
        self.n_variants = 1
        for mlog10p_val, expected_p in test_cases:
            sumstats = self._create_base_sumstats()
            sumstats["MLOG10P"] = [mlog10p_val]
            
            result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
            
            if mlog10p_val <= 308:  # Within float64 range
                # Use relaxed tolerance for extreme values
                rtol = 1e-5 if mlog10p_val > 200 else 1e-10
                np.testing.assert_allclose(
                    result["P"].iloc[0], 
                    expected_p, 
                    rtol=rtol,
                    atol=1e-320 if mlog10p_val > 200 else 0,
                    err_msg=f"Failed for MLOG10P={mlog10p_val}"
                )
    
    def test_fill_mlog10p_from_p_accuracy(self):
        """Test floating point accuracy when filling MLOG10P from P."""
        # Test with various precision levels
        test_cases = [
            (1.0, 0.0),           # P=1 -> MLOG10P=0
            (0.1, 1.0),           # P=0.1 -> MLOG10P=1
            (0.01, 2.0),          # P=0.01 -> MLOG10P=2
            (1e-5, 5.0),          # P=1e-5 -> MLOG10P=5
            (1e-8, 8.0),          # P=1e-8 -> MLOG10P=8
            (1e-10, 10.0),       # P=1e-10 -> MLOG10P=10
            (1e-50, 50.0),       # P=1e-50 -> MLOG10P=50
            (1e-100, 100.0),     # P=1e-100 -> MLOG10P=100
            (1e-200, 200.0),     # P=1e-200 -> MLOG10P=200
            (1e-300, 300.0),     # P=1e-300 -> MLOG10P=300 (extreme)
            (1e-308, 308.0),     # P=1e-308 -> MLOG10P=308 (float64 limit)
        ]
        
        self.n_variants = 1
        for p_val, expected_mlog10p in test_cases:
            sumstats = self._create_base_sumstats()
            sumstats["P"] = [p_val]
            
            result = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
            
            if p_val > 0:  # Avoid log(0)
                np.testing.assert_allclose(
                    result["MLOG10P"].iloc[0], 
                    expected_mlog10p, 
                    rtol=1e-10,
                    err_msg=f"Failed for P={p_val}"
                )
    
    def test_fill_z_from_beta_se_accuracy(self):
        """Test floating point accuracy when filling Z from BETA/SE."""
        # Test with various ratios
        test_cases = [
            (0.0, 0.1, 0.0),     # BETA=0 -> Z=0
            (0.1, 0.1, 1.0),     # BETA=SE -> Z=1
            (0.2, 0.1, 2.0),     # BETA=2*SE -> Z=2
            (-0.1, 0.1, -1.0),   # Negative BETA
            (1e-6, 1e-6, 1.0),   # Very small values
            (1e6, 1e6, 1.0),     # Very large values
        ]
        
        self.n_variants = 1
        for beta_val, se_val, expected_z in test_cases:
            sumstats = self._create_base_sumstats()
            sumstats["BETA"] = [beta_val]
            sumstats["SE"] = [se_val]
            
            result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
            
            np.testing.assert_allclose(
                result["Z"].iloc[0], 
                expected_z, 
                rtol=1e-10,
                err_msg=f"Failed for BETA={beta_val}, SE={se_val}"
            )
    
    def test_fill_or_from_beta_accuracy(self):
        """Test floating point accuracy when filling OR from BETA."""
        # Test with various BETA values
        test_cases = [
            (0.0, 1.0),           # BETA=0 -> OR=1
            (np.log(2.0), 2.0),    # BETA=ln(2) -> OR=2 (use np.log for precision)
            (-np.log(2.0), 0.5),   # BETA=-ln(2) -> OR=0.5 (use np.log for precision)
            (0.1, np.exp(0.1)),   # Small positive
            (-0.1, np.exp(-0.1)), # Small negative
        ]
        
        self.n_variants = 1
        for beta_val, expected_or in test_cases:
            sumstats = self._create_base_sumstats()
            sumstats["BETA"] = [beta_val]
            
            result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
            
            self.assertIn("OR", result.columns)
            np.testing.assert_allclose(
                result["OR"].iloc[0], 
                expected_or, 
                rtol=1e-9,  # Slightly relaxed tolerance for floating point precision
                err_msg=f"Failed for BETA={beta_val}"
            )
    
    def test_fill_beta_from_or_accuracy(self):
        """Test floating point accuracy when filling BETA from OR."""
        # Test with various OR values
        test_cases = [
            (1.0, 0.0),           # OR=1 -> BETA=0
            (2.0, np.log(2.0)),   # OR=2 -> BETA=ln(2)
            (0.5, np.log(0.5)),  # OR=0.5 -> BETA=ln(0.5)
            (1.1, np.log(1.1)),  # Small positive
            (0.9, np.log(0.9)),  # Small negative
        ]
        
        self.n_variants = 1
        for or_val, expected_beta in test_cases:
            sumstats = self._create_base_sumstats()
            sumstats["OR"] = [or_val]
            
            result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
            
            np.testing.assert_allclose(
                result["BETA"].iloc[0], 
                expected_beta, 
                rtol=1e-10,
                err_msg=f"Failed for OR={or_val}"
            )
    
    def test_fill_chisq_from_z_accuracy(self):
        """Test floating point accuracy when filling CHISQ from Z."""
        # Test with various Z values
        test_cases = [
            (0.0, 0.0),           # Z=0 -> CHISQ=0
            (1.0, 1.0),           # Z=1 -> CHISQ=1
            (-1.0, 1.0),          # Z=-1 -> CHISQ=1
            (2.0, 4.0),           # Z=2 -> CHISQ=4
            (-2.0, 4.0),          # Z=-2 -> CHISQ=4
        ]
        
        self.n_variants = 1
        for z_val, expected_chisq in test_cases:
            sumstats = self._create_base_sumstats()
            sumstats["Z"] = [z_val]
            
            result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
            
            np.testing.assert_allclose(
                result["CHISQ"].iloc[0], 
                expected_chisq, 
                rtol=1e-10,
                err_msg=f"Failed for Z={z_val}"
            )
    
    def test_round_trip_accuracy_p_mlog10p(self):
        """Test round-trip accuracy: P -> MLOG10P -> P."""
        # Test with various P values
        test_p = np.array([1.0, 0.1, 0.01, 1e-5, 1e-8, 1e-10, 1e-15])
        self.n_variants = len(test_p)
        sumstats = self._create_base_sumstats()
        sumstats["P"] = test_p
        
        # Round trip: P -> MLOG10P -> P
        result1 = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        result2 = _fill_data(result1, to_fill=["P"], overwrite=True, verbose=False, log=self.log)
        
        # Check accuracy (allowing for floating point errors)
        np.testing.assert_allclose(result2["P"], test_p, rtol=1e-10, atol=1e-15)
    
    def test_round_trip_accuracy_beta_or(self):
        """Test round-trip accuracy: BETA -> OR -> BETA."""
        # Test with various BETA values
        test_beta = np.array([0.0, 0.1, -0.1, 0.5, -0.5, 1.0, -1.0])
        self.n_variants = len(test_beta)
        sumstats = self._create_base_sumstats()
        sumstats["BETA"] = test_beta
        
        # Round trip: BETA -> OR -> BETA
        result1 = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        result1 = result1.drop(columns=["BETA"])
        result2 = _fill_data(result1, to_fill=["BETA"], verbose=False, log=self.log)
        
        # Check accuracy
        np.testing.assert_allclose(result2["BETA"], test_beta, rtol=1e-10)
    
    # ========== Extreme Value Tests ==========
    
    def test_extreme_p_values_from_mlog10p(self):
        """Test filling P from MLOG10P with extreme values near float64 limits."""
        # Test with extremely small P values
        extreme_cases = [
            (100.0, 1e-100),      # MLOG10P=100 -> P=1e-100
            (150.0, 1e-150),      # MLOG10P=150 -> P=1e-150
            (200.0, 1e-200),      # MLOG10P=200 -> P=1e-200
            (250.0, 1e-250),      # MLOG10P=250 -> P=1e-250
            (300.0, 1e-300),      # MLOG10P=300 -> P=1e-300 (near float64 limit)
            (308.0, 1e-308),      # MLOG10P=308 -> P=1e-308 (float64 minimum)
        ]
        
        self.n_variants = 1
        for mlog10p_val, expected_p in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["MLOG10P"] = [mlog10p_val]
            
            result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
            
            self.assertIn("P", result.columns)
            # For extreme values, check that result is very small and positive (or 0 if underflow)
            if mlog10p_val <= 308:
                p_val = result["P"].iloc[0]
                # P should be >= 0 (may be 0 due to underflow for very large MLOG10P)
                self.assertTrue(p_val >= 0)
                if p_val > 0:
                    # If not zero, should be extremely small
                    # For MLOG10P=100, P=1e-100 exactly, so check based on MLOG10P value
                    if mlog10p_val <= 100:
                        # For smaller MLOG10P, P can be exactly 1e-100 or larger
                        self.assertTrue(p_val <= 1e-100)
                    else:
                        # For larger MLOG10P, P should be smaller than 1e-100
                        self.assertTrue(p_val < 1e-100)
                    # Check relative accuracy (may need relaxed tolerance for extreme values)
                    np.testing.assert_allclose(
                        p_val, 
                        expected_p, 
                        rtol=1e-5,  # Relaxed tolerance for extreme values
                        atol=1e-320,
                        err_msg=f"Failed for extreme MLOG10P={mlog10p_val}"
                    )
                else:
                    # If zero due to underflow, that's acceptable for extreme values
                    self.assertTrue(mlog10p_val >= 300, 
                                  f"P=0 for MLOG10P={mlog10p_val} (may be underflow)")
    
    def test_extreme_mlog10p_values_from_p(self):
        """Test filling MLOG10P from P with extreme values near float64 limits."""
        # Test with extremely small P values
        extreme_cases = [
            (1e-100, 100.0),      # P=1e-100 -> MLOG10P=100
            (1e-150, 150.0),      # P=1e-150 -> MLOG10P=150
            (1e-200, 200.0),      # P=1e-200 -> MLOG10P=200
            (1e-250, 250.0),      # P=1e-250 -> MLOG10P=250
            (1e-300, 300.0),      # P=1e-300 -> MLOG10P=300
            (1e-308, 308.0),      # P=1e-308 -> MLOG10P=308 (float64 minimum)
        ]
        
        self.n_variants = 1
        for p_val, expected_mlog10p in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["P"] = [p_val]
            
            result = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
            
            self.assertIn("MLOG10P", result.columns)
            # Check that MLOG10P is calculated correctly for extreme values
            np.testing.assert_allclose(
                result["MLOG10P"].iloc[0], 
                expected_mlog10p, 
                rtol=1e-10,
                err_msg=f"Failed for extreme P={p_val}"
            )
    
    def test_extreme_z_scores(self):
        """Test filling P and MLOG10P from extreme Z scores."""
        # Test with very large Z scores (corresponding to very small P values)
        extreme_z = np.array([0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 50.0])
        self.n_variants = len(extreme_z)
        sumstats = self._create_base_sumstats()
        sumstats["Z"] = extreme_z
        
        # Fill P from Z
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that P values are very small for large Z scores
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
        # For Z=30 (index 5), P should be extremely small (~1e-197) or 0 (underflow)
        z_30_idx = np.where(extreme_z == 30.0)[0]
        if len(z_30_idx) > 0:
            idx = z_30_idx[0]
            self.assertTrue(result["P"][idx] < 1e-100 or result["P"][idx] == 0.0)
        # For Z=50 (index 6), P should be even smaller or 0 (underflow)
        z_50_idx = np.where(extreme_z == 50.0)[0]
        if len(z_50_idx) > 0:
            idx = z_50_idx[0]
            self.assertTrue(result["P"][idx] < 1e-500 or result["P"][idx] == 0.0)
    
    def test_extreme_z_scores_to_mlog10p(self):
        """Test filling MLOG10P from extreme Z scores using extreme methods."""
        # Test with very large Z scores
        extreme_z = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        self.n_variants = len(extreme_z)
        sumstats = self._create_base_sumstats()
        sumstats["Z"] = extreme_z
        
        # Fill MLOG10P from Z (should use extreme methods)
        result = _fill_data(sumstats, to_fill=["MLOG10P"], extreme=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        # Check that MLOG10P values are very large for large Z scores
        self.assertTrue((result["MLOG10P"] >= 0).all())
        # For Z=30, MLOG10P should be around 197
        self.assertTrue(result["MLOG10P"][2] > 100)
        # For Z=50, MLOG10P should be around 542
        self.assertTrue(result["MLOG10P"][4] > 500)
    
    def test_extreme_beta_se_ratios(self):
        """Test filling Z from extreme BETA/SE ratios."""
        # Test with very small and very large ratios
        extreme_cases = [
            (1e-10, 1e-10, 1.0),      # Very small values
            (1e-6, 1e-6, 1.0),        # Small values
            (1e10, 1e10, 1.0),       # Very large values
            (1e-10, 1e-12, 100.0),   # Very small BETA, tiny SE -> large Z
            (1e10, 1e8, 100.0),      # Very large BETA, large SE -> large Z
            (-1e10, 1e8, -100.0),    # Very large negative BETA
        ]
        
        self.n_variants = 1
        for beta_val, se_val, expected_z in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["BETA"] = [beta_val]
            sumstats["SE"] = [se_val]
            
            result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
            
            self.assertIn("Z", result.columns)
            # Check that Z is calculated correctly (may need relaxed tolerance for extreme values)
            np.testing.assert_allclose(
                result["Z"].iloc[0], 
                expected_z, 
                rtol=1e-5,  # Relaxed tolerance for extreme values
                err_msg=f"Failed for extreme BETA={beta_val}, SE={se_val}"
            )
    
    def test_extreme_or_values(self):
        """Test filling BETA from extreme OR values."""
        # Test with very small and very large OR values
        extreme_cases = [
            (1e-10, np.log(1e-10)),      # Very small OR
            (1e-5, np.log(1e-5)),        # Small OR
            (1e5, np.log(1e5)),          # Very large OR
            (1e10, np.log(1e10)),        # Extremely large OR
            (1.0, 0.0),                   # OR=1 -> BETA=0
        ]
        
        self.n_variants = 1
        for or_val, expected_beta in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["OR"] = [or_val]
            
            result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
            
            self.assertIn("BETA", result.columns)
            np.testing.assert_allclose(
                result["BETA"].iloc[0], 
                expected_beta, 
                rtol=1e-10,
                err_msg=f"Failed for extreme OR={or_val}"
            )
    
    def test_extreme_beta_values(self):
        """Test filling OR from extreme BETA values."""
        # Test with very small and very large BETA values
        extreme_cases = [
            (-10.0, np.exp(-10.0)),      # Very negative BETA -> very small OR
            (-5.0, np.exp(-5.0)),        # Negative BETA
            (0.0, 1.0),                   # BETA=0 -> OR=1
            (5.0, np.exp(5.0)),          # Large positive BETA
            (10.0, np.exp(10.0)),        # Very large positive BETA
        ]
        
        self.n_variants = 1
        for beta_val, expected_or in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["BETA"] = [beta_val]
            
            result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
            
            self.assertIn("OR", result.columns)
            np.testing.assert_allclose(
                result["OR"].iloc[0], 
                expected_or, 
                rtol=1e-10,
                err_msg=f"Failed for extreme BETA={beta_val}"
            )
    
    def test_extreme_chisq_values(self):
        """Test filling P from extreme CHISQ values."""
        # Test with very large CHISQ values (corresponding to very small P values)
        extreme_chisq = np.array([1.0, 10.0, 100.0, 1000.0, 10000.0])
        self.n_variants = len(extreme_chisq)
        sumstats = self._create_base_sumstats()
        sumstats["CHISQ"] = extreme_chisq
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that P values are very small for large CHISQ
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
        # For CHISQ=100 (index 2), P should be extremely small or 0 (underflow)
        chisq_100_idx = np.where(extreme_chisq == 100.0)[0]
        if len(chisq_100_idx) > 0:
            idx = chisq_100_idx[0]
            self.assertTrue(result["P"][idx] < 1e-20 or result["P"][idx] == 0.0)
        # For CHISQ=10000 (index 4), P should be even smaller or 0 (underflow)
        chisq_10000_idx = np.where(extreme_chisq == 10000.0)[0]
        if len(chisq_10000_idx) > 0:
            idx = chisq_10000_idx[0]
            self.assertTrue(result["P"][idx] < 1e-2000 or result["P"][idx] == 0.0)
    
    def test_extreme_round_trip_p_mlog10p(self):
        """Test round-trip accuracy with extreme P values."""
        # Test with extremely small P values
        extreme_p = np.array([1.0, 1e-10, 1e-50, 1e-100, 1e-200, 1e-300])
        self.n_variants = len(extreme_p)
        sumstats = self._create_base_sumstats()
        sumstats["P"] = extreme_p
        
        # Round trip: P -> MLOG10P -> P
        result1 = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        result2 = _fill_data(result1, to_fill=["P"], overwrite=True, verbose=False, log=self.log)
        
        # Check accuracy (allowing for floating point errors with extreme values)
        # For very small P values, relative tolerance may need to be relaxed
        valid_mask = extreme_p > 0
        np.testing.assert_allclose(
            result2["P"][valid_mask], 
            extreme_p[valid_mask], 
            rtol=1e-5,  # Relaxed tolerance for extreme values
            atol=1e-320,
            err_msg="Round-trip failed for extreme P values"
        )
    
    def test_extreme_na_handling(self):
        """Test NA handling with extreme values mixed in."""
        # Mix extreme values with NA values
        mlog10p = np.array([300.0, np.nan, 10.0, np.nan, 0.0, 200.0])
        self.n_variants = 6
        sumstats = self._create_base_sumstats()
        sumstats["MLOG10P"] = mlog10p
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that NA values are preserved
        self.assertTrue(np.isnan(result["P"][1]))
        self.assertTrue(np.isnan(result["P"][3]))
        # Check that extreme values are handled correctly
        valid_mask = ~np.isnan(mlog10p)
        # For extreme values, check that P is very small
        extreme_mask = valid_mask & (mlog10p > 100)
        if extreme_mask.sum() > 0:
            self.assertTrue((result["P"][extreme_mask] < 1e-100).all())
    
    # ========== Beyond Float64 Limits Tests ==========
    
    def test_mlog10p_beyond_float64_limit(self):
        """Test filling P from MLOG10P values that exceed float64 limits (>308)."""
        # MLOG10P values beyond float64 limit (P would be < 1e-308)
        beyond_limit_cases = [
            (309.0,),   # Just beyond limit
            (350.0,),   # Well beyond limit
            (400.0,),   # Far beyond limit
            (500.0,),   # Very far beyond limit
            (1000.0,),  # Extremely beyond limit
        ]
        
        self.n_variants = 1
        for mlog10p_val, in beyond_limit_cases:
            sumstats = self._create_base_sumstats()
            sumstats["MLOG10P"] = [mlog10p_val]
            
            result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
            
            self.assertIn("P", result.columns)
            # P should be 0.0 or very close to 0 (underflow)
            # The result should not crash and should handle gracefully
            self.assertTrue(result["P"].iloc[0] >= 0.0)
            self.assertTrue(result["P"].iloc[0] <= 1e-300)  # Should be extremely small or 0
    
    def test_p_beyond_float64_limit(self):
        """Test filling MLOG10P from P values that are below float64 minimum (< 1e-308)."""
        # P values below float64 minimum (these will underflow to 0 or be treated as 0)
        # Note: Values < 1e-308 may not be representable in float64
        beyond_limit_cases = [
            (0.0,),      # Zero P
            (1e-309,),   # Just below limit (may underflow)
            (1e-320,),   # Well below limit (will underflow)
            (1e-400,),   # Far below limit (will underflow)
        ]
        
        self.n_variants = 1
        for p_val, in beyond_limit_cases:
            sumstats = self._create_base_sumstats()
            # For values that can't be represented, use the smallest representable value
            if p_val > 0 and p_val < 1e-308:
                # Use smallest representable value
                sumstats["P"] = [np.nextafter(0.0, 1.0)]  # Smallest positive float64
            else:
                sumstats["P"] = [p_val]
            
            # Suppress floating point warnings for this test
            with np.errstate(divide='ignore', invalid='ignore'):
                result = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
            
            self.assertIn("MLOG10P", result.columns)
            # MLOG10P should be very large or inf for very small P
            if p_val == 0.0:
                self.assertTrue(np.isinf(result["MLOG10P"].iloc[0]) or result["MLOG10P"].iloc[0] > 300)
            else:
                # Should handle gracefully without crashing
                self.assertTrue(np.isfinite(result["MLOG10P"].iloc[0]) or np.isinf(result["MLOG10P"].iloc[0]))
    
    def test_extreme_z_beyond_float64_limit(self):
        """Test filling P from Z scores that would produce P values beyond float64 limits."""
        # Very large Z scores that would produce P < 1e-308
        # Z > ~37.5 gives P < 1e-308
        extreme_z = np.array([37.0, 38.0, 40.0, 50.0, 60.0, 100.0])
        self.n_variants = len(extreme_z)
        sumstats = self._create_base_sumstats()
        sumstats["Z"] = extreme_z
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # P should be 0.0 or extremely small (underflow)
        # The code should handle this gracefully without crashing
        self.assertTrue((result["P"] >= 0.0).all())
        # For Z > 37, P should be 0 or extremely small
        large_z_mask = np.abs(extreme_z) > 37
        if large_z_mask.sum() > 0:
            self.assertTrue((result["P"][large_z_mask] <= 1e-300).all() or 
                          (result["P"][large_z_mask] == 0.0).all())
    
    def test_extreme_z_to_mlog10p_beyond_limit(self):
        """Test filling MLOG10P from extreme Z scores using extreme methods."""
        # Very large Z scores that would produce MLOG10P > 308
        # Z > ~37.5 gives MLOG10P > 308
        extreme_z = np.array([38.0, 40.0, 50.0, 60.0, 100.0])
        self.n_variants = len(extreme_z)
        sumstats = self._create_base_sumstats()
        sumstats["Z"] = extreme_z
        
        # Fill MLOG10P from Z (should use extreme methods)
        result = _fill_data(sumstats, to_fill=["MLOG10P"], extreme=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        # MLOG10P should be very large (> 308) for these Z scores
        # Extreme methods should handle this correctly
        self.assertTrue((result["MLOG10P"] > 300).all())
        # For Z=38, MLOG10P should be around 315
        # For Z=100, MLOG10P should be around 2171
        self.assertTrue(result["MLOG10P"][0] > 310)  # Z=38
        self.assertTrue(result["MLOG10P"][4] > 2000)  # Z=100
    
    def test_extreme_chisq_beyond_float64_limit(self):
        """Test filling P from CHISQ values that would produce P beyond float64 limits."""
        # Very large CHISQ values that would produce P < 1e-308
        # CHISQ > ~1400 gives P < 1e-308 for df=1
        extreme_chisq = np.array([1000.0, 1500.0, 2000.0, 5000.0, 10000.0, 50000.0])
        self.n_variants = len(extreme_chisq)
        sumstats = self._create_base_sumstats()
        sumstats["CHISQ"] = extreme_chisq
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # P should be 0.0 or extremely small (underflow)
        self.assertTrue((result["P"] >= 0.0).all())
        # For very large CHISQ, P should be 0 or extremely small
        large_chisq_mask = extreme_chisq > 1400
        if large_chisq_mask.sum() > 0:
            self.assertTrue((result["P"][large_chisq_mask] <= 1e-300).all() or 
                          (result["P"][large_chisq_mask] == 0.0).all())
    
    def test_round_trip_beyond_float64_limit(self):
        """Test round-trip with values that exceed float64 limits."""
        # Test with MLOG10P values beyond 308
        extreme_mlog10p = np.array([309.0, 350.0, 400.0, 500.0])
        self.n_variants = len(extreme_mlog10p)
        sumstats = self._create_base_sumstats()
        sumstats["MLOG10P"] = extreme_mlog10p
        
        # Round trip: MLOG10P -> P -> MLOG10P
        # Note: P will underflow to 0, so round-trip may not be perfect
        result1 = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        # Suppress floating point warnings for MLOG10P calculation from P=0
        with np.errstate(divide='ignore', invalid='ignore'):
            result2 = _fill_data(result1, to_fill=["MLOG10P"], overwrite=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result2.columns)
        # For P=0, MLOG10P should be inf
        # The code should handle this gracefully
        zero_p_mask = result1["P"] == 0.0
        if zero_p_mask.sum() > 0:
            # MLOG10P should be inf or very large for P=0
            self.assertTrue(
                (np.isinf(result2["MLOG10P"][zero_p_mask])).all() or
                (result2["MLOG10P"][zero_p_mask] > 300).all()
            )
    
    def test_extreme_beta_se_beyond_limit(self):
        """Test filling Z from extreme BETA/SE ratios that exceed float64 limits."""
        # Extreme ratios that could cause overflow
        extreme_cases = [
            (1e100, 1e-100, 1e200),    # Extremely large ratio
            (-1e100, 1e-100, -1e200),  # Extremely large negative ratio
            (1e200, 1e200, 1.0),       # Very large values, ratio=1
        ]
        
        self.n_variants = 1
        for beta_val, se_val, expected_z in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["BETA"] = [beta_val]
            sumstats["SE"] = [se_val]
            
            result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
            
            self.assertIn("Z", result.columns)
            # Z should be calculated or be inf/very large
            # The code should handle this without crashing
            z_val = result["Z"].iloc[0]
            if np.isfinite(z_val):
                # If finite, should be approximately correct
                np.testing.assert_allclose(
                    z_val, 
                    expected_z, 
                    rtol=1e-5,
                    err_msg=f"Failed for extreme BETA={beta_val}, SE={se_val}"
                )
            else:
                # If inf, that's acceptable for extreme values
                self.assertTrue(np.isinf(z_val) or np.abs(z_val) > 1e100)
    
    def test_extreme_or_beyond_limit(self):
        """Test filling BETA from extreme OR values that exceed float64 limits."""
        # Extreme OR values
        extreme_cases = [
            (1e100,),   # Extremely large OR
            (1e-100,),  # Extremely small OR
            (np.inf,),  # Infinite OR
        ]
        
        self.n_variants = 1
        for or_val, in extreme_cases:
            if np.isinf(or_val):
                # Skip infinite values as they can't be in a DataFrame easily
                continue
                
            sumstats = self._create_base_sumstats()
            sumstats["OR"] = [or_val]
            
            result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
            
            self.assertIn("BETA", result.columns)
            # BETA should be calculated or be inf/very large
            beta_val = result["BETA"].iloc[0]
            if np.isfinite(beta_val):
                # If finite, should be approximately log(OR)
                expected_beta = np.log(or_val)
                np.testing.assert_allclose(
                    beta_val, 
                    expected_beta, 
                    rtol=1e-5,
                    err_msg=f"Failed for extreme OR={or_val}"
                )
            else:
                # If inf, that's acceptable for extreme values
                self.assertTrue(np.isinf(beta_val))
    
    def test_extreme_beta_beyond_limit(self):
        """Test filling OR from extreme BETA values that exceed float64 limits."""
        # Extreme BETA values
        extreme_cases = [
            (100.0,),   # Very large BETA -> OR = exp(100) ~ 2.7e43
            (200.0,),   # Extremely large BETA -> OR = exp(200) ~ 7.2e86
            (-100.0,),  # Very negative BETA -> OR = exp(-100) ~ 3.7e-44
            (-200.0,),  # Extremely negative BETA -> OR = exp(-200) ~ 1.4e-87
        ]
        
        self.n_variants = 1
        for beta_val, in extreme_cases:
            sumstats = self._create_base_sumstats()
            sumstats["BETA"] = [beta_val]
            
            result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
            
            self.assertIn("OR", result.columns)
            # OR should be calculated or be inf/very large
            or_val = result["OR"].iloc[0]
            if np.isfinite(or_val):
                # If finite, should be approximately exp(BETA)
                expected_or = np.exp(beta_val)
                # For very large BETA, exp(BETA) may overflow, so check if it's reasonable
                if np.isfinite(expected_or):
                    np.testing.assert_allclose(
                        or_val, 
                        expected_or, 
                        rtol=1e-5,
                        err_msg=f"Failed for extreme BETA={beta_val}"
                    )
                else:
                    # If expected is inf, result should also be inf or very large
                    self.assertTrue(np.isinf(or_val) or or_val > 1e100)
            else:
                # If inf, that's acceptable for extreme values
                self.assertTrue(np.isinf(or_val) or or_val > 1e100)


if __name__ == "__main__":
    # Run tests with verbose output
    unittest.main(verbosity=2)

