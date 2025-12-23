"""
Comprehensive test suite for SSF format validators.
"""
import os
import sys
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
import pytest

# Add src to path
ROOT = Path(__file__).parent.parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from gwaslab.extension.gwas_sumstats_tools.validate_ssf import validate_ssf_file
from gwaslab.info.g_Log import Log


def create_test_file(data: pd.DataFrame, base_path: Path, compress=False) -> Path:
    """Create test SSF file."""
    if compress:
        output = base_path.with_suffix('.tsv.gz')
        data.to_csv(output, sep='\t', index=False, compression='gzip')
    else:
        output = base_path.with_suffix('.tsv')
        data.to_csv(output, sep='\t', index=False)
    return output


class TestSSFValidator:
    """Test suite for SSF validator."""
    
    def test_valid_ssf_file_basic(self):
        """Test 1: Valid SSF file with all required columns in correct order."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "valid_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == True, f"Should pass validation: {message}"
            assert errors is None, "Should have no errors"
    
    def test_valid_ssf_file_with_odds_ratio(self):
        """Test 2: Valid SSF file with odds_ratio instead of beta."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'odds_ratio': np.random.uniform(0.1, 10, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "valid_or_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == True, f"Should pass validation: {message}"
    
    def test_wrong_column_order(self):
        """Test 3: Valid data but wrong column order - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'base_pair_location': np.random.randint(1, 250000000, n_rows),  # Wrong position
            'chromosome': np.random.randint(1, 23, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "wrong_order_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on wrong column order"
            assert "order" in message.lower() or (errors and any("order" in str(e).lower() for e in errors))
    
    def test_missing_required_column(self):
        """Test 4: Missing required column - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            # Missing 'other_allele'
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "missing_col_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on missing required column"
            assert "missing" in message.lower() or (errors and any("missing" in str(e).lower() for e in errors))
    
    def test_insufficient_rows(self):
        """Test 5: File with fewer than minimum rows - should fail."""
        n_rows = 50000  # Below minimum
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "insufficient_rows_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on insufficient rows"
            assert "minimum" in message.lower() or "rows" in message.lower()
    
    def test_invalid_chromosome_values(self):
        """Test 6: Invalid chromosome values (outside 1-25) - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': [1, 2, 3, 26, 27, 0, -1] + [np.random.randint(1, 23) for _ in range(n_rows - 7)],
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid_chr_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on invalid chromosome values"
            assert "chromosome" in message.lower() or (errors and any("chromosome" in str(e).lower() for e in errors))
    
    def test_invalid_p_value_range(self):
        """Test 7: Invalid p-value (outside [0,1]) - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': [1.5, 2.0, -0.1] + list(np.random.uniform(1e-10, 1, n_rows - 3)),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid_pval_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on invalid p-value"
            assert "p_value" in message.lower() or (errors and any("p_value" in str(e).lower() or "p-value" in str(e).lower() for e in errors))
    
    def test_very_small_p_value_scientific_notation(self):
        """Test 8: Very small p-values in scientific notation - should pass."""
        n_rows = 150000
        # Create data with string p-values in scientific notation
        p_values = ['1.23e-300', '5.67e-200'] + [str(x) for x in np.random.uniform(1e-10, 1, n_rows - 2)]
        
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': p_values,
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "small_pval_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            # Should pass because mantissa/exponent splitting handles this
            assert is_valid == True, f"Should pass with very small p-values: {message}"
    
    def test_missing_autosomes(self):
        """Test 9: Missing autosomes (only has chr 1-10) - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 11, n_rows),  # Only 1-10
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "missing_autosomes_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on missing autosomes"
            assert "missing" in message.lower() or (errors and any("missing" in str(e).lower() for e in errors))
    
    def test_invalid_allele(self):
        """Test 10: Invalid allele (non-ACTG) - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': ['X'] + list(np.random.choice(['A', 'T', 'G', 'C'], n_rows - 1)),  # Invalid
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid_allele_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on invalid allele"
            assert "allele" in message.lower() or (errors and any("allele" in str(e).lower() for e in errors))
    
    def test_negative_odds_ratio(self):
        """Test 11: Negative odds_ratio - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'odds_ratio': [-0.5] + list(np.random.uniform(0.1, 10, n_rows - 1)),  # Invalid
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "negative_or_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on negative odds_ratio"
    
    def test_invalid_file_extension(self):
        """Test 12: Invalid file extension - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid_ext_test.csv"  # Wrong extension
            data.to_csv(test_file, sep='\t', index=False)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                test_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on invalid file extension"
            assert "extension" in message.lower() or (errors and any("extension" in str(e).lower() for e in errors))
    
    def test_valid_compressed_file(self):
        """Test 13: Valid compressed .tsv.gz file - should pass."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "valid_compressed_test"
            ssf_file = create_test_file(data, test_file, compress=True)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == True, f"Should pass compressed file: {message}"
    
    def test_invalid_effect_allele_frequency(self):
        """Test 14: Invalid effect_allele_frequency (outside [0,1]) - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': [1.5, -0.1] + list(np.random.uniform(0, 1, n_rows - 2)),  # Invalid
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid_eaf_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on invalid effect_allele_frequency"
            assert "effect_allele_frequency" in message.lower() or (errors and any("effect_allele_frequency" in str(e).lower() for e in errors))
    
    def test_invalid_rsid_format(self):
        """Test 15: Invalid rsid format - should fail."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': np.random.randint(1, 23, n_rows),
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
            'rsid': ['invalid_rsid', 'rs123'] + ['rs' + str(i) for i in range(100000, 100000 + n_rows - 2)],
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid_rsid_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            assert is_valid == False, "Should fail on invalid rsid format"
            assert "rsid" in message.lower() or (errors and any("rsid" in str(e).lower() for e in errors))
    
    def test_chromosome_x_only(self):
        """Test 16: File with only chromosome X (23) - should pass (special case)."""
        n_rows = 150000
        data = pd.DataFrame({
            'chromosome': [23] * n_rows,  # Only X
            'base_pair_location': np.random.randint(1, 250000000, n_rows),
            'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n_rows),
            'beta': np.random.normal(0, 0.1, n_rows),
            'standard_error': np.abs(np.random.normal(0, 0.05, n_rows)),
            'effect_allele_frequency': np.random.uniform(0, 1, n_rows),
            'p_value': np.random.uniform(1e-10, 1, n_rows),
        })
        
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "chr_x_only_test"
            ssf_file = create_test_file(data, test_file)
            
            log = Log()
            is_valid, message, errors = validate_ssf_file(
                ssf_file, log=log, verbose=False, minimum_rows=100000
            )
            
            # Special case: only chromosome X is allowed
            assert is_valid == True, f"Should pass with only chromosome X: {message}"




if __name__ == "__main__":
    pytest.main([__file__, "-v"])

