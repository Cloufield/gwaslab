"""
Test suite for SSF validator.
"""
import os
import sys
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from gwaslab.extension.gwas_sumstats_tools.validate_ssf import validate_ssf_file
from gwaslab.info.g_Log import Log


def create_test_ssf_file(
    data: pd.DataFrame,
    filename: Path,
    compress: bool = False
) -> Path:
    """Create a test SSF file."""
    if compress:
        output_file = filename.with_suffix('.tsv.gz')
        data.to_csv(output_file, sep='\t', index=False, compression='gzip')
    else:
        output_file = filename.with_suffix('.tsv')
        data.to_csv(output_file, sep='\t', index=False)
    return output_file



