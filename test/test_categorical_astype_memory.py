"""Regression tests for #222: Categorical.astype(str) memory blow-up on allele columns."""

import os
import re
import sys
import tracemalloc
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm.hm_harmonize_sumstats import is_palindromic, check_status
from gwaslab.hm.hm_assign_rsid import _expand_multiallelic_fast
from gwaslab.util.util_ex_ldsc import _munge_sumstats
from gwaslab.util.util_in_filter_value import _filter_values, _filter_palindromic, _filter_indel, _filter_snp
from gwaslab.qc.qc_check_datatype import (
    categorical_str_len,
    categorical_str_upper,
    categorical_str_contains,
)
from gwaslab.extension.gwas_sumstats_tools.validate_ssf import _validate_dataframe
from gwaslab.io.io_to_formats import _check_indel
from gwaslab.qc.qc_fix_sumstats import _fix_allele, _flip_allele_stats
from gwaslab.info.g_Log import Log

N_ROWS = 100_000
LONG_LEN = 662
MAX_PEAK_MB = 100.0
GWASLAB_SRC = Path(ROOT) / "src" / "gwaslab"
UNSAFE_ALLELE_ASTYPE_STR = re.compile(
    r'\[(["\'])(EA|NEA|REF|ALT)\1\]\.astype\(str\)'
)
TO_NUMPY_UNICODE_ASTYPE = re.compile(
    r'\.to_numpy\(\)\.astype\(f?[\'"]<U'
)
UNGUARDED_ALLELE_STR = re.compile(
    r'\[(["\'])(EA|NEA|REF|ALT)\1\]\.str\.'
)
FULL_COLUMN_CATEGORICAL_SAFE_STR = re.compile(
    r'categorical_safe_str\(sumstats\[(["\'])(EA|NEA|REF|ALT)\1\]\)'
)
SUMSTATS_SCAN_PARTS = {
    "qc", "hm", "io", "util", "info", "extension", "bd", "viz", "algorithm",
}
HM_QC_IO_SCAN_PARTS = {"qc", "hm", "io"}
TO_NUMPY_UNICODE_ALLOWLIST = {
    "hm_harmonize_sumstats.py",
}
FULL_COLUMN_CATEGORICAL_SAFE_STR_ALLOWLIST = {
    "hm_casting.py",  # _ASET construction requires full-column allele strings
}
UNGUARDED_ALLELE_STR_ALLOWLIST: set[str] = set()


def _allowed_astype_str_line(line: str, filepath: Path) -> bool:
    stripped = line.strip()
    if stripped.startswith("#"):
        return True
    if ".index.astype(str)" in line or ".columns.astype(str)" in line:
        return True
    if "dtype=str" in line or "dtype = str" in line:
        return True
    if "astype(str)" in line and (
        '"""' in line
        or "'''" in line
        or "``" in line
        or "Series.astype(str)" in line
        or "Categorical.astype(str)" in line
    ):
        return True
    return False


def _line_has_safe_str_context(lines: list[str], line_idx: int, window: int = 8) -> bool:
    start = max(0, line_idx - window)
    end = min(len(lines), line_idx + window + 1)
    context = "\n".join(lines[start:end])
    return (
        "categorical_safe_str" in context
        or "categorical_str_len" in context
        or "categorical_str_upper" in context
        or "categorical_str_contains" in context
    )


def _peak_mb(fn):
    tracemalloc.start()
    tracemalloc.reset_peak()
    fn()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak / 1e6


def _make_categorical_allele_df(include_long_allele: bool = True) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    alleles = ["A", "C", "G", "T"]
    long_allele = "A" * LONG_LEN
    ea_codes = rng.integers(0, 4, size=N_ROWS)
    nea_codes = rng.integers(0, 4, size=N_ROWS)
    if include_long_allele:
        ea_codes[0] = 4
    categories = alleles + ([long_allele] if include_long_allele else [])
    cat_index = pd.Index(categories)
    return pd.DataFrame(
        {
            "CHR": np.ones(N_ROWS, dtype=int),
            "POS": np.arange(N_ROWS, dtype=np.int64) + 1,
            "EA": pd.Categorical.from_codes(ea_codes, categories=cat_index),
            "NEA": pd.Categorical.from_codes(nea_codes, categories=cat_index),
            "P": rng.uniform(0.001, 0.5, size=N_ROWS),
            "BETA": rng.normal(size=N_ROWS),
        }
    )


def _make_categorical_allele_df_with_status(include_long_allele: bool = True) -> pd.DataFrame:
    df = _make_categorical_allele_df(include_long_allele=include_long_allele)
    df["STATUS"] = pd.array([9909999] * N_ROWS, dtype="Int64")
    df.loc[0, "STATUS"] = 9909440  # digit 6 = 4 -> reverse-complement branch in _flip_allele_stats
    return df


class TestCategoricalAstypeMemory(unittest.TestCase):
    def test_hm_casting_allele_set_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)
        mold = sumstats.iloc[:10].copy()

        peak = _peak_mb(
            lambda: _merge_mold_with_sumstats_by_chrpos(mold, sumstats, verbose=False)
        )
        self.assertLess(peak, MAX_PEAK_MB, f"hm_casting merge peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_ldsc_allele_filter_long_allele(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)
        log = Log()

        peak = _peak_mb(
            lambda: _munge_sumstats(
                sumstats,
                log,
                info=0.0,
                maf=0.0,
                nopalindromic=False,
                exclude_hla=False,
                exclude_sexchr=False,
                verbose=False,
            )
        )
        self.assertLess(peak, MAX_PEAK_MB, f"LDSC munge peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_filter_value_ea_str_len(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)

        peak = _peak_mb(
            lambda: _filter_values(sumstats, "EA.str.len() > 0", verbose=False)
        )
        self.assertLess(peak, MAX_PEAK_MB, f"filter_value peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_validate_ssf_categorical_alleles(self):
        df = _make_categorical_allele_df(include_long_allele=True)
        df = df.rename(columns={"EA": "effect_allele", "NEA": "other_allele", "CHR": "chromosome", "POS": "base_pair_location"})
        df["standard_error"] = 0.01
        df["effect_allele_frequency"] = 0.3
        df["p_value"] = 0.05
        df["beta"] = 0.1

        peak = _peak_mb(
            lambda: _validate_dataframe(
                df,
                header=list(df.columns),
                effect_field="beta",
                p_value_field="p_value",
                pval_zero=False,
                log=Log(),
                verbose=False,
            )
        )
        self.assertLess(peak, MAX_PEAK_MB, f"validate_ssf peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_is_palindromic_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)

        peak = _peak_mb(lambda: is_palindromic(sumstats))
        self.assertLess(peak, MAX_PEAK_MB, f"is_palindromic peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_filter_palindromic_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)

        peak = _peak_mb(
            lambda: _filter_palindromic(sumstats, mode="out", verbose=False, log=Log())
        )
        self.assertLess(peak, MAX_PEAK_MB, f"filter_palindromic peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_check_status_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df_with_status(include_long_allele=True)[
            ["CHR", "POS", "EA", "NEA", "STATUS"]
        ]
        record = np.full(500, 2, dtype=np.uint8)
        starting_positions_dict = {1: 0}
        records_len_dict = {1: len(record)}

        peak = _peak_mb(
            lambda: check_status(
                sumstats,
                record=record,
                starting_positions_dict=starting_positions_dict,
                records_len_dict=records_len_dict,
                verbose=False,
            )
        )
        self.assertLess(peak, MAX_PEAK_MB, f"check_status peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_fix_allele_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df_with_status(include_long_allele=True)

        peak = _peak_mb(
            lambda: _fix_allele(sumstats, remove=False, verbose=False, log=Log())
        )
        self.assertLess(peak, MAX_PEAK_MB, f"_fix_allele peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_flip_allele_stats_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df_with_status(include_long_allele=True)

        peak = _peak_mb(
            lambda: _flip_allele_stats(sumstats, verbose=False, log=Log())
        )
        self.assertLess(peak, MAX_PEAK_MB, f"_flip_allele_stats peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_check_indel_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)
        log = Log()

        peak = _peak_mb(
            lambda: _check_indel(sumstats, log, verbose=False)
        )
        self.assertLess(peak, MAX_PEAK_MB, f"_check_indel peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_expand_multiallelic_fast_no_memory_blowup(self):
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["A,T", "G"] + ["C"] * (N_ROWS - 2)),
                "NEA": pd.Categorical(["G"] * N_ROWS),
            }
        )

        peak = _peak_mb(
            lambda: _expand_multiallelic_fast(df, "EA", "NEA")
        )
        self.assertLess(peak, MAX_PEAK_MB, f"_expand_multiallelic_fast peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_filter_indel_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)
        peak = _peak_mb(
            lambda: _filter_indel(sumstats, mode="in", verbose=False, log=Log())
        )
        self.assertLess(peak, MAX_PEAK_MB, f"_filter_indel peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")

    def test_filter_snp_no_memory_blowup(self):
        sumstats = _make_categorical_allele_df(include_long_allele=True)
        peak = _peak_mb(
            lambda: _filter_snp(sumstats, mode="in", verbose=False, log=Log())
        )
        self.assertLess(peak, MAX_PEAK_MB, f"_filter_snp peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")


class TestCategoricalStrHelpers(unittest.TestCase):
    def test_categorical_str_len_long_category_no_blowup(self):
        long_allele = "A" * LONG_LEN
        ser = pd.Series(
            pd.Categorical.from_codes(
                [0, 1, 2, 3, 4],
                categories=pd.Index(["A", "C", "G", "T", long_allele]),
            )
        )
        peak = _peak_mb(lambda: categorical_str_len(ser))
        self.assertLess(peak, MAX_PEAK_MB, f"categorical_str_len peak {peak:.1f} MB exceeds {MAX_PEAK_MB} MB")
        lens = categorical_str_len(ser)
        self.assertEqual(int(lens.iloc[4]), LONG_LEN)
        self.assertEqual(int(lens.iloc[0]), 1)

    def test_categorical_str_upper_preserves_categorical(self):
        ser = pd.Categorical(["a", "c", "g", "t"])
        upper = categorical_str_upper(ser)
        self.assertIsInstance(upper.dtype, pd.CategoricalDtype)
        self.assertTrue(set(upper.cat.categories) <= {"A", "C", "G", "T"})

    def test_categorical_str_upper_merges_case_collisions(self):
        ser = pd.Categorical(["a", "A", "c"])
        upper = categorical_str_upper(ser)
        self.assertIsInstance(upper.dtype, pd.CategoricalDtype)
        self.assertTrue(set(upper.cat.categories) <= {"A", "C"})
        self.assertEqual(list(upper), ["A", "A", "C"])

    def test_categorical_str_upper_indel_collision(self):
        ser = pd.Categorical(["At", "at"])
        upper = categorical_str_upper(ser)
        self.assertIsInstance(upper.dtype, pd.CategoricalDtype)
        self.assertEqual(list(upper.cat.categories), ["AT"])
        self.assertEqual(list(upper), ["AT", "AT"])

    def test_fix_allele_mixed_case_categorical(self):
        df = pd.DataFrame(
            {
                "EA": pd.Categorical(["a", "A", "At", "G"]),
                "NEA": pd.Categorical(["g", "T", "A", "C"]),
                "STATUS": pd.array([9909999, 9909999, 9909999, 9909999], dtype="Int64"),
            }
        )
        _fix_allele(df, remove=False, verbose=False, log=Log())

    def test_categorical_str_contains_on_categories(self):
        ser = pd.Categorical(["A", "C", "XY", "T"])
        bad = categorical_str_contains(ser, "[^ATCG]", na=True, regex=True)
        self.assertTrue(bad.iloc[2])
        self.assertFalse(bad.iloc[0])


class TestNoUnsafeAstypeStrOnAlleles(unittest.TestCase):
    """CI guard: ban full-column .astype(str) on reserved allele headers."""

    def test_no_unsafe_astype_str_on_allele_columns(self):
        violations = []
        for py in GWASLAB_SRC.rglob("*.py"):
            for line_no, line in enumerate(py.read_text(encoding="utf-8").splitlines(), 1):
                if UNSAFE_ALLELE_ASTYPE_STR.search(line):
                    violations.append(f"{py.relative_to(ROOT)}:{line_no}: {line.strip()}")
        self.assertEqual(
            violations,
            [],
            "Found unsafe .astype(str) on EA/NEA/REF/ALT:\n" + "\n".join(violations),
        )

    def test_prefer_pandas_string_in_sumstats_modules(self):
        violations = []
        for py in GWASLAB_SRC.rglob("*.py"):
            if py.name == "qc_check_datatype.py":
                continue
            if not SUMSTATS_SCAN_PARTS.intersection(py.parts):
                continue
            for line_no, line in enumerate(py.read_text(encoding="utf-8").splitlines(), 1):
                if ".astype(str)" not in line:
                    continue
                if _allowed_astype_str_line(line, py):
                    continue
                violations.append(f"{py.relative_to(ROOT)}:{line_no}: {line.strip()}")
        self.assertEqual(
            violations,
            [],
            "Use astype('string') or categorical_safe_str instead of .astype(str):\n"
            + "\n".join(violations),
        )


class TestDtypeMemoryPatterns(unittest.TestCase):
    """Static scan for dtype conversion patterns that can cause memory spikes."""

    def test_to_numpy_unicode_astype_is_allowlisted(self):
        """Tier-1 sites using fixed-width unicode must be explicitly allowlisted."""
        found = {}
        for py in GWASLAB_SRC.rglob("*.py"):
            if not HM_QC_IO_SCAN_PARTS.intersection(py.parts):
                continue
            for line_no, line in enumerate(py.read_text(encoding="utf-8").splitlines(), 1):
                if TO_NUMPY_UNICODE_ASTYPE.search(line):
                    found.setdefault(py.name, []).append(f"{py.relative_to(ROOT)}:{line_no}: {line.strip()}")
        unexpected = {k: v for k, v in found.items() if k not in TO_NUMPY_UNICODE_ALLOWLIST}
        self.assertEqual(
            unexpected,
            {},
            "New .to_numpy().astype('<U...') sites must be reviewed and allowlisted:\n"
            + "\n".join(f"{k}:\n  " + "\n  ".join(v) for k, v in unexpected.items()),
        )
        self.assertIn(
            "hm_harmonize_sumstats.py",
            found,
            "Expected intentional fixed-width unicode path in hm_harmonize_sumstats.py",
        )

    def test_allele_str_ops_use_safe_conversion(self):
        violations = []
        for py in GWASLAB_SRC.rglob("*.py"):
            if py.name in UNGUARDED_ALLELE_STR_ALLOWLIST:
                continue
            if not HM_QC_IO_SCAN_PARTS.intersection(py.parts):
                continue
            lines = py.read_text(encoding="utf-8").splitlines()
            for line_no, line in enumerate(lines, 1):
                if not UNGUARDED_ALLELE_STR.search(line):
                    continue
                if _line_has_safe_str_context(lines, line_no - 1):
                    continue
                if line.strip().startswith("#"):
                    continue
                violations.append(f"{py.relative_to(ROOT)}:{line_no}: {line.strip()}")
        self.assertEqual(
            violations,
            [],
            "Wrap allele columns with categorical_safe_str before .str ops:\n"
            + "\n".join(violations),
        )

    def test_no_full_column_categorical_safe_str_on_alleles(self):
        violations = []
        for py in GWASLAB_SRC.rglob("*.py"):
            if py.name in FULL_COLUMN_CATEGORICAL_SAFE_STR_ALLOWLIST:
                continue
            if not {"qc", "hm"}.intersection(py.parts):
                continue
            for line_no, line in enumerate(py.read_text(encoding="utf-8").splitlines(), 1):
                if FULL_COLUMN_CATEGORICAL_SAFE_STR.search(line):
                    violations.append(f"{py.relative_to(ROOT)}:{line_no}: {line.strip()}")
        self.assertEqual(
            violations,
            [],
            "Avoid full-column categorical_safe_str(sumstats['EA']) materialization:\n"
            + "\n".join(violations),
        )


if __name__ == "__main__":
    unittest.main()
