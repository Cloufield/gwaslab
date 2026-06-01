"""Tests for builtin HapMap3 core EAF panels."""

from pathlib import Path

import pandas as pd
import pytest

from gwaslab.bd.bd_ancestry_ref import builtin_eaf_path, resolve_ancestry_af
from gwaslab.util.util_ex_infer_ancestry import _infer_ancestry, calculate_fst

DATA = Path(__file__).resolve().parents[1] / "src" / "gwaslab" / "data" / "hapmap3_EAF"
HG19_CORE = DATA / "PAN.hapmap3.hg19.EAF.core.tsv.gz"
HG38_CORE = DATA / "PAN.hapmap3.hg38.EAF.core.tsv.gz"

META = ["SNPID", "CHR", "POS", "REF", "ALT"]
POPS = [
    "GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV",
    "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB",
    "JPT", "LWK", "ASW", "MXL", "TSI", "GIH",
    "EUR", "EAS", "AMR", "SAS", "AFR",
]
EXPECTED_COLS = META + POPS
MAX_TOTAL_BYTES = 5_000_000


@pytest.fixture(scope="module")
def hg19_core_df():
    return pd.read_csv(HG19_CORE, sep="\t", nrows=500)


def test_core_files_exist():
    assert HG19_CORE.is_file(), f"missing {HG19_CORE}"
    assert HG38_CORE.is_file(), f"missing {HG38_CORE}"


def test_combined_size_under_5mb():
    total = HG19_CORE.stat().st_size + HG38_CORE.stat().st_size
    assert total < MAX_TOTAL_BYTES, f"combined gzip size {total} exceeds 5 MB"


def test_core_header_columns(hg19_core_df):
    assert list(hg19_core_df.columns) == EXPECTED_COLS


def test_core_variant_count():
    import gzip
    with gzip.open(HG19_CORE, "rt") as f:
        n19 = sum(1 for _ in f) - 1
    with gzip.open(HG38_CORE, "rt") as f:
        n38 = sum(1 for _ in f) - 1
    assert 10_000 <= n19 <= 35_000
    assert n19 == n38


def test_builtin_eaf_path():
    assert builtin_eaf_path("19").name == "PAN.hapmap3.hg19.EAF.core.tsv.gz"
    assert builtin_eaf_path("38").name == "PAN.hapmap3.hg38.EAF.core.tsv.gz"


def test_resolve_builtin_when_not_downloaded(monkeypatch):
    monkeypatch.setattr(
        "gwaslab.bd.bd_ancestry_ref.get_path",
        lambda *a, **k: False,
    )
    p = resolve_ancestry_af(None, "19", verbose=False)
    assert p.endswith("PAN.hapmap3.hg19.EAF.core.tsv.gz")


def test_infer_ancestry_runs_on_core_panel():
    ref = pd.read_csv(HG19_CORE, sep="\t", nrows=150)
    sumstats = pd.DataFrame({
        "CHR": ref["CHR"].values,
        "POS": ref["POS"].values,
        "EA": ref["ALT"].values,
        "NEA": ref["REF"].values,
        "EAF": ref["EUR"].values,
    })
    label = _infer_ancestry(sumstats, ancestry_af=str(HG19_CORE), build="19", verbose=False)
    assert label in POPS
    assert len(label) >= 2


def test_calculate_fst_bounds():
    assert 0 <= calculate_fst(0.8, 0.2) <= 1
