"""Validate NumPy docstring style across gwaslab."""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest
from docstring_parser import parse as parse_docstring

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "gwaslab"

_FORBIDDEN_HEADERS = re.compile(r"^\s*(Args|Parameters|Returns|Raises|Yields):\s*$", re.MULTILINE)
_DEFAULT_EQUALS = re.compile(r"^\s*\w+\s*:\s*.+, default=", re.MULTILINE)


def _python_files() -> list[Path]:
    return sorted(SRC.rglob("*.py"))


def _docstrings_in_file(path: Path) -> list[str]:
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError:
        return []
    docs: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Module)):
            doc = ast.get_docstring(node, clean=False)
            if doc:
                docs.append(doc)
    return docs


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: p.relative_to(ROOT).as_posix())
def test_no_google_or_colon_section_headers(path: Path):
    for doc in _docstrings_in_file(path):
        assert not _FORBIDDEN_HEADERS.search(doc), f"{path}: docstring uses non-NumPy section header"


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: p.relative_to(ROOT).as_posix())
def test_no_default_equals_in_param_lines(path: Path):
    for doc in _docstrings_in_file(path):
        if "Parameters" not in doc:
            continue
        assert not _DEFAULT_EQUALS.search(doc), f"{path}: parameter line uses 'default=' instead of 'default '"


def test_g_numpy_doc_helpers():
    from gwaslab.info.g_numpy_doc import format_parameter, infer_type_annotation, normalize_sections

    assert infer_type_annotation(False) == "bool, default False"
    assert infer_type_annotation(None, has_default=True) == "Any, optional"
    lines = format_parameter("verbose", "Print progress", True)
    assert lines[0] == "verbose : bool, default True"
    normalized = normalize_sections("Parameters:\n-----------\n")
    assert "Parameters:" not in normalized
    assert "Parameters\n----------" in normalized


def test_sumstats_fix_id_has_parameters():
    from gwaslab.g_Sumstats import Sumstats

    doc = Sumstats.fix_id.__doc__ or ""
    assert "Parameters" in doc
    parsed = parse_docstring(doc)
    names = {p.arg_name for p in parsed.params}
    assert "fixprefix" in names
    assert "fixchrpos" in names


def test_sumstats_adapted_methods_hide_log():
    from gwaslab.g_Sumstats import Sumstats

    for name in ("fix_id", "filter_value", "remove_dup"):
        params = {p.arg_name for p in parse_docstring(getattr(Sumstats, name).__doc__ or "").params}
        assert "log" not in params


def test_docstring_parser_parses_sumstats_basic_check_parameters():
    from gwaslab.g_Sumstats import Sumstats

    parsed = parse_docstring(Sumstats.basic_check.__doc__ or "")
    names = {p.arg_name for p in parsed.params}
    assert "remove" in names
    assert "threads" in names


def test_docstring_parser_parses_plot_region_parameters():
    from gwaslab.g_Sumstats import Sumstats

    parsed = parse_docstring(Sumstats.plot_region.__doc__ or "")
    names = {p.arg_name for p in parsed.params}
    assert "region" in names
    assert "sig_level" in names
    sig = next(p for p in parsed.params if p.arg_name == "sig_level")
    assert sig.type_name is not None


@pytest.fixture(scope="module")
def gwaslab_module():
    import griffe

    return griffe.load("gwaslab", search_paths=[str(ROOT / "src")])


def test_griffe_loads_sumstats_class(gwaslab_module):
    assert "Sumstats" in gwaslab_module.classes


def test_normalize_docstrings_idempotent():
    from gwaslab.info.g_numpy_doc import normalize_sections

    sample = '''"""
Summary.

Parameters
----------
x : bool, default False
    Flag.
"""
'''
    once = normalize_sections(sample.strip())
    twice = normalize_sections(once)
    assert once == twice
