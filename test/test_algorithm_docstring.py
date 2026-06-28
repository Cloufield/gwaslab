"""Validate NumPy docstring style under gwaslab.algorithm."""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ALGORITHM = ROOT / "src" / "gwaslab" / "algorithm"

_PROSE_DEFAULT = re.compile(r"Default is ", re.IGNORECASE)
_RETURNS_OVER_INDENT = re.compile(
    r"^Returns\n-+\n    \S",
    re.MULTILINE,
)
_PARAM_SECTION = re.compile(r"^\s*Parameters\s*$", re.MULTILINE)


def _algorithm_python_files() -> list[Path]:
    return sorted(
        p for p in ALGORITHM.rglob("*.py") if p.name != "__init__.py"
    )


def _public_functions(path: Path) -> list[tuple[str, ast.FunctionDef | ast.AsyncFunctionDef]]:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    out: list[tuple[str, ast.FunctionDef | ast.AsyncFunctionDef]] = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
            out.append((node.name, node))
    return out


@pytest.mark.parametrize("path", _algorithm_python_files(), ids=lambda p: p.relative_to(ROOT).as_posix())
def test_algorithm_no_prose_default(path: Path):
    text = path.read_text(encoding="utf-8")
    assert not _PROSE_DEFAULT.search(text), f"{path}: docstring uses 'Default is' prose"


@pytest.mark.parametrize("path", _algorithm_python_files(), ids=lambda p: p.relative_to(ROOT).as_posix())
def test_algorithm_returns_not_over_indented(path: Path):
    for _, node in _public_functions(path):
        doc = ast.get_docstring(node, clean=False)
        if not doc or "Returns" not in doc:
            continue
        assert not _RETURNS_OVER_INDENT.search(doc), (
            f"{path}:{node.name}: Returns type line is over-indented"
        )


@pytest.mark.parametrize("path", _algorithm_python_files(), ids=lambda p: p.relative_to(ROOT).as_posix())
def test_algorithm_public_functions_with_args_have_parameters(path: Path):
    for name, node in _public_functions(path):
        args = [
            a.arg
            for a in node.args.args
            if a.arg not in {"self", "cls"}
        ]
        if not args:
            continue
        doc = ast.get_docstring(node) or ""
        assert _PARAM_SECTION.search(doc), (
            f"{path}:{name}: public function with args missing Parameters section"
        )


def test_fix_returns_indent_idempotent():
    from gwaslab.info.g_numpy_doc import normalize_sections

    sample = '''"""Summary.

Parameters
----------
x : float
    Input.

Returns
-------
    float
        Output.
"""
'''
    once = normalize_sections(sample.strip())
    twice = normalize_sections(once)
    assert once == twice
    assert "Returns\n-------\nfloat\n" in once
    assert "    float\n" not in once.split("Returns\n-------\n")[1].split("\n\n")[0]
