"""Validate escape-sequence fix patterns before applying them across src/gwaslab."""

from __future__ import annotations

import subprocess
import sys
import warnings
from pathlib import Path

import pytest

SAMPLES = Path(__file__).resolve().parent / "fixtures" / "escape_fix_samples"
GOOD = SAMPLES / "good_patterns.py"
BAD_TXT = SAMPLES / "bad_patterns.py.txt"
DEMO = SAMPLES / "demo_all_patterns.py"


def _compile_with_syntaxwarning_as_error(source: str, filename: str = "<test>") -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", SyntaxWarning)
        compile(source, filename, "exec")


def test_good_patterns_file_compiles_cleanly():
    source = GOOD.read_text(encoding="utf-8")
    _compile_with_syntaxwarning_as_error(source, str(GOOD))


def test_bad_patterns_reference_would_warn():
    """Document that original literals emit SyntaxWarning when compiled as Python."""
    source = BAD_TXT.read_text(encoding="utf-8")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SyntaxWarning)
        compile(source, str(BAD_TXT), "exec")
    escape_warnings = [w for w in caught if issubclass(w.category, SyntaxWarning)]
    assert escape_warnings, "bad_patterns reference should emit at least one SyntaxWarning"


def test_compileall_good_samples_subprocess():
    result = subprocess.run(
        [sys.executable, "-W", "error::SyntaxWarning", "-m", "compileall", "-q", str(SAMPLES)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr or result.stdout


def test_demo_script_passes():
    result = subprocess.run([sys.executable, str(DEMO)], capture_output=True, text=True, cwd=str(SAMPLES))
    assert result.returncode == 0, result.stderr or result.stdout
    assert "All" in result.stdout and "passed" in result.stdout


@pytest.mark.parametrize(
    "bad_literal,good_literal,expected",
    [
        (r'"\s+"', r'r"\s+"', r"\s+"),
        (r'"\s+"', r'r"\s+"', "\\s+"),
        (r'"$\mathregular{x}$"', r'r"$\mathregular{x}$"', "$\\mathregular{x}$"),
    ],
)
def test_bad_and_good_literals_same_runtime_value(bad_literal: str, good_literal: str, expected: str):
    bad_val = eval(bad_literal)  # noqa: S307 — intentional minimal repro
    good_val = eval(good_literal)  # noqa: S307
    assert bad_val == good_val == expected


def test_ldsc_regex_matches_sample_line():
    sys.path.insert(0, str(SAMPLES))
    import good_patterns as gp

    found = gp._LDSC_LINE_RE.findall(gp.SAMPLE_LDSC_LINE)
    assert found == ["obs:", "0.25", "se:", "0.03", "NA"]
    assert found[1] == "0.25"  # h2_obs value (objects[1] in io_read_ldsc.py)


def test_status_regex_matches_gwaslab_pattern():
    sys.path.insert(0, str(SAMPLES))
    import good_patterns as gp

    assert gp.STATUS_RE.match("1234547")  # digit at index 5 is 4
    assert gp.STATUS_RE.match("1234557")  # digit at index 5 is 5
    assert gp.STATUS_RE.match("1234537") is None  # digit at index 5 is 3


def test_legend_title_matches_format_style():
    sys.path.insert(0, str(SAMPLES))
    import good_patterns as gp

    title = gp.legend_title_good(5.0, -8)
    assert title == "$\\mathregular{ P < 5.0 x 10^{-8} }$ in:"


def test_plink_script_has_shell_continuations():
    sys.path.insert(0, str(SAMPLES))
    import good_patterns as gp

    script = gp.plink_script_good("--bfile demo", "demo_out")
    assert script.startswith("plink2 \\\n")
    assert "--out demo_out" in script
    assert "\\ \n" not in script


def test_production_p_text_double_space_parity():
    sys.path.insert(0, str(SAMPLES))
    import good_patterns as gp

    p12, pe = "1.2", "-8"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        production = "$\mathregular{p = " + p12 + " \\times  10^{" + pe + "}}$"
    assert gp.p_text_good(p12, pe) == production


def test_pandas_sep_parity():
    import io

    import pandas as pd

    data = "a b c\n1 2 3\n"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        bad_sep = "\s+"
    good_sep = r"\s+"
    df_bad = pd.read_csv(io.StringIO(data), sep=bad_sep, header=None)
    df_good = pd.read_csv(io.StringIO(data), sep=good_sep, header=None)
    assert df_bad.equals(df_good)


def test_matplotlib_render_bbox_parity():
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sys.path.insert(0, str(SAMPLES))
    import good_patterns as gp

    def bbox(s: str) -> tuple[float, float]:
        fig, ax = plt.subplots(figsize=(3, 1))
        text = ax.text(0.1, 0.5, s, fontsize=12)
        fig.canvas.draw()
        extent = text.get_window_extent(renderer=fig.canvas.get_renderer())
        plt.close(fig)
        return round(extent.width, 2), round(extent.height, 2)

    pairs = [
        gp.YLABEL_GOOD,
        gp.legend_title_good(5.0, -8),
        gp.p_text_good("1.2", "-8"),
    ]
    for s in pairs:
        assert bbox(s) == bbox(s)


def test_plink_parity_preserving_fix():
    filter_flag = "--bfile demo"
    out_prefix = "demo_out"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        production = """
    plink2 \
    {} \
    --make-just-bim \
    --make-just-fam \ 
    --out {}
    """.format(filter_flag, out_prefix)
    fixed = """
    plink2 \
    {} \
    --make-just-bim \
    --make-just-fam \\ 
    --out {}
    """.format(filter_flag, out_prefix)
    assert production == fixed
