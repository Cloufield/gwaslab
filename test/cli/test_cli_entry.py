"""Smoke tests for the `gwaslab` CLI entry point (Python API)."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd
import pytest

from gwaslab_cli.main import main

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_main_version_exits_cleanly() -> None:
    """`gwaslab version` fast path should not raise."""
    main(["version"])


def test_cli_help_does_not_crash() -> None:
    """`--help` should print and exit via SystemExit(0)."""
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0


def test_cli_qc_writes_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Minimal QC run with output file (matches public CLI examples)."""
    raw = REPO_ROOT / "test" / "raw" / "dirty_sumstats.tsv"
    if not raw.is_file():
        pytest.skip(f"missing test data: {raw}")
    out = tmp_path / "cleaned.tsv"
    monkeypatch.chdir(tmp_path)
    main(
        [
            "--input",
            str(raw),
            "--qc",
            "--output",
            str(out),
            "--quiet",
        ]
    )
    assert out.is_file() or any(tmp_path.glob("cleaned*"))


def test_cli_formatbook_list() -> None:
    """Subcommand `formatbook list` should run."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        main(["formatbook", "list"])
    text = buf.getvalue()
    assert text.strip()


def test_cli_plot_qq_routes_to_plot_mqq(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """`--plot qq` must call `plot_mqq(mode='qq')` with the requested save path."""
    import gwaslab_cli.cli as cli

    called: dict[str, object] = {}

    class _DummyLog:
        def write(self, *args, **kwargs) -> None:
            return None

    class _DummySumstats:
        def __init__(self) -> None:
            self.data = pd.DataFrame({"CHR": [1], "POS": [12345], "P": [0.5]})
            self.log = _DummyLog()

        def fix_chr(self, **kwargs) -> None:
            return None

        def fix_pos(self, **kwargs) -> None:
            return None

        def plot_mqq(self, **kwargs) -> None:
            called.update(kwargs)

        def plot_qq(self, **kwargs) -> None:
            raise AssertionError("CLI qq plotting should use plot_mqq(mode='qq'), not plot_qq()")

    dummy = _DummySumstats()
    monkeypatch.setattr(cli, "emit_cli_mode_banner", lambda *args, **kwargs: None)
    monkeypatch.setattr(cli, "load_sumstats", lambda *args, **kwargs: dummy)

    out = tmp_path / "qq.png"
    cli.main(["--input", "dummy.tsv", "--fmt", "auto", "--plot", "qq", "--out", str(out), "--quiet"])

    assert called.get("mode") == "qq"
    assert called.get("save") == str(out)


def test_gwaslab_dunder_version_matches_version_info() -> None:
    """Top-level package should expose `__version__` for standard tooling."""
    import gwaslab
    from gwaslab.info.g_version import gwaslab_info

    assert hasattr(gwaslab, "__version__")
    assert gwaslab.__version__ == gwaslab_info()["version"]
