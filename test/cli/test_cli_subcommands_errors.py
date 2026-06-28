"""Subcommands, parser errors, launcher edge cases, and extra pipeline coverage."""

from __future__ import annotations

import io
import json
import re
import sys
from contextlib import redirect_stdout
from importlib.metadata import PackageNotFoundError
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW = REPO_ROOT / "test" / "raw" / "dirty_sumstats.tsv"


@pytest.fixture(autouse=True)
def _ensure_src_path() -> None:
    src = str(REPO_ROOT / "src")
    if src not in sys.path:
        sys.path.insert(0, src)


def test_cli_requires_input_or_subcommand() -> None:
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main([])
    assert exc.value.code == 2


def test_subcommand_version_via_full_parser() -> None:
    """`gwaslab_cli.cli.main` path uses handler `run_version` (not launcher fast path)."""
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["version"])
    assert buf.getvalue().strip()


def test_launcher_main_uses_sys_argv_when_argv_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from gwaslab_cli.main import main as launcher_main

    monkeypatch.setattr(sys, "argv", ["gwaslab", "version"])
    buf = io.StringIO()
    with redirect_stdout(buf):
        launcher_main()
    assert buf.getvalue().strip()


def test_launcher_version_when_metadata_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    from gwaslab_cli.main import main as launcher_main
    import gwaslab_cli.main as launcher_mod

    def _boom(_name: str) -> str:
        raise PackageNotFoundError

    monkeypatch.setattr(launcher_mod, "version", _boom)
    buf = io.StringIO()
    with redirect_stdout(buf):
        launcher_main(["version"])
    assert "unavailable" in buf.getvalue()


def test_config_subcommand_text() -> None:
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["config"])
    text = buf.getvalue()
    assert "paths" in text.lower() or ":" in text


def test_config_subcommand_json() -> None:
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["config", "--json"])
    data = json.loads(buf.getvalue())
    assert "paths" in data


def test_config_set_data_directory(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    import gwaslab as gl

    target = tmp_path / "refs"
    target.mkdir()
    buf = io.StringIO()
    with redirect_stdout(buf):
        from gwaslab_cli.cli import main as cli_main

        cli_main(["config", "set", "data_directory", str(target)])
    assert gl.options.paths["data_directory"].startswith(str(target))
    assert "data_directory:" in buf.getvalue()


def test_init_subcommand_invokes_scan(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock_scan = MagicMock(return_value=True)
    mock_set_dir = MagicMock()
    mock_ensure = MagicMock()
    monkeypatch.setattr(gwaslab, "scan_downloaded_files", mock_scan)
    monkeypatch.setattr(gwaslab, "set_default_directory", mock_set_dir)
    monkeypatch.setattr("gwaslab.bd.bd_config.ensure_user_layout", mock_ensure)
    from gwaslab_cli.cli import main as cli_main

    ref_dir = tmp_path / "refs"
    ref_dir.mkdir()
    cli_main(["init", "--dir", str(ref_dir)])
    mock_ensure.assert_called_once()
    mock_set_dir.assert_called_once()
    mock_scan.assert_called_once()
    _, kwargs = mock_scan.call_args
    assert kwargs["verbose"] is True
    assert Path(kwargs["directory"]) == ref_dir.resolve()


def test_init_subcommand_quiet(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock_scan = MagicMock(return_value=True)
    mock_ensure = MagicMock()
    monkeypatch.setattr(gwaslab, "scan_downloaded_files", mock_scan)
    monkeypatch.setattr("gwaslab.bd.bd_config.ensure_user_layout", mock_ensure)
    from gwaslab_cli.cli import main as cli_main

    cli_main(["init", "--quiet"])
    _, kwargs = mock_scan.call_args
    assert kwargs["verbose"] is False
    assert kwargs["directory"] == gwaslab.get_default_directory()


def test_config_show_config_json() -> None:
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["config", "show", "config"])
    payload = json.loads(buf.getvalue())
    assert isinstance(payload, dict)


def test_path_config_resolves() -> None:
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["path", "config"])
    line = buf.getvalue().strip().splitlines()[-1]
    assert Path(line).suffix == ".json" or "config" in line.lower()


def test_formatbook_update_invokes_library(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock = MagicMock()
    monkeypatch.setattr(gwaslab, "update_formatbook", mock)
    from gwaslab_cli.cli import main as cli_main

    cli_main(["formatbook", "update"])
    mock.assert_called_once()


def test_download_sumstats_invokes_library(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock = MagicMock()
    monkeypatch.setattr(gwaslab, "download_sumstats", mock)
    from gwaslab_cli.cli import main as cli_main

    cli_main(["download-sumstats", "GCST000000", "--output-dir", "/tmp/gwaslab_cli_test_dl"])
    mock.assert_called_once()
    call_kw = mock.call_args.kwargs
    assert call_kw.get("gcst_id") == "GCST000000"
    assert call_kw.get("output_dir") == "/tmp/gwaslab_cli_test_dl"


def test_download_ref_group_invokes_library(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock = MagicMock()
    monkeypatch.setattr(gwaslab, "download_ref", mock)
    from gwaslab_cli.cli import main as cli_main

    cli_main(["download", "ref", "1kg_eas_hg19", "--directory", "/tmp/gwaslab_cli_ref", "--overwrite"])
    mock.assert_called_once()
    call_kw = mock.call_args.kwargs
    assert call_kw.get("name") == "1kg_eas_hg19"
    assert call_kw.get("directory") == "/tmp/gwaslab_cli_ref"
    assert call_kw.get("overwrite") is True


def test_download_sumstats_group_invokes_library(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock = MagicMock()
    monkeypatch.setattr(gwaslab, "download_sumstats", mock)
    from gwaslab_cli.cli import main as cli_main

    cli_main(["download", "sumstats", "GCST000000", "-d", "/tmp/gwaslab_cli_test_dl"])
    mock.assert_called_once()
    call_kw = mock.call_args.kwargs
    assert call_kw.get("gcst_id") == "GCST000000"
    assert call_kw.get("output_dir") == "/tmp/gwaslab_cli_test_dl"


def test_download_sumstats_directory_alias(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock = MagicMock()
    monkeypatch.setattr(gwaslab, "download_sumstats", mock)
    from gwaslab_cli.cli import main as cli_main

    cli_main(["download-sumstats", "GCST000000", "--directory", "/tmp/gwaslab_cli_test_dl2"])
    mock.assert_called_once()
    assert mock.call_args.kwargs.get("output_dir") == "/tmp/gwaslab_cli_test_dl2"


def test_list_ref_json(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    monkeypatch.setattr(
        gwaslab,
        "check_available_ref",
        MagicMock(return_value={"ref_a": {"description": "test"}}),
    )
    monkeypatch.setattr(gwaslab, "check_downloaded_ref", MagicMock(return_value={}))
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["list", "ref", "--available", "--json", "-q"])
    data = json.loads(buf.getvalue())
    assert "available" in data
    assert "ref_a" in data["available"]
    assert "downloaded" not in data


def test_download_ref_invokes_library(monkeypatch: pytest.MonkeyPatch) -> None:
    from unittest.mock import MagicMock

    import gwaslab

    mock = MagicMock()
    monkeypatch.setattr(gwaslab, "download_ref", mock)
    from gwaslab_cli.cli import main as cli_main

    cli_main(
        [
            "download-ref",
            "1kg_eas_hg19",
            "--directory",
            "/tmp/gwaslab_cli_ref",
            "--local-filename",
            "1kg_eas_hg19.vcf.gz",
            "--overwrite",
        ]
    )
    mock.assert_called_once()
    call_kw = mock.call_args.kwargs
    assert call_kw.get("name") == "1kg_eas_hg19"
    assert call_kw.get("directory") == "/tmp/gwaslab_cli_ref"
    assert call_kw.get("local_filename") == "1kg_eas_hg19.vcf.gz"
    assert call_kw.get("overwrite") is True


def test_formatbook_list_json() -> None:
    from gwaslab_cli.cli import main as cli_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["formatbook", "list", "--json"])
    data = json.loads(buf.getvalue())
    assert "formats" in data
    assert data["formats"]


def test_formatbook_show_gwaslab(capsys: pytest.CaptureFixture[str]) -> None:
    from gwaslab_cli.cli import main as cli_main

    cli_main(["formatbook", "show", "gwaslab"])
    captured = capsys.readouterr().out
    dec = json.JSONDecoder()
    data = None
    for i, ch in enumerate(captured):
        if ch != "{":
            continue
        try:
            data, _ = dec.raw_decode(captured[i:])
            break
        except json.JSONDecodeError:
            continue
    assert data is not None
    assert "gwaslab" in data


def test_config_show_unknown_keyword_exits(capsys: pytest.CaptureFixture[str]) -> None:
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(["config", "show", "__not_a_real_gwaslab_path_keyword__xyz__"])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "not found" in captured.err.lower()


def test_path_unknown_keyword_exits(capsys: pytest.CaptureFixture[str]) -> None:
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(["path", "__not_a_real_gwaslab_path_keyword__xyz__"])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "not found" in (captured.err + captured.out).lower()


def test_plot_miami_errors() -> None:
    if not RAW.is_file():
        pytest.skip(f"missing {RAW}")
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(
            [
                "--input",
                str(RAW),
                "--plot",
                "miami",
                "--output",
                str(RAW.parent / "_miami_skip.png"),
                "--quiet",
            ]
        )
    assert exc.value.code == 2


def test_plot_regional_missing_region_errors() -> None:
    if not RAW.is_file():
        pytest.skip(f"missing {RAW}")
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(
            [
                "--input",
                str(RAW),
                "--plot",
                "regional",
                "--output",
                str(RAW.parent / "_regional_skip.png"),
                "--quiet",
            ]
        )
    assert exc.value.code == 2


def test_extract_proxy_errors() -> None:
    if not RAW.is_file():
        pytest.skip(f"missing {RAW}")
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(
            [
                "--input",
                str(RAW),
                "--get",
                "proxy",
                "--output",
                str(RAW.parent / "_proxy_skip.tsv"),
                "--quiet",
            ]
        )
    assert exc.value.code == 2


def test_extract_lead_writes_tsv(tmp_path: Path) -> None:
    if not RAW.is_file():
        pytest.skip(f"missing {RAW}")
    from gwaslab_cli.cli import main as cli_main

    out = tmp_path / "leads.tsv"
    cli_main(
        [
            "--input",
            str(RAW),
            "--get",
            "lead",
            "--output",
            str(out),
            "--nrows",
            "200",
            "--quiet",
        ]
    )
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert len(text.splitlines()) >= 1


def test_plot_qq_cli_calls_plot_mqq(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from unittest.mock import MagicMock

    from gwaslab_cli import cli as cli_mod

    fake = MagicMock()
    monkeypatch.setattr(cli_mod, "load_sumstats", lambda *a, **k: fake)
    png = tmp_path / "qq.png"
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--plot",
            "qq",
            "--output",
            str(png),
            "--quiet",
        ]
    )
    fake.plot_mqq.assert_called_once()
    kw = fake.plot_mqq.call_args.kwargs
    assert kw.get("mode") == "qq"
    assert kw.get("save") == str(png)


def test_liftover_writes_output(tmp_path: Path) -> None:
    if not RAW.is_file():
        pytest.skip(f"missing {RAW}")
    from gwaslab_cli.cli import main as cli_main

    out = tmp_path / "lifted.tsv"
    cli_main(
        [
            "--input",
            str(RAW),
            "--liftover",
            "19",
            "38",
            "--output",
            str(out),
            "--nrows",
            "50",
            "--quiet",
        ]
    )
    assert out.with_suffix(out.suffix + ".gz").is_file() or out.is_file() or any(
        tmp_path.glob("lifted*")
    )


def test_infer_build_and_liftover_writes_output(tmp_path: Path) -> None:
    if not RAW.is_file():
        pytest.skip(f"missing {RAW}")
    from gwaslab_cli.cli import main as cli_main

    out = tmp_path / "inferred.tsv"
    cli_main(
        [
            "--input",
            str(RAW),
            "--qc",
            "--infer-build",
            "--liftover",
            "19",
            "38",
            "--output",
            str(out),
            "--nrows",
            "80",
            "--quiet",
        ]
    )
    matches = list(tmp_path.glob("inferred*"))
    assert matches


def test_plot_manhattan_cli_calls_plot_mqq(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from unittest.mock import MagicMock

    from gwaslab_cli import cli as cli_mod

    fake = MagicMock()
    monkeypatch.setattr(cli_mod, "load_sumstats", lambda *a, **k: fake)
    out = tmp_path / "m.png"
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--plot",
            "manhattan",
            "--output",
            str(out),
            "--quiet",
        ]
    )
    fake.plot_mqq.assert_called_once()
    kw = fake.plot_mqq.call_args.kwargs
    assert kw.get("mode") == "m"
    assert kw.get("save") == str(out)


def test_plot_regional_cli_calls_plot_mqq(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from unittest.mock import MagicMock

    from gwaslab_cli import cli as cli_mod

    fake = MagicMock()
    monkeypatch.setattr(cli_mod, "load_sumstats", lambda *a, **k: fake)
    out = tmp_path / "regional.png"
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--plot",
            "regional",
            "--chr",
            "1",
            "--start",
            "100",
            "--end",
            "200000",
            "--output",
            str(out),
            "--quiet",
        ]
    )
    fake.plot_mqq.assert_called_once()
    kw = fake.plot_mqq.call_args.kwargs
    assert kw.get("mode") == "r"
    assert kw.get("region") == (1, 100, 200000)


def test_plot_mqq_combined_cli_calls_plot_mqq(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from unittest.mock import MagicMock

    from gwaslab_cli import cli as cli_mod

    fake = MagicMock()
    monkeypatch.setattr(cli_mod, "load_sumstats", lambda *a, **k: fake)
    out = tmp_path / "mqq.png"
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--plot",
            "mqq",
            "--output",
            str(out),
            "--quiet",
        ]
    )
    fake.plot_mqq.assert_called_once()
    assert fake.plot_mqq.call_args.kwargs.get("mode") == "mqq"


def test_plot_forest_not_supported_as_cli_choice(tmp_path: Path) -> None:
    """Forest plots use gl.plot_forest(); --plot forest is not accepted."""
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(
            [
                "--input",
                str(tmp_path / "in.tsv"),
                "--plot",
                "forest",
                "--output",
                str(tmp_path / "f.png"),
            ]
        )
    assert exc.value.code == 2


def test_extract_novel_calls_get_novel(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from unittest.mock import MagicMock

    import pandas as pd

    from gwaslab_cli import cli as cli_mod

    fake = MagicMock()
    fake.get_novel.return_value = pd.DataFrame({"SNPID": ["x"]})
    monkeypatch.setattr(cli_mod, "load_sumstats", lambda *a, **k: fake)
    out = tmp_path / "novel.tsv"
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--get",
            "novel",
            "--output",
            str(out),
            "--quiet",
        ]
    )
    fake.get_novel.assert_called_once()
    assert out.is_file()


def test_ref_infer_passes_ref_alt_freq_to_harmonize(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Covers ref_alt_freq default vs explicit when --ref-infer is set."""
    from unittest.mock import MagicMock

    from gwaslab_cli import cli as cli_mod

    fake = MagicMock()
    monkeypatch.setattr(cli_mod, "load_sumstats", lambda *a, **k: fake)
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--harmonize",
            "--ref-infer",
            "/tmp/ref.vcf",
            "--quiet",
        ]
    )
    assert fake.harmonize.call_count == 1
    kw = fake.harmonize.call_args.kwargs
    assert kw["ref_alt_freq"] == "AF"

    fake.reset_mock()
    cli_mod.main(
        [
            "--input",
            "/dev/null/fake.tsv",
            "--harmonize",
            "--ref-infer",
            "/tmp/ref.vcf",
            "--ref-alt-freq",
            "AF_ALT",
            "--quiet",
        ]
    )
    kw = fake.harmonize.call_args.kwargs
    assert kw["ref_alt_freq"] == "AF_ALT"


def test_subcommand_version_matches_semver_line() -> None:
    from gwaslab_cli.main import main as launcher_main

    buf = io.StringIO()
    with redirect_stdout(buf):
        launcher_main(["version"])
    line = buf.getvalue().strip()
    assert re.search(r"\d+\.\d+", line)


def test_report_subcommand_rejects_non_html_pdf_output(tmp_path: Path) -> None:
    from gwaslab_cli.cli import main as cli_main

    with pytest.raises(SystemExit) as exc:
        cli_main(
            [
                "report",
                "-i",
                str(RAW),
                "-o",
                str(tmp_path / "out.txt"),
                "-q",
            ]
        )
    assert exc.value.code == 1


def test_report_subcommand_writes_html(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Smoke-test report CLI: HTML output only (no regional plots)."""
    monkeypatch.setenv("MPLBACKEND", "Agg")
    from gwaslab_cli.cli import main as cli_main

    inp = REPO_ROOT / "test" / "raw" / "realistic_sumstats.tsv"
    outp = tmp_path / "qc.html"
    cli_main(
        [
            "report",
            "-i",
            str(inp),
            "-o",
            str(outp),
            "--fmt",
            "auto",
            "--build",
            "19",
            "--sig-level",
            "1e-50",
            "--mqq-dpi",
            "72",
            "-q",
        ]
    )
    assert outp.is_file()
    head = outp.read_text(encoding="utf-8")[:200].lower()
    assert "<!doctype html>" in head or "<html" in head
