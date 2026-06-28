"""CLI tests for ref add/remove and list --kind."""

from __future__ import annotations

import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def _ensure_src_path() -> None:
    src = str(REPO_ROOT / "src")
    if src not in sys.path:
        sys.path.insert(0, src)


def test_ref_add_and_remove(tmp_path, monkeypatch):
    import gwaslab as gl
    from gwaslab_cli.cli import main as cli_main

    data_file = tmp_path / "panel.vcf.gz"
    data_file.write_bytes(b"vcf")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"downloaded": {}}), encoding="utf-8")

    monkeypatch.setitem(gl.options.paths, "config", str(config))
    monkeypatch.setitem(gl.options.paths, "data_directory", str(tmp_path) + "/")

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["ref", "add", "test_panel", str(data_file)])
    assert gl.get_path("test_panel", verbose=False) == str(data_file.resolve())

    meta = json.loads(config.read_text())["downloaded"]["test_panel"]
    assert meta["kind"] == "ref"
    assert meta["source"] == "local"

    cli_main(["ref", "remove", "test_panel"])
    payload = json.loads(config.read_text())
    assert "test_panel" not in payload["downloaded"]
    assert data_file.exists()


def test_ref_remove_delete_file(tmp_path, monkeypatch):
    import gwaslab as gl
    from gwaslab_cli.cli import main as cli_main

    data_file = tmp_path / "panel.vcf.gz"
    data_file.write_bytes(b"vcf")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"downloaded": {}}), encoding="utf-8")

    monkeypatch.setitem(gl.options.paths, "config", str(config))
    monkeypatch.setitem(gl.options.paths, "data_directory", str(tmp_path) + "/")

    cli_main(["ref", "add", "test_panel", str(data_file)])
    cli_main(["ref", "remove", "test_panel", "--delete-file"])
    assert not data_file.exists()


def test_list_ref_kind_filter(monkeypatch):
    import gwaslab as gl
    from gwaslab_cli.cli import main as cli_main

    monkeypatch.setattr(
        gl,
        "check_downloaded_ref",
        lambda **k: {
            "1kg_eas_hg19": {"local_path": "/a", "kind": "ref"},
            "GCST1": {"local_path": "/b", "kind": "sumstats", "source": "gwas_catalog"},
        },
    )
    monkeypatch.setattr(gl, "check_available_ref", lambda **k: {})

    buf = io.StringIO()
    with redirect_stdout(buf):
        cli_main(["list", "ref", "--downloaded", "--kind", "sumstats", "--json"])
    data = json.loads(buf.getvalue())
    assert list(data["downloaded"].keys()) == ["GCST1"]
