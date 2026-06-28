"""Tests for registry kind, filtering, get_path, check_and_download, scan."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gwaslab.bd.bd_download import (
    backfill_registry_kinds,
    check_and_download,
    filter_downloaded_registry,
    get_path,
    infer_registry_kind,
    scan_downloaded_files,
)


def test_infer_registry_kind():
    assert infer_registry_kind("1kg_eas_hg19") == "ref"
    assert infer_registry_kind("GCST90270926") == "sumstats"
    assert infer_registry_kind("GCST90270926", {"source": "gwas_catalog"}) == "sumstats"
    assert infer_registry_kind("custom", {"kind": "ref"}) == "ref"


def test_filter_downloaded_registry():
    entries = {
        "1kg_eas_hg19": {"local_path": "/a", "kind": "ref", "source": "catalog"},
        "GCST123": {"local_path": "/b", "kind": "sumstats", "source": "gwas_catalog"},
    }
    assert len(filter_downloaded_registry(entries, kind="ref")) == 1
    assert len(filter_downloaded_registry(entries, source="gwas_catalog")) == 1


def test_backfill_registry_kinds():
    cfg = {"downloaded": {"GCST999": {"local_path": "/x"}}}
    assert backfill_registry_kinds(cfg) is True
    assert cfg["downloaded"]["GCST999"]["kind"] == "sumstats"


def test_get_path_raise_on_missing(tmp_path, monkeypatch):
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"downloaded": {}}), encoding="utf-8")
    monkeypatch.setitem(__import__("gwaslab.bd.bd_config", fromlist=["options"]).options.paths, "config", str(config))

    assert get_path("missing_kw", verbose=False) is False
    with pytest.raises(KeyError):
        get_path("missing_kw", verbose=False, raise_on_missing=True)


def test_check_and_download_raises_on_failure(monkeypatch):
    monkeypatch.setattr(
        "gwaslab.bd.bd_download.get_path",
        lambda *a, **k: False,
    )
    monkeypatch.setattr(
        "gwaslab.bd.bd_download.get_default_directory",
        lambda: "/tmp/gwaslab_test/",
    )
    monkeypatch.setattr(
        "gwaslab.bd.bd_download.download_ref",
        lambda *a, **k: None,
    )
    with pytest.raises(RuntimeError, match="Failed to download"):
        check_and_download("demo_ref", log=MagicMock())


def test_scan_registers_source_and_kind(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    vcf = data_dir / "demo.vcf.gz"
    vcf.write_bytes(b"x")

    config = tmp_path / "config.json"
    config.write_text(json.dumps({"downloaded": {}}), encoding="utf-8")

    import gwaslab.bd.bd_config as bd_config

    monkeypatch.setitem(bd_config.options.paths, "config", str(config))
    monkeypatch.setitem(bd_config.options.paths, "data_directory", str(data_dir) + "/")

    catalog = {
        "demo_ref": {
            "url": "http://example.com/demo.vcf.gz",
            "description": "demo",
        }
    }
    log = MagicMock()

    with patch("gwaslab.bd.bd_download.check_available_ref", return_value=catalog):
        with patch("gwaslab.bd.bd_download.check_downloaded_ref", return_value={}):
            assert scan_downloaded_files(log=log, verbose=False, directory=str(data_dir)) is True

    payload = json.loads(config.read_text(encoding="utf-8"))
    assert payload["downloaded"]["demo_ref"]["source"] == "catalog"
    assert payload["downloaded"]["demo_ref"]["kind"] == "ref"
    assert payload["downloaded"]["demo_ref"]["local_path"] == str(vcf)
