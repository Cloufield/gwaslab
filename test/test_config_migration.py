"""Tests for config path resolution and legacy registry migration."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gwaslab.bd import bd_config


def test_resolve_default_paths_uses_home_config(tmp_path, monkeypatch):
    data_dir = tmp_path / "gwaslab_home"
    monkeypatch.setenv("GWASLAB_DATA_DIR", str(data_dir) + "/")
    monkeypatch.delenv("GWASLAB_CONFIG", raising=False)
    paths = bd_config._resolve_default_paths()
    assert paths["data_directory"] == str(data_dir) + "/"
    assert paths["config"] == str(data_dir / "config.json")


def test_migrate_legacy_config(tmp_path, monkeypatch):
    legacy_path = tmp_path / "legacy_config.json"
    legacy_path.write_text(
        json.dumps(
            {"downloaded": {"test_ref": {"local_path": "/tmp/test_ref.vcf.gz"}}},
            indent=4,
        ),
        encoding="utf-8",
    )
    new_config = tmp_path / "new" / "config.json"
    monkeypatch.setattr(bd_config, "_LEGACY_CONFIG", legacy_path)

    migrated = bd_config._migrate_legacy_config(str(new_config))
    assert migrated is True
    payload = json.loads(new_config.read_text(encoding="utf-8"))
    assert "test_ref" in payload["downloaded"]


def test_settings_persist_data_directory(tmp_path, monkeypatch):
    data_dir = tmp_path / "home"
    data_dir.mkdir()
    custom = str(tmp_path / "other_refs") + "/"
    bd_config.save_settings({"data_directory": custom}, data_dir=str(data_dir))
    loaded = bd_config.load_settings(str(data_dir))
    assert loaded.get("data_directory") == custom


def test_ensure_user_layout_creates_config(tmp_path, monkeypatch):
    data_dir = tmp_path / "layout_home"
    config_path = data_dir / "config.json"
    monkeypatch.setitem(bd_config.options.paths, "data_directory", str(data_dir) + "/")
    monkeypatch.setitem(bd_config.options.paths, "config", str(config_path))
    if config_path.exists():
        config_path.unlink()
    bd_config.ensure_user_layout()
    assert config_path.is_file()
