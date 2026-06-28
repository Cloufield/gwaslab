"""Tests for shared bd_io download utilities."""

from __future__ import annotations

import hashlib
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gwaslab.bd.bd_io import stream_download, verify_md5


def test_stream_download_atomic_replace(tmp_path):
    dest = tmp_path / "file.bin"
    part = tmp_path / "file.bin.part"

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size=0):
            yield b"hello"

    with patch("gwaslab.bd.bd_io.requests.get", return_value=FakeResponse()):
        out = stream_download("http://example.com/file.bin", str(dest))

    assert out == str(dest)
    assert dest.read_bytes() == b"hello"
    assert not part.exists()


def test_verify_md5(tmp_path):
    target = tmp_path / "data.txt"
    target.write_text("abc", encoding="utf-8")
    digest = hashlib.md5(b"abc").hexdigest()
    log = MagicMock()
    assert verify_md5(str(target), digest, log=log) is True
    assert verify_md5(str(target), "deadbeef", log=log) is False
