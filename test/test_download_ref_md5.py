"""Tests for download_ref md5sum verification."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from gwaslab.bd.bd_download import download_ref


def test_download_ref_checks_md5sum_key():
    catalog = {
        "demo_ref": {
            "url": "http://example.com/demo.txt.gz",
            "md5sum": "abc123",
        }
    }
    log = MagicMock()

    with patch("gwaslab.bd.bd_download.check_available_ref", return_value=catalog):
        with patch("gwaslab.bd.bd_download.download_file"):
            with patch("gwaslab.bd.bd_download.search_local", return_value=False):
                with patch("gwaslab.bd.bd_download.check_file_integrity", return_value=0) as verify:
                    with patch("gwaslab.bd.bd_download.update_record") as update:
                        result = download_ref("demo_ref", directory="/tmp/", log=log)
                        verify.assert_called_once()
                        assert verify.call_args.kwargs["md5sum"] == "abc123"
                        update.assert_not_called()
                        assert result is None
