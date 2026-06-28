"""Shared download and integrity-check utilities for reference and sumstats modules.
"""

from __future__ import annotations

import hashlib
import os
from typing import Optional

import requests

from gwaslab.info.g_Log import Log


def ensure_directory(dir_path: str) -> str:
    """Create directory if missing and return the path.
"""
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


def stream_download(
    url: str,
    output_path: str,
    timeout: int = 60,
    chunk_size: int = 1024 * 1024,
    overwrite: bool = True,
) -> str:
    """Download URL content to output_path with atomic rename via a .part temp file.

    Returns the output_path on success.
"""
    if not overwrite and os.path.exists(output_path):
        return output_path

    tmp_path = output_path + ".part"
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with open(tmp_path, "wb") as handle:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    handle.write(chunk)
    os.replace(tmp_path, output_path)
    return output_path


def verify_md5(local_path: str, expected_md5: str, log: Optional[Log] = None) -> bool:
    """Return True when file MD5 matches expected_md5.
"""
    md5_hash = hashlib.md5()
    with open(local_path, "rb") as handle:
        for byte_block in iter(lambda: handle.read(4096 * 1000), b""):
            md5_hash.update(byte_block)
    digest = md5_hash.hexdigest()
    if log is not None:
        log.write(f" -File path: {local_path}")
        log.write(f" -MD5 check: {digest}")
    if digest == expected_md5:
        if log is not None:
            log.write(" -MD5 verified.")
        return True
    if log is not None:
        log.warning("-MD5 VERIFICATION FAILED!")
    return False
