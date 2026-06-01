"""Builtin and downloaded paths for HapMap3 ancestry EAF panels."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from gwaslab.bd.bd_download import get_path
from gwaslab.info.g_Log import Log

_EAF_KEYWORDS = {
    "19": "1kg_hm3_hg19_eaf",
    "38": "1kg_hm3_hg38_eaf",
}


def builtin_eaf_path(build: str) -> Path:
    """Return packaged core EAF path for genome build 19 or 38."""
    b = str(build).replace("hg", "").replace("GRCh", "")
    if b not in ("19", "38"):
        raise ValueError(f"build must be '19' or '38', got {build!r}")
    base = Path(__file__).resolve().parents[1] / "data" / "hapmap3_EAF"
    return base / f"PAN.hapmap3.hg{b}.EAF.core.tsv.gz"


def resolve_ancestry_af(
    ancestry_af: Optional[str],
    build: Optional[str],
    log: Log = Log(),
    verbose: bool = True,
) -> str:
    """
    Resolve allele-frequency table path for infer_ancestry.

    Priority:
    1. Explicit file path (existing file)
    2. Keyword or None with build: downloaded full PAN if present, else builtin core
    """
    if build is None:
        raise ValueError("build is required when ancestry_af is None or a download keyword")

    b = str(build).replace("hg", "").replace("GRCh", "")
    if b not in ("19", "38"):
        raise ValueError(f"build must be '19' or '38', got {build!r}")

    keyword = _EAF_KEYWORDS[b]

    if ancestry_af is not None and ancestry_af not in (keyword, f"1kg_hm3_hg{b}_eaf"):
        path = Path(ancestry_af)
        if path.is_file():
            return str(path.resolve())
        raise ValueError(f"ancestry_af file not found: {ancestry_af}")

    downloaded = get_path(keyword, log=log, verbose=verbose)
    if downloaded is not False and downloaded is not None:
        log.write(f"  -Using downloaded reference: {keyword}", verbose=verbose)
        return downloaded

    core = builtin_eaf_path(b)
    if not core.is_file():
        raise ValueError(
            f"Builtin core EAF not found at {core}. "
            f"Reinstall gwaslab or download full panel: download_ref('{keyword}')"
        )
    log.write(
        f"  -Using builtin core EAF panel ({core.name}; ~30k SNPs). "
        f"For full HapMap3 panel: download_ref('{keyword}')",
        verbose=verbose,
    )
    return str(core)
