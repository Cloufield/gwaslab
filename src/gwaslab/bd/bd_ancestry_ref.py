"""Resolve HapMap3 / 1KG effective-allele-frequency (EAF) panels for ancestry inference.

Used by :func:`gwaslab.util.util_ex_infer_ancestry._infer_ancestry` to pick a reference
table before Fst-based ancestry assignment. Two panel tiers exist:

- **Core** (~30k SNPs): shipped under ``gwaslab/data/hapmap3_EAF/`` as
  ``PAN.hapmap3.hg{19,38}.EAF.core.tsv.gz``.
- **Full** (~1.2M SNPs): optional download via ``download_ref("1kg_hm3_hg{19,38}_eaf")``,
  registered in the local gwaslab config and resolved through :func:`get_path`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from gwaslab.bd.bd_download import get_path
from gwaslab.info.g_Log import Log

# Maps normalized build ("19" / "38") to download_ref / get_path keywords.
_EAF_KEYWORDS = {
    "19": "1kg_hm3_hg19_eaf",
    "38": "1kg_hm3_hg38_eaf",
}


def builtin_eaf_path(build: str) -> Path:
    """Path to the packaged core HapMap3 EAF panel for a genome build.

    Parameters
    ----------
    build : str
        Genome build: ``"19"``, ``"38"``, or common aliases (``hg19``, ``GRCh38``, …).

    Returns
    -------
    pathlib.Path
        ``.../data/hapmap3_EAF/PAN.hapmap3.hg{b}.EAF.core.tsv.gz`` (existence not checked).

    Raises
    ------
    ValueError
        If ``build`` cannot be normalized to 19 or 38.
    """
    b = str(build).replace("hg", "").replace("GRCh", "")
    if b not in ("19", "38"):
        raise ValueError(f"build must be '19' or '38', got {build!r}")
    base = Path(__file__).resolve().parents[1] / "data" / "hapmap3_EAF"
    return base / f"PAN.hapmap3.hg{b}.EAF.core.tsv.gz"


def resolve_ancestry_af(
    ancestry_af: Optional[str],
    build: Optional[str],
    _core: bool = False,
    log: Log = Log(),
    verbose: bool = True,
) -> str:
    """Resolve the EAF TSV path used by :func:`~gwaslab.util.util_ex_infer_ancestry._infer_ancestry`.

    Resolution order when ``ancestry_af`` is ``None`` or a download keyword
    (``1kg_hm3_hg19_eaf`` / ``1kg_hm3_hg38_eaf``):

    1. Downloaded full panel in ``~/.gwaslab/`` if :func:`get_path` finds it.
    2. Builtin core panel under ``gwaslab/data/hapmap3_EAF/``.

    Any other non-keyword string is treated as a user-supplied file path and must exist.

    Parameters
    ----------
    ancestry_af : str or None
        Explicit path, a download keyword, or ``None`` (keyword implied from ``build``).
    build : str or None
        Required unless ``ancestry_af`` is an existing file path. Normalized like
        :func:`builtin_eaf_path`.
    _core : bool, optional
        Internal flag. If ``True``, skip downloaded panel lookup and force builtin core.
    log : Log, optional
        Logger for which panel tier was chosen.
    verbose : bool, optional
        Whether to write resolution messages to ``log``.

    Returns
    -------
    str
        Resolved path to a tab-separated EAF file.

    Raises
    ------
    ValueError
        Missing ``build``, invalid build, missing custom file, or missing builtin core.
    """
    if build is None:
        raise ValueError("build is required when ancestry_af is None or a download keyword")

    b = str(build).replace("hg", "").replace("GRCh", "")
    if b not in ("19", "38"):
        raise ValueError(f"build must be '19' or '38', got {build!r}")

    keyword = _EAF_KEYWORDS[b]

    # Custom path: not None and not the canonical keyword for this build.
    if ancestry_af is not None and ancestry_af not in (keyword, f"1kg_hm3_hg{b}_eaf"):
        path = Path(ancestry_af)
        if path.is_file():
            return str(path.resolve())
        raise ValueError(f"ancestry_af file not found: {ancestry_af}")

    # Prefer user-downloaded full PAN unless _core explicitly forces builtin core.
    if not _core:
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
