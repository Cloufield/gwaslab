"""Shared demo data for panel visualization tests."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from gwaslab.viz.viz_aux_track_bundle import ContactBundle

_FIXTURES_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_DIR = os.path.join(_FIXTURES_DIR, "data")
_REPO_ROOT = os.path.abspath(os.path.join(_FIXTURES_DIR, "..", ".."))
_REF_DIR = os.path.join(_REPO_ROOT, "test", "ref")

REGION = (7, 126253550, 128253550)
_VCF_PREFIX = os.path.join(_REF_DIR, "1kg_eas_hg19.chr7_126253550_128253550")

_DEMO: dict | None = None


def _make_sumstats_from_bim(n: int = 80) -> pd.DataFrame:
    bim_path = f"{_VCF_PREFIX}.bim"
    rows = []
    with open(bim_path, encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            if i >= n:
                break
            chrom, snpid, _, pos, ea, nea = line.rstrip("\n").split("\t")[:6]
            pos_i = int(pos)
            if pos_i < REGION[1] or pos_i > REGION[2]:
                continue
            p = 10 ** (-(4 + (i % 6)))
            rows.append(
                {
                    "CHR": int(chrom),
                    "POS": pos_i,
                    "SNPID": snpid,
                    "rsID": snpid,
                    "EA": ea,
                    "NEA": nea,
                    "P": p,
                    "MLOG10P": -np.log10(p),
                    "BETA": 0.05 * ((i % 3) - 1),
                    "SE": 0.01,
                }
            )
    if not rows:
        raise RuntimeError(f"No variants from {_VCF_PREFIX}.bim overlap {REGION}")
    return pd.DataFrame(rows)


def _make_pipcs(sumstats: pd.DataFrame) -> pd.DataFrame:
    pipcs = sumstats[["CHR", "POS"]].copy()
    pipcs["PIP"] = np.linspace(0.2, 0.99, len(pipcs))
    pipcs["CREDIBLE_SET_INDEX"] = np.where(pipcs["PIP"] > 0.8, 1, 0)
    pipcs["CS_CATEGORY"] = np.where(
        pipcs["CREDIBLE_SET_INDEX"] > 0, "CS1", "None"
    )
    return pipcs


def _make_ag_contact() -> ContactBundle:
    resolution = 10_000
    start, end = REGION[1], REGION[2]
    n_bins = max(2, (end - start) // resolution)
    rng = np.random.RandomState(0)
    values = rng.random((n_bins, n_bins, 1)) * 0.5
    for i in range(n_bins):
        values[i, i, 0] = 1.0
    metadata = pd.DataFrame({"name": ["Hi-C"]})
    return ContactBundle(
        values=values,
        metadata=metadata,
        resolution=resolution,
        region=REGION,
    )


def setup_demo() -> dict:
    """Return paths and tables for stacked-panel visualization tests."""
    global _DEMO
    if _DEMO is not None:
        return _DEMO

    sumstats = _make_sumstats_from_bim()
    _DEMO = {
        "region": REGION,
        "gtf": os.path.join(_DATA_DIR, "demo_panel.gtf"),
        "bedpe": os.path.join(_DATA_DIR, "demo_panel.bedpe"),
        "chromatin_bed": os.path.join(_DATA_DIR, "demo_panel_chromatin.bed"),
        "pipcs": _make_pipcs(sumstats),
        "sumstats": sumstats,
        "ld_sumstats": sumstats.copy(),
        "vcf": f"{_VCF_PREFIX}.vcf.gz",
        "ag_contact": _make_ag_contact(),
    }
    return _DEMO
