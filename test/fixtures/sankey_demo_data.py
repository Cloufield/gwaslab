"""Offline fixtures for the plot_sankey development tutorial."""

from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd

OVERALL_SIGNAL_ORDER = [
    "Overall GW significant",
    "Overall suggestive",
    "Overall non-significant",
]

SUBTYPE_SIGNAL_ORDER = [
    "Multiple subtypes",
    "Subtype A only",
    "Subtype B only",
    "Subtype C only",
    "Subtype suggestive",
    "No subtype signal",
]


def simulate_disease_subtype_sumstats(
    n_variants: int = 6000,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Simulate a multi-trait GWAS table for Sankey demos.

    Each variant has an overall-disease P-value, subtype P-values (A/B/C),
    and derived categorical labels ``overall_signal`` and ``subtype_signal``.
    Use with ``columns=["overall_signal", "subtype_signal", "MAF"]``.
    """
    rng = np.random.default_rng(seed)

    eaf = rng.beta(0.35, 3.5, size=n_variants)
    beta = rng.normal(0, 0.12, size=n_variants)
    p_overall = 10 ** (-rng.uniform(4, 11, size=n_variants))
    p_subtypes: Dict[str, np.ndarray] = {
        subtype: 10 ** (-rng.uniform(4, 11, size=n_variants))
        for subtype in ("A", "B", "C")
    }

    strata = rng.choice(
        ["shared", "overall_only", "subtype_only", "null"],
        size=n_variants,
        p=[0.04, 0.03, 0.05, 0.88],
    )

    for i, label in enumerate(strata):
        if label == "shared":
            p_overall[i] = 10 ** (-rng.uniform(8.5, 14))
            for subtype in p_subtypes:
                if rng.random() < 0.7:
                    p_subtypes[subtype][i] = 10 ** (-rng.uniform(8.0, 13))
        elif label == "overall_only":
            p_overall[i] = 10 ** (-rng.uniform(8.5, 14))
        elif label == "subtype_only":
            hit = rng.choice(["A", "B", "C"])
            p_subtypes[hit][i] = 10 ** (-rng.uniform(8.5, 14))

    overall_signal = [_classify_overall_signal(p) for p in p_overall]
    subtype_signal = [
        _classify_subtype_signal({name: p_subtypes[name][i] for name in p_subtypes})
        for i in range(n_variants)
    ]

    return pd.DataFrame(
        {
            "SNPID": [f"rs{i}" for i in range(n_variants)],
            "CHR": rng.integers(1, 23, size=n_variants),
            "POS": rng.integers(1_000_000, 100_000_000, size=n_variants),
            "EA": "A",
            "NEA": "G",
            "EAF": eaf,
            "BETA": beta,
            "P": p_overall,
            "P_SUB_A": p_subtypes["A"],
            "P_SUB_B": p_subtypes["B"],
            "P_SUB_C": p_subtypes["C"],
            "overall_signal": overall_signal,
            "subtype_signal": subtype_signal,
        }
    )


def _classify_overall_signal(p_value: float) -> str:
    if p_value < 5e-8:
        return "Overall GW significant"
    if p_value < 5e-6:
        return "Overall suggestive"
    return "Overall non-significant"


def _classify_subtype_signal(p_values: Dict[str, float]) -> str:
    hits: List[str] = [name for name, p in p_values.items() if p < 5e-8]
    if len(hits) >= 2:
        return "Multiple subtypes"
    if len(hits) == 1:
        return f"Subtype {hits[0]} only"
    if any(p < 5e-6 for p in p_values.values()):
        return "Subtype suggestive"
    return "No subtype signal"


def overall_signal_colors() -> Dict[str, str]:
    return {
        "Overall GW significant": "#CB181D",
        "Overall suggestive": "#FDAE6B",
        "Overall non-significant": "#D9D9D9",
    }
