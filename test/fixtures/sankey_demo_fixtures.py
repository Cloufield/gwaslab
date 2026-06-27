"""Demo data for disease-subtype Sankey visualization tests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gwaslab.viz.viz_aux_sankey import _categorize_maf


def overall_signal_colors() -> dict[str, str]:
    return {
        "Overall GW significant": "#E51819",
        "Not overall GW significant": "#BBBBBB",
        "Subtype A signal": "#1f77b4",
        "Subtype B signal": "#2ca02c",
        "No subtype signal": "#DDDDDD",
    }


def simulate_disease_subtype_sumstats(
    n_variants: int = 500,
    seed: int = 0,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    eaf = rng.uniform(0.0005, 0.5, size=n_variants)
    p_overall = 10 ** (-rng.uniform(2, 12, size=n_variants))
    p_subtype = 10 ** (-rng.uniform(2, 10, size=n_variants))

    overall_sig = p_overall < 5e-8
    if not overall_sig.any():
        overall_sig[: max(5, n_variants // 20)] = True
        p_overall[: max(5, n_variants // 20)] = 1e-10

    subtype_a = (p_subtype < 1e-5) & overall_sig
    subtype_b = (p_subtype < 1e-4) & overall_sig & ~subtype_a

    overall_signal = np.where(
        overall_sig,
        "Overall GW significant",
        "Not overall GW significant",
    )
    subtype_signal = np.where(
        subtype_a,
        "Subtype A signal",
        np.where(subtype_b, "Subtype B signal", "No subtype signal"),
    )

    df = pd.DataFrame(
        {
            "EAF": eaf,
            "P": p_overall,
            "P_SUBTYPE": p_subtype,
            "BETA": rng.normal(0, 0.2, size=n_variants),
            "overall_signal": overall_signal,
            "subtype_signal": subtype_signal,
        }
    )
    df["MAF"] = _categorize_maf(df["EAF"])
    return df
