"""Smoke tests for plot_associations."""

import matplotlib
matplotlib.use("Agg")

import pandas as pd

from gwaslab.g_Sumstats import Sumstats


def test_plot_associations_smoke():
    df = pd.DataFrame({
        "rsID": ["rs1", "rs1", "rs2"],
        "reported_trait_GCV2": ["TraitA", "TraitB", "TraitA"],
        "BETA_GCV2": [0.1, -0.2, 0.05],
        "P_GCV2": [1e-8, 1e-6, 1e-5],
    })
    ss = Sumstats(df[["rsID"]], fmt="gwaslab")
    ss.associations = df
    fig, ax = ss.plot_associations(verbose=False)
    assert fig is not None
    assert ax is not None
