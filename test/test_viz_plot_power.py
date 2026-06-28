"""Smoke tests for plot_power registry integration."""

import matplotlib
matplotlib.use("Agg")

import gwaslab as gl


def test_plot_power_quantitative_smoke():
    fig = gl.plot_power(mode="q", ns=500, n_matrix=50, verbose=False)
    assert fig is not None


def test_plot_power_doc_has_parameters():
    doc = gl.plot_power.__doc__ or ""
    assert "Parameters" in doc
    assert "ns" in doc
