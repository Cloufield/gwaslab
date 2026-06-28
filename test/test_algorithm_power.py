"""Unit tests for gwaslab.algorithm.core.power helpers."""

from __future__ import annotations

import time

import numpy as np
import pytest

from gwaslab.algorithm.core import power as ap
from gwaslab.util.util_in_calculate_power import get_beta, get_beta_binary, get_power

# Regression baselines captured before util→algorithm delegation refactor.
_REGRESSION = {
    "q_power": 0.9999999999999705,
    "b_power": 0.488850823667201,
    "b_power_or_to_rr": 0.1862354472609722,
    "q_beta_first": (0.47959591836734694, 0.20417959183673468),
    "beta_rows": 48,
}


def test_resolve_genotype_rr_or_equals_rr_when_no_conversion():
    rr = ap.resolve_genotype_rr(genotype_or=1.5, prevalence=0.15, use_or_to_rr=False)
    assert rr == pytest.approx(1.5)


def test_resolve_genotype_rr_zhang_kai_conversion():
    or_val = 2.0
    prev = 0.15
    rr = ap.resolve_genotype_rr(genotype_or=or_val, prevalence=prev, use_or_to_rr=True)
    assert rr == pytest.approx(ap.or_to_rr(or_val, prev))
    assert rr != pytest.approx(or_val)


def test_resolve_genotype_rr_from_beta():
    beta = 0.2
    rr = ap.resolve_genotype_rr(beta=beta, prevalence=0.15, use_or_to_rr=False)
    assert rr == pytest.approx(np.exp(beta))


def test_power_quantitative_smoke_and_monotonicity():
    p_low = ap.power_quantitative(beta=0.1, eaf=0.3, n=5000)
    p_high = ap.power_quantitative(beta=0.3, eaf=0.3, n=5000)
    assert 0.0 <= p_low <= 1.0
    assert p_high > p_low

    p_small_n = ap.power_quantitative(beta=0.2, eaf=0.3, n=2000)
    p_large_n = ap.power_quantitative(beta=0.2, eaf=0.3, n=20000)
    assert p_large_n > p_small_n


def test_power_binary_smoke_and_monotonicity():
    kwargs = dict(daf=0.2, ncase=2000, ncontrol=15000, prevalence=0.15)
    p_low = ap.power_binary(genotype_rr=1.1, **kwargs)
    p_high = ap.power_binary(genotype_rr=1.5, **kwargs)
    assert 0.0 <= p_low <= 1.0
    assert p_high > p_low


def test_power_contour_at_target_quantitative():
    grid, eafs, betas = ap.power_grid_quantitative(
        eaf_range=(0.0001, 0.5),
        beta_range=(0.0001, 3.0),
        n=10000,
        n_matrix=50,
    )
    contour = ap.power_contour_at_target(grid, eafs, betas, target_power=0.8)
    assert len(contour) > 0
    assert all(len(pair) == 2 for pair in contour)


def test_infer_maf_and_beta_ranges():
    maf_range = ap.infer_maf_range(0.001)
    assert maf_range[0] <= 1e-3
    assert maf_range[1] == 0.5

    assert ap.infer_beta_range(2.0) == (0.0001, 3.0)
    assert ap.infer_beta_range(4.0) == (0.0001, 4.0)


def test_get_power_regression_quantitative():
    val = get_power(mode="q", beta=0.2, eaf=0.3, n=10000, verbose=False)
    assert val == pytest.approx(_REGRESSION["q_power"])


def test_get_power_regression_binary():
    kwargs = dict(
        mode="b",
        beta=0.2,
        daf=0.2,
        ncase=2000,
        ncontrol=15000,
        prevalence=0.15,
        verbose=False,
    )
    assert get_power(**kwargs) == pytest.approx(_REGRESSION["b_power"])
    assert get_power(or_to_rr=True, **kwargs) == pytest.approx(_REGRESSION["b_power_or_to_rr"])


def test_get_beta_regression_quantitative():
    df = get_beta(mode="q", t=0.8, n=10000, n_matrix=50, verbose=False)
    assert len(df) == _REGRESSION["beta_rows"]
    first = tuple(df.iloc[0].tolist())
    assert first[0] == pytest.approx(_REGRESSION["q_beta_first"][0], rel=1e-3)
    assert first[1] == pytest.approx(_REGRESSION["q_beta_first"][1], rel=1e-3)


def test_get_beta_binary_regression():
    df = get_beta_binary(
        t=0.8,
        ncase=2000,
        ncontrol=15000,
        prevalence=0.15,
        n_matrix=50,
        verbose=False,
    )
    assert len(df) == _REGRESSION["beta_rows"]


def test_power_grid_binary_matches_scalar_power_binary():
    grid, eafs, betas = ap.power_grid_binary(
        eaf_range=(0.0001, 0.5),
        beta_range=(0.0001, 3.0),
        ncase=2000,
        ncontrol=15000,
        prevalence=0.15,
        n_matrix=20,
    )
    for i, daf in enumerate(eafs):
        for j, beta in enumerate(betas):
            rr = ap.resolve_genotype_rr(beta=float(beta), prevalence=0.15, use_or_to_rr=False)
            expected = ap.power_binary(
                genotype_rr=rr,
                daf=float(daf),
                ncase=2000,
                ncontrol=15000,
                prevalence=0.15,
            )
            assert grid[i, j] == pytest.approx(expected)


def test_power_grid_binary_performance():
    start = time.perf_counter()
    ap.power_grid_binary(
        eaf_range=(0.0001, 0.5),
        beta_range=(0.0001, 3.0),
        ncase=2000,
        ncontrol=15000,
        prevalence=0.15,
        n_matrix=200,
    )
    elapsed = time.perf_counter() - start
    assert elapsed < 2.0, f"power_grid_binary(n_matrix=200) took {elapsed:.2f}s"
