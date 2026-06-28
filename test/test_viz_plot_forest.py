"""Smoke tests for plot_forest registry filtering."""

import pandas as pd
import pytest

import gwaslab as gl


@pytest.fixture
def forest_df():
    return pd.DataFrame({
        "study": ["A", "B"],
        "beta": [0.1, 0.2],
        "se": [0.05, 0.06],
        "group": ["G1", "G1"],
    })


def test_plot_forest_accepts_column_mapping(forest_df):
    fig, _ = gl.plot_forest(
        forest_df,
        study_col="study",
        beta_col="beta",
        se_col="se",
        group_col="group",
        verbose=False,
    )
    assert fig is not None
