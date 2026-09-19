"""Tests for the validation harness (signatures_over_axis, run_blind)."""
from __future__ import annotations

import pandas as pd

from mls_emergence.validation.harness import run_blind, signatures_over_axis
from mls_emergence.validation.mechanisms import (
    T_DEFAULT,
    gen_aggregated_signaling,
    gen_drift_space,
    gen_group_emergence,
    gen_patchiness,
)

SIGNATURE_COLUMNS = ["neutral_departure", "seriability", "fst", "spatial_boundary"]


def test_signatures_over_axis_shape():
    slices, coords = gen_group_emergence(seed=1)
    panel = signatures_over_axis(slices, coords)
    assert isinstance(panel, pd.DataFrame)
    assert list(panel.columns) == SIGNATURE_COLUMNS
    assert len(panel) == T_DEFAULT
    assert list(panel.index) == list(range(T_DEFAULT))
    assert panel.notna().all().all()


def test_run_blind_one_panel_per_mechanism():
    generators = {
        "group_emergence": gen_group_emergence,
        "aggregated_signaling": gen_aggregated_signaling,
        "patchiness": gen_patchiness,
        "drift_space": gen_drift_space,
    }
    panels = run_blind(generators, seed=7)
    assert set(panels) == set(generators)
    for name, panel in panels.items():
        assert isinstance(panel, pd.DataFrame)
        assert list(panel.columns) == SIGNATURE_COLUMNS
        assert len(panel) == T_DEFAULT
