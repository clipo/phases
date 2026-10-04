"""Mechanism-faithful contracts for the four generators.

Each test encodes the qualitative behaviour the generator MUST exhibit for the
validation to mean anything. If a mechanism cannot meet its contract honestly,
that is a finding to report, not a contract to relax.
"""
from __future__ import annotations

import numpy as np

from mls_emergence.signatures import neutral, variance
from mls_emergence.signatures.assortativity import spatial_assortativity
from mls_emergence.validation.mechanisms import (
    gen_aggregated_signaling,
    gen_drift_space,
    gen_group_emergence,
    gen_patchiness,
    simulate_slice,
)


def _within_conformity_departure(slice_counts: np.ndarray) -> float:
    """Mean within-group signed conformity departure (1 - theta_f/theta_e).

    Positive => conformist (fewer effective types than Ewens neutrality predicts).
    Measured within groups and averaged: conformity is a within-group bias, and
    pooling divergent groups cancels it (see mechanisms.py rationale).
    """
    vals = []
    for row in slice_counts:
        tf = neutral.theta_f(row)
        te = neutral.theta_e(row)
        if np.isfinite(tf) and te > 0:
            vals.append(1.0 - tf / te)
    return float(np.mean(vals)) if vals else 0.0


def test_simulate_slice_shape_and_counts():
    rng = np.random.default_rng(0)
    coords = rng.random((6, 2))
    m = simulate_slice(
        n_groups=6,
        n_per_group=200,
        n_types=8,
        between_divergence=0.5,
        within_conformity=0.3,
        spatial_rule="bounded",
        coords=coords,
        rng=rng,
    )
    assert m.shape == (6, 8)
    assert m.dtype.kind in "iu"
    # Each group draws exactly n_per_group items.
    assert np.all(m.sum(axis=1) == 200)


def test_simulate_slice_panmictic_low_fst():
    """d=0 => all groups share one pool => F_ST near zero."""
    rng = np.random.default_rng(1)
    coords = rng.random((12, 2))
    m = simulate_slice(
        n_groups=12,
        n_per_group=300,
        n_types=10,
        between_divergence=0.0,
        within_conformity=0.0,
        spatial_rule="none",
        coords=coords,
        rng=rng,
    )
    assert variance.cultural_fst(m) < 0.05


def test_group_emergence_contract():
    slices, coords = gen_group_emergence(seed=42)
    fst0 = variance.cultural_fst(slices[0])
    fstT = variance.cultural_fst(slices[-1])
    assert fstT > fst0 + 0.1
    assert _within_conformity_departure(slices[-1]) > 0.1


def test_aggregated_signaling_contract():
    slices, coords = gen_aggregated_signaling(seed=42)
    max_fst = max(variance.cultural_fst(s) for s in slices)
    assert max_fst < 0.05
    assert _within_conformity_departure(slices[-1]) > 0.1


def test_patchiness_contract():
    slices, coords = gen_patchiness(seed=42)
    fsts = [variance.cultural_fst(s) for s in slices]
    # Flat in expectation (fixed latent process); the spread is sampling noise at
    # N=300/group and is far below the emergence span (~0.0 -> ~0.99).
    assert (max(fsts) - min(fsts)) < 0.05  # flat over time
    r0 = spatial_assortativity(slices[0], coords)["mantel_r"]
    assert r0 < -0.1  # some spatial structure present


def test_drift_space_contract():
    slices, coords = gen_drift_space(seed=42)
    for s in slices:
        assert _within_conformity_departure(s) < 0.4  # near-neutral throughout
    r0 = spatial_assortativity(slices[0], coords)["mantel_r"]
    assert r0 < -0.1  # smooth decay present
