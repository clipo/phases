"""Tests for the transmission-layer coupling of ceramic-style assortment to the
parent monument-mls emergence phi.

The chosen bistable operating point is (sigma=0.5, lambda_W=0.5), whose interior
saddle phi_star ~= 0.4299 lies safely inside (0,1) so that phi_star +/- 0.1 are
both valid initial conditions on opposite sides of the saddle.
"""
from __future__ import annotations

import numpy as np

from mls_emergence.signatures import convergence
from mls_emergence.transmission.model import (
    SIGMA,
    LAMBDA_W,
    coupling_robustness,
    emit_signatures,
    phi_trajectory,
    simulate_copying,
)

from signaling.emergence import phi_star


def _coords(n_groups: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    # A few well-separated clusters so the bounded spatial rule has structure.
    offsets = np.array([[0, 0], [10, 10], [0, 10], [10, 0]], float)
    labels = np.repeat(np.arange(4), n_groups // 4 + 1)[:n_groups]
    return offsets[labels] + rng.normal(0, 0.4, size=(n_groups, 2))


def test_phi_trajectory_bistable():
    ps = phi_star(SIGMA, LAMBDA_W)
    assert np.isfinite(ps) and 0.0 < ps < 1.0

    up = phi_trajectory(SIGMA, LAMBDA_W, ps + 0.1)
    dn = phi_trajectory(SIGMA, LAMBDA_W, ps - 0.1)

    assert up[-1] > 0.5, f"above-saddle start should rise toward 1, got {up[-1]}"
    assert dn[-1] < 0.5, f"below-saddle start should fall toward 0, got {dn[-1]}"


def test_coupling_zero_flat_one_rises():
    ps = phi_star(SIGMA, LAMBDA_W)
    phi_t = phi_trajectory(SIGMA, LAMBDA_W, ps + 0.1)  # rising trajectory
    coords = _coords(12)

    def trend(coupling):
        slices = simulate_copying(
            phi_t,
            coupling=coupling,
            n_groups=12,
            n_per_group=300,
            n_types=10,
            coords=coords,
            seed=7,
        )
        panel = emit_signatures(slices, coords)
        return convergence.time_derivative(convergence.convergence_score(panel))

    t0 = trend(0.0)
    t1 = trend(1.0)

    # coupling=0: assortment is not driven by phi, so the convergence score
    # should not show a meaningful rising ordinal trend.
    assert abs(t0) < 0.15, f"coupling=0 trend should be ~flat, got {t0}"
    # coupling=1: phi drives assortment, signatures co-rise.
    assert t1 > 0.15, f"coupling=1 trend should be clearly positive, got {t1}"
    assert t1 - t0 > 0.15, f"coupling=1 should exceed coupling=0 by a clear margin: {t1} vs {t0}"


def test_coupling_robustness_monotone_ish():
    ps = phi_star(SIGMA, LAMBDA_W)
    coords = _coords(12)
    df = coupling_robustness(
        SIGMA,
        LAMBDA_W,
        phi_0=ps + 0.1,
        couplings=[0.0, 0.25, 0.5, 0.75, 1.0],
        n_groups=12,
        n_per_group=300,
        n_types=10,
        coords=coords,
        seed=7,
    )
    assert list(df["coupling"]) == [0.0, 0.25, 0.5, 0.75, 1.0]
    top = df.loc[df["coupling"] == 1.0, "convergence_trend"].iloc[0]
    bot = df.loc[df["coupling"] == 0.0, "convergence_trend"].iloc[0]
    assert top > bot, f"top coupling trend {top} should exceed bottom {bot}"
    # Generally increasing: the trend at the top half should exceed the bottom half.
    lower = df[df["coupling"] <= 0.25]["convergence_trend"].mean()
    upper = df[df["coupling"] >= 0.75]["convergence_trend"].mean()
    assert upper > lower
