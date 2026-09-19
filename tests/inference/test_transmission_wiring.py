import importlib

import numpy as np
import pytest


@pytest.fixture(scope="module")
def wiring():
    return importlib.import_module("analyses.38_abc_smc_transmission")


def test_prior_bounds_and_logpdf(wiring):
    assert wiring.PRIOR_LO.shape == (4,)
    assert wiring.PRIOR_HI.shape == (4,)
    inside = (wiring.PRIOR_LO + wiring.PRIOR_HI) / 2
    assert np.isfinite(wiring.prior_logpdf(inside))
    outside = wiring.PRIOR_HI + 1.0
    assert wiring.prior_logpdf(outside) == -np.inf


def test_fit_target_runs_on_synthetic_obs(wiring):
    rng = np.random.default_rng(0)
    K = 4
    # small synthetic 6-bin count matrix
    obs = rng.integers(5, 40, size=(wiring.__dict__["N_BINS_LOCAL"], K)).astype(float)
    res, adj_b, adj_joint = wiring.fit_target(obs, slice(None), n_particles=60,
                                              n_rounds=2, seed=1)
    assert adj_b.shape == (60,)
    assert np.isclose(res.weights.sum(), 1.0)
    # b is column 1 of theta; adjusted b should stay within a sane range
    assert np.all(adj_b > -2.0) and np.all(adj_b < 2.0)

    # The joint posterior is returned alongside the b marginal (added
    # 2026-09-02 so the neutrality posterior predictive check can draw the
    # population size). Every adjusted particle must be inside the prior box:
    # the unbounded local-linear adjustment used to put 300/800 innovation
    # rates and 150/800 population sizes below zero (F21).
    assert adj_joint.shape == (60, 4)
    assert np.all(adj_joint >= wiring.PRIOR_LO) and np.all(adj_joint <= wiring.PRIOR_HI), (
        "adjusted particles left the prior box: "
        f"{adj_joint.min(0)}, {adj_joint.max(0)}")
    assert np.allclose(adj_joint[:, 1], adj_b), "column 1 must be the b marginal"
