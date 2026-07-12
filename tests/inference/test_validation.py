import importlib

import numpy as np


def test_posterior_rank_bounds():
    val = importlib.import_module("analyses.39_abc_smc_validation")
    samples = np.linspace(0, 1, 100)
    w = np.ones(100)
    assert 0.0 <= val.posterior_rank(0.5, samples, w) <= 1.0
    assert val.posterior_rank(-1.0, samples, w) == 0.0
    assert val.posterior_rank(2.0, samples, w) == 1.0


def test_sbc_ranks_shape_small():
    val = importlib.import_module("analyses.39_abc_smc_validation")
    ranks = val.sbc_ranks(n_sbc=3, K=4, n_particles=50, n_rounds=2, seed=0,
                          n_workers=1)
    assert ranks.shape == (3, 4)
    assert np.all((ranks >= 0.0) & (ranks <= 1.0))


def test_sbc_ranks_order_independent():
    """Per-index seeding makes the ranks independent of worker scheduling:
    the same n_sbc must give identical ranks whether run serially or pooled."""
    val = importlib.import_module("analyses.39_abc_smc_validation")
    a = val.sbc_ranks(n_sbc=4, K=4, n_particles=50, n_rounds=2, seed=0, n_workers=1)
    b = val.sbc_ranks(n_sbc=4, K=4, n_particles=50, n_rounds=2, seed=0, n_workers=2)
    assert np.allclose(a, b)
