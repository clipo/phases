"""Pair statistics for the rim R1 readings (spec 4.2, 4.3): total variation, Spearman, partial Spearman after linear
adjustment on ranks, with the definition witnesses Codex RIM-1 and RIM-6 asked for."""
import numpy as np
import pytest

from mls_emergence.signatures import assortativity as A


def test_total_variation_is_zero_for_identical_compositions_and_reduces_to_abs_diff_for_binary():
    P = np.array([[0.5, 0.5], [0.5, 0.5], [1.0, 0.0], [0.9, 0.1]])
    tv = A.total_variation(P)
    assert tv[0, 1] == 0.0 and np.isclose(tv[2, 3], 0.1)                      # identical rows 0; (1,0) against (.9,.1) is .1
    mismatch = 1 - P @ P.T                                                     # the rejected sherd-mismatch form
    assert np.isclose(mismatch[0, 1], 0.5) and np.isclose(mismatch[2, 3], 0.1)   # it scores identical mixtures above different ones
    b = np.array([[0.2, 0.8], [0.7, 0.3], [0.5, 0.5]])
    assert np.allclose(A.total_variation(b), np.abs(b[:, [0]] - b[:, [0]].T))   # binary: |p_i - p_j|
    with pytest.raises(ValueError):
        A.total_variation(np.array([[0.6, 0.6], [0.5, 0.5]]))


def test_partial_spearman_direction_guards_and_the_quadratic_distance_witness():
    rng = np.random.default_rng(3); n = 171
    d = rng.uniform(0, 100, n)
    own = rng.normal(0, 10, n)
    assert A.partial_spearman(d + own, own + rng.normal(0, 3, n), d) > 0.6      # shared non-distance structure: positive
    assert A.partial_spearman(d - own, own + rng.normal(0, 3, n), d) < -0.6     # opposite: negative
    assert abs(A.partial_spearman(d + rng.normal(0, 10, n), d + rng.normal(0, 10, n), d)) < 0.2   # distance alone: near zero
    assert np.isnan(A.partial_spearman(np.ones(n), own, d))                     # a constant vector is undefined, not 0
    assert np.isnan(A.partial_spearman(d, own, d))                              # zero residual after removing distance: undefined
    # The witness of the statistic's limit: both readouts are a quadratic function of distance alone, yet the linear rank
    # adjustment leaves a residual association of 1. The statistic is association after LINEAR adjustment on ranks, not
    # "at a given distance".
    q = (np.linspace(0, 1, n) - 0.5) ** 2; dd = np.linspace(0, 1, n)
    assert np.isclose(A.partial_spearman(q, q, dd), 1.0)
    assert np.isclose(A.spearman(np.arange(10), np.arange(10) ** 2), 1.0) and np.isnan(A.spearman(np.ones(5), np.arange(5)))
