"""Analytic pins for the A1 readouts (D15, D14, D18; R2 items 4 and 5)."""
import numpy as np
import pytest

from mls_emergence.probe import readouts as ro

TRUTH = np.repeat([0, 1, 2], [12, 11, 5])


def test_delta_contrast_cases():
    pi = np.array([0.5, 0.3, 0.2])
    assert ro.delta_contrast(np.ones((3, 9)) * 0.7, pi) == 0.0                      # all shifts equal
    d = np.random.default_rng(0).normal(size=(3, 9))
    assert ro.delta_contrast(d, np.array([1.0, 0.0, 0.0])) == pytest.approx(0.0)      # one occupied component
    g = 1.3; d2 = np.zeros((2, 9)); d2[0, 3] = g / 2; d2[1, 3] = -g / 2
    assert ro.delta_contrast(d2, np.array([0.5, 0.5])) == pytest.approx(g / 2)        # two equal groups, gap g in one class
    # Invariant to a constant added to all shifts in a class.
    assert ro.delta_contrast(d + 3.0, pi) == pytest.approx(ro.delta_contrast(d, pi))
    # Prior identity, Monte Carlo: E[Delta^2 | pi, sigma_g] = 9 sigma_g^2 (1 - sum pi^2) for the sum-to-zero projection.
    rng = np.random.default_rng(1); sg = 1.7; n = 20000
    acc = 0.0
    for _ in range(n):
        u = rng.normal(0, sg, size=(3, 9)); u -= u.mean(axis=0, keepdims=True)
        acc += ro.delta_contrast(u, pi) ** 2
    expect = 9 * sg ** 2 * (1 - np.sum(pi ** 2))
    assert abs(acc / n - expect) / expect < 0.03
    with pytest.raises(ValueError):
        ro.delta_contrast(d, np.array([0.5, 0.6, 0.2]))


def test_comembership_and_coherence_cases():
    # Hard labels on compact groups along a path graph: positive coherence.
    z = np.array([0, 0, 0, 1, 1, 1])
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
    q = ro.hard_comembership(z)
    assert ro.coherence(q, edges) == pytest.approx(4 / 5 - 6 / 15)
    # Alternating labels: neighbors never share a group, coherence negative.
    assert ro.coherence(ro.hard_comembership(np.array([0, 1, 0, 1, 0, 1])), edges) == pytest.approx(0 - 6 / 15)
    # One group holds every site: exactly zero.
    assert ro.coherence(ro.hard_comembership(np.zeros(6, int)), edges) == 0.0
    # Soft: responsibilities p_i = (1/2, 1/2) everywhere give q = 1/2 off the diagonal and coherence 0.
    resp = np.full((6, 2), 0.5)
    q = ro.comembership(resp)
    assert np.all(np.diag(q) == 1) and np.allclose(q[np.triu_indices(6, 1)], 0.5)
    assert ro.coherence(q, edges) == pytest.approx(0.0)
    # The Codex enumeration (R2 technical notes, section 3): rows (.8,.2), (.7,.3), (.3,.7), (.1,.9) on a path give soft C = 0.12.
    resp = np.array([[.8, .2], [.7, .3], [.3, .7], [.1, .9]])
    assert ro.coherence(ro.comembership(resp), np.array([[0, 1], [1, 2], [2, 3]])) == pytest.approx(0.12, abs=1e-12)


def test_separation_margins_on_the_three_named_partitions():
    perfect = ro.hard_comembership(TRUTH)
    assert np.allclose(ro.separation_margins(perfect, TRUTH), [1, 1, 1])
    all_one = ro.hard_comembership(np.zeros(28, int))
    assert np.allclose(ro.separation_margins(all_one, TRUTH), [0, 0, 0])
    merged = ro.hard_comembership(np.where(TRUTH == 2, 1, TRUTH))
    assert np.allclose(ro.separation_margins(merged, TRUTH), [1, 0, 0])
    B = ro.block_matrix(perfect, TRUTH)
    assert np.allclose(B, np.eye(3))
