"""Regional drift of the innovation target (regional-drift spec, sections 3.1 and 3.3; tests 7.1 to 7.5).

The target is redrawn from Dirichlet(kappa * target) after every post-burn-in
learner draw. With target_kappa=None the model is the paper's model, exactly.
"""
import numpy as np
import pytest

from mls_emergence.transmission import spatial as sp

W_PRIMARY = 2336
GRID = (0.01, 0.03, 0.1, 0.3)
KAPPA_PRIMARY = (232430.0, 76692.0, 22171.0, 6549.0)   # Codex's inversion at W = 2336 (draft-2 technical notes, section 1)
POOLED = np.array([7702, 4504, 1224, 150, 215, 138, 49, 67, 34, 18], float) / 14101


def _weights(n=3, seed=0):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0, 30, (n, 2))
    d = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    return sp.copying_weights(d, 6.0)


def test_no_kappa_is_the_current_model():
    """target_kappa=None consumes no randomness: identical output to a call without the argument."""
    w = _weights()
    a = sp.drift_record(w, k=4, n_ind=200, mixing=.1, innovation=.01, steps=60, n_records=60, burnin=20, seed=11,
                        target=np.array([.5, .3, .15, .05]))
    b = sp.drift_record(w, k=4, n_ind=200, mixing=.1, innovation=.01, steps=60, n_records=60, burnin=20, seed=11,
                        target=np.array([.5, .3, .15, .05]), target_kappa=None)
    assert np.array_equal(a, b)


def test_step_target_moments_and_zeros():
    rng = np.random.default_rng(1)
    p = np.array([.5, .3, .15, .05]); kappa = 1000.0; n = 200_000
    draws = np.array([sp.step_target(p, kappa, rng) for _ in range(n)])
    cov_true = (np.diag(p) - np.outer(p, p)) / (kappa + 1)
    se_mean = np.sqrt(np.diag(cov_true) / n)
    assert np.all(np.abs(draws.mean(0) - p) < 4 * se_mean)          # a biased step would mean the recursion is not mean-preserving
    cov = np.cov(draws.T)
    for i in range(4):
        for j in range(4):
            tol = max(0.05 * abs(cov_true[i, j]), 4 * abs(cov_true[i, j]) * np.sqrt(2.0 / n))
            assert abs(cov[i, j] - cov_true[i, j]) < tol, (i, j, cov[i, j], cov_true[i, j])
    z = sp.step_target(np.array([0.0, 0.5, 0.5]), kappa, rng)
    assert z[0] == 0.0 and abs(z.sum() - 1) < 1e-12                  # a zero component stays exactly zero (numerical absorbing face)
    one = sp.step_target(np.array([0.0, 1.0, 0.0]), kappa, rng)
    assert np.array_equal(one, np.array([0.0, 1.0, 0.0]))            # the single-survivor profile is returned unchanged
    near = sp.step_target(p, 1e9, rng)
    assert np.max(np.abs(near - p)) < 1e-4                           # probe point: a huge kappa is the identity to 1e-4
    for bad in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError):
            sp.step_target(p, bad, rng)
    with pytest.raises(ValueError):
        sp.step_target(np.tile(p, (2, 1)), kappa, rng)
    with pytest.raises(ValueError):
        sp.step_target(np.array([1.0]), kappa, rng)


def test_kappa_for_window_f_closed_form():
    for F, k_exp in zip(GRID, KAPPA_PRIMARY):
        k = sp.kappa_for_window_f(F, W_PRIMARY)
        assert np.isclose(k, k_exp, rtol=5e-4), (F, k)
        assert np.isclose(1 - (k / (k + 1)) ** W_PRIMARY, F, rtol=1e-9)
    assert sp.kappa_for_window_f(0.0, W_PRIMARY) is None
    for bad in (-0.1, 1.0, 1.5, np.nan, np.inf):
        with pytest.raises(ValueError):
            sp.kappa_for_window_f(bad, W_PRIMARY)
    for badW in (0, -3, 2.5, True):
        with pytest.raises(ValueError):
            sp.kappa_for_window_f(0.1, badW)


def test_regional_f_accumulates_to_window_f():
    """E[p(1-p)] after W updates equals (kappa/(kappa+1))^W times its start value: ratio 0.9 at F_W = 0.1 (Codex witness 0.901, 0.899, 0.894)."""
    rng = np.random.default_rng(2); W = 300; F = 0.1; kappa = sp.kappa_for_window_f(F, W); n = 5000
    ends = np.empty((n, len(POOLED)))
    for r in range(n):
        p = POOLED.copy()
        for _ in range(W):
            p = sp.step_target(p, kappa, rng)
        ends[r] = p
    ratio = (ends * (1 - ends)).mean(0) / (POOLED * (1 - POOLED))
    assert np.all(np.abs(ratio[:3] - (1 - F)) < 0.02 * (1 - F)), ratio[:3]
    assert sp.kappa_for_window_f(0.0, W) is None                     # probe point: F_W = 0 leaves p(1-p) exactly unchanged (no update)


def test_update_runs_on_every_step_not_only_saved_ones():
    w = _weights(); tgt = np.array([.5, .3, .15, .05]); kw = {"k": 4, "n_ind": 200, "mixing": .1, "innovation": .01, "burnin": 10, "seed": 5,
                                                            "target": tgt, "target_kappa": 500.0, "return_targets": True}
    rec_d, tg_d = sp.drift_record(w, steps=400, n_records=400, **kw)
    rec_s, tg_s = sp.drift_record(w, steps=400, n_records=100, **kw)
    saved = sorted(set(np.linspace(0, 399, 100).astype(int)))
    assert tg_d.shape == (400, 4) and tg_s.shape == (100, 4)
    assert np.array_equal(tg_d[saved], tg_s) and np.array_equal(rec_d[saved], rec_s)
    assert np.array_equal(tg_d[0], tgt)                              # the first recorded state used the unwalked target
    assert not np.array_equal(tg_d[1], tgt)                          # and the target has moved by the second
    with pytest.raises(ValueError):
        sp.drift_record(w, steps=20, n_records=20, k=4, n_ind=50, mixing=.1, innovation=.01, burnin=2, seed=1,
                        target=np.tile(tgt, (3, 1)), target_kappa=500.0)
    with pytest.raises(ValueError):
        sp.drift_record(w, steps=20, n_records=20, k=1, n_ind=50, mixing=.1, innovation=.01, burnin=2, seed=1,
                        target=np.array([1.0]), target_kappa=500.0)


def test_domain_guards_reject_strings_denormals_and_bools():
    for bad in ("0.1", 5e-324):
        with pytest.raises(ValueError):
            sp.kappa_for_window_f(bad, 2336)
    w = _weights()
    with pytest.raises(ValueError):
        sp.drift_record(w, k=4, n_ind=50, mixing=.1, innovation=.01, steps=20, n_records=20, burnin=2, seed=1,
                        target=np.array([.5, .3, .15, .05]), target_kappa=True)
