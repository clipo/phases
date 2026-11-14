"""Generator C4 (type-origination spec 2026-11-14, draft 1, section 8). Expected values are by hand or by a separate
calculation; the rival each test rejects is named (rule 3)."""
import numpy as np
import pytest

from mls_emergence.transmission import origination as og


def test_largest_remainder_by_hand():
    # shares 0.5, 0.3, 0.15, 0.05 at 7 draws: 3.5, 2.1, 1.05, 0.35 -> floors 3, 2, 1, 0 (sum 6); the largest remainder
    # is 0.5 (class 0) -> 4, 2, 1, 0. A class below half a draw is absent.
    assert og.largest_remainder([0.5, 0.3, 0.15, 0.05], 7).tolist() == [4, 2, 1, 0]
    assert og.largest_remainder([1, 1], 3).tolist() == [2, 1]        # tie to the lower index
    with pytest.raises(ValueError):
        og.largest_remainder([0.5, -0.1], 5)


def test_pool_step_moments_match_the_dirichlet_redraw_it_replaces():
    # One step from counts (60, 30, 10) of N = 100 at mu = 0: mean share unchanged, variance p (1 - p) / N.
    rng = np.random.default_rng(1); c = np.array([60, 30, 10])
    x = np.array([og.pool_step(c, 0.0, rng)[0] for _ in range(40000)]) / 100
    np.testing.assert_allclose(x.mean(0), [0.6, 0.3, 0.1], atol=0.002)
    np.testing.assert_allclose(x.var(0), [0.6 * 0.4 / 100, 0.3 * 0.7 / 100, 0.1 * 0.9 / 100], rtol=0.03)
    # Rival rejected: a pool that resamples only half its draws per step (variance halved).


def test_accumulated_f_over_a_window():
    # Expected heterozygosity decays by (1 - 1/N) per step, so F over W steps is 1 - (1 - 1/N)^W: N = 50, W = 40 -> 0.5543.
    rng = np.random.default_rng(2); h = []
    for _ in range(4000):
        c = np.array([25, 15, 10])
        for _s in range(40):
            c, _ = og.pool_step(c, 0.0, rng)
        h.append(1 - ((c / 50) ** 2).sum())
    h0 = 1 - (0.25 + 0.09 + 0.04)
    assert 1 - np.mean(h) / h0 == pytest.approx(1 - (1 - 1 / 50) ** 40, abs=0.015)
    # Rival rejected: N = kappa in place of kappa + 1 is not separable here; a per-step factor (1 - 2/N) gives 0.805.


def test_a_lost_type_never_returns_and_mu_zero_originates_nothing():
    rng = np.random.default_rng(3); c = np.array([0, 7, 3])
    for _ in range(200):
        c, n_new = og.pool_step(c, 0.0, rng)
        assert c[0] == 0 and n_new == 0 and c.sum() == 10


def test_stationary_richness_matches_the_ewens_expectation():
    # N = 200, theta = 2 (mu = theta / (2 N) = 0.005). Expected number of types in a sample of 50 draws:
    # sum_{i=0}^{49} 2 / (2 + i) = 7.04. Rival rejected: mu = theta / N, which behaves as theta = 4 (expected 10.9).
    rng = np.random.default_rng(4); N, theta, n = 200, 2.0, 50
    expected = sum(theta / (theta + i) for i in range(n))
    assert expected == pytest.approx(7.04, abs=0.01)
    counts = {0: N}; nxt = 1; ks = []
    for step in range(30000):
        ids = list(counts); c, n_new = og.pool_step(np.array([counts[i] for i in ids]), theta / (2 * N), rng)
        counts = {i: int(v) for i, v in zip(ids, c) if v > 0}
        for _ in range(n_new):
            counts[nxt] = 1; nxt += 1
        if step >= 3000 and step % 25 == 0:
            pool = np.repeat(list(counts), list(counts.values()))
            ks.append(len(set(rng.choice(pool, n, replace=False).tolist())))
    assert np.mean(ks) == pytest.approx(expected, rel=0.06), np.mean(ks)


def _small(theta, trace=True, seed=11, pool_seed=12, max_slots=4096):
    n = 4; w = np.full((n, n), 1 / n)
    spans = {"a": (np.array([29, 20, 12, 29]), np.array([5, 1, 13, 30])), "b": (np.array([4, 29, 15, 9]), np.array([5, 3, 1, 10]))}
    return og.origination_profiles(w, n_ind=60, mixing=0.3, innovation=0.2, steps=30, burnin=20, seed=seed, pool_seed=pool_seed,
                                   target=[0.6, 0.3, 0.1], pool_size=20, theta=theta, spans=spans, trace=trace, max_slots=max_slots), spans


def test_profiles_equal_the_brute_force_average_over_each_span():
    # The accumulators and the freeing of dead types are checked against a direct average of the traced shares per type.
    out, spans = _small(theta=6.0)
    assert out["info"]["originated"] > 20                         # new types did enter, so slots were freed and reused
    for name, (ends, lengths) in spans.items():
        ids = out["type_id"][name]; col = {int(i): j for j, i in enumerate(ids)}; brute = np.zeros((4, len(ids)))
        for site in range(4):
            for t in range(ends[site] - lengths[site] + 1, ends[site] + 1):
                tid, pa = out["trace"][t]
                for j, i in enumerate(tid):
                    if pa[site, j] > 0:
                        brute[site, col[int(i)]] += pa[site, j]
            brute[site] /= lengths[site]
        np.testing.assert_allclose(out["profiles"][name], brute, atol=1e-12)
        assert (out["origin"][name][ids < 3] == -1).all() and (out["origin"][name][ids >= 3] >= 0).all()
    # Rival rejected: an accumulator that drops a type's mass when its slot is freed (rows would not match).


def test_theta_zero_has_only_the_starting_classes_and_the_run_is_reproducible():
    a, _ = _small(theta=0.0); b, _ = _small(theta=0.0)
    assert a["info"]["originated"] == 0 and a["profiles"]["a"].shape[1] <= 3
    np.testing.assert_array_equal(a["profiles"]["a"], b["profiles"]["a"])
    c, _ = _small(theta=0.0, pool_seed=99)
    assert not np.array_equal(a["profiles"]["a"], c["profiles"]["a"])   # the pool drifts on its own stream


def test_guards():
    with pytest.raises(RuntimeError):
        _small(theta=30.0, max_slots=6)                             # slot exhaustion raises
    for kw in ({"theta": -1.0}, {"pool_size": 1}, {"pool_size": 20.0}, {"theta": 100.0}):
        args = {"n_ind": 10, "mixing": 0.1, "innovation": 0.1, "steps": 5, "burnin": 1, "seed": 0, "pool_seed": 0, "target": [0.5, 0.5],
                "pool_size": 20, "theta": 1.0, "spans": {"a": (np.array([4, 4]), np.array([1, 1]))}}
        args.update(kw)
        with pytest.raises(ValueError):
            og.origination_profiles(np.full((2, 2), 0.5), **args)
    with pytest.raises(ValueError):                                 # a span that starts before the stored steps
        og.origination_profiles(np.full((2, 2), 0.5), n_ind=10, mixing=0.1, innovation=0.1, steps=5, burnin=1, seed=0, pool_seed=0,
                                target=[0.5, 0.5], pool_size=20, theta=1.0, spans={"a": (np.array([2, 4]), np.array([4, 1]))})
