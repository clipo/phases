"""Rim R1 statistics (spec sections 4, 6, 7): posterior draws, geography and association statistics, prior predictives,
the recovery program. Tests written before the module."""
import numpy as np
import pytest

from mls_emergence.rims import stats as st


def _sites(seed=0, n=19):
    rng = np.random.default_rng(seed)
    lat = rng.uniform(34.7, 35.6, n); lon = rng.uniform(-90.8, -90.1, n)
    xy = np.column_stack([(lon + 90.5) * 101.0, (lat - 35.0) * 111.32]); D = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))
    N = rng.integers(27, 278, n)
    return lat, lon, D, N


def test_frequency_draws_have_the_conjugate_mean_and_prior_weight():
    rng = np.random.default_rng(1)
    k = np.array([3, 30, 0]); N = np.array([27, 100, 50])
    d = st.frequency_draws(k, N, rng, draws=200000, alpha=0.5)
    assert d.shape == (200000, 3)
    assert np.allclose(d.mean(0), (k + 0.5) / (N + 1), atol=3e-3)              # one pseudo-observation in total: weight 1/(N+1)
    assert np.allclose(st.prior_weight(N), 1 / (N + 1)) and np.isclose(st.prior_weight(27), 1 / 28)
    with pytest.raises(ValueError):
        st.frequency_draws(np.array([5]), np.array([4]), rng, draws=10)        # k > N
    c = st.composition_draws(np.array([[10, 5, 0], [2, 2, 2]]), rng, draws=20000, alpha=0.5)
    assert c.shape == (20000, 2, 3) and np.allclose(c.sum(-1), 1) and np.allclose(c.mean(0)[0], (np.array([10, 5, 0]) + 0.5) / 16.5, atol=5e-3)


def test_geography_statistics_recover_constructed_directions():
    lat, lon, D, _N = _sites()
    p_north = 0.05 + 0.6 * (lat - lat.min()) / np.ptp(lat)                     # frequency rises to the north
    g = st.geography_stats(p_north, D, lat, lon)
    assert g["rho_latitude"] > 0.95 and g["rho_pair_distance"] > 0.3            # pairs farther apart differ more
    g0 = st.geography_stats(np.full(19, 0.3), D, lat, lon)
    assert all(np.isnan(g0[k]) for k in ("rho_latitude", "rho_longitude", "rho_pair_distance"))   # constant: undefined, not 0
    south = st.geography_stats(1 - p_north, D, lat, lon)
    assert south["rho_latitude"] < -0.95                                        # the sign follows the direction: falls to the north


def test_pair_association_uses_total_variation_and_linear_distance_adjustment():
    _lat, _lon, D, _N = _sites()
    rng = np.random.default_rng(5)
    shared = rng.normal(size=19)                                                 # a non-geographic pattern both traits carry
    p = 0.5 + 0.3 * np.tanh(shared)
    comp = np.column_stack([0.4 + 0.3 * np.tanh(shared), 0.3 - 0.15 * np.tanh(shared), 0.3 - 0.15 * np.tanh(shared)])
    a = st.pair_association(p, comp, D)
    assert a["partial_spearman"] > 0.8 and a["spearman"] > 0.8
    comp_un = np.tile([0.5, 0.3, 0.2], (19, 1))
    assert np.isnan(st.pair_association(p, comp_un, D)["partial_spearman"])      # identical compositions: zero differences, undefined


def test_posterior_summary_reports_median_interval_and_monte_carlo_error():
    rng = np.random.default_rng(2)
    lat, lon, D, N = _sites(); k = rng.binomial(N, 0.2)
    s = st.posterior_summary(lambda p: st.geography_stats(p, D, lat, lon)["rho_latitude"], k, N, draws=2000, seeds=(10, 11))
    assert set(s) >= {"observed", "median", "q05", "q95", "mc_error_median"} and s["q05"] <= s["median"] <= s["q95"]
    assert s["mc_error_median"] < 0.1


def test_prior_predictive_counts_are_beta_binomial_at_the_actual_totals():
    N = np.array([27, 100, 277])
    pp = st.prior_predictive_counts(N, alpha=0.5)
    assert pp.shape == (3, 3) and np.allclose(pp[:, 1] / N, 0.5, atol=0.05)    # median share near 1/2 under Beta(1/2, 1/2)
    assert (pp[:, 0] < N * 0.05).all() and (pp[:, 2] > N * 0.95).all()         # 5th and 95th percentiles near the ends


def test_recovery_program_covers_truth_and_recovers_signs():
    """Spec section 7: 20 replicates at three true configurations on actual-sized totals; 90 percent intervals cover the
    true rho_latitude in a share consistent with 0.9 (binomial band for 60 trials: at least 0.80); signs recovered where
    the truth is non-zero."""
    lat, lon, D, N = _sites()
    out = st.recovery_program(D, lat, lon, N, replicates=20, draws=1000, seed=7)
    assert set(out) == {"positive", "zero", "opposite"}
    cover = np.mean([out[c]["coverage"] for c in out])
    assert cover >= 0.80, out
    assert out["positive"]["sign_recovered"] >= 0.9 and out["opposite"]["sign_recovered"] >= 0.9
