"""ACF estimator, doubling rule and stationarity (spec 3.3; Finding 8).

Expected values come from the AR(1) autocorrelation phi^h, not from the
functions under test.
"""

import numpy as np

from mls_emergence.probe import timescales as ts


def ar1(phi, S, rng, n_series=1):
    x = np.zeros((S, n_series)); e = rng.standard_normal((S, n_series))
    for t in range(1, S):
        x[t] = phi * x[t - 1] + e[t]
    return x


def test_acf_matches_ar1_theory_and_first_crossing_is_near_minus_one_over_log_phi():
    rng = np.random.default_rng(3); phi = 0.95; S = 200_000
    x = ar1(phi, S, rng)[:, 0]
    r = ts.acf(x, 200)
    assert abs(r[0] - 1) < 1e-12
    for h in (1, 10, 40):
        assert abs(r[h] - phi ** h) < 0.02
    # Theoretical crossing of 1/e: h = -1 / log(phi) = 19.5; estimator gives 19 or 20.
    assert ts.first_crossing(r) in (19, 20)
    # A faster process crosses sooner (discriminates against a constant answer).
    y = ar1(0.5, S, rng)[:, 0]
    assert ts.first_crossing(ts.acf(y, 50)) in (1, 2)


def test_local_and_regional_curves_exclude_flat_series_and_average_the_rest():
    rng = np.random.default_rng(4); S, N, K = 20_000, 4, 3
    freq = np.zeros((S, N, K))
    for i in range(N):
        freq[:, i, :2] = ar1(0.9, S, rng, 2)
    # Class 3 is flat everywhere: excluded from the local curve (N series) and from the regional one (1 series).
    curve, n_in, n_ex = ts.mean_local_curve(freq, max_lag=100)
    assert n_in == N * 2 and n_ex == N
    assert abs(curve[10] - 0.9 ** 10) < 0.05
    rc, _, n_exr = ts.regional_curve(freq, max_lag=100)
    assert n_exr == 1 and rc is not None
    # All flat: no curve.
    c2, n2, e2 = ts.mean_local_curve(np.zeros((100, 2, 2)))
    assert c2 is None and n2 == 0 and e2 == 4


def test_biased_acf_of_a_smooth_series_crosses_by_a_fixed_fraction_of_the_window():
    # The reason for the S/10 acceptance rule: a linear trend crosses 1/e at
    # about 0.22 S for every S, so "no crossing by S/4" could never trigger.
    for S in (100, 400, 1600):
        c = ts.first_crossing(ts.acf(np.arange(S, dtype=float), S // 4))
        assert c is not None and 0.20 * S <= c <= 0.24 * S


def test_measure_time_scales_doubles_until_an_accepted_crossing_and_reports_status():
    # Stand-in for a process slow relative to the window: a linear trend,
    # whose crossing (about 0.22 S) is never within S/10, so the window must
    # double; once the window reaches 400 the stand-in returns an AR(1) with
    # phi = 0.95, whose theory crossing -1/log(0.95) = 19.5 lies within
    # 400/10. The function must double twice (100 -> 200 -> 400) and report
    # the AR(1) crossing (the biased estimator at S = 400 sits a little low).
    def run(seed, S):
        x = np.arange(S, dtype=float) if S < 400 else ar1(0.95, S, np.random.default_rng(seed))[:, 0]
        return x[:, None, None]
    out = ts.measure_time_scales(run, seeds=[1, 2], S0=100, max_doublings=3)
    assert out["status"] == "ok" and out["doublings"] == 2 and out["S"] == 400
    assert 12 <= out["L_d"] <= 22 and 12 <= out["L_R"] <= 22
    # With doublings capped below what the crossing needs, the status says so.
    out2 = ts.measure_time_scales(run, seeds=[1], S0=100, max_doublings=1)
    assert out2["status"] == "no_crossing" and out2["S"] == 200
    # All-flat input reports all_excluded.
    out3 = ts.measure_time_scales(lambda seed, S: np.zeros((S, 2, 2)), seeds=[1], S0=50)
    assert out3["status"] == "all_excluded"


def test_stationarity_rule_passes_a_stationary_run_and_fails_a_drifting_one():
    rng = np.random.default_rng(6); S, N, K = 4000, 6, 3
    stat = [rng.standard_normal((S, N, K)) * 0.1 for _ in range(5)]
    assert ts.stationarity(stat)["stationary"]
    drift = [rng.standard_normal((S, N, K)) * (0.1 + 0.5 * np.linspace(0, 1, S))[:, None, None] for _ in range(5)]
    assert not ts.stationarity(drift)["stationary"]
    d = ts.between_site_dispersion(stat[0])
    assert d.shape == (S,) and abs(d.mean() - 0.01) < 0.002   # variance of N(0, 0.1^2) across sites is 0.01
