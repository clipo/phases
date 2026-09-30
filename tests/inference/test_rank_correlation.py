"""Tests for the Bayesian rank correlation that replaced Spearman-plus-p-value.

This module had no tests at all. It is not incidental code: it is what several
manuscript claims now rest on after the rule-18 conversion, including the
"no resolved divergence trend" statement in the main text.

Rule 3 governs what is written here. Every expected value is derived by hand or
from a property that holds independently of the implementation, never by calling
the code under test. Where a property is asserted, there is a probe point at
which the nearest rival account (a Pearson correlation on the raw values)
predicts something different, so a test can actually fail.
"""
from __future__ import annotations

import numpy as np
import pytest

from mls_emergence.inference.rank_correlation import (
    bayesian_rank_correlation,
    format_result,
    normal_scores,
)


# --------------------------------------------------------------------------- #
# normal_scores
# --------------------------------------------------------------------------- #

def test_normal_scores_are_antisymmetric_and_centered():
    """Phi^-1((rank - 0.5)/n) is odd about the middle rank, so the scores sum to
    zero and mirror each other. Hand-checked, not read back from the function."""
    z = normal_scores([10.0, 20.0, 30.0, 40.0, 50.0])
    assert z.sum() == pytest.approx(0.0, abs=1e-12)
    assert z[0] == pytest.approx(-z[-1], abs=1e-12)
    assert z[1] == pytest.approx(-z[-2], abs=1e-12)
    assert z[2] == pytest.approx(0.0, abs=1e-12)
    assert np.all(np.diff(z) > 0)  # strictly increasing in the input


def test_normal_scores_depend_only_on_order_not_on_spacing():
    """The discriminating property. A rank-based score is invariant under ANY
    strictly monotone transformation; a correlation on the raw values is not.
    The second assertion is the probe: it fails for the rival account."""
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    wild = np.exp(10.0 * x)  # same order, wildly different spacing
    assert np.allclose(normal_scores(x), normal_scores(wild))

    # The rival: Pearson on the raw values does move under the same monotone
    # transform, so the invariance asserted above is a real property and not a
    # tautology. Kept to a range that does not overflow.
    y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert not np.isclose(np.corrcoef(x, y)[0, 1],
                          np.corrcoef(x, np.exp(2.0 * y))[0, 1])


def test_normal_scores_average_ties():
    """Tied inputs must receive identical scores, and the untied neighbours must
    straddle them. Averaging ranks is the documented tie rule."""
    z = normal_scores([1.0, 5.0, 5.0, 9.0])
    assert z[1] == pytest.approx(z[2], abs=1e-12)
    assert z[0] < z[1] < z[3]
    assert z.sum() == pytest.approx(0.0, abs=1e-12)


def test_normal_scores_ignore_non_finite_entries():
    """NaN propagates to its own slot but must not shift the finite scores,
    which are computed over the finite subset only."""
    z_all = normal_scores([1.0, 2.0, 3.0])
    z_gap = normal_scores([1.0, 2.0, np.nan, 3.0])
    assert np.isnan(z_gap[2])
    assert np.allclose(z_gap[[0, 1, 3]], z_all)


@pytest.mark.parametrize(
    "bad, match",
    [
        (np.zeros((2, 3)), "1-D"),
        ([1.0, 2.0], "at least 3"),
        ([1.0, np.nan, np.nan, np.nan], "at least 3"),
    ],
)
def test_normal_scores_guards_its_domain(bad, match):
    """Rule 5: out-of-domain input raises with a reason rather than returning
    something plausible-looking."""
    with pytest.raises(ValueError, match=match):
        normal_scores(bad)


# --------------------------------------------------------------------------- #
# bayesian_rank_correlation
# --------------------------------------------------------------------------- #

def test_recovers_a_known_correlation_from_bivariate_normal_draws():
    """Recovery against a truth the prior does not favour. r ~ Uniform(-1, 1)
    puts only 10 percent of its mass above 0.8, so recovering 0.8 is a real
    check and not the prior speaking.

    The assertion is a band on the posterior median, not coverage of the truth
    by one 95 percent interval. A single interval misses its own truth 5 percent
    of the time by construction, so a coverage assertion at one seed is a test
    that fails one run in twenty for the right reason, which is worse than
    useless. The band is wide enough to survive that and narrow enough to fail
    if the estimator were wrong: it excludes zero, excludes one, and excludes
    the Spearman-transform value 0.786 being mistaken for something far off.

    The estimand is the generating Pearson r itself, checked separately rather
    than assumed: the normal-scores correlation converges to it (0.7999 at
    n = 20,000 against a true 0.80), with mild small-sample attenuation
    (0.7947 at n = 250, 0.7827 at n = 50).
    """
    rng = np.random.default_rng(11)
    true_r = 0.8
    cov = np.array([[1.0, true_r], [true_r, 1.0]])
    xy = rng.multivariate_normal([0.0, 0.0], cov, size=250)

    res = bayesian_rank_correlation(xy[:, 0], xy[:, 1], draws=400, tune=400,
                                    chains=2, random_seed=3)
    assert abs(res["r_median"] - true_r) < 0.1, res["hdi95"]
    lo, hi = res["hdi95"]
    assert lo > 0.6, "posterior should resolve a strong positive correlation"
    assert hi < 0.95, "posterior should not run up against the boundary"
    assert res["p_positive"] > 0.99


def test_monotone_but_strongly_nonlinear_data_reads_as_near_perfect():
    """The discriminating case for a RANK method. y = x**3 over a range that
    includes negatives is perfectly monotone, so a rank correlation should sit
    near +1, while Pearson on the raw values is visibly below 1. If this module
    were secretly correlating raw values, the first assertion would fail."""
    x = np.linspace(-3.0, 3.0, 40)
    y = x ** 3
    pearson_raw = float(np.corrcoef(x, y)[0, 1])
    assert pearson_raw < 0.95, "probe is only meaningful if Pearson is attenuated"

    res = bayesian_rank_correlation(x, y, draws=400, tune=400, chains=2,
                                    random_seed=4)
    assert res["r_median"] > 0.97
    assert res["hdi95"][0] > 0.9


def test_reversed_order_gives_a_resolved_negative_correlation():
    """Sign check with a one-sided assertion. The other side would mean the
    posterior had placed a decreasing relationship on the increasing side of
    zero, which is a sign error rather than a precision problem."""
    x = np.arange(30, dtype=float)
    res = bayesian_rank_correlation(x, -x, draws=400, tune=400, chains=2,
                                    random_seed=5)
    assert res["r_median"] < -0.9
    assert res["p_positive"] < 0.01


def test_small_n_widens_the_interval_rather_than_reporting_absence():
    """The module's stated reason for existing: with six bins the posterior is
    wide, where the frequentist version said 'not significant' and read as
    evidence of absence. Same underlying association, fewer points."""
    rng = np.random.default_rng(7)
    cov = np.array([[1.0, 0.6], [0.6, 1.0]])
    big = rng.multivariate_normal([0.0, 0.0], cov, size=200)

    wide = bayesian_rank_correlation(big[:6, 0], big[:6, 1], draws=400,
                                     tune=400, chains=2, random_seed=6)
    narrow = bayesian_rank_correlation(big[:, 0], big[:, 1], draws=400,
                                       tune=400, chains=2, random_seed=6)
    w_width = wide["hdi95"][1] - wide["hdi95"][0]
    n_width = narrow["hdi95"][1] - narrow["hdi95"][0]
    assert w_width > 3 * n_width
    assert wide["n"] == 6 and narrow["n"] == 200


def test_every_rule_16_diagnostic_is_present_and_finite():
    """Rule 16: diagnostics are attached by construction, so a caller cannot
    report this correlation without them. Asserting the keys exist is not
    enough; a NaN R-hat would satisfy that and mean nothing."""
    x = np.arange(25, dtype=float)
    rng = np.random.default_rng(2)
    res = bayesian_rank_correlation(x, x + rng.normal(0, 5, x.size), draws=400,
                                    tune=400, chains=2, random_seed=8)
    for key in ("rhat", "ess_bulk", "ess_tail", "n_div", "n_total",
                "pct_div", "n_treedepth", "ebfmi"):
        assert key in res, f"missing diagnostic {key}"
        assert np.isfinite(res[key]), f"non-finite diagnostic {key}"
    assert res["rhat"] < 1.05
    assert res["ess_bulk"] > 100
    assert res["n_total"] == 800  # 2 chains x 400 draws, not the tuning draws
    assert format_result(res, "trend ").startswith("trend r = ")
    assert "R-hat" in format_result(res)


def test_guards_too_few_complete_pairs():
    """Two complete pairs cannot support a correlation and must raise rather
    than return a confident-looking number off a degenerate fit."""
    with pytest.raises(ValueError, match="at least 3"):
        bayesian_rank_correlation([1.0, 2.0, 3.0], [1.0, np.nan, np.nan])
