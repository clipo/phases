import numpy as np

from mls_emergence.inference import (
    abc_smc,
    weighted_mean,
    weighted_quantile,
    resample,
)


# ---- analytic conjugate problem: Normal mean, known variance -------------
def _normal_mean_problem(seed=0, n_obs=30, sigma=1.0, m0=0.0, s0=5.0):
    """Return (callbacks, s_obs, analytic_mean, analytic_sd) for inferring the
    mean of a Normal(mu, sigma) with a Normal(m0, s0) prior from the sample mean
    of n_obs draws. Posterior is analytic (conjugate)."""
    rng = np.random.default_rng(seed)
    mu_true = 2.0
    y = rng.normal(mu_true, sigma, n_obs)
    y_bar = float(y.mean())

    def prior_sampler(r):
        return np.array([r.normal(m0, s0)])

    def prior_logpdf(theta):
        z = (theta[0] - m0) / s0
        return -0.5 * z * z  # up to constant; proper for weighting

    def simulator(theta, r):
        return r.normal(theta[0], sigma, n_obs)

    def summary(sim):
        return np.array([sim.mean()])

    def distance(s, s_obs):
        return float(abs(s[0] - s_obs[0]))

    s_obs = np.array([y_bar])
    post_var = 1.0 / (1.0 / s0**2 + n_obs / sigma**2)
    post_mean = post_var * (m0 / s0**2 + n_obs * y_bar / sigma**2)
    return (prior_sampler, prior_logpdf, simulator, summary, distance,
            s_obs, post_mean, post_var**0.5)


def test_recovers_conjugate_normal_posterior():
    (ps, pl, sim, summ, dist, s_obs, post_mean, post_sd) = _normal_mean_problem()
    res = abc_smc(ps, pl, sim, summ, dist, s_obs,
                  n_particles=1500, n_rounds=8, quantile=0.4, seed=1)
    est_mean = weighted_mean(res.thetas[:, 0], res.weights)
    lo = weighted_quantile(res.thetas[:, 0], res.weights, 0.025)
    hi = weighted_quantile(res.thetas[:, 0], res.weights, 0.975)
    est_sd = (hi - lo) / (2 * 1.959963985)
    assert abs(est_mean - post_mean) < 0.15
    assert abs(est_sd - post_sd) < 0.10


def test_tolerance_non_increasing():
    (ps, pl, sim, summ, dist, s_obs, *_) = _normal_mean_problem()
    res = abc_smc(ps, pl, sim, summ, dist, s_obs,
                  n_particles=400, n_rounds=6, seed=2)
    eps = np.array(res.eps_schedule)
    assert np.all(np.diff(eps) <= 1e-9)


def test_particles_in_support_and_weights_normalized():
    (ps, pl, sim, summ, dist, s_obs, *_) = _normal_mean_problem()
    res = abc_smc(ps, pl, sim, summ, dist, s_obs,
                  n_particles=400, n_rounds=5, seed=3)
    assert np.isclose(res.weights.sum(), 1.0)
    assert np.all(np.isfinite([pl(t) for t in res.thetas]))


def test_seed_determinism():
    (ps, pl, sim, summ, dist, s_obs, *_) = _normal_mean_problem()
    r1 = abc_smc(ps, pl, sim, summ, dist, s_obs, n_particles=300, n_rounds=4, seed=7)
    r2 = abc_smc(ps, pl, sim, summ, dist, s_obs, n_particles=300, n_rounds=4, seed=7)
    assert np.array_equal(r1.thetas, r2.thetas)
    assert np.array_equal(r1.weights, r2.weights)


def test_resample_and_weighted_quantile():
    vals = np.array([0.0, 1.0, 2.0, 3.0])
    w = np.array([0.1, 0.2, 0.3, 0.4])
    assert abs(weighted_quantile(vals, w, 0.5) - 2.0) < 0.6
    rng = np.random.default_rng(0)
    draws = resample(vals, w, 10000, rng)
    # heavier weight on larger values -> mean above the unweighted 1.5
    assert draws.mean() > 1.5


def test_bounded_prior_keeps_particles_in_support():
    """A prior bounded to [0, 1] with truth near the upper edge forces some
    perturbed proposals out of support, exercising the prior_logpdf rejection
    branch. Every accepted particle must stay in support and recovery must hold."""
    rng0 = np.random.default_rng(0)
    lo, hi = 0.0, 1.0
    p_true = 0.9
    y = (rng0.uniform(0, 1, 50) < p_true).astype(float)
    y_bar = float(y.mean())

    def prior_sampler(r):
        return np.array([r.uniform(lo, hi)])

    def prior_logpdf(theta):
        return 0.0 if lo <= theta[0] <= hi else -np.inf

    def simulator(theta, r):
        return (r.uniform(0, 1, 50) < theta[0]).astype(float)

    def summary(sim):
        return np.array([sim.mean()])

    def distance(s, s_obs):
        return float(abs(s[0] - s_obs[0]))

    res = abc_smc(prior_sampler, prior_logpdf, simulator, summary, distance,
                  np.array([y_bar]), n_particles=600, n_rounds=5, seed=1)
    # every accepted particle lies inside the bounded support
    assert np.all((res.thetas[:, 0] >= lo) & (res.thetas[:, 0] <= hi))
    # posterior mean near the observed rate
    est = weighted_mean(res.thetas[:, 0], res.weights)
    assert abs(est - y_bar) < 0.15
    # n_sims now counts simulator calls only (finite, at least one per round)
    assert all(n >= 1 for n in res.n_sims)


def test_weights_raise_rather_than_silently_going_uniform():
    """F7. The linear-scale weight recursion underflowed to zero when the prior
    log-density was strongly negative, and the guard replaced an all-zero vector
    with UNIFORM weights: a plausible-looking result that has discarded the
    importance-sampling correction entirely.

    Probe: an ordinary prior must still produce non-uniform weights, or the test
    would pass on a sampler that had stopped weighting at all."""
    import numpy as np
    import pytest
    # NOTE: mls_emergence.inference re-exports the abc_smc FUNCTION, so the
    # module has to be imported by its full path.
    from mls_emergence.inference.abc_smc import abc_smc as run_smc

    rng = np.random.default_rng(0)

    def prior_sampler(r):
        return r.uniform(-1.0, 1.0, size=2)

    def simulator(theta, r):
        return theta + r.normal(0, 0.1, size=2)

    summary = lambda x: np.asarray(x, float)
    distance = lambda a, b: float(np.sqrt(np.sum((a - b) ** 2)))
    s_obs = np.array([0.0, 0.0])

    # ordinary run: weights must vary, otherwise the probe below proves nothing
    res = run_smc(prior_sampler, lambda th: 0.0, simulator, summary,
                  distance, s_obs, n_particles=40, n_rounds=3, seed=1)
    assert np.isclose(res.weights.sum(), 1.0)
    assert res.weights.std() > 0, "weights are uniform; the probe cannot discriminate"

    # a prior returning -inf everywhere after round 0 makes every weight
    # non-finite; that must raise, not quietly become uniform
    calls = {"n": 0}

    def collapsing_logpdf(theta):
        calls["n"] += 1
        return 0.0 if calls["n"] <= 40 else -np.inf

    with pytest.raises((FloatingPointError, RuntimeError)):
        run_smc(prior_sampler, collapsing_logpdf, simulator, summary,
                distance, s_obs, n_particles=40, n_rounds=3, seed=2,
                max_sims_per_round=2000)
