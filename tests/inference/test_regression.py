import numpy as np

from mls_emergence.inference import abc_smc, regression_adjust, weighted_mean


def _binomial_rate_problem(seed=0, n=40, a=2.0, b=2.0):
    """Infer Binomial rate p with a Beta(a,b) prior; k is sufficient so the
    analytic posterior is Beta(a+k, b+n-k)."""
    rng = np.random.default_rng(seed)
    p_true = 0.7
    k_obs = int(rng.binomial(n, p_true))

    def prior_sampler(r):
        return np.array([r.beta(a, b)])

    def prior_logpdf(theta):
        p = theta[0]
        if p <= 0.0 or p >= 1.0:
            return -np.inf
        return (a - 1) * np.log(p) + (b - 1) * np.log(1 - p)

    def simulator(theta, r):
        return np.array([r.binomial(n, theta[0])])

    def summary(sim):
        return np.array([sim[0] / n])

    def distance(s, s_obs):
        return float(abs(s[0] - s_obs[0]))

    s_obs = np.array([k_obs / n])
    post_mean = (a + k_obs) / (a + b + n)
    return (prior_sampler, prior_logpdf, simulator, summary, distance,
            s_obs, post_mean)


def test_regression_adjustment_reduces_bias():
    (ps, pl, sim, summ, dist, s_obs, post_mean) = _binomial_rate_problem()
    # Deliberately loose tolerance schedule so the unadjusted posterior retains
    # ABC bias that the local-linear adjustment should remove.
    res = abc_smc(ps, pl, sim, summ, dist, s_obs,
                  n_particles=1200, n_rounds=4, quantile=0.6, seed=5)
    unadj_mean = weighted_mean(res.thetas[:, 0], res.weights)
    adj = regression_adjust(res, s_obs)
    adj_mean = weighted_mean(adj[:, 0], res.weights)
    err_unadj = abs(unadj_mean - post_mean)
    err_adj = abs(adj_mean - post_mean)
    assert err_adj <= err_unadj
    assert err_adj < 0.03


def test_regression_adjust_shapes_and_subset():
    (ps, pl, sim, summ, dist, s_obs, _) = _binomial_rate_problem()
    res = abc_smc(ps, pl, sim, summ, dist, s_obs,
                  n_particles=300, n_rounds=3, seed=6)
    full = regression_adjust(res, s_obs)
    assert full.shape == res.thetas.shape
    sub = regression_adjust(res, s_obs, param_indices=[0])
    assert sub.shape == (res.thetas.shape[0], 1)


def test_bounded_adjustment_stays_inside_the_prior_box():
    """The support fix (F21). Measured on the real ABC-SMC fit, the unbounded
    local-linear adjustment put 300/800 innovation rates and 150/800 population
    sizes below zero. The discriminating probe is the pair: the SAME particles
    adjusted without bounds must leave the box, or this test is not testing the
    thing it is named for."""
    import numpy as np
    from mls_emergence.inference.abc_smc import SMCResult
    from mls_emergence.inference.regression import regression_adjust

    rng = np.random.default_rng(0)
    n, lo, hi = 300, np.array([0.001, 50.0]), np.array([0.05, 1000.0])
    thetas = rng.uniform(lo, hi, size=(n, 2))
    # summaries strongly correlated with theta and offset far from s_obs, which
    # is what makes the local-linear shift a long extrapolation
    summ = np.column_stack([thetas[:, 0] * 40.0, thetas[:, 1] / 400.0]) + \
        rng.normal(0, 0.02, (n, 2)) + 3.0
    res = SMCResult(thetas, np.full(n, 1.0 / n), summ, rng.random(n) + 0.1,
                    [1.0], [0.5], [n])
    s_obs = np.array([0.0, 0.0])

    unbounded = regression_adjust(res, s_obs)
    bounded = regression_adjust(res, s_obs, bounds=(lo, hi))

    n_out = int(((unbounded < lo) | (unbounded > hi)).sum())
    assert n_out > 0, "probe failed: unbounded adjustment stayed in the box, so this test cannot detect the defect"
    assert np.all(bounded >= lo) and np.all(bounded <= hi), \
        f"bounded adjustment left the prior box: {bounded.min(0)}, {bounded.max(0)}"
