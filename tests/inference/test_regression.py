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
