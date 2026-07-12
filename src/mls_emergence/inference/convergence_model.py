"""Hierarchical Bayesian convergence model (PyMC).

Four emergence slopes with partial pooling and a heterogeneous likelihood: the
three CA-panel signatures contribute per-slice linear trajectories observed with
known measurement error; the seriation signature contributes its fragmentation
slope directly. The posterior yields P(all four slopes > 0) - the convergence
criterion as a single coherent quantity.
"""
from __future__ import annotations

import arviz as az
import numpy as np
import pymc as pm


def build_convergence_model(y, se, t, b_ser_obs, se_ser, *,
                            mu_sd=1.0, tau_sd=1.0, sigma_sd=0.5, a_sd=1.0):
    y = np.asarray(y, float)          # (n_panel, T)
    se = np.asarray(se, float)        # (n_panel, T)
    t = np.asarray(t, float)          # (T,)
    n_panel = y.shape[0]
    with pm.Model() as model:
        mu = pm.Normal("mu", 0.0, mu_sd)
        tau = pm.HalfNormal("tau", tau_sd)
        z_b = pm.Normal("z_b", 0.0, 1.0, shape=n_panel + 1)
        b = pm.Deterministic("b", mu + tau * z_b)        # last = seriation
        a = pm.Normal("a", 0.0, a_sd, shape=n_panel)
        sigma = pm.HalfNormal("sigma", sigma_sd)
        mu_panel = a[:, None] + b[:n_panel][:, None] * t[None, :]   # (n_panel, T)
        obs_sd = pm.math.sqrt(se ** 2 + sigma ** 2)
        pm.Normal("y_obs", mu=mu_panel, sigma=obs_sd, observed=y)
        pm.Normal("b_ser_obs", mu=b[n_panel], sigma=se_ser, observed=float(b_ser_obs))
    return model


def sample_convergence(y, se, t, b_ser_obs, se_ser, *, draws=2000, tune=2000,
                       chains=4, target_accept=0.9, random_seed=0, **priors):
    model = build_convergence_model(y, se, t, b_ser_obs, se_ser, **priors)
    with model:
        idata = pm.sample(draws=draws, tune=tune, chains=chains,
                          target_accept=target_accept, random_seed=random_seed,
                          progressbar=False)
    return idata


def convergence_summary(idata) -> dict:
    post = idata.posterior
    b = post["b"].stack(s=("chain", "draw")).values      # (n_panel+1, S)
    mu = post["mu"].stack(s=("chain", "draw")).values     # (S,)
    cscore = b.mean(axis=0)                               # mean slope per draw
    hdi = np.array([az.hdi(b[i], prob=0.95) for i in range(b.shape[0])])  # (n_panel+1, 2)
    lo = hdi[:, 0]
    hi = hdi[:, 1]
    return {
        "p_convergence": float(np.mean(np.all(b > 0, axis=0))),
        "b_mean": b.mean(axis=1).tolist(),
        "b_p_pos": (b > 0).mean(axis=1).tolist(),
        "b_hdi95": [[float(lo[i]), float(hi[i])] for i in range(b.shape[0])],
        "mu_mean": float(mu.mean()),
        "mu_p_pos": float((mu > 0).mean()),
        "cscore_slope_mean": float(cscore.mean()),
        "cscore_slope_p_pos": float((cscore > 0).mean()),
    }
