"""Bayesian credible interval on the cultural F_ST (Gini-Simpson estimator).

A hierarchical Dirichlet-multinomial generative model on group sherd counts:
each group's latent type-frequency vector is drawn from a Dirichlet centered on a
shared ancestral distribution with a concentration parameter, and the observed
counts are multinomial. For each posterior draw of the group frequencies we
compute the SAME estimator the manuscript reports, (H_T - H_S)/H_T on
Gini-Simpson diversity, so the posterior is a credible interval on the reported
F_ST with small-sample bias correction. This quantifies estimation uncertainty;
it does not replace the separate stochastic-drift-null comparison.
"""
from __future__ import annotations

import numpy as np
import pymc as pm


def fst_from_frequencies(p, sizes) -> float:
    """Gini-Simpson F_ST for one (G, K) frequency array with group sizes (G,)."""
    p = np.asarray(p, float)
    sizes = np.asarray(sizes, float)
    w = sizes / sizes.sum()
    h_within = 1.0 - np.sum(p ** 2, axis=1)          # (G,)
    h_s = float(np.sum(h_within * w))
    p_pool = np.sum(p * w[:, None], axis=0)          # (K,)
    h_t = float(1.0 - np.sum(p_pool ** 2))
    return (h_t - h_s) / h_t if h_t > 0 else 0.0


def build_fst_model(group_counts, *, conc_mu: float = 3.0, conc_sd: float = 2.0):
    gc = np.asarray(group_counts)
    counts = gc.astype(int)
    G, K = counts.shape
    N = counts.sum(axis=1)
    with pm.Model() as model:
        p_anc = pm.Dirichlet("p_anc", np.ones(K))
        log_conc = pm.Normal("log_conc", conc_mu, conc_sd)     # concentration on log scale
        conc = pm.Deterministic("conc", pm.math.exp(log_conc))
        p = pm.Dirichlet("p", conc * p_anc, shape=(G, K))
        pm.Multinomial("counts", n=N, p=p, observed=counts)
    return model


def sample_fst(group_counts, *, draws: int = 2000, tune: int = 2000, chains: int = 4,
               target_accept: float = 0.9, random_seed: int = 0, **priors):
    model = build_fst_model(group_counts, **priors)
    with model:
        idata = pm.sample(draws=draws, tune=tune, chains=chains,
                          target_accept=target_accept, random_seed=random_seed,
                          progressbar=False)
    return idata


def fst_posterior(idata, sizes) -> np.ndarray:
    """Posterior sample of the Gini-Simpson F_ST from the group-frequency draws."""
    sizes = np.asarray(sizes, float)
    w = sizes / sizes.sum()
    p = idata.posterior["p"].values                  # (chain, draw, G, K)
    p = p.reshape(-1, p.shape[-2], p.shape[-1])       # (S, G, K)
    h_within = 1.0 - np.sum(p ** 2, axis=2)           # (S, G)
    h_s = np.sum(h_within * w[None, :], axis=1)       # (S,)
    p_pool = np.sum(p * w[None, :, None], axis=1)      # (S, K)
    h_t = 1.0 - np.sum(p_pool ** 2, axis=1)           # (S,)
    return np.where(h_t > 0, (h_t - h_s) / h_t, 0.0)


def fst_summary(idata, sizes, group_counts) -> dict:
    import arviz as az
    from mls_emergence.signatures.variance import cultural_fst
    fst = fst_posterior(idata, sizes)
    try:
        hdi = az.hdi(fst, prob=0.95)                    # arviz 1.2: kwarg is `prob`
        hdi = [float(hdi[0]), float(hdi[1])]
    except TypeError:
        hdi = [float(v) for v in np.percentile(fst, [2.5, 97.5])]
    return {
        "fst_mean": float(fst.mean()),
        "fst_hdi95": hdi,
        "p_fst_gt0": float((fst > 0).mean()),
        "plugin_fst": float(cultural_fst(np.asarray(group_counts, float))),
        "fst_samples": fst,
    }
