"""Bayesian cultural F_ST for the decorated-ceramic clusters (Balding-Nichols).

This mirrors the hyperlocality project's ``src/bayes.py`` so the two papers use
the same basic model, adapted to the mls-emergence data. It replaces the earlier
non-centered Dirichlet-multinomial (which put an informative prior on the
concentration and sampled the full per-group simplex, producing a low-F_ST
funnel) with the standard population-genetics parameterization.

Model (Balding & Nichols 1995; the cultural reading follows the same
classes-as-alleles / assemblages-as-demes move as the frequentist ``variance``
signature). For one locus with K decorated classes across G spatial clusters:

    F   ~ Uniform(0, 1)                     cultural F_ST (the BN theta); a=(1-F)/F
    pi  ~ Dirichlet(1, ..., 1)              basin-wide class frequencies
    x_g ~ DirichletMultinomial(n_g, a*pi)   cluster g's observed counts

F_ST is a model PARAMETER estimated directly, with a FLAT prior on F_ST itself
(so the prior does not push the estimate toward the drift level the paper reports).
The per-cluster frequencies are marginalized analytically inside the
Dirichlet-multinomial, so the sampler sees only F and pi: this is well
conditioned and needs no funnel-taming ``target_accept``.

Two readouts come out of the one fit:

  * ``fst_summary``            -- posterior of the BN F_ST parameter (median, 95%
                                  HDI, R-hat, ESS).
  * ``gini_simpson_summary``   -- posterior of the SAME Gini-Simpson estimator the
                                  manuscript reports, (H_T-H_S)/H_T, reconstructed
                                  conjugately for the OBSERVED clusters
                                  (p_g | x_g ~ Dirichlet(a*pi + x_g)). This keeps a
                                  credible interval on the exact reported quantity.

Structure vs panmixia is tested by a Bayes factor (``bayes_factor_structure``):
panmixia (M0) is the F -> 0 limit where every cluster shares one pi and the
counts are plain Multinomial. Because M0 is nested in M1, the SMC marginal
likelihoods are directly comparable and 2 ln BF10 is read on the Kass & Raftery
(1995) scale. This quantifies estimation uncertainty and evidence for structure;
it does not replace the separate stochastic-drift-null comparison, which asks
whether spatial drift alone could produce the value.
"""
from __future__ import annotations

import contextlib
import io
import logging

import numpy as np
import pymc as pm

# PyMC/SMC emit sampling chatter through the logging module (bypassing stdout
# redirection); quiet it so analysis-script output stays readable. Diagnostics
# (R-hat, ESS) are computed explicitly instead.
for _name in ("pymc", "pymc.sampling", "pymc.smc"):
    logging.getLogger(_name).setLevel(logging.ERROR)


def fst_from_frequencies(p, sizes) -> float:
    """Gini-Simpson F_ST for one (G, K) frequency array with group sizes (G,).

    F_ST = (H_T - H_S) / H_T on Gini-Simpson diversity H = 1 - sum p^2, with H_S
    the size-weighted mean within-group diversity and H_T the pooled diversity.
    This is the exact estimator ``signatures.variance.cultural_fst`` computes on
    counts; kept here for the validation battery and the conjugate readout.
    """
    p = np.asarray(p, float)
    sizes = np.asarray(sizes, float)
    w = sizes / sizes.sum()
    h_within = 1.0 - np.sum(p ** 2, axis=1)          # (G,)
    h_s = float(np.sum(h_within * w))
    p_pool = np.sum(p * w[:, None], axis=0)          # (K,)
    h_t = float(1.0 - np.sum(p_pool ** 2))
    return (h_t - h_s) / h_t if h_t > 0 else 0.0


def _hdi95(samples, hdi_prob: float = 0.95):
    """95% HDI of a 1-D sample array, robust across arviz keyword conventions."""
    import arviz as az
    x = np.asarray(samples, float).ravel()
    for kw in ("prob", "hdi_prob"):
        try:
            hdi = az.hdi(x, **{kw: hdi_prob})
            return [float(hdi[0]), float(hdi[1])]
        except TypeError:
            continue
    q = [100 * (1 - hdi_prob) / 2, 100 * (1 + hdi_prob) / 2]
    lo, hi = np.percentile(x, q)
    return [float(lo), float(hi)]


@contextlib.contextmanager
def _quiet():
    """Silence PyMC sampling chatter so run-script output stays readable."""
    import sys
    sys.stdout.flush()
    sys.stderr.flush()
    with contextlib.redirect_stdout(io.StringIO()), \
            contextlib.redirect_stderr(io.StringIO()):
        yield


def _prep(counts):
    """Drop empty clusters, return integer (G, K) counts and row sums (G,)."""
    counts = np.asarray(counts, float)
    counts = counts[counts.sum(1) > 0]
    counts = np.rint(counts).astype(int)
    return counts, counts.sum(1)


# --------------------------------------------------------------------------- #
#  Balding-Nichols model (structured M1) and the panmixia limit (M0)          #
# --------------------------------------------------------------------------- #
def build_fst_model(group_counts, *, f_prior=("uniform",)):
    """Structured (M1) Balding-Nichols model. F is the cultural F_ST parameter.

    ``f_prior`` is ('uniform',) or ('beta', a, b) for the prior-sensitivity check.
    The per-cluster frequencies are marginalized inside the Dirichlet-multinomial.
    """
    counts, n_g = _prep(group_counts)
    K = counts.shape[1]
    model = pm.Model()
    with model:
        if f_prior[0] == "uniform":
            F = pm.Uniform("F", 0.0, 1.0)
        elif f_prior[0] == "beta":
            F = pm.Beta("F", alpha=f_prior[1], beta=f_prior[2])
        else:
            raise ValueError(f"unknown F prior {f_prior!r}")
        alpha = (1.0 - F) / F
        pi = pm.Dirichlet("pi", a=np.ones(K))
        pm.DirichletMultinomial("x", n=n_g, a=alpha * pi, observed=counts)
    return model


def build_panmixia_model(group_counts):
    """Panmixia (M0): all clusters share one pi (the F -> 0 limit)."""
    counts, n_g = _prep(group_counts)
    K = counts.shape[1]
    model = pm.Model()
    with model:
        pi = pm.Dirichlet("pi", a=np.ones(K))
        pm.Multinomial("x", n=n_g, p=pi, observed=counts)
    return model


def sample_fst(group_counts, *, draws: int = 2000, tune: int = 2000, chains: int = 4,
               target_accept: float = 0.9, random_seed: int = 0, f_prior=("uniform",)):
    """Posterior for the cultural F_ST parameter F. Returns an ArviZ InferenceData.

    The marginalized model is well conditioned, so the default ``target_accept``
    of 0.9 suffices (no funnel-taming needed).
    """
    model = build_fst_model(group_counts, f_prior=f_prior)
    with model, _quiet():
        idata = pm.sample(draws=draws, tune=tune, chains=chains, cores=1,
                          target_accept=target_accept, random_seed=random_seed,
                          progressbar=False, compute_convergence_checks=False)
    return idata


def fst_summary(idata, *, hdi_prob: float = 0.95) -> dict:
    """Posterior of the BN cultural F_ST parameter: median, mean, 95% HDI, R-hat, ESS."""
    import arviz as az
    f = np.asarray(idata.posterior["F"].values).ravel()
    ess = float(np.asarray(az.ess(idata, var_names=["F"])["F"].values))
    rhat = float(np.asarray(az.rhat(idata, var_names=["F"])["F"].values))
    return {
        "fst_median": float(np.median(f)),
        "fst_mean": float(np.mean(f)),
        "fst_hdi95": _hdi95(f, hdi_prob),
        "rhat": rhat,
        "ess": ess,
        "fst_samples": f,
    }


# --------------------------------------------------------------------------- #
#  Gini-Simpson readout: credible interval on the exact manuscript estimator  #
# --------------------------------------------------------------------------- #
def gini_simpson_posterior(idata, group_counts, *, seed: int = 0) -> np.ndarray:
    """Posterior of the Gini-Simpson F_ST for the OBSERVED clusters.

    Reconstructs each cluster's frequency conjugately from the fitted BN
    posterior: with likelihood DirichletMultinomial(n_g, a*pi) and the latent
    p_g ~ Dirichlet(a*pi), the posterior of the observed cluster's frequency is
    p_g | x_g, F, pi ~ Dirichlet(a*pi + x_g). For each posterior draw of (F, pi)
    we draw p_g for every cluster and compute the size-weighted Gini-Simpson
    F_ST = (H_T - H_S)/H_T. This yields a credible interval on the exact estimator
    the manuscript reports, from the same well-conditioned fit.
    """
    counts, sizes = _prep(group_counts)
    K = counts.shape[1]
    w = sizes / sizes.sum()

    F = np.asarray(idata.posterior["F"].values).ravel()          # (S,)
    pi = idata.posterior["pi"].values.reshape(-1, K)              # (S, K)
    alpha = (1.0 - F) / F                                         # (S,)

    rng = np.random.default_rng(seed)
    conc = alpha[:, None, None] * pi[:, None, :] + counts[None, :, :]   # (S, G, K)
    gam = rng.standard_gamma(conc)                               # (S, G, K)
    p = gam / gam.sum(axis=2, keepdims=True)                     # (S, G, K)

    h_within = 1.0 - np.sum(p ** 2, axis=2)                      # (S, G)
    h_s = np.sum(h_within * w[None, :], axis=1)                  # (S,)
    p_pool = np.sum(p * w[None, :, None], axis=1)                # (S, K)
    h_t = 1.0 - np.sum(p_pool ** 2, axis=1)                      # (S,)
    return np.where(h_t > 0, (h_t - h_s) / h_t, 0.0)


def gini_simpson_summary(idata, group_counts, *, seed: int = 0,
                         hdi_prob: float = 0.95) -> dict:
    """Summary of the Gini-Simpson F_ST posterior plus the frequentist plug-in."""
    from mls_emergence.signatures.variance import cultural_fst
    gst = gini_simpson_posterior(idata, group_counts, seed=seed)
    return {
        "gst_median": float(np.median(gst)),
        "gst_mean": float(np.mean(gst)),
        "gst_hdi95": _hdi95(gst, hdi_prob),
        "plugin_fst": float(cultural_fst(np.asarray(group_counts, float))),
        "gst_samples": gst,
    }


# --------------------------------------------------------------------------- #
#  Bayes factor: structure (M1) vs panmixia (M0) -- the drift/null-free test  #
# --------------------------------------------------------------------------- #
def _final_logml(idata) -> np.ndarray:
    """Per-chain full-data log marginal likelihood from an SMC InferenceData."""
    arr = np.asarray(idata.sample_stats["log_marginal_likelihood"].values).ravel()
    vals = []
    for entry in arr:
        seq = np.asarray(entry, dtype=float).ravel()
        seq = seq[~np.isnan(seq)]
        if seq.size:
            vals.append(seq[-1])
    return np.asarray(vals, dtype=float)


def _smc_logml(group_counts, *, structured: bool, draws: int, chains: int,
               seed: int, f_prior=("uniform",)) -> np.ndarray:
    model = (build_fst_model(group_counts, f_prior=f_prior) if structured
             else build_panmixia_model(group_counts))
    with model, _quiet():
        idata = pm.sample_smc(draws=draws, chains=chains, cores=1,
                              random_seed=seed, progressbar=False)
    return _final_logml(idata)


def _kass_raftery(two_ln_bf: float) -> str:
    """Kass & Raftery (1995) evidence label for 2 ln BF10 (structure vs panmixia)."""
    a = abs(two_ln_bf)
    band = ("not worth more than a bare mention" if a < 2 else
            "positive" if a < 6 else
            "strong" if a < 10 else
            "very strong")
    direction = "for structure" if two_ln_bf > 0 else "for panmixia"
    return f"{band} {direction}"


def bayes_factor_structure(group_counts, *, draws: int = 2000, chains: int = 4,
                           seed: int = 0, f_prior=("uniform",)) -> dict:
    """Bayes factor for between-cluster structure (M1) vs panmixia (M0).

    Fits both models with SMC and returns log marginal likelihoods, log10 BF10,
    2 ln BF10, the Kass-Raftery label, and the across-chain spread of 2 ln BF10 as
    a stability check.
    """
    lml1 = _smc_logml(group_counts, structured=True, draws=draws, chains=chains,
                      seed=seed, f_prior=f_prior)
    lml0 = _smc_logml(group_counts, structured=False, draws=draws, chains=chains,
                      seed=seed)
    m1, m0 = float(lml1.mean()), float(lml0.mean())
    ln_bf = m1 - m0
    two_ln_bf = 2.0 * ln_bf
    chain_sd = float(2.0 * np.sqrt(lml1.var() / len(lml1) + lml0.var() / len(lml0)))
    return {
        "logml_structure": m1,
        "logml_panmixia": m0,
        "log10_bf": float(ln_bf / np.log(10)),
        "two_ln_bf": two_ln_bf,
        "two_ln_bf_chain_sd": chain_sd,
        "evidence": _kass_raftery(two_ln_bf),
    }
