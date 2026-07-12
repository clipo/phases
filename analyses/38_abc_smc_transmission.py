"""38_abc_smc_transmission.py — ABC-SMC posterior of the transmission bias.

Same generative copying-with-time-averaging model as 19_abc_transmission.py,
fit with sequential Monte-Carlo ABC and local-linear regression adjustment
instead of plain rejection ABC. Runs three targets (whole sequence, early half,
late half) and writes output/abc_smc_transmission.md + output/abc_smc_posterior.npz.
19_abc_transmission.py is left untouched as the cross-check baseline.

Usage:
    .venv/bin/python analyses/38_abc_smc_transmission.py          # full run
    .venv/bin/python analyses/38_abc_smc_transmission.py --fast   # smoke config
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.inference import (  # noqa: E402
    abc_smc, regression_adjust, weighted_mean, weighted_quantile, resample,
)

abc19 = importlib.import_module("19_abc_transmission")
import make_figures as mf  # noqa: E402
res = importlib.import_module("17_basin_results")

N_BINS_LOCAL = abc19.N_BINS
OUT_MD = ROOT / "output" / "abc_smc_transmission.md"
OUT_NPZ = ROOT / "output" / "abc_smc_posterior.npz"

# priors match 19_abc_transmission.py exactly; theta order = (mu, b, N, w)
PRIOR_LO = np.array([0.001, -0.5, 50.0, 5.0])
PRIOR_HI = np.array([0.05, 0.5, 1000.0, 30.0])

# full-run and smoke settings
FULL = dict(n_particles=800, n_rounds=6)
FAST = dict(n_particles=60, n_rounds=2)


def prior_sampler(rng):
    return rng.uniform(PRIOR_LO, PRIOR_HI)


def prior_logpdf(theta):
    theta = np.asarray(theta, float)
    if np.all(theta >= PRIOR_LO) and np.all(theta <= PRIOR_HI):
        return 0.0
    return -np.inf


def make_simulator(K):
    def simulator(theta, rng):
        mu, b, N, w = theta
        return abc19.simulate(mu, b, int(round(N)), w, K, rng)
    return simulator


def target_summary(sl):
    def summ(sim):
        return abc19.summary(sim[sl])
    return summ


def distance(s, s_obs):
    return float(np.sqrt(np.sum((np.asarray(s) - np.asarray(s_obs)) ** 2)))


def load_obs():
    counts, coords = mf._load_curated()
    ca, _ = res.oriented_ca(counts)
    rank = ca.rank().to_numpy()
    bins = pd.qcut(pd.Series(rank, index=counts.index), N_BINS_LOCAL,
                   labels=False, duplicates="drop")
    obs = np.array([counts.loc[bins.index[bins == bb]].to_numpy(float).sum(0)
                    for bb in sorted(pd.Series(bins).dropna().unique())])
    return obs, counts


def fit_target(obs, sl, n_particles, n_rounds, seed):
    """Fit one target (a bin slice). Returns (SMCResult, adjusted b samples)."""
    K = obs.shape[1]
    summ = target_summary(sl)
    s_obs = summ(obs)
    result = abc_smc(prior_sampler, prior_logpdf, make_simulator(K),
                     summ, distance, s_obs,
                     n_particles=n_particles, n_rounds=n_rounds, seed=seed)
    adj = regression_adjust(result, s_obs)   # (n_particles, 4)
    return result, adj[:, 1]                  # column 1 = b


def _summarize(result, adj_b):
    lo = weighted_quantile(adj_b, result.weights, 0.025)
    hi = weighted_quantile(adj_b, result.weights, 0.975)
    mean = weighted_mean(adj_b, result.weights)
    # posterior SD from an equally-weighted resample of the adjusted particles
    rng = np.random.default_rng(0)
    draws = resample(adj_b, result.weights, 20000, rng)
    return mean, lo, hi, float(np.mean(draws > 0)), float(np.std(draws)), draws


def main(fast=False):
    cfg = FAST if fast else FULL
    obs, counts = load_obs()
    K = obs.shape[1]

    r_all, b_all = fit_target(obs, slice(None), seed=101, **cfg)
    r_e, b_e = fit_target(obs, slice(0, 3), seed=102, **cfg)
    r_l, b_l = fit_target(obs, slice(3, 6), seed=103, **cfg)

    mean, lo, hi, p_pos, post_sd, draws_all = _summarize(r_all, b_all)
    _, e_lo, e_hi, _, _, draws_e = _summarize(r_e, b_e)
    _, l_lo, l_hi, _, _, draws_l = _summarize(r_l, b_l)
    e_mean = weighted_mean(b_e, r_e.weights)
    l_mean = weighted_mean(b_l, r_l.weights)
    prior_sd = (PRIOR_HI[1] - PRIOR_LO[1]) / np.sqrt(12)

    L = ["# ABC-SMC inference of the transmission bias", "",
         f"Basin curated set (n = {counts.shape[0]}), {K} decorated types, "
         f"{N_BINS_LOCAL} ordinal bins. ABC-SMC with local-linear regression "
         f"adjustment, {cfg['n_particles']} particles, {cfg['n_rounds']} rounds "
         f"({'SMOKE' if fast else 'full'} config). Same model and summary "
         "statistics as 19_abc_transmission.py.", "",
         "## Posterior of the transmission-bias parameter b", "",
         f"- Posterior mean b = {mean:+.3f}, 95% interval [{lo:+.3f}, {hi:+.3f}].",
         f"- The 95% interval {'INCLUDES' if lo <= 0 <= hi else 'EXCLUDES'} the "
         "neutral value b = 0.",
         f"- P(b > 0) = {p_pos:.2f} (0.5 = no directional information).",
         f"- Posterior SD {post_sd:.3f} vs prior SD {prior_sd:.3f} "
         f"(ratio {post_sd/prior_sd:.2f}).", "",
         "## Early vs late halves", "",
         f"- Early-half b = {e_mean:+.3f} [{e_lo:+.3f}, {e_hi:+.3f}].",
         f"- Late-half b = {l_mean:+.3f} [{l_lo:+.3f}, {l_hi:+.3f}].",
         f"- Shift early->late: {l_mean - e_mean:+.3f}.", "",
         "## Round diagnostics (whole-sequence target)", "",
         f"- Tolerance schedule: {[round(e,4) for e in r_all.eps_schedule]}.",
         f"- Acceptance rates: {[round(a,4) for a in r_all.accept_rates]}.",
         f"- Simulations per round: {r_all.n_sims}.", ""]

    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    np.savez(OUT_NPZ, post_b=draws_all, early=draws_e, late=draws_l)
    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
