"""39_abc_smc_validation.py — validation battery for the ABC-SMC posterior.

Three checks, per the design spec:
  1. Simulation-based calibration (SBC): rank statistics of the true parameter
     within its posterior must be uniform under the prior.
  2. Coverage / power on b: 95% intervals over a grid of true b must cover truth
     ~95% of the time; interval width is the resolvable-effect-size story.
  3. Cross-check: the ABC-SMC posterior of b on the real data must agree with
     the rejection-ABC posterior from 19_abc_transmission.py.

Writes output/abc_smc_validation.md and figures/abc_smc_sbc,
figures/abc_smc_crosscheck.

Usage:
    .venv/bin/python analyses/39_abc_smc_validation.py          # full battery
    .venv/bin/python analyses/39_abc_smc_validation.py --fast   # smoke config
"""
from __future__ import annotations

import importlib
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from figstyle import OI_BLUE, OI_VERMIL, save  # noqa: E402

from mls_emergence.inference import (  # noqa: E402
    abc_smc, regression_adjust, weighted_quantile,
)

W = importlib.import_module("38_abc_smc_transmission")

OUT_MD = ROOT / "output" / "abc_smc_validation.md"

FULL = dict(n_sbc=200, sbc_particles=300, sbc_rounds=5,
            grid=np.linspace(-0.4, 0.4, 9), reps=20,
            cov_particles=300, cov_rounds=5)
FAST = dict(n_sbc=15, sbc_particles=50, sbc_rounds=2,
            grid=np.linspace(-0.4, 0.4, 3), reps=3,
            cov_particles=50, cov_rounds=2)


def posterior_rank(value, samples, weights):
    """Weighted fraction of posterior mass below ``value`` (in [0, 1])."""
    samples = np.asarray(samples, float)
    weights = np.asarray(weights, float)
    below = weights[samples < value].sum()
    return float(below / weights.sum())


def _fit_synthetic(theta_true, K, n_particles, n_rounds, seed_seq):
    """Simulate a dataset from theta_true and fit the whole-sequence target.

    seed_seq is a numpy SeedSequence; the observation RNG and the ABC-SMC
    sampler get independent child streams so no two replicates share entropy.
    """
    obs_ss, smc_ss = seed_seq.spawn(2)
    rng = np.random.default_rng(obs_ss)
    sim = W.make_simulator(K)
    summ = W.target_summary(slice(None))
    s_obs = summ(sim(theta_true, rng))
    result = abc_smc(W.prior_sampler, W.prior_logpdf, sim, summ, W.distance,
                     s_obs, n_particles=n_particles, n_rounds=n_rounds,
                     seed=smc_ss)
    adj = regression_adjust(result, s_obs)
    return result, adj


def _sbc_one(args):
    """One SBC replicate, fully determined by its index (parallel-safe)."""
    i, K, n_particles, n_rounds, base = args
    theta_ss, fit_ss = np.random.SeedSequence([int(base), 0, int(i)]).spawn(2)
    theta_true = W.prior_sampler(np.random.default_rng(theta_ss))
    result, adj = _fit_synthetic(theta_true, K, n_particles, n_rounds, fit_ss)
    return [posterior_rank(theta_true[j], adj[:, j], result.weights)
            for j in range(4)]


def sbc_ranks(n_sbc, K, n_particles, n_rounds, seed=11, n_workers=None):
    """Rank of each true parameter within its (regression-adjusted) posterior.

    Each replicate is independent and index-seeded, so the result is identical
    for any ``n_workers``. ``n_workers=1`` runs serially (no process pool).
    """
    args = [(i, K, n_particles, n_rounds, seed) for i in range(n_sbc)]
    if n_workers == 1:
        rows = [_sbc_one(a) for a in args]
    else:
        with Pool(processes=n_workers) as pool:
            rows = pool.map(_sbc_one, args)
    return np.array(rows)


def _cov_one(args):
    """One coverage replicate at a fixed true b, index-seeded (parallel-safe)."""
    gi, r, b_true, K, n_particles, n_rounds, base = args
    theta_ss, fit_ss = np.random.SeedSequence(
        [int(base), 1, int(gi), int(r)]).spawn(2)
    theta = W.prior_sampler(np.random.default_rng(theta_ss))
    theta[1] = b_true  # fix b, keep nuisance drawn from prior
    result, adj = _fit_synthetic(theta, K, n_particles, n_rounds, fit_ss)
    lo = weighted_quantile(adj[:, 1], result.weights, 0.025)
    hi = weighted_quantile(adj[:, 1], result.weights, 0.975)
    return (gi, int(lo <= b_true <= hi), hi - lo)


def coverage_on_b(b_grid, reps, K, n_particles, n_rounds, seed=22,
                  n_workers=None):
    """Fraction of 95% intervals for b that cover truth, and mean width.

    Index-seeded and parallel-safe like :func:`sbc_ranks`.
    """
    args = [(gi, r, float(b), K, n_particles, n_rounds, seed)
            for gi, b in enumerate(b_grid) for r in range(reps)]
    if n_workers == 1:
        out = [_cov_one(a) for a in args]
    else:
        with Pool(processes=n_workers) as pool:
            out = pool.map(_cov_one, args)
    n = len(b_grid)
    cover = np.zeros(n)
    width = np.zeros(n)
    cnt = np.zeros(n)
    for gi, hit, wd in out:
        cover[gi] += hit
        width[gi] += wd
        cnt[gi] += 1
    return {"grid": np.asarray(b_grid), "coverage": cover / cnt,
            "mean_width": width / cnt}


def cross_check():
    """Compare ABC-SMC and rejection-ABC posteriors of b on the real data."""
    smc = np.load(ROOT / "output" / "abc_smc_posterior.npz")["post_b"]
    rej = np.load(ROOT / "output" / "abc_posterior.npz")["post_b"]
    return {
        "smc_mean": float(smc.mean()),
        "smc_ci": (float(np.percentile(smc, 2.5)), float(np.percentile(smc, 97.5))),
        "rej_mean": float(rej.mean()),
        "rej_ci": (float(np.percentile(rej, 2.5)), float(np.percentile(rej, 97.5))),
        "mean_gap": float(abs(smc.mean() - rej.mean())),
    }


def main(fast=False, n_workers=None):
    cfg = FAST if fast else FULL
    # K must match the real inference (analysis 38): derive the decorated-type
    # count from the curated basin set rather than hard-coding it, so SBC and
    # coverage validate the same-sized problem we actually solved.
    obs, _ = W.load_obs()
    K = obs.shape[1]
    if n_workers is None:
        n_workers = os.cpu_count()

    # 1. SBC
    ranks = sbc_ranks(cfg["n_sbc"], K, cfg["sbc_particles"], cfg["sbc_rounds"],
                      seed=11, n_workers=n_workers)
    # persist the ranks so the manuscript figure (Figure S6) can show the b-rank
    # calibration histogram without recomputing the (expensive) SBC.
    np.savez(ROOT / "output" / "abc_smc_sbc_ranks.npz", ranks=ranks)
    names = ["mu", "b", "N", "w"]
    # chi-square uniformity test on the b ranks
    from scipy.stats import kstest
    ks_b = kstest(ranks[:, 1], "uniform")

    fig, axes = plt.subplots(1, 4, figsize=(7, 2))
    for j, ax in enumerate(axes):
        ax.hist(ranks[:, j], bins=10, range=(0, 1), color=OI_BLUE)
        ax.set_xlabel(f"SBC rank ({names[j]})")
        ax.set_yticks([])
    save(fig, "abc_smc_sbc")
    plt.close(fig)

    # 2. Coverage on b
    cov = coverage_on_b(cfg["grid"], cfg["reps"], K, cfg["cov_particles"],
                        cfg["cov_rounds"], seed=22, n_workers=n_workers)

    # 3. Cross-check (only if both posterior files exist)
    cc = None
    if (ROOT / "output" / "abc_smc_posterior.npz").exists() and \
       (ROOT / "output" / "abc_posterior.npz").exists():
        cc = cross_check()
        fig, ax = plt.subplots(figsize=(3.5, 3))
        smc = np.load(ROOT / "output" / "abc_smc_posterior.npz")["post_b"]
        rej = np.load(ROOT / "output" / "abc_posterior.npz")["post_b"]
        ax.hist(rej, bins=40, density=True, alpha=0.5, color=OI_VERMIL,
                label="rejection ABC (19)")
        ax.hist(smc, bins=40, density=True, alpha=0.5, color=OI_BLUE,
                label="ABC-SMC (38)")
        ax.set_xlabel("transmission bias b")
        ax.legend()
        save(fig, "abc_smc_crosscheck")
        plt.close(fig)

    L = ["# ABC-SMC validation battery", "",
         f"Config: {'SMOKE' if fast else 'full'}, {n_workers} worker process(es). "
         f"SBC n={cfg['n_sbc']} (particles {cfg['sbc_particles']}, rounds "
         f"{cfg['sbc_rounds']}); coverage grid {list(np.round(cfg['grid'],2))} "
         f"x {cfg['reps']} reps.", "",
         "## 1. Simulation-based calibration", "",
         f"- b rank KS vs uniform: D = {ks_b.statistic:.3f}, p = {ks_b.pvalue:.3f} "
         "(p > 0.05 is consistent with calibration).",
         "- Rank histograms: figures/abc_smc_sbc (should be flat).", "",
         "## 2. Coverage / power on b", ""]
    for g, c, wd in zip(cov["grid"], cov["coverage"], cov["mean_width"]):
        L.append(f"- true b = {g:+.2f}: 95% coverage {c:.2f}, mean width {wd:.3f}.")
    L += ["", "## 3. Cross-check vs rejection ABC (analysis 19)", ""]
    if cc is None:
        L.append("- SKIPPED: run analyses/19 and analyses/38 first to produce "
                 "both posterior files.")
    else:
        L += [f"- ABC-SMC b = {cc['smc_mean']:+.3f} "
              f"[{cc['smc_ci'][0]:+.3f}, {cc['smc_ci'][1]:+.3f}].",
              f"- Rejection ABC b = {cc['rej_mean']:+.3f} "
              f"[{cc['rej_ci'][0]:+.3f}, {cc['rej_ci'][1]:+.3f}].",
              f"- Posterior-mean gap = {cc['mean_gap']:.3f} "
              "(should be within Monte-Carlo error).",
              "- Overlay: figures/abc_smc_crosscheck."]

    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
