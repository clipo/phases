"""44_bayesian_fst_validation.py - Validation battery for the BN F_ST fit.

Checks on the Balding-Nichols cultural-F_ST model (analysis 43):

1. Coverage of the BN F_ST parameter: at known true F, does the 95% credible
   interval contain it ~95% of the time over simulated datasets?
2. Coverage of the Gini-Simpson readout: does its 95% interval contain the
   REALIZED Gini-Simpson F_ST of each simulated dataset?
3. Small-sample bias: at small group sizes, does the posterior track truth better
   than the frequentist plug-in Gini-Simpson estimator?
4. Simulation-based calibration (SBC): over datasets drawn from the prior, is the
   rank of the true F within its posterior uniform (Talts et al. 2018)?
5. MCMC health on the real St. Francis basin fit: max R-hat, min ESS, divergences.
6. Bayes-factor sanity: panmixia-generated data favors panmixia; strongly
   structured data favors structure.

The BN model is well conditioned (marginalized frequencies), so no funnel-taming
target_accept is needed. Datasets are simulated from the BN generative process
itself, so coverage/SBC test the estimand the model targets.

Usage: .venv/bin/python analyses/44_bayesian_fst_validation.py [--fast] [--serial]
"""
from __future__ import annotations

import importlib
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
from figstyle import OI_BLUE, OI_VERMIL, OI_GREEN, save  # noqa: E402

OUT_MD = ROOT / "output" / "bayesian_fst_validation.md"

FULL = dict(f_grid=(0.02, 0.05, 0.10, 0.20, 0.35, 0.50), reps=40,
            draws=1000, tune=1000, chains=2, n_sbc=200, sbc_bins=20,
            small_sizes=(15, 15, 15))
FAST = dict(f_grid=(0.05, 0.30), reps=4,
            draws=300, tune=400, chains=2, n_sbc=12, sbc_bins=4,
            small_sizes=(15, 15, 15))

P_ANC = np.array([0.4, 0.35, 0.25])
SIZES = np.array([200, 200, 200])


def simulate_bn(f_true, p_anc, sizes, rng):
    """Draw (G, K) counts from the Balding-Nichols process at true F_ST f_true.

    alpha = (1 - f_true) / f_true; for each group p_g ~ Dirichlet(alpha * p_anc),
    counts_g ~ Multinomial(sizes[g], p_g). Returns the counts, the true F, and the
    realized Gini-Simpson F_ST of the drawn frequencies.
    """
    from mls_emergence.inference import fst_from_frequencies
    f_true = float(np.clip(f_true, 1e-3, 1 - 1e-3))
    p_anc = np.asarray(p_anc, float)
    sizes = np.asarray(sizes, int)
    alpha = (1.0 - f_true) / f_true
    G, K = len(sizes), len(p_anc)
    p = np.zeros((G, K))
    counts = np.zeros((G, K), dtype=int)
    for g in range(G):
        p[g] = rng.dirichlet(alpha * p_anc)
        counts[g] = rng.multinomial(sizes[g], p[g])
    return counts, f_true, fst_from_frequencies(p, sizes)


# --------------------------------------------------------------------------- #
#  Coverage + small-sample bias (parallel, index-seeded)                      #
# --------------------------------------------------------------------------- #
def _cov_one(args):
    """One coverage replicate at fixed true F, index-seeded (parallel-safe)."""
    gi, r, f_true, p_anc, sizes, draws, tune, chains, base = args
    from mls_emergence.inference import (
        sample_fst, fst_summary, gini_simpson_summary,
    )
    obs_ss, fit_ss = np.random.SeedSequence([int(base), int(gi), int(r)]).spawn(2)
    rng = np.random.default_rng(obs_ss)
    gc, f_tru, g_realized = simulate_bn(f_true, p_anc, np.asarray(sizes), rng)
    idata = sample_fst(gc, draws=draws, tune=tune, chains=chains,
                       random_seed=int(fit_ss.generate_state(1)[0]))
    s = fst_summary(idata)
    g = gini_simpson_summary(idata, gc, seed=0)
    cov_f = int(s["fst_hdi95"][0] <= f_tru <= s["fst_hdi95"][1])
    cov_g = int(g["gst_hdi95"][0] <= g_realized <= g["gst_hdi95"][1])
    return {
        "gi": gi,
        "cov_f": cov_f,
        "cov_g": cov_g,
        "realized": g_realized,
        "plugin_bias": g["plugin_fst"] - g_realized,
        "post_bias": g["gst_median"] - g_realized,
        "abs_plugin_bias": abs(g["plugin_fst"] - g_realized),
        "abs_post_bias": abs(g["gst_median"] - g_realized),
    }


def coverage(f_grid, reps, p_anc, sizes, *, draws, tune, chains, seed, n_workers):
    args = [(gi, r, float(f), tuple(p_anc), tuple(int(x) for x in sizes),
             draws, tune, chains, seed)
            for gi, f in enumerate(f_grid) for r in range(reps)]
    rows = ([_cov_one(a) for a in args] if n_workers == 1
            else Pool(processes=n_workers).map(_cov_one, args))
    out = []
    for gi, f in enumerate(f_grid):
        rs = [x for x in rows if x["gi"] == gi]
        out.append({
            "f_true": float(f),
            "reps": len(rs),
            "mean_realized": float(np.mean([x["realized"] for x in rs])),
            "coverage_f": float(np.mean([x["cov_f"] for x in rs])),
            "coverage_g": float(np.mean([x["cov_g"] for x in rs])),
            "mean_plugin_bias": float(np.mean([x["plugin_bias"] for x in rs])),
            "mean_post_bias": float(np.mean([x["post_bias"] for x in rs])),
            "abs_plugin_bias": float(np.mean([x["abs_plugin_bias"] for x in rs])),
            "abs_post_bias": float(np.mean([x["abs_post_bias"] for x in rs])),
        })
    return out


# --------------------------------------------------------------------------- #
#  SBC on the BN F parameter (parallel, index-seeded)                         #
# --------------------------------------------------------------------------- #
def _sbc_one(args):
    """One SBC replicate, drawing the FULL prior (F and pi), then ranking F.

    Proper simulation-based calibration requires the simulated data to come from
    the model's joint prior predictive. The model has F ~ Uniform(0,1) and
    pi ~ Dirichlet(1,...,1), so both are drawn fresh here (holding pi fixed would
    break rank uniformity).
    """
    i, K, sizes, draws, tune, chains, base = args
    from mls_emergence.inference import sample_fst
    f_ss, pi_ss, obs_ss, fit_ss = np.random.SeedSequence([int(base), int(i)]).spawn(4)
    f_true = float(np.random.default_rng(f_ss).uniform(0.0, 1.0))
    pi_true = np.random.default_rng(pi_ss).dirichlet(np.ones(K))
    rng = np.random.default_rng(obs_ss)
    gc, f_tru, _ = simulate_bn(f_true, pi_true, np.asarray(sizes), rng)
    idata = sample_fst(gc, draws=draws, tune=tune, chains=chains,
                       random_seed=int(fit_ss.generate_state(1)[0]))
    f_samples = np.asarray(idata.posterior["F"].values).ravel()
    return float((f_samples < f_tru).mean())


def sbc_ranks(n_sbc, K, sizes, *, draws, tune, chains, seed, n_workers):
    args = [(i, int(K), tuple(int(x) for x in sizes), draws, tune, chains, seed)
            for i in range(n_sbc)]
    if n_workers == 1:
        return np.array([_sbc_one(a) for a in args])
    return np.array(Pool(processes=n_workers).map(_sbc_one, args))


def sbc_uniformity(ranks, n_bins):
    """Chi-square test that SBC ranks are uniform on [0,1]."""
    from scipy import stats
    counts, _ = np.histogram(ranks, bins=n_bins, range=(0, 1))
    expected = len(ranks) / n_bins
    chi2 = float(((counts - expected) ** 2 / expected).sum())
    p = float(stats.chi2.sf(chi2, n_bins - 1))
    ks = stats.kstest(ranks, "uniform")
    return {"chi2": chi2, "chi2_p": p, "ks_stat": float(ks.statistic),
            "ks_p": float(ks.pvalue), "counts": counts.tolist(), "expected": expected}


# --------------------------------------------------------------------------- #
#  Real-basin MCMC health + Bayes-factor sanity                               #
# --------------------------------------------------------------------------- #
def real_fit_mcmc_health(*, draws, tune, chains, seed=0):
    import arviz as az
    from mls_emergence.inference import sample_fst, fst_summary, gini_simpson_summary
    from mls_emergence.signatures.variance import cultural_fst

    a43 = importlib.import_module("43_bayesian_fst")
    a07 = importlib.import_module("07_refined_empirical")
    inp = a07.prepare_inputs()
    gc, sizes = a43.basin_group_counts(inp)

    idata = sample_fst(gc, draws=draws, tune=tune, chains=chains, random_seed=seed)
    s = fst_summary(idata)
    g = gini_simpson_summary(idata, gc, seed=0)

    def _flat(diag):
        ds = diag.dataset if hasattr(diag, "dataset") else diag
        return np.concatenate([np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()])
    return {
        "n_groups": int(gc.shape[0]),
        "n_types": int(gc.shape[1]),
        "fst_median": s["fst_median"],
        "fst_hdi95": s["fst_hdi95"],
        "gst_median": g["gst_median"],
        "gst_hdi95": g["gst_hdi95"],
        "plugin_fst": float(cultural_fst(gc.astype(float))),
        "max_rhat": float(np.max(_flat(az.rhat(idata)))),
        "min_ess": float(np.min(_flat(az.ess(idata)))),
        "n_divergences": int(idata.sample_stats["diverging"].sum()),
    }


def bayes_factor_sanity(*, draws, chains, seed=0):
    """BF on panmixia-generated and strongly-structured synthetic data."""
    from mls_emergence.inference import bayes_factor_structure
    rng = np.random.default_rng(seed)
    # panmixia: two clusters from one shared distribution
    pan, _, _ = simulate_bn(1e-3, P_ANC, SIZES, rng)
    # structure: strong differentiation
    struct, _, _ = simulate_bn(0.30, P_ANC, SIZES, rng)
    return {
        "panmixia": bayes_factor_structure(pan, draws=draws, chains=chains, seed=1),
        "structure": bayes_factor_structure(struct, draws=draws, chains=chains, seed=2),
    }


def _make_figures(cov_rows, ranks, n_bins):
    fig, axes = plt.subplots(1, 2, figsize=(7, 3))
    ax = axes[0]
    f_true = [r["f_true"] for r in cov_rows]
    ax.plot(f_true, [r["coverage_f"] for r in cov_rows], marker="o", color=OI_BLUE,
            lw=1.2, label="BN F parameter")
    ax.plot(f_true, [r["coverage_g"] for r in cov_rows], marker="s", color=OI_GREEN,
            lw=1.2, label="Gini-Simpson readout")
    ax.axhline(0.95, color=OI_VERMIL, lw=1.0, ls="--", label="nominal 95%")
    ax.set_xlabel("true cultural $F_{ST}$")
    ax.set_ylabel("95% CI coverage")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=7)

    ax2 = axes[1]
    ax2.hist(ranks, bins=n_bins, range=(0, 1), color=OI_BLUE, alpha=0.75)
    ax2.axhline(len(ranks) / n_bins, color=OI_VERMIL, lw=1.0, ls="--",
                label="uniform expectation")
    ax2.set_xlabel("SBC rank of true $F_{ST}$")
    ax2.set_ylabel("count")
    ax2.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "bayesian_fst_coverage")
    plt.close(fig)


def main(fast=False, serial=False):
    cfg = FAST if fast else FULL
    n_workers = 1 if serial else None

    cov_rows = coverage(cfg["f_grid"], cfg["reps"], P_ANC, SIZES,
                        draws=cfg["draws"], tune=cfg["tune"], chains=cfg["chains"],
                        seed=100, n_workers=n_workers)
    small_bias = coverage(cfg["f_grid"], cfg["reps"], P_ANC, cfg["small_sizes"],
                          draws=cfg["draws"], tune=cfg["tune"], chains=cfg["chains"],
                          seed=200, n_workers=n_workers)
    ranks = sbc_ranks(cfg["n_sbc"], len(P_ANC), SIZES,
                      draws=cfg["draws"], tune=cfg["tune"], chains=cfg["chains"],
                      seed=300, n_workers=n_workers)
    unif = sbc_uniformity(ranks, cfg["sbc_bins"])
    health = real_fit_mcmc_health(draws=300 if fast else 2000,
                                  tune=500 if fast else 2000,
                                  chains=2 if fast else 4)
    bf = bayes_factor_sanity(draws=500 if fast else 1500, chains=2)

    abs_plugin = np.mean([r["abs_plugin_bias"] for r in small_bias])
    abs_post = np.mean([r["abs_post_bias"] for r in small_bias])

    L = ["# Bayesian F_ST validation battery (Balding-Nichols model)", "",
         f"Config: {'FAST' if fast else 'FULL'} (f_grid={list(cfg['f_grid'])}, "
         f"reps={cfg['reps']}, draws={cfg['draws']} x tune={cfg['tune']} x "
         f"chains={cfg['chains']}; large sizes={SIZES.tolist()}, small sizes="
         f"{list(cfg['small_sizes'])}; SBC n={cfg['n_sbc']}).", "",
         "## 1-2. Coverage (large sizes)", "",
         "Coverage of the BN F parameter is checked against the generating F; "
         "coverage of the Gini-Simpson readout is checked against the REALIZED "
         "Gini-Simpson F_ST of each simulated dataset.", "",
         "| true F_ST | mean realized G-S F_ST | reps | coverage (BN F) | coverage (G-S readout) |",
         "|---|---|---|---|---|"]
    for r in cov_rows:
        L.append(f"| {r['f_true']:.3f} | {r['mean_realized']:.4f} | {r['reps']} | "
                 f"{r['coverage_f']:.3f} | {r['coverage_g']:.3f} |")
    L += ["", "## 3. Small-sample bias vs the realized Gini-Simpson F_ST "
          f"(sizes = {list(cfg['small_sizes'])})", "",
          "| true F_ST | mean realized | reps | plug-in bias | posterior bias |",
          "|---|---|---|---|---|"]
    for r in small_bias:
        L.append(f"| {r['f_true']:.3f} | {r['mean_realized']:.4f} | {r['reps']} | "
                 f"{r['mean_plugin_bias']:+.4f} | {r['mean_post_bias']:+.4f} |")
    L += ["", f"Mean absolute bias at small sizes: plug-in = {abs_plugin:.4f}, "
          f"posterior median = {abs_post:.4f} "
          f"({'posterior closer to truth' if abs_post < abs_plugin else 'comparable / plug-in closer'}).",
          "", "## 4. Simulation-based calibration (BN F parameter)", "",
          f"- {cfg['n_sbc']} datasets drawn from the prior; rank of true F within "
          f"its posterior binned into {cfg['sbc_bins']} bins.",
          f"- Chi-square uniformity: chi2 = {unif['chi2']:.2f}, p = "
          f"{unif['chi2_p']:.3f} (want p > 0.05, i.e. not distinguishable from uniform).",
          f"- KS uniformity: D = {unif['ks_stat']:.3f}, p = {unif['ks_p']:.3f}.", "",
          "## 5. MCMC health on the real St. Francis basin fit", "",
          f"- {health['n_groups']} spatial clusters, {health['n_types']} decorated "
          f"types; default target_accept (no funnel-taming needed).",
          f"- max R-hat = {health['max_rhat']:.4f} (want < 1.01); min ESS = "
          f"{health['min_ess']:.0f}; divergences = {health['n_divergences']}.",
          f"- BN F_ST median = {health['fst_median']:.4f}, 95% CI "
          f"[{health['fst_hdi95'][0]:.4f}, {health['fst_hdi95'][1]:.4f}].",
          f"- Gini-Simpson readout median = {health['gst_median']:.4f}, 95% CI "
          f"[{health['gst_hdi95'][0]:.4f}, {health['gst_hdi95'][1]:.4f}]; plug-in = "
          f"{health['plugin_fst']:.4f} "
          f"({'inside' if health['gst_hdi95'][0] <= health['plugin_fst'] <= health['gst_hdi95'][1] else 'outside'} the interval).",
          "", "## 6. Bayes-factor sanity", "",
          f"- Panmixia-generated data: 2 ln BF10 = {bf['panmixia']['two_ln_bf']:.2f} "
          f"({bf['panmixia']['evidence']}).",
          f"- Strongly-structured data (F=0.30): 2 ln BF10 = "
          f"{bf['structure']['two_ln_bf']:.2f} ({bf['structure']['evidence']}).", ""]

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    _make_figures(cov_rows, ranks, cfg["sbc_bins"])

    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv, serial="--serial" in sys.argv)
