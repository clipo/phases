"""44_bayesian_fst_validation.py - Validation battery for the Bayesian F_ST fit.

Four checks on the hierarchical Dirichlet-multinomial credible interval
(mls_emergence.inference.sample_fst / fst_summary, analysis 43):

1. Coverage: at known true F_ST values, does the nominal 95% credible interval
   contain the truth ~95% of the time over repeated simulated datasets?
2. Small-sample bias: at small group sizes, does the posterior mean track truth
   better than the frequentist plug-in Gini-Simpson estimator?
3. MCMC health on the real St. Francis basin fit (analysis 43's
   basin_group_counts): max R-hat, min ESS, divergence count. This script only
   reports health; it does not attempt to fix divergences (that is analysis 43's
   / a later task's job).
4. Frequentist cross-read: posterior 95% interval vs variance.cultural_fst on
   the same basin data.

Usage: .venv/bin/python analyses/44_bayesian_fst_validation.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
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

# Full battery config vs a fast smoke config. --fast STATES this in the report;
# it does not silently down-scale without saying so.
FULL = dict(
    f_grid=(0.02, 0.05, 0.10, 0.20, 0.35, 0.50),
    reps=40,
    draws=1000, tune=1000, chains=4,
    small_sizes=np.array([15, 15, 15]),
)
FAST = dict(
    f_grid=(0.05, 0.30),
    reps=4,
    draws=200, tune=300, chains=2,
    small_sizes=np.array([15, 15, 15]),
)

P_ANC = np.array([0.4, 0.35, 0.25])
SIZES = np.array([200, 200, 200])


def simulate_groups(f_true, p_anc, sizes, rng):
    """Draw (G, K) counts from the F-model at true F_ST f_true.

    conc = (1 - f_true) / f_true; for each group draw p_g ~ Dirichlet(conc *
    p_anc), then counts_g ~ Multinomial(sizes[g], p_g). f_true is guarded away
    from 0/1 so the concentration stays finite and positive.
    """
    from mls_emergence.inference import fst_from_frequencies
    f_true = float(np.clip(f_true, 1e-3, 1 - 1e-3))
    p_anc = np.asarray(p_anc, float)
    sizes = np.asarray(sizes, int)
    conc = (1.0 - f_true) / f_true
    G = len(sizes)
    K = len(p_anc)
    counts = np.zeros((G, K), dtype=int)
    p = np.zeros((G, K))
    for g in range(G):
        p[g] = rng.dirichlet(conc * p_anc)
        counts[g] = rng.multinomial(sizes[g], p[g])
    # The model estimates the Gini-Simpson F_ST of the drawn groups, NOT the
    # F-model concentration parameter f_true (a different quantity). Return the
    # realized Gini-Simpson F_ST of the actual generating frequencies so coverage
    # and bias are assessed against what the model actually targets.
    g_true = fst_from_frequencies(p, sizes)
    return counts, g_true


def coverage(f_grid, reps, p_anc, sizes, *, draws, tune, chains, seed):
    """Simulation-based calibration of the 95% credible interval on F_ST.

    For each true F_ST in f_grid, simulate `reps` datasets, fit each with
    sample_fst, and record whether the 95% HDI covers the truth, plus the
    plug-in and posterior-mean values (bias comparison). Returns one dict per
    grid point with the coverage fraction and mean plug-in / posterior mean.
    """
    from mls_emergence.inference import sample_fst, fst_summary

    rng = np.random.default_rng(seed)
    results = []
    for f_true in f_grid:
        covered = []
        realized = []
        plugin_bias = []
        post_bias = []
        for r in range(reps):
            gc, g_true = simulate_groups(f_true, p_anc, sizes, rng)
            idata = sample_fst(gc, draws=draws, tune=tune, chains=chains,
                                target_accept=0.99,
                                random_seed=int(rng.integers(0, 2**31 - 1)))
            s = fst_summary(idata, sizes=gc.sum(axis=1), group_counts=gc)
            lo, hi = s["fst_hdi95"]
            # target the REALIZED Gini-Simpson F_ST (what the model estimates)
            covered.append(lo <= g_true <= hi)
            realized.append(g_true)
            plugin_bias.append(s["plugin_fst"] - g_true)
            post_bias.append(s["fst_mean"] - g_true)
        results.append({
            "f_true": float(f_true),
            "mean_realized": float(np.mean(realized)),
            "reps": int(reps),
            "coverage": float(np.mean(covered)),
            "mean_plugin_bias": float(np.mean(plugin_bias)),
            "mean_post_bias": float(np.mean(post_bias)),
            "abs_plugin_bias": float(np.mean(np.abs(plugin_bias))),
            "abs_post_bias": float(np.mean(np.abs(post_bias))),
        })
    return results


def _flat_diag(diag):
    """Flatten an arviz rhat/ess result across possible DataTree/Dataset shapes."""
    ds = diag.dataset if hasattr(diag, "dataset") else diag
    return np.concatenate([np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()])


def real_fit_mcmc_health(*, draws, tune, chains, seed=0):
    """Fit the real St. Francis basin data and report MCMC health + the
    frequentist cross-read. Reuses analysis 43's basin_group_counts."""
    import arviz as az
    from mls_emergence.inference import sample_fst, fst_summary
    from mls_emergence.signatures.variance import cultural_fst

    a43 = importlib.import_module("43_bayesian_fst")
    a07 = importlib.import_module("07_refined_empirical")

    inp = a07.prepare_inputs()
    gc, sizes = a43.basin_group_counts(inp)

    idata = sample_fst(gc, draws=draws, tune=tune, chains=chains,
                       target_accept=0.99, random_seed=seed)
    s = fst_summary(idata, sizes=sizes, group_counts=gc)

    rhat = float(np.max(_flat_diag(az.rhat(idata))))
    ess = float(np.min(_flat_diag(az.ess(idata))))
    ndiv = int(idata.sample_stats["diverging"].sum())

    plugin_direct = float(cultural_fst(gc.astype(float)))

    return {
        "n_groups": int(gc.shape[0]),
        "n_types": int(gc.shape[1]),
        "fst_mean": s["fst_mean"],
        "fst_hdi95": s["fst_hdi95"],
        "plugin_fst": s["plugin_fst"],
        "plugin_fst_direct": plugin_direct,
        "max_rhat": rhat,
        "min_ess": ess,
        "n_divergences": ndiv,
    }


def _make_coverage_figure(cov_rows, small_bias):
    fig, axes = plt.subplots(1, 2, figsize=(7, 3))

    ax = axes[0]
    f_true = [r["f_true"] for r in cov_rows]
    cov = [r["coverage"] for r in cov_rows]
    ax.plot(f_true, cov, marker="o", color=OI_BLUE, lw=1.2)
    ax.axhline(0.95, color=OI_VERMIL, lw=1.0, ls="--", label="nominal 95%")
    ax.set_xlabel("true cultural $F_{ST}$")
    ax.set_ylabel("95% CI coverage")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=7)

    ax2 = axes[1]
    xs = np.arange(len(small_bias))
    width = 0.35
    truths = [b["f_true"] for b in small_bias]
    plug_bias = [b["mean_plugin_bias"] for b in small_bias]
    post_bias = [b["mean_post_bias"] for b in small_bias]
    ax2.bar(xs - width / 2, plug_bias, width, label="plug-in bias", color=OI_GREEN)
    ax2.bar(xs + width / 2, post_bias, width, label="posterior-mean bias", color=OI_BLUE)
    ax2.axhline(0.0, color=OI_VERMIL, lw=1.0)
    ax2.set_xticks(xs)
    ax2.set_xticklabels([f"{t:.2f}" for t in truths])
    ax2.set_xlabel("F-model $F_{ST}$ (small sizes)")
    ax2.set_ylabel("bias vs realized $F_{ST}$")
    ax2.legend(frameon=False, fontsize=7)

    fig.tight_layout()
    save(fig, "bayesian_fst_coverage")
    plt.close(fig)


def main(fast=False):
    cfg = FAST if fast else FULL

    cov_rows = coverage(
        cfg["f_grid"], cfg["reps"], P_ANC, SIZES,
        draws=cfg["draws"], tune=cfg["tune"], chains=cfg["chains"], seed=0,
    )

    # Small-sample bias: reuse coverage() at small group sizes so we get both
    # the coverage fraction and the plug-in/posterior-mean bias comparison at
    # a size regime where the plug-in Gini-Simpson estimator is known to be
    # biased upward.
    small_bias = coverage(
        cfg["f_grid"], cfg["reps"], P_ANC, cfg["small_sizes"],
        draws=cfg["draws"], tune=cfg["tune"], chains=cfg["chains"], seed=1,
    )

    # The real-fit health check uses analysis 43's PRODUCTION config (larger
    # budget + target_accept=0.99) so it reports the health of the headline fit,
    # not the smaller per-replicate coverage budget.
    health = real_fit_mcmc_health(
        draws=300 if fast else 3000, tune=500 if fast else 3000,
        chains=2 if fast else 4,
    )

    L = ["# Bayesian F_ST validation battery", "",
         f"Config: {'FAST' if fast else 'FULL'} "
         f"(f_grid={list(cfg['f_grid'])}, reps={cfg['reps']}, "
         f"draws={cfg['draws']} x tune={cfg['tune']} x chains={cfg['chains']}; "
         f"large sizes={SIZES.tolist()}, small sizes={cfg['small_sizes'].tolist()}).",
         "", "## 1. Coverage of the 95% credible interval (large sizes)", "",
         "Coverage is assessed against the REALIZED Gini-Simpson F_ST of each "
         "simulated dataset (the quantity the model estimates), not the F-model "
         "concentration parameter used to generate it.", "",
         "| F-model F_ST | mean realized F_ST | reps | coverage |",
         "|---|---|---|---|"]
    for r in cov_rows:
        L.append(f"| {r['f_true']:.3f} | {r['mean_realized']:.4f} | {r['reps']} | "
                  f"{r['coverage']:.3f} |")
    L += ["", "Nominal target is 0.95 coverage at each grid point. With few "
          "reps (fast mode) coverage estimates are noisy; the full battery "
          "should be run for a publication-grade calibration check.", "",
          "## 2. Small-sample bias vs the realized F_ST (sizes = "
          f"{cfg['small_sizes'].tolist()})", "",
          "| F-model F_ST | mean realized F_ST | reps | coverage | "
          "plug-in bias | posterior bias |",
          "|---|---|---|---|---|---|"]
    for r in small_bias:
        L.append(f"| {r['f_true']:.3f} | {r['mean_realized']:.4f} | {r['reps']} | "
                  f"{r['coverage']:.3f} | "
                  f"{r['mean_plugin_bias']:+.4f} | "
                  f"{r['mean_post_bias']:+.4f} |")
    abs_plugin_bias = np.mean([r["abs_plugin_bias"] for r in small_bias])
    abs_post_bias = np.mean([r["abs_post_bias"] for r in small_bias])
    L += ["", f"Mean absolute bias at small sizes (vs the realized F_ST): plug-in "
          f"= {abs_plugin_bias:.4f}, posterior mean = {abs_post_bias:.4f} "
          f"({'posterior closer to truth' if abs_post_bias < abs_plugin_bias else 'comparable / plug-in closer'}).",
          "", "## 3. MCMC health on the real St. Francis basin fit", "",
          f"- {health['n_groups']} spatial clusters, {health['n_types']} decorated "
          f"types, draws={cfg['draws']} x tune={cfg['tune']} x chains={cfg['chains']}.",
          f"- max R-hat = {health['max_rhat']:.4f} (want < 1.01).",
          f"- min ESS = {health['min_ess']:.0f}.",
          f"- divergences = {health['n_divergences']}.",
          "", "Divergences, if any, are not addressed here: reparameterization "
          "or prior tuning is a separate follow-on task; this script only "
          "reports the diagnostic.", "",
          "## 4. Frequentist cross-read on the real basin data", "",
          f"- Posterior mean F_ST = {health['fst_mean']:.4f}, 95% credible "
          f"interval [{health['fst_hdi95'][0]:.4f}, {health['fst_hdi95'][1]:.4f}].",
          f"- Frequentist plug-in (variance.cultural_fst, direct) = "
          f"{health['plugin_fst_direct']:.4f}.",
          f"- Plug-in falls "
          f"{'inside' if health['fst_hdi95'][0] <= health['plugin_fst_direct'] <= health['fst_hdi95'][1] else 'outside'} "
          "the 95% credible interval.", ""]

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    _make_coverage_figure(cov_rows, small_bias)

    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
