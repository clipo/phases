"""50_perbin_bayesian_fst.py - per-bin Bayesian F_ST, and the trend as a posterior.

Converts two frequentist objects at once (docs/FREQUENTIST_INVENTORY.md, items
B and H), under the rule-18 ruling of 2026-08-31:

  B  The main text reports the size-controlled F_ST trajectory as a Spearman
     rank correlation with a RAREFACTION RESAMPLING interval (400 replicates,
     analyses/21_signal_recovery.py).
  H  analyses/40_hierarchical_convergence.py feeds its Bayesian trend model
     per-cell standard errors taken from 800 ASSEMBLAGE BOOTSTRAP rebuilds
     (`se_raw = np.nanstd(boot, axis=0)`), which enter the likelihood as if
     known.

Both are replaced here by the same thing: a Balding-Nichols posterior for each
seriation bin, whose spread is generative rather than resampled. The trend is
then estimated over those posteriors.

WHAT THIS IS AND IS NOT. This is a two-stage ("cut") analysis: each bin's F_ST
posterior is summarised and the summary is carried into the trend model. It is
NOT a single joint model of bins and trend. The honest consequence is that
uncertainty in the trend is propagated through a normal approximation to each
bin's posterior rather than exactly. That is a strict improvement on a bootstrap
standard deviation, which has no generative status at all, but it is not the
fully joint model and should not be described as one.

Rule 20: fitted under BOTH the uniform prior and Beta(1,3), with the movement
reported whether or not it is small.

Usage: .venv/bin/python analyses/50_perbin_bayesian_fst.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

OUT_MD = ROOT / "output" / "findings" / "perbin_bayesian_fst.md"
N_BINS = 6
# target_accept 0.95 rather than sample_fst's default 0.9: the full run on
# 2026-08-31 produced 7 divergent transitions in bin 4 under the uniform prior.
# Bin 4 pools 14,224 sherds, by far the most, so its posterior is the tightest
# and the default step size is too coarse for it. Rule 17: divergences are a
# validity signal and get fixed, not reported around.
FULL = dict(draws=2000, tune=2000, chains=4, target_accept=0.95)
FAST = dict(draws=400, tune=600, chains=2, target_accept=0.95)
# Beta(1,10) is PRIMARY from 2026-09-02; see analyses/57_fst_prior_predictive.py
# and the note in analyses/43. Uniform(0,1) is retained as the rule-20(b)
# sensitivity refit, deliberately leaning the other way.
PRIORS = [("Beta(1,10)", ("beta", 1.0, 10.0)), ("uniform", ("uniform",))]


def per_bin_group_counts(inp, n_bins=N_BINS):
    """(bin -> (n_clusters, K) counts) on the BASIN, clusters re-selected on the
    basin's own coordinates. Mirrors 43_bayesian_fst.basin_group_counts."""
    import pandas as pd
    import make_figures as mf
    a43 = importlib.import_module("43_bayesian_fst")

    members = mf._basin_members("curated")
    ids_all = [str(a) for a in inp.have_coords_ids]
    keep = [i for i, a in enumerate(ids_all) if a in members]
    ids = [ids_all[i] for i in keep]
    coords_b = np.asarray(inp.cc_c)[keep]
    labels = a43._basin_cluster_labels(coords_b)
    counts_b = inp.counts.loc[ids].to_numpy(float)

    ca = inp.ca.reindex(ids).to_numpy(float)
    bins = pd.qcut(pd.Series(ca, index=ids), q=n_bins, labels=False,
                   duplicates="drop").to_numpy()

    out = {}
    for b in sorted(set(bins[~np.isnan(bins)])):
        m = bins == b
        lab_b = labels[m]
        uniq = np.unique(lab_b)
        if len(uniq) < 2:            # F_ST undefined with one group present
            out[int(b)] = None
            continue
        gc = np.array([counts_b[m][lab_b == c].sum(axis=0) for c in uniq])
        gc = gc[gc.sum(1) > 0]
        out[int(b)] = gc if gc.shape[0] >= 2 else None
    return out, len(ids), len(np.unique(labels))


def main(fast=False):
    import arviz as az
    from mls_emergence.inference import sample_fst, gini_simpson_summary
    a07 = importlib.import_module("07_refined_empirical")
    cfg = FAST if fast else FULL

    inp = a07.prepare_inputs()
    per_bin, n_ids, k = per_bin_group_counts(inp)
    print(f"basin: {n_ids} assemblages, {k} spatial clusters, {N_BINS} bins")

    results = {}
    for pname, prior in PRIORS:
        rows = []
        for b, gc in sorted(per_bin.items()):
            if gc is None:
                rows.append(dict(bin=b, ok=False))
                print(f"  [{pname}] bin {b}: fewer than two clusters present, skipped")
                continue
            idata = sample_fst(gc, random_seed=0, f_prior=prior, **cfg)
            g = gini_simpson_summary(idata, gc, seed=0)
            s = np.asarray(g["gst_samples"], float)
            ndiv = int(idata.sample_stats["diverging"].sum())
            rhat = float(np.max(np.concatenate(
                [np.atleast_1d(v.values).ravel()
                 for v in az.rhat(idata).data_vars.values()])))
            ess = float(np.min(np.concatenate(
                [np.atleast_1d(v.values).ravel()
                 for v in az.ess(idata).data_vars.values()])))
            rows.append(dict(bin=b, ok=True, med=float(np.median(s)),
                             sd=float(np.std(s)),
                             lo=float(np.percentile(s, 2.5)),
                             hi=float(np.percentile(s, 97.5)),
                             plugin=g["plugin_fst"], rhat=rhat, ess=ess, ndiv=ndiv,
                             n=int(gc.sum())))
            print(f"  [{pname}] bin {b}: F_ST {rows[-1]['med']:.4f} "
                  f"[{rows[-1]['lo']:.4f}, {rows[-1]['hi']:.4f}]  "
                  f"n={rows[-1]['n']}  R-hat {rhat:.4f} div {ndiv}")
        results[pname] = rows

    # ---- the trend, as a posterior over per-bin posteriors ----------------- #
    import pymc as pm
    trend = {}
    for pname, rows in results.items():
        ok = [r for r in rows if r["ok"]]
        y = np.array([r["med"] for r in ok])
        se = np.array([r["sd"] for r in ok])
        x = np.array([r["bin"] for r in ok], float)
        x = (x - x.mean()) / (x.std() if x.std() > 0 else 1.0)
        with pm.Model():
            a = pm.Normal("a", 0.0, 0.5)
            bslope = pm.Normal("b", 0.0, 0.5)
            extra = pm.HalfNormal("extra", 0.05)
            pm.Normal("obs", mu=a + bslope * x,
                      sigma=pm.math.sqrt(se ** 2 + extra ** 2), observed=y)
            idt = pm.sample(draws=cfg["draws"], tune=cfg["tune"],
                            chains=cfg["chains"], random_seed=0,
                            target_accept=0.95, progressbar=False,
                            compute_convergence_checks=False)
        bs = np.asarray(idt.posterior["b"].values).ravel()
        trend[pname] = dict(
            med=float(np.median(bs)), lo=float(np.percentile(bs, 2.5)),
            hi=float(np.percentile(bs, 97.5)), p_pos=float((bs > 0).mean()),
            ndiv=int(idt.sample_stats["diverging"].sum()),
            rhat=float(np.max(np.concatenate(
                [np.atleast_1d(v.values).ravel()
                 for v in az.rhat(idt).data_vars.values()]))),
            ess=float(np.min(np.concatenate(
                [np.atleast_1d(v.values).ravel()
                 for v in az.ess(idt).data_vars.values()]))))
        print(f"  [{pname}] TREND slope {trend[pname]['med']:+.4f} "
              f"[{trend[pname]['lo']:+.4f}, {trend[pname]['hi']:+.4f}]  "
              f"P(>0) = {trend[pname]['p_pos']:.3f}")

    L = ["# Per-bin Bayesian F_ST, and the trend as a posterior", "",
         f"Produced by `analyses/50_perbin_bayesian_fst.py` "
         f"({'FAST' if fast else 'full'}). St. Francis basin, {n_ids} "
         f"assemblages, {k} spatial clusters, {N_BINS} seriation bins.", "",
         "Replaces two frequentist objects with one generative one "
         "(`docs/FREQUENTIST_INVENTORY.md` items B and H): the rarefaction "
         "resampling interval on the reported trajectory, and the 800-draw "
         "assemblage-bootstrap standard errors that "
         "`analyses/40_hierarchical_convergence.py` feeds into its likelihood "
         "as if known.", "",
         "## Per-bin Gini-Simpson F_ST posteriors", "",
         "| bin | prior | median | 95% CI | posterior SD | plug-in | sherds | R-hat | min ESS | div |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for pname, rows in results.items():
        for r in rows:
            if not r["ok"]:
                L.append(f"| {r['bin']} | {pname} | (fewer than two clusters) | - | - | - | - | - | - | - |")
                continue
            L.append(f"| {r['bin']} | {pname} | {r['med']:.4f} | "
                     f"[{r['lo']:.4f}, {r['hi']:.4f}] | {r['sd']:.4f} | "
                     f"{r['plugin']:.4f} | {r['n']} | {r['rhat']:.4f} | "
                     f"{r['ess']:.0f} | {r['ndiv']} |")
    L += ["", "## The trend, as a posterior slope", "",
          "Slope of per-bin F_ST against standardized bin position. This is the "
          "quantity the manuscript currently reports as a Spearman rank "
          "correlation with a resampling interval.", "",
          "| prior | slope median | 95% CI | P(slope > 0) | R-hat | min ESS | div |",
          "|---|---|---|---|---|---|---|"]
    for pname, tr in trend.items():
        L.append(f"| {pname} | {tr['med']:+.4f} | [{tr['lo']:+.4f}, "
                 f"{tr['hi']:+.4f}] | {tr['p_pos']:.3f} | {tr['rhat']:.4f} | "
                 f"{tr['ess']:.0f} | {tr['ndiv']} |")
    # Reference the priors BY POSITION in PRIORS, never by hard-coded name.
    # These lines used to read trend["uniform"], trend["Beta(1,3)"], which
    # silently stopped matching when the primary prior changed on 2026-09-02:
    # the script then crashed AFTER printing its results and BEFORE writing its
    # finding, so the console looked right while the committed file went stale
    # (rule 4, no silent parameters; rule 11, verify the write).
    primary_name, alt_name = PRIORS[0][0], PRIORS[1][0]
    u, bta = trend[primary_name], trend[alt_name]
    L += ["", "## Prior sensitivity (rule 20b)", "",
          f"The slope moves from {u['med']:+.4f} under the {primary_name} prior to "
          f"{bta['med']:+.4f} under {alt_name}, a shift of "
          f"{abs(u['med'] - bta['med']):.4f}. P(slope > 0) moves from "
          f"{u['p_pos']:.3f} to {bta['p_pos']:.3f}. Reported whether or not it "
          f"is small, as the rule requires.", "",
          "## Limitation, stated rather than buried", "",
          "This is a two-stage analysis. Each bin's posterior is summarised by "
          "its median and standard deviation, and those summaries are carried "
          "into the trend model, so trend uncertainty is propagated through a "
          "normal approximation to each bin's posterior rather than exactly. "
          "A single joint model of bins and trend would be better and is not "
          "what this is. The improvement over what it replaces is nonetheless "
          "categorical: a posterior standard deviation has generative status, "
          "and a bootstrap standard deviation does not.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"\nwrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
