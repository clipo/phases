"""41_hierarchical_convergence_validation.py - validation battery for the
hierarchical Bayesian convergence model (analysis 40 / mls_emergence.inference
.convergence_model).

Four checks:
1. p_convergence calibration: synthetic panels with known true slope signs
   should recover a high P(all four > 0) when all four are truly positive and
   a low one when all four are truly negative (and something intermediate for
   a mixed case).
2. Prior sensitivity: refit the REAL panel under narrower/wider mu_sd / tau_sd
   / sigma_sd priors; the convergence verdict should not flip.
3. Posterior-predictive check on the real fit: do observed panel cells fall
   within the model's predictive interval?
4. Cross-read vs analysis 07: recompute 07's three panel-signature OLS slopes
   directly (a07.panel_for_bins + a07.ols_slope, not by parsing the report
   text) and confirm the model's per-signature slope posteriors carry the
   same sign / story.

Usage:
    .venv/bin/python analyses/41_hierarchical_convergence_validation.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

OUT_MD = ROOT / "output" / "hierarchical_convergence_validation.md"

FULL = dict(draws=2000, tune=2000, chains=4, n_boot=800,
            calib_draws=1000, calib_tune=1000, calib_chains=4,
            prior_draws=1000, prior_tune=1000, prior_chains=4)
FAST = dict(draws=300, tune=300, chains=2, n_boot=80,
            calib_draws=300, calib_tune=300, calib_chains=2,
            prior_draws=300, prior_tune=300, prior_chains=2)

FULL_CASES = [
    {"label": "all_up", "slope_panel": 1.2, "slope_ser": 1.2},
    {"label": "mixed", "slope_panel": 1.2, "slope_ser": -1.2},
    {"label": "all_down", "slope_panel": -1.2, "slope_ser": -1.2},
]
FAST_CASES = [
    {"label": "all_up", "slope_panel": 1.2, "slope_ser": 1.2},
    {"label": "all_down", "slope_panel": -1.2, "slope_ser": -1.2},
]

PRIOR_GRIDS = [
    {"mu_sd": 0.5, "tau_sd": 0.5, "sigma_sd": 0.25},   # narrower
    {"mu_sd": 1.0, "tau_sd": 1.0, "sigma_sd": 0.5},    # baseline (model defaults)
    {"mu_sd": 2.0, "tau_sd": 2.0, "sigma_sd": 1.0},    # wider
]


def _synthetic_panel(slope_panel, slope_ser, T=6, n_panel=3, noise=0.15, seed=0):
    """Panel z-values with a common linear slope per signature plus noise.

    Same construction as the Task-3 unit test's ``_synthetic_panel``
    (tests/inference/test_convergence_model.py): a shared true slope per
    signature drives all n_panel rows, plus independent Gaussian noise, with a
    matching fixed per-cell SE and a separately supplied seriation slope/SE.
    """
    rng = np.random.default_rng(seed)
    t = (np.arange(T) - (T - 1) / 2)
    t = t / t.std()
    y = np.array([slope_panel * t + rng.normal(0, noise, T) for _ in range(n_panel)])
    se = np.full((n_panel, T), noise)
    b_ser_obs = slope_ser
    se_ser = 0.15
    return y, se, t, b_ser_obs, se_ser


def _mcmc_health(idata):
    """max R-hat, min ESS, divergence count. arviz >=1.x returns a DataTree
    (not a plain Dataset) from az.rhat/az.ess; pull the flat array of
    per-variable diagnostic values via .dataset (mirrors analysis 40)."""
    import arviz as az

    def _extremum(diag_result, reducer):
        ds = diag_result.dataset if hasattr(diag_result, "dataset") else diag_result
        vals = np.concatenate(
            [np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()]
        )
        return float(reducer(vals))

    rhat = _extremum(az.rhat(idata), np.max)
    ess = _extremum(az.ess(idata), np.min)
    ndiv = int(idata.sample_stats["diverging"].sum())
    return rhat, ess, ndiv


def pconv_calibration(cases, *, T=6, draws=1000, tune=1000, chains=4, seed=0):
    """For each case (dict with label, slope_panel, slope_ser), synthesize a
    panel with that panel slope + seriation slope, fit the hierarchical
    convergence model, and record P(all four slopes > 0). Validates that
    p_convergence is high when all four true slopes are positive and low when
    they are not."""
    from mls_emergence.inference import sample_convergence, convergence_summary

    out = []
    for i, case in enumerate(cases):
        y, se, t, bso, ses = _synthetic_panel(
            case["slope_panel"], case["slope_ser"], T=T, seed=seed + i
        )
        idata = sample_convergence(y, se, t, bso, ses, draws=draws, tune=tune,
                                   chains=chains, target_accept=0.99,
                                   random_seed=seed + i)
        s = convergence_summary(idata)
        ndiv = int(idata.sample_stats["diverging"].sum())
        true_all_positive = bool(case["slope_panel"] > 0 and case["slope_ser"] > 0)
        out.append({
            "label": case["label"],
            "true_all_positive": true_all_positive,
            "p_convergence": s["p_convergence"],
            "ndiv": ndiv,
        })
    return out


def prior_sensitivity(y, se, t, b_ser_obs, se_ser, grids, *, draws=1000, tune=1000,
                      chains=4, seed=0):
    """Refit the REAL panel under each prior setting in ``grids`` (dicts of
    mu_sd / tau_sd / sigma_sd overrides); return p_convergence per setting so
    the verdict's stability can be checked."""
    from mls_emergence.inference import sample_convergence, convergence_summary

    out = []
    for i, grid in enumerate(grids):
        idata = sample_convergence(y, se, t, b_ser_obs, se_ser, draws=draws,
                                   tune=tune, chains=chains, random_seed=seed + i,
                                   target_accept=0.99, **grid)
        s = convergence_summary(idata)
        out.append({**grid, "p_convergence": s["p_convergence"]})
    return out


def posterior_predictive_summary(model, idata, y, seed=0):
    """Sample the posterior predictive for y_obs and report the fraction of
    observed panel cells falling within the 95% predictive interval."""
    import pymc as pm

    with model:
        ppc = pm.sample_posterior_predictive(idata, var_names=["y_obs"],
                                              random_seed=seed, progressbar=False)
    pred = ppc.posterior_predictive["y_obs"].stack(s=("chain", "draw")).values
    # pred: (n_panel, T, S)
    lo = np.percentile(pred, 2.5, axis=-1)
    hi = np.percentile(pred, 97.5, axis=-1)
    within = (y >= lo) & (y <= hi)
    return {
        "n_cells": int(within.size),
        "n_within": int(within.sum()),
        "frac_within": float(within.mean()),
    }


def cross_read_vs_analysis07(inp, a07, s_real, labels, n_bins=6):
    """Recompute analysis 07's three panel-signature OLS slopes directly
    (panel_for_bins + ols_slope, not by parsing output/empirical_refined.md)
    and compare their sign to the model's per-signature slope posteriors.
    Also confirms the seriation sign is self-consistent (it is a direct model
    input, so this is a sanity check, not an independent test)."""
    panel = a07.panel_for_bins(inp, n_bins)
    axis = panel.index.to_numpy(float)
    rows = []
    for i, sig in enumerate(labels):
        slope_07 = a07.ols_slope(axis, panel[sig].to_numpy(float))
        b_model = s_real["b_mean"][i]
        rows.append({
            "signature": sig,
            "slope_07": float(slope_07),
            "b_model": float(b_model),
            "same_sign": bool(np.sign(slope_07) == np.sign(b_model)),
        })
    return rows


def _bar_figure(calib_results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figstyle import OI_BLUE, OI_VERMIL, save

    labels = [r["label"] for r in calib_results]
    vals = [r["p_convergence"] for r in calib_results]
    colors = [OI_BLUE if r["true_all_positive"] else OI_VERMIL for r in calib_results]

    fig, ax = plt.subplots(figsize=(3.5, 3))
    ax.bar(labels, vals, color=colors)
    ax.axhline(0.5, color="0.5", lw=0.8, ls="--")
    ax.set_ylabel("P(all four slopes > 0)")
    ax.set_ylim(0, 1)
    save(fig, "hierarchical_convergence_calibration")
    plt.close(fig)


def main(fast=False):
    import matplotlib
    matplotlib.use("Agg")

    from mls_emergence.inference import build_convergence_model, convergence_summary

    a07 = importlib.import_module("07_refined_empirical")
    a40 = importlib.import_module("40_hierarchical_convergence")

    cfg = FULL if not fast else FAST
    cases = FULL_CASES if not fast else FAST_CASES

    L = ["# Hierarchical convergence model: validation battery", "",
         f"Config: {'FAST' if fast else 'full'} — real-panel fit "
         f"{cfg['draws']} draws x {cfg['chains']} chains "
         f"({cfg['n_boot']} bootstrap reps for panel/seriation SEs); "
         f"calibration fits {cfg['calib_draws']} draws x {cfg['calib_chains']} "
         f"chains over {len(cases)} synthetic cases "
         f"({'reduced from 3 to 2 cases under --fast' if fast else '3 cases'}); "
         f"prior-sensitivity refits {cfg['prior_draws']} draws x "
         f"{cfg['prior_chains']} chains over {len(PRIOR_GRIDS)} prior settings.",
         ""]

    # ------------------------------------------------------------------
    # 1. p_convergence calibration on synthetic panels
    # ------------------------------------------------------------------
    calib = pconv_calibration(cases, T=6, draws=cfg["calib_draws"],
                              tune=cfg["calib_tune"], chains=cfg["calib_chains"],
                              seed=0)
    L += ["## 1. p_convergence calibration (synthetic panels)", "",
          "| case | true all-positive | P(all four > 0) | divergences |",
          "|---|---|---|---|"]
    for r in calib:
        L.append(f"| {r['label']} | {r['true_all_positive']} | "
                 f"{r['p_convergence']:.3f} | {r['ndiv']} |")
    L.append("")
    _bar_figure(calib)

    # ------------------------------------------------------------------
    # Real panel: build once, sample once, reuse for prior sensitivity's
    # baseline row, the PPC, and the cross-read.
    # ------------------------------------------------------------------
    inp = a07.prepare_inputs()
    y, se, t, labels = a40.build_panel_and_ses(inp, n_bins=6, n_boot=cfg["n_boot"],
                                               seed=0)
    b_ser_obs, se_ser = a40.seriation_slope_and_se(inp, n_boot=cfg["n_boot"], seed=1)

    model = build_convergence_model(y, se, t, b_ser_obs, se_ser)
    import pymc as pm
    with model:
        idata_real = pm.sample(draws=cfg["draws"], tune=cfg["tune"],
                               chains=cfg["chains"], target_accept=0.99,
                               random_seed=0, progressbar=False)
    s_real = convergence_summary(idata_real)
    rhat, ess, ndiv = _mcmc_health(idata_real)

    L += ["## Real-panel fit (baseline)", "",
          f"- P(all four slopes > 0) = {s_real['p_convergence']:.3f}.",
          f"- MCMC health: max R-hat = {rhat:.4f}, min ESS = {ess:.0f}, "
          f"divergences = {ndiv}.", ""]

    # ------------------------------------------------------------------
    # 2. Prior sensitivity on the real panel
    # ------------------------------------------------------------------
    prior = prior_sensitivity(y, se, t, b_ser_obs, se_ser, PRIOR_GRIDS,
                              draws=cfg["prior_draws"], tune=cfg["prior_tune"],
                              chains=cfg["prior_chains"], seed=10)
    pvals = [r["p_convergence"] for r in prior]
    stable = (max(pvals) - min(pvals)) < 0.3 and \
        all((p > 0.5) == (pvals[0] > 0.5) for p in pvals)
    L += ["## 2. Prior sensitivity (real panel)", "",
          "| mu_sd | tau_sd | sigma_sd | P(all four > 0) |",
          "|---|---|---|---|"]
    for r in prior:
        L.append(f"| {r['mu_sd']} | {r['tau_sd']} | {r['sigma_sd']} | "
                 f"{r['p_convergence']:.3f} |")
    L += ["", f"- Verdict {'STABLE' if stable else 'UNSTABLE'} across priors "
          f"(range {min(pvals):.3f}-{max(pvals):.3f}).", ""]

    # ------------------------------------------------------------------
    # 3. Posterior-predictive check
    # ------------------------------------------------------------------
    ppc = posterior_predictive_summary(model, idata_real, y, seed=0)
    L += ["## 3. Posterior-predictive check (real panel)", "",
          f"- {ppc['n_within']}/{ppc['n_cells']} observed panel cells "
          f"({ppc['frac_within']:.1%}) fall within the model's 95% posterior-"
          f"predictive interval.", ""]

    # ------------------------------------------------------------------
    # 4. Cross-read vs analysis 07
    # ------------------------------------------------------------------
    cross = cross_read_vs_analysis07(inp, a07, s_real, labels, n_bins=6)
    L += ["## 4. Cross-read vs analysis 07 (recomputed, not parsed from "
          "output/empirical_refined.md)", "",
          "| signature | 07 recomputed OLS slope | model posterior mean slope | "
          "same sign |",
          "|---|---|---|---|"]
    all_match = True
    for r in cross:
        L.append(f"| {r['signature']} | {r['slope_07']:+.5f} | "
                 f"{r['b_model']:+.3f} | {r['same_sign']} |")
        all_match = all_match and r["same_sign"]
    ser_sign_match = np.sign(b_ser_obs) == np.sign(s_real["b_mean"][-1])
    L.append(f"| seriation | {b_ser_obs:+.5f} | {s_real['b_mean'][-1]:+.3f} | "
             f"{bool(ser_sign_match)} |")
    L += ["", f"- All panel-signature signs match between the recomputed "
          f"07 slopes and the model posterior: {all_match}.", ""]

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
