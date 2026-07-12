"""40_hierarchical_convergence.py - Bayesian convergence posterior.

Builds the three-signature CA panel and per-cell bootstrap SEs plus the
seriation fragmentation slope (all reusing analysis 07's prepare_inputs), fits
the hierarchical convergence model, and reports P(all four emergence slopes > 0).
Analysis 07 is unchanged. Usage:
    .venv/bin/python analyses/40_hierarchical_convergence.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

OUT_MD = ROOT / "output" / "hierarchical_convergence.md"
OUT_NC = ROOT / "output" / "hierarchical_convergence.nc"
PANEL_SIGS = ["neutral_departure", "fst", "spatial_boundary"]

FULL = dict(draws=2000, tune=2000, chains=4, n_boot=800)
FAST = dict(draws=300, tune=300, chains=2, n_boot=80)


def _zstd(x):
    x = np.asarray(x, float)
    m = np.nanmean(x)
    sd = np.nanstd(x)
    sd = sd if sd > 0 else 1.0
    return (x - m) / sd, float(m), float(sd)


def build_panel_and_ses(inp, n_bins=6, n_boot=800, seed=0):
    """Return z-standardized panel (3, n_bins), z-scale per-cell SE, std axis t.

    Point panel from a07.panel_for_bins(inp, n_bins). Per-cell SE from n_boot
    assemblage-bootstrap rebuilds of the panel (SD of each cell). Each signature
    row is z-standardized across slices; SEs divided by that row's SD; the axis
    (bin index or bin mean CA coordinate) standardized.
    """
    a07 = importlib.import_module("07_refined_empirical")
    panel = a07.panel_for_bins(inp, n_bins)          # DataFrame, cols=PANEL_SIGS
    y_raw = np.array([panel[s].to_numpy(float) for s in PANEL_SIGS])   # (3, n_bins)
    axis = panel.index.to_numpy(float)
    t, _, _ = _zstd(axis)

    # bootstrap per-cell SE: resample assemblages, rebuild panel, collect cells
    rng = np.random.default_rng(seed)
    boot = []  # list of (3, n_bins) arrays
    for _ in range(n_boot):
        p = a07.bootstrap_panel_once(inp, n_bins, rng)   # 07's exact bootstrap
        if p is None or list(p.index) != list(panel.index):
            continue
        boot.append(np.array([p[s].to_numpy(float) for s in PANEL_SIGS]))
    boot = np.array(boot)                              # (B, 3, n_bins)
    se_raw = np.nanstd(boot, axis=0)                   # (3, n_bins)

    # z-standardize each signature row; scale SE by the same row SD
    y = np.empty_like(y_raw)
    se = np.empty_like(se_raw)
    for i in range(len(PANEL_SIGS)):
        y[i], _, sd = _zstd(y_raw[i])
        se[i] = se_raw[i] / sd
    return y, se, t, PANEL_SIGS


def seriation_slope_and_se(inp, n_boot=800, seed=1):
    """Emergence-oriented seriation slope (-slope_frag), standardized, + boot SD."""
    a07 = importlib.import_module("07_refined_empirical")
    zc, _, _ = _zstd(inp.ca_vals)
    zn, _, _ = _zstd(inp.nmemb)
    b_ser = -a07.ols_slope(zc, zn)                    # emergence-oriented, standardized
    rng = np.random.default_rng(seed)
    boot = []
    n = len(inp.idx)
    for _ in range(n_boot):
        s = rng.integers(0, n, n)
        zcs, _, _ = _zstd(inp.ca_vals[s])
        zns, _, _ = _zstd(inp.nmemb[s])
        val = -a07.ols_slope(zcs, zns)
        if np.isfinite(val):
            boot.append(val)
    return float(b_ser), float(np.std(boot))


def main(fast=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figstyle import OI_BLUE, OI_VERMIL, save

    from mls_emergence.inference import sample_convergence, convergence_summary

    a07 = importlib.import_module("07_refined_empirical")

    cfg = FAST if fast else FULL
    inp = a07.prepare_inputs()
    y, se, t, labels = build_panel_and_ses(inp, n_bins=6, n_boot=cfg["n_boot"], seed=0)
    b_ser_obs, se_ser = seriation_slope_and_se(inp, n_boot=cfg["n_boot"], seed=1)

    idata = sample_convergence(y, se, t, b_ser_obs, se_ser,
                               draws=cfg["draws"], tune=cfg["tune"],
                               chains=cfg["chains"], random_seed=0,
                               target_accept=0.99)
    s = convergence_summary(idata)

    import arviz as az

    def _diag_extremum(diag_result, reducer):
        # arviz >=1.x returns a DataTree (not a plain Dataset) from az.rhat/az.ess;
        # pull the flat array of per-variable diagnostic values via .dataset.
        ds = diag_result.dataset if hasattr(diag_result, "dataset") else diag_result
        vals = np.concatenate(
            [np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()]
        )
        return float(reducer(vals))

    rhat = _diag_extremum(az.rhat(idata), np.max)
    ess = _diag_extremum(az.ess(idata), np.min)
    ndiv = int(idata.sample_stats["diverging"].sum())

    sig_names = labels + ["seriation"]
    L = ["# Hierarchical Bayesian convergence posterior", "",
         f"Four emergence slopes (3 CA-panel signatures + seriation fragmentation), "
         f"partial pooling, measurement-error likelihood. {'FAST' if fast else 'full'} "
         f"config: {cfg['draws']} draws x {cfg['chains']} chains.", "",
         "## Convergence", "",
         f"- **P(all four slopes > 0) = {s['p_convergence']:.3f}**.",
         f"- Convergence-score slope mean {s['cscore_slope_mean']:+.3f}, "
         f"P(>0) = {s['cscore_slope_p_pos']:.3f}.",
         f"- Shared emergence slope mu: mean {s['mu_mean']:+.3f}, "
         f"P(>0) = {s['mu_p_pos']:.3f}.", "",
         "## Per-signature slopes (emergence-oriented, standardized)", "",
         "| signature | posterior mean | 95% interval | P(>0) |",
         "|---|---|---|---|"]
    for i, nm in enumerate(sig_names):
        L.append(f"| {nm} | {s['b_mean'][i]:+.3f} | "
                 f"[{s['b_hdi95'][i][0]:+.3f}, {s['b_hdi95'][i][1]:+.3f}] | "
                 f"{s['b_p_pos'][i]:.3f} |")
    L += ["", "## MCMC diagnostics", "",
          f"- max R-hat = {rhat:.4f} (want < 1.01); min ESS = {ess:.0f}; "
          f"divergences = {ndiv}.", ""]

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    idata.to_netcdf(OUT_NC)

    # forest plot of the four slope posteriors
    fig, ax = plt.subplots(figsize=(3.5, 3))
    means = s["b_mean"]
    for i, nm in enumerate(sig_names):
        lo, hi = s["b_hdi95"][i]
        ax.plot([lo, hi], [i, i], color=OI_BLUE)
        ax.plot(means[i], i, "o", color=OI_VERMIL)
    ax.axvline(0.0, color="0.5", lw=0.8)
    ax.set_yticks(range(len(sig_names)))
    ax.set_yticklabels(sig_names)
    ax.set_xlabel("emergence slope (standardized)")
    save(fig, "hierarchical_convergence_slopes")
    plt.close(fig)

    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
