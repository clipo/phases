"""29_concept_figure.py — the conceptual figure: what we are trying to do.

Illustrates the test graphically, from synthetic data matched to the record:
  Row 1: what the assemblages SHOULD look like if the Parkin phase were a bounded
         interaction group (sharp clusters), what spatially structured drift
         produces instead (a continuous blur), and what the OBSERVED assemblages
         look like (a blur, like drift).
  Row 2: why the determination cannot be made. Cultural F_ST recovered from
         synthetic assemblages as the strength of group closure rises from 0
         (drift) to 1 (strong groups), with the detection threshold s* and the
         observed value marked. Below s* a weak group and pure drift produce the
         same F_ST at this record's size and time-averaging, so the observed
         value (in the shaded region) is consistent with both.

Read-only on the manuscript. Writes figures/fig2_concept.png.

Usage: .venv/bin/python analyses/29_concept_figure.py
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
import pandas as pd  # noqa: E402
import make_figures as mf  # noqa: E402
demo = importlib.import_module("25_drift_vs_groups_demo")
sd = importlib.import_module("23_phases_vs_spatial_drift")
r21 = importlib.import_module("21_signal_recovery")
res17 = importlib.import_module("17_basin_results")

OUT = ROOT / "figures" / "fig2_concept.png"
C0 = "#0072B2"   # spatial cluster 0
C1 = "#D55E00"   # spatial cluster 1
LEN = 24.0       # calibrated interaction range (km)
MB = 0.02        # low between-node mixing


def fst_groups(counts, labels):
    pools = np.array([counts[labels == c].sum(0) for c in np.unique(labels)
                      if counts[labels == c].sum() > 0])
    return mf.cultural_fst(pools) if len(pools) >= 2 else np.nan


def main():
    counts_df, coords_df = mf._load_curated()
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    K = counts.shape[1]
    totals = counts.sum(1)
    cc = coords - coords.mean(0)
    glab = mf._kmeans_labels(cc, 2, seed=7)

    Wd = demo.drift_weights(coords, length_km=LEN)
    Wg = demo.group_weights(coords, glab, length_km=LEN, leak=0.03)
    Mg = demo.simulate(Wg, K=K, total_per_node=totals, m_between=MB, seed=11)
    Md = demo.simulate(Wd, K=K, total_per_node=totals, m_between=MB, seed=11)

    # ---- Row 2 recovery curve (same generator as the recovery experiment) - #
    # F_ST trend (Spearman rho vs ordinal bin) recovered from synthetic
    # assemblages as the injected group-closure strength s rises 0 -> 1, matched
    # to the record's clusters, bins, sizes, and rarefaction (script 21).
    ca, _ = res17.oriented_ca(counts_df)
    cdf = coords_df.dropna()
    have = list(cdf.index)
    cch = cdf[["Latitude", "Longitude"]].to_numpy(float)
    cch_c = cch - cch.mean(0)
    sil = {kk: mf.silhouette_mean(cch_c, mf._kmeans_labels(cch_c, kk, seed=7))
           for kk in range(2, 7)}
    kk = max(sil, key=sil.get)
    clh = mf._kmeans_labels(cch_c, kk, seed=7)
    binsh_s = pd.Series(pd.qcut(ca.reindex(have).to_numpy(), r21.N_BINS,
                                labels=False, duplicates="drop"))
    ok = (~binsh_s.isna()).to_numpy()
    have = list(np.array(have)[ok])
    clh = clh[ok]
    cch_c = cch_c[ok]
    binsh = np.asarray(binsh_s[ok].astype(int))
    Narr = counts_df.reindex(have).to_numpy(float).sum(1)
    real_have = counts_df.reindex(have).to_numpy(float)
    Khave = counts.shape[1]
    s_grid = r21.S_GRID
    means, los, his = [], [], []
    for s in s_grid:
        vals = []
        for seed in range(40):
            rg = np.random.default_rng(1000 + seed)
            grid = r21.emergence_profiles(kk, Khave, float(s), rg)
            M = r21.gen_synth_M(grid, clh, binsh, Narr, rg)
            Mr = r21.rarefy(M, r21.NRARE, rg)
            rho = r21.sig_rhos(Mr, clh, binsh, cch_c, which=["fst"])["fst"]
            if np.isfinite(rho):
                vals.append(rho)
        vals = np.array(vals)
        means.append(vals.mean()); los.append(np.percentile(vals, 2.5)); his.append(np.percentile(vals, 97.5))
    means, los, his = map(np.array, (means, los, his))
    obs_vals = []
    for b in range(200):
        rg = np.random.default_rng(5000 + b)
        Mr = r21.rarefy(real_have, r21.NRARE, rg)
        rho = r21.sig_rhos(Mr, clh, binsh, cch_c, which=["fst"])["fst"]
        if np.isfinite(rho):
            obs_vals.append(rho)
    obs_fst = float(np.mean(obs_vals))
    s_star = 0.5     # detection threshold from the calibrated recovery (Figure 4)

    # ---- figure ----------------------------------------------------------- #
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8})
    fig = plt.figure(figsize=(7.0, 5.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95], hspace=0.5, wspace=0.3)

    panels = [
        ("If the phase were a bounded group", Mg, "what we would expect"),
        ("Spatially structured drift", Md, "continuous interaction"),
        ("Observed assemblages", counts, "the real data"),
    ]
    for j, (title, M, sub) in enumerate(panels):
        ax = fig.add_subplot(gs[0, j])
        xy = demo.mds2(M)
        for gi, col in zip((0, 1), (C0, C1)):
            m = glab == gi
            ax.scatter(xy[m, 0], xy[m, 1], s=22, c=col, edgecolor="white", linewidth=0.4)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=8)
        ax.set_xlabel(sub, fontsize=7, color="#555555")
        for spi in ax.spines.values():
            spi.set_edgecolor("#bbbbbb")

    # recovery panel spanning the bottom row
    axr = fig.add_subplot(gs[1, :])
    ymin = float(min(los.min(), obs_fst)) - 0.05
    ymax = float(max(his.max(), 0.0)) + 0.20          # headroom for the annotation
    axr.set_ylim(ymin, ymax)
    axr.set_xlim(0, 1)
    axr.axvspan(0, s_star, color="#cccccc", alpha=0.40, zorder=0)
    axr.fill_between(s_grid, los, his, color=C1, alpha=0.20, zorder=1)
    axr.plot(s_grid, means, color=C1, lw=1.8, zorder=2,
             label="F_ST signal recovered from synthetic groups")
    axr.axhline(obs_fst, color="black", lw=1.4, ls="--", zorder=3,
                label=f"observed F_ST signal ({obs_fst:+.02f})")
    axr.axvline(s_star, color="#444444", lw=1.0, ls=":", zorder=3)
    # annotation high in the shaded (below-resolution) region, clear of the curve
    axr.text(s_star / 2, ymax - 0.02,
             "below resolution:\nweak group and drift\nare indistinguishable",
             ha="center", va="top", fontsize=6.5, color="#333333")
    # threshold label just right of the s* line, near the bottom, clear of the curve
    axr.text(s_star + 0.015, ymin + 0.03, f"detection threshold s* = {s_star:.2f}",
             ha="left", va="bottom", fontsize=6.5, color="#333333")
    axr.set_xlabel("strength of group closure injected into synthetic assemblages "
                   "(0 = drift, 1 = strong bounded groups)")
    axr.set_ylabel("recovered F_ST trend (Spearman rho)")
    # legend below the panel so it cannot overlap the curve or the annotations
    axr.legend(fontsize=6.5, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, -0.28), ncol=2)
    for spi in ("top", "right"):
        axr.spines[spi].set_visible(False)

    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"s* = {s_star:.2f}; observed F_ST trend = {obs_fst:+.3f}; "
          f"recovered at s=0.5 = {means[np.argmin(abs(s_grid - 0.5))]:+.3f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
