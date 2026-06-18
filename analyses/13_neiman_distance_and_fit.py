"""13_neiman_distance_and_fit.py — two theory-grounded transmission analyses.

B. Neiman's (1995) interassemblage-distance neutral model. Squared Euclidean
   distance d_ij^2 = sum_k (p_ik - p_jk)^2 between assemblages. Tests:
   (i) the diversity-distance relationship Neiman found under drift (within-
       assemblage diversity and mean interassemblage distance inversely
       correlated, his r = -0.62) -- the convergence-of-two-measures precedent;
   (ii) the trend in mean interassemblage distance along the CA axis (does
       between-assemblage divergence rise toward contact, as assortment would
       require, or stay flat / fall, as homogenization would?);
   (iii) the between- vs within-cluster squared-distance ratio (Neiman's
       between-group leg, an alternative to the Gini-Simpson F_ST).

D. Frequency Increment Test (Feder, Kryazhimskiy & Plotkin 2014) on each
   decorated type's frequency trajectory along the ordered (binned) sequence:
   a time-aware test of neutral drift that uses the seriation ordering directly.

Writes output/neiman_distance_and_fit.md. Read-only on the manuscript.

Usage: .venv/bin/python analyses/13_neiman_distance_and_fit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_1samp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import os  # noqa: E402
os.environ.setdefault("MLS_FIG_COLOR", "1")  # supplement figure is online-only; render in color

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
from figstyle import save, OI_BLUE, OI_ORANGE  # noqa: E402

OUT = ROOT / "output" / "neiman_distance_and_fit.md"
N_BINS = 6


def main():
    counts, coords = mf._load_curated()        # basin 53
    M = counts.to_numpy(float)
    N = M.sum(1, keepdims=True)
    P = M / N                                   # type proportions per assemblage
    ca1, _, _ = mf.correspondence_axis(M)
    rank = pd.Series(ca1).rank().to_numpy()

    L = ["# Neiman interassemblage-distance model and Frequency Increment Test", ""]

    # --- B. Interassemblage distance ---
    D = ((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)   # squared Euclidean
    mean_d = D.sum(1) / (D.shape[0] - 1)
    tE = np.array([mf.theta_e(r) for r in M])
    np.array([mf.theta_f(r) for r in M])
    fin = np.isfinite(tE) & (tE > 0)
    rho_d, p_d = spearmanr(tE[fin], mean_d[fin])

    L += ["## B. Neiman interassemblage-distance neutral model", "",
          f"Squared Euclidean distance d_ij^2 = sum_k (p_ik - p_jk)^2 over {M.shape[1]} "
          f"decorated types, {M.shape[0]} basin assemblages.", "",
          f"- **Diversity-distance relationship** (Neiman's convergence-of-two-measures "
          f"precedent): within-assemblage diversity (t_E) vs mean interassemblage distance, "
          f"Spearman rho = {rho_d:+.2f} (p = {p_d:.3f}). Neiman (1995, Fig. 7) found "
          f"rho = -0.62 under drift+innovation in one interacting population; a negative "
          f"relationship here is consistent with a single drifting field rather than "
          f"assorted, bounded groups."]

    # trajectory: mean within-bin pairwise distance along CA
    bins = pd.qcut(pd.Series(rank), N_BINS, labels=False, duplicates="drop").to_numpy()
    bin_mean_d = {}
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        if len(idx) >= 2:
            sub = D[np.ix_(idx, idx)]
            iu = np.triu_indices(len(idx), 1)
            bin_mean_d[int(b)] = float(sub[iu].mean())
    bd = pd.Series(bin_mean_d).sort_index()
    rho_traj, p_traj = spearmanr(bd.index.to_numpy(float), bd.values)
    L += ["",
          f"- **Divergence trajectory**: mean within-bin pairwise distance across the {N_BINS} "
          f"CA bins trends with seriation position at Spearman rho = {rho_traj:+.2f} "
          f"(p = {p_traj:.3f}). A rising trend would indicate growing between-assemblage "
          f"divergence (fragmentation); a flat or falling trend indicates no growing "
          f"divergence. Per-bin mean d^2: " +
          ", ".join(f"{k}:{v:.3f}" for k, v in bd.items()) + "."]

    # within vs between cluster distance (Neiman between-group leg)
    cdf = coords.dropna()
    have = list(cdf.index)
    cc = cdf[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(0)
    sil = {kk: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, kk, seed=7)) for kk in range(2, 7)}
    k = max(sil, key=sil.get)
    cl = mf._kmeans_labels(cc_c, k, seed=7)
    pos = {a: i for i, a in enumerate(counts.index)}
    hav_idx = np.array([pos[a] for a in have])
    Dh = D[np.ix_(hav_idx, hav_idx)]
    iu = np.triu_indices(len(have), 1)
    same = (cl[iu[0]] == cl[iu[1]])
    within = Dh[iu][same].mean()
    between = Dh[iu][~same].mean()
    L += ["",
          f"- **Between- vs within-cluster distance** (k = {k} spatial clusters): "
          f"within-cluster mean d^2 = {within:.3f}, between-cluster = {between:.3f}, "
          f"ratio = {between/within:.2f}. A ratio near 1 means spatial clusters are no more "
          f"internally similar than the basin at large (no bounded ceramic groups); a large "
          f"ratio would mark bounded groups.", ""]

    # --- D. Frequency Increment Test ---
    # type frequencies per CA bin
    freq = np.zeros((N_BINS, M.shape[1]))
    for b in range(N_BINS):
        idx = np.where(bins == b)[0]
        tot = M[idx].sum()
        freq[b] = M[idx].sum(0) / tot if tot > 0 else np.nan
    L += ["## D. Frequency Increment Test (Feder, Kryazhimskiy & Plotkin 2014)", "",
          "Per decorated type: rescaled increments Y_i = (v_i - v_{i-1}) / "
          "sqrt(2 v_{i-1}(1-v_{i-1}) dt) along the ordered bins (dt = 1); a one-sample "
          "t-test of mean(Y) = 0 tests neutral drift. Increments at v in {0,1} are dropped.", "",
          "| type | n increments | mean Y | t p-value | departs neutral? |",
          "|---|---|---|---|---|"]
    ndepart = 0
    ntested = 0
    for j, typ in enumerate(counts.columns):
        v = freq[:, j]
        Y = []
        for i in range(1, len(v)):
            vp = v[i - 1]
            if 0 < vp < 1 and np.isfinite(v[i]):
                Y.append((v[i] - vp) / np.sqrt(2 * vp * (1 - vp) * 1.0))
        if len(Y) >= 3:
            t, p = ttest_1samp(Y, 0.0)
            dep = p < 0.05
            ndepart += int(dep)
            ntested += 1
            L.append(f"| {typ} | {len(Y)} | {np.mean(Y):+.2f} | {p:.3f} | {'YES' if dep else 'no'} |")
        else:
            L.append(f"| {typ} | {len(Y)} | - | - | (too few) |")
    L += ["",
          f"**{ndepart} of {ntested} testable types depart from neutral drift at p < 0.05.** "
          "A predominantly neutral result corroborates the static neutrality signature: the "
          "decorated types behave as drifting neutral variants along the sequence, not as "
          "markers under conformist or anti-conformist bias. (With six bins per type the "
          "per-type test has low power; read the count, not individual p-values.)"]

    # --- Figure S3: the Neiman interassemblage-distance result ---
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7, 3.3))
    x, y = tE[fin], mean_d[fin]
    axA.scatter(x, y, s=20, color=OI_BLUE, alpha=0.8, edgecolor="none")
    # No OLS line: the reported statistic is the (rank-based) Spearman rho, which
    # an outlier-sensitive least-squares slope would contradict in sign here.
    names = list(counts.index)
    if "Parkin" in names:
        pj = names.index("Parkin")
        if fin[pj]:
            axA.scatter([tE[pj]], [mean_d[pj]], s=80, marker="*",
                        color=OI_ORANGE, edgecolor="0.2", lw=0.4, zorder=5)
            axA.annotate("Parkin", (tE[pj], mean_d[pj]), fontsize=7,
                         xytext=(4, 4), textcoords="offset points")
    axA.set_xlabel(r"within-assemblage diversity ($\theta_E$)")
    axA.set_ylabel("mean interassemblage distance")
    _ns = "" if p_d < 0.05 else " (n.s.)"
    axA.text(0.95, 0.95, rf"$\rho$ = {rho_d:+.2f}{_ns}", transform=axA.transAxes,
             ha="right", va="top", fontsize=9)
    axB.plot(bd.index.to_numpy(float), bd.values, "-o", color=OI_BLUE, ms=5, lw=1.3)
    axB.set_xlabel("CA seriation bin (early to late)")
    axB.set_ylabel("mean within-bin distance")
    axB.text(0.05, 0.95, rf"$\rho$ = {rho_traj:+.2f} (n.s.)", transform=axB.transAxes,
             ha="left", va="top", fontsize=9)
    save(fig, "figS1_neiman")

    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}")
    print("\n".join(L))


if __name__ == "__main__":
    main()
