"""29_concept_figure.py — the conceptual figure: what we are trying to do.

Illustrates the test graphically, from synthetic data matched to the record:
  Row 1: what the assemblages SHOULD look like if the Parkin phase were a bounded
         interaction group (sharp clusters), what spatially structured drift
         produces instead (a continuous blur), and what the OBSERVED assemblages
         look like (a blur, like drift).
  Row 2: why the determination cannot be made. Cultural F_ST recovered from
         synthetic assemblages as the strength of group closure rises from 0
         (drift) to 1 (strong groups), with the observed value and the
         posterior over injected closure strength marked. The shaded band is
         the 95 percent credible interval for that strength: it reaches from
         near zero to moderate closure, so weak closure and pure drift are both
         consistent with the observation at this record's size and
         time-averaging. An earlier version drew a hard detection threshold
         s* = 0.5 here and shaded everything below it; that was a frequentist
         decision rule, withdrawn under rule 18 and replaced by the posterior
         that script 53 computes.

Read-only on the manuscript. Writes figures/fig3_concept.png.

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
from mls_emergence.signatures.sampling import rarefy  # noqa: E402
res17 = importlib.import_module("17_basin_results")

OUT = ROOT / "figures" / "fig3_concept.png"
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
    Khave = counts.shape[1]
    s_grid = r21.S_GRID
    means, los, his = [], [], []
    for s in s_grid:
        vals = []
        for seed in range(40):
            rg = np.random.default_rng(1000 + seed)
            grid = r21.emergence_profiles(kk, Khave, float(s), rg)
            M = r21.gen_synth_M(grid, clh, binsh, Narr, rg)
            Mr = rarefy(M, r21.NRARE, rg)
            rho = r21.sig_rhos(Mr, clh, binsh, cch_c, which=["fst"])["fst"]
            if np.isfinite(rho):
                vals.append(rho)
        vals = np.array(vals)
        means.append(vals.mean()); los.append(np.percentile(vals, 2.5)); his.append(np.percentile(vals, 97.5))
    means, los, his = map(np.array, (means, los, his))
    # The observed trend comes from the one shared definition in script 21
    # rather than being recomputed here. This panel used to run its own 200
    # rarefactions from a different seed and print -0.02 while the main text
    # printed +0.01: the same statistic, on opposite sides of zero, because its
    # Monte Carlo error (single-draw SD about 0.26) swamps its magnitude.
    trend = r21.empirical_fst_trend()
    # The RAREFACTION spread, not the Monte Carlo interval of the mean. This
    # panel used to shade the latter, which is about thirty times narrower and
    # made an unresolved trend look resolved; CLAUDE.md rule 6 is explicit that
    # the spread across draws, not the precision of their average, is what
    # belongs next to this number.
    obs_fst, obs_lo, obs_hi = trend["mean"], trend["p_lo"], trend["p_hi"]

    # The withdrawn apparatus drew a hard "detection threshold s* = 0.50" here
    # and shaded everything below it as unresolvable. That is a frequentist
    # decision rule and rule 18 removed it; script 53 replaced it with a
    # posterior over the injected closure strength, which is what this panel
    # now shows.
    # Same staleness guard as Figure 4; see 53.load_posterior.
    post = importlib.import_module("53_closure_strength_posterior").load_posterior()
    s_fine, p_w1 = post["s_grid"], post["post_w1"]
    cdf_w1 = np.cumsum(p_w1) / p_w1.sum()
    s_med = float(np.interp(0.5, cdf_w1, s_fine))
    s_lo = float(np.interp(0.025, cdf_w1, s_fine))
    s_hi = float(np.interp(0.975, cdf_w1, s_fine))

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
        # MDS coordinates are dimensionless (arbitrary scale/sign/rotation), so
        # the axes carry no ticks; label them as the ordination dimensions. The
        # scenario subtitle moves into the title so nothing is lost.
        ax.set_title(f"{title}\n({sub})", fontsize=8)
        ax.set_xlabel("MDS 1", fontsize=7, labelpad=2)
        ax.set_ylabel("MDS 2", fontsize=7, labelpad=2)
        for spi in ax.spines.values():
            spi.set_edgecolor("#bbbbbb")

    # recovery panel spanning the bottom row
    axr = fig.add_subplot(gs[1, :])
    # The observed band is the rarefaction spread now, which is far wider than
    # the Monte Carlo interval it replaced, so the limits have to admit it or
    # the band is silently clipped.
    ymin = float(min(los.min(), obs_lo)) - 0.05
    ymax = float(max(his.max(), obs_hi, 0.0)) + 0.20  # headroom for the annotation
    axr.set_ylim(ymin, ymax)
    axr.set_xlim(0, 1)
    # Credible interval for the injected closure strength given the observation.
    axr.axvspan(s_lo, s_hi, color="#cccccc", alpha=0.40, zorder=0)
    axr.fill_between(s_grid, los, his, color=C1, alpha=0.20, zorder=1)
    axr.plot(s_grid, means, color=C1, lw=1.8, zorder=2,
             label="$F_{ST}$ trend recovered from synthetic groups")
    axr.axhline(obs_fst, color="black", lw=1.4, ls="--", zorder=3,
                label=f"observed trend ({obs_fst:+.03f}; shaded band is the "
                      f"2.5th-97.5th rarefaction percentiles)")
    axr.fill_between([0, 1], obs_lo, obs_hi, color="black", alpha=0.12, zorder=2,
                     linewidth=0)
    axr.axvline(s_med, color="#444444", lw=1.0, ls=":", zorder=3)
    axr.text((s_lo + s_hi) / 2, ymax - 0.02,
             "closure strengths consistent\nwith the observation\n"
             f"(95% credible interval {s_lo:.2f}-{s_hi:.2f})",
             ha="center", va="top", fontsize=6.5, color="#333333")
    axr.text(s_med + 0.015, ymin + 0.03, f"posterior median s = {s_med:.2f}",
             ha="left", va="bottom", fontsize=6.5, color="#333333")
    axr.set_xlabel("strength of group closure injected into synthetic assemblages "
                   "(0 = drift, 1 = strong bounded groups)")
    axr.set_ylabel("recovered $F_{ST}$ trend (Spearman rho)")
    # legend below the panel so it cannot overlap the curve or the annotations
    axr.legend(fontsize=6.5, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, -0.28), ncol=2)
    for spi in ("top", "right"):
        axr.spines[spi].set_visible(False)

    from figstyle import save_all
    save_all(fig, OUT)
    print(f"observed F_ST trend = {obs_fst:+.4f} "
          f"[{obs_lo:+.4f}, {obs_hi:+.4f}] over {trend['n_draws']} rarefactions "
          f"(single-draw SD {trend['sd']:.3f}); closure-strength posterior median "
          f"{s_med:.2f}, 95% interval [{s_lo:.2f}, {s_hi:.2f}]; "
          f"recovered at s=0.5 = {means[np.argmin(abs(s_grid - 0.5))]:+.3f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
