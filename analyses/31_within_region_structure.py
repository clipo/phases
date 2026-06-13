"""31_within_region_structure.py — within-region interaction structure: CMV vs LMV.

Demonstrates that NEITHER region's decorated-ceramic assemblages form bounded
interaction groups: within each region the structure is drift-consistent, the
between-cluster cultural F_ST sitting at the neutral-drift level and far below the
bounded-groups expectation. The CMV is analyzed on Mississippian types only, so
the comparison is synchronic-to-synchronic and not driven by the Woodland-to-
Mississippian composition gradient.

Top row: multidimensional scaling of decorated Brainerd-Robinson dissimilarity
for the LMV St. Francis basin and the CMV southeast-Missouri assemblages, points
colored by within-region spatial cluster, annotated with the observed F_ST.
Bottom: observed between-cluster F_ST for each region against the 95% nulls from a
neutral spatial-drift model and a bounded-groups model run on the real
coordinates (time-averaged, rarefied to the observed counts).

Read-only on the manuscript. Writes figures/fig9_within_region.png.

Usage: .venv/bin/python analyses/31_within_region_structure.py
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

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
demo = importlib.import_module("25_drift_vs_groups_demo")
g26 = importlib.import_module("26_cmv_phase_groupness")

OUT = ROOT / "figures" / "fig9_within_region.png"
C_DRIFT = "#0072B2"
C_GROUP = "#D55E00"
LEN, MB = 24.0, 0.02


def fst_by(counts, lab):
    pools = np.array([counts[lab == c].sum(0) for c in np.unique(lab)
                      if counts[lab == c].sum() > 0])
    return mf.cultural_fst(pools) if len(pools) >= 2 else np.nan


def region_nulls(counts, coords, lab, n_seed=80):
    K = counts.shape[1]
    tot = counts.sum(1)
    Wd = demo.drift_weights(coords, length_km=LEN)
    Wg = demo.group_weights(coords, lab, length_km=LEN, leak=0.03)
    drift, grp = [], []
    for s in range(n_seed):
        Md = demo.simulate(Wd, K=K, total_per_node=tot, m_between=MB, seed=200 + s)
        Mg = demo.simulate(Wg, K=K, total_per_node=tot, m_between=MB, seed=200 + s)
        fd, fg = fst_by(Md, lab), fst_by(Mg, lab)
        if np.isfinite(fd):
            drift.append(fd)
        if np.isfinite(fg):
            grp.append(fg)
    return np.array(drift), np.array(grp)


def load_lmv():
    cdf, codf = mf._load_curated()
    counts = cdf.to_numpy(float)
    coords = codf[["Latitude", "Longitude"]].to_numpy(float)
    return counts, coords


def load_cmv_miss():
    M, _ = g26.load_cmv()
    miss = [c for c in M.columns if c not in g26.WOODLAND]
    Mm = M[miss]
    Mm = Mm[Mm.sum(1) >= 30]
    cc = pd.read_csv(ROOT / "data" / "processed" / "williams1954_cmv_coords.tsv",
                     sep="\t").set_index("site_name")
    sites = [a.split(" [")[0] for a in Mm.index]
    lat = np.array([cc.loc[s, "lat"] for s in sites], float)
    lon = np.array([cc.loc[s, "lon"] for s in sites], float)
    return Mm.to_numpy(float), np.column_stack([lat, lon])


def main():
    regions = []
    for name, (counts, coords) in [("LMV St. Francis basin", load_lmv()),
                                   ("CMV southeast Missouri (Mississippian)", load_cmv_miss())]:
        cc_c = coords - coords.mean(0)
        # pick k by silhouette over 2..4
        sil = {k: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, k, seed=7))
               for k in range(2, 5)}
        k = max(sil, key=sil.get)
        lab = mf._kmeans_labels(cc_c, k, seed=7)
        obs = fst_by(counts, lab)
        drift, grp = region_nulls(counts, coords, lab)
        regions.append(dict(name=name, counts=counts, lab=lab, k=k, obs=obs,
                            drift=drift, grp=grp, n=len(counts)))

    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8})
    fig = plt.figure(figsize=(7.0, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85], hspace=0.42, wspace=0.28)
    pal = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]

    for j, R in enumerate(regions):
        ax = fig.add_subplot(gs[0, j])
        xy = demo.mds2(R["counts"])
        for c in np.unique(R["lab"]):
            m = R["lab"] == c
            ax.scatter(xy[m, 0], xy[m, 1], s=24, c=pal[c % len(pal)],
                       edgecolor="white", linewidth=0.4)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{R['name']}\n(n={R['n']}, {R['k']} spatial clusters, "
                     f"F_ST={R['obs']:+.03f})", fontsize=7.5)
        for sp in ax.spines.values():
            sp.set_edgecolor("#bbbbbb")

    # bottom: observed F_ST vs drift / bounded-groups nulls, both regions
    axb = fig.add_subplot(gs[1, :])
    ypos = [1, 0]
    for R, y in zip(regions, ypos):
        dlo, dhi = np.percentile(R["drift"], [2.5, 97.5])
        glo, ghi = np.percentile(R["grp"], [2.5, 97.5])
        axb.plot([dlo, dhi], [y + 0.10, y + 0.10], color=C_DRIFT, lw=6, alpha=0.5,
                 solid_capstyle="butt")
        axb.plot([glo, ghi], [y - 0.10, y - 0.10], color=C_GROUP, lw=6, alpha=0.5,
                 solid_capstyle="butt")
        axb.plot([R["obs"]], [y], "D", ms=8, color="black", zorder=5)
        axb.text(R["obs"], y + 0.22, f"observed {R['obs']:+.03f}", ha="center",
                 fontsize=6.5)
    axb.set_yticks(ypos)
    axb.set_yticklabels(["LMV basin", "CMV (Miss.)"], fontsize=8)
    axb.set_ylim(-0.6, 1.6)
    axb.set_xlabel("between-cluster cultural F_ST")
    from matplotlib.patches import Patch
    axb.legend(handles=[Patch(color=C_DRIFT, alpha=0.5, label="neutral spatial-drift null (95%)"),
                        Patch(color=C_GROUP, alpha=0.5, label="bounded-groups null (95%)"),
                        plt.Line2D([], [], marker="D", color="black", ls="", label="observed")],
               fontsize=6.5, frameon=False, loc="upper right", ncol=1)
    for sp in ("top", "right", "left"):
        axb.spines[sp].set_visible(False)

    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    for R in regions:
        print(f"{R['name']}: obs F_ST={R['obs']:+.3f} drift95=[{np.percentile(R['drift'],2.5):+.3f},"
              f"{np.percentile(R['drift'],97.5):+.3f}] groups95=[{np.percentile(R['grp'],2.5):+.3f},"
              f"{np.percentile(R['grp'],97.5):+.3f}]")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
