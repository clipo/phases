#!/usr/bin/env python3
"""The phase tests, replicated on Mainfort's (2003) table under his own phases.

The primary analysis tests three phases on the matrix built from the survey's
counts and Lipo's compilation. Mainfort (2003, Table 1) assembled a different
tally of the same region, 39 sites of 800 or more sherds, and assigned them to
phases under a DIFFERENT CONFIGURATION: Smith's (1990) Horseshoe Lake phase
split out of Walls, Hollywood and Commerce moved from Walls to Kent, Castile
Landing to Parkin, Carson Lake to Walls, Jeter to Walls. Nine of the 36 sites
in both schemes carry different assignments. Two workers, the same sites,
different lines.

That makes his table worth two things. First, a replication of the three
tests that carry the paper's account, on a second tally and a second scheme.
Second, a question the primary analysis cannot ask: where the two schemes put
a site in different phases, does its pottery favour either assignment beyond
what its position predicts? If the pottery is indifferent between two
experts' lines, the lines are a reading of a gradient.

WHAT IS RUN, at a minimum of 100 decorated sherds (stricter than the primary
75, because a replication should not be the looser test), on his phase labels,
with Tipton's one surviving site not carried as a phase of one:

  1. What recovers his scheme: site map alone, map plus composition at
     several weights, composition alone (as analyses/76_phase_recovery.py).
  2. His partition against size-matched partitions with boundaries placed at
     random and drawn to be compact (as analyses/74_phase_partition_test.py).
  3. The contested sites. For each of the nine, the chi-square distance from
     the site's profile to the pooled profile of each candidate phase (the
     site itself left out of both), and the same for the site's straight-line
     distance to each phase's centroid. A site whose pottery is closer to the
     phase it is NOT assigned to under one scheme is a line the pottery does
     not support.

WHAT THIS IS NOT. His counts are not independent of the survey's for the 21
shared sites (11 identical, 10 larger), so this is a different sampling of the
same collections plus new Tennessee material, not an independent record. The
drift comparison is not repeated here: the model is calibrated to the primary
matrix, and recalibrating it to his would be a second paper.

Usage:
    python analyses/83_mainfort_replication.py [--min 100] [--seeds 40] [--alt 1000]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "mainfort_replication.md"
FIG = "fig17_mainfort_replication"
WEIGHTS = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, np.inf]
MIN_PHASE_MEMBERS = 2
# The nine sites the two schemes assign differently, with the Figure 1 phase.
CONTESTED = {"Beck": "Walls", "Belle Meade": "Walls", "Mound Place": "Walls", "Young": "Walls",
             "Hollywood": "Walls", "Commerce": "Walls", "Castile Lg.": "Kent",
             "Carson Lake": "Nodena", "Jeter": "Tipton"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min", type=int, default=100)
    ap.add_argument("--seeds", type=int, default=40)
    ap.add_argument("--alt", type=int, default=1000)
    args = ap.parse_args()

    import matplotlib.pyplot as plt
    from mls_emergence.dataio.matrix import read_mainfort_replication
    fs = importlib.import_module("figstyle")
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    t74 = importlib.import_module("74_phase_partition_test")
    t76 = importlib.import_module("76_phase_recovery")

    counts, coords, phases = read_mainfort_replication(args.min)
    sizes = phases.value_counts()
    small = sizes[sizes < MIN_PHASE_MEMBERS].index.tolist()
    if small:
        print(f"not carried as a phase (fewer than {MIN_PHASE_MEMBERS} sites): "
              + ", ".join(f"{p} ({', '.join(phases.index[phases == p])})" for p in small))
    keep = ~phases.isin(small)
    counts, coords, phases = counts[keep], coords[keep], phases[keep]
    names = list(counts.index)
    m = counts.to_numpy(float)
    xy = coords.to_numpy(float)
    plist = sorted(set(phases))
    pidx = np.array([plist.index(p) for p in phases])
    k = len(plist)
    pts = t74.km_xy(xy)
    G = t76.unit_scale(pts)
    p = m / m.sum(1, keepdims=True)
    C = t76.unit_scale(t76.composition_features(p, "chisq"))

    # 1. What recovers his scheme.
    rec = []
    for w in WEIGHTS:
        F = G if w == 0 else (C if not np.isfinite(w) else np.column_stack([G, np.sqrt(w) * C]))
        a = [t76.ari(pidx, mf._kmeans_labels(F, k, seed=s)) for s in range(args.seeds)]
        rec.append((w, float(np.median(a)), float(np.min(a)), float(np.max(a))))

    # 2. His partition against size-matched alternatives.
    fst_phase = float(rev.fst_by(m.astype(int), pidx))
    sz = np.bincount(pidx)
    slot = np.repeat(np.arange(k), sz)
    rng = np.random.default_rng(83)
    alt, opt, exact = [], [], 0
    for _ in range(args.alt):
        seeds = pts[rng.choice(len(pts), size=k, replace=False)]
        a = t74.assign_exact(pts, seeds, slot)
        alt.append(rev.fst_by(m.astype(int), a))
        b = t74.local_search(pts, a, sz)
        opt.append(rev.fst_by(m.astype(int), b))
        exact += int(t74.ari(pidx, b) > 0.9999)
    alt, opt = np.array(alt), np.array(opt)
    pct_alt = 100 * float(np.mean(alt < fst_phase))
    pct_opt = 100 * float(np.mean(opt < fst_phase))
    rnd_ari = [t76.ari(pidx, t74.assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot))
               for _ in range(400)]

    # 3. The contested sites.
    cm = p.mean(0)
    w = np.sqrt(np.where(cm > 0, cm, 1.0))
    contested = []
    for site, fig1 in CONTESTED.items():
        if site not in names:
            contested.append((site, "under the minimum, not tested", "", "", "", "", ""))
            continue
        i = names.index(site)
        his = phases.iloc[i]
        row = [site, his, fig1]
        for cand in (his, fig1):
            members = [j for j in range(len(names)) if phases.iloc[j] == cand and j != i]
            if not members:
                row += ["no other member", ""]
                continue
            pooled = m[members].sum(0) / m[members].sum()
            cd = float(np.sqrt((((p[i] - pooled) / w) ** 2).sum()))
            gd = float(np.linalg.norm(pts[i] - pts[members].mean(0)))
            row += [f"{cd:.3f}", f"{gd:.0f}"]
        contested.append(tuple(row))

    # Figure: recovery curve and the contested sites.
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = fig.add_gridspec(1, 2, wspace=0.35)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(WEIGHTS))
    ax.plot(x, [r[1] for r in rec], color="0.1", marker="o", ms=3.5, lw=1.4, label="Mainfort's scheme")
    ax.fill_between(x, [r[2] for r in rec], [r[3] for r in rec], color="0.85", label="range over initializations")
    ax.axhline(np.median(rnd_ari), color="0.55", ls="-.", lw=1, label="same sizes, boundaries at random")
    ax.set_xticks(x)
    ax.set_xticklabels(["map\nalone" if w == 0 else ("pottery\nalone" if not np.isfinite(w) else f"{w:g}")
                        for w in WEIGHTS], fontsize=6)
    ax.set_xlabel("weight on composition relative to geography")
    ax.set_ylabel("agreement with Mainfort's phases\n(adjusted Rand index)")
    ax.legend(fontsize=5.5, frameon=False, loc="lower left")
    fs.panel_label(ax, "A")
    axB = fig.add_subplot(gs[0, 1])
    tested = [c for c in contested if c[1] != "under the minimum, not tested" and c[3] != "no other member" and c[5] != "no other member"]
    ys = np.arange(len(tested))
    axB.barh(ys - 0.18, [float(c[3]) for c in tested], height=0.36, color="0.25", label="to Mainfort's phase")
    axB.barh(ys + 0.18, [float(c[5]) for c in tested], height=0.36, color="0.7", label="to the Figure 1 phase")
    axB.set_yticks(ys)
    axB.set_yticklabels([f"{c[0]}\n({c[1]} / {c[2]})" for c in tested], fontsize=5.5)
    axB.invert_yaxis()
    axB.set_xlabel("ceramic distance to the phase's pooled profile\n(site left out)")
    axB.legend(fontsize=5.5, frameon=False, loc="lower right")
    fs.panel_label(axB, "B")
    png = fs.save_all(fig, FIG, close=True)

    L = ["# The phase tests on Mainfort's (2003) table, under his phases", "",
         f"{len(names)} sites with at least {args.min} decorated sherds, "
         f"{int(m.sum()):,} decorated sherds, {k} phases "
         f"({', '.join(f'{pl} {int((pidx == i).sum())}' for i, pl in enumerate(plist))})"
         + (f"; not carried as a phase: {', '.join(small)}." if small else "."), "",
         "## 1. What recovers his scheme", "",
         "| weight on composition | ARI, median over initializations | range |", "|---|---|---|"]
    for w_, med, lo, hi in rec:
        L.append(f"| {'map alone' if w_ == 0 else ('pottery alone' if not np.isfinite(w_) else f'{w_:g}')} | "
                 f"**{med:.3f}** | {lo:.3f}-{hi:.3f} |")
    L += [f"| same sizes, boundaries at random (400 draws) | {np.median(rnd_ari):.3f} | |", "",
          "## 2. His partition against size-matched alternatives", "",
          f"Cultural F_ST under his phases: **{fst_phase:.4f}**. Against {args.alt} partitions with his "
          f"group sizes: seeds at random, median {np.median(alt):.4f}, his sits at the "
          f"**{t74.ordinal(pct_alt)} percentile**; drawn to be compact, median {np.median(opt):.4f}, "
          f"the **{t74.ordinal(pct_opt)} percentile**. The compactness search reproduced his partition "
          f"exactly {exact} times in {args.alt}.", "",
          "## 3. The nine sites the two schemes assign differently", "",
          "Chi-square distance from the site's profile to each candidate phase's pooled profile with the "
          "site left out, and straight-line km to that phase's centroid.", "",
          "| site | Mainfort 2003 | Figure 1 | ceramic, to his | km, to his | ceramic, to Figure 1's | km, to Figure 1's |",
          "|---|---|---|---|---|---|---|"]
    for c in contested:
        L.append("| " + " | ".join(str(x) for x in c) + " |")
    L += ["", "A site whose pottery is closer to the phase it is NOT assigned to under one of the schemes "
              "is a line that",
          "scheme drew where the pottery does not go. Read the ceramic columns against each other, and "
          "against the km",
          "columns, which say what proximity alone would predict.", "",
          "Not repeated here: the drift comparison, which is calibrated to the primary matrix. His counts "
          "are not",
          "independent of the survey's for the 21 shared sites, so this is a second sampling of the same "
          "collections",
          "plus new Tennessee material, not an independent record.", "",
          f"Figure written to {png.name} and its siblings.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[2:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
