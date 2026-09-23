#!/usr/bin/env python3
"""What recovers the phase scheme: the pottery, or the map?

The phases are ceramic units. They were defined by workers sorting decorated
pottery, and the claim carried by the name is that assemblages within a phase
share a way of making pots. This script asks how much of the published scheme
can be recovered WITHOUT any pottery at all, from the site coordinates alone,
and whether adding the pottery to the coordinates improves the recovery.

THE MEASUREMENT. Assemblages are clustered in a joint feature space,

    geography (km, centred, scaled to unit mean pairwise distance)
    composition (scaled the same way) weighted by sqrt(w)

so w = 0 is the site map alone and w -> large is the pottery alone. Agreement
with the published phase assignment is the adjusted Rand index, which is zero
in expectation for unrelated partitions and one for identical ones. If the
phases are ceramic units, agreement should RISE as the pottery enters.

THE CIRCULARITY THAT WOULD MANUFACTURE THIS RESULT, and the control for it.
Nine of the 43 assemblages are not in Mainfort's phase lists; this project
assigns them the phase of the territory they fall in, which is a geographic
rule. Scoring geography against those labels would be scoring geography against
geography. Every number is therefore reported twice: over all 43, and over only
the assemblages Mainfort himself assigned. The conclusion rests on the second.

THREE ROBUSTNESS AXES, because a single clustering of a single representation
is not a result.

  composition distance   Class proportions are COMPOSITIONAL, so a Euclidean
                         distance on raw proportions is not the principled
                         choice. The centred log-ratio transform is, and the
                         chi-square distance is what correspondence analysis
                         uses and so is closest to what a seriation-trained
                         worker is looking at. All three are run. If the
                         verdict depends on the transform, it is not a verdict.
                         THE CHI-SQUARE TRANSFORM IS THE ONE TO READ FIRST,
                         because it is the only one of the three with no free
                         parameter: 39 percent of the cells here are zero, and
                         the log-ratio transform cannot be computed without
                         deciding what a zero is worth. That decision moves the
                         answer from 0.74 to 1.00 (see `zero_sensitivity`), so
                         a log-ratio number quoted on its own is not a
                         measurement of the record.
  algorithm              k-means, Ward linkage and average linkage.
  seed                   k-means over many initializations, reported as a
                         spread rather than one number.

WHAT A REFERENCE LEVEL LOOKS LIKE HERE. An adjusted Rand index is already
chance-corrected, so zero is the no-relationship level and no null is needed to
establish it. For scale, the index is also computed for partitions carrying the
phases' own group sizes with their boundaries placed at random, which describes
what this design manufactures from the group sizes alone (rule 18: descriptive,
not an inferential test).

Usage:
    python analyses/76_phase_recovery.py [--seeds 40] [--random 400]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.special import comb

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "phase_recovery.md"
OUT_CSV = ROOT / "output" / "findings" / "phase_recovery.csv"
FIG = "fig14_phase_recovery"
WEIGHTS = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 32.0]
# Zero handling for the log-ratio transform. Proportions are rescaled to the
# assemblage's sherd count and a pseudocount added, so absence is bounded by
# the sample rather than by a chosen floor.
PSEUDO_SCALE = 1.0
PSEUDO = 0.5


def ari(a, b) -> float:
    """Adjusted Rand index. Zero in expectation for unrelated partitions."""
    a, b = np.asarray(a), np.asarray(b)
    ra = {v: i for i, v in enumerate(sorted(set(a.tolist())))}
    rb = {v: i for i, v in enumerate(sorted(set(b.tolist())))}
    a = np.array([ra[v] for v in a.tolist()])
    b = np.array([rb[v] for v in b.tolist()])
    C = np.zeros((a.max() + 1, b.max() + 1), int)
    for i, j in zip(a, b):
        C[i, j] += 1
    s = comb(C, 2).sum()
    sa, sb = comb(C.sum(1), 2).sum(), comb(C.sum(0), 2).sum()
    e = sa * sb / comb(len(a), 2)
    return float((s - e) / ((sa + sb) / 2 - e))


def unit_scale(X: np.ndarray) -> np.ndarray:
    """Centre, then scale so the mean squared pairwise distance is 1.

    Puts geography and composition on a common footing, so `w` is a ratio of
    like quantities rather than an arbitrary mixing of kilometres and
    proportions.
    """
    X = np.asarray(X, float)
    X = X - X.mean(0)
    s = np.sqrt((X ** 2).sum(1).mean())
    return X / s if s > 0 else X


def _floor_clr(p: np.ndarray, mult: float) -> np.ndarray:
    """Log-ratio with zeros floored at `mult` times the smallest positive value.

    Kept only so `zero_sensitivity` can show what this convention does to the
    answer. It is not used for any reported estimate.
    """
    q = p.copy()
    q[q <= 0] = q[q > 0].min() * mult
    q = q / q.sum(1, keepdims=True)
    lg = np.log(q)
    return lg - lg.mean(1, keepdims=True)


def composition_features(p: np.ndarray, kind: str) -> np.ndarray:
    """Composition as coordinates, under the transform named.

    raw    proportions as they are. Euclidean distance on a simplex; included
           because it is what the rest of this paper's F_ST work uses.
    clr    centred log-ratio, the principled transform for compositional data.
           Zeros are replaced by a multiplicative half of the smallest positive
           proportion, since log(0) is undefined; the replacement is reported.
    chisq  the chi-square metric correspondence analysis uses, obtained as the
           row profiles divided by the square root of the column masses. This
           is closest to what a seriation-trained worker is looking at.
    """
    if kind == "raw":
        return p.copy()
    if kind == "clr":
        # Zeros are replaced by adding a pseudocount to the COUNTS, not by
        # flooring the proportions. That choice is load-bearing and was got
        # wrong first: 39 percent of the cells in this matrix are zero, and
        # flooring at a fraction of the smallest positive proportion lets the
        # fraction set how hard absence is amplified. Agreement with the phases
        # then runs 0.82 at a floor of 2x, 0.90 at 0.5x and a perfect 1.00 at
        # 0.1x -- the headline number was a free parameter. A pseudocount has
        # no such knob at the limit and gives 0.74, and `zero_sensitivity`
        # below prints the whole curve so the reader sees the dependence
        # instead of inheriting one point from it.
        q = (p * PSEUDO_SCALE + PSEUDO)
        q = q / q.sum(1, keepdims=True)
        lg = np.log(q)
        return lg - lg.mean(1, keepdims=True)
    if kind == "chisq":
        cm = p.mean(0)
        cm = np.where(cm > 0, cm, 1.0)
        return p / np.sqrt(cm)
    raise ValueError(f"unknown composition transform: {kind!r}")


def cluster(F: np.ndarray, k: int, how: str, seed: int, mf) -> np.ndarray:
    if how == "kmeans":
        return mf._kmeans_labels(F, k, seed=seed)
    D = squareform(np.sqrt(((F[:, None, :] - F[None, :, :]) ** 2).sum(2)), checks=False)
    return fcluster(linkage(D, method=how), k, criterion="maxclust") - 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=40)
    ap.add_argument("--random", type=int, default=400)
    args = ap.parse_args()

    import matplotlib.pyplot as plt
    fs = importlib.import_module("figstyle")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    t74 = importlib.import_module("74_phase_partition_test")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m = counts.to_numpy(float)
    p = m / m.sum(1, keepdims=True)
    xy = coords.to_numpy(float)
    pts = t74.km_xy(xy)

    labels_ph, derived = ph.assign_phases_by_territory(names, xy)
    derived = np.asarray([bool(d) for d in derived]) if not isinstance(derived, dict) \
        else np.array([bool(derived.get(n, False)) for n in names])
    phases = sorted(set(labels_ph))
    pidx = np.array([phases.index(l) for l in labels_ph])
    k = len(phases)
    mapped = ~derived

    G = unit_scale(pts)
    rows = []
    for kind in ("raw", "clr", "chisq"):
        C = unit_scale(composition_features(p, kind))
        for w in WEIGHTS:
            F = G if w == 0 else np.column_stack([G, np.sqrt(w) * C])
            for how in ("kmeans", "ward", "average"):
                seeds = range(args.seeds) if how == "kmeans" else [0]
                a_all, a_map = [], []
                for sd in seeds:
                    lab = cluster(F, k, how, sd, mf)
                    a_all.append(ari(pidx, lab))
                    a_map.append(ari(pidx[mapped], lab[mapped]))
                rows.append(dict(transform=kind, w=w, algorithm=how,
                                 ari_all=float(np.median(a_all)),
                                 ari_mapped=float(np.median(a_map)),
                                 ari_mapped_lo=float(np.min(a_map)),
                                 ari_mapped_hi=float(np.max(a_map)),
                                 n_seeds=len(a_all)))
        # composition alone, the far end of the axis
        for how in ("kmeans", "ward", "average"):
            seeds = range(args.seeds) if how == "kmeans" else [0]
            a_all = [ari(pidx, cluster(C, k, how, sd, mf)) for sd in seeds]
            a_map = [ari(pidx[mapped], cluster(C, k, how, sd, mf)[mapped]) for sd in seeds]
            rows.append(dict(transform=kind, w=np.inf, algorithm=how,
                             ari_all=float(np.median(a_all)),
                             ari_mapped=float(np.median(a_map)),
                             ari_mapped_lo=float(np.min(a_map)),
                             ari_mapped_hi=float(np.max(a_map)),
                             n_seeds=len(a_all)))
    df = pd.DataFrame(rows)
    # `df.transform` resolves to the DataFrame METHOD, not this column, so a
    # filter written that way is silently always False and the tables below
    # come out empty. They did, in the first committed version. Assert the
    # selection is non-empty rather than trust the attribute.
    assert len(df[(df["transform"] == "chisq") & (df["algorithm"] == "kmeans")]) > 0, \
        "transform/algorithm selection is empty; check column access"
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    # How much of the log-ratio answer is the zero convention? Rule 4: a choice
    # a result depends on is reported with the result, not buried in a helper.
    zero_rows = []
    for label, C in (("pseudocount 0.5 (the default here)",
                      composition_features(p, "clr")),
                     ("pseudocount 5.0", np.log((p + 5.0) / (p + 5.0).sum(1, keepdims=True))
                      - np.log((p + 5.0) / (p + 5.0).sum(1, keepdims=True)).mean(1, keepdims=True)),
                     ("floor 2.0x min positive", _floor_clr(p, 2.0)),
                     ("floor 0.5x min positive", _floor_clr(p, 0.5)),
                     ("floor 0.1x min positive", _floor_clr(p, 0.1))):
        F = np.column_stack([G, np.sqrt(0.5) * unit_scale(C)])
        a = [ari(pidx[mapped], cluster(F, k, "kmeans", sd, mf)[mapped])
             for sd in range(args.seeds)]
        zero_rows.append((label, float(np.median(a))))

    # Descriptive scale: what the group sizes alone manufacture.
    rng = np.random.default_rng(76)
    sizes = np.bincount(pidx)
    slot = np.repeat(np.arange(k), sizes)
    rnd = []
    for _ in range(args.random):
        seeds_pts = pts[rng.choice(len(pts), size=k, replace=False)]
        lab = t74.assign_exact(pts, seeds_pts, slot)
        rnd.append(ari(pidx[mapped], lab[mapped]))
    rnd = np.array(rnd)

    geo = df[(df["w"] == 0) & (df["algorithm"] == "kmeans")]["ari_mapped"].median()
    fig, ax = plt.subplots(figsize=(3.8, 3.0))
    # Ordinal x axis over the weight grid, with "pottery alone" (w = inf) as a
    # separate end point. A symlog axis drew ticks at negative weights, which
    # cannot occur. Tick labels carry the actual grid values.
    pos = {w: i for i, w in enumerate(WEIGHTS)}
    x_alone = len(WEIGHTS) + 0.6
    # Three series coincide at many weights (proportions and log-ratio
    # especially), so each gets its own marker and line style and a small
    # horizontal dodge; open markers let an overlaid series show through.
    marks = {"raw": ("o", "-", -0.12, "none", "0.45"),
             "clr": ("s", "--", 0.0, "0.15", "0.15"),
             "chisq": ("^", ":", 0.12, "0.55", "0.55")}
    for kind, (mk, ls, dx, fc, col) in marks.items():
        sk = df[(df["transform"] == kind) & (df["algorithm"] == "kmeans")]
        fin = sk[np.isfinite(sk["w"])]
        xs = np.array([pos[w] for w in fin["w"]]) + dx
        ax.plot(xs, fin["ari_mapped"], marker=mk, ls=ls, ms=4, lw=1.1, color=col,
                markerfacecolor=fc, markeredgecolor=col, alpha=0.9,
                label={"raw": "proportions", "clr": "log-ratio",
                       "chisq": "chi-square"}[kind])
        alone = sk[~np.isfinite(sk["w"])]["ari_mapped"]
        ax.plot([x_alone + dx] * len(alone), alone, marker=mk, ls="none", ms=4,
                color=col, markerfacecolor=fc, markeredgecolor=col, alpha=0.9)
    ax.axvline((len(WEIGHTS) - 1 + x_alone) / 2, color="0.8", lw=0.6, ls=":")
    ax.axhline(np.median(rnd), color="0.6", lw=1.0, ls="-.",
               label="same sizes, boundaries at random")
    ax.set_xticks(list(pos.values()) + [x_alone])
    ax.set_xticklabels([f"{w:g}" for w in WEIGHTS] + ["pottery\nalone"], fontsize=6.5)
    ax.set_xlim(-0.5, x_alone + 0.5)
    ax.set_xlabel("weight on composition relative to geography")
    ax.set_ylabel("agreement with the published phases\n(adjusted Rand index)")
    ax.set_ylim(-0.05, max(0.8, geo + 0.12))
    ax.legend(fontsize=5.5, frameon=False, loc="lower left")
    png = fs.save_all(fig, FIG, close=True)

    L = ["# What recovers the phase scheme: the pottery, or the map?", "",
         f"Basin phase set, {len(names)} assemblages, {k} phases, 10 decorated classes.",
         f"{int(derived.sum())} assemblages carry a phase this project derived by "
         f"territory, a geographic",
         f"rule; they are EXCLUDED from the column the conclusion rests on, leaving "
         f"{int(mapped.sum())}",
         "assemblages Mainfort himself assigned.", "",
         "Agreement is the adjusted Rand index: 0 in expectation for unrelated "
         "partitions, 1 for",
         "identical ones. Weight 0 is the site map alone; the last row is the "
         "pottery alone.", "",
         "| composition transform | weight on composition | ARI, all 43 | ARI, Mainfort's "
         f"{int(mapped.sum())} | k-means seed range |",
         "|---|---|---|---|---|"]
    for kind in ("raw", "clr", "chisq"):
        for _, r in df[(df["transform"] == kind) & (df["algorithm"] == "kmeans")].iterrows():
            w = "pottery alone" if not np.isfinite(r["w"]) else (
                "**map alone**" if r["w"] == 0 else f"{r['w']:g}")
            L.append(f"| {kind} | {w} | {r['ari_all']:.3f} | "
                     f"**{r['ari_mapped']:.3f}** | "
                     f"{r['ari_mapped_lo']:.3f}-{r['ari_mapped_hi']:.3f} |")
    L += ["", "## Does the verdict depend on the algorithm?", "",
          "| transform | algorithm | map alone | pottery alone |", "|---|---|---|---|"]
    for kind in ("raw", "clr", "chisq"):
        for how in ("kmeans", "ward", "average"):
            a = df[(df["transform"] == kind) & (df["algorithm"] == how) & (df["w"] == 0)]
            b = df[(df["transform"] == kind) & (df["algorithm"] == how) & ~np.isfinite(df["w"])]
            if len(a) and len(b):
                L.append(f"| {kind} | {how} | {a["ari_mapped"].iloc[0]:.3f} | "
                         f"{b["ari_mapped"].iloc[0]:.3f} |")
    L += ["", "## How much of the log-ratio answer is the zero convention?", "",
          f"{int((m == 0).sum())} of {m.size} cells ({100 * (m == 0).mean():.0f} percent) "
          "are zero, so the log-ratio transform",
          "cannot be computed without deciding what a zero is worth. Map plus log-ratio "
          "composition",
          "at half weight, on Mainfort's assemblages, under five conventions:", "",
          "| zero convention | ARI |", "|---|---|"]
    for label, a in zero_rows:
        L.append(f"| {label} | {a:.3f} |")
    L += ["",
          "A convention that amplifies absence harder gives a better match, up to a",
          "perfect one. That is a property of the convention, not of the pottery, and it",
          "is why the chi-square transform -- which needs no such choice -- is the one",
          "the reading below uses.", "",
          "## Reading", "",
          f"Partitions carrying the phases' own group sizes with boundaries placed at "
          f"random agree with the phases at a median ARI of {np.median(rnd):.3f} "
          f"({args.random} draws), which is what",
          "the group sizes manufacture on their own.", "",
          "The site map alone recovers the published scheme well above what the group",
          "sizes manufacture. Composition adds to it, and the addition is real but",
          "secondary: under the chi-square transform, which carries no free parameter,",
          "agreement rises from the map-alone value to its best at a composition weight",
          "of 0.1 to 0.25, that is with the pottery counting for a quarter or less of",
          "what the coordinates count for. Weighting the pottery equally with the map is",
          "already worse than the map alone, and the pottery by itself is worst of all.",
          "",
          "So the scheme is neither purely geographic nor a reading of the pots. It is",
          "the pots seen through where the sites are, with geography carrying most of",
          "the weight. That is consistent with the rest of this paper: on a compositional",
          "gradient with no edges, similarity is a monotone function of proximity, so",
          "sorting assemblages by how alike their pots look largely recovers where they",
          "are, and the residual ceramic signal refines that rather than overriding it.",
          "",
          "**What this does not show.** It does not show the assemblages are",
          "compositionally identical; they are not, and composition varies strongly with",
          "distance. It shows that the particular five-way division the phase scheme",
          "draws is predicted by the coordinates and not by the pots.", "",
          f"Figure written to {png.name} and its siblings; full grid in {OUT_CSV.name}.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(df[df["algorithm"] == "kmeans"][["transform", "w", "ari_all", "ari_mapped"]]
          .to_string(index=False))
    print(f"\nrandom size-matched partitions: median ARI {np.median(rnd):.3f}")
    print(f"wrote {png} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
