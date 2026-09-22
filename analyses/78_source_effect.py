#!/usr/bin/env python3
"""Does the drift shortfall depend on which tally an assemblage's counts come from?

Every assemblage in the analysis matrix takes its counts from ONE source
(scripts/build_analysis_matrix.py): PFGData.xlsx, the survey's counts alone, or Lipo's (2001)
compilation where his tab holds the assemblage. (An earlier version of this
docstring and its labels said "Mainfort's table" for the first kind; the matrix
has held no Mainfort rows since 2026-09-21.) The two are
not independent collections, since Mainfort's numbers are Lipo's, which are
the Phillips-Ford-Griffin counts plus additions. What differs between them is
completeness, and with it sample size. If the rows from one source sit in
particular places, that difference can read as spatial differentiation that no
drift model would reproduce.

Two questions, and they need different measurements.

1. IS SOURCE CONFOUNDED WITH PLACE? A cross-tabulation of source against the
   five spatial clusters.

2. DO THE SOURCES DIFFER ONCE PLACE IS HELD? Within each cluster holding at
   least two rows of each kind, the cultural F_ST between the two kinds of row,
   compared with the same statistic for random splits of that cluster into
   groups of the same sizes, which describes what a split of that size produces
   from within-cluster variation alone (rule 18: descriptive, not a test).

Then the consequence that matters: the between-cluster F_ST and its excess over
calibrated drift, recomputed on each kind of row separately, on the clusters
both kinds cover.

HISTORY. Written on 2026-09-21 against the old matrix, which SUMMED the two
sources at 19 assemblages and doubled Parkin Punctated in one of them. It then
contrasted "Mainfort-only" rows with "merged" rows and found the shortfall
12.8-fold in the first and 4.3-fold in the second. Those numbers described the
defective matrix and are superseded.

Usage:
    python analyses/78_source_effect.py [--splits 2000] [--reps 200]
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

OUT_MD = ROOT / "output" / "findings" / "source_effect.md"
SOURCES = ROOT / "data" / "processed" / "analysis_matrix_provenance.csv"
MODEL = "pooled"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--splits", type=int, default=2000)
    ap.add_argument("--reps", type=int, default=200)
    args = ap.parse_args()

    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m = counts.to_numpy(int)
    xy = coords.to_numpy(float)
    labels = mf._kmeans_labels(xy - xy.mean(0), args.k, seed=7)

    src = pd.read_csv(SOURCES).set_index("assemblage")
    absent = [n for n in names if n not in src.index]
    if absent:
        raise ValueError(f"no source recorded for: {absent}; rerun "
                         "scripts/build_analysis_matrix.py")
    # Since 2026-09-21 no row is a sum. `merged` keeps its name for the code
    # below but now means "this row was taken from Lipo's compilation", the
    # larger tally at that assemblage, as against Mainfort's table.
    merged = np.array([str(src.loc[n, "source_used"]).startswith("Lipo") for n in names])

    L = ["# Is the drift shortfall an artefact of merging data sources?", "",
         f"Basin phase set, {len(names)} assemblages: {int((~merged).sum())} rows carrying the "
         f"survey's counts alone (PFGData.xlsx) and {int(merged.sum())} carrying Lipo's (2001) compilation, which holds "
         f"the Phillips-Ford-Griffin counts plus his 1996-97 field collections. The two are "
         f"tallies of the same material, so each assemblage uses one of them, the larger, and "
         f"none is a sum (scripts/build_analysis_matrix.py).", "",
         "## 1. Is source confounded with place?", "",
         "| cluster | assemblages | survey alone | Lipo compilation | sherds, survey alone | "
         "sherds, Lipo compilation |", "|---|---|---|---|---|---|"]
    for g in np.unique(labels):
        s = labels == g
        L.append(f"| {int(g)} ({', '.join([n for n, x in zip(names, s) if x][:3])}...) | "
                 f"{int(s.sum())} | {int((s & ~merged).sum())} | {int((s & merged).sum())} | "
                 f"{int(m[s & ~merged].sum()):,} | {int(m[s & merged].sum()):,} |")

    L += ["", "## 2. Do the sources differ once place is held?", "",
          "| cluster | F_ST, survey-alone rows against compilation rows | random splits of the same "
          "sizes, median (5th-95th) | position |", "|---|---|---|---|"]
    rng = np.random.default_rng(78)
    for g in np.unique(labels):
        s = np.where(labels == g)[0]
        a = merged[s]
        if a.sum() < 2 or (~a).sum() < 2:
            L.append(f"| {int(g)} | fewer than two rows of one kind | | not testable |")
            continue
        obs = float(rev.fst_by(m[s], a.astype(int)))
        ref = np.array([rev.fst_by(m[s], rng.permutation(a).astype(int))
                        for _ in range(args.splits)])
        L.append(f"| {int(g)} | {obs:.4f} | {np.median(ref):.4f} "
                 f"({np.percentile(ref, 5):.4f}-{np.percentile(ref, 95):.4f}) | "
                 f"above {100 * np.mean(ref < obs):.0f} percent of {args.splits:,} |")

    # 3. The excess over drift, on each source subset separately.
    data = rev.load_sets()["basin"]
    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == MODEL)].iloc[0]
    rates = dict(innovation=float(row.innovation), mixing=float(row.mixing),
                 n_ind=int(row.n_ind))
    totals = data["m"].sum(1)
    # The merged rows reach only the clusters Lipo's compilation covers, so the
    # raw Mainfort-against-merged contrast compares five clusters with four. The
    # "shared clusters" rows hold the clusters fixed, which is the comparison
    # to read.
    shared = np.isin(labels, [g for g in np.unique(labels)
                              if (merged & (labels == g)).any()
                              and (~merged & (labels == g)).sum() >= 1
                              and (merged & (labels == g)).sum() >= 1])
    shared &= np.isin(labels, np.unique(labels[merged]))
    subsets = {"all 43": np.ones(len(names), bool), "rows carrying the survey's counts alone": ~merged,
               "rows carrying Lipo's compilation": merged,
               "rows carrying the survey's counts alone, shared clusters only": ~merged & shared,
               "rows carrying Lipo's compilation, shared clusters only": merged & shared}
    sims = {k: [] for k in subsets}
    for seed in range(args.reps):
        f = rev.simulate(data, 99000 + seed, MODEL, **rates)
        mm = rev.sample_record(f, data["ranks"], totals, np.random.default_rng(99500 + seed))
        for k, sel in subsets.items():
            if len(np.unique(labels[sel])) >= 2:
                sims[k].append(rev.fst_by(mm[sel], labels[sel]))
        if (seed + 1) % 50 == 0:
            print(f"  {seed + 1}/{args.reps}", flush=True)
    L += ["", "## 3. Does the excess over drift survive in a single source?", "",
          f"Between-cluster F_ST at k = {args.k}, the same clusters throughout, against "
          f"{args.reps} calibrated drift realizations scored on the same subset.", "",
          "| subset | assemblages | clusters represented | median sherds per assemblage | "
          "observed F_ST | drift median (95 percent) | shortfall |",
          "|---|---|---|---|---|---|---|"]
    for k, sel in subsets.items():
        d = np.array(sims[k], float)
        obs = float(rev.fst_by(m[sel], labels[sel]))
        L.append(f"| {k} | {int(sel.sum())} | {len(np.unique(labels[sel]))} | "
                 f"{int(np.median(m[sel].sum(1)))} | {obs:.4f} | "
                 f"{np.median(d):.4f} ({np.percentile(d, 2.5):.4f}-"
                 f"{np.percentile(d, 97.5):.4f}) | **{obs / np.median(d):.1f}x** |")
    L += ["", "**Source and sample size are not separable here.** The rows carrying the survey's counts alone "
          "are also the small",
          "ones, so a larger shortfall among them may be an analyst or collection-regime "
          "effect, or",
          "overdispersion in small surface collections, or both. Either is a property of "
          "the record",
          "rather than of past interaction, which is the distinction that matters for the "
          "residual.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[4:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
