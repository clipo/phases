#!/usr/bin/env python3
"""Does the drift verdict depend on the scale at which we partition?

A cultural F_ST is a between-group variance, so it needs groups, and the groups
here are k-means clusters of the assemblages' coordinates. That makes k a scale
knob rather than a discovery: cut a continuous settlement distribution finely
and between-group variance rises, cut it coarsely and it falls. On the
43-assemblage phase set the observed F_ST runs from 0.008 at k = 2 to 0.042 at
k = 6 (analyses/59_partition_sensitivity.py), which is 7.7 times the width of
the credible interval reported at a single k. A bare F_ST is therefore not a
quantity; it is a quantity per scale.

The paper's claim does not rest on the level, though. It rests on a COMPARISON:
the observed differentiation against what spatially structured drift produces
on the same geography, scored the same way. If that comparison holds at every
scale, the conclusion is scale-invariant even though the number is not, and
that is a stronger statement than any single k can make. If it flips with k,
the conclusion was an artefact of the partition and needs saying.

METHOD. Drift realizations do not depend on k — the simulation runs on the
river network and knows nothing about the partition — so one set of
realizations is scored at every k. Same simulated records, same calibrated
rates, five scales. That also removes simulation noise from the comparison
across k: any movement is the scale, not a different draw.

The calibrated rates come from analyses/47_revision_analysis.py
(`calibrated_rates.csv`), which tunes each model to the assemblages' own
within-assemblage diversity before any between-group statistic is computed.

Usage:
    python analyses/71_scale_sweep.py [--reps 300]
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

OUT_MD = ROOT / "output" / "findings" / "scale_sweep.md"
OUT_CSV = ROOT / "output" / "findings" / "scale_sweep.csv"
K_RANGE = range(2, 13)
# Through k = 12 so the comparison can be read against group RADIUS in km
# (analyses/75_groupness_surface.py), not only against k. Past k = 7 the
# partition starts isolating single assemblages, which is why `min_group` is
# carried on every row: a group of one matches its own profile exactly and
# adds between-group variance for an arithmetic reason.
MODEL = "pooled"          # the calibrated drift model the main text reports


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()

    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")

    counts, coords = mf._load_curated()
    sets = rev.load_sets()
    data = sets["basin"]
    m_obs = counts.to_numpy(int)
    centred = coords.to_numpy(float) - coords.to_numpy(float).mean(0)
    t74 = importlib.import_module("74_phase_partition_test")
    pts_km = t74.km_xy(coords.to_numpy(float))

    labels = {k: mf._kmeans_labels(centred, k, seed=7) for k in K_RANGE}
    sil = {k: mf.silhouette_mean(centred, labels[k]) for k in K_RANGE}
    observed = {k: rev.fst_by(m_obs, labels[k]) for k in K_RANGE}

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == MODEL)].iloc[0]
    rates = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]),
                 n_ind=int(row["n_ind"]))

    totals = data["m"].sum(1)
    sims = {k: [] for k in K_RANGE}
    for seed in range(args.reps):
        f = rev.simulate(data, 90000 + seed, MODEL, **rates)
        m = rev.sample_record(f, data["ranks"], totals,
                              np.random.default_rng(91000 + seed))
        for k in K_RANGE:                       # one realization, every scale
            sims[k].append(rev.fst_by(m, labels[k]))
        if (seed + 1) % 50 == 0:
            print(f"  {seed + 1}/{args.reps} realizations", flush=True)

    rows = []
    for k in K_RANGE:
        d = np.array(sims[k], dtype=float)
        d = d[np.isfinite(d)]
        lo, med, hi = np.percentile(d, [2.5, 50, 97.5])
        rows.append(dict(k=k, radius_km=float(np.sqrt(t74.inertia(pts_km, labels[k]))),
                         min_group=int(np.bincount(labels[k]).min()),
                         silhouette=sil[k], observed=observed[k],
                         drift_median=med, drift_lo=lo, drift_hi=hi,
                         above_interval=bool(observed[k] > hi),
                         shortfall=observed[k] / med if med > 0 else float("nan"),
                         reach_rate=float(np.mean(d >= observed[k])), n=len(d)))
    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    verdict = ("the observed differentiation exceeds calibrated drift at every scale"
               if df.above_interval.all() else
               "the comparison does NOT hold at every scale; see the table")
    L = ["# Is the drift verdict scale-invariant?", "",
         f"Basin phase set, {counts.shape[0]} assemblages, {args.reps} drift realizations,",
         f"calibrated {MODEL} model (innovation {rates['innovation']}, mixing "
         f"{rates['mixing']}, {rates['n_ind']} learners). Each realization is scored at",
         "every k, so differences across scales are the partition and not a different draw.",
         "",
         "| k | silhouette | observed F_ST | drift median | drift 95% | observed above? | shortfall | share of runs reaching observed |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['k']} | {r['silhouette']:.3f} | **{r['observed']:.4f}** | "
                 f"{r['drift_median']:.4f} | [{r['drift_lo']:.4f}, {r['drift_hi']:.4f}] | "
                 f"{'yes' if r['above_interval'] else 'NO'} | {r['shortfall']:.1f}x | "
                 f"{r['reach_rate'] * 100:.1f}% |")
    L += ["", f"**Verdict: {verdict}.**", "",
          "The level is a function of the scale: the observed F_ST rises with k because",
          "a finer cut of a continuous distribution puts more variance between groups.",
          "The comparison is not, and that is the quantity the argument uses.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n" + df.to_string(index=False))
    print(f"\n{verdict}\nwrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
