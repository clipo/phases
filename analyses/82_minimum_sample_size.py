#!/usr/bin/env python3
"""How many decorated sherds does an assemblage need before F_ST means anything?

The membership rule of 2026-09-21 drops assemblages with fewer than 75
decorated sherds. This script is the reason for the number.

The Gini-Simpson cultural F_ST is a plug-in estimator and it is biased upward
by sampling: two assemblages drawn from one and the same profile look
different because each is a finite sample, and the spurious between-assemblage
variance scales roughly as 1/n. The question is at what n that artifact stops
being able to masquerade as the differentiation the basin shows. It is
measured here by drawing pairs of multinomial samples of size n from the
basin's pooled class profile and computing F_ST between them: everything that
comes out is sampling, since the two share one profile by construction.

Reported against the observed between-cluster F_ST at five spatial clusters,
so the reader can see the artifact as a share of the signal. The drift null is
drawn at the observed sherd totals and so carries the same bias; the minimum
protects the per-cluster decomposition and any reading of a single small
cluster, not the calibrated comparison.

Usage:
    python analyses/82_minimum_sample_size.py [--draws 4000]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "minimum_sample_size.md"
SIZES = (25, 50, 75, 100, 150, 200, 300, 500, 1000)
THRESHOLD = 75


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=4000)
    args = ap.parse_args()
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")

    counts, coords = mf._load_curated()
    m = counts.to_numpy(float)
    pooled = m.sum(0) / m.sum()
    xy = coords.to_numpy(float)
    obs = float(rev.fst_by(m.astype(int), mf._kmeans_labels(xy - xy.mean(0), 5, seed=7)))
    rng = np.random.default_rng(82)
    rows = []
    for N in SIZES:
        f = np.array([rev.fst_by(rng.multinomial(N, pooled, size=2), np.array([0, 1]))
                      for _ in range(args.draws)])
        rows.append((N, float(np.median(f)), float(np.percentile(f, 95))))

    L = ["# How many decorated sherds does an assemblage need?", "",
         f"Pairs of multinomial samples of size n from the basin's pooled profile "
         f"(Gini-Simpson diversity {1 - (pooled ** 2).sum():.3f}), {args.draws:,} draws per "
         "size. The two share one profile, so all of the F_ST between them is sampling.", "",
         f"Observed between-cluster F_ST at five clusters, {len(counts)} assemblages: "
         f"**{obs:.4f}**.", "",
         "| decorated sherds per assemblage | spurious F_ST, median | 95th percentile | "
         "median as a share of the observed |", "|---|---|---|---|"]
    for N, med, hi in rows:
        mark = "**" if N == THRESHOLD else ""
        L.append(f"| {mark}{N}{mark} | {med:.4f} | {hi:.4f} | {100 * med / obs:.0f}% |")
    t_med, t_hi = next((med, hi) for N, med, hi in rows if N == THRESHOLD)
    L += ["", f"At {THRESHOLD} sherds the artifact's median is {100 * t_med / obs:.0f} percent of the observed "
              "differentiation and its 95th",
          f"percentile ({t_hi:.4f}) {'falls below' if t_hi < obs else 'does not fall below'} it; "
          "above that the bias shrinks slowly. Mainfort (2003) "
          "kept sites with 800",
          "or more sherds in all, and decorated sherds are 13 percent of his table, so his "
          "rule is about 100",
          "decorated sherds. The analytic expectation for two equal samples from one "
          "profile is close to 1/n.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[6:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
