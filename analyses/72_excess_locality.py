#!/usr/bin/env python3
"""Is the drift shortfall concentrated in particular clusters, or spread?

analyses/71_scale_sweep.py shows the observed cultural F_ST exceeding calibrated
drift at every scale the test can discriminate, and the shortfall growing with
grain: 2.9x at k = 2, 8.4x at k = 6. Two accounts fit that shape and they say
different things.

  CLOSURE.  Some part of the basin really did narrow its copying. Then the
            excess should be LOCAL: one or two clusters carry it, and removing
            them should bring the rest close to what drift produces.
  GRAIN.    The simulated field is smooth and the real settlement pattern is
            lumpy, so a finer cut isolates neighbours that resemble each other
            more than a distance-decay model predicts. Then the excess should
            be EVEN: every cluster contributes in proportion to its size, and
            dropping any one leaves the picture unchanged.

THE MEASUREMENT. Cultural F_ST is a ratio of between-group to total diversity,
so each cluster contributes a term to the numerator: how far its own class
profile sits from the pooled profile, weighted by its sherds. This script
computes that contribution per cluster for the observed data and for each drift
realization, then asks whether the observed contributions are uniformly above
the drift ones or concentrated in a few clusters.

Two readouts, because one number would hide the shape:

  per-cluster excess   observed contribution divided by the drift median for
                       the same cluster. Closure predicts a few large values;
                       grain predicts a flat profile near the overall shortfall.
  leave-one-out        the overall shortfall recomputed with each cluster
                       dropped. Closure predicts one cluster whose removal
                       collapses it; grain predicts little movement whichever
                       cluster goes.

The same drift realizations are used throughout, so differences between
clusters are the data, not separate draws.

Usage:
    python analyses/72_excess_locality.py [--k 5] [--reps 300]
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

OUT_MD = ROOT / "output" / "findings" / "excess_locality.md"
OUT_CSV = ROOT / "output" / "findings" / "excess_locality.csv"
MODEL = "pooled"


def cluster_terms(m: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Each cluster's contribution to the between-group term of cultural F_ST.

    For Gini-Simpson diversity the between-group term decomposes exactly into
    non-negative per-cluster pieces:

        D_total - D_within = sum_g w_g * sum_i (p_gi - p_i)^2

    with w_g the cluster's share of sherds, p_gi its class proportions and p_i
    the pooled ones. Dividing by D_total gives terms that sum to F_ST itself,
    so the parts add to the whole and none can be negative.

    An earlier version of this function used w_g * (D_total - D_g), which sums
    to the same total but goes NEGATIVE for any cluster more diverse than the
    pooled assemblage. Ratios of those terms to their drift counterparts are
    meaningless, and the script's own verdict was computed from them. The
    identity above is the one to use.
    """
    tot = m.sum()
    if tot == 0:
        return np.zeros(len(np.unique(labels)))
    pooled = m.sum(0) / tot
    div_pooled = 1.0 - float((pooled ** 2).sum())
    out = []
    for g in np.unique(labels):
        sel = labels == g
        gm = m[sel].sum(0)
        n = gm.sum()
        if n == 0:
            out.append(0.0)
            continue
        p = gm / n
        out.append((n / tot) * float(((p - pooled) ** 2).sum()))
    return np.array(out) / div_pooled if div_pooled > 0 else np.zeros(len(out))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()

    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")

    counts, coords = mf._load_curated()
    data = rev.load_sets()["basin"]
    m_obs = counts.to_numpy(int)
    xy = coords.to_numpy(float)
    labels = mf._kmeans_labels(xy - xy.mean(0), args.k, seed=7)
    names = list(counts.index)

    obs_terms = cluster_terms(m_obs, labels)
    obs_fst = float(obs_terms.sum())

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == MODEL)].iloc[0]
    rates = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]),
                 n_ind=int(row["n_ind"]))
    totals = data["m"].sum(1)

    sim_terms, sim_fst = [], []
    for seed in range(args.reps):
        f = rev.simulate(data, 95000 + seed, MODEL, **rates)
        m = rev.sample_record(f, data["ranks"], totals, np.random.default_rng(96000 + seed))
        sim_terms.append(cluster_terms(m, labels))
        sim_fst.append(float(cluster_terms(m, labels).sum()))
        if (seed + 1) % 50 == 0:
            print(f"  {seed + 1}/{args.reps}", flush=True)
    sim_terms = np.array(sim_terms)
    sim_fst = np.array(sim_fst)

    rows = []
    for i, g in enumerate(np.unique(labels)):
        members = [n for n, l in zip(names, labels) if l == g]
        drift_med = float(np.median(sim_terms[:, i]))
        rows.append(dict(
            cluster=int(g), n_assemblages=len(members),
            sherds=int(m_obs[labels == g].sum()),
            observed_term=float(obs_terms[i]),
            share_of_observed=float(obs_terms[i] / obs_fst) if obs_fst else np.nan,
            drift_median_term=drift_med,
            excess=float(obs_terms[i] / drift_med) if drift_med > 0 else np.nan,
            members=", ".join(members)))
    df = pd.DataFrame(rows).sort_values("excess", ascending=False)

    # Leave one cluster out, observed and drift alike, and recompute the shortfall.
    loo = []
    for i, g in enumerate(np.unique(labels)):
        keep = labels != g
        o = float(cluster_terms(m_obs[keep], labels[keep]).sum())
        d = np.array([float(cluster_terms(mm, labels[keep]).sum())
                      for mm in [rev.sample_record(
                          rev.simulate(data, 95000 + s, MODEL, **rates),
                          data["ranks"], totals, np.random.default_rng(96000 + s))[keep]
                          for s in range(min(args.reps, 120))]])
        loo.append(dict(dropped=int(g), observed=o, drift_median=float(np.median(d)),
                        shortfall=o / float(np.median(d)) if np.median(d) > 0 else np.nan))
    loo_df = pd.DataFrame(loo)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    overall = obs_fst / float(np.median(sim_fst))
    ex = df.excess.to_numpy(dtype=float)
    spread = np.nanmax(ex) / np.nanmin(ex) if np.nanmin(ex) > 0 else np.nan
    verdict = ("CONCENTRATED: one or two clusters carry the excess"
               if spread >= 3 else
               "EVEN: every cluster is above drift by a similar factor, which is what "
               "grain predicts and closure does not")
    L = ["# Where does the drift shortfall live?", "",
         f"Basin phase set, {counts.shape[0]} assemblages, k = {args.k}, "
         f"{args.reps} drift realizations, calibrated {MODEL} model.",
         f"Observed F_ST {obs_fst:.4f} against a drift median of "
         f"{float(np.median(sim_fst)):.4f}, an overall shortfall of {overall:.1f}x.", "",
         "## Per-cluster contribution", "",
         "| cluster | n | sherds | observed term | share of observed F_ST | drift median term | excess |",
         "|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        L.append(f"| {r['cluster']} | {r['n_assemblages']} | {r['sherds']:,} | "
                 f"{r['observed_term']:.4f} | {r['share_of_observed'] * 100:.0f}% | "
                 f"{r['drift_median_term']:.4f} | **{r['excess']:.1f}x** |")
    L += ["", f"Excess ranges from {np.nanmin(ex):.1f}x to {np.nanmax(ex):.1f}x, "
              f"a ratio of {spread:.1f}.", "",
          f"**{verdict}.**", "",
          "## Leave one cluster out", "",
          "| dropped | observed F_ST | drift median | shortfall |",
          "|---|---|---|---|"]
    for _, r in loo_df.iterrows():
        L.append(f"| {r['dropped']} | {r['observed']:.4f} | {r['drift_median']:.4f} | "
                 f"{r['shortfall']:.1f}x |")
    L += ["", "If the excess were one cluster's closure, dropping that cluster would "
              "collapse the shortfall toward 1. Read the column for that.", "",
          "## Cluster membership", ""]
    for _, r in df.iterrows():
        L.append(f"- **cluster {r['cluster']}** ({r['n_assemblages']}): {r['members']}")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n" + df[["cluster", "n_assemblages", "sherds", "observed_term",
                     "drift_median_term", "excess"]].to_string(index=False))
    print("\nleave-one-out:\n" + loo_df.to_string(index=False))
    print(f"\noverall shortfall {overall:.1f}x; {verdict}\nwrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
