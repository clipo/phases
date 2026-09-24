#!/usr/bin/env python3
"""Does letting copying scale with connectivity close the drift shortfall?

analyses/72_excess_locality.py found the shortfall concentrated where the model
is least able to represent the place: per-cluster excess over calibrated drift
correlates with raw connectivity at Spearman -0.80, from 82x in the least
connected cluster to 0.45x in the best connected one. Distance itself does not
predict it (-0.10), nor the effective number of partners (+0.10). What predicts
it is the total influence available before the kernel is normalised.

THE ASSUMPTION UNDER TEST. `copying_weights` returns a row-stochastic matrix,
and `drift_record` applies one mixing rate to every node:

    q = (1 - mixing) * p + mixing * (w @ p)

Because each row of `w` sums to 1, every assemblage takes in exactly `mixing`
worth of other people's pots per step, however isolated it is. The kernel sets
WHO you copy and never HOW MUCH. A site at the edge of the network with a
quarter of the neighbours is modelled as fully participating, so the model
cannot produce the extra divergence isolation would cause, and the difference
lands in the residual we have been calling a shortfall.

THE RELAXATION. Give each node its own mixing rate, proportional to the total
pre-normalisation kernel weight it can reach:

    mixing_i = mixing * conn_i / mean(conn),  clipped to [0, 1]

scaled so the mean over nodes stays at the calibrated value. Nothing else
changes: same geography, same innovation, same populations, same calibration
target. Only the assumption that participation is uniform is dropped.

WHAT WOULD COUNT AS AN ANSWER. If the uniform-connectivity assumption is doing
the work, the shortfall should fall, and fall furthest in the clusters that
carry it now. If the shortfall survives, then isolation as modelled here does
not explain the excess and the residual is about something else.

This is the fifth relaxed assumption, to sit alongside the four in Table 2, and
it is reported the same way: it counts only if it raises differentiation while
the assemblages' own diversity stays matched.

Usage:
    python analyses/73_connectivity_mixing.py [--k 5] [--reps 200]
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

OUT_MD = ROOT / "output" / "findings" / "connectivity_mixing.md"
MODEL = "pooled"
LENGTH = 24.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--reps", type=int, default=200)
    args = ap.parse_args()

    from mls_emergence.transmission.spatial import copying_weights, drift_record
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    loc = importlib.import_module("72_excess_locality")

    counts, coords = mf._load_curated()
    data = rev.load_sets()["basin"]
    m_obs = counts.to_numpy(int)
    xy = coords.to_numpy(float)
    labels = mf._kmeans_labels(xy - xy.mean(0), args.k, seed=7)

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == MODEL)].iloc[0]
    innovation, mixing, n_ind = (float(row["innovation"]), float(row["mixing"]),
                                 int(row["n_ind"]))

    d = data["d"]
    raw = np.exp(-d / LENGTH)
    np.fill_diagonal(raw, 0.0)
    conn = raw.sum(1)
    per_node = np.clip(mixing * conn / conn.mean(), 0.0, 1.0)
    per_node *= mixing / per_node.mean()          # keep the calibrated average
    per_node = np.clip(per_node, 0.0, 1.0)

    w = copying_weights(d, LENGTH, labels=data["labels"], leak=1.0)
    totals = data["m"].sum(1)
    obs_terms = loc.cluster_terms(m_obs, labels)
    obs_fst = float(obs_terms.sum())

    def run(mix) -> tuple[np.ndarray, np.ndarray, list]:
        fst, terms, divs = [], [], []
        for seed in range(args.reps):
            f = drift_record(w, k=data["m"].shape[1], n_ind=n_ind, seed=90000 + seed,
                             innovation=innovation, mixing=mix,
                             target=data["pooled"])
            m = rev.sample_record(f, data["ranks"], totals,
                                  np.random.default_rng(91000 + seed))
            t = loc.cluster_terms(m, labels)
            terms.append(t)
            fst.append(float(t.sum()))
            divs.append(rev.diversity(m))
            if (seed + 1) % 50 == 0:
                print(f"    {seed + 1}/{args.reps}", flush=True)
        return np.array(fst), np.array(terms), divs

    print("  uniform participation (the model as published)")
    fst_u, terms_u, div_u = run(mixing)
    print("  participation scaled by connectivity")
    fst_c, terms_c, div_c = run(per_node)

    obs_div = rev.diversity(m_obs)
    def med(ds, key):
        return float(np.median([x[key] for x in ds if np.isfinite(x.get(key, np.nan))]))

    rows = []
    for i, g in enumerate(np.unique(labels)):
        rows.append(dict(
            cluster=int(g), n=int((labels == g).sum()),
            conn=float(conn[labels == g].mean()),
            mixing_scaled=float(per_node[labels == g].mean()),
            observed=float(obs_terms[i]),
            excess_uniform=float(obs_terms[i] / np.median(terms_u[:, i]))
            if np.median(terms_u[:, i]) > 0 else np.nan,
            excess_scaled=float(obs_terms[i] / np.median(terms_c[:, i]))
            if np.median(terms_c[:, i]) > 0 else np.nan))
    df = pd.DataFrame(rows).sort_values("excess_uniform", ascending=False)

    sf_u = obs_fst / float(np.median(fst_u))
    sf_c = obs_fst / float(np.median(fst_c))
    keys = [k for k in obs_div if isinstance(obs_div[k], (int, float))]
    L = ["# Does connectivity-scaled participation close the shortfall?", "",
         f"Basin phase set, {counts.shape[0]} assemblages, k = {args.k}, "
         f"{args.reps} realizations each, calibrated {MODEL} model "
         f"(innovation {innovation}, mixing {mixing}, {n_ind} learners).",
         "Per-node mixing is proportional to pre-normalisation connectivity, rescaled",
         f"so the mean stays at {mixing}. Range {per_node.min():.4f} to {per_node.max():.4f}.",
         "",
         "| model | observed F_ST | drift median | shortfall |",
         "|---|---|---|---|",
         f"| uniform participation | {obs_fst:.4f} | {float(np.median(fst_u)):.4f} | **{sf_u:.1f}x** |",
         f"| scaled by connectivity | {obs_fst:.4f} | {float(np.median(fst_c)):.4f} | **{sf_c:.1f}x** |",
         "",
         "## Per cluster", "",
         "| cluster | n | raw connectivity | its mixing rate | excess, uniform | excess, scaled |",
         "|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        L.append(f"| {r['cluster']} | {r['n']} | {r['conn']:.2f} | {r['mixing_scaled']:.4f} | "
                 f"{r['excess_uniform']:.1f}x | **{r['excess_scaled']:.1f}x** |")
    L += ["", "## Did the diversity match survive?", "",
          "The relaxation counts only if the assemblages' own diversity stays matched;",
          "a model that reaches the observed differentiation by abandoning the",
          "calibration has fitted one quantity at the other's expense.", "",
          "| summary | observed | uniform | scaled |", "|---|---|---|---|"]
    for key in keys:
        L.append(f"| {key} | {obs_div[key]:.4f} | {med(div_u, key):.4f} | {med(div_c, key):.4f} |")
    L.append("")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n" + df.to_string(index=False))
    print(f"\nshortfall: uniform {sf_u:.1f}x -> connectivity-scaled {sf_c:.1f}x")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
