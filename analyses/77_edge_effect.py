#!/usr/bin/env python3
"""Is the drift shortfall an edge effect of simulating a closed basin?

analyses/72_excess_locality.py found the shortfall concentrated in two small
clusters, and they are not arbitrary ones. Dundee, Parchman, Salomon and West
Mounds are the Parchman phase, on the southern edge across the Mississippi;
Carson Lake, Notgrass and Upper Nodena are the Nodena phase, on the northern
edge. Both are fragments of phases whose bulk lies outside the 43-assemblage
set. The interior clusters carry little excess (the Parkin sites 3.9-fold) or
none (the Walls sites, below drift).

THE ASSUMPTION UNDER TEST. The drift model treats the 43 assemblages as the
whole world. An edge node has nobody beyond it, so all its copying points
inward and the model pulls it toward the basin's interior. A real edge
settlement also copied from neighbours outside the set, which would hold it
apart from the interior. A closed model therefore UNDERSTATES how distinct
edge clusters should be under drift alone, and the difference lands in the
residual.

THE TEST. Run the same calibrated drift on the larger 55-assemblage network,
which adds twelve assemblages outside the phase set, then score only the 43
basin assemblages, with the same clusters, sherd totals and innovation profile.
Nothing changes but whether the edge has neighbours. The same seeds are used in
both arms so the difference is the geometry, not the draw.

WHAT IT CAN AND CANNOT SHOW. The twelve added assemblages are the Tennessee
sites east and north-east of the basin. They open the edge beside the Nodena
and Walls clusters. NOTHING in the curated matrix lies south of the Parchman
cluster, so this test cannot open that edge, and a null result for Parchman is
uninformative rather than negative. The script prints, for each cluster, how
many outside assemblages fall within one interaction length of it, so the
reader can see which edges the test reaches.

Usage:
    python analyses/77_edge_effect.py [--k 5] [--reps 300]
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

OUT_MD = ROOT / "output" / "findings" / "edge_effect.md"
MODEL = "pooled"
LENGTH = 24.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()

    from mls_emergence.transmission.spatial import copying_weights, drift_record
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    loc = importlib.import_module("72_excess_locality")

    sets = rev.load_sets()
    basin, wide = sets["basin"], sets["valley"]
    bn = [str(x) for x in basin["names"]]
    wn = [str(x) for x in wide["names"]]
    missing = [n for n in bn if n not in wn]
    if missing:
        raise ValueError(f"basin assemblages absent from the wide set: {missing}")
    idx = np.array([wn.index(n) for n in bn])          # basin rows inside the wide set
    outside = np.array([i for i in range(len(wn)) if i not in set(idx)], dtype=int)
    if len(outside) == 0:
        # Since the matrix was rebuilt from PFGData and Lipo's compilation
        # (2026-09-21) every assemblage in it lies inside the basin: the ones
        # outside were Mainfort's alone, and his table is now a separate
        # comparison. There is nothing to open the edge WITH, so the test is
        # reported as not runnable rather than run on an empty set.
        OUT_MD.write_text("\n".join([
            "# Is the drift shortfall an edge effect of a closed basin?", "",
            f"NOT RUNNABLE on the primary matrix. All {len(wn)} assemblages in it lie inside "
            "the basin, so there are no",
            "outside neighbours to add. The assemblages beyond the basin's edge are known "
            "only from Mainfort's (2003)",
            "table, which this project keeps as a separate comparison and does not mix into "
            "the matrix the drift model",
            "is calibrated on. The question stays open. Answering it would need counts for "
            "sites beyond the Nodena and",
            "Parchman edges from a source independent of his, or a ruling that his counts "
            "may stand in as boundary",
            "conditions only (never scored).", "",
            "An earlier run, on the matrix that summed the two sources, found a partial "
            "effect at the one edge it could",
            "reach (Nodena, 15.8-fold to 10.3-fold). That number described the defective "
            "matrix and is superseded.", ""]), encoding="utf-8")
        print(f"edge-effect test not runnable: no assemblages outside the basin; wrote {OUT_MD}")
        return 0

    counts, coords = mf._load_curated()
    if [str(i) for i in counts.index] != bn:
        raise ValueError("curated order differs from the simulation set's order")
    m_obs = counts.to_numpy(int)
    xy = coords.to_numpy(float)
    labels = mf._kmeans_labels(xy - xy.mean(0), args.k, seed=7)

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == MODEL)].iloc[0]
    innovation, mixing, n_ind = float(row.innovation), float(row.mixing), int(row.n_ind)

    # Which edges does the wide set actually open?
    d_wide = wide["d"]
    reach = {}
    for g in np.unique(labels):
        members = idx[labels == g]
        near = (d_wide[np.ix_(members, outside)] <= LENGTH).any(0).sum()
        reach[int(g)] = int(near)

    w_closed = copying_weights(basin["d"], LENGTH, labels=basin["labels"], leak=1.0)
    w_open = copying_weights(d_wide, LENGTH, labels=wide["labels"], leak=1.0)
    k_cls = basin["m"].shape[1]
    totals_b = basin["m"].sum(1)
    totals_w = wide["m"].sum(1)
    target = basin["pooled"]                           # one innovation profile, both arms

    # Third arm. With one innovation profile everywhere, outside neighbours are
    # more of the same and cannot hold an edge cluster apart, so the "open" arm
    # alone cannot show an edge effect even if one operates. Here the outside
    # nodes innovate from THEIR OWN observed pooled profile. That is not
    # circular: the outside assemblages are never scored, so this conditions on
    # what the world beyond the basin looked like and asks how distinct the
    # basin's edge should then be under drift.
    out_pool = wide["m"][outside].sum(0).astype(float)
    out_pool = out_pool / out_pool.sum()
    target_open = np.tile(np.asarray(target, float), (len(wn), 1))
    target_open[outside] = out_pool
    l1 = float(np.abs(out_pool - np.asarray(target, float)).sum())

    res = {"closed": ([], [], []), "open": ([], [], []), "open, distinct outside": ([], [], [])}
    for seed in range(args.reps):
        f = drift_record(w_closed, k=k_cls, n_ind=n_ind, seed=97000 + seed,
                         innovation=innovation, mixing=mixing, target=target)
        m = rev.sample_record(f, basin["ranks"], totals_b,
                              np.random.default_rng(98000 + seed))
        t = loc.cluster_terms(m, labels)
        res["closed"][0].append(t); res["closed"][1].append(t.sum())
        res["closed"][2].append(rev.diversity(m))

        f = drift_record(w_open, k=k_cls, n_ind=n_ind, seed=97000 + seed,
                         innovation=innovation, mixing=mixing, target=target)
        m = rev.sample_record(f, wide["ranks"], totals_w,
                              np.random.default_rng(98000 + seed))[idx]
        t = loc.cluster_terms(m, labels)
        res["open"][0].append(t); res["open"][1].append(t.sum())
        res["open"][2].append(rev.diversity(m))
        f = drift_record(w_open, k=k_cls, n_ind=n_ind, seed=97000 + seed,
                         innovation=innovation, mixing=mixing, target=target_open)
        m = rev.sample_record(f, wide["ranks"], totals_w,
                              np.random.default_rng(98000 + seed))[idx]
        t = loc.cluster_terms(m, labels)
        a3 = "open, distinct outside"
        res[a3][0].append(t); res[a3][1].append(t.sum()); res[a3][2].append(rev.diversity(m))
        if (seed + 1) % 50 == 0:
            print(f"  {seed + 1}/{args.reps}", flush=True)

    obs_terms = loc.cluster_terms(m_obs, labels)
    obs_fst = float(obs_terms.sum())
    obs_div = rev.diversity(m_obs)
    med = {a: np.median(np.array(res[a][0]), axis=0) for a in res}
    fst = {a: float(np.median(res[a][1])) for a in res}
    div = {a: {k: float(np.median([x[k] for x in res[a][2]])) for k in ("hs", "rich", "ht")}
           for a in res}

    L = ["# Is the drift shortfall an edge effect of a closed basin?", "",
         f"Basin phase set, {len(bn)} assemblages scored; the open arm simulates "
         f"{len(wn)} ({len(outside)} outside the set).",
         f"k = {args.k}, {args.reps} realizations per arm on shared seeds, calibrated "
         f"{MODEL} cell ({n_ind} learners, innovation {innovation}, mixing {mixing}), "
         f"interaction length {LENGTH:g} river-km.", "",
         "| arm | observed F_ST | drift median | shortfall | H_S | richness | H_T |",
         "|---|---|---|---|---|---|---|",
         f"| observed | {obs_fst:.4f} | | | {obs_div['hs']:.3f} | {obs_div['rich']:.2f} | "
         f"{obs_div['ht']:.3f} |"]
    for a in res:
        L.append(f"| {a} | {obs_fst:.4f} | {fst[a]:.4f} | **{obs_fst / fst[a]:.1f}x** | "
                 f"{div[a]['hs']:.3f} | {div[a]['rich']:.2f} | {div[a]['ht']:.3f} |")
    L += ["", "## Per cluster", "",
          "| cluster | n | outside assemblages within one interaction length | "
          "observed term | excess, closed | excess, open | excess, open with distinct "
          "outside |", "|---|---|---|---|---|---|---|"]
    names_by = {int(g): [n for n, l in zip(bn, labels) if l == g] for g in np.unique(labels)}
    for i, g in enumerate(np.unique(labels)):
        ec = obs_terms[i] / med["closed"][i] if med["closed"][i] > 0 else np.nan
        eo = obs_terms[i] / med["open"][i] if med["open"][i] > 0 else np.nan
        m3 = med["open, distinct outside"][i]
        e3 = obs_terms[i] / m3 if m3 > 0 else np.nan
        L.append(f"| {int(g)} ({', '.join(names_by[int(g)][:3])}...) | "
                 f"{int((labels == g).sum())} | {reach[int(g)]} | {obs_terms[i]:.4f} | "
                 f"{ec:.1f}x | {eo:.1f}x | **{e3:.1f}x** |")
    L += ["", f"The outside assemblages' pooled profile differs from the basin's by an "
              f"L1 distance of {l1:.3f}", "(0 identical, 2 disjoint).",
          "", "A cluster with no outside assemblage within one interaction length has "
              "an edge this test",
          "cannot open; its row says nothing about whether an edge effect operates "
          "there.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[5:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
