"""85_excess_robustness.py - is the finer-scale excess over drift a sampling effect?

Analysis 72 decomposed the between-cluster cultural F_ST at five spatial
clusters into per-cluster terms and found the excess over calibrated drift
concentrated in the two most thinly sampled clusters (765 and 918 sherds).
Three checks here, all at the same five clusters and calibrated cell:

(a) SAMPLING ALONE. Every assemblage is redrawn from the basin's pooled class
    profile at its own sherd count, 2,000 times, and each cluster's term is
    computed. This is what the term would be with no spatial structure and no
    drift at all, only finite samples.
(b) UNCERTAINTY IN THE OBSERVED COUNTS. Each assemblage's class proportions
    are drawn from their posterior, Dirichlet(counts + 1/2) (the Jeffreys
    prior, which leans toward even profiles and so slightly toward LARGER
    terms than the counts; named per rule 20), 2,000 times, and each
    cluster's term recomputed with the observed sherd weights. The ratio of
    each draw to the drift median gives a posterior for the excess.
(c) ONE ASSEMBLAGE AT A TIME. Within each cluster, every assemblage is dropped
    from the observed record and from every drift run alike, and the
    cluster's excess recomputed, to show whether one collection carries it.

Drift runs reuse analysis 72's seed families (95000 simulate, 96000 sample)
so (c) is comparable with the published decomposition.

Outputs: output/findings/excess_robustness.md
Usage: .venv/bin/python analyses/85_excess_robustness.py [--reps 300] [--draws 2000]
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
OUT_MD = ROOT / "output" / "findings" / "excess_robustness.md"
K = 5


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--draws", type=int, default=2000)
    args = ap.parse_args()

    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    a72 = importlib.import_module("72_excess_locality")
    terms = a72.cluster_terms

    counts, coords = mf._load_curated()
    data = rev.load_sets()["basin"]
    m_obs = counts.to_numpy(int)
    names = list(counts.index)
    xy = coords.to_numpy(float)
    labels = mf._kmeans_labels(xy - xy.mean(0), K, seed=7)
    groups = np.unique(labels)
    n_a = m_obs.sum(1)
    pooled = m_obs.sum(0) / m_obs.sum()

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == a72.MODEL)].iloc[0]
    rates = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))

    obs_terms = terms(m_obs, labels)
    sims = [rev.sample_record(rev.simulate(data, 95000 + s, a72.MODEL, **rates),
                              data["ranks"], n_a, np.random.default_rng(96000 + s))
            for s in range(args.reps)]
    sim_terms = np.array([terms(mm, labels) for mm in sims])
    drift_med = np.median(sim_terms, axis=0)

    rng = np.random.default_rng(85001)
    samp = np.array([terms(np.array([rng.multinomial(n, pooled) for n in n_a]), labels)
                     for _ in range(args.draws)])
    rng = np.random.default_rng(85002)
    post = []
    for _ in range(args.draws):
        p = np.array([rng.dirichlet(row + 0.5) for row in m_obs])
        post.append(terms(p * n_a[:, None], labels))
    post = np.array(post)

    L = ["# Is the finer-scale excess a sampling effect?", "",
         f"Produced by `analyses/85_excess_robustness.py`. Basin, {len(names)} assemblages, k = {K} "
         f"spatial clusters (the division of analysis 72), calibrated pooled-profile cell "
         f"({rates['n_ind']} learners, innovation {rates['innovation']}, mixing {rates['mixing']}), "
         f"{args.reps} drift runs, {args.draws} draws for (a) and (b).", "",
         "| cluster | assemblages | sherds | observed term | (a) sampling alone, median [95%] | drift median | "
         "(b) excess over drift, posterior median [95%] | P(excess > 1) |", "|---|---|---|---|---|---|---|---|"]
    order = np.argsort(-(obs_terms / np.where(drift_med > 0, drift_med, np.nan)))
    for i in order:
        g = groups[i]; sel = labels == g
        ratio = post[:, i] / drift_med[i]
        L.append(f"| {g} ({', '.join(n for n, l in zip(names, labels) if l == g)}) | {int(sel.sum())} | "
                 f"{int(n_a[sel].sum()):,} | {obs_terms[i]:.4f} | {np.median(samp[:, i]):.5f} "
                 f"[{np.percentile(samp[:, i], 2.5):.5f}, {np.percentile(samp[:, i], 97.5):.5f}] | "
                 f"{drift_med[i]:.4f} | {np.median(ratio):.1f} [{np.percentile(ratio, 2.5):.1f}, "
                 f"{np.percentile(ratio, 97.5):.1f}] | {float((ratio > 1).mean()):.2f} |")
    L += ["", "## (c) Dropping one assemblage at a time", "",
          "| cluster | assemblage dropped | observed term | drift median | excess |", "|---|---|---|---|---|"]
    for i in order:
        g = groups[i]
        members = np.flatnonzero(labels == g)
        if len(members) < 3:
            continue
        for j in members:
            keep = np.ones(len(names), bool); keep[j] = False
            lab = labels[keep]
            gi = list(np.unique(lab)).index(g)
            o = terms(m_obs[keep], lab)[gi]
            dm = float(np.median([terms(mm[keep], lab)[gi] for mm in sims]))
            L.append(f"| {g} | {names[j]} | {o:.4f} | {dm:.4f} | {o / dm:.1f} |")
    L += ["", "## Reading", "",
          "(a) says how large each cluster's term would be from finite samples alone. (b) carries the "
          "uncertainty in the observed counts into the excess; a posterior that stays above 1 says the "
          "excess is not an accident of which sherds were collected. (c) says whether the excess belongs "
          "to the cluster or to one collection in it.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
