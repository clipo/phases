"""Why the simulated two-cluster F_ST is below the observed value: span variants at the operating point.

Latent-group probe, pilot step 5 follow-up. Exploratory design validation
(rule 1). The restricted screen found the operating point's two-cluster
F_ST median at 0.0009 against the observed 0.0062 (no cell above 0.0032).
This script measures the same summaries at the operating point under four
sampling variants, 20 seeds each, to see what the paired-span time averaging
and the end jitter do to the between-cluster contrast: (a) T1 as in the
screen; (b) T0 (same lengths, common end); (c) a single stored step, no
averaging; (d) the old eight-snapshot window (slots 20 steps apart, as in
Carl's runs) at a common end. Same functions as the screen.

Usage: .venv/bin/python analyses/probe_latent_groups/resemblance_variants.py
Output: output/findings/resemblance_variants.md
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses")); sys.path.insert(0, str(Path(__file__).resolve().parent))

from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.probe import resemblance as rs
from mls_emergence.probe import spans
from mls_emergence.transmission.spatial import copying_weights, drift_record

OUTD = ROOT / "output" / "probe_latent_groups"
FIND = ROOT / "output" / "findings" / "resemblance_variants.md"
OP = {"n_ind": 2000, "innovation": 0.001, "mixing": 0.02}
N_SEEDS = 20


def main() -> int:
    m47 = importlib.import_module("47_revision_analysis")
    from a0_predictive_summaries import spatial_clusters
    counts, coords = load_curated(); xy = coords.to_numpy(float)
    tsj = json.loads((OUTD / "time_scales.json").read_text()); L_d, S, burnin = int(tsj["L_d"]), int(tsj["probe_stored_length"]), int(tsj["probe_burnin"])
    N, K = counts.shape; totals = counts.to_numpy(int).sum(1)
    pooled = counts.to_numpy(float).sum(0); pooled /= pooled.sum()
    labels = spatial_clusters(xy); w = copying_weights(geo_km(xy), length=24.0)
    obs = rs.summaries(counts.to_numpy(int), labels, m47.diversity)
    rows = []
    for j in range(N_SEEDS):
        ss = np.random.SeedSequence(2026101660 + j); kids = ss.spawn(4)
        rec = drift_record(w, k=K, steps=S, n_records=S, burnin=burnin, seed=int(np.random.default_rng(kids[0]).integers(0, 2**31 - 1)),
                           initial="uniform", target=pooled, **OP)
        lengths = spans.draw_lengths(np.random.default_rng(kids[1]), N, L_d)
        variants = {
            "T1 (screen): lengths 1..L_d, ends over last L_d": (spans.ends_t1(np.random.default_rng(kids[2]), N, S, L_d), lengths),
            "T0: same lengths, common end": (spans.ends_t0(N, S), lengths),
            "single step, no averaging": (spans.ends_t0(N, S), np.ones(N, int)),
        }
        sub = np.random.default_rng(kids[3])
        for name, (ends, w_) in variants.items():
            y = spans.sample_spans(rec, ends, w_, totals, sub)
            rows.append({"variant": name, "seed": j, **rs.summaries(y, labels, m47.diversity)})
        # Old eight-snapshot window: 8 slots 20 steps apart ending at the last step, common end.
        slots = S - 1 - 20 * np.arange(8)
        p = rec[slots].mean(0); p /= p.sum(1, keepdims=True)
        y = np.array([sub.multinomial(int(n), q) for n, q in zip(totals, p)])
        rows.append({"variant": "eight snapshots 20 steps apart (old window), common end", "seed": j, **rs.summaries(y, labels, m47.diversity)})
        # Spatial field alone: the stored frequencies at the last step without sherd sampling noise (plug-in on expected counts).
        p = rec[-1]; y = np.rint(p * totals[:, None]).astype(int)
        rows.append({"variant": "last step, expected counts (no multinomial noise)", "seed": j, **rs.summaries(y, labels, m47.diversity)})
    df = pd.DataFrame(rows)
    g = df.groupby("variant", sort=False).agg(hs=("hs", "mean"), rich=("rich", "mean"), ht=("ht", "mean"),
                                              fst_median=("fst_clusters", "median"), fst_q10=("fst_clusters", lambda v: v.quantile(.1)),
                                              fst_q90=("fst_clusters", lambda v: v.quantile(.9)))
    L = ["# Span variants at the operating point: where the between-cluster contrast goes (pilot step 5)", "",
         ("Produced by `analyses/probe_latent_groups/resemblance_variants.py`, 20 seeds per variant, same functions as the "
          "screen. Exploratory design validation, not evidence about the record."), "",
         (f"Observed: hs {obs['hs']:.4f}, rich {obs['rich']:.3f}, ht {obs['ht']:.4f}, two-cluster F_ST {obs['fst_clusters']:.4f}. "
         f"Operating point {OP}, stored length {S}, L_d {L_d}, burn-in {burnin}."), "",
         "| variant | mean hs | mean rich | mean ht | F_ST median [10%, 90%] |", "|---|---|---|---|---|",
         *[f"| {k} | {r.hs:.4f} | {r.rich:.3f} | {r.ht:.4f} | {r.fst_median:.4f} [{r.fst_q10:.4f}, {r.fst_q90:.4f}] |" for k, r in g.iterrows()]]
    FIND.write_text("\n".join(L) + "\n"); df.to_csv(OUTD / "resemblance_variants.csv", index=False)
    print("\n".join(L[4:])); print(f"\nwrote {FIND}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
