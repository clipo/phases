"""Resemblance summaries of A0's prior- and posterior-predictive counts (rule 20a).

Pilot step 4 of the latent-group probe (spec 3.3, 4.1a). Reads the replicated
counts written by refit_a0_basin.R and the observed counts through the basin
loader, and computes the four summaries the operating-point screen used, with
the screen's own functions: `diversity()` from analyses/47_revision_analysis.py
(mean within-assemblage Gini-Simpson, mean class richness, pooled Gini-Simpson)
and the plug-in `cultural_fst` on counts aggregated by the screen's spatial
clusters (k-means on centered coordinates, k chosen by silhouette among 2 to 4,
seed 7, as in 47.load_sets). The cluster labels are written with their hash.

Exploratory design validation; not evidence about the record.
Usage: .venv/bin/python analyses/probe_latent_groups/a0_predictive_summaries.py
"""
from __future__ import annotations

import hashlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.dataio.basin import data_fingerprint, load_curated
from mls_emergence.signatures.variance import cultural_fst

OUTD = ROOT / "output" / "probe_latent_groups"
FIND = ROOT / "output" / "findings" / "a0_predictive_summaries.md"
SUMMARIES = ("hs", "rich", "ht", "fst_clusters")
LABEL = {"hs": "mean within-assemblage Gini-Simpson", "rich": "mean class richness",
         "ht": "pooled Gini-Simpson", "fst_clusters": "plug-in F_ST between the spatial clusters"}


def spatial_clusters(coords: np.ndarray) -> np.ndarray:
    """The screen's rule (47.load_sets): k-means on centered coordinates, k by silhouette among 2..4, seed 7."""
    mf = importlib.import_module("make_figures")
    centered = coords - coords.mean(0)
    labs = {k: mf._kmeans_labels(centered, k, seed=7) for k in range(2, 5)}
    k = max(labs, key=lambda k: mf.silhouette_mean(centered, labs[k]))
    return np.asarray(labs[k])


def summaries(m: np.ndarray, labels: np.ndarray, diversity) -> dict:
    d = diversity(m)
    d["fst_clusters"] = float(cultural_fst(np.array([m[labels == g].sum(0) for g in np.unique(labels)])))
    return d


def read_yrep(path: Path, N: int, K: int) -> np.ndarray:
    df = pd.read_csv(path)
    n_draws = int(df["draw"].max())
    arr = np.zeros((n_draws, N, K), dtype=int)
    arr[df["draw"].to_numpy() - 1, df["i"].to_numpy() - 1, df["k"].to_numpy() - 1] = df["count"].to_numpy()
    return arr


def main() -> int:
    m47 = importlib.import_module("47_revision_analysis")
    counts, coords = load_curated()
    y = counts.to_numpy(int); xy = coords.to_numpy(float)
    fp = data_fingerprint(counts, coords)
    labels = spatial_clusters(xy)
    lab_json = {"assemblage": list(map(str, counts.index)), "cluster": [int(v) for v in labels],
                "rule": "k-means on centered lat/lon, k by silhouette among 2..4, seed 7 (47.load_sets)",
                "sha256": hashlib.sha256(",".join(map(str, labels)).encode()).hexdigest()}
    OUTD.mkdir(parents=True, exist_ok=True)
    (OUTD / "spatial_clusters.json").write_text(json.dumps(lab_json, indent=1))
    obs = summaries(y, labels, m47.diversity)

    rows = []
    for kind in ("prior", "posterior"):
        path = OUTD / f"a0_{kind}_predictive_yrep.csv"
        if not path.is_file():
            print(f"missing {path}; run refit_a0_basin.R first"); return 1
        arr = read_yrep(path, *y.shape)
        vals = {s: [] for s in SUMMARIES}
        for d in range(arr.shape[0]):
            sm = summaries(arr[d], labels, m47.diversity)
            for s in SUMMARIES:
                vals[s].append(sm[s])
        for s in SUMMARIES:
            v = np.asarray(vals[s], float); v = v[np.isfinite(v)]
            rows.append({"kind": kind, "summary": s, "n_draws": int(len(v)), "median": float(np.median(v)),
                         "q5": float(np.quantile(v, .05)), "q95": float(np.quantile(v, .95)),
                         "observed": obs[s], "frac_below_obs": float(np.mean(v < obs[s]))})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUTD / "a0_predictive_summaries.csv", index=False)

    intro = ("Produced by `analyses/probe_latent_groups/a0_predictive_summaries.py` from the counts written by "
             "`refit_a0_basin.R` (200 draws each). Summaries computed with the operating-point screen's own functions "
             "(`47_revision_analysis.diversity`, plug-in `cultural_fst` on the screen's spatial clusters). "
             "Exploratory design validation, not evidence about the record.")
    data_line = (f"Observed counts: {y.shape[0]} assemblages, {int(y.sum())} sherds, counts sha256 `{fp['sha256_counts']}`. "
                 f"Spatial clusters: {len(np.unique(labels))} (sizes {', '.join(str(int(v)) for v in np.bincount(labels))}), "
                 f"labels sha256 `{lab_json['sha256']}`, in `output/probe_latent_groups/spatial_clusters.json`.")
    L = ["# A0 prior- and posterior-predictive resemblance summaries (rule 20a)", "", intro, "", data_line, "",
         "| summary | observed | prior predictive median [5%, 95%] | fraction of prior draws below observed | posterior predictive median [5%, 95%] | fraction of posterior draws below observed |",
         "|---|---|---|---|---|---|"]
    for s in SUMMARIES:
        p = tab[(tab.kind == "prior") & (tab.summary == s)].iloc[0]; q = tab[(tab.kind == "posterior") & (tab.summary == s)].iloc[0]
        cell = lambda r: f"{r['median']:.4f} [{r['q5']:.4f}, {r['q95']:.4f}] | {r['frac_below_obs']:.2f}"
        L.append(f"| {LABEL[s]} | {obs[s]:.4f} | {cell(p)} | {cell(q)} |")
    reading = ("Reading: the prior-predictive column says what the priors alone generate on this geometry and these sherd totals; "
               "a fraction near 0 or 1 means the observed value sits in a tail of the prior predictive. The posterior-predictive "
               "column is the ordinary predictive check of the fitted model on the same four summaries.")
    L += ["", reading]
    FIND.parent.mkdir(parents=True, exist_ok=True); FIND.write_text("\n".join(L) + "\n")
    print("\n".join(L[6:]))
    print(f"\nwrote {FIND}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
