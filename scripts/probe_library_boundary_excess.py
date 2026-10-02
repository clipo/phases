#!/usr/bin/env python3
"""What would a minimum-pairs rule change in the library boundary excess?

`signatures.assortativity.boundary_excess` averages a within-minus-between
similarity gap over four distance bins and counts a bin as soon as it holds one
pair of each kind. The analysis-level statistic at the phase lines
(`23_phases_vs_spatial_drift.boundary_excess_labeled`) had the same rule until
2026-10-02, when one pair in one bin was found to decide its sign, and it now
needs five pairs of each kind. This script measures what the same rule would
do to the LIBRARY function, which the per-period trajectories call on four to
ten assemblages at a time.

It runs one analysis script with the library function wrapped: the wrapper
returns the ORIGINAL value to the caller (so nothing the script writes
changes) and records, for every call, the pairs in each counted bin and the
value a five-pair rule would have returned.

MEASURED 2026-10-02 (basin, after the Phillips-label rerun):

  17_basin_results.py          6 calls, 5 with a bin under five pairs, 1 sign change
  12_sensitivity_grid.py     120 calls, 75 with a bin under five pairs, 30 sign changes
  05_empirical_application.py  6 calls, 3, 1
  06_empirical_two_level.py    6 calls, 3, 1

WHAT FOLLOWS. On windows of four to ten assemblages almost every counted bin
holds one or two pairs, so a five-pair rule empties the bins and the function
falls back to its raw within-minus-between gap, which is a different statistic,
not a repaired one. The per-period spatial measure is therefore unreliable at
this resolution whichever rule is used. That is what the supplement already
says of it ("cannot discriminate at this resolution ... shown for
completeness", Figure S3), and no claim in either manuscript rests on it. The
library function is left as it is, with this limit recorded in its docstring;
the phase-line statistic, which the paper does report, carries the five-pair
rule.

Usage: PYTHONPATH=src .venv/bin/python scripts/probe_library_boundary_excess.py 17_basin_results.py
"""
import importlib, sys, runpy, io, contextlib
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "analyses", ROOT / "src"):
    sys.path.insert(0, str(p))
import mls_emergence.signatures.assortativity as A
orig = A.boundary_excess
log = []
def variant(counts, coords, n_clusters=None, n_bins=4, seed=0, min_pairs=5):
    counts = np.asarray(counts, float); coords = np.asarray(coords, float)
    n = counts.shape[0]
    if n < 4: return 0.0, None
    k = n_clusters if n_clusters is not None else A._select_k(coords, seed=seed)
    labels = A._kmeans_labels(coords, k, seed=seed)
    sim = A.similarity_matrix(counts); dist = A.geo_distance(coords)
    iu = np.triu_indices(n, k=1); d = dist[iu]; s = sim[iu]
    same = labels[iu[0]] == labels[iu[1]]
    if same.sum() == 0 or (~same).sum() == 0: return 0.0, None
    if d.max() == d.min(): return float(s[same].mean() - s[~same].mean()), None
    edges = np.linspace(d.min(), d.max() + 1e-9, n_bins + 1)
    gaps, pairs = [], []
    for b in range(n_bins):
        inb = (d >= edges[b]) & (d < edges[b + 1])
        w, bt = inb & same, inb & ~same
        if w.sum() == 0 or bt.sum() == 0: continue
        pairs.append((int(w.sum()), int(bt.sum())))
        if w.sum() >= min_pairs and bt.sum() >= min_pairs:
            gaps.append(float(s[w].mean() - s[bt].mean()))
    raw = float(s[same].mean() - s[~same].mean())
    return (float(np.mean(gaps)) if gaps else raw), pairs
def wrapped(counts, coords, *a, **k):
    v0 = orig(counts, coords, *a, **k)
    v5, pairs = variant(counts, coords, *a, **k)
    log.append((np.asarray(counts).shape[0], v0, v5, pairs))
    return v0
A.boundary_excess = wrapped
import mls_emergence.signatures as SG
for mod in list(sys.modules.values()):
    if mod and getattr(mod, "boundary_excess", None) is orig:
        mod.boundary_excess = wrapped
script = sys.argv[1]
mf = importlib.import_module("make_figures")
if getattr(mf, "boundary_excess", None) is orig: mf.boundary_excess = wrapped
with contextlib.redirect_stdout(io.StringIO()):
    try:
        runpy.run_path(str(ROOT / "analyses" / script), run_name="__main__")
    except SystemExit:
        pass
small = [r for r in log if r[3] and any(min(p) < 5 for p in r[3])]
print(script, "calls", len(log), "| calls with a bin under 5 pairs:", len(small))
diffs = [(n, v0, v5, pr) for n, v0, v5, pr in log if abs(v0 - v5) > 1e-9]
print("  calls whose value would change:", len(diffs), "| sign flips:", sum(1 for n, v0, v5, pr in diffs if v0 * v5 < 0))
for n, v0, v5, pr in diffs[:12]:
    print(f"   n={n:3d}  now {v0:+7.2f}  five-pair {v5:+7.2f}  bins(within,between)={pr}")
