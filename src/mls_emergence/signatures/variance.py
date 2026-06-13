from __future__ import annotations
import numpy as np

def gini_simpson(counts: np.ndarray) -> float:
    c = np.asarray(counts, float); N = c.sum()
    if N == 0: return 0.0
    p = c / N
    return float(1.0 - np.sum(p ** 2))

def cultural_fst(group_counts: np.ndarray) -> float:
    """F_ST = (H_T - H_S)/H_T via Gini-Simpson. rows=groups, cols=types.
    H_S = size-weighted mean within-group diversity; H_T = total-pool diversity."""
    g = np.asarray(group_counts, float)
    sizes = g.sum(axis=1)
    if sizes.sum() == 0: return 0.0
    H_within = np.array([gini_simpson(row) for row in g])
    H_S = float(np.average(H_within, weights=sizes))
    H_T = gini_simpson(g.sum(axis=0))
    if H_T == 0: return 0.0
    return float((H_T - H_S) / H_T)
