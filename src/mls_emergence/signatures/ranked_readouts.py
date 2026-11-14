"""Readouts on classes ranked by pooled count (type-origination spec 2026-11-14, draft 1, section 5).

Under origination a simulated class has no counterpart among the record's
named classes, so every statistic is defined by rank: classes are ordered by
pooled count, largest first, ties broken by column order, and classes with a
pooled count of zero are dropped. The same function is applied to the record
and to simulated records. Signed readouts reuse `spatial_readouts`
(positive lat_rank: share higher to the north; positive D_dist: farther
sites differ more).
"""
from __future__ import annotations

import numpy as np

from mls_emergence.signatures import spatial_readouts as sr
from mls_emergence.signatures.variance import moment_f

N_COMMON = 3
N_RANKED = 10


def rank_order(y) -> np.ndarray:
    """Column indices of the classes present, by pooled count descending, ties in column order."""
    y = np.asarray(y)
    if y.ndim != 2 or (y < 0).any():
        raise ValueError("counts must be a non-negative sites-by-classes matrix")
    tot = y.sum(0); order = np.argsort(-tot, kind="stable")
    return order[tot[order] > 0]


def ranked_readouts(y, latitude, D) -> dict:
    """F and |latitude rank correlation| of ranks 1 to 3, occupancy of ranks 4 to 10, richness, share beyond rank 10, D_dist.

    A record with fewer than three classes has NaN at the missing F and
    latitude ranks; one with fewer than ten has occupancy 0 at the missing
    ranks (the class is present at no site).
    """
    y = np.asarray(y); order = rank_order(y); yr = y[:, order].astype(float); K = yr.shape[1]
    f = moment_f(yr) if K else np.array([]); lat, _why = sr.lat_rank(yr, latitude) if K else (np.array([]), [])
    out = {"K_obs": int(K), "beyond10": float(yr[:, N_RANKED:].sum() / yr.sum())}
    for r in range(N_COMMON):
        out[f"F_r{r + 1}"] = float(f[r]) if r < K else float("nan")
        out[f"lat_r{r + 1}"] = float(lat[r]) if r < K else float("nan")
    for r in range(N_COMMON, N_RANKED):
        out[f"occ_r{r + 1}"] = int((yr[:, r] > 0).sum()) if r < K else 0
    out["D_dist"] = sr.distance_decay_rank(yr, D)[0]
    return out
