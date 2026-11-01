"""Barrier placements from geography only (A3 design, section 3; D2).

G1: the two spatial clusters the screen uses (k-means on centered
coordinates, k = 2, seed 7; a partition from site locations alone). G2: a
declared candidate set of straight lines across the basin, each splitting
the sites by side; the placement index is a discrete latent variable in
stan/spatial/composition_gp_barrier.stan. Phillips's areas (G3) are not
produced here; they enter only as the descriptive comparison of D5 part 2.
"""
from __future__ import annotations

import numpy as np

from mls_emergence.dataio.basin import km_xy

MIN_SIDE = 3


def kmeans2_regions(xy_latlon: np.ndarray, seed: int = 7, n_init: int = 50) -> np.ndarray:
    """Two spatial clusters by k-means on centered coordinates (the screen's rule at k = 2), labels 1 and 2.

    Deterministic restarts from the given seed; the lower-mean-easting cluster is label 1.
    """
    p = km_xy(np.asarray(xy_latlon, float))
    rng = np.random.default_rng(seed)
    best, best_inertia = None, np.inf
    for _ in range(n_init):
        c = p[rng.choice(len(p), 2, replace=False)]
        for _ in range(100):
            lab = np.argmin(((p[:, None, :] - c[None, :, :]) ** 2).sum(-1), axis=1)
            new = np.array([p[lab == k].mean(0) if (lab == k).any() else c[k] for k in range(2)])
            if np.allclose(new, c):
                break
            c = new
        inertia = ((p - c[lab]) ** 2).sum()
        if inertia < best_inertia - 1e-9:
            best, best_inertia = lab.copy(), inertia
    order = np.argsort([p[best == k, 0].mean() for k in range(2)])
    out = np.empty(len(p), int)
    for new_label, k in enumerate(order, start=1):
        out[best == k] = new_label
    return out


def line_regions(xy_latlon: np.ndarray, angle_deg: float, offset_km: float) -> np.ndarray:
    """Side of a line: region 1 where x cos(a) + y sin(a) < offset (centered km coordinates), else 2."""
    p = km_xy(np.asarray(xy_latlon, float))
    a = np.deg2rad(angle_deg)
    proj = p[:, 0] * np.cos(a) + p[:, 1] * np.sin(a)
    return np.where(proj < offset_km, 1, 2)


def candidate_lines(xy_latlon: np.ndarray, angles=(0, 30, 60, 90, 120, 150), n_offsets: int = 5,
                    min_side: int = MIN_SIDE) -> list[dict]:
    """The declared candidate set: for each angle, offsets at evenly spaced quantiles of the projections.

    A candidate is kept when both sides hold at least `min_side` sites and
    its partition is not identical to an earlier candidate's. Returns dicts
    with angle, offset, regions (labels 1 and 2) and sides (counts).
    """
    p = km_xy(np.asarray(xy_latlon, float))
    seen = set(); out = []
    qs = np.linspace(0, 1, n_offsets + 2)[1:-1]
    for ang in angles:
        a = np.deg2rad(ang); proj = p[:, 0] * np.cos(a) + p[:, 1] * np.sin(a)
        for q in qs:
            off = float(np.quantile(proj, q))
            reg = np.where(proj < off, 1, 2)
            if min(np.sum(reg == 1), np.sum(reg == 2)) < min_side:
                continue
            key = tuple(reg.tolist())
            if key in seen:
                continue
            seen.add(key)
            out.append({"angle": float(ang), "offset": off, "regions": reg, "sides": (int(np.sum(reg == 1)), int(np.sum(reg == 2)))})
    return out


def best_line_for(partition: np.ndarray, candidates: list[dict]) -> dict:
    """The candidate line whose split agrees best with a given two-way partition (Rand-style agreement).

    Used only to pick a TEST placement for simulated records; never for inference.
    """
    t = np.asarray(partition)
    best, score = None, -1.0
    for c in candidates:
        r = c["regions"]
        agree = np.mean((t[:, None] == t[None, :]) == (r[:, None] == r[None, :]))
        if agree > score:
            best, score = c, agree
    return dict(best, agreement=float(score))
