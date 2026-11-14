"""Two plug-in spatial readouts of a counts matrix (drift-extension spec 2026-11-14, draft 3, section 4).

Both are computed on observed shares p = y / n per site, by one function each,
applied identically to the record and to simulated records. They are not the
posterior quantities the main text reports for the record.

Sign conventions (rule 21a; direction tests in tests/signatures/test_sign_conventions.py):
  lat_rank             positive means the class's share is higher to the north.
  distance_decay_rank  positive means sites farther apart differ more in composition.

A class absent from a site has share zero there and the site is kept. A
correlation with a constant input is undefined: the value is NaN and the
reason string is "constant". Out-of-domain input raises (rule 5).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata

CONSTANT = "constant"


def _shares(y) -> np.ndarray:
    y = np.asarray(y, float)
    if y.ndim != 2 or y.shape[0] < 2 or y.shape[1] < 1:
        raise ValueError("counts must be a sites-by-classes matrix with at least two sites")
    if not np.isfinite(y).all() or (y < 0).any():
        raise ValueError("counts must be finite and non-negative")
    n = y.sum(1)
    if (n <= 0).any():
        raise ValueError("a site with zero total has no shares")
    return y / n[:, None]


def _rank_corr(a: np.ndarray, b: np.ndarray) -> tuple[float, str]:
    """Spearman correlation with average ranks for ties: (value, "") or (NaN, "constant")."""
    ra, rb = rankdata(a, method="average"), rankdata(b, method="average")
    da, db = ra - ra.mean(), rb - rb.mean()
    saa, sbb = float((da * da).sum()), float((db * db).sum())
    if saa == 0.0 or sbb == 0.0:
        return float("nan"), CONSTANT
    return float((da * db).sum() / np.sqrt(saa * sbb)), ""


def lat_rank(y, latitude) -> tuple[np.ndarray, list[str]]:
    """Per class, the Spearman correlation of the class's share with latitude over the sites.

    Positive means the share is higher to the north. Returns the values and one
    reason string per class ("" or "constant").
    """
    p = _shares(y); lat = np.asarray(latitude, float)
    if lat.shape != (p.shape[0],):
        raise ValueError("latitude must have one entry per site")
    if not np.isfinite(lat).all() or np.all(lat == lat[0]):
        raise ValueError("latitudes must be finite and not all equal")
    res = [_rank_corr(p[:, k], lat) for k in range(p.shape[1])]
    return np.array([v for v, _ in res]), [why for _, why in res]


def a_lat(y, latitude) -> np.ndarray:
    """|lat_rank| per class: the strength of a latitude trend whatever its direction. NaN where lat_rank is."""
    return np.abs(lat_rank(y, latitude)[0])


def tv_pairs(y) -> np.ndarray:
    """Total-variation distance 0.5 * sum_k |p_ik - p_jk| for the unique pairs i < j, in row-major order."""
    p = _shares(y); i, j = np.triu_indices(p.shape[0], k=1)
    return 0.5 * np.abs(p[i] - p[j]).sum(1)


def distance_decay_rank(y, D) -> tuple[float, str]:
    """Spearman correlation of pairwise total-variation distance with geographic distance over the unique pairs.

    Positive means sites farther apart differ more in composition. Returns the
    value and a reason string ("" or "constant").
    """
    tv = tv_pairs(y); D = np.asarray(D, float); n = np.asarray(y).shape[0]
    if D.shape != (n, n):
        raise ValueError("the distance matrix must be sites by sites")
    if not np.isfinite(D).all() or (D < 0).any() or not np.array_equal(D, D.T):
        raise ValueError("the distance matrix must be finite, non-negative and symmetric")
    i, j = np.triu_indices(n, k=1)
    return _rank_corr(tv, D[i, j])
