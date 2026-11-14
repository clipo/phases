"""Rival-cause generators for the per-type step (D30; spec 6.2 and 6.4).

Each is a switch on the D22 neutral simulator's output, applied by
analyses/per_type/per_class_table.py: C2 ties the record's totals to span
length; R1 relabels sherds at a low rate; R2 draws sherds in clumps of one
class each (vessels counted as several sherds). None has an anchor in the
repository for its rate; every rate is swept (D28).
"""
from __future__ import annotations

import numpy as np


def totals_by_span_length(totals: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """C2: the record's totals permuted so the longest span gets the largest total (ties by index)."""
    totals = np.asarray(totals); lengths = np.asarray(lengths)
    if totals.shape != lengths.shape or totals.ndim != 1:
        raise ValueError("totals and lengths must be one-dimensional and the same length")
    order = np.lexsort((np.arange(len(lengths)), lengths))     # ascending length, ties by index
    out = np.empty_like(totals); out[order] = np.sort(totals)
    return out


def misclassify(y: np.ndarray, rate: float, pooled: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """R1: each sherd relabeled with probability `rate` to a class drawn from `pooled`; row sums unchanged."""
    if not 0.0 <= rate <= 1.0:
        raise ValueError("rate must be a probability")
    y = np.asarray(y, int); pooled = np.asarray(pooled, float); pooled = pooled / pooled.sum()
    if y.shape[1] != len(pooled):
        raise ValueError("pooled must have one entry per class")
    out = np.zeros_like(y)
    for i in range(y.shape[0]):
        moved = rng.binomial(y[i], rate)
        out[i] = y[i] - moved + rng.multinomial(int(moved.sum()), pooled)
    return out


def clumped_counts(p: np.ndarray, totals: np.ndarray, mean_clump: float, rng: np.random.Generator) -> np.ndarray:
    """R2: per site, clumps of one class each (class from p[i], size geometric with mean `mean_clump`), truncated to the total."""
    if mean_clump < 1.0:
        raise ValueError("mean_clump must be at least 1")
    p = np.asarray(p, float); totals = np.asarray(totals, int)
    if p.shape[0] != len(totals):
        raise ValueError("one profile per site is required")
    out = np.zeros((len(totals), p.shape[1]), int)
    for i, n in enumerate(totals):
        q = p[i] / p[i].sum(); left = int(n)
        while left > 0:
            size = min(int(rng.geometric(1.0 / mean_clump)), left)
            out[i, rng.choice(len(q), p=q)] += size; left -= size
    return out


def plant_marker(prof: np.ndarray, labels: np.ndarray, cluster, c: int, m: float) -> np.ndarray:
    """M (regional-drift spec 4.1): class c's unnormalized weight times m in every row whose label equals `cluster`, renormalized.

    m is an odds multiplier: the planted share is m p / (1 + (m - 1) p). Other
    rows are untouched; the input is not modified. m must be finite and
    positive (m = 1 returns a copy of the input); the cluster label must occur;
    c must index a class; rows must sum to 1 within 1e-9.
    """
    prof = np.asarray(prof, float); labels = np.asarray(labels)
    if prof.ndim != 2 or labels.shape != (prof.shape[0],):
        raise ValueError("prof must be N x k with one label per row")
    if not np.isfinite(m) or m <= 0:
        raise ValueError("m must be finite and positive")
    if not np.any(labels == cluster):
        raise ValueError(f"cluster label {cluster!r} does not occur in labels")
    if isinstance(c, (bool, np.bool_)) or not (isinstance(c, (int, np.integer)) and 0 <= c < prof.shape[1]):
        raise ValueError("c must index a class")
    if (prof < 0).any() or not np.allclose(prof.sum(1), 1.0, atol=1e-9, rtol=0):
        raise ValueError("profile rows must be non-negative and sum to 1")
    out = prof.copy()
    if m == 1:
        return out                      # the identity, bit for bit: renormalizing rows that sum to 1 within an ulp would move them
    rows = labels == cluster
    out[rows, c] *= float(m)
    out[rows] /= out[rows].sum(1, keepdims=True)
    return out
