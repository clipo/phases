"""Sampling of finite archaeological count assemblages."""
import numpy as np


def rarefy(counts, depth, rng):
    """Draw without replacement; retain rows smaller than ``depth`` in full."""
    a = np.asarray(counts)
    if a.ndim != 2 or not np.isfinite(a).all() or (a < 0).any():
        raise ValueError("counts must be a finite, nonnegative matrix")
    if not np.equal(a, np.floor(a)).all():
        raise ValueError("rarefaction requires integer counts")
    if int(depth) != depth or depth < 0:
        raise ValueError("depth must be a nonnegative integer")
    a = a.astype(np.int64)
    out = np.empty_like(a)
    for i, row in enumerate(a):
        n = min(int(depth), int(row.sum()))
        out[i] = row if n == row.sum() else rng.multivariate_hypergeometric(row, n)
    return out
