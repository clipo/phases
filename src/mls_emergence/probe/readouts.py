"""Label-invariant readouts of an A1 fit, in numpy (spec 4.3; D13, D14, D15, D18; R2 items 4, 5).

Each is defined once here and mirrored in the generated quantities of
stan/spatial/composition_gp_groups.stan; tests/probe/test_readouts.py pins
them on constructed cases with analytic answers, and the validation script
compares the Stan values with these on the same draws (rule 3).
"""
from __future__ import annotations

import numpy as np


def delta_contrast(delta: np.ndarray, weights: np.ndarray) -> float:
    """sqrt( sum_c sum_g w_g (delta_gc - mean_c)^2 ), mean_c = sum_g w_g delta_gc (D15).

    With weights = pi this is Delta; with the mean responsibility per group
    it is Delta_w. Weights must be non-negative and sum to one.
    """
    d = np.asarray(delta, float); w = np.asarray(weights, float)
    if d.ndim != 2 or w.shape != (d.shape[0],) or (w < 0).any() or not np.isclose(w.sum(), 1.0):
        raise ValueError("delta must be (G, C) and weights a probability vector over its G rows")
    m = w @ d
    return float(np.sqrt(np.sum(w[:, None] * (d - m[None, :]) ** 2)))


def comembership(resp: np.ndarray) -> np.ndarray:
    """sum_g resp[i,g] resp[j,g] for i != j, diagonal 1 (D13; the diagonal is a constant by definition)."""
    r = np.asarray(resp, float)
    q = r @ r.T
    np.fill_diagonal(q, 1.0)
    return q


def coherence(q: np.ndarray, edges: np.ndarray) -> float:
    """Mean co-membership over the graph edges minus the mean over all unordered pairs (D14, D18).

    Positive: neighbors share a group more often than arbitrary pairs do.
    Zero when one group holds every site. Works for a soft co-membership
    matrix (conditional-mean coherence) and for a hard one (0/1 entries).
    """
    q = np.asarray(q, float); e = np.asarray(edges, int)
    n = len(q)
    iu = np.triu_indices(n, 1)
    all_mean = float(q[iu].mean())
    if len(e) == 0:
        return 0.0 - all_mean
    return float(q[e[:, 0], e[:, 1]].mean()) - all_mean


def hard_comembership(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z)
    return (z[:, None] == z[None, :]).astype(float)


def block_matrix(q: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Mean co-membership within and between true groups: entry (a, b) averages q over pairs
    (i in a, j in b, i != j). Diagonal entries are within-group means over pairs of members."""
    q = np.asarray(q, float); t = np.asarray(truth)
    groups = sorted(set(t.tolist()))
    B = np.full((len(groups), len(groups)), np.nan)
    for a, ga in enumerate(groups):
        ia = np.where(t == ga)[0]
        for b, gb in enumerate(groups):
            ib = np.where(t == gb)[0]
            if a == b:
                if len(ia) < 2:
                    continue
                sub = q[np.ix_(ia, ia)]
                B[a, b] = sub[np.triu_indices(len(ia), 1)].mean()
            else:
                B[a, b] = q[np.ix_(ia, ib)].mean()
    return B


def separation_margins(q: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Per true group: within-group mean minus the largest between-group mean (within minus between, rule 21).

    Zero is the merge case; one is perfect recovery. Expected on the 12/11/5
    truth: all-in-one (0, 0, 0); 5-site group merged into the 11-site group
    (1, 0, 0); perfect (1, 1, 1) (R2 item 4).
    """
    B = block_matrix(q, truth)
    out = np.empty(len(B))
    for a in range(len(B)):
        between = np.delete(B[a], a)
        out[a] = B[a, a] - np.nanmax(between)
    return out
