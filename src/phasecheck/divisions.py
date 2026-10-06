"""Alternative divisions of the same assemblages into groups of the same sizes.

Two kinds, as in the paper. A division "around random centers" picks one
assemblage per group as a center and assigns every assemblage to a center at
least total squared distance, keeping the group sizes exact. A "compact"
division starts from one of those and swaps pairs between groups until no swap
lowers the within-group spread.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def spread(pts: np.ndarray, labels: np.ndarray) -> float:
    """Mean squared distance to the group center, per assemblage (km^2)."""
    return float(sum(((pts[labels == g] - pts[labels == g].mean(0)) ** 2).sum()
                     for g in np.unique(labels))) / len(pts)


def assign_exact(pts: np.ndarray, centers: np.ndarray, slot: np.ndarray) -> np.ndarray:
    """Least-cost assignment of assemblages to centers under exact group sizes."""
    cost = ((pts[:, None, :] - centers[None, :, :]) ** 2).sum(2)[:, slot]
    r, c = linear_sum_assignment(cost)
    out = np.empty(len(pts), int)
    out[r] = slot[c]
    return out


def make_compact(pts: np.ndarray, labels: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Lower the spread by swapping pairs between groups until no swap helps."""
    lab = labels.copy()
    k = len(sizes)
    while True:
        sums = np.array([pts[lab == g].sum(0) for g in range(k)])
        gain, move = 1e-12, None
        for a in range(k):
            ia = np.where(lab == a)[0]
            for b in range(a + 1, k):
                ib = np.where(lab == b)[0]
                pa, pb = pts[ia][:, None, :], pts[ib][None, :, :]
                sa, sb = sums[a] - pa + pb, sums[b] - pb + pa
                d = (((sa ** 2).sum(-1) - (sums[a] ** 2).sum()) / sizes[a]
                     + ((sb ** 2).sum(-1) - (sums[b] ** 2).sum()) / sizes[b])
                i, j = np.unravel_index(np.argmax(d), d.shape)
                if d[i, j] > gain:
                    gain, move = float(d[i, j]), (ia[i], ib[j])
        if move is None:
            return lab
        i, j = move
        lab[i], lab[j] = lab[j], lab[i]


def alternatives(pts: np.ndarray, group_index: np.ndarray, n_alt: int, rng: np.random.Generator):
    """`n_alt` random-center divisions and their compact versions, at the given group sizes."""
    if n_alt < 1:
        raise ValueError("n_alt must be at least 1")
    k = len(np.unique(group_index))
    sizes = np.bincount(group_index)
    slot = np.repeat(np.arange(k), sizes)
    around, compact = [], []
    for _ in range(n_alt):
        a = assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)
        around.append(a)
        compact.append(make_compact(pts, a, sizes))
    return around, compact
