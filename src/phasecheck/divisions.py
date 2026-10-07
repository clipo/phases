"""Alternative divisions of the same assemblages into groups of the same sizes.

Two kinds, as in the paper. A division "around random centers" picks one
assemblage per group as a center and assigns every assemblage to a center at
least total squared distance, keeping the group sizes exact. A "compact"
division starts from one of those and swaps pairs between groups until no swap
lowers the within-group spread.

These are the "other lines on the same map" against which the phase lines are
judged in questions 2 and 3 (`questions.compare_division`) and 9
(`copying.power`). Both kinds keep each group's size exactly equal to a
phase's size, so a comparison between them and the phases differs only in
where the lines fall. Points are flat east/north coordinates in km from
`geo.planar_km`; spread is in km^2. The random source is a numpy Generator
passed in by the caller, which owns the seed.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def spread(pts: np.ndarray, labels: np.ndarray) -> float:
    """Mean squared distance to the group center, per assemblage (km^2).

    Parameters
    ----------
    pts : numpy.ndarray, shape (n, 2)
        Flat coordinates in km.
    labels : array_like, shape (n,)
        Group of each assemblage (any hashable labels).

    Returns
    -------
    float
        The within-group sum of squares divided by n. Lower is tighter.
    """
    return float(sum(((pts[labels == g] - pts[labels == g].mean(0)) ** 2).sum()
                     for g in np.unique(labels))) / len(pts)


def assign_exact(pts: np.ndarray, centers: np.ndarray, slot: np.ndarray) -> np.ndarray:
    """Least-cost assignment of assemblages to centers under exact group sizes.

    Parameters
    ----------
    pts : numpy.ndarray, shape (n, 2)
        Flat coordinates in km.
    centers : numpy.ndarray, shape (k, 2)
        One center per group.
    slot : numpy.ndarray, shape (n,)
        The group of each of n slots, group g repeated as many times as its
        size (so `slot` fixes the sizes).

    Returns
    -------
    numpy.ndarray of int, shape (n,)
        The group of each assemblage, minimizing total squared distance to
        its center subject to the sizes in `slot`.
    """
    # An n x n assignment problem: one column per slot, a slot costing the
    # squared distance to its group's center. Solving it fills every slot once.
    cost = ((pts[:, None, :] - centers[None, :, :]) ** 2).sum(2)[:, slot]
    r, c = linear_sum_assignment(cost)
    out = np.empty(len(pts), int)
    out[r] = slot[c]
    return out


def make_compact(pts: np.ndarray, labels: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Lower the spread by swapping pairs between groups until no swap helps.

    Each pass makes the single best swap of one assemblage in group a for one
    in group b, over every pair of groups, and stops when no swap lowers the
    spread. Swaps keep the group sizes; the result is a local optimum, so
    different starts give different compact divisions.

    Parameters
    ----------
    pts : numpy.ndarray, shape (n, 2)
        Flat coordinates in km.
    labels : numpy.ndarray of int, shape (n,)
        Starting division, groups numbered 0..k-1. Not modified.
    sizes : numpy.ndarray of int, shape (k,)
        Size of each group (`np.bincount(labels)`).

    Returns
    -------
    numpy.ndarray of int, shape (n,)
        The compact division, with the same group sizes.
    """
    lab = labels.copy()
    k = len(sizes)
    while True:
        sums = np.array([pts[lab == g].sum(0) for g in range(k)])
        # With sizes fixed, within-group sum of squares = sum|x|^2 - sum_g |S_g|^2 / n_g,
        # so a swap lowers the spread by the rise in sum_g |S_g|^2 / n_g, which is `d`.
        # A gain must exceed 1e-12 to count, so rounding cannot loop forever.
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
    """`n_alt` random-center divisions and their compact versions, at the given group sizes.

    Parameters
    ----------
    pts : numpy.ndarray, shape (n, 2)
        Flat coordinates in km.
    group_index : numpy.ndarray of int, shape (n,)
        The tested division, groups numbered 0..k-1 with none empty; only
        its group sizes are used.
    n_alt : int
        Number of divisions of each kind, at least 1. Duplicates are not
        removed (callers report how many are distinct).
    rng : numpy.random.Generator
        Source of the random centers.

    Returns
    -------
    (list, list)
        `around`: `n_alt` label arrays, each from k distinct assemblages
        drawn as centers and an exact-size least-cost assignment to them.
        `compact`: the same divisions after `make_compact`, in the same order.

    Raises
    ------
    ValueError
        If `n_alt` is below 1.
    """
    if n_alt < 1:
        raise ValueError("n_alt must be at least 1")
    k = len(np.unique(group_index))
    sizes = np.bincount(group_index)
    slot = np.repeat(np.arange(k), sizes)          # group g repeated sizes[g] times
    around, compact = [], []
    for _ in range(n_alt):
        a = assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)
        around.append(a)
        compact.append(make_compact(pts, a, sizes))
    return around, compact
