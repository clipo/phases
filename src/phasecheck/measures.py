"""The measures: cultural F_ST, the boundary excess, agreement, and profile difference.

Every question reads its numbers from here:

  cultural F_ST        questions 2, 3, 7 and 9: how much groups differ in class
                       mix, (H_T - H_S) / H_T on Gini-Simpson diversity. 0 when
                       the groups have the same mix, 1 when they share no class.
  boundary excess      questions 2, 3 and 7: the step in Brainerd-Robinson
                       similarity at a set of lines, at matched distance.
  adjusted Rand index  questions 1 and 8: agreement between two divisions.
  profile difference   questions 5 and 7: the share of sherds (percent) that
                       would have to change class to make two profiles alike.

Sign conventions, defined once here and used oriented everywhere:

  boundary excess      within-group similarity MINUS between-group similarity,
                       at matched distance. Positive means assemblages on the
                       same side of a line are more alike than distance alone
                       predicts, so the lines separate assemblages; negative
                       means they are more alike ACROSS the lines. Units are
                       Brainerd-Robinson points (the scale runs 0 to 200).
  profile difference   always non-negative: 0 identical, 100 nothing in common.

The direction of each is pinned by `tests/phasecheck/test_phasecheck.py`
(`test_boundary_excess_direction`).
"""
from __future__ import annotations

import numpy as np
from scipy.special import comb

# A distance bin enters the boundary excess only when it holds at least this
# many pairs of each kind (same group, different groups); fewer would let one
# or two pairs set the bin's gap.
BOUNDARY_EXCESS_MIN_PAIRS = 5


def gini_simpson(counts: np.ndarray) -> float:
    """Gini-Simpson diversity, 1 - sum(p^2), of one row of counts (plug-in).

    Parameters
    ----------
    counts : array_like, shape (k,)
        Non-negative counts (or posterior expected counts) by class.

    Returns
    -------
    float
        Between 0 (one class) and 1 - 1/k (all classes equal): the chance
        that two sherds drawn with replacement differ in class.

    Raises
    ------
    ValueError
        If the row sums to zero or less.
    """
    c = np.asarray(counts, float)
    t = c.sum()
    if t <= 0:
        raise ValueError("diversity is undefined for an empty row of counts")
    p = c / t
    return float(1.0 - (p ** 2).sum())


def cultural_fst(group_counts: np.ndarray) -> float:
    """F_ST = (H_T - H_S) / H_T on Gini-Simpson diversity. Rows are groups.

    H_S is the size-weighted mean diversity within groups and H_T the diversity
    of all groups pooled. Zero means the groups have the same class mix; larger
    means they differ more. Refuses fewer than two non-empty groups and a pool
    holding a single class, where the quantity is undefined.

    Parameters
    ----------
    group_counts : array_like, shape (groups, classes)
        Pooled counts of each group. Empty groups are dropped first.

    Returns
    -------
    float
        F_ST, 0 to 1. Size weighting means a group with more sherds counts
        for more in H_S.

    Raises
    ------
    ValueError
        If the array is not 2-D, holds a negative count, has fewer than two
        non-empty groups, or its pooled counts hold a single class.
    """
    g = np.asarray(group_counts, float)
    if g.ndim != 2:
        raise ValueError(f"expected a groups-by-classes array, got shape {g.shape}")
    if np.any(g < 0):
        raise ValueError("counts must be non-negative")
    sizes = g.sum(1)
    g, sizes = g[sizes > 0], sizes[sizes > 0]
    if g.shape[0] < 2:
        raise ValueError(f"cultural F_ST needs at least two non-empty groups, got {g.shape[0]}")
    h_t = gini_simpson(g.sum(0))
    if h_t == 0:
        raise ValueError("cultural F_ST is undefined when the pooled counts hold a single class")
    h_s = float(np.average([gini_simpson(r) for r in g], weights=sizes))
    return float((h_t - h_s) / h_t)


def fst_by(counts: np.ndarray, labels: np.ndarray) -> float:
    """Cultural F_ST between the groups that `labels` assigns the assemblages to.

    Parameters
    ----------
    counts : array_like, shape (n, classes)
        Counts of each assemblage.
    labels : array_like, shape (n,)
        Group of each assemblage. Each group's assemblages are pooled
        (summed) before `cultural_fst`.

    Returns
    -------
    float
        As `cultural_fst`, which raises on the same conditions.
    """
    counts, labels = np.asarray(counts, float), np.asarray(labels)
    return cultural_fst(np.array([counts[labels == g].sum(0) for g in np.unique(labels)]))


def similarity_matrix(counts: np.ndarray) -> np.ndarray:
    """Brainerd-Robinson similarity between every pair of assemblages, 0 to 200.

    200 minus the summed absolute difference of the two percentage profiles:
    200 for identical profiles, 0 for profiles with no class in common.

    Parameters
    ----------
    counts : array_like, shape (n, classes)

    Returns
    -------
    numpy.ndarray, shape (n, n)
        Symmetric, 200 on the diagonal.

    Raises
    ------
    ValueError
        If any assemblage has no sherds.
    """
    c = np.asarray(counts, float)
    t = c.sum(1, keepdims=True)
    if np.any(t <= 0):
        raise ValueError("every assemblage needs at least one sherd")
    p = 100.0 * c / t
    return 200.0 - np.abs(p[:, None, :] - p[None, :, :]).sum(-1)


def boundary_excess(counts: np.ndarray, dist: np.ndarray, labels: np.ndarray, n_bins: int = 4,
                    min_pairs: int = BOUNDARY_EXCESS_MIN_PAIRS) -> float:
    """Within-group minus between-group similarity, at matched distance.

    The value alone; see `boundary_excess_detail` for the definition, the
    arguments and the case with no distance control.
    """
    return boundary_excess_detail(counts, dist, labels, n_bins, min_pairs)[0]


def boundary_excess_detail(counts: np.ndarray, dist: np.ndarray, labels: np.ndarray, n_bins: int = 4,
                           min_pairs: int = BOUNDARY_EXCESS_MIN_PAIRS) -> tuple:
    """Within-group minus between-group similarity, and the number of distance bins it rests on.

    Returns (value, bins_used). `bins_used` of 0 means no bin held enough pairs
    of each kind, and the value is the gap over all pairs with NO control for
    distance; a caller must say so.

    Pairs of assemblages are sorted into `n_bins` equal-width distance bins. In
    each bin holding at least `min_pairs` pairs of each kind, the mean
    Brainerd-Robinson similarity of same-group pairs minus that of
    different-group pairs is taken, and the bins are averaged. Where similarity
    only falls with distance the result is near zero; a positive value is a
    step in similarity at the lines beyond what distance predicts. If no bin
    qualifies, the gap over all pairs is returned without the distance control.

    Parameters
    ----------
    counts : array_like, shape (n, classes)
        Counts of each assemblage (observed, or one posterior draw).
    dist : array_like, shape (n, n)
        Pairwise distances in km (straight-line or supplied).
    labels : array_like, shape (n,)
        Group of each assemblage.
    n_bins : int
        Number of equal-width distance bins between the shortest and the
        longest pair distance. The report assumes the default, 4.
    min_pairs : int
        Pairs of each kind a bin needs to enter the average, at least 1.

    Returns
    -------
    (float, int)
        The excess in Brainerd-Robinson points (positive: the lines separate
        assemblages beyond what distance predicts), and the number of bins
        averaged. With no pair of one kind (a single group, or every
        assemblage alone in its group) the result is (0.0, 0); with every
        pair at the same distance it is the raw gap and 0 bins.

    Raises
    ------
    ValueError
        If `min_pairs` is below 1, or an assemblage has no sherds.
    """
    if min_pairs < 1:
        raise ValueError("min_pairs must be at least 1")
    labels = np.asarray(labels)
    s_full = similarity_matrix(counts)
    iu = np.triu_indices(s_full.shape[0], k=1)     # each unordered pair once
    d, s = np.asarray(dist, float)[iu], s_full[iu]
    same = labels[iu[0]] == labels[iu[1]]
    if same.sum() == 0 or (~same).sum() == 0:
        return 0.0, 0
    raw = float(s[same].mean() - s[~same].mean())
    if d.max() == d.min():
        return raw, 0
    # Bins are half-open [lo, hi); the top edge is nudged up so the longest pair
    # falls in the last bin rather than outside every bin.
    edges = np.linspace(d.min(), d.max() + 1e-9, n_bins + 1)
    gaps = []
    for b in range(n_bins):
        in_bin = (d >= edges[b]) & (d < edges[b + 1])
        w, btw = in_bin & same, in_bin & ~same
        if w.sum() < min_pairs or btw.sum() < min_pairs:
            continue
        gaps.append(float(s[w].mean() - s[btw].mean()))
    return (float(np.mean(gaps)), len(gaps)) if gaps else (raw, 0)


def adjusted_rand(a, b) -> float:
    """Adjusted Rand index: 1 for identical divisions, about 0 for unrelated ones.

    Group labels need not match between the two divisions; only who is
    grouped with whom counts. Negative values (agreement below chance) are
    possible.

    Parameters
    ----------
    a, b : array_like, shape (n,)
        Two divisions of the same assemblages, in the same order.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If the divisions differ in length, or the index is undefined (both
        a single group, or both every assemblage alone).
    """
    a, b = np.asarray(a), np.asarray(b)
    if len(a) != len(b):
        raise ValueError("the two divisions must label the same assemblages")
    ia = {v: i for i, v in enumerate(sorted(set(a.tolist())))}
    ib = {v: i for i, v in enumerate(sorted(set(b.tolist())))}
    c = np.zeros((len(ia), len(ib)), int)          # contingency table of the two divisions
    for x, y in zip(a.tolist(), b.tolist()):
        c[ia[x], ib[y]] += 1
    s = comb(c, 2).sum()
    sa, sb = comb(c.sum(1), 2).sum(), comb(c.sum(0), 2).sum()
    e = sa * sb / comb(len(a), 2)
    denom = (sa + sb) / 2 - e
    if denom == 0:
        raise ValueError("the adjusted Rand index is undefined when both divisions are a single group "
                         "or both put every assemblage alone")
    return float((s - e) / denom)


def pooled_profiles(counts: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Percent of each group's pooled sherds in each class (groups x classes).

    Parameters
    ----------
    counts : array_like, shape (n, classes)
    labels : array_like, shape (n,)

    Returns
    -------
    numpy.ndarray, shape (groups, classes)
        Rows in the sorted order of the distinct labels, each summing to 100.
    """
    counts, labels = np.asarray(counts, float), np.asarray(labels)
    pooled = np.array([counts[labels == g].sum(0) for g in np.unique(labels)])
    return 100.0 * pooled / pooled.sum(1, keepdims=True)


def profile_difference(p: np.ndarray, q: np.ndarray) -> float:
    """Half the summed absolute difference of two percentage profiles (0 to 100).

    The share of sherds that would have to change class to make the two
    profiles identical.

    Parameters
    ----------
    p, q : array_like, shape (classes,)
        Two percentage profiles, each summing to 100.

    Returns
    -------
    float
        0 for identical profiles, 100 for profiles with no class in common.
    """
    return float(0.5 * np.abs(np.asarray(p, float) - np.asarray(q, float)).sum())
