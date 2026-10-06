"""The measures: cultural F_ST, the boundary excess, agreement, and profile difference.

Sign conventions, defined once here and used oriented everywhere:

  boundary excess   within-group similarity MINUS between-group similarity, at
                    matched distance. Positive means assemblages on the same
                    side of a line are more alike than assemblages across it.
  profile difference   always non-negative: the share of sherds that would have
                    to change class to make two profiles identical.
"""
from __future__ import annotations

import numpy as np
from scipy.special import comb

BOUNDARY_EXCESS_MIN_PAIRS = 5


def gini_simpson(counts: np.ndarray) -> float:
    """Gini-Simpson diversity, 1 - sum(p^2), of one row of counts (plug-in)."""
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
    """Cultural F_ST between the groups that `labels` assigns the assemblages to."""
    counts, labels = np.asarray(counts, float), np.asarray(labels)
    return cultural_fst(np.array([counts[labels == g].sum(0) for g in np.unique(labels)]))


def similarity_matrix(counts: np.ndarray) -> np.ndarray:
    """Brainerd-Robinson similarity between every pair of assemblages, 0 to 200."""
    c = np.asarray(counts, float)
    t = c.sum(1, keepdims=True)
    if np.any(t <= 0):
        raise ValueError("every assemblage needs at least one sherd")
    p = 100.0 * c / t
    return 200.0 - np.abs(p[:, None, :] - p[None, :, :]).sum(-1)


def boundary_excess(counts: np.ndarray, dist: np.ndarray, labels: np.ndarray, n_bins: int = 4,
                    min_pairs: int = BOUNDARY_EXCESS_MIN_PAIRS) -> float:
    """Within-group minus between-group similarity, at matched distance. See `boundary_excess_detail`."""
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
    """
    if min_pairs < 1:
        raise ValueError("min_pairs must be at least 1")
    labels = np.asarray(labels)
    s_full = similarity_matrix(counts)
    iu = np.triu_indices(s_full.shape[0], k=1)
    d, s = np.asarray(dist, float)[iu], s_full[iu]
    same = labels[iu[0]] == labels[iu[1]]
    if same.sum() == 0 or (~same).sum() == 0:
        return 0.0, 0
    raw = float(s[same].mean() - s[~same].mean())
    if d.max() == d.min():
        return raw, 0
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
    """Adjusted Rand index: 1 for identical divisions, about 0 for unrelated ones."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) != len(b):
        raise ValueError("the two divisions must label the same assemblages")
    ia = {v: i for i, v in enumerate(sorted(set(a.tolist())))}
    ib = {v: i for i, v in enumerate(sorted(set(b.tolist())))}
    c = np.zeros((len(ia), len(ib)), int)
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
    """Percent of each group's pooled sherds in each class (groups x classes)."""
    counts, labels = np.asarray(counts, float), np.asarray(labels)
    pooled = np.array([counts[labels == g].sum(0) for g in np.unique(labels)])
    return 100.0 * pooled / pooled.sum(1, keepdims=True)


def profile_difference(p: np.ndarray, q: np.ndarray) -> float:
    """Half the summed absolute difference of two percentage profiles (0 to 100).

    The share of sherds that would have to change class to make the two
    profiles identical.
    """
    return float(0.5 * np.abs(np.asarray(p, float) - np.asarray(q, float)).sum())
