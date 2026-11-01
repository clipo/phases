"""Signature 4, spatial assortativity of ceramic similarity.

Provides Brainerd-Robinson similarity (``brainerd_robinson``,
``similarity_matrix``), pairwise Euclidean distance (``geo_distance``), the
Mantel correlation between the two (``mantel``, ``spatial_assortativity``), the
k-means spatial clustering used to define groups (``_kmeans_labels``,
``_select_k``) and the distance-controlled ``boundary_excess``. Used by the
empirical and phase analyses (05, 06, 07, 09, 15, 23, 25, 26, 28, 43, 47, 59, 60),
``make_figures.py`` and ``validation.harness``.
"""
from __future__ import annotations

import numpy as np


def brainerd_robinson(a: np.ndarray, b: np.ndarray) -> float:
    """BR similarity in [0,200]: 200 - sum|%a-%b|."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    pa = 100 * a / a.sum() if a.sum() else a * 0
    pb = 100 * b / b.sum() if b.sum() else b * 0
    return float(200.0 - np.sum(np.abs(pa - pb)))

def similarity_matrix(counts: np.ndarray) -> np.ndarray:
    """Pairwise Brainerd-Robinson similarity of the rows of ``counts``.

    ``counts`` is an assemblages x types array. Returns a symmetric n x n array
    in [0, 200]; 200 means identical percentage profiles.
    """
    n = counts.shape[0]
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            S[i, j] = brainerd_robinson(counts[i], counts[j])
    return S

def geo_distance(coords: np.ndarray) -> np.ndarray:
    """Pairwise Euclidean distance between the rows of ``coords`` (n x 2).

    Returns an n x n array in the units of ``coords`` (meters for projected
    UTM input; convert before passing if km are wanted). No geodesic
    correction is made, so latitude/longitude input gives degrees, not a
    distance.
    """
    d = coords[:, None, :] - coords[None, :, :]
    return np.sqrt((d ** 2).sum(-1))

def mantel(sim: np.ndarray, dist: np.ndarray, n_perm: int = 999, seed: int = 0):
    """Mantel test between a similarity matrix and a distance matrix."""
    iu = np.triu_indices_from(sim, k=1)
    s, d = sim[iu], dist[iu]
    r_obs = np.corrcoef(s, d)[0, 1]
    rng = np.random.default_rng(seed)
    n = sim.shape[0]; count = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        sp = sim[p][:, p][iu]
        if abs(np.corrcoef(sp, d)[0, 1]) >= abs(r_obs): count += 1
    return float(r_obs), float((count + 1) / (n_perm + 1))

def spatial_assortativity(counts: np.ndarray, coords: np.ndarray, n_perm: int = 999) -> dict:
    """Mantel r vs geography; boundary_residual = 1 - r^2 (structure beyond IBD)."""
    sim, dist = similarity_matrix(counts), geo_distance(coords)
    r, p = mantel(sim, dist, n_perm)
    return {"mantel_r": r, "mantel_p": p, "boundary_residual": float(1 - r ** 2)}


def _kmeans_once(coords: np.ndarray, k: int, rng: np.random.Generator):
    n = coords.shape[0]
    centers = coords[rng.choice(n, size=k, replace=False)].copy()
    labels = np.full(n, -1)
    for _ in range(50):
        d = ((coords[:, None, :] - centers[None, :, :]) ** 2).sum(-1)
        new = d.argmin(axis=1)
        if np.array_equal(new, labels):
            break
        labels = new
        for j in range(k):
            members = coords[labels == j]
            if len(members):
                centers[j] = members.mean(axis=0)
    inertia = float(
        sum(((coords[labels == j] - centers[j]) ** 2).sum() for j in range(k))
    )
    return labels, inertia


def _kmeans_labels(coords: np.ndarray, k: int, seed: int = 0, n_init: int = 500) -> np.ndarray:
    """k-means (Lloyd) with multiple restarts; returns the lowest-inertia labeling.

    Restarts make the spatial partition stable: with a single fixed init, k-means
    on a clustered layout can land in an orthogonal local optimum and invert the
    within/between-cluster similarity gap. Taking the best of several inits
    removes that instability.

    n_init was 12 until 2026-09-23. A blind re-derivation found that 12 is
    enough on the 28 basin assemblages at k = 2 to 4 but not at k >= 5, where
    the returned partition was a worse local optimum (inertia 0.2924 against
    0.2768 at k = 5, changing the between-cluster F_ST from 0.0307 to 0.0206).
    500 restarts on 28 points is cheap and returns the best-of-3,000 optimum
    at k = 2 to 7 (tests/signatures/test_kmeans_restarts.py).
    """
    coords = np.asarray(coords, float)
    n = coords.shape[0]
    k = max(1, min(k, n))
    rng = np.random.default_rng(seed)
    best_labels, best_inertia = None, np.inf
    for _ in range(n_init):
        labels, inertia = _kmeans_once(coords, k, rng)
        if inertia < best_inertia:
            best_inertia, best_labels = inertia, labels
    return best_labels


def _select_k(coords: np.ndarray, k_range=(2, 3, 4), seed: int = 0) -> int:
    """Pick the number of spatial clusters by a simple silhouette-like score.

    For each candidate k, score = (mean nearest-other-cluster centroid distance)
    relative to (mean within-cluster spread). The k with the cleanest separation
    wins. This keeps boundary_excess from depending on a hard-coded cluster count
    (e.g. k=2 for a two-blob layout, k=3 for three clusters).
    """
    coords = np.asarray(coords, float)
    n = coords.shape[0]
    best_k, best_score = k_range[0], -np.inf
    for k in k_range:
        if k >= n:
            continue
        labels = _kmeans_labels(coords, k, seed=seed)
        if len(np.unique(labels)) < k:
            continue
        centers = np.array([coords[labels == j].mean(axis=0) for j in range(k)])
        within = np.mean(
            [np.linalg.norm(coords[i] - centers[labels[i]]) for i in range(n)]
        )
        cc = geo_distance(centers)
        np.fill_diagonal(cc, np.inf)
        between = cc.min(axis=1).mean()
        score = between / (within + 1e-9)
        if score > best_score:
            best_score, best_k = score, k
    return best_k


def boundary_excess(
    counts: np.ndarray,
    coords: np.ndarray,
    n_clusters: int | None = None,
    n_bins: int = 4,
    seed: int = 0,
) -> float:
    """Distance-controlled within- minus between-cluster similarity (BR units).

    A plain Mantel r cannot separate a SHARP spatial boundary from SMOOTH
    isolation-by-distance: both yield negative r. This metric does. It clusters
    the coordinates, then within each pairwise-distance bin compares the mean
    Brainerd-Robinson similarity of WITHIN-cluster pairs against BETWEEN-cluster
    pairs at the SAME distance. Under pure IBD, similarity depends only on
    distance, so within and between are comparable and the excess is ~0. A sharp
    boundary makes within-cluster pairs far more similar than between-cluster
    pairs at matched distance, yielding a large positive excess.

    Returns the distance-bin-averaged (within - between) BR similarity gap.

    LIMIT (measured 2026-10-02, scripts/probe_library_boundary_excess.py). A
    bin counts as soon as it holds one within-cluster and one between-cluster
    pair. On the per-period windows of four to ten assemblages this function
    is called on, most counted bins hold one or two pairs, and requiring five
    would change the sign of one call in four to six. The value on such a
    window is not a resolved estimate. A minimum-pairs rule is not added here
    because it would empty the bins and return the raw gap below, a different
    statistic; the phase-line statistic the paper reports
    (`23_phases_vs_spatial_drift.boundary_excess_labeled`) has the rule.
    """
    counts = np.asarray(counts, float)
    coords = np.asarray(coords, float)
    n = counts.shape[0]
    if n < 4:
        return 0.0

    if n_clusters is None:
        n_clusters = _select_k(coords, seed=seed)
    labels = _kmeans_labels(coords, n_clusters, seed=seed)
    sim = similarity_matrix(counts)
    dist = geo_distance(coords)

    iu = np.triu_indices(n, k=1)
    d = dist[iu]
    s = sim[iu]
    same = labels[iu[0]] == labels[iu[1]]

    if same.sum() == 0 or (~same).sum() == 0:
        return 0.0
    if d.max() == d.min():
        return float(s[same].mean() - s[~same].mean())

    edges = np.linspace(d.min(), d.max() + 1e-9, n_bins + 1)
    gaps = []
    n_overlap = 0
    for b in range(n_bins):
        in_bin = (d >= edges[b]) & (d < edges[b + 1])
        w = in_bin & same
        btw = in_bin & ~same
        if w.sum() == 0 or btw.sum() == 0:
            continue
        n_overlap += 1
        gaps.append(float(s[w].mean() - s[btw].mean()))

    raw_gap = float(s[same].mean() - s[~same].mean())
    if not gaps:
        # No distance bin contains both within- and between-cluster pairs: the
        # clusters are separated in distance, which is itself the signature of a
        # sharp spatial boundary. Report the raw within-between gap.
        return raw_gap
    # When clusters share distance ranges (as under a smooth IBD gradient), the
    # distance-controlled gap collapses toward zero even though the raw gap is
    # large; that is exactly the sharp-vs-smooth discrimination we want.
    return float(np.mean(gaps))


# ---------------------------------------------------------------------------
# Pair statistics for the rim R1 readings (spec 2026-10-19-rim-attributes-r1, 4.2 and 4.3)
# ---------------------------------------------------------------------------

def total_variation(P: np.ndarray) -> np.ndarray:
    """Pairwise total variation between composition rows: 0.5 * sum_c |p_ic - p_jc|.

    Zero for identical compositions and equal to |p_i - p_j| for a binary
    trait. Not 1 - p.q, which adds half the within-site diversities of the
    two rows to half their squared difference and can reorder pairs (Codex
    RIM-1). Rows must be probability vectors.
    """
    P = np.asarray(P, float)
    if P.ndim != 2 or not np.isfinite(P).all() or (P < 0).any() or not np.allclose(P.sum(1), 1.0, atol=1e-9, rtol=0):
        raise ValueError("total_variation takes rows that are probability vectors")
    return 0.5 * np.abs(P[:, None, :] - P[None, :, :]).sum(-1)


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation; NaN when either vector is constant."""
    from scipy.stats import rankdata
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("spearman takes two vectors of the same length")
    ra, rb = rankdata(a), rankdata(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def partial_spearman(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Association of a with b after LINEAR adjustment on ranks for c: the correlation of the rank residuals of a and b
    after least-squares removal of the centered rank of c.

    This is the form analysis 100 used. It is not "at a given distance": a
    nonlinear dependence on c survives the linear adjustment (two quadratic
    functions of c return 1; tests). NaN when a vector is constant or a
    residual vanishes; when c is constant the plain Spearman is returned.
    """
    from scipy.stats import rankdata
    a = np.asarray(a, float); b = np.asarray(b, float); c = np.asarray(c, float)
    if not (a.shape == b.shape == c.shape) or a.ndim != 1:
        raise ValueError("partial_spearman takes three vectors of the same length")
    ra, rb, rc = (rankdata(v) - (len(v) + 1) / 2.0 for v in (a, b, c))
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    if rc.std() == 0:
        return spearman(a, b)
    ea = ra - rc * (ra @ rc) / (rc @ rc)
    eb = rb - rc * (rb @ rc) / (rc @ rc)
    if ea.std() < 1e-12 or eb.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(ea, eb)[0, 1])
