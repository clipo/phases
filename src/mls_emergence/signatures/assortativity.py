"""Spatial-assortativity signature: does decorated similarity follow social
rather than geographic distance?

Provides Brainerd-Robinson similarity, a Mantel test against geographic
distance, and the distance-controlled boundary-excess statistic that separates a
sharp social boundary from smooth isolation by distance (a plain Mantel r cannot
tell the two apart). Spatial clusters for the boundary excess are found by
k-means with multiple restarts.
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
    """Pairwise Brainerd-Robinson similarity matrix (n x n) over assemblages."""
    n = counts.shape[0]
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            S[i, j] = brainerd_robinson(counts[i], counts[j])
    return S

def geo_distance(coords: np.ndarray) -> np.ndarray:
    """Pairwise Euclidean distance matrix for the given coordinates."""
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
    """One Lloyd k-means run from a random seeding; returns (labels, inertia)."""
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


def _kmeans_labels(coords: np.ndarray, k: int, seed: int = 0, n_init: int = 12) -> np.ndarray:
    """k-means (Lloyd) with multiple restarts; returns the lowest-inertia labeling.

    Restarts make the spatial partition stable: with a single fixed init, k-means
    on a clustered layout can land in an orthogonal local optimum and invert the
    within/between-cluster similarity gap. Taking the best of several inits
    removes that instability.
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
