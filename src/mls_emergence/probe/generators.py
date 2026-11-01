"""Copying-network families for the probe's generator (D21, 2026-10-16).

The inherited generator weights copying between sites by a distance kernel
with a free interaction length (exponential at 24 km, a number the project's
decision register says "is not derived from anything", D-18). Under rule 19
the length is re-screened rather than inherited, and a second family is
added: copying restricted to Gabriel neighbors (two sites are neighbors if no
third site lies inside the circle on their connecting segment), which has no
free length and adapts to settlement density; it is the adjacency the
coherence readout uses (D14). Both families produce row-stochastic weight
matrices that `drift_record` accepts, and a barrier multiplies cross-label
weights by `leak` exactly as `copying_weights` does.
"""
from __future__ import annotations

import numpy as np

from mls_emergence.transmission.spatial import copying_weights


def gabriel_edges(D: np.ndarray) -> np.ndarray:
    """Gabriel graph edges from a pairwise distance matrix, as an (m, 2) array of index pairs (i < j).

    i and j are neighbors when no k has d(i,k)^2 + d(j,k)^2 < d(i,j)^2.
    Computed on the distance matrix so the model and the generator share one
    geometry (D10).
    """
    D = np.asarray(D, float)
    n = len(D)
    if D.ndim != 2 or D.shape != (n, n) or not np.allclose(D, D.T) or np.any(D < 0):
        raise ValueError("D must be a symmetric non-negative distance matrix")
    edges = []
    D2 = D ** 2
    for i in range(n):
        for j in range(i + 1, n):
            s = D2[i] + D2[j]                         # d(i,k)^2 + d(j,k)^2 for every k
            s[i] = s[j] = np.inf
            if not np.any(s < D2[i, j]):
                edges.append((i, j))
    return np.array(edges, dtype=int).reshape(-1, 2)


def graph_weights(edges: np.ndarray, n: int, labels=None, leak: float = 1.0) -> np.ndarray:
    """Row-stochastic copying weights with equal weight on each graph edge.

    A site copies from each of its graph neighbors with equal weight; a
    barrier multiplies cross-label edges by `leak` (0 <= leak <= 1). A site
    whose every edge is cut keeps itself as its only source, as
    `copying_weights` does for an isolated node.
    """
    edges = np.asarray(edges, int)
    if edges.ndim != 2 or edges.shape[1] != 2 or (edges < 0).any() or (edges >= n).any():
        raise ValueError("edges must be an (m, 2) array of indices below n")
    if not 0 <= leak <= 1:
        raise ValueError("leak must lie in [0, 1]")
    w = np.zeros((n, n))
    for i, j in edges:
        w[i, j] = w[j, i] = 1.0
    if labels is not None:
        labels = np.asarray(labels)
        if labels.shape != (n,):
            raise ValueError("one label is required per node")
        w *= np.where(labels[:, None] == labels[None, :], 1.0, leak)
    sums = w.sum(1)
    isolated = sums == 0
    w[np.where(isolated)[0], np.where(isolated)[0]] = 1.0
    return w / w.sum(1, keepdims=True)


FAMILIES = (
    ("exponential", 6.0), ("exponential", 12.0), ("exponential", 24.0),
    ("gaussian", 6.0), ("gaussian", 12.0), ("gaussian", 24.0),
    ("power", 6.0), ("power", 12.0), ("power", 24.0),
    ("gabriel", None),
)


def family_name(kernel: str, length) -> str:
    return "gabriel" if kernel == "gabriel" else f"{kernel}-{length:g}km"


def family_weights(kernel: str, length, D: np.ndarray, labels=None, leak: float = 1.0) -> np.ndarray:
    """Weights for one family: a distance kernel at `length` km, or the Gabriel graph."""
    if kernel == "gabriel":
        return graph_weights(gabriel_edges(D), len(D), labels=labels, leak=leak)
    return copying_weights(D, length=float(length), labels=labels, leak=leak, kernel=kernel)
