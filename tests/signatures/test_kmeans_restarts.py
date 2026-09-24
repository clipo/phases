"""The spatial k-means must return the global optimum on the basin layout.

With 12 restarts it did not at k >= 5 (blind re-derivation, 2026-09-23). The
reference optimum is computed with scipy's k-means2, an implementation
independent of the one under test, over many restarts (rule 3).
"""
from pathlib import Path
import numpy as np
import pytest

XY = Path("data/raw/mainfort-pfg-cplXY.txt")


@pytest.mark.skipif(not XY.exists(), reason="coordinate file absent")
@pytest.mark.parametrize("k", [2, 3, 4, 5, 6])
def test_kmeans_labels_reach_global_optimum(k):
    import importlib, sys
    sys.path.insert(0, "analyses")
    from scipy.cluster.vq import kmeans2
    from mls_emergence.signatures.assortativity import _kmeans_labels
    mf = importlib.import_module("make_figures")
    _, coords = mf._load_curated()
    c = coords.to_numpy(float); c = c - c.mean(0)

    def inertia(lab):
        return sum(((c[lab == g] - c[lab == g].mean(0)) ** 2).sum() for g in np.unique(lab))

    best = np.inf
    for s in range(3000):
        _, lab = kmeans2(c, k, minit="++", seed=s)
        if len(np.unique(lab)) == k:
            best = min(best, inertia(lab))
    got = inertia(_kmeans_labels(c, k, seed=7))
    # got > best would mean a local optimum was returned (the defect);
    # got < best would only mean scipy's restarts missed it, not a failure.
    assert got <= best + 1e-9
