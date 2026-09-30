import numpy as np
from mls_emergence.signatures import assortativity as A

def test_brainerd_robinson_identical_is_200():
    a = np.array([50, 30, 20]); b = np.array([50, 30, 20])
    assert abs(A.brainerd_robinson(a, b) - 200.0) < 1e-9

def test_brainerd_robinson_disjoint_is_0():
    a = np.array([100, 0]); b = np.array([0, 100])
    assert abs(A.brainerd_robinson(a, b)) < 1e-9

def test_mantel_detects_distance_decay():
    rng = np.random.default_rng(1)
    coords = np.vstack([rng.normal([0,0],1,(10,2)), rng.normal([100,100],1,(10,2))])
    counts = np.vstack([np.tile([80,20,0],(10,1)), np.tile([0,20,80],(10,1))]).astype(float)
    sim = A.similarity_matrix(counts)
    r, p = A.mantel(sim, A.geo_distance(coords))
    assert r < 0 and p < 0.05

def test_boundary_excess_positive_for_sharp_clusters():
    """Sharp spatial clusters: high within-cluster similarity, low between =>
    large positive boundary excess controlling for distance."""
    rng = np.random.default_rng(2)
    # Two tight spatial clusters with distinct, internally-uniform assemblages.
    coords = np.vstack([rng.normal([0, 0], 0.3, (8, 2)),
                        rng.normal([20, 20], 0.3, (8, 2))])
    counts = np.vstack([np.tile([90, 10, 0, 0], (8, 1)),
                        np.tile([0, 0, 10, 90], (8, 1))]).astype(float)
    excess = A.boundary_excess(counts, coords)
    assert excess > 20.0  # strong sharp boundary on the 0-200 BR scale


def test_boundary_excess_near_zero_for_smooth_ibd():
    """Smooth isolation-by-distance gradient along a transect has no sharp edge,
    so within- vs between-cluster similarity at matched distance is comparable
    and the boundary excess stays small relative to a sharp boundary."""
    n = 16
    coords = np.column_stack([np.linspace(0, 20, n), np.zeros(n)])
    # Profiles interpolate smoothly between two endpoints; no discontinuity.
    frac = np.linspace(0, 1, n)
    counts = np.column_stack([
        100 * (1 - frac), 100 * frac, np.full(n, 20.0)
    ])
    sharp_coords = np.vstack([np.zeros((8, 2)), np.full((8, 2), 20.0)])
    sharp_counts = np.vstack([np.tile([90, 0, 20], (8, 1)),
                            np.tile([0, 90, 20], (8, 1))]).astype(float)
    ibd_excess = A.boundary_excess(counts, coords)
    sharp_excess = A.boundary_excess(sharp_counts, sharp_coords)
    assert ibd_excess < sharp_excess
