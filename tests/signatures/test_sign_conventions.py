"""Every signed statistic the paper reports, fed a record whose direction is known.

`scripts/check_signs.py` compares the signs in the manuscripts with the signs in
the outputs. It cannot see a sign flipped INSIDE a calculation, because that
flip reaches the outputs and the text together. These tests are the other
half: each builds the smallest record on which the direction of a statistic is
not in doubt and asserts the sign the code returns. They pin conventions, not
values. If a statistic is ever redefined with the opposite orientation, the
test for it must change in the same commit as every sentence that quotes it.

Rule 3: each test names what the opposite sign would mean.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from mls_emergence.signatures import assortativity as A  # noqa: E402
from mls_emergence.signatures import variance as V  # noqa: E402

sd = importlib.import_module("23_phases_vs_spatial_drift")


def _gradient_record(n=12):
    """Sites along a line whose composition shifts steadily from one class to
    another: similarity falls with distance, with no groups anywhere."""
    x = np.linspace(0.0, 110.0, n)
    share = np.linspace(0.9, 0.1, n)
    counts = np.column_stack([share, 1 - share]) * 200.0
    coords = np.column_stack([x, np.zeros(n)])
    return counts, coords, np.abs(x[:, None] - x[None, :])


def test_similarity_is_high_for_alike_and_low_for_unalike():
    # Brainerd-Robinson is a SIMILARITY: 200 for identical profiles, 0 for
    # disjoint ones. Read as a distance, every correlation below flips.
    assert A.brainerd_robinson(np.array([30.0, 70.0]), np.array([3.0, 7.0])) == pytest.approx(200.0)
    assert A.brainerd_robinson(np.array([10.0, 0.0]), np.array([0.0, 10.0])) == pytest.approx(0.0)


def test_distance_decay_correlation_is_negative_on_a_gradient():
    # The paper quotes the correlation of SIMILARITY with distance, so decay
    # is a negative number (-0.20 on the basin). A positive value here would
    # mean a dissimilarity was passed, or the fade was reported as its opposite.
    counts, _, dist = _gradient_record()
    r, _ = A.mantel(A.similarity_matrix(counts), dist, n_perm=1)
    assert r < -0.9
    assert sd.stats_for_matrix(counts, dist, seed=0)["dd_r"] == pytest.approx(r)


def test_cultural_fst_rises_with_differentiation_and_is_zero_without_it():
    # F_ST is the share of diversity BETWEEN groups: 0 when groups share one
    # profile, larger as they diverge. A measure that fell as groups diverged
    # would be the within-group share instead.
    same = np.array([[60.0, 40.0], [120.0, 80.0]])
    mild = np.array([[65.0, 35.0], [55.0, 45.0]])
    strong = np.array([[90.0, 10.0], [10.0, 90.0]])
    assert V.cultural_fst(same) == pytest.approx(0.0, abs=1e-12)
    assert 0 < V.cultural_fst(mild) < V.cultural_fst(strong)


def test_boundary_excess_is_positive_when_groups_are_more_alike_inside():
    # Within-group similarity minus between-group similarity at matched
    # distance. Negative would mean pairs across the line are MORE alike than
    # pairs on one side of it. (The record on which one far pair used to set
    # this sign is pinned in test_boundary_excess_labeled.py.)
    x = np.array([0.0, 10, 20, 30, 40, 50, 5, 15, 25, 35, 45, 55])
    labels = np.array([0] * 6 + [1] * 6)
    counts = np.where(labels[:, None] == 0, [80.0, 20.0], [60.0, 40.0])
    dist = np.abs(x[:, None] - x[None, :])
    assert sd.boundary_excess_labeled(counts, dist, labels) > 0
    assert sd.boundary_excess_labeled(counts, dist, 1 - labels) > 0, "relabeling the groups cannot change the sign"
    # ...and it is about zero when composition follows distance alone.
    g_counts, _, g_dist = _gradient_record()
    halves = (np.arange(len(g_counts)) >= len(g_counts) // 2).astype(int)
    interleaved = np.arange(len(g_counts)) % 2
    assert abs(sd.boundary_excess_labeled(g_counts, g_dist, interleaved)) < \
        sd.boundary_excess_labeled(counts, dist, labels)
    assert np.isfinite(sd.boundary_excess_labeled(g_counts, g_dist, halves))


@pytest.mark.data
def test_seriation_axis_is_oriented_later_is_higher_and_north_is_higher():
    # The axis's own sign is arbitrary; the project fixes it so that the
    # pooled radiocarbon medians of the dated assemblages do not DECREASE
    # along it, and under that orientation its rank correlation with latitude
    # is the +0.77 the main text quotes. A negative value here means the axis
    # was used unoriented, and every "along the sequence" trend has the
    # opposite sign from the text.
    mf = importlib.import_module("make_figures")
    basin = importlib.import_module("17_basin_results")
    counts, coords = mf._load_curated()
    ca, n_anchors = basin.oriented_ca(counts)
    assert n_anchors >= 3, "too few dated assemblages to orient the axis"
    lat = coords.reindex(ca.index).to_numpy(float)[:, 0]
    rho = spearmanr(ca.to_numpy(float), lat)[0]
    assert rho == pytest.approx(0.77, abs=0.01)
