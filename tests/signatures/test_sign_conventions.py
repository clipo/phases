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

from mls_emergence.signatures import assortativity as A
from mls_emergence.signatures import variance as V

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


def test_rim_partial_rank_correlation_keeps_direction_and_removes_a_shared_cause():
    """Analysis 100, part D: decoration difference against attribute difference, distance removed.

    Positive means pairs that differ more in the attribute also differ more in
    decoration at a given distance. A negative value would mean they differ
    LESS, and reading the paper's +0.06 with the other sign would turn "age
    does not go with the difference" into "older pairs are more alike".
    """
    ra = importlib.import_module("100_rim_attributes")
    rng = np.random.default_rng(100)
    dist = rng.uniform(0, 100, 400)
    own = rng.normal(0, 30, 400)
    # the attribute carries its own signal into decoration: positive with distance removed
    attr = dist + own
    assert ra.partial_spearman(dist + own + rng.normal(0, 5, 400), attr, dist) > 0.5
    # reversed signal: negative
    assert ra.partial_spearman(dist - own + rng.normal(0, 5, 400), attr, dist) < -0.5
    # both follow distance and nothing else: positive as it stands, near zero with distance removed.
    # A rival that did not remove distance would report the first number for the second.
    a, b = dist + rng.normal(0, 10, 400), dist + rng.normal(0, 10, 400)
    assert ra.spearman(a, b) > 0.7
    assert abs(ra.partial_spearman(a, b, dist)) < 0.15


def test_per_class_share_size_correlation_sign():
    # Per-type step (spec 4.4 item 5; findings per_type_generators_*.md): rho_c is the Spearman correlation of a
    # class's share with the site's sherd total. Positive means the share RISES with collection size (the record's
    # Parkin_Punctated, +0.42); negative means it falls (Painted, -0.57). A flipped sign would read the common class
    # as concentrated in small collections.
    for _q in (ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups"):
        if str(_q) not in sys.path:
            sys.path.insert(0, str(_q))
    pct = importlib.import_module("per_class_table")
    n = np.array([100, 200, 400, 800]); share = np.array([0.1, 0.2, 0.3, 0.4])
    y = np.column_stack([np.rint(n * share), n - np.rint(n * share)]).astype(int)
    row = pct.per_class_rows(y, labels=np.array([0, 0, 1, 1]))
    assert row["rho_0"] > 0 and row["rho_1"] < 0


def test_collection_log_odds_ratio_sign():
    # Composition by collection (D36 closing step; findings composition_by_collection.md): the log odds ratio is the Lipo
    # addition MINUS the PFG part, so positive means the class is a larger share of Lipo's 1996-97 collection than of the
    # 1940 collection at the same site. A flipped sign would attribute to the survey what belongs to the later collection.
    for _q in (ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups", ROOT / "scripts"):
        if str(_q) not in sys.path:
            sys.path.insert(0, str(_q))
    cbc = importlib.import_module("composition_by_collection")
    r = cbc.log_odds_ratio(np.array([80, 20]), np.array([20, 80]), np.random.default_rng(1), draws=4000)
    assert r.loc[0, "median"] > 0 and r.loc[1, "median"] < 0
