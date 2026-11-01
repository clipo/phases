"""Rival-cause generators for the per-type step (spec 6.2, 6.4): C2 totals by span length, R1 misclassification, R2 vessel clumps."""
import numpy as np
import pytest

from mls_emergence.probe import rivals


def test_totals_by_span_length_is_a_permutation_ordered_by_length():
    tot = np.array([79, 1709, 251, 300]); L = np.array([5, 50, 20, 20])
    out = rivals.totals_by_span_length(tot, L)
    assert sorted(out) == sorted(tot) and out[1] == 1709 and out[0] == 79


def test_misclassify_preserves_row_sums_and_moves_mass_at_rate_one():
    y = np.array([[10, 0], [0, 10]]); pooled = np.array([0.5, 0.5])
    out = rivals.misclassify(y, 1.0, pooled, np.random.default_rng(1))
    assert (out.sum(1) == 10).all() and out[0, 1] > 0
    big = np.array([[100000, 0]]); moved = rivals.misclassify(big, 0.01, pooled, np.random.default_rng(3))
    assert 400 < moved[0, 1] < 600          # rate 0.01 of 100,000 sherds, half of which land in class 2
    assert (rivals.misclassify(y, 0.0, pooled, np.random.default_rng(1)) == y).all()
    with pytest.raises(ValueError):
        rivals.misclassify(y, 1.5, pooled, np.random.default_rng(1))


def test_clumped_counts_meets_totals_exactly_and_truncates():
    p = np.array([[0.9, 0.1]])
    out = rivals.clumped_counts(p, np.array([3]), mean_clump=10.0, rng=np.random.default_rng(2))
    assert out.sum() == 3 and (out > 0).sum() == 1      # one clump, one class, truncated to the total
    many = rivals.clumped_counts(np.array([[0.5, 0.5]]), np.array([100000]), mean_clump=5.0, rng=np.random.default_rng(4))
    assert many.sum() == 100000 and abs(many[0, 0] - 50000) < 2500   # clumps of mean 5 inflate the SD about sqrt(5)-fold
    with pytest.raises(ValueError):
        rivals.clumped_counts(p, np.array([3]), mean_clump=0.5, rng=np.random.default_rng(2))


def test_plant_marker_touches_only_the_cluster():
    """M (regional-drift spec 4.1): class c's odds times m in the marked cluster's rows only; m = 1 is the identity."""
    prof = np.array([[.5, .3, .2], [.1, .6, .3], [.25, .25, .5], [.7, .2, .1]])
    labels = np.array([0, 1, 1, 0])
    out = rivals.plant_marker(prof, labels, cluster=1, c=1, m=3.0)
    assert np.array_equal(out[[0, 3]], prof[[0, 3]])                                   # unmarked rows untouched
    for i in (1, 2):
        p = prof[i, 1]
        assert np.isclose(out[i, 1], 3 * p / (1 + 2 * p)) and np.isclose(out[i].sum(), 1.0)   # odds multiplier, renormalized
        assert np.allclose(out[i, [0, 2]], prof[i, [0, 2]] / (1 + 2 * p))
    assert np.array_equal(rivals.plant_marker(prof, labels, 1, 1, 1.0), prof)
    assert np.array_equal(prof, np.array([[.5, .3, .2], [.1, .6, .3], [.25, .25, .5], [.7, .2, .1]]))   # input not modified
    for bad_m in (0.0, -2.0, np.nan, np.inf):
        with pytest.raises(ValueError):
            rivals.plant_marker(prof, labels, 1, 1, bad_m)
    with pytest.raises(ValueError):
        rivals.plant_marker(prof, labels, 7, 1, 2.0)        # absent cluster label
    with pytest.raises(ValueError):
        rivals.plant_marker(prof, labels, 1, 3, 2.0)        # class out of range
    with pytest.raises(ValueError):
        rivals.plant_marker(prof * 1.01, labels, 1, 1, 2.0)  # rows do not sum to 1


def test_plant_marker_at_m_one_is_bitwise_identity_on_normalized_rows():
    """Rows normalized by division sum to 1 only within an ulp; m = 1 must return them untouched (no renormalization)."""
    rng = np.random.default_rng(11)
    raw = rng.dirichlet(np.ones(10) * 0.3, size=5000); prof = raw / raw.sum(1, keepdims=True)
    labels = rng.integers(0, 2, size=5000)
    out = rivals.plant_marker(prof, labels, 1, 3, 1.0)
    assert np.array_equal(out, prof)                         # a renormalizing implementation changes about 90 percent of rows by an ulp
    with pytest.raises(ValueError):
        rivals.plant_marker(prof[:4], labels[:4], 1, True, 2.0)   # a bool is not a class index
