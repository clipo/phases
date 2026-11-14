"""The two spatial readouts of the drift-extension design (spec 2026-11-14, draft 3, sections 4 and 8).

Expected values are worked by hand in the comments. Direction tests are in
tests/signatures/test_sign_conventions.py (rule 21a).
"""
import numpy as np
import pytest
from scipy.stats import pearsonr

from mls_emergence.signatures import spatial_readouts as sr

# Four sites, ten classes, ten sherds each. Shares: s0 = all class 0; s1 = half class 0, half
# class 1; s2 = all class 2; s3 = half class 0, half class 2.
Y4 = np.zeros((4, 10), int)
Y4[0, 0] = 10; Y4[1, 0] = 5; Y4[1, 1] = 5; Y4[2, 2] = 10; Y4[3, 0] = 5; Y4[3, 2] = 5


def test_tv_pairs_by_hand_with_the_factor_of_one_half_and_unique_pairs():
    # TV = 0.5 * sum |p_i - p_j|, pairs in the order (0,1) (0,2) (0,3) (1,2) (1,3) (2,3):
    # (0,1) 0.5(0.5 + 0.5) = 0.5; (0,2) 0.5(1 + 1) = 1; (0,3) 0.5(0.5 + 0.5) = 0.5;
    # (1,2) 0.5(0.5 + 0.5 + 1) = 1; (1,3) 0.5(0 + 0.5 + 0.5) = 0.5; (2,3) 0.5(0.5 + 0.5) = 0.5.
    tv = sr.tv_pairs(Y4)
    assert tv.shape == (6,)                              # n(n-1)/2 unique pairs, not n^2 or n(n-1)
    np.testing.assert_allclose(tv, [0.5, 1.0, 0.5, 1.0, 0.5, 0.5])
    # Rivals rejected: the sum without the half (values doubled); all ordered pairs (length 12).
    # A correlation cannot detect either, which is why the array is tested directly.


def test_lat_rank_is_spearman_not_pearson():
    # Share of class 0 rises with latitude monotonically but not linearly: 0.1, 0.2, 0.3, 0.9.
    y = np.array([[1, 9], [2, 8], [3, 7], [9, 1]])
    lat = np.array([1.0, 2.0, 3.0, 4.0])
    val, reason = sr.lat_rank(y, lat)
    assert val[0] == pytest.approx(1.0) and val[1] == pytest.approx(-1.0) and reason == ["", ""]
    assert pearsonr(y[:, 0] / 10, lat)[0] < 0.95           # the Pearson rival does not give 1


def test_lat_rank_average_ranks_for_ties_by_hand():
    # Shares of class 0: [0.1, 0.1, 0.3, 0.5] -> average ranks [1.5, 1.5, 3, 4]; latitude ranks
    # [1, 2, 3, 4]. Deviations from 2.5: x = [-1, -1, 0.5, 1.5], y = [-1.5, -0.5, 0.5, 1.5].
    # Sum xy = 1.5 + 0.5 + 0.25 + 2.25 = 4.5; sum xx = 4.5; sum yy = 5.
    # r = 4.5 / sqrt(4.5 * 5) = 0.9486833.
    y = np.array([[1, 9], [1, 9], [3, 7], [5, 5]])
    val, _ = sr.lat_rank(y, np.array([1.0, 2.0, 3.0, 4.0]))
    assert val[0] == pytest.approx(0.9486833, abs=1e-6)
    assert val[1] == pytest.approx(-0.9486833, abs=1e-6)


def test_absent_at_one_site_keeps_the_site_with_share_zero():
    # Class 0 shares [0, 0.2, 0.4, 0.6] rise strictly with latitude: +1 only if the zero is kept.
    y = np.array([[0, 10], [2, 8], [4, 6], [6, 4]])
    val, reason = sr.lat_rank(y, np.array([1.0, 2.0, 3.0, 4.0]))
    assert val[0] == pytest.approx(1.0) and reason[0] == ""


def test_absent_at_an_interior_site_changes_the_value_if_the_site_is_dropped():
    # The monotone fixture above gives +1 with or without the absent site, so it cannot reject a
    # function that drops sites where the class is absent (Codex build audit 2026-11-14, F5).
    # Here the absent site is interior. Class 0 shares [0.3, 0, 0.2, 0.6] at latitudes 1, 2, 3, 4.
    # Kept: share ranks [3, 1, 2, 4] against [1, 2, 3, 4]; d = [2, -1, -1, 0], sum d^2 = 6;
    # rho = 1 - 6 * 6 / (4 * 15) = 0.4.
    # Dropped: shares [0.3, 0.2, 0.6] at latitudes [1, 3, 4]; ranks [2, 1, 3] against [1, 2, 3];
    # sum d^2 = 2; rho = 1 - 6 * 2 / (3 * 8) = 0.5.
    y = np.array([[3, 7], [0, 10], [2, 8], [6, 4]])
    val, reason = sr.lat_rank(y, np.array([1.0, 2.0, 3.0, 4.0]))
    assert val[0] == pytest.approx(0.4, abs=1e-9) and reason[0] == ""


def test_constant_inputs_return_nan_with_reason():
    lat = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([[0, 5, 5], [0, 5, 5], [0, 5, 5], [0, 5, 5]])   # class 0 absent everywhere; all sites identical
    val, reason = sr.lat_rank(y, lat)
    assert np.isnan(val).all() and reason == ["constant"] * 3
    a = sr.a_lat(y, lat)
    assert np.isnan(a).all()
    D = np.abs(lat[:, None] - lat[None, :])
    d, why = sr.distance_decay_rank(y, D)
    assert np.isnan(d) and why == "constant"


def test_a_lat_is_the_absolute_value_of_lat_rank():
    y = np.array([[1, 9], [2, 8], [3, 7], [9, 1]])
    lat = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(sr.a_lat(y, lat), [1.0, 1.0])
    assert sr.lat_rank(y, lat)[0][1] < 0                      # the signed value is kept signed


def test_distance_decay_uses_unique_pairs_of_the_distance_matrix():
    # Y4 with distances chosen so that the two TV = 1 pairs, (0,2) and (1,2), are the two longest.
    D = np.array([[0, 1, 5, 2], [1, 0, 6, 3], [5, 6, 0, 4], [2, 3, 4, 0]], float)
    # TV = [0.5, 1, 0.5, 1, 0.5, 0.5] -> average ranks [2.5, 5.5, 2.5, 5.5, 2.5, 2.5];
    # d  = [1, 5, 2, 6, 3, 4] -> ranks [1, 5, 2, 6, 3, 4]. Deviations from 3.5:
    # x = [-1, 2, -1, 2, -1, -1], y = [-2.5, 1.5, -1.5, 2.5, -0.5, 0.5].
    # Sum xy = 2.5 + 3 + 1.5 + 5 + 0.5 - 0.5 = 12; sum xx = 12; sum yy = 17.5.
    # r = 12 / sqrt(12 * 17.5) = 0.8280787.
    val, reason = sr.distance_decay_rank(Y4, D)
    assert val == pytest.approx(0.8280787, abs=1e-6) and reason == ""
    # Two rivals this fixture cannot reject (Codex build audit 2026-11-14, F5). Pearson: the
    # distances 1..6 are their own ranks and TV takes two values, so Pearson equals Spearman here;
    # the next test separates them. Counting every pair twice: duplicating all pairs leaves any
    # correlation unchanged, so no correlation value can detect it. The pair count is tested on
    # the array in test_tv_pairs_by_hand_with_the_factor_of_one_half_and_unique_pairs.


def test_distance_decay_is_spearman_not_pearson():
    # The same record with pair (1,2) moved from 6 to 100 km. The distance ranks are unchanged
    # ([1, 5, 2, 6, 3, 4]), so the Spearman value is still 12 / sqrt(12 * 17.5) = 0.8280787.
    # Pearson on TV = [0.5, 1, 0.5, 1, 0.5, 0.5] and d = [1, 5, 2, 100, 3, 4]: means 2/3 and 115/6;
    # sum xy = (1/6)(66.667) + (1/3)(66.667) = 33.333; sum xx = 1/3; sum yy = 7850.83;
    # r = 33.333 / sqrt(7850.83 / 3) = 0.6516.
    D = np.array([[0, 1, 5, 2], [1, 0, 100, 3], [5, 100, 0, 4], [2, 3, 4, 0]], float)
    val, reason = sr.distance_decay_rank(Y4, D)
    assert val == pytest.approx(0.8280787, abs=1e-6) and reason == ""
    assert pearsonr([0.5, 1, 0.5, 1, 0.5, 0.5], [1, 5, 2, 100, 3, 4])[0] == pytest.approx(0.6516, abs=1e-4)


@pytest.mark.parametrize("y", [np.array([[0, 0], [1, 2]]), np.array([[1, -1], [1, 2]]), np.array([[1.0, np.nan], [1, 2]]), np.array([1, 2, 3])])
def test_readouts_reject_invalid_counts(y):
    with pytest.raises(ValueError):
        sr.tv_pairs(y)
    with pytest.raises(ValueError):
        sr.lat_rank(y, np.array([1.0, 2.0]))


def test_readouts_reject_mismatched_shapes():
    y = np.array([[1, 2], [3, 4], [5, 6]])
    with pytest.raises(ValueError):
        sr.lat_rank(y, np.array([1.0, 2.0]))
    with pytest.raises(ValueError):
        sr.distance_decay_rank(y, np.zeros((2, 2)))
    with pytest.raises(ValueError):
        sr.distance_decay_rank(y, np.array([[0, 1, 2], [9, 0, 3], [2, 3, 0]], float))   # not symmetric
