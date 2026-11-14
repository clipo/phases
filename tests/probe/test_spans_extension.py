"""Timing functions of the drift-extension design (spec 2026-11-14, draft 3, sections 3.3, 5 and 8).

Each of the three functions is a deterministic function of its arguments, with
no random generator inside, so every expected value here is worked by hand in
the comment beside it and none is computed through the function under test.
Rule 3: each test names the wrong implementation it rejects.
"""
import numpy as np
import pytest
from scipy.stats import spearmanr

from mls_emergence.probe import spans

# Five sites. Latitude ranks (1 = southernmost): [3, 1, 5, 2, 4].
LAT5 = np.array([30.0, 10.0, 50.0, 20.0, 40.0])
U5 = np.array([0.9, 0.8, 0.1, 0.2, 0.1])
ENDS5 = np.array([3, 7, 1, 7, 3])          # tied end times; sorted: [1, 3, 3, 7, 7]


def test_lengths_from_uniforms_exact():
    # 1 + floor(v L) at L = 46: floor(0) = 0, floor(1.38) = 1, floor(23.0) = 23, floor(45.9954) = 45.
    v = np.array([0.0, 0.03, 0.5, 0.9999])
    got = spans.lengths_from_uniforms(v, 46)
    assert got.tolist() == [1, 2, 24, 46]
    # Rivals rejected: all-one lengths; endpoint-only lengths (1 or 46), which have the uniform mean
    # 23.5 and so pass a mean check; a 0-based mapping floor(v L).
    assert np.issubdtype(got.dtype, np.integer)


def test_lengths_limit_one_is_all_ones():
    assert spans.lengths_from_uniforms(np.array([0.0, 0.5, 0.9999]), 1).tolist() == [1, 1, 1]


@pytest.mark.parametrize("v, L", [(np.array([1.0]), 5), (np.array([-0.1]), 5), (np.array([0.5]), 0),
                                  (np.array([0.5]), 2.5), (np.array([np.nan]), 5), (np.array([0.5]), True)])
def test_lengths_guards(v, L):
    with pytest.raises(ValueError):
        spans.lengths_from_uniforms(v, L)


def test_ends_coupled_interior_c_north_late_by_hand():
    # q = (r - 0.5) / 5 = [0.5, 0.1, 0.9, 0.3, 0.7]; z = 0.5 q + 0.5 u = [0.70, 0.45, 0.50, 0.25, 0.40].
    # Ascending z: site 3, site 4, site 1, site 2, site 0. Sorted ends [1, 3, 3, 7, 7] go to them in
    # that order: site3 <- 1, site4 <- 3, site1 <- 3, site2 <- 7, site0 <- 7.
    got = spans.ends_coupled(ENDS5, LAT5, 0.5, True, U5)
    assert got.tolist() == [7, 3, 7, 1, 3]
    # Rivals rejected: ignoring the interior c and ordering by latitude alone gives [3, 1, 7, 3, 7];
    # ordering by u alone gives [7, 7, 1, 3, 3].


def test_ends_coupled_interior_c_south_late_by_hand():
    # q = 1 - (r - 0.5) / 5 = [0.5, 0.9, 0.1, 0.7, 0.3]; z = [0.70, 0.85, 0.10, 0.45, 0.20].
    # Ascending z: site 2, site 4, site 3, site 0, site 1 -> site2 <- 1, site4 <- 3, site3 <- 3,
    # site0 <- 7, site1 <- 7.
    got = spans.ends_coupled(ENDS5, LAT5, 0.5, False, U5)
    assert got.tolist() == [7, 7, 1, 3, 3]
    # Rival rejected: a south-late design that applies the north-late assignment ([7, 3, 7, 1, 3]).


# Not testable, and recorded so that it is not mistaken for a gap: the offset in q = (r - 0.5) / n cannot be
# distinguished from q = r / n by any assignment. The two differ by the constant c / (2 n) in every site's score,
# which leaves the order of the scores unchanged for every c, u and orientation. A mutation pass on 2026-11-14
# confirmed that replacing one with the other changes no output. The offset matters only to the stated expected
# correlation of latitude with the score, which is also invariant to it.


def test_ends_coupled_full_coupling_orders_ends_by_latitude_with_ties():
    north = spans.ends_coupled(ENDS5, LAT5, 1.0, True, U5)
    south = spans.ends_coupled(ENDS5, LAT5, 1.0, False, U5)
    by_lat = np.argsort(LAT5)                      # south to north: sites 1, 3, 0, 4, 2
    assert north[by_lat].tolist() == [1, 3, 3, 7, 7]   # nondecreasing northward
    assert south[by_lat].tolist() == [7, 7, 3, 3, 1]   # nonincreasing northward
    # The opposite order would mean the later ends go to the other end of the basin.


def test_ends_coupled_zero_coupling_uses_u_and_breaks_ties_by_site_order():
    # At c = 0 the score is u exactly. u ties at sites 2 and 4 (0.1): the stable sort keeps site 2
    # first. Ascending: site 2, site 4, site 3, site 1, site 0.
    got = spans.ends_coupled(ENDS5, LAT5, 0.0, True, U5)
    assert got.tolist() == [7, 7, 1, 3, 3]


@pytest.mark.parametrize("c", [0.0, 0.25, 0.5, 1.0])
@pytest.mark.parametrize("north", [True, False])
def test_ends_coupled_preserves_the_multiset_of_ends(c, north):
    # Rival rejected: a function that redraws or interpolates end times.
    got = spans.ends_coupled(ENDS5, LAT5, c, north, U5)
    assert sorted(got.tolist()) == sorted(ENDS5.tolist())


def test_ends_coupled_signed_expected_correlation_both_orientations():
    # At c = 0.5 and 28 sites the expected rank correlation of latitude with the score is
    # +549/784 = +0.7003 north-late and -0.7003 south-late (repair check RC-F3). The score is not
    # returned, so the test uses distinct end times 0..27, for which the assigned end's rank equals
    # the score's rank. Positive means the later ends lie to the north.
    n = 28
    lat = np.linspace(34.0, 36.0, n)[np.random.default_rng(5).permutation(n)]
    ends = np.arange(n)
    rng = np.random.default_rng(20261114)
    north, south = [], []
    for _ in range(2000):
        u = rng.random(n)
        north.append(spearmanr(lat, spans.ends_coupled(ends, lat, 0.5, True, u))[0])
        south.append(spearmanr(lat, spans.ends_coupled(ends, lat, 0.5, False, u))[0])
    assert abs(np.mean(north) - 0.700) < 0.01
    assert abs(np.mean(south) + 0.700) < 0.01


@pytest.mark.parametrize("kw", [{"c": -0.1}, {"c": 1.1}, {"c": np.nan}, {"north_late": 1}, {"latitude": np.array([1.0, 1.0, 2.0, 3.0, 4.0])},
                                {"u": np.array([0.1, 0.2, 0.3, 0.4, 1.0])}, {"u": np.array([0.1, 0.2])}, {"ends": np.array([1, 2])}])
def test_ends_coupled_guards(kw):
    args = {"ends": ENDS5, "latitude": LAT5, "c": 0.5, "north_late": True, "u": U5}; args.update(kw)
    with pytest.raises(ValueError):
        spans.ends_coupled(**args)


LABELS6 = np.array([0, 0, 1, 1, 0, 1])
W6 = np.array([0.0, 0.99, 0.0, 0.5, 0.34, 0.99])


def test_ends_step_exact_odd_window_late_label_one():
    # S = 20, W = 7: window steps 13..19. Early half: floor(7/2) = 3 steps, 13..15. Late half: 4
    # steps, 16..19. End = h0 + floor(w h).
    # site0 early 13 + floor(0) = 13; site1 early 13 + floor(2.97) = 15; site2 late 16 + 0 = 16;
    # site3 late 16 + floor(2.0) = 18; site4 early 13 + floor(1.02) = 14; site5 late 16 + floor(3.96) = 19.
    got = spans.ends_step(W6, LABELS6, 1, 20, 7)
    assert got.tolist() == [13, 15, 16, 18, 14, 19]
    # Rivals rejected: half-window midpoints (14 and 17 or 18 everywhere); endpoint-only draws,
    # which cannot return 14 or 18; a split at ceil(W/2).


def test_ends_step_exact_odd_window_late_label_zero():
    # Roles swapped: site0 late 16; site1 late 16 + floor(3.96) = 19; site2 early 13; site3 early
    # 13 + floor(1.5) = 14; site4 late 16 + floor(1.36) = 17; site5 early 13 + floor(2.97) = 15.
    got = spans.ends_step(W6, LABELS6, 0, 20, 7)
    assert got.tolist() == [16, 19, 13, 14, 17, 15]


def test_ends_step_separates_the_clusters_and_stays_in_the_window():
    rng = np.random.default_rng(3)
    labels = np.array([0] * 17 + [1] * 11)
    for W in (80, 81, 2336):
        e = spans.ends_step(rng.random(28), labels, 0, 2600, W)
        assert e[labels == 0].min() > e[labels == 1].max()   # the late cluster is strictly later
        assert e.min() >= 2600 - W and e.max() <= 2599
        assert e[labels == 1].max() <= 2600 - W + W // 2 - 1  # the split is at floor(W/2)


@pytest.mark.parametrize("kw", [{"labels": np.array([0, 0, 0, 0, 0, 0])}, {"labels": np.array([0, 1, 2, 0, 1, 2])}, {"late_label": 5},
                                {"W": 1}, {"W": 21}, {"W": 7.0}, {"w": np.array([0.0, 0.5, 0.5, 0.5, 0.5, 1.0])}, {"w": np.array([0.1, 0.2])}])
def test_ends_step_guards(kw):
    args = {"w": W6, "labels": LABELS6, "late_label": 1, "S": 20, "W": 7}; args.update(kw)
    with pytest.raises(ValueError):
        spans.ends_step(**args)
