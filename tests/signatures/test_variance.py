import numpy as np
from mls_emergence.signatures import variance

def test_fst_zero_when_groups_identical():
    g = np.array([[10, 10, 10], [10, 10, 10], [10, 10, 10]])
    assert abs(variance.cultural_fst(g)) < 1e-9

def test_fst_high_when_groups_disjoint():
    g = np.array([[30, 0, 0], [0, 30, 0], [0, 0, 30]])
    assert variance.cultural_fst(g) > 0.6


def test_degenerate_input_raises_instead_of_returning_zero():
    """F4. The old code returned 0.0 for a degenerate pool, which is exactly the
    value 'no differentiation' takes, so an undefined case was indistinguishable
    from a real null and would have entered a trajectory as evidence for the
    paper's own conclusion.

    The probe that makes this discriminate: a pool with TWO classes and genuinely
    equal frequencies also gives F_ST = 0.0, and must keep doing so. If the guard
    were implemented by refusing everything that returns zero, this assertion
    fails."""
    import numpy as np
    import pytest
    from mls_emergence.signatures.variance import cultural_fst, gini_simpson

    # genuinely zero differentiation: identical groups, two classes present
    real_null = np.array([[10.0, 10.0], [20.0, 20.0]])
    assert cultural_fst(real_null) == pytest.approx(0.0, abs=1e-12)

    # degenerate: only one class in the whole pool -> undefined, not zero
    with pytest.raises(ValueError, match="single class"):
        cultural_fst(np.array([[5.0, 0.0], [3.0, 0.0]]))

    # fewer than two non-empty groups
    with pytest.raises(ValueError, match="at least two non-empty groups"):
        cultural_fst(np.array([[5.0, 2.0], [0.0, 0.0]]))

    # empty assemblage
    with pytest.raises(ValueError, match="empty assemblage"):
        gini_simpson(np.array([0.0, 0.0]))

    # negative counts are out of domain
    with pytest.raises(ValueError, match="non-negative"):
        cultural_fst(np.array([[-1.0, 2.0], [3.0, 4.0]]))


# --------------------------------------------------------------------------- #
# Value tests. Added 2026-09-04 after a mutation sweep showed the tests above
# could not tell cultural_fst from a wrong implementation: p**2 -> p**3,
# dividing by H_S instead of H_T, and dropping the size weighting all SURVIVED.
# The reason is that the two behavioural tests pin only the endpoints, and every
# one of those mutants still returns 0 for identical groups and something large
# for disjoint ones. Nothing pinned the arithmetic in between.
#
# Every expected value below is worked out by hand in exact fractions and
# written as the arithmetic, not as a decimal read back from the code (rule 3).
# --------------------------------------------------------------------------- #

def test_gini_simpson_exact_values():
    """1 - sum(p^2), pinned at two points. Kills any change to the exponent:
    with p**3 the second value would be 1 - (0.512 + 0.001 + 0.001) = 0.486."""
    import pytest
    assert variance.gini_simpson(np.array([1.0, 1.0])) == pytest.approx(0.5, abs=1e-12)
    #   p = (0.8, 0.1, 0.1); sum p^2 = 0.64 + 0.01 + 0.01 = 0.66
    assert variance.gini_simpson(np.array([8.0, 1.0, 1.0])) == pytest.approx(
        1.0 - 0.66, abs=1e-12)


def test_cultural_fst_exact_value_on_asymmetric_groups():
    """Pins the whole computation on a case with a known answer.

    groups (80,10,10) and (10,10,80), 100 sherds each:
        each group  p = (0.8, 0.1, 0.1) or its mirror, sum p^2 = 0.66, GS = 0.34
        H_S         = 0.34            (equal sizes)
        pooled      (90, 20, 90), p = (0.45, 0.10, 0.45), sum p^2 = 0.415
        H_T         = 0.585
        F_ST        = (0.585 - 0.34) / 0.585 = 0.245 / 0.585

    Dividing by H_S instead of H_T would give 0.245 / 0.34 = 0.72, which this
    assertion rejects.
    """
    import pytest
    g = np.array([[80, 10, 10], [10, 10, 80]])
    assert variance.cultural_fst(g) == pytest.approx(0.245 / 0.585, abs=1e-12)


def test_within_group_diversity_is_weighted_by_group_size():
    """The discriminating probe for the size weighting.

    Groups of very different size AND different diversity:
        (50, 50): 100 sherds, GS = 0.5
        (9, 1):    10 sherds, GS = 0.18
        size-weighted H_S = (0.5*100 + 0.18*10) / 110 = 51.8 / 110
        unweighted   H_S = (0.5 + 0.18) / 2          = 0.34
        pooled (59, 51): H_T = 1 - (59^2 + 51^2) / 110^2

    The two answers are 0.053 and 0.316, a factor of six apart, so an unweighted
    mean cannot hide here. The endpoint tests above cannot see this difference
    at all, which is why this case is the one that was missing.
    """
    import pytest
    g = np.array([[50, 50], [9, 1]])
    H_T = 1.0 - (59.0 ** 2 + 51.0 ** 2) / 110.0 ** 2
    H_S = (0.5 * 100.0 + 0.18 * 10.0) / 110.0
    assert variance.cultural_fst(g) == pytest.approx((H_T - H_S) / H_T, abs=1e-12)
    # and it is far from what the unweighted mean would give
    assert abs(variance.cultural_fst(g) - (H_T - 0.34) / H_T) > 0.25


def test_fst_rises_monotonically_as_groups_separate():
    """Ordering rather than endpoints: mixing two groups from identical toward
    disjoint must increase F_ST at every step. A monotonicity claim needs more
    than the two ends, since the mutants reproduce both ends correctly."""
    prev = -1.0
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        a = 50.0 + 50.0 * frac
        g = np.array([[a, 100.0 - a], [100.0 - a, a]])
        v = variance.cultural_fst(g)
        assert v > prev, f"F_ST not increasing at separation {frac}"
        prev = v
