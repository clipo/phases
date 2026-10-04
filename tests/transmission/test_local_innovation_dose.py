"""How many areas a given "strength" of local innovation actually changes.

`65_other_departures.group_targets` gives each area a reordered copy of the
regional profile. A class is drawn with probability `strength`, and the drawn
classes are permuted among themselves, so an area's profile changes only when
two or more classes are drawn AND the permutation is not the identity. The
share of areas changed is therefore well below the strength's face value, and
analysis 92's site-by-site rows have to be read at the dose they deliver.

Rule 3: the expected value is worked from the binomial and the count of
non-identity permutations, not through the function under test, and the two
nearest wrong readings predict clearly different shares at strength 0.1 with
ten classes: "any drawn class changes the profile" gives 0.65, and "two or more
drawn classes change it" gives 0.26. The function gives 0.16.
"""
import importlib
import sys
from math import comb, factorial
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

a65 = importlib.import_module("65_other_departures")

# Ten classes with distinct frequencies, skewed like the record's.
POOLED = np.array([0.52, 0.38, 0.03, 0.02, 0.015, 0.012, 0.009, 0.007, 0.005, 0.002])


def expected_changed(strength: float, k: int) -> float:
    """P(an area's profile differs): m classes drawn, and not the identity among them."""
    return sum(comb(k, m) * strength ** m * (1 - strength) ** (k - m) * (1 - 1 / factorial(m))
               for m in range(2, k + 1))


@pytest.mark.parametrize("strength", [0.1, 0.2, 0.25])
def test_share_of_areas_changed_matches_the_binomial(strength):
    n = 6000
    tgt = a65.group_targets(POOLED, np.arange(n), strength, np.random.default_rng(11))
    got = float((np.abs(tgt - POOLED).sum(1) > 0).mean())
    want = expected_changed(strength, POOLED.size)
    # three binomial standard errors
    assert abs(got - want) < 3 * np.sqrt(want * (1 - want) / n)


def test_strength_point_one_changes_about_one_area_in_six():
    want = expected_changed(0.1, 10)
    assert want == pytest.approx(0.157, abs=0.001)
    # The rival readings: any drawn class (0.65), or two or more drawn (0.26).
    assert 1 - 0.9 ** 10 == pytest.approx(0.651, abs=0.001)
    assert 1 - 0.9 ** 10 - 10 * 0.1 * 0.9 ** 9 == pytest.approx(0.264, abs=0.001)
    tgt = a65.group_targets(POOLED, np.arange(6000), 0.1, np.random.default_rng(3))
    got = float((np.abs(tgt - POOLED).sum(1) > 0).mean())
    assert 0.13 < got < 0.19   # far from both rivals


def test_reordering_keeps_the_frequencies_and_zero_strength_is_the_regional_profile():
    tgt = a65.group_targets(POOLED, np.arange(500), 0.3, np.random.default_rng(5))
    assert np.allclose(np.sort(tgt, axis=1), np.sort(POOLED))
    assert np.array_equal(a65.group_targets(POOLED, np.arange(5), 0.0, np.random.default_rng(5)), POOLED)
