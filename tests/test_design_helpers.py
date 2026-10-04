"""Helpers of the design analysis (97) and the collections inventory (98).

Rule 3: each assertion includes a case where the nearest wrong implementation
gives a different answer.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "analyses", ROOT / "scripts", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

a97 = importlib.import_module("97_design_to_separate")
a98 = importlib.import_module("98_small_collections")


def test_split_profile_adds_back_to_the_original_classes():
    pooled = np.array([0.6, 0.3, 0.1])
    fine = a97.split_profile(pooled, 4)
    assert fine.size == 12
    # class by class, not just in total: an even split of the grand total would
    # also sum to one and would fail here
    assert np.allclose(fine.reshape(3, 4).sum(1), pooled)
    # shares within a class fall by halves, so each split has rarer members
    assert np.allclose(fine[:4] / fine[0], [1, 0.5, 0.25, 0.125])
    assert np.array_equal(a97.split_profile(pooled, 1), pooled)
    with pytest.raises(ValueError):
        a97.split_profile(pooled, 0)


def test_base_number_strips_collection_units_and_rejects_other_ids():
    assert a98.base_number("13-N-4/B,C,D,E") == "13-N-4"
    assert a98.base_number("13-N-4/C") == "13-N-4"
    assert a98.base_number(" 12-N-3/A&B") == "12-N-3"
    # 13-N-4 and 13-N-14 are different sites; a prefix match would merge them
    assert a98.base_number("13-N-14") != a98.base_number("13-N-4")
    assert a98.base_number("40LA7") is None
    assert a98.base_number("22Pa528") is None


def test_every_analyzed_assemblage_has_a_recorded_site_number():
    mf = importlib.import_module("make_figures")
    counts, _ = mf._load_curated()
    names = [str(i) for i in counts.index]
    assert set(names) <= set(a98.ANALYZED_NUMBER)
    numbers = [v for k, v in a98.ANALYZED_NUMBER.items() if v is not None]
    assert len(numbers) == len(set(numbers))   # two assemblages on one number would hide a candidate
