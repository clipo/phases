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


def test_phase_profile_difference_is_the_share_of_sherds_to_move():
    a99 = importlib.import_module("99_phase_profiles")
    # two phases of two assemblages each; phase 0 is 80/20, phase 1 is 50/50
    m = np.array([[40, 10], [40, 10], [10, 10], [40, 40]])
    phase = np.array([0, 0, 1, 1])
    prof = a99.phase_profiles(m, phase)
    assert np.allclose(prof, [[80, 20], [50, 50]])     # pooled by phase, not a mean of assemblage percentages
    # moving 30 percent of phase 0's sherds from class 0 to class 1 makes the profiles alike;
    # the summed absolute difference (60) would double-count each moved sherd
    assert a99.mean_pair_difference(prof) == pytest.approx(30.0)
    assert a99.pair_differences(prof) == {(0, 1): pytest.approx(30.0)}
    assert a99.mean_pair_difference(np.array([[100.0, 0.0], [0.0, 100.0]])) == pytest.approx(100.0)


def test_rim_count_table_matches_the_published_counts():
    """The count table must carry Mainfort's (2003b) printed counts, not a rounding of percentages."""
    import pandas as pd
    built = pd.read_csv(ROOT / "data" / "processed" / "mainfort_rim_counts.csv").set_index("site")
    pub = pd.read_csv(ROOT / "data" / "raw" / "mainfort2003_rim_table1_counts.csv", comment="#").set_index("site")
    assert len(built) == len(pub) == 55 and int(pub["N"].sum()) == 7195       # the paper: 7,195 rims, 55 sites
    assert (built.loc[pub.index, "N"] == pub["N"]).all()
    # Upper Nodena as printed: 26 beaded, 27 notched exterior, 29 beveled of 103 (Mainfort et al. 2007: Table 8
    # gives the same sample as 25.2, 26.2 and 28.2 percent)
    un = built.loc["Upper Nodena"]
    assert (un.beaded, un.notched_exterior, un.beveled, un.N) == (26, 27, 29, 103)
    assert un.bevel_sharp + un.bevel_moderate + un.bevel_slight == un.beveled
    assert un.bevel_none == un.N - un.beveled and un.lip_flat + un.lip_round == un.N
    # Author's ruling: the printed counts are the record, even the one that contradicts the paper's own
    # percentage table (Pecan Point's beaded count, 10 of 30, where Table 2 gives 10 percent).
    six = ["beaded", "noded", "wide_nodes", "notched_top", "notched_exterior", "beveled"]
    assert (built.loc[pub.index, six].to_numpy() == pub[six].to_numpy()).all()
    assert built.loc["Pecan Point", "beaded"] == 10


def test_largest_remainder_keeps_the_total():
    sys.path.insert(0, str(ROOT / "scripts"))
    b = importlib.import_module("build_mainfort_rim_counts")
    out = b.largest_remainder(np.array([8.0, 15.0, 6.0]), 29)
    assert out.sum() == 29 and list(out) == [8, 15, 6]
    # plain rounding of 1/3 each of 10 gives 3 + 3 + 3 = 9 and loses a rim
    assert b.largest_remainder(np.array([1.0, 1.0, 1.0]), 10).sum() == 10
    assert list(b.largest_remainder(np.array([0.0, 0.0, 0.0]), 0)) == [0, 0, 0]
