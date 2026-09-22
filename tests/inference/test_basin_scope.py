"""The Bayesian F_ST fit must run on the basin, and prove it.

Regression test for docs/CODE_REVIEW_2026-08-31.md F18: until 2026-08-31
`43_bayesian_fst.basin_group_counts` ignored the canonical drainage-basin
membership and silently fitted the whole 55-assemblage curated set while being
named for the basin.

Rule 3 (tests must discriminate): the assertions below are paired with the
`scope="region"` probe, where the rival account -- "it does not matter, the sets
are effectively the same" -- predicts no difference and is wrong. Expected values
come from `data/processed/basin_members_curated.txt`, not from the code path
under test.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for p in (str(ROOT), str(ROOT / "analyses")):
    if p not in sys.path:
        sys.path.insert(0, p)


@pytest.fixture(scope="module")
def inp():
    return importlib.import_module("07_refined_empirical").prepare_inputs()


@pytest.fixture(scope="module")
def a43():
    return importlib.import_module("43_bayesian_fst")


def _membership_from_file() -> set:
    """Read the expected ids from the committed artifact, not from the code."""
    txt = ROOT / "data" / "processed" / "basin_members_curated.txt"
    return set(txt.read_text().split())


@pytest.mark.data
def test_basin_fit_uses_exactly_the_canonical_membership(inp, a43):
    fitted = set(a43.fitted_basin_ids(inp))
    expected = _membership_from_file() & {str(a) for a in inp.have_coords_ids}
    assert fitted == expected, (
        f"fitted set differs from the membership file: "
        f"{sorted(fitted ^ expected)}")
    # The membership file is 43 assemblages; all of them carry coordinates.
    # It went 29 -> 30 on 2026-09-19 (Beck's coordinate corrected) and 30 -> 43
    # on 2026-09-20, when membership became phase membership rather than a
    # drainage corridor. The count is pinned deliberately: an accidental change
    # of scope should fail here rather than quietly alter every basin number.
    # 29 -> 30 (Beck's coordinate, 2026-09-19) -> 43 (phase membership,
    # 2026-09-20) -> 38 (2026-09-21, matrix rebuilt from PFGData and Lipo's
    # compilation; five basin assemblages known only from Mainfort's table
    # left) -> 28 (2026-09-21, author rulings: Parchman out on geography, a
    # 75-decorated-sherd minimum, and Nodena's one survivor not carried as a
    # phase of one).
    assert len(fitted) == 28


@pytest.mark.data
def test_region_scope_is_the_matrix_and_basin_scope_applies_the_rules(inp, a43):
    """Two probes. The matrix must not contain the assemblages known only
    from Mainfort's table (a loader falling back to the raw workbook would put
    them back), and the basin scope must exclude, by name, the ten assemblages
    the membership rules remove while the region scope keeps them. The sherd
    totals are read off the processed matrix, not off either scope.
    """
    gone = {"40LA007", "40TP026", "Bishop", "Fullen", "Graves_Lake", "Hatchie",
            "Jeter", "Jones_Bayou", "Porter", "Rast", "Richardsons_Landing",
            "Wilder", "Chuccalissa", "Soudan", "Wall", "West_Mounds", "Young"}
    present = {str(a) for a in inp.counts.index}
    assert not (gone & present), (
        f"assemblages known only from Mainfort's table are back in the matrix: "
        f"{sorted(gone & present)}; is a reader using the raw workbook?")
    # Excluded by the membership rules: Parchman's three on geography, six
    # under 75 decorated sherds, and Nodena's lone survivor. They stay in the
    # matrix file and must be absent from the INPUTS, since 07 now applies the
    # membership before fitting the correspondence axis (2026-09-22).
    excluded = {"Dundee", "Parchman", "Salomon", "Cheatham", "Connor", "Norfolk",
                "Notgrass", "Pouncey", "Upper_Nodena", "Carson_Lake"}
    import pandas as pd
    mat = pd.read_csv(ROOT / "data" / "processed" / "analysis_matrix.csv")
    assert excluded <= set(mat["Assemblages"].astype(str)), "excluded rows should stay in the matrix file"
    assert not (excluded & present), sorted(excluded & present)
    fitted = set(a43.fitted_basin_ids(inp))
    assert not (excluded & fitted), sorted(excluded & fitted)
    gc_b, _ = a43.basin_group_counts(inp, scope="basin")
    gc_r, _ = a43.basin_group_counts(inp, scope="region")
    # With the membership applied at load, the two scopes coincide, and the
    # total is pinned against the matrix file restricted by the members list.
    members = set((ROOT / "data" / "processed" / "basin_members_curated.txt").read_text().split())
    pinned = int(mat[mat["Assemblages"].isin(members)].iloc[:, 1:].to_numpy().sum())
    assert int(gc_b.sum()) == int(gc_r.sum()) == pinned == 14101


@pytest.mark.data
def test_default_scope_is_the_basin(inp, a43):
    """A default that every caller must override is not a default (rule 4)."""
    gc_default, _ = a43.basin_group_counts(inp)
    gc_basin, _ = a43.basin_group_counts(inp, scope="basin")
    assert np.array_equal(gc_default, gc_basin)


def test_unknown_scope_raises(inp, a43):
    """Out-of-domain input raises with a reason; it never returns plausible
    garbage (rule 5)."""
    with pytest.raises(ValueError, match="scope must be"):
        a43.basin_group_counts(inp, scope="whole-lmv")
