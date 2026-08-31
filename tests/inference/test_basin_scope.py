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
    # The membership file is 29 assemblages; all of them carry coordinates.
    assert len(fitted) == 29


@pytest.mark.data
def test_region_scope_is_genuinely_different(inp, a43):
    """The probe. If basin and region agreed, F18 would have been harmless."""
    gc_b, _ = a43.basin_group_counts(inp, scope="basin")
    gc_r, _ = a43.basin_group_counts(inp, scope="region")

    # Cluster counts differ: k is re-selected on the basin's own coordinates.
    assert gc_b.shape[0] == 3, "basin should give three spatial clusters"
    assert gc_r.shape[0] == 5, "region should give five"
    # And the region pools strictly more sherds, because it pools more sites.
    assert gc_r.sum() > gc_b.sum()
    # Same type vocabulary either way; only the grouping and scope change.
    assert gc_b.shape[1] == gc_r.shape[1]


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
