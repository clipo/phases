"""The painted class by type: row selection and totals of analyses/per_type/painted_by_type.py."""
import importlib
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "per_type", ROOT / "scripts", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
pbt = importlib.import_module("painted_by_type")
TYPES = ["Old Town Red", "Nodena Red and White", "Avenue Polychrome"]


@pytest.fixture(scope="module")
def table():
    return pbt.by_type()[0]


@pytest.fixture(scope="module")
def sheet():
    """PFGData read here directly, not through the script, so expected values do not pass through the code under test."""
    d = pd.read_excel(ROOT / "data" / "raw" / "PFGData.xlsx", sheet_name="SherdData")
    d.columns = [str(c).strip() for c in d.columns]
    d["Site Number"] = d["Site Number"].astype(str).str.strip()
    return d


def test_kent_place_takes_the_ruled_row_and_does_not_count_collection_c_twice(table, sheet):
    """PFGData lists Kent Place as 13-N-4/B,C,D,E and again as 13-N-4/C. The ruled 1940 part is the first row alone."""
    first = sheet[sheet["Site Number"] == "13-N-4/B,C,D,E"][TYPES].apply(pd.to_numeric, errors="coerce").fillna(0).sum()
    both = sheet[sheet["Site Number"].str.startswith("13-N-4")][TYPES].apply(pd.to_numeric, errors="coerce").fillna(0).sum()
    assert int(both.sum()) > int(first.sum())            # the rival reading (summing both rows) gives a larger count
    assert [int(table.loc["Kent_Place", t]) for t in TYPES] == [int(first[t]) for t in TYPES]


def test_sites_keyed_by_grid_number_are_found(table, sheet):
    """Nickel is spelled 'Nickle' and split over two rows in PFGData; a match on the site name alone finds nothing."""
    assert not sheet["Site Name"].astype(str).str.strip().str.lower().eq("nickel").any()
    rows = sheet[sheet["Site Number"].str.split("/").str[0] == "13-N-15"]
    assert len(rows) == 2
    assert int(table.loc["Nickel", "Old Town Red"]) == int(pd.to_numeric(rows["Old Town Red"], errors="coerce").fillna(0).sum()) > 0
    assert table.loc["Holden_Lake", TYPES].sum() == 0 and table.loc["Holden_Lake", "decorated, 1940"] == 0   # Lipo's alone


def test_types_and_additions_account_for_every_painted_sherd_of_the_matrix(table):
    assert ((table[TYPES].sum(axis=1) + table["painted, Lipo addition"]) == table["painted, matrix"]).all()
    s = pbt.summary(table)
    assert s["n_sites"] == 27 and s["painted_1940"] + s["addition_painted"] == s["matrix_painted"]
    assert s["sites_named"]["Nodena Red and White"] == ["Hollywood"]
    # The sentence in the paper: nine-tenths of the 1940 painted sherds are Old Town Red. Below one half would mean
    # the pooled class is mainly something else and its latitude trend could not be read as Old Town Red's.
    assert 0.89 < s["share"]["Old Town Red"] < 0.91


def test_the_two_further_painted_types_are_counted_at_hollywood_and_left_out_at_williamson(table, sheet):
    """Lipo (2001: Table 7.4, note a) pools Carson Red-and-Buff and Hollywood White into Painted. Expected counts are
    read from the sheet here, not through the script."""
    def extra(site_number):
        r = sheet[sheet["Site Number"] == site_number][["Carson Red-on-Buff", "Hollywood White Slipped"]]
        return int(r.apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy().sum())
    s = pbt.summary(table)
    assert s["extra_counted"] == {"Hollywood": extra("13-O-10")} and extra("13-O-10") > 0
    assert s["extra_left_out"] == {"Williamson": extra("11-N-13")} and extra("11-N-13") > 0
    # If the matrix pooled three types everywhere, Hollywood's painted count would equal its three-type 1940 count
    # and extra_counted would be empty; if it pooled five everywhere, Williamson would not be left out.
    assert int(table.loc["Hollywood", "painted, matrix"]) == int(table.loc["Hollywood", TYPES].sum()) + extra("13-O-10")
    assert int(table.loc["Williamson", "painted, matrix"]) == int(table.loc["Williamson", TYPES].sum())
    assert s["addition_painted_own"] + sum(s["extra_counted"].values()) == s["addition_painted"]
