"""Rim data plumbing (rim plan, 2026-10-19): matching to the basin set, published-count agreement, sums, attribute summaries."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from mls_emergence.dataio.basin import load_curated
from mls_emergence.rims import data as rd


def test_nineteen_basin_assemblages_have_rims_and_kent_place_is_kent():
    counts, _ = load_curated()
    matches, unmatched = rd.match_basin(list(counts.index))
    assert len(matches) == 19 and matches["Kent_Place"] == "Kent"
    assert sorted(unmatched) == ["Commerce", "Cramor_Place", "Cummins", "Grant", "Holden_Lake", "Hollywood", "Lake_Cormorant", "Nickel", "Turnbow"]
    with pytest.raises(ValueError):
        rd.match_basin(["Kent_Place", "Kent"])                      # two names onto one rim site is refused


def test_basin_rim_table_matches_the_published_counts_and_sums():
    counts, _ = load_curated()
    t, unmatched = rd.basin_rim_table(list(counts.index))
    assert len(t) == 19 and len(unmatched) == 9 and list(t.index) == [n for n in counts.index if n in set(t.index)]
    assert int(t.loc["Kent_Place", "N"]) == 153 and int(t.loc["Parkin", "N"]) == 67 and int(t.loc["Irby", "N"]) == 277
    assert t["N"].min() == 27 and t["N"].max() == 277                 # Woodlyn to Irby
    for cols in (rd.BEVEL_DEGREES, rd.LIP, rd.PASTE):
        assert (t[list(cols)].sum(axis=1) == t["N"]).all()
    assert (t[list(rd.PUBLISHED)].le(t["N"], axis=0)).all().all()
    assert int(t["N"].sum()) == sum(int(x) for x in t["N"])          # 19-site rim total is the sum of the published Ns


def test_attribute_summary_reports_rarity_without_deciding_anything():
    counts, _ = load_curated()
    t, _ = rd.basin_rim_table(list(counts.index))
    s = rd.attribute_summary(t).set_index("attribute")
    # The 19 matched sites hold 2,065 rims. Measured 2026-10-19 and pinned as data facts, not as a rule about what is usable.
    assert int(s.rims_examined.iloc[0]) == 2065
    assert {a: (int(s.loc[a, "rims_with"]), int(s.loc[a, "sites_with_any"])) for a in rd.PUBLISHED} == {
        "beaded": (28, 10), "noded": (13, 7), "wide_nodes": (10, 8), "notched_top": (23, 13), "notched_exterior": (160, 14), "beveled": (248, 18)}
    assert (int(s.loc["lip_flat", "rims_with"]), int(s.loc["paste_bell", "rims_with"]), int(s.loc["paste_other", "rims_with"])) == (708, 641, 550)
    assert not bool(s.loc["beveled", "tabulated_as_1951_type"]) and bool(s.loc["paste_bell", "tabulated_as_1951_type"])
    assert np.isclose(s["pooled_frequency"], s["rims_with"] / s["rims_examined"]).all()
