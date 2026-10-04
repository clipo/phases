"""The Phillips (1970) phase outlines must sit where the source figure puts them.

`data/raw/phillips1970_phase_outlines.csv` is read out of the vector artwork of
Lipo (2001: Figure 2.6) and georeferenced from the figure's quadrangle frame
(`scripts/build_phillips1970_outlines.py`). Manuscript Figure 1 draws it. The
failure that would matter is a wrong georeference: a flipped axis, or the frame
taken one quadrangle off, still yields four tidy outlines inside a plausible
box.

Rule 3: the probes are sites whose phase area the source figure shows
unambiguously, with coordinates written here from the project's assemblage
table rather than computed through the builder. Several sit near an outline's
edge on purpose (Barton Ranch 5.7 km inside Parkin's line, Clay Hill 1.2 km
inside Kent's, both as written to output/findings/phillips1970_outline_membership.md),
so a shift of a quarter-degree quadrangle, or a mirrored axis, moves them out:
each of those six mutations fails at least seven of the nine probes. The
builder names each outline by one type site; the probes below are different
sites except where noted.
"""
from pathlib import Path

import pandas as pd
import pytest
from shapely.geometry import Point, Polygon

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "data" / "raw" / "phillips1970_phase_outlines.csv"

# (longitude, latitude), from data/raw/mainfort-pfg-cplXY.txt as corrected.
INSIDE = {
    "Nodena": {"Carson_Lake": (-90.0470, 35.5970), "Notgrass": (-90.1035, 35.5384)},
    "Parkin": {"Barton_Ranch": (-90.3895, 35.3791), "Castile_Landing": (-90.6778, 34.9777),
               "Parkin": (-90.5560, 35.2760)},   # Parkin is also the builder's type site
    "Walls": {"Irby": (-90.1564, 34.9179), "Mound_Place": (-90.2390, 35.0790)},
    "Kent": {"Clay_Hill": (-90.7242, 34.8770), "Starkley": (-90.5844, 34.7706)},
}
# In the frame but in no outline: a point 15.5 km north of the Parkin area. It
# is the coordinate the assemblage file carried for Cummins until 2026-10-01,
# which belongs to a different site; Cummins itself (11-O-4) is inside the
# Parkin area. (Holden Lake and Cramor Place are outside by 1.0 km, which is
# inside the source's drafting error, so they are not used as probes.)
OUTSIDE_ALL = {"north of the Parkin area": (-90.6333, 35.6193)}


def _outlines():
    df = pd.read_csv(CSV, comment="#")
    return df, {ph: Polygon(g.sort_values("vertex")[["longitude", "latitude"]].to_numpy())
                for ph, g in df.groupby("phase")}


@pytest.mark.data
def test_four_simple_outlines_inside_the_frame():
    df, poly = _outlines()
    assert sorted(poly) == ["Kent", "Nodena", "Parkin", "Walls"]
    assert all(p.is_valid for p in poly.values()), "an outline crosses itself"
    # The source figure's frame, LMS columns N-Q by rows 10-14.
    assert df["longitude"].between(-90.75, -89.75).all()
    assert df["latitude"].between(34.50, 35.75).all()
    # Pinned so a re-extraction that picks up different paths is noticed.
    assert df.groupby("phase").size().to_dict() == {
        "Kent": 70, "Nodena": 265, "Parkin": 86, "Walls": 54}


@pytest.mark.data
def test_each_probe_site_is_in_its_own_outline_and_no_other():
    _, poly = _outlines()
    for phase, sites in INSIDE.items():
        for name, (lon, lat) in sites.items():
            hit = sorted(ph for ph, p in poly.items() if p.contains(Point(lon, lat)))
            assert hit == [phase], f"{name} falls in {hit or 'no outline'}, expected {phase}"


@pytest.mark.data
def test_sites_the_source_leaves_outside_are_outside():
    # One-sided on purpose: a site inside an outline here would mean the
    # outlines had grown or shifted, not that the site had moved.
    _, poly = _outlines()
    for name, (lon, lat) in OUTSIDE_ALL.items():
        hit = [ph for ph, p in poly.items() if p.contains(Point(lon, lat))]
        assert not hit, f"{name} falls inside {hit}"


@pytest.mark.data
def test_outlines_do_not_overlap():
    # The four areas are separate in the source; Figure 1 fills them one gray
    # on that basis, and the membership table counts no site in two.
    _, poly = _outlines()
    names = sorted(poly)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            assert poly[a].intersection(poly[b]).area == 0, f"{a} and {b} overlap"


def _m36():
    import importlib
    import sys
    for p in (ROOT, ROOT / "analyses", ROOT / "src"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    return importlib.import_module("36_canonical_phase_map")


@pytest.mark.data
def test_phillips_assignment_differs_from_the_1996_one_where_the_sources_say():
    # Rule 3: the rival account is the 1996 assignment the pipeline uses, which
    # puts Beck and Hollywood in Walls and Castile Landing in Kent. Mainfort
    # (2003:176-177) names Belle Meade, Beck, Hollywood and Commerce as
    # Phillips's Kent phase sites and Mound Place as his Walls; his Table 1,
    # which follows Phillips there, has Castile Landing in Parkin.
    m36 = _m36()
    sites = {  # (latitude, longitude)
        "Beck": (34.9125, -90.3939), "Hollywood": (34.7774, -90.3707),
        "Castile_Landing": (34.9777, -90.6778), "Mound_Place": (35.0790, -90.2390),
        "north of Parkin": (35.6193, -90.6333),   # in the frame, in no outline
        "Cummins": (35.458177, -90.477238),       # 11-O-4, as corrected 2026-10-01
        "Parchman": (34.3586, -90.5657),     # south of the frame
    }
    names = list(sites)
    labels, derived = m36.assign_phases_phillips(names, list(sites.values()))
    got = dict(zip(names, zip(labels, derived)))
    assert got["Beck"] == ("Kent", False)
    assert got["Hollywood"] == ("Kent", False)
    assert got["Castile_Landing"] == ("Parkin", False)
    assert got["Mound_Place"] == ("Walls", False)
    assert got["north of Parkin"] == ("Parkin", True), "nearest outline, flagged as derived"
    assert got["Cummins"] == ("Parkin", False)
    assert got["Parchman"] == ("unassigned", True), "outside the frame is not assigned"
    old = dict(zip(names, m36.assign_phases(names)))
    assert (old["Beck"], old["Hollywood"], old["Castile_Landing"]) == ("Walls", "Walls", "Kent")
