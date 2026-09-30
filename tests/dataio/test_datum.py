"""The NAD27 conversion moves only what it can demonstrate, by the right amount.

Expected positions are computed with pyproj directly from PFG's published UTM,
not through the loader under test. The discriminating probes: Parkin, whose
file coordinate is 41 m from PFG's UTM read as NAD83, must NOT move (a rule of
"convert every PFG site" would move it); Grant has no PFG UTM and must not
move; a settlement row from another compiler must not move.
"""
from pathlib import Path

import numpy as np
import pytest

from mls_emergence.dataio import coords as C
from mls_emergence.dataio import settlement as S

XY = Path("data/raw/mainfort-pfg-cplXY.txt")
LMV = Path("data/LMVData.xlsx")
SPATIAL = Path("data/raw/mainfort-spatial-data.xlsx")
needs_xy = pytest.mark.skipif(not (XY.exists() and SPATIAL.exists()), reason="raw data absent")
needs_lmv = pytest.mark.skipif(not LMV.exists(), reason="LMVData gitignored/absent")

# PFG site table, Nickel 13-N-15: 724350 E, 3871500 N, zone 15 (From_Lipo sheet).
NICKEL_UTM = (724350.0, 3871500.0, 15)


def _nad27_reading(e, n, z):
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(f"EPSG:267{z}", "EPSG:4326", always_xy=True).transform(e, n)
    return lat, lon


def _m(lat1, lon1, lat2, lon2):
    return float(np.hypot((lat1 - lat2) * 111320,
                          (lon1 - lon2) * 111320 * np.cos(np.radians(lat1))))


@needs_xy
def test_demonstrated_misreading_moves_to_the_nad27_reading():
    raw = C.read_assemblage_xy(XY, verbose=False)
    # a loader without the datum step: same file, corrections only
    import pandas as pd
    plain = pd.read_csv(XY, sep="\t")
    plain["Assemblages"] = plain["Assemblages"].str.strip()
    r = raw.set_index("Assemblages"); p = plain.set_index("Assemblages")
    lat, lon = _nad27_reading(*NICKEL_UTM)
    assert _m(r.loc["Nickel", "Latitude"], r.loc["Nickel", "Longitude"], lat, lon) < 1.0
    shift = _m(p.loc["Nickel", "Latitude"], p.loc["Nickel", "Longitude"],
               r.loc["Nickel", "Latitude"], r.loc["Nickel", "Longitude"])
    assert 200 < shift < 220                      # the NAD27->NAD83 offset here
    assert r.loc["Nickel", "Latitude"] > p.loc["Nickel", "Latitude"]   # northward


@needs_xy
def test_near_and_unmatched_coordinates_stand():
    raw = C.read_assemblage_xy(XY, verbose=False).set_index("Assemblages")
    import pandas as pd
    plain = pd.read_csv(XY, sep="\t")
    plain["Assemblages"] = plain["Assemblages"].str.strip()
    p = plain.set_index("Assemblages")
    for name in ["Parkin", "Grant", "Kent_Place"]:
        assert _m(raw.loc[name, "Latitude"], raw.loc[name, "Longitude"],
                  p.loc[name, "Latitude"], p.loc[name, "Longitude"]) < 0.01
    t = C.pfg_datum_table(plain)
    assert t.loc["Parkin", "klass"] == "near"          # 41 m off: not demonstrated
    assert t.loc["Grant", "klass"] == "no_pfg_utm"
    assert t.loc["Nickel", "klass"] == "converted"


@needs_xy
def test_correction_transcribed_from_pfg_utm_is_converted_too():
    raw = C.read_assemblage_xy(XY, verbose=False).set_index("Assemblages")
    lat, lon = _nad27_reading(721080.0, 3850050.0, 15)   # PFG site table, 13-N-16
    assert _m(raw.loc["Starkley", "Latitude"], raw.loc["Starkley", "Longitude"], lat, lon) < 1.0


@needs_lmv
def test_settlement_rows_convert_by_source():
    import pandas as pd
    frames = [pd.read_excel(LMV, sheet_name=s) for s in S.ZONE_SHEETS]
    plain = S.apply_site_coordinate_corrections(pd.concat(frames, ignore_index=True), verbose=False)
    conv = S.load_lmv(LMV)
    norm = lambda v: S.normalize_grid(v) if isinstance(v, str) else ""
    key = plain["Number"].map(norm)
    ckey = conv["Number"].map(norm)
    b = plain[key == "13-N-15"].iloc[0]; a = conv[ckey == "13-N-15"].iloc[0]
    assert str(b["Source"]).startswith("PFG")
    assert 200 < a["Northing"] - b["Northing"] < 220 and abs(a["Easting"] - b["Easting"]) < 40
    assert a["Datum"].startswith("NAD27")
    other = conv[~conv["Source"].astype(str).str.startswith("PFG")]
    assert (other["Datum"] == "as recorded").all()
    o = other.dropna(subset=["Easting"]).iloc[0]
    same = plain[plain["Number"].astype(str) == str(o["Number"])].iloc[0]
    assert o["Easting"] == same["Easting"] and o["Northing"] == same["Northing"]
