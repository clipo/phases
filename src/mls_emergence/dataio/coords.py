"""Assemblage coordinates, with published corrections applied on read.

The assemblage positions in `data/raw/mainfort-pfg-cplXY.txt` drive everything
spatial in this project: the drainage-basin membership rule, the k-means
clusters the cultural F_ST partitions over, the river-distance matrix, and the
maps. Two of them were wrong by kilometres, found on 2026-09-19 by comparing
the file against the Phillips, Ford and Griffin site table and then asking the
lidar which position had anything at it:

  Belle Meade  the file's point sits 4.81 km west of the site, on ground with
               no compact rise over 1 m. NE1/4 SE1/4 S30 T4N R7E, the legal
               description, matches the section PFG records and has a 2.73 m
               rise beside it.
  Starkley     the file's point falls outside section 02N/05E/21, which PFG
               records for the site; the coordinate transcribed from PFG's own
               UTM falls inside it, 2.46 km away.

Corrections live in `data/raw/coordinate_corrections.csv`, one row per
assemblage naming its source, rather than as edits to the raw file: the raw
file is what Mainfort and Lipo compiled, and a silent edit would leave no trace
of which points we moved or why.

Not every assemblage is a PFG site, and those have no township-and-range
description to check against. Their coordinates stand as recorded (author
ruling, 2026-09-19); the four in the basin set are Grant, Holden Lake, Soudan,
and West Mounds.

For the record, correcting both changed nothing in the result on the set then
in use (2026-09-19): the silhouette still picked three clusters, no assemblage
changed cluster, and the between-cluster F_ST stayed 0.017942. That is a
robustness finding, not a reason to skip the correction.

THE DATUM (2026-09-22). Phillips, Ford and Griffin's UTMs are NAD27; the survey
predates NAD83 by decades, and PFG's UTM for Nickel lands in the author's
quarter-quarter only when read as NAD27 (docs/METHODS_DECISIONS.md, "The
analysis matrix and its sources"). The XY file's compilers read those UTMs as
NAD83, which puts each such point about 209 m south of where PFG put it. The
conversion is applied here, on read, and ONLY where the derivation is
demonstrated: the file's coordinate must coincide, within PFG_DATUM_TOL_M, with
PFG's UTM for that assemblage read as NAD83. A file coordinate that is merely
near a PFG UTM (tens of metres off) was derived some other way, and moving it
by 209 m would be a guess; those are printed and left as recorded, as are
assemblages with no PFG UTM at all. `scripts/measure_datum.py` writes the
ledger of which assemblages fall in each class.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

CORRECTIONS = "coordinate_corrections.csv"
PFG_SPATIAL = "mainfort-spatial-data.xlsx"   # sheet From_Lipo: PFG site table UTMs
PFG_DATUM_TOL_M = 10.0   # PFG UTMs are rounded to 10 m; half a cell's diagonal is 7 m


def _planar_m(lat1, lon1, lat2, lon2):
    import math
    return math.hypot((lat1 - lat2) * 111320,
                      (lon1 - lon2) * 111320 * math.cos(math.radians(lat1)))


def load_pfg_utm(path: str | Path | None = None) -> pd.DataFrame:
    """PFG site-table UTMs by assemblage name, from the From_Lipo sheet.

    One row per assemblage (the sheet repeats a few rows verbatim); rows with
    no UTM are dropped. Columns: easting, northing, zone, site_number.
    """
    if path is None:
        path = Path(__file__).resolve().parents[3] / "data" / "raw" / PFG_SPATIAL
    df = pd.read_excel(Path(path), sheet_name="From_Lipo")
    df.columns = [str(c).strip() for c in df.columns]
    # The sheet also carries lower-case copies of these columns; take the
    # originals by name before renaming so no name is duplicated.
    df = df[["Site_Name", "Site Number", "UTM E", "UTM North", "Zone"]]
    df = df.rename(columns={"Site_Name": "assemblage", "UTM E": "easting",
                            "UTM North": "northing", "Zone": "zone",
                            "Site Number": "site_number"})
    df = df.dropna(subset=["assemblage", "easting", "northing", "zone"])
    df["assemblage"] = df["assemblage"].astype(str).str.strip()
    dup = df[df.duplicated("assemblage", keep=False)]
    for name, g in dup.groupby("assemblage"):
        if g[["easting", "northing", "zone"]].nunique().max() > 1:
            raise ValueError(f"From_Lipo sheet gives {name} more than one UTM")
    df = df.drop_duplicates("assemblage").set_index("assemblage")
    bad = sorted(set(df["zone"].astype(int)) - {15, 16})
    if bad:
        raise ValueError(f"PFG UTM zones outside 15/16: {bad}")
    return df[["easting", "northing", "zone", "site_number"]]


def pfg_datum_table(xy: pd.DataFrame, utm: pd.DataFrame | None = None,
                    tol_m: float = PFG_DATUM_TOL_M) -> pd.DataFrame:
    """Classify each assemblage's file coordinate against PFG's UTM.

    For every assemblage that has a PFG UTM, the distance from the file's
    coordinate to that UTM read as NAD83 and read as NAD27 is measured. Class
    "converted" means the file coordinate IS the NAD83 reading (within tol_m)
    and the NAD27 reading is what the assemblage should carry; "near" means
    within 100 m of the NAD83 reading but not within tol_m, so the derivation
    is not demonstrated and the coordinate stands; "other" is everything else
    with a UTM; "no_pfg_utm" has nothing to compare against. Returned with the
    NAD27 reading for every row that has one, so a caller can apply or merely
    report.
    """
    from pyproj import Transformer
    if utm is None:
        utm = load_pfg_utm()
    tr = {}
    rows = []
    for _, r in xy.iterrows():
        name = str(r["Assemblages"]).strip()
        lat, lon = float(r["Latitude"]), float(r["Longitude"])
        if name not in utm.index:
            rows.append(dict(assemblage=name, site_number=None, klass="no_pfg_utm",
                             d_nad83_m=float("nan"), d_nad27_m=float("nan"),
                             lat_nad27=float("nan"), lon_nad27=float("nan")))
            continue
        u = utm.loc[name]
        z = int(u["zone"])
        if z not in tr:
            tr[z] = (Transformer.from_crs(f"EPSG:269{z}", "EPSG:4326", always_xy=True),
                     Transformer.from_crs(f"EPSG:267{z}", "EPSG:4326", always_xy=True))
        lo83, la83 = tr[z][0].transform(float(u["easting"]), float(u["northing"]))
        lo27, la27 = tr[z][1].transform(float(u["easting"]), float(u["northing"]))
        d83 = _planar_m(lat, lon, la83, lo83)
        d27 = _planar_m(lat, lon, la27, lo27)
        klass = "converted" if d83 <= tol_m else ("near" if d83 <= 100 else "other")
        rows.append(dict(assemblage=name, site_number=u["site_number"], klass=klass,
                         d_nad83_m=d83, d_nad27_m=d27, lat_nad27=la27, lon_nad27=lo27))
    return pd.DataFrame(rows).set_index("assemblage")


def apply_pfg_datum(xy: pd.DataFrame, verbose: bool = True,
                    utm_path: str | Path | None = None) -> pd.DataFrame:
    """Move demonstrated NAD83 misreadings of PFG UTMs to their NAD27 reading."""
    table = pfg_datum_table(xy, load_pfg_utm(utm_path))
    out = xy.copy()
    conv = table[table["klass"] == "converted"]
    for name, row in conv.iterrows():
        sel = out["Assemblages"] == name
        before = (float(out.loc[sel, "Latitude"].iloc[0]), float(out.loc[sel, "Longitude"].iloc[0]))
        out.loc[sel, "Latitude"] = row["lat_nad27"]
        out.loc[sel, "Longitude"] = row["lon_nad27"]
        if verbose:
            d = _planar_m(before[0], before[1], row["lat_nad27"], row["lon_nad27"])
            print(f"datum: {name} read PFG UTM as NAD27, moved {d:.0f} m "
                  f"(was {row['d_nad83_m']:.1f} m from the NAD83 reading)")
    if verbose:
        near = table[table["klass"] == "near"]
        for name, row in near.iterrows():
            print(f"datum: {name} left as recorded; {row['d_nad83_m']:.0f} m from PFG's UTM "
                  f"read as NAD83, so its derivation is not demonstrated")
        print(f"datum: {len(conv)} converted, {len(near)} near and left, "
              f"{int((table['klass'] == 'other').sum())} other, "
              f"{int((table['klass'] == 'no_pfg_utm').sum())} without a PFG UTM")
    return out


def load_corrections(path: str | Path | None = None) -> pd.DataFrame:
    """The correction table, indexed by assemblage name as the XY file spells it."""
    if path is None:
        path = Path(__file__).resolve().parents[3] / "data" / "raw" / CORRECTIONS
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["latitude", "longitude", "source", "note"])
    df = pd.read_csv(path)
    df["assemblage"] = df["assemblage"].astype(str).str.strip()
    return df.drop_duplicates("assemblage").set_index("assemblage")


def read_assemblage_xy(path: str | Path, corrections: str | Path | None = None,
                       verbose: bool = True) -> pd.DataFrame:
    """The assemblage coordinate table with corrections applied.

    Returns the file's own columns (Assemblages, Latitude, Longitude) so that
    callers keep their existing handling; only the values change. Applied
    corrections are printed, never applied silently, and a correction naming an
    assemblage the file does not contain raises rather than passing unnoticed —
    a typo in the correction table would otherwise look like a clean run.
    """
    xy = pd.read_csv(Path(path), sep="\t")
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    corr = load_corrections(corrections)
    if corr.empty:
        return apply_pfg_datum(xy, verbose=verbose)

    known = set(xy["Assemblages"])
    missing = [a for a in corr.index if a not in known]
    if missing:
        raise ValueError(
            f"coordinate corrections name assemblages absent from {Path(path).name}: "
            f"{missing}")

    for name, row in corr.iterrows():
        sel = xy["Assemblages"] == name
        before = (float(xy.loc[sel, "Latitude"].iloc[0]),
                  float(xy.loc[sel, "Longitude"].iloc[0]))
        xy.loc[sel, "Latitude"] = float(row["latitude"])
        xy.loc[sel, "Longitude"] = float(row["longitude"])
        if verbose:
            import math
            d = math.hypot((float(row["latitude"]) - before[0]) * 111320,
                           (float(row["longitude"]) - before[1]) * 111320
                           * math.cos(math.radians(before[0])))
            print(f"coordinate correction: {name} moved {d / 1000:.2f} km "
                  f"({before[0]:.5f},{before[1]:.5f} -> "
                  f"{row['latitude']:.5f},{row['longitude']:.5f}) [{row['source']}]")
    # Corrections first, then the datum: a correction transcribed from PFG's
    # own UTM (Starkley) is itself a NAD83 reading and is converted like the rest.
    return apply_pfg_datum(xy, verbose=verbose)
