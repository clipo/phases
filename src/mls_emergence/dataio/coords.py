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

For the record, correcting both changed nothing in the result: the silhouette
still picks three clusters, no assemblage changes cluster, and the
between-cluster F_ST stays 0.017942. That is a robustness finding, not a reason
to skip the correction.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

CORRECTIONS = "coordinate_corrections.csv"


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
        return xy

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
    return xy
