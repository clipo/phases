"""build_mainfort2003_sites.py - coordinates for the 39 sites of Mainfort (2003) Table 1.

Writes `data/processed/mainfort2003_sites.csv`, which `dataio/matrix.py`
reads for the replication analysis (83). The file was first assembled by hand
on 2026-09-21; this script reproduces it from the three sources it named, row
by row, so that it can be regenerated (rule 15) and so that the four rows
taken from the assemblage coordinate file pass through `read_assemblage_xy`
and receive the published corrections and the NAD27 conversion like every
other use of that file.

Sources, in the order tried for each site:
  1. an explicit assignment in ASSEMBLAGE_FILE (four sites whose coordinate
     the 2026-09-21 assembly took from `mainfort-pfg-cplXY.txt`);
  2. the `MainfortXY` sheet of `data/raw/mainfort-pfg-cpl.xlsx`, Mainfort's
     own coordinates, by the name alias in MAINFORT_XY_NAME;
  3. the `Spatial Data` sheet (settlement compilation), for Sweat.
A site found nowhere raises. Precision follows the source: Mainfort's own
points carry three decimals, and are written as he gave them.

Usage: .venv/bin/python scripts/build_mainfort2003_sites.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from mls_emergence.dataio.coords import read_assemblage_xy  # noqa: E402

TABLE = ROOT / "data" / "raw" / "mainfort2003_table1.csv"
WORKBOOK = ROOT / "data" / "raw" / "mainfort-pfg-cpl.xlsx"
XY = ROOT / "data" / "raw" / "mainfort-pfg-cplXY.txt"
OUT = ROOT / "data" / "processed" / "mainfort2003_sites.csv"

# Table 1 name -> assemblage name in mainfort-pfg-cplXY.txt
ASSEMBLAGE_FILE = {"Cramor": "Cramor_Place", "Williamson": "Williamson",
                   "Barton Ranch": "Barton_Ranch", "Kent": "Kent_Place"}
# Table 1 name -> Site_Name in the MainfortXY sheet, where they differ
MAINFORT_XY_NAME = {"40LA7": "40LA007", "Richardson's Lg.": "Richardsons_Landing",
                    "Jones Bayou": "Jones_Bayou", "Lake Cormorant": "Lake_Cormorant",
                    "Upper Nodena": "Upper_Nodena", "Castile Lg.": "Castile_Landing",
                    "Rose Mound": "Rose_Mound", "Neeley's Ferry": "Neeleys_Ferry",
                    "Carson Lake": "Carson_Lake", "Mound Place": "Mound_Place",
                    "Chucalissa": "Chuccalissa", "Graves Lake": "Graves_Lake",
                    "Clay Hill": "Clay_Hill", "Belle Meade": "Belle_Meade"}
SPATIAL_DATA = {"Sweat": "Sweat"}

SRC_XY = "mainfort-pfg-cplXY.txt (with published corrections and the NAD27 conversion applied on read)"
SRC_MXY = "MainfortXY sheet of mainfort-pfg-cpl.xlsx"
SRC_SD = "Spatial Data sheet of mainfort-pfg-cpl.xlsx (settlement compilation)"


def main():
    table = pd.read_csv(TABLE, comment="#")
    xy = read_assemblage_xy(XY, verbose=False).set_index("Assemblages")
    mxy = pd.read_excel(WORKBOOK, sheet_name="MainfortXY")
    mxy["Site_Name"] = mxy["Site_Name"].astype(str).str.strip()
    mxy = mxy.drop_duplicates("Site_Name").set_index("Site_Name")
    sd = pd.read_excel(WORKBOOK, sheet_name="Spatial Data")
    sd["Name"] = sd["Name"].astype(str).str.strip()
    sd = sd.drop_duplicates("Name").set_index("Name")

    # Sites whose coordinate the project has corrected (coordinate_corrections.csv,
    # keyed by assemblage-file name) take the corrected, datum-converted value
    # from the assemblage file rather than Mainfort's own point (2026-09-23;
    # until then Beck, Belle Meade and Starkley carried ruled-wrong points here).
    from mls_emergence.dataio.coords import load_corrections
    corrected = {n.replace("_", " "): n for n in load_corrections().index}
    rows = []
    for _, r in table.iterrows():
        site = str(r["site"]).strip()
        if site in corrected and site not in ASSEMBLAGE_FILE:
            a = corrected[site]
            rows.append(dict(site=site, phase_mainfort2003=r["phase"], latitude=float(xy.loc[a, "Latitude"]),
                             longitude=float(xy.loc[a, "Longitude"]),
                             coordinate_source=SRC_XY + "; project coordinate correction"))
            continue
        if site in ASSEMBLAGE_FILE:
            a = ASSEMBLAGE_FILE[site]
            lat, lon, src = float(xy.loc[a, "Latitude"]), float(xy.loc[a, "Longitude"]), SRC_XY
        elif site in SPATIAL_DATA:
            lat, lon, src = float(sd.loc[SPATIAL_DATA[site], "latitude"]), float(sd.loc[SPATIAL_DATA[site], "longitude"]), SRC_SD
        else:
            key = MAINFORT_XY_NAME.get(site, site)
            if key not in mxy.index:
                raise KeyError(f"{site!r}: no coordinate in any source (MainfortXY key {key!r})")
            lat, lon, src = float(mxy.loc[key, "Latitude"]), float(mxy.loc[key, "Longitude"]), SRC_MXY
        rows.append(dict(site=site, phase_mainfort2003=r["phase"], latitude=lat, longitude=lon,
                         coordinate_source=src))
    out = pd.DataFrame(rows)
    header = ("# Coordinates for the 39 sites of Mainfort (2003) Table 1, written by scripts/build_mainfort2003_sites.py.\n"
              "# Each row names where its coordinate came from. The MainfortXY sheet is Mainfort's own; the other two sources are this project's.\n"
              "# Not independently verified against legal descriptions except where the assemblage file carries a published correction.\n")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(header + out.to_csv(index=False), encoding="utf-8")
    print(f"wrote {OUT} ({len(out)} sites): "
          + ", ".join(f"{k} {v}" for k, v in out["coordinate_source"].str.split(" ").str[0].value_counts().items()))


if __name__ == "__main__":
    main()
