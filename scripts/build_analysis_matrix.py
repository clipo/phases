#!/usr/bin/env python3
"""Build the decorated-class analysis matrix from the two definitive sources.

THE SOURCES, by the author's ruling of 2026-09-21.

  data/raw/PFGData.xlsx, sheet SherdData
      The definitive Phillips-Ford-Griffin counts. Parkin carries 865 Parkin
      Punctated here, not the 1,298 of Mainfort's table.
  data/raw/mainfort-pfg-cpl.xlsx, sheet cpl-pfg
      Lipo's (2001) compilation: the PFG counts plus his own collections. It
      adds four substantial field collections (Beck +1,154 decorated sherds,
      Belle Meade +1,455, Cramor +286, Castile Landing +220), a few sherds at
      three other sites, and one new assemblage, Holden Lake.

Mainfort's (2003) table is NOT a source for this matrix. His counts exceed the
PFG counts at ten sites and may include collections of his own, so they are
written out separately, as published, to data/processed/mainfort2003_matrix.csv
and used as a comparison, never mixed in.

THE RULE. One source per assemblage, never a sum. An assemblage in the cpl-pfg
tab takes that row. Otherwise it takes its PFGData row, matched on the exact
normalized site name. An assemblage in neither is left out of the primary
matrix and listed, since its only counts are Mainfort's.

THE CLASS MAPPING is recovered, not assumed. Collapsing PFGData as

    Barton/Kent/MPI = Barton Incised + Kent Incised + Mound Place Incised (T)
    Painted         = Old Town Red + Nodena Red and White + Avenue Polychrome

reproduces nine of the sixteen cpl-pfg rows that have a PFGData row exactly,
and leaves no negative cell in the other seven, which is what "PFG plus
additions" requires. Every other candidate tried reproduces fewer. The script
re-checks this on every run and refuses to write if a negative cell appears.

ASSEMBLAGE IDENTITIES RULED ON BY THE AUTHOR, 2026-09-21.

  Wall     is Walls with a missing "s": a second row for one site in the old
           Mainfort sheet, carrying Walls' coordinates to five decimals. It is
           not in this matrix; Walls comes from the cpl-pfg tab (13-P-1).
  Soudan   is PFG 13-N-1 ("Sudan Cemetery" in PFGData). The PFG collection there
           holds 5 decorated sherds (3 Parkin Punctated, 2 Barton/Kent/MPI), so
           it is left out on sample size, not for want of a source. The 1,311
           sherds Mainfort reports for Soudan are later work.
  Nickel   is 13-N-15, which is how the cpl-pfg tab keys it. PFGData splits it
           across "13-N-15/B" and "13-N-15/C,D,E" and spells it "Nickle".

Composite site numbers matter when reading PFGData. Kent Place appears as
"13-N-4/B,C,D,E" and again as "13-N-4/C", so collection C is in both rows and
summing them double-counts it; the cpl-pfg tab's 241 equals the first row
alone. Matched on base site number, the tab's additions to PFGData are seven
assemblages, the number Lipo (2001, Ch. 5) reports: Belle Meade +1,455
decorated sherds, Beck +1,154, Nickel +1,120, Rose Mound +738, Cramor +286,
Castile Landing +220, and Holden Lake, which PFG never collected.

One location is NOT settled. The PFG site table gives 13-N-1 and 13-N-15 the
identical section and UTM (T4N R5E S11; 3871500 N, 724350 E), so one of those
rows is a copy of the other. The settlement compilation resolves it as Nickel
there and Soudan 15 km south-west at House's (1991) location, which is what
the coordinate file uses. Soudan is not in this matrix, so nothing here
depends on it, but Nickel's position does rest on that reading.

WHAT THIS REPLACES. The `pfg-cpl-mainfort` sheet, which every result rested on
until 2026-09-21, summed the cpl-pfg tab with a "Mainfort" sheet whose Parkin
Punctated column was doubled (scripts/check_matrix_against_mainfort2003.py).
The same PFG sherds were counted up to three times over.

Usage:
    python scripts/build_analysis_matrix.py [--check]
"""
from __future__ import annotations

import contextlib
import io
import re
import runpy
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_CSV = ROOT / "data" / "raw" / "mainfort-pfg-cpl.csv"
TABLE = ROOT / "data" / "raw" / "mainfort2003_table1.csv"
OUT = ROOT / "data" / "processed" / "analysis_matrix.csv"
PROV = ROOT / "data" / "processed" / "analysis_matrix_provenance.csv"
ALIAS = {"40la7": "40la007", "richardsonslg": "richardsonslanding",
         "castilelg": "castilelanding", "bartonranch": "barton",
         "chucalissa": "chuccalissa"}
UNTALLIED_BY_MAINFORT = ("Wallace_Incised", "Hull_Engraved")


def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


PFG_XLSX = ROOT / "data" / "raw" / "PFGData.xlsx"
WORKBOOK = ROOT / "data" / "raw" / "mainfort-pfg-cpl.xlsx"
MAINFORT_OUT = ROOT / "data" / "processed" / "mainfort2003_matrix.csv"
BK = ["Barton Incised", "Kent Incised", "Mound Place Incised (T)"]
PAINTED = ["Old Town Red", "Nodena Red and White", "Avenue Polychrome"]
ONE_TO_ONE = {"Parkin_Punctated": "Parkin Punctated", "Fortune_Noded": "Fortune Noded",
              "Ranch_Incised": "Ranch Incised", "Walls_Engraved": "Walls Engraved",
              "Wallace_Incised": "Wallace Incised", "Rhodes_Incised": "Rhodes Incised",
              "Vernon_Paul_Applique": "Vernon Paul Applique", "Hull_Engraved": "Hull Engraved"}


def pfg_table(TY) -> pd.DataFrame:
    """PFGData collapsed into the analysis classes, one row per site number."""
    d = pd.read_excel(PFG_XLSX, sheet_name="SherdData")
    d.columns = [str(c).strip() for c in d.columns]
    num = lambda c: pd.to_numeric(d[c], errors="coerce").fillna(0)
    out = pd.DataFrame({"id": d["Site Number"].astype(str).str.strip(),
                        "name": d["Site Name"].astype(str).str.strip()})
    for cls, col in ONE_TO_ONE.items():
        out[cls] = num(col)
    out["Barton/Kent/MPI"] = sum(num(c) for c in BK)
    out["Painted"] = sum(num(c) for c in PAINTED)
    return out[["id", "name"] + TY]


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    import os
    here = os.getcwd()
    os.chdir(ROOT)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            g = runpy.run_path(str(ROOT / "scripts" / "reconcile_ceramic_matrix.py"))
    finally:
        os.chdir(here)
    u, TY, cpl = g["u"], g["TY"], g["cpl"]
    P = pfg_table(TY)

    # The mapping is a claim about the data, so it is tested every time.
    C = pd.read_excel(WORKBOOK, sheet_name="cpl-pfg")
    C["id"] = C["Assemblages"].astype(str).str.strip()
    C = C.drop_duplicates("id").set_index("id")[TY].apply(pd.to_numeric, errors="coerce").fillna(0)
    Pid = P.groupby("id")[TY].sum()
    shared = [i for i in C.index if i in Pid.index]
    diff = C.loc[shared] - Pid.loc[shared]
    if (diff < 0).any().any():
        bad = diff[(diff < 0).any(axis=1)]
        raise ValueError("cpl-pfg holds FEWER sherds than PFGData in some class, so the "
                         f"class mapping no longer holds:\n{bad}")
    exact = int((diff.abs().sum(axis=1) == 0).sum())

    by_name = P.assign(k=P["name"].map(norm)).groupby("k")[TY].sum()
    # Assemblages PFGData lists under another name. Each entry is verified by
    # location, not by resemblance: "West" and "West Mounds" are different
    # sites in different states, and that mistake has been made here before.
    #   Grant -> 13-N-11 "Grant Place", 0.1 km from the assemblage's coordinate.
    PFG_NAME = {"grant": "grantplace"}
    rows, prov, left_out = [], [], []
    for _, r in u.iterrows():
        name, k = str(r["Assemblages"]).strip(), r["k"]
        Cv = cpl(k).astype(float)
        if Cv.sum() > 0:
            row, src = Cv, "Lipo 2001 compilation (cpl-pfg tab: PFG plus his collections)"
        elif PFG_NAME.get(k, k) in by_name.index and by_name.loc[PFG_NAME.get(k, k)].sum() > 0:
            row, src = by_name.loc[PFG_NAME.get(k, k)].to_numpy(float), "PFGData.xlsx"
        else:
            left_out.append(name)
            continue
        rows.append([name] + [int(v) for v in row])
        prov.append(dict(assemblage=name, source_used=src,
                         previous_total=int(r[TY].astype(float).sum()),
                         corrected_total=int(row.sum())))
    prov = pd.DataFrame(prov)
    prov.attrs["left_out"] = left_out
    prov.attrs["mapping_exact"] = (exact, len(shared))
    return pd.DataFrame(rows, columns=["Assemblages"] + TY), prov


def mainfort_comparison() -> pd.DataFrame:
    """Mainfort (2003) Table 1 exactly as published, for comparison only."""
    return pd.read_csv(TABLE, comment="#")


def main() -> int:
    mat, prov = build()
    if "--check" in sys.argv:
        ok = OUT.exists() and pd.read_csv(OUT).equals(mat)
        print("analysis matrix is current." if ok else "analysis matrix is STALE.")
        return 0 if ok else 1
    OUT.parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(OUT, index=False)
    prov.to_csv(PROV, index=False)
    mainfort_comparison().to_csv(MAINFORT_OUT, index=False)
    e, n = prov.attrs["mapping_exact"]
    print(f"class mapping: {e} of {n} cpl-pfg rows reproduce PFGData exactly, no negative cells")
    print(prov.source_used.value_counts().to_string())
    print(f"\nleft out of the primary matrix (Mainfort is their only source): "
          f"{len(prov.attrs['left_out'])}\n  " + ", ".join(prov.attrs["left_out"]))
    print(f"\nsherds: {prov.previous_total.sum():,} in these rows before -> "
          f"{prov.corrected_total.sum():,}")
    print(f"wrote {OUT.relative_to(ROOT)}, {PROV.relative_to(ROOT)} and "
          f"{MAINFORT_OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
