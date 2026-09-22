#!/usr/bin/env python3
"""Confront the workbook's "Mainfort" sheet with Mainfort's published table.

`data/raw/mainfort-pfg-cpl.xlsx` sheet `mainfort-data-collapsed` is one of the
two components the analysis matrix sums. This script compares it, cell by cell,
with Table 1 of Mainfort (2003, Southeastern Archaeology 22(2):178), transcribed
into `data/raw/mainfort2003_table1.csv`.

What it found on 2026-09-21:

  - Ranch Incised, Rhodes Incised, Walls Engraved, Vernon Paul Applique and
    Fortune Noded match the published table in 38 of 38 shared sites. The sheet
    IS Mainfort's table, and the transcription is sound.
  - PARKIN PUNCTATED IS DOUBLED. The sheet holds exactly twice the published
    count in 34 of 38 sites and twice plus one to nine sherds in the other
    four. Eight of the sheet's sites that are not in the 2003 table are exactly
    twice the PFG sheet's count as well, so the doubling runs through the whole
    sheet. The PFG sheet carries both a "Parkin Punctated" and a "Parkin
    Incised" column; a collapse meant to add the two that pointed at the first
    twice would produce exactly this.
  - Sixteen sites in the sheet are absent from the 2003 table, including all
    four Parchman-cluster sites. Mainfort (2003) keeps only sites with 800 or
    more sherds and says his data are those of Mainfort (1999), so these rows
    come from somewhere else and their source is not yet recorded.

And one thing the paper itself establishes (rule 12, read from the PDF): his
counts "were transcribed from files of the Lower Mississippi Valley Survey",
which is the Phillips-Ford-Griffin material Lipo's compilation also carries.
The two components are NOT independent, and the matrix sums them.

Usage:
    python scripts/check_matrix_against_mainfort2003.py
Exits non-zero while the doubling is present, so it can gate a rerun.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / "raw" / "mainfort-pfg-cpl.xlsx"
TABLE = ROOT / "data" / "raw" / "mainfort2003_table1.csv"
OUT = ROOT / "output" / "findings" / "matrix_integrity.md"
ALIAS = {"40la7": "40la007", "richardsonslg": "richardsonslanding",
         "castilelg": "castilelanding", "bartonranch": "barton",
         "chucalissa": "chuccalissa"}
SAME = [("ranch_inc", "Ranch_Incised"), ("rhodes_inc", "Rhodes_Incised"),
        ("walls_eng", "Walls_Engraved"), ("vernon_pal", "Vernon_Paul_Applique"),
        ("fortune_nd", "Fortune_Noded")]


def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main() -> int:
    pub = pd.read_csv(TABLE, comment="#")
    pub["k"] = pub.site.map(norm).replace(ALIAS)
    wb = pd.read_excel(XLSX, sheet_name="mainfort-data-collapsed")
    wb["k"] = wb.site_name.map(norm)
    b = pub.merge(wb, on="k", suffixes=("_pub", "_wb"))
    only_wb = sorted(set(wb.site_name[~wb.k.isin(pub.k)]))
    only_pub = sorted(set(pub.site[~pub.k.isin(wb.k)]))
    ratio = (b.Parkin_Punctated_wb / b.Parkin_Punctated_pub)
    exact2 = int((b.Parkin_Punctated_wb == 2 * b.Parkin_Punctated_pub).sum())
    agree = int((b.Parkin_Punctated_wb == b.Parkin_Punctated_pub).sum())

    typed = pub[list(pub.columns[2:-2])].sum(axis=1)
    over = list(pub.site[typed > pub.Total])

    L = ["# The workbook's Mainfort sheet against Mainfort (2003), Table 1", "",
         f"{len(b)} sites appear in both; {len(only_wb)} only in the workbook sheet; "
         f"{len(only_pub)} only in the published table.", "",
         "| class | sites identical to the published table |", "|---|---|"]
    for w, p in SAME:
        L.append(f"| {p} | {int((b[w] == b[p]).sum())} of {len(b)} |")
    L += [f"| **Parkin Punctated** | **{agree} of {len(b)}**; exactly twice the published "
          f"count in {exact2}, ratio {ratio.min():.2f} to {ratio.max():.2f} in all |", "",
          "Workbook-only sites, whose source is not yet recorded: "
          + ", ".join(only_wb) + ".", "",
          "Transcription self-check: rows whose typed columns exceed the printed row total "
          "(Mainfort's totals include sherds he does not list, so columns may fall short "
          f"but never exceed): {over if over else 'none'}.", ""]
    # Mainfort's published table against the definitive PFG counts, site by site.
    pfg = pd.read_excel(ROOT / "data" / "raw" / "PFGData.xlsx", sheet_name="SherdData")
    pfg.columns = [str(c).strip() for c in pfg.columns]
    pfg["k"] = pfg["Site Name"].map(norm)
    g = pfg.groupby("k").agg(pp=("Parkin Punctated", "sum"), tot=("Sherd Total", "sum"))
    pk = pub.assign(k=pub.site.map(norm).replace({"castilelg": "castile"}))
    j = pk.merge(g, left_on="k", right_index=True)
    same = int((j.Parkin_Punctated == j.pp).sum())
    more = j[j.Parkin_Punctated > j.pp]
    fewer = int((j.Parkin_Punctated < j.pp).sum())
    dec = ["Parkin_Punctated", "Barton_Incised", "Kent_Incised", "Ranch_Incised",
           "Old_Town_Red", "Nodena_Red_and_White", "Rhodes_Incised", "Walls_Engraved",
           "Vernon_Paul_Applique", "Fortune_Noded"]
    plain = int(pub.Mississippi_Plain.sum() + pub.Bell_Plain.sum())
    decorated = int(pub[dec].sum().sum())
    L += ["## Mainfort's published counts against PFGData.xlsx", "",
          f"{len(j)} of his sites match a PFGData site by exact name. Parkin Punctated is "
          f"identical at {same}, larger in his table at {len(more)}, and smaller at {fewer}.",
          "Sites where his count exceeds the survey's: "
          + ", ".join(f"{r.site} ({int(r.Parkin_Punctated)} against {int(r.pp)})"
                      for r in more.itertuples()) + ".", "",
          f"Plainware in his table: {plain:,} Mississippi Plain and Bell Plain sherds against "
          f"{decorated:,} decorated, that is {100 * plain / (plain + decorated):.0f} percent "
          "plain.", ""]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    if exact2 or agree < len(b):
        print("DEFECT PRESENT: Parkin Punctated in the workbook sheet does not match the "
              "published table.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
