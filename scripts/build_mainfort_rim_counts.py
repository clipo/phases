#!/usr/bin/env python3
"""Turn Mainfort's rim-attribute file into a table of counts, and say what each column is.

THE TWO SOURCES.
  data/raw/mainfort_miss_rims.xls
      Sent by Robert C. Mainfort, Jr. to Carl Lipo, to use as he wishes (author,
      2026-10-04; acknowledge him). 55 assemblages, 13 attribute columns, a
      total, and rough coordinates. Every attribute value is a PERCENTAGE.
  data/raw/mainfort2003_rim_table1_counts.csv
      Table 1 of Mainfort (2003b), "Late Period Ceramic Rim Attribute
      Variation in the Central Mississippi Valley", Southeastern Archaeology
      22(1):33-46: the published COUNTS of six rim attributes for the same 55
      sites and the same 7,195 rims. Transcribed from the page image.

WHAT THE COLUMNS ARE, established by comparing the two (this script repeats the
comparison every time and stops if it fails).
  total      N, the number of rims analyzed. Equal to the published N at all 55
             sites.
  beadrim    beaded applique strip      } each is the percentage of the N rims
  noderim    band of nodes              } that carry the attribute. They are
  widenode   widely spaced nodes        } separate yes-or-no observations, not
  notchtop   notched top of the lip     } parts of one whole, so they do not
  notchext   notched exterior           } sum to 100 with one another.
  sharpbev, modbev, slitebev, nobev
             degree of interior bevel. These four ARE parts of one whole and
             sum to 100; the first three together are the paper's "Beveled".
             The paper publishes only beveled or not; the degrees are in the
             file alone (recorded as in Mainfort and Moore 1998).
  flatlip, roundlip
             lip form. Parts of one whole, summing to 100. Not in the paper.
  nf, bell   the percentage of sherds in each paste, Neeley's Ferry and Bell.
             Not in the paper. The author's reading (2026-10-04): the paste of
             the rim sherds on which the attributes were recorded, so not
             comparable with the Phillips-Ford-Griffin type counts, which
             needed only a surface and were made on body sherds too. The file
             agrees: at all but two sites each value is a whole number of rims
             out of N, and the values are far from the sherd percentages of
             the same sites in Mainfort (2003a: Table 1). They do NOT sum to
             100 (20 to 102, median 89); what the remaining rims are is not
             established.

THE COUNTS WRITTEN.
  six attributes   the PUBLISHED counts, as printed, which are the record
                   (author, 2026-10-04: use the paper's table; the file is
                   older and may have been revised since). Where a printed
                   count and the file disagree the printed count is still
                   used, and the disagreement is listed. One printed value
                   is probably a misprint: Pecan Point's beaded count is 10 of
                   30 rims, where the paper's own Table 2 and the file give
                   10 percent, which is 3. It is kept as printed and flagged.
  bevel degrees    the published beveled count divided among sharp, moderate
                   and slight in the file's proportions (largest remainders);
                   none = N - beveled.
  lip form         flat = file percentage x N / 100, rounded; round = N - flat.
  plainwares       file percentage x N / 100, rounded, with the remainder as
                   "other". Marked unconfirmed.

Outputs: data/processed/mainfort_rim_counts.csv, output/findings/rim_table_structure.md
Usage: .venv/bin/python scripts/build_mainfort_rim_counts.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLS = ROOT / "data" / "raw" / "mainfort_miss_rims.xls"
PUBLISHED = ROOT / "data" / "raw" / "mainfort2003_rim_table1_counts.csv"
OUT_CSV = ROOT / "data" / "processed" / "mainfort_rim_counts.csv"
OUT_MD = ROOT / "output" / "findings" / "rim_table_structure.md"

# file name -> name in the published table
TO_PUBLISHED = {"Walls/Harris": "Walls", "CastileLndg": "Castile Landing", "RoseMound": "Rose Mound",
                "NeeleysFerry": "Neeleys Ferry", "CT40": "Neely (3CT40)", "CT245": "3CT245",
                "CT7": "Bradley (3CT7)", "Shawnee Vill": "Shawnee Village", "Richardson": "Richardson's Landing",
                "MidNodena": "Middle Nodena", "GravesLake": "Graves Lake", "UpNodena": "Upper Nodena",
                "LA17": "40LA17", "Knappenberg": "Knappenberger", "Denton Mnds": "Denton Mounds"}
MODES = [("beaded", "beadrim"), ("noded", "noderim"), ("wide_nodes", "widenode"),
         ("notched_top", "notchtop"), ("notched_exterior", "notchext")]
TOLERANCE = 1.5      # percentage points between a published count / N and the file's percentage


def largest_remainder(shares: np.ndarray, total: int) -> np.ndarray:
    """Integers summing to `total`, in proportion to `shares`."""
    if total == 0 or shares.sum() == 0:
        return np.zeros(len(shares), int)
    raw = shares / shares.sum() * total
    out = np.floor(raw).astype(int)
    for i in np.argsort(-(raw - out))[: total - out.sum()]:
        out[i] += 1
    return out


def paste_check(m: pd.DataFrame) -> list[str]:
    """Evidence on what nf and bell are percentages of."""
    n = m["N"].to_numpy(int)

    def whole(col):
        ok = 0
        for p, k in zip(m[col].to_numpy(int), n):
            v = 100.0 * np.arange(k + 1) / k
            ok += bool(((np.floor(v + 0.5) == p) | (np.round(v) == p)).any())
        return ok

    L = ["", "## What the paste percentages are percentages of", "",
         ("The author reads `nf` and `bell` as the paste of the rim sherds on which the attributes were "
         "recorded. Two checks agree that the base is the N rims and not the whole sherd assemblage."), "",
         (f"- Each value is a whole number of rims out of N at {whole('nf')} of {len(m)} sites for `nf` and "
         f"{whole('bell')} for `bell`.")]
    sherds = ROOT / "data" / "raw" / "mainfort2003_table1.csv"
    if sherds.exists():
        t = pd.read_csv(sherds, comment="#")
        alias = {"Castile Landing": "Castile Lg.", "Neeleys Ferry": "Neeley's Ferry",
                 "Richardson's Landing": "Richardson's Lg."}
        j = m.assign(k=m["published_name"].map(lambda s: alias.get(s, s))).merge(
            t, left_on="k", right_on="site", suffixes=("", "_sherd"))
        bell_all = 100.0 * j["Bell_Plain"] / j["Total"]
        pk = j["k"] == "Parkin"
        L.append(f"- At the {len(j)} sites also in Mainfort (2003a: Table 1), Bell Plain's share of all sherds "
                 f"there differs from this file's `bell` by a mean of {np.abs(j['bell'] - bell_all).mean():.0f} "
                 f"points. Parkin is {int(j.loc[pk, 'bell'].iloc[0])} percent Bell here and "
                 f"{bell_all[pk].iloc[0]:.0f} percent there. Rims and body sherds are different samples of "
                 "the pottery, so the two need not agree.")
    return L


def main() -> int:
    x = pd.read_excel(XLS)
    x.columns = [str(c).strip() for c in x.columns]
    x["site"] = x["site"].astype(str).str.strip()
    x["published_name"] = x["site"].map(lambda s: TO_PUBLISHED.get(s, s))
    t = pd.read_csv(PUBLISHED, comment="#")
    m = x.merge(t, left_on="published_name", right_on="site", how="outer", suffixes=("", "_pub"), indicator=True)
    if not (m["_merge"] == "both").all():
        raise SystemExit("the file and the published table do not list the same sites: "
                         + str(m[m["_merge"] != "both"][["site", "site_pub"]].to_dict("records")))
    if not (m["total"] == m["N"]).all():
        raise SystemExit("the file's totals are not the published N at every site")

    notes = []
    rows = []
    for r in m.itertuples(index=False):
        n = int(r.N)
        row = {"site": r.published_name, "file_name": r.site, "N": n}
        for pub, col in MODES:
            printed = int(getattr(r, pub))
            pct = float(getattr(r, col))
            if abs(100.0 * printed / n - pct) > TOLERANCE:
                notes.append(f"{r.published_name}, {pub.replace('_', ' ')}: printed count {printed} of {n} "
                             f"({100.0 * printed / n:.1f} percent) against {pct:.0f} percent in the file "
                             f"(which would be {int(round(pct * n / 100.0))}); the printed count is used")
            row[pub] = printed
        bev_printed = int(r.beveled)
        bev_pct = float(r.sharpbev + r.modbev + r.slitebev)
        four = float(bev_pct + r.nobev)
        if abs(four - 100.0) > 2.0:
            notes.append(f"{r.published_name}, bevel degrees: the file's four bevel percentages sum to {four:.0f}, "
                         f"not 100; the published beveled count ({bev_printed} of {n}) is used and divided in "
                         f"the file's proportions")
        elif abs(100.0 * bev_printed / n - bev_pct) > TOLERANCE:
            notes.append(f"{r.published_name}, beveled: published {bev_printed} of {n} "
                         f"({100.0 * bev_printed / n:.1f} percent) against {bev_pct:.0f} percent in the file; "
                         f"the published count is used")
        deg = largest_remainder(np.array([r.sharpbev, r.modbev, r.slitebev], float), bev_printed)
        row.update(beveled=bev_printed, bevel_sharp=int(deg[0]), bevel_moderate=int(deg[1]),
                   bevel_slight=int(deg[2]), bevel_none=n - bev_printed)
        flat = int(round(float(r.flatlip) / float(r.flatlip + r.roundlip) * n))
        row.update(lip_flat=flat, lip_round=n - flat)
        nf, bell = int(round(r.nf * n / 100.0)), int(round(r.bell * n / 100.0))
        if nf + bell > n:                       # rounding (Rast sums to 102 percent)
            nf = n - bell
        row.update(nf_percent=int(r.nf), bell_percent=int(r.bell), paste_nf=nf, paste_bell=bell,
                   paste_other=n - nf - bell, latitude_rough=r.latitude, longitude_rough=r.longitude)
        rows.append(row)
    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    lip = m["flatlip"] + m["roundlip"]
    bev4 = m[["sharpbev", "modbev", "slitebev", "nobev"]].sum(axis=1)
    paste = m["nf"] + m["bell"]
    L = ["# What Mainfort's rim-attribute file holds, as counts", "",
         ("Written by `scripts/build_mainfort_rim_counts.py`. The file "
         "(`data/raw/mainfort_miss_rims.xls`, sent by R. C. Mainfort, Jr.) is the data behind Mainfort (2003b), "
         "\"Late Period Ceramic Rim Attribute Variation in the Central Mississippi Valley\", *Southeastern "
         "Archaeology* 22(1):33-46. That paper's Table 1 prints the counts of six attributes for the same "
         f"{len(m)} sites and {int(m['N'].sum()):,} rims; it is transcribed in "
         "`data/raw/mainfort2003_rim_table1_counts.csv` and compared with the file here."), "",
         "## The columns", "",
         "| column(s) | what it is | how the values relate |",
         "|---|---|---|",
         f"| `total` | N, the number of rims analyzed | equals the published N at all {len(m)} sites |",
         ("| `beadrim`, `noderim`, `widenode`, `notchtop`, `notchext` | beaded applique strip; band of nodes; "
         "widely spaced nodes; notched top of the lip; notched exterior | each is the percentage of the N rims "
         "with that attribute. Separate yes-or-no observations: they do not sum to 100 with one another |"),
         ("| `sharpbev`, `modbev`, `slitebev`, `nobev` | degree of interior bevel | parts of one whole: sum "
         f"{int(bev4.min())} to {int(bev4.max())} (100 within rounding at {int(((bev4 >= 98) & (bev4 <= 102)).sum())} "
         f"of {len(m)} sites). The first three together are the paper's \"Beveled\" |"),
         f"| `flatlip`, `roundlip` | lip form | parts of one whole: sum {int(lip.min())} to {int(lip.max())} |",
         (f"| `nf`, `bell` | Neeley's Ferry Plain and Bell Plain | NOT parts of one whole: sum {int(paste.min())} to "
         f"{int(paste.max())}, median {int(paste.median())}. The paste of the rim sherds (author), so percentages "
         "of N; what the remaining rims are is not established. Not in the paper |"),
         ("| `latitude`, `longitude` | rough locations | some are wrong (Davis is at latitude 24.86); the "
         "project's own coordinates are used for matched sites |"), "",
         "## The check against the published counts", "",
         ("For each of the six published attributes, the published count divided by N is compared with the "
         f"file's percentage. They agree within {TOLERANCE:g} points except where listed, so the file's "
         "percentages are percentages of N."), ""]
    L += [f"- {n}" for n in notes] or ["- no disagreements"]
    L += ["", ("The paper's Table 2 (percentages) lists Parkin with N = 111, where Table 1 and the file have 67; "
          "the file's Parkin percentages are those of 67 rims."), "",
          "## The table as counts", "",
          (f"`data/processed/mainfort_rim_counts.csv`, {len(out)} rows. The six attributes are the published "
          "counts as printed; bevel degrees divide the published beveled count in the file's proportions; lip "
          "form and paste are the file's percentages times N and are not in the paper."), "",
          ("| site | N | beaded | noded | wide nodes | notched top | notched exterior | bevel: sharp / moderate / "
          "slight / none | lip: flat / round | Neeley's Ferry / Bell / other |"),
          "|---|---:|---:|---:|---:|---:|---:|---|---|---|"]
    for r in out.itertuples(index=False):
        L.append(f"| {r.site} | {r.N} | {r.beaded} | {r.noded} | {r.wide_nodes} | {r.notched_top} | "
                 f"{r.notched_exterior} | {r.bevel_sharp} / {r.bevel_moderate} / {r.bevel_slight} / {r.bevel_none} | "
                 f"{r.lip_flat} / {r.lip_round} | {r.paste_nf} / {r.paste_bell} / {r.paste_other} |")
    L += paste_check(m)
    L += ["", "## To confirm with Mainfort", "",
          ("1. What the rims in neither paste are. `nf` and `bell` are the paste of the rim sherds (see above), "
          "and the remainder is large where decorated jars are common "
          f"(Fortune {int(m.loc[m['site'] == 'Fortune', 'nf'].iloc[0] + m.loc[m['site'] == 'Fortune', 'bell'].iloc[0])} "
          "percent in total)."),
          ("2. Whether lip form and bevel degree were recorded on all N rims. A few percentages cannot be a whole "
          "number of rims out of N (for example Bishop's lip values), which suggests a slightly smaller base "
          "at some sites."),
          "3. The rows where the file and the published counts disagree, listed above.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:26 + len(notes)]))
    print(f"wrote {OUT_CSV} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
