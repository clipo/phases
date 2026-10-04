"""98_small_collections.py - which existing collections would a lower sherd minimum admit?

The primary analysis uses 28 assemblages with at least 75 decorated sherds. The
archaeological record will not supply new assemblages, but the
Phillips-Ford-Griffin table holds many collections from the same ground that
were never candidates, because the analysis matrix was built from a list of
named assemblages and then cut at 75. This script lists them: every PFG
collection in Phillips's Kent, Parkin and Walls areas that is not already
analyzed, with its decorated-sherd count in the ten analysis classes, its
location, and its phase.

It changes nothing in the analysis. It is the inventory a decision about the
minimum needs, and `97_design_to_separate.py` reads its table to ask what
adding these collections would buy.

HOW A COLLECTION IS COUNTED.
  classes   PFGData.xlsx collapsed into the ten analysis classes by the
            mapping `scripts/build_analysis_matrix.py` verifies on every run.
  one row   PFGData lists some sites more than once, as overlapping groups of
            collection units ("13-N-4/B,C,D,E" and "13-N-4/C"). Summing them
            would count sherds twice, so each site number takes its single
            largest row.
  location  the PFG site number's coordinate in the settlement compilation
            (data/LMVData.xlsx, with its recorded corrections).
  phase     the Phillips (1970) area the site lies in, or the nearest one
            inside the source figure's frame, exactly as the analyzed
            assemblages are assigned; a site in or nearest the Nodena area, or
            outside the frame, is not a candidate.
  analyzed  an analyzed assemblage is recognized by its PFG site number and
            left out. A different site number within 1 km of one is kept and
            flagged.

Outputs: output/findings/small_collections.md, output/findings/small_collections.csv
Usage: PYTHONPATH=src .venv/bin/python analyses/98_small_collections.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "small_collections.md"
OUT_CSV = ROOT / "output" / "findings" / "small_collections.csv"
PHASES = ("Kent", "Parkin", "Walls")
SAME_SITE_KM = 1.0
CURRENT_MINIMUM = 75
MINIMA = (75, 50, 25, 10)
PLAIN = ("Neeley's Ferry Plain", "Bell Plain")
# The PFG site number of each analyzed assemblage. By number, not by distance:
# several analyzed assemblages carry coordinates corrected by the author, which
# sit kilometers from the settlement compilation's point for the same site, and
# matching on distance alone listed Belle Meade and Lake Cormorant as new
# collections. Names are PFGData's; Holden Lake is Lipo's (2001) collection and
# has no PFG number.
ANALYZED_NUMBER = {
    "Barton_Ranch": "11-O-10", "Beck": "13-O-7", "Belle_Meade": "13-O-5", "Big_Eddy": "12-N-4",
    "Castile_Landing": "13-N-21", "Clay_Hill": "13-N-7", "Commerce": "13-O-11",
    "Cramor_Place": "12-O-5", "Cummins": "11-O-4", "Davis": "13-N-5", "Fortune": "11-N-15",
    "Grant": "13-N-11", "Holden_Lake": None, "Hollywood": "13-O-10", "Irby": "13-P-10",
    "Kent_Place": "13-N-4", "Lake_Cormorant": "13-P-8", "Mound_Place": "12-P-1",
    "Neeleys_Ferry": "11-N-4", "Nickel": "13-N-15", "Parkin": "11-N-1", "Rose_Mound": "12-N-3",
    "Starkley": "13-N-16", "Turnbow": "11-N-12", "Vernon_Paul": "11-N-9", "Walls": "13-P-1",
    "Williamson": "11-N-13", "Woodlyn": "13-P-11",
}


def base_number(site_id: str) -> str | None:
    """'13-N-4/B,C,D,E' -> '13-N-4'; None when the id is not a PFG site number."""
    m = re.match(r"^\s*(\d+)\s*-\s*([A-Z])\s*-\s*(\d+)", str(site_id))
    return f"{int(m.group(1))}-{m.group(2)}-{int(m.group(3))}" if m else None


def main() -> int:
    mf = importlib.import_module("make_figures")
    mm = importlib.import_module("make_map")
    ph = importlib.import_module("36_canonical_phase_map")
    bam = importlib.import_module("build_analysis_matrix")
    from mls_emergence.dataio.settlement import load_lmv
    from pyproj import Transformer

    counts, coords = mf._load_curated()
    names28 = [str(i) for i in counts.index]
    classes = list(counts.columns)
    xy28 = coords.to_numpy(float)
    lab28, _ = ph.assign_primary_phases(names28, xy28)
    lab28 = np.asarray(lab28)
    tr = Transformer.from_crs("EPSG:4326", mm.UTM15N, always_xy=True)
    E28, N28 = (np.asarray(v, float) for v in tr.transform(xy28[:, 1], xy28[:, 0]))

    # ---- PFG collections, one row per site number
    P = bam.pfg_table(classes)
    P["number"] = P["id"].map(base_number)
    P = P[P["number"].notna()].copy()
    P["decorated"] = P[classes].sum(axis=1)
    raw = pd.read_excel(bam.PFG_XLSX, sheet_name="SherdData")
    raw.columns = [str(c).strip() for c in raw.columns]
    plain_cols = [c for c in raw.columns if c in PLAIN]
    if len(plain_cols) != len(PLAIN):
        raise RuntimeError(f"plainware columns not found: have {plain_cols}")
    P["plain"] = raw.loc[P.index, plain_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    n_rows = len(P)
    P = P.sort_values("decorated", ascending=False).drop_duplicates("number")
    n_sites = len(P)

    # ---- locations
    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    lmv = lmv[pd.to_numeric(lmv["Zone"], errors="coerce") == 15].copy()
    lmv["number"] = lmv["Number"].map(base_number)
    lmv["E"] = pd.to_numeric(lmv["Easting"], errors="coerce")
    lmv["N"] = pd.to_numeric(lmv["Northing"], errors="coerce")
    loc = lmv[lmv["number"].notna() & lmv["E"].notna() & lmv["N"].notna()].drop_duplicates("number")
    P = P.merge(loc[["number", "E", "N"]], on="number", how="left")
    n_no_coord = int(P["E"].isna().sum())
    P = P[P["E"].notna()].copy()
    back = Transformer.from_crs(mm.UTM15N, "EPSG:4326", always_xy=True)
    lon, lat = back.transform(P["E"].to_numpy(float), P["N"].to_numpy(float))
    P["lat"], P["lon"] = lat, lon

    # ---- phase, and whether it is already an analyzed assemblage
    lab, derived = ph.assign_phases_phillips(list(P["number"]), P[["lat", "lon"]].to_numpy(float))
    P["phase"], P["by_nearest"] = lab, derived
    d28 = np.hypot(P["E"].to_numpy(float)[:, None] - E28[None, :],
                   P["N"].to_numpy(float)[:, None] - N28[None, :]) / 1000.0
    P["km_to_analyzed"] = d28.min(1)
    missing = [n for n in names28 if n not in ANALYZED_NUMBER]
    if missing:
        raise RuntimeError(f"no PFG site number recorded for analyzed assemblages {missing}; "
                           "add them to ANALYZED_NUMBER before listing candidates")
    by_number = {v: k for k, v in ANALYZED_NUMBER.items() if v is not None and k in names28}
    P["analyzed_as"] = P["number"].map(by_number).fillna("")
    # A different site number within SAME_SITE_KM of an analyzed assemblage is
    # kept as a candidate but flagged: it may be a second collection at one place.
    P["flag_close"] = (P["analyzed_as"] == "") & (P["km_to_analyzed"] <= SAME_SITE_KM)
    in_scheme = P["phase"].isin(PHASES)
    n_out_of_scheme = int((~in_scheme).sum())
    cand = P[in_scheme & (P["analyzed_as"] == "") & (P["decorated"] > 0)].copy()
    # How close to a phase line: distance to the nearest analyzed assemblage of another phase.
    Ec, Nc = cand["E"].to_numpy(float), cand["N"].to_numpy(float)
    other = np.array([np.hypot(Ec[i] - E28[lab28 != p], Nc[i] - N28[lab28 != p]).min() / 1000.0
                      for i, p in enumerate(cand["phase"])]) if len(cand) else np.array([])
    cand["km_to_other_phase"] = other
    cand = cand.sort_values("decorated", ascending=False)
    cols = ["number", "name", "phase", "by_nearest", "flag_close", "decorated", "plain", "lat", "lon",
            "km_to_analyzed", "km_to_other_phase"] + classes
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    cand[cols].to_csv(OUT_CSV, index=False)

    cur = {p: int((lab28 == p).sum()) for p in PHASES}
    L = ["# Which existing collections would a lower sherd minimum admit?", "",
         "Produced by `analyses/98_small_collections.py`. PFGData.xlsx holds "
         f"{n_rows} collection rows with a PFG site number, {n_sites} distinct site numbers "
         f"(each takes its largest row, since some rows overlap); {n_no_coord} have no coordinate in the "
         f"settlement compilation; {n_out_of_scheme} of the rest lie in or nearest the Nodena area or "
         "outside the frame of Phillips's figure. Of the sites in the Kent, Parkin and Walls areas, "
         f"{int((in_scheme & (P['analyzed_as'] != '')).sum())} are analyzed assemblages, recognized by site "
         f"number; {int((in_scheme & (P['analyzed_as'] == '') & (P['decorated'] == 0)).sum())} "
         f"hold no sherd of the ten classes. **{len(cand)} collections remain as candidates.**", "",
         f"The analysis now holds {sum(cur.values())} assemblages (" + ", ".join(f"{p} {n}" for p, n in cur.items())
         + f") at a minimum of {CURRENT_MINIMUM} decorated sherds.", "",
         "| minimum decorated sherds | collections added | Kent / Parkin / Walls added | assemblages in all | "
         "added within 10 km of an analyzed assemblage of another phase | median decorated sherds of those added |",
         "|---|---|---|---|---|---|"]
    for mn in MINIMA:
        sub = cand[cand["decorated"] >= mn]
        L.append(f"| {mn} | {len(sub)} | " + " / ".join(str(int((sub['phase'] == p).sum())) for p in PHASES)
                 + f" | {sum(cur.values()) + len(sub)} | {int((sub['km_to_other_phase'] <= 10).sum())} | "
                 + (f"{sub['decorated'].median():.0f}" if len(sub) else "-") + " |")
    L += ["",
          f"A candidate with {CURRENT_MINIMUM} or more decorated sherds was left out because the analysis matrix "
          "was built from a list of named assemblages, not because of its count.", "",
          "## The candidates", "",
          "| site | name | phase | decorated sherds | plainware sherds | km to nearest analyzed assemblage | "
          "km to nearest analyzed assemblage of another phase |",
          "|---|---|---|---|---|---|---|"]
    for r in cand.itertuples(index=False):
        L.append(f"| {r.number} | {r.name}{' (within 1 km of an analyzed assemblage)' if r.flag_close else ''} | "
                 f"{r.phase}{' (nearest area)' if r.by_nearest else ''} | "
                 f"{int(r.decorated)} | {int(r.plain)} | {r.km_to_analyzed:.1f} | {r.km_to_other_phase:.1f} |")
    L += ["", "## Reading", "",
          "Counts are of the ten analysis classes only. Plainware counts (Neeley's Ferry Plain and Bell Plain) "
          "are listed because most collections hold far more plain than decorated sherds; they are paste and "
          "temper categories, not decoration, and are not part of the analysis.", "",
          "Phase is by position in Phillips's areas, as for the analyzed assemblages. The coordinates are the "
          "settlement compilation's, which carry the PFG site table's precision and the corrections recorded "
          "there; none of these sites has been checked individually as the analyzed ones were.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:14]))
    print(f"wrote {OUT_MD} and {OUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
