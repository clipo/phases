"""46_radiocarbon_table.py - the manuscript's radiocarbon determination list.

Emits the list of radiocarbon determinations that underlie the St. Francis basin
chronology (analysis 11 + the analysis 45 Alvey robustness augmentation), so the
exact dates used are documented for the Supplemental. Each row is one
determination: site, trinomial, lab number, conventional 14C age +/- error,
calibrated median (AD, IntCal20), source, and which component of the analysis it
feeds.

Writes:
  - docs/manuscript/radiocarbon_dates_used.csv  (supplemental data file)
  - output/radiocarbon_dates_used.md            (formatted table for review)

Read-only on the manuscript. Usage: .venv/bin/python analyses/46_radiocarbon_table.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

a11 = importlib.import_module("11_chronology_14c")
a45 = importlib.import_module("45_alvey_14c_robustness")

CSV = ROOT / "docs" / "manuscript" / "radiocarbon_dates_used.csv"
OUT_MD = ROOT / "output" / "radiocarbon_dates_used.md"

# Basin site name -> trinomial (shared-lab crosswalk; Nodena loci hand-confirmed).
PROV2TRI = {
    "Parkin": "3CS29", "Hazel": "3PO6", "Kent": "3LE8", "Clay Hill": "3LE11",
    "Upper Nodena": "3MS4", "Neeley's Ferry": "3CS24", "Callahan-Thompson": "23MI71",
}
GRID_AD = np.arange(1000, 1750)


def _cal_median(bp, err):
    d = pd.DataFrame([{"bp": bp, "err": err}])
    return a11.median_ad(d, GRID_AD)


def main() -> None:
    rows = []

    # --- Mainfort (2001) basin determinations: the primary compilation ---
    mf = a11.parse_dates()
    mf["lab"] = mf["Sample ID"].astype(str).str.strip()
    mf_basin = mf[mf["prov"].isin(a11.BASIN_PROV)]
    for _, r in mf_basin.iterrows():
        rows.append({
            "site": r["prov"],
            "trinomial": PROV2TRI.get(r["prov"], ""),
            "lab_number": r["lab"],
            "c14_age_bp": int(r["bp"]),
            "error": int(r["err"]),
            "cal_median_ad": round(_cal_median(r["bp"], r["err"])),
            "source": "Mainfort 2001",
            "used_in": "basin SPD; seriation correlation (curated sites)",
        })

    # --- Alvey (pers. comm.) net-new Mississippian determinations: robustness ---
    mlabs = set(mf["Sample ID"].apply(a45._norm_lab))
    new = a45.alvey_new_basin_dates(mlabs)
    for _, r in new.iterrows():
        rows.append({
            "site": r["prov"],
            "trinomial": PROV2TRI.get(r["prov"], ""),
            "lab_number": str(r["lab"]),
            "c14_age_bp": int(r["bp"]),
            "error": int(r["err"]),
            "cal_median_ad": round(_cal_median(r["bp"], r["err"])),
            "source": "Alvey (personal communication, 2022)",
            "used_in": "robustness augmentation (Supplemental)",
        })

    df = pd.DataFrame(rows)
    # one row per determination: drop duplicate lab numbers within a site (the
    # Mainfort compilation lists TX-848 twice)
    df = df.drop_duplicates(subset=["site", "lab_number"], keep="first")
    df = df.sort_values(["site", "c14_age_bp"]).reset_index(drop=True)
    CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(CSV, index=False)

    n_mf = int((df["source"] == "Mainfort 2001").sum())
    n_al = len(df) - n_mf
    L = ["# Radiocarbon determinations used in the St. Francis basin chronology", "",
         f"{len(df)} determinations: {n_mf} from Mainfort (2001) (primary basin SPD and "
         f"seriation-vs-calendar test) and {n_al} net-new Mississippian dates from the "
         f"Alvey lower-valley compilation (personal communication, 2022; robustness check "
         f"only). Calibrated medians are point calibrations against IntCal20.", "",
         "Lab-number prefixes denote the radiocarbon laboratory, not a location: "
         "Beta = Beta Analytic; TX = University of Texas at Austin; M = University "
         "of Michigan; SMU = Southern Methodist University; UGa = University of "
         "Georgia; DRI = Desert Research Institute. Sites are keyed by Smithsonian "
         "trinomial (state-county-site).", "",
         "| Site | Trinomial | Lab number | 14C age BP | ± | Cal median AD | Source | Used in |",
         "|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        L.append(f"| {r['site']} | {r['trinomial']} | {r['lab_number']} | {r['c14_age_bp']} "
                 f"| {r['error']} | {r['cal_median_ad']} | {r['source']} | {r['used_in']} |")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    print(f"wrote {CSV} ({len(df)} determinations: {n_mf} Mainfort, {n_al} Alvey net-new)")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
