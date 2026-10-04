"""45_alvey_14c_robustness.py - does adding the Alvey LMV 14C dates move the basin chronology?

Non-destructive robustness check. Analysis 11 builds the St. Francis basin
chronology (SPD span, contact-truncation mass, and the seriation-axis-vs-
calendar-age Spearman) from the Mainfort (2001) compilation alone. The Alvey LMV
14C Database (data/Alvey LMV 14C Database_220714.xlsx) is a wider compilation
that already contains 96 of our 108 Mainfort determinations plus many more.

Here we (1) map Alvey trinomials to our basin site names via shared lab numbers
(authoritative) plus the two hand-confirmed Nodena assignments (Upper Nodena =
3MS4, Middle Nodena = 3MS3), (2) pull the Alvey determinations NOT already in
Mainfort, screened to the Mississippian window (normalized age <= 800 BP), and
(3) re-run analysis 11's SPD and seriation-correlation on the merged set,
reporting BEFORE (Mainfort only) vs AFTER (Mainfort + Alvey net-new).

Read-only on the manuscript; writes output/alvey_14c_robustness.md.

Usage: .venv/bin/python analyses/45_alvey_14c_robustness.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

a11 = importlib.import_module("11_chronology_14c")
import make_figures as mf  # noqa: E402

DATA = ROOT / "data"
ALVEY = DATA / "Alvey LMV 14C Database_220714.xlsx"
OUT = ROOT / "output" / "alvey_14c_robustness.md"

# Basin site name -> Alvey trinomial. Established by shared Mainfort lab numbers,
# with the two Nodena loci confirmed by hand (Upper = 3MS4, Middle = 3MS3).
# Neeley's Ferry (3CS24) is absent from Alvey; Callahan-Thompson (23MI71) adds nothing new.
PROV2TRI = {
    "Parkin": "3CS0029",
    "Hazel": "3PO0006",
    "Kent": "3LE0008",
    "Clay Hill": "3LE0011",
    "Upper Nodena": "3MS0004",
}
MISSISSIPPIAN_MAX_BP = 800.0   # ~ post AD 1250; screens out Archaic/Woodland dates


def _norm_lab(s: str) -> str:
    return re.sub(r"[\s\-_]", "", str(s)).upper()


def alvey_new_basin_dates(mainfort_labs: set) -> pd.DataFrame:
    """Alvey determinations at our basin sites that are NOT in Mainfort and fall
    in the Mississippian window. Returns rows with bp/err/prov like parse_dates()."""
    alv = pd.read_excel(ALVEY)
    alv["lab"] = alv["LABNR"].apply(_norm_lab)
    # Use the normalized age where present, else the measured age: older TX/SMU
    # determinations in Alvey carry a MEAS_AGE but a blank NORM_AGE, and Mainfort's
    # own ages are uncorrected, so this keeps the two compilations comparable.
    alv["age"] = alv["NORM_AGE"].fillna(alv["MEAS_AGE"])
    alv["sig"] = alv["NA_SIGMA"].fillna(alv["MA_SIGMA"])
    rows = []
    for prov, tri in PROV2TRI.items():
        sub = alv[(alv["SITE_IDENTIFIER"] == tri)
                  & (~alv["lab"].isin(mainfort_labs))
                  & (alv["age"].notna())
                  & (alv["age"] <= MISSISSIPPIAN_MAX_BP)]
        for _, r in sub.iterrows():
            err = r["sig"] if pd.notna(r["sig"]) else 40.0
            rows.append({"prov": prov, "bp": float(r["age"]),
                         "err": max(float(err), 20.0), "lab": r["LABNR"]})
    return pd.DataFrame(rows)


def spd_summary(dates: pd.DataFrame, grid_ad, grid_bp) -> dict:
    d = a11.spd(dates, grid_bp)
    return {
        "n": len(dates),
        "median_ad": a11.median_ad(dates, grid_ad),
        "post_1541": float(d[grid_ad >= a11.CONTACT_AD].sum()),
        "post_1600": float(d[grid_ad >= 1600].sum()),
    }


def seriation_corr(df: pd.DataFrame, grid_ad):
    """Replicate analysis 11's CA1-vs-pooled-median-age Spearman on the given
    date table. Returns (anchors DataFrame, rho, p)."""
    counts, _ = mf._load_curated()
    M = counts.to_numpy(float)
    ca1, _, _ = mf.correspondence_axis(M)
    ca = pd.Series(ca1, index=counts.index)
    cur_norm = {mf.norm_name(a): a for a in counts.index}

    def match(prov):
        pn = mf.norm_name(prov)
        hit = cur_norm.get(pn)
        if hit is None:
            for k, v in cur_norm.items():
                if pn and (pn in k or k in pn):
                    return v
        return hit

    # orient CA so increasing = later (same convention as analysis 11)
    tmp = {}
    for prov, g in df.groupby("prov"):
        h = match(prov)
        if h is not None:
            tmp[h] = a11.median_ad(g, grid_ad)
    if len(tmp) >= 3:
        rr, _ = spearmanr([ca[h] for h in tmp], [tmp[h] for h in tmp])
        if np.isfinite(rr) and rr < 0:
            ca = -ca

    anchors = []
    for prov, g in df.groupby("prov"):
        h = match(prov)
        if h is not None:
            anchors.append({"prov": prov, "assem": h, "n": len(g),
                            "median_ad": a11.median_ad(g, grid_ad),
                            "ca1": float(ca[h])})
    A = pd.DataFrame(anchors).drop_duplicates("assem")
    if len(A) >= 3:
        rho, p = spearmanr(A["ca1"], A["median_ad"])
    else:
        rho, p = np.nan, np.nan
    return A, rho, p


def main() -> None:
    grid_ad = np.arange(1000, 1750)
    grid_bp = 1950 - grid_ad

    mf_all = a11.parse_dates()
    mf_all["lab"] = mf_all["Sample ID"].apply(_norm_lab)
    mainfort_labs = set(mf_all["lab"])
    mf_basin = mf_all[mf_all["prov"].isin(a11.BASIN_PROV)].copy()

    new = alvey_new_basin_dates(mainfort_labs)
    merged_basin = pd.concat([mf_basin[["prov", "bp", "err"]],
                              new[["prov", "bp", "err"]]], ignore_index=True)

    before = spd_summary(mf_basin, grid_ad, grid_bp)
    after = spd_summary(merged_basin, grid_ad, grid_bp)
    A_before, rho_b, p_b = seriation_corr(mf_all, grid_ad)
    # merged full table (all Mainfort + Alvey net-new at basin sites) for the correlation
    merged_full = pd.concat([mf_all[["prov", "bp", "err"]],
                             new[["prov", "bp", "err"]]], ignore_index=True)
    A_after, rho_a, p_a = seriation_corr(merged_full, grid_ad)

    L = ["# Alvey 14C robustness check for the St. Francis basin chronology", "",
         "Does adding the Alvey LMV net-new determinations move analysis 11's "
         "basin chronology? BEFORE = Mainfort (2001) only; AFTER = Mainfort + "
         "Alvey determinations at our basin sites that are not already in Mainfort "
         f"and fall in the Mississippian window (normalized age <= {MISSISSIPPIAN_MAX_BP:.0f} BP).", "",
         "## Net-new Alvey dates added (by basin site)", "",
         "| site | Alvey trinomial | new Mississippian dates | ages (BP) |",
         "|---|---|---|---|"]
    for prov, tri in PROV2TRI.items():
        sub = new[new["prov"] == prov]
        ages = ", ".join(f"{int(b)}" for b in sorted(sub["bp"])) if len(sub) else "-"
        L.append(f"| {prov} | {tri} | {len(sub)} | {ages} |")
    L.append(f"| **total** | | **{len(new)}** | |")
    L += ["",
          "Neeley's Ferry (3CS24) is absent from Alvey; Callahan-Thompson (23MI71) "
          "and Parkin/Kent contribute no net-new Mississippian dates.", "",
          "## Basin SPD: before vs after", "",
          "| quantity | Mainfort only | + Alvey net-new |",
          "|---|---|---|",
          f"| basin determinations | {before['n']} | {after['n']} |",
          f"| basin SPD median (AD) | {before['median_ad']:.0f} | {after['median_ad']:.0f} |",
          f"| mass after AD 1541 | {before['post_1541']:.2f} | {after['post_1541']:.2f} |",
          f"| mass after AD 1600 | {before['post_1600']:.2f} | {after['post_1600']:.2f} |", "",
          "## Seriation axis vs calendar age: before vs after", ""]

    def anchor_block(A, rho, p, label):
        out = [f"**{label}** (n = {len(A)} anchors, {int(A['n'].sum())} dates pooled): "
               f"Spearman rho = {rho:+.2f}, p = {p:.3f}", "",
               "| assemblage | n dates | median cal AD | CA1 |", "|---|---|---|---|"]
        for _, r in A.sort_values("median_ad").iterrows():
            out.append(f"| {r['assem']} | {int(r['n'])} | {r['median_ad']:.0f} | {r['ca1']:+.3f} |")
        out.append("")
        return out

    L += anchor_block(A_before, rho_b, p_b, "BEFORE (Mainfort only)")
    L += anchor_block(A_after, rho_a, p_a, "AFTER (+ Alvey net-new)")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
