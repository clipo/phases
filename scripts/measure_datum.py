"""measure_datum.py - the ledger of the NAD27 conversion, and its measured effect.

Phillips, Ford and Griffin's UTMs are NAD27. `dataio/coords.py` converts the
assemblage coordinates that are demonstrably those UTMs read as NAD83, and
`dataio/settlement.py` converts every PFG-sourced row of the settlement table.
This script records which assemblage fell in which class and what the
conversion did to the quantities that depend on position: basin membership,
the k-means partitions at k = 2 to 6, and the plug-in cultural F_ST on each.
It writes `output/findings/datum_conversion.md`. Rule 1: the numbers in
docs/METHODS_DECISIONS.md under the datum ruling come from here.

Usage: .venv/bin/python scripts/measure_datum.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.dataio import coords as C  # noqa: E402
from mls_emergence.dataio.settlement import load_lmv  # noqa: E402
from mls_emergence.dataio.matrix import read_analysis_matrix  # noqa: E402
from mls_emergence.signatures.assortativity import _kmeans_labels  # noqa: E402
from mls_emergence.signatures.variance import cultural_fst  # noqa: E402

XY = ROOT / "data" / "raw" / "mainfort-pfg-cplXY.txt"
OUT = ROOT / "output" / "findings" / "datum_conversion.md"
MEMBERS = ROOT / "data" / "processed" / "basin_members_curated.txt"
K_RANGE = range(2, 7)


def _group(counts, labels):
    return np.array([counts[labels == c].sum(axis=0) for c in np.unique(labels)])


def _to_utm(df):
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", "EPSG:26915", always_xy=True)
    e, n = tr.transform(df["Longitude"].to_numpy(float), df["Latitude"].to_numpy(float))
    return np.column_stack([e, n]) / 1000.0


def main():
    # --- the assemblage ledger -------------------------------------------
    raw = pd.read_csv(XY, sep="\t")
    raw["Assemblages"] = raw["Assemblages"].astype(str).str.strip()
    corr = C.load_corrections()
    for name, row in corr.iterrows():
        sel = raw["Assemblages"] == name
        raw.loc[sel, "Latitude"] = float(row["latitude"])
        raw.loc[sel, "Longitude"] = float(row["longitude"])
    table = C.pfg_datum_table(raw)
    mx = pd.read_excel(ROOT / "data" / "raw" / "mainfort-spatial-data.xlsx", sheet_name="mainfortXY")
    mx["Site_Name"] = mx["Site_Name"].astype(str).str.strip()
    mx = mx.drop_duplicates("Site_Name").set_index("Site_Name")
    rr = raw.set_index("Assemblages")
    own = {}
    for name in table.index:
        if name in mx.index:
            own[name] = C._planar_m(rr.loc[name, "Latitude"], rr.loc[name, "Longitude"],
                                    float(mx.loc[name, "Latitude"]), float(mx.loc[name, "Longitude"]))
    table["d_mainfort_own_m"] = pd.Series(own)
    members = [l.strip() for l in MEMBERS.read_text().splitlines() if l.strip()]
    table["basin_member"] = table.index.isin(members)

    converted = C.apply_pfg_datum(raw, verbose=False)

    # --- effect on the basin set ------------------------------------------
    mat = read_analysis_matrix()
    name_col = next(c for c in mat.columns if str(c).lower() in ("assemblage", "assemblages", "name"))
    counts = mat.set_index(name_col).reindex(members).select_dtypes("number")
    if counts.isna().any().any():
        raise RuntimeError("basin member absent from the analysis matrix")
    before = raw.set_index("Assemblages").reindex(members)
    after = converted.set_index("Assemblages").reindex(members)
    xb, xa = _to_utm(before), _to_utm(after)
    moved = np.hypot(*(xa - xb).T * 1000)
    rows = []
    for k in K_RANGE:
        lb = _kmeans_labels(xb, k, seed=7)
        la = _kmeans_labels(xa, k, seed=7)
        # label permutation-invariant comparison
        same = len(set(zip(lb, la))) == k
        fb = cultural_fst(_group(counts.to_numpy(float), lb))
        fa = cultural_fst(_group(counts.to_numpy(float), la))
        rows.append((k, same, fb, fa))

    # --- settlement table ------------------------------------------------
    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    n_conv = int((lmv["Datum"] != "as recorded").sum())

    L = ["# The NAD27 conversion: ledger and measured effect", "",
         f"Produced by `scripts/measure_datum.py`. Tolerance {C.PFG_DATUM_TOL_M:.0f} m "
         f"(PFG UTMs are rounded to 10 m). Basin set of {len(members)} assemblages "
         f"from `basin_members_curated.txt`.", "",
         "## Assemblage coordinates (`mainfort-pfg-cplXY.txt`, after the published corrections)", "",
         "| class | all assemblages | basin members | meaning |", "|---|---|---|---|"]
    meaning = {"converted": "file coordinate is PFG's UTM read as NAD83; moved to the NAD27 reading",
               "near": "within 100 m of that reading but not within tolerance; derivation not demonstrated; left as recorded",
               "other": "more than 100 m from PFG's UTM either way; left as recorded",
               "no_pfg_utm": "no PFG UTM to compare against; left as recorded"}
    for kl in ["converted", "near", "other", "no_pfg_utm"]:
        s = table[table["klass"] == kl]
        L.append(f"| {kl} | {len(s)} | {int(s['basin_member'].sum())} | {meaning[kl]} |")
    L += ["", "### Basin members, row by row", "",
          "| assemblage | site | class | m from NAD83 reading | m from NAD27 reading | m from Mainfort's own XY |",
          "|---|---|---|---|---|---|"]
    for name, r in table[table["basin_member"]].sort_values(["klass", "d_nad83_m"]).iterrows():
        f = lambda v: "-" if pd.isna(v) else f"{v:.0f}"
        L.append(f"| {name} | {r['site_number'] or '-'} | {r['klass']} | {f(r['d_nad83_m'])} | "
                 f"{f(r['d_nad27_m'])} | {f(r['d_mainfort_own_m'])} |")
    near_own = table[(table["klass"] == "near") & (table["d_mainfort_own_m"] <= C.PFG_DATUM_TOL_M)]
    L += ["", f"Of the {int((table['klass']=='near').sum())} near cases, "
          f"{len(near_own)} coincide with Mainfort's own coordinate for the site "
          f"(`mainfortXY` sheet) within tolerance; those are his points, whose datum is "
          f"not recorded, and they stand as recorded under the 2026-09-19 ruling on "
          f"non-PFG coordinates.", "",
          "## Settlement table (`LMVData.xlsx`)", "",
          f"{n_conv} rows whose Source begins with \"PFG\", or that the settlement correction "
          f"table replaced with PFG UTMs, converted from EPSG:267zz to EPSG:269zz; all others "
          f"as recorded. For the basin's assemblage sites the compilation's UTMs are identical "
          f"to PFG's site table digit for digit, which is why the Source label is trusted.", "",
          "## Effect on the basin set", "",
          f"Assemblages moved: {int((moved > 1).sum())} of {len(members)}, by "
          f"{moved[moved > 1].min():.0f} to {moved[moved > 1].max():.0f} m.", "",
          "| k | same partition? | plug-in F_ST before | after |", "|---|---|---|---|"]
    for k, same, fb, fa in rows:
        L.append(f"| {k} | {'yes' if same else 'NO'} | {fb:.6f} | {fa:.6f} |")
    L += ["", "The partition is the same at every k listed above and the plug-in F_ST is "
          "unchanged to the printed precision" if all(r[1] and abs(r[2]-r[3]) < 5e-7 for r in rows)
          else "**The partition or the F_ST changed at some k; every position-dependent "
          "result must be rerun and re-read.**", ""]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
