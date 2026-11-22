#!/usr/bin/env python3
"""The painted class by type: Old Town Red, Nodena Red and White, and Avenue Polychrome in the 1940 collections.

The analysis matrix pools three painted types into one class where a row is built from PFGData
(scripts/build_analysis_matrix.py, PAINTED). Rows taken from Lipo's (2001) table follow his Painted category, which
his Table 7.4, note a, gives as five types: those three with Carson Red-and-Buff and Hollywood White. This script
states what the class consists of at the 28 basin assemblages, from the type-level counts of data/raw/PFGData.xlsx,
and reports where the two further types are counted in the matrix and where they are not (`extra_account`).

Scope. This script splits only the 1940 Phillips-Ford-Griffin collections by type. The data files used here record
Lipo's 1996 to 1997 collections as the pooled class. Lipo (2001: Table 6.1) publishes those collections by type, and
that table is not read here.

Rows. The 1940 part of each assemblage is taken from the same PFGData rows as
`composition_by_collection.decompose` uses: for an assemblage in the cpl-pfg tab, the rows of its grid id (with the
ruled override for Kent Place, whose collection C appears in two rows); for any other assemblage, the PFGData rows of
its normalized site name, which is how the matrix was built. The selection is checked on every run: collapsed into
the ten analysis classes, the selected rows must equal the 1940 part that `decompose` returns, class by class, and the
three types must sum to that part's painted count. The script raises if either fails.

Descriptive throughout (rule 18): counts, shares, and a rank correlation of each type's share of the 1940 decorated
sherds with latitude (positive: a larger share to the north). No interval is attached.

Usage: .venv/bin/python analyses/per_type/painted_by_type.py
Output: output/findings/painted_by_type.md
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "per_type", ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from mls_emergence.dataio.basin import load_curated  # noqa: E402

OUT = ROOT / "output" / "findings" / "painted_by_type.md"
PFG_XLSX = ROOT / "data" / "raw" / "PFGData.xlsx"
TYPES = ("Old Town Red", "Nodena Red and White", "Avenue Polychrome")
# Two further painted types of the 1940 classification. Lipo (2001: Table 7.4, note a) includes them in his Painted
# category; scripts/build_analysis_matrix.py (PAINTED) does not. They are counted here so the difference is on record.
EXTRA = ("Carson Red-on-Buff", "Hollywood White Slipped")
# PFGData site names that differ from the assemblage name (scripts/build_analysis_matrix.py, PFG_NAME, and
# scripts/reconcile_ceramic_matrix.py, ALIAS). Each is confirmed below by the ten-class equality, not by resemblance.
NAME_ALIAS = {"grant": "grantplace", "bartonranch": "barton"}


def raw_rows() -> pd.DataFrame:
    """PFGData with the three painted types as integers, the site number, its base, and the normalized name."""
    cbc = importlib.import_module("composition_by_collection")
    d = pd.read_excel(PFG_XLSX, sheet_name="SherdData")
    d.columns = [str(c).strip() for c in d.columns]
    out = pd.DataFrame({"id": d["Site Number"].astype(str).str.strip(), "name": d["Site Name"].astype(str).str.strip()})
    out["base"] = out["id"].str.split("/").str[0]
    out["key"] = out["name"].map(cbc.norm)
    for t in TYPES + EXTRA:
        v = pd.to_numeric(d[t], errors="coerce").fillna(0)
        if (v < 0).any() or not np.equal(v, np.rint(v)).all():
            raise ValueError(f"{t}: PFGData holds a negative or non-integer count")
        out[t] = v.astype(int)
    return out


def rows_for(name: str, raw: pd.DataFrame, ids_all: list[str], g2n: dict) -> pd.DataFrame:
    """The PFGData rows that make up the 1940 part of one assemblage (see the module docstring)."""
    cbc = importlib.import_module("composition_by_collection")
    if cbc.tab_ids(name, ids_all, g2n):
        gid = cbc.one_tab_id(name, ids_all, g2n)
        if gid in cbc.PFG_ROW_OVERRIDE:
            return raw[raw["id"] == cbc.PFG_ROW_OVERRIDE[gid]]
        return raw[raw["base"] == gid]
    k = cbc.norm(name)
    return raw[raw["key"] == NAME_ALIAS.get(k, k)]


def by_type() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per assemblage: the three types in the 1940 part, that part's decorated total, the pooled painted count of the
    matrix row and of Lipo's addition, and latitude. Second value: the provenance table of `decompose`."""
    cbc = importlib.import_module("composition_by_collection")
    bam = importlib.import_module("build_analysis_matrix")
    counts, coords = load_curated()
    TY = list(counts.columns)
    _c, pfg_part, addition, cls, prov = cbc.decompose()
    _tab, ids_all, g2n = cbc.load_sources()
    raw = raw_rows()
    P = bam.pfg_table(TY).reset_index(drop=True)          # the same file collapsed to the ten classes, row for row
    if len(P) != len(raw) or not (P["id"].to_numpy() == raw["id"].to_numpy()).all():
        raise ValueError("the type-level and class-level readings of PFGData are not row-aligned")
    out = []
    for name in counts.index:
        rows = rows_for(name, raw, ids_all, g2n)
        ten = P.loc[rows.index, TY].sum().to_numpy(float)
        if not np.array_equal(ten, pfg_part.loc[name].to_numpy(float)):
            raise ValueError(f"{name}: the selected PFGData rows {rows['id'].tolist()} do not reproduce the 1940 part of the matrix")
        t = rows[list(TYPES)].sum()
        if int(t.sum()) != int(pfg_part.loc[name, "Painted"]):
            raise ValueError(f"{name}: the three types sum to {int(t.sum())}, the 1940 painted count is {int(pfg_part.loc[name, 'Painted'])}")
        out.append({"assemblage": name, "latitude": float(coords.loc[name, "Latitude"]), "collection": cls[name],
                    "PFG rows": "; ".join(rows["id"]) or "none",
                    **{ty: int(t[ty]) for ty in TYPES}, **{ty: int(rows[ty].sum()) for ty in EXTRA},
                    "decorated, 1940": int(pfg_part.loc[name].sum()),
                    "painted, matrix": int(counts.loc[name, "Painted"]), "painted, Lipo addition": int(addition.loc[name, "Painted"])})
    return pd.DataFrame(out).set_index("assemblage"), prov


def summary(t: pd.DataFrame) -> dict:
    """Totals, shares, sites present, and each type's rank correlation with latitude, over assemblages with a 1940 part."""
    has = t[t["decorated, 1940"] > 0]
    tot = {ty: int(has[ty].sum()) for ty in TYPES}
    painted = sum(tot.values())
    rho = {}
    for ty in TYPES:
        share = has[ty] / has["decorated, 1940"]
        rho[ty] = float(spearmanr(has["latitude"], share)[0]) if (share > 0).sum() >= 3 else float("nan")
    return {"n_sites": len(has), "totals": tot, "painted_1940": painted, "share": {ty: tot[ty] / painted for ty in TYPES},
            "sites_present": {ty: int((has[ty] > 0).sum()) for ty in TYPES},
            "sites_named": {ty: list(has.index[has[ty] > 0]) for ty in TYPES}, "rho_latitude": rho,
            "addition_painted": int(t["painted, Lipo addition"].sum()), "matrix_painted": int(t["painted, matrix"].sum()),
            **extra_account(t)}


def extra_account(t: pd.DataFrame) -> dict:
    """Where the two further painted types stand in the matrix. An assemblage whose painted "addition" equals its
    1940 count of the two types, with that count above zero, has them counted in its Painted class (the row comes
    from Lipo's table). An assemblage with such sherds and no painted addition has them left out."""
    ex = t[list(EXTRA)].sum(axis=1)
    counted = [a for a in t.index if ex[a] > 0 and t.loc[a, "painted, Lipo addition"] == ex[a]]
    left_out = [a for a in t.index if ex[a] > 0 and t.loc[a, "painted, Lipo addition"] == 0]
    other = [a for a in t.index if ex[a] > 0 and a not in counted and a not in left_out]
    if other:
        raise ValueError(f"{other}: the two further painted types are neither wholly counted nor wholly left out; read these rows by hand")
    return {"extra_totals": {ty: int(t[ty].sum()) for ty in EXTRA}, "extra_counted": {a: int(ex[a]) for a in counted},
            "extra_left_out": {a: int(ex[a]) for a in left_out},
            "addition_painted_own": int(t["painted, Lipo addition"].sum()) - int(sum(ex[a] for a in counted))}


def main() -> int:
    t, _prov = by_type()
    s = summary(t)
    L = ["# The painted class by type, 1940 collections", "",
         ("Produced by `analyses/per_type/painted_by_type.py` from `data/raw/PFGData.xlsx` (sheet SherdData) and the analysis matrix. "
          "The analysis class \"Painted\" pools Old Town Red, Nodena Red and White, and Avenue Polychrome in rows built from PFGData, and follows "
          "Lipo's five-type category in rows taken from his table (next paragraph but one). Only the 1940 Phillips-Ford-Griffin collections "
          "are split by type here; the data files record Lipo's 1996 to 1997 collections as the pooled class. "
          "Descriptive counts and shares; no interval is attached."), "",
         f"Assemblages with a 1940 part: {s['n_sites']} of {len(t)}. Painted sherds in those parts: {s['painted_1940']}. "
         f"Painted sherds in the matrix above those three-type counts, not split by type: {s['addition_painted']}. The two sum to the matrix's "
         f"{s['matrix_painted']} painted sherds.", "",
         ("**Two further painted types.** Lipo (2001: Table 7.4, note a) states that his Painted category includes Nodena "
          "Red-and-White, Carson Red-and-Buff, Old Town Red, Hollywood White, and Avenue Polychrome. The 1940 rows of the 28 "
          f"assemblages hold {s['extra_totals'][EXTRA[0]]} {EXTRA[0]} and {s['extra_totals'][EXTRA[1]]} {EXTRA[1]} sherds. "
          "Counted in the matrix's Painted class, because the row comes from his table: "
          + ("; ".join(f"{a} {n}" for a, n in s["extra_counted"].items()) or "none") + ". Not counted, because the row was built from "
          "PFGData by name with the three types: " + ("; ".join(f"{a} {n}" for a, n in s["extra_left_out"].items()) or "none") + ". "
          f"Of the {s['addition_painted']} painted sherds above the three-type 1940 counts, {s['addition_painted'] - s['addition_painted_own']} "
          f"are therefore 1940 sherds of the two further types and {s['addition_painted_own']} are in Lipo's own collections, which he "
          "recorded under the five-type category. The table below covers the three types only."), "",
         "| type | sherds, 1940 | share of the painted class, 1940 | sites present | rank correlation of the type's share of the 1940 decorated sherds with latitude (positive: larger to the north) |",
         "|---|---|---|---|---|"]
    for ty in TYPES:
        r = s["rho_latitude"][ty]
        rtxt = f"{r:+.2f}" if np.isfinite(r) else f"not computed (present at {s['sites_present'][ty]} site{'s' if s['sites_present'][ty] != 1 else ''})"
        L.append(f"| {ty} | {s['totals'][ty]} | {100 * s['share'][ty]:.1f} percent | {s['sites_present'][ty]} of {s['n_sites']} | {rtxt} |")
    nrw = s["sites_named"]["Nodena Red and White"]
    L += ["", f"Nodena Red and White occurs at: {', '.join(nrw) if nrw else 'no assemblage'}.", "",
          "## By assemblage, north to south", "",
          "| assemblage | latitude | collection | PFGData rows | Old Town Red | Nodena Red and White | Avenue Polychrome | Carson Red-on-Buff | Hollywood White Slipped | decorated sherds, 1940 (ten classes) | painted in the matrix above the three types |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, r in t.sort_values("latitude", ascending=False).iterrows():
        L.append(f"| {name} | {r['latitude']:.3f} | {r['collection']} | {r['PFG rows']} | {r[TYPES[0]]} | {r[TYPES[1]]} | {r[TYPES[2]]} | {r[EXTRA[0]]} | {r[EXTRA[1]]} | "
                 f"{r['decorated, 1940']} | {r['painted, Lipo addition']} |")
    L += ["", ("Checks made on this run: for every assemblage, the selected PFGData rows collapsed into the ten analysis classes equal the "
               "1940 part returned by `composition_by_collection.decompose`, and the three types sum to that part's painted count."), ""]
    OUT.write_text("\n".join(L))
    print(f"wrote {OUT}")
    print({k: s[k] for k in ("n_sites", "totals", "painted_1940", "sites_present", "rho_latitude", "addition_painted", "matrix_painted")})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
