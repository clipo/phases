"""Rim attribute counts for the basin assemblages (rim plan R1 plumbing, 2026-10-19).

Source: `data/processed/mainfort_rim_counts.csv`, built by
`scripts/build_mainfort_rim_counts.py` from Mainfort (2003b, Table 1), the
published counts of six rim attributes on 7,195 rims at 55 sites, and the
file Mainfort sent Carl, which adds the degree of interior bevel, lip form
and paste as percentages (converted to counts there, with the published N).
That script documents what each column is and its reading of the paste
columns (unconfirmed; the paste percentages do not sum to 100 in the file,
so `paste_other` is a remainder).

This module does three things and states no result: it matches the rim
sites to the 28 basin assemblages by normalized name with the one alias
the two sources need (Kent for Kent_Place), it re-checks the published
counts against the raw transcription at the matched sites, and it gives the
attribute counts, rim totals and pooled frequencies the R1 readings use.
Which attributes were tabulated as types in the 1951 typology (paste: Bell
Plain, Neeley's Ferry Plain) and which were not (bevel, notching, lip form)
is recorded in the constants below. Not tabulated as a type is not the same
as independent of the phase definitions: Phillips, Ford and Griffin knew
interior beveling and the rim fillet (Mainfort 2003b:34 cites 1951:116 and
123), and whether any attribute entered Phillips's 1970 phase definitions
is UNVERIFIED (Codex RIM-3, 2026-10-19). The reading scripts keep paste
separate.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RIM_CSV = ROOT / "data" / "processed" / "mainfort_rim_counts.csv"
PUBLISHED_CSV = ROOT / "data" / "raw" / "mainfort2003_rim_table1_counts.csv"

PUBLISHED = ("beaded", "noded", "wide_nodes", "notched_top", "notched_exterior", "beveled")
BEVEL_DEGREES = ("bevel_sharp", "bevel_moderate", "bevel_slight", "bevel_none")
LIP = ("lip_flat", "lip_round")
PASTE = ("paste_nf", "paste_bell", "paste_other")
# Observations not tabulated as types in the 1951 typology (their relation to the phase definitions is unverified).
FORM_ATTRIBUTES = PUBLISHED + BEVEL_DEGREES + LIP
# Basin assemblage name -> rim-table site name, where normalization alone does not match.
ALIASES = {"Kent_Place": "Kent"}


def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def load_rim_counts() -> pd.DataFrame:
    """The 55-site count table, indexed by the rim table's site name, with its sums checked."""
    r = pd.read_csv(RIM_CSV).set_index("site")
    for cols in (BEVEL_DEGREES, LIP, PASTE):
        bad = r.index[(r[list(cols)].sum(axis=1) != r["N"])]
        if len(bad):
            raise ValueError(f"{cols} do not sum to N at {bad.tolist()}")
    if ((r[["bevel_sharp", "bevel_moderate", "bevel_slight"]].sum(axis=1)) != r["beveled"]).any():
        raise ValueError("bevel degrees do not sum to the published beveled count")
    return r


def match_basin(basin_names) -> tuple[dict, list]:
    """Map each basin assemblage to its rim-table site (exactly one match, by normalized name or alias); return (matches, unmatched)."""
    r = load_rim_counts(); by_norm: dict[str, list[str]] = {}
    for site in r.index:
        by_norm.setdefault(norm(site), []).append(site)
    matches, unmatched = {}, []
    for name in basin_names:
        key = norm(ALIASES.get(name, name)); hits = by_norm.get(key, [])
        if len(hits) > 1:
            raise ValueError(f"{name}: {len(hits)} rim sites normalize to {key!r}: {hits}")
        if hits:
            matches[name] = hits[0]
        else:
            unmatched.append(name)
    if len(set(matches.values())) != len(matches):
        raise ValueError("two basin assemblages matched the same rim site")
    return matches, unmatched


def basin_rim_table(basin_names) -> tuple[pd.DataFrame, list]:
    """Rim counts for the matched basin assemblages, indexed by basin name in the basin order, with the published counts re-checked."""
    r = load_rim_counts(); matches, unmatched = match_basin(basin_names)
    t = r.loc[[matches[n] for n in basin_names if n in matches]].copy()
    t.insert(0, "rim_site", t.index); t.index = pd.Index([n for n in basin_names if n in matches], name="assemblage")
    pub = pd.read_csv(PUBLISHED_CSV, comment="#").set_index("site")
    for n, site in matches.items():
        for c in PUBLISHED + ("N",):
            if int(pub.loc[site, c]) != int(t.loc[n, c]):
                raise ValueError(f"{n} ({site}): {c} is {t.loc[n, c]} in the count table and {pub.loc[site, c]} as published")
    return t, unmatched


def attribute_summary(t: pd.DataFrame, attributes=FORM_ATTRIBUTES + PASTE) -> pd.DataFrame:
    """Per attribute over the matched sites: rims carrying it, rims examined, pooled frequency, sites with any, sites with all, and whether it is a 1951 type-level distinction (paste)."""
    rows = []
    for a in attributes:
        k = t[a].to_numpy(); n = t["N"].to_numpy()
        rows.append({"attribute": a, "rims_with": int(k.sum()), "rims_examined": int(n.sum()), "pooled_frequency": float(k.sum() / n.sum()),
                     "sites_with_any": int((k > 0).sum()), "sites_with_all": int((k == n).sum()), "sites": int(len(t)),
                     "min_rims": int(n.min()), "median_rims": float(pd.Series(n).median()), "max_rims": int(n.max()),
                     "tabulated_as_1951_type": a in PASTE})
    return pd.DataFrame(rows)
