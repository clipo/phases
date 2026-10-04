"""63_cmv_site_coords.py - geocode the Williams 1954 CMV sites.

Writes `data/processed/williams1954_cmv_coords.tsv`, which analyses 28, 30 and
31 read, and which Figures 10 and 11 therefore depend on.

WHY THIS EXISTS. That table was in the repository as a committed artifact with
no generator: nothing in the tree wrote it, and when it went missing the two
figures could not be rebuilt at all. This is the same shape as the
LMVHydrology provenance defect, an artifact outliving the chain that produced
it, and rule 15 wants the chain committed rather than the memory of it.

HOW THE GEOCODING WORKS. The Williams (1954) assemblages carry an LMS grid
designation (`site_desig`, e.g. "5-S-4"). The digitized Phillips, Ford and
Griffin point files under `data/CMV/` carry the same designation in `NAME1_`
and a site name in `NAME2_`, so the join is:

  1. normalise the designation (strip, upper-case, remove spaces) and look it
     up in the pooled point table; record `match = "grid"`.
  2. failing that, match on the site NAME instead, and record `match = "name"`.

The fallback is not decoration. Beckwith's Fort is "6-T-1" in Williams and
"S-T-1" in PHILLIPS.SHP, a 6/S transcription difference of exactly the kind the
LMS grid designations attract, so it is the one site of 28 that the grid lookup
misses and the name lookup recovers. The `match` column is written out so a
reader can see which sites rest on which join rather than having to trust that
they all rest on the same one.

The three point layers are pooled in the order PHILLIPS, GEOFILE, PFGOUT, and
the first hit for a designation wins, matching `28_macro_boundary.grid_coords`.
Points are restricted to the study window (lat 33-38, lon -92 to -88) before
the join, so a stray digitized point elsewhere cannot capture a designation.

Verified 2026-09-04 to reproduce the previously committed table exactly, all 28
sites and both match types.

Usage: .venv/bin/python analyses/63_cmv_site_coords.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CMV = ROOT / "data" / "CMV"
COUNTS = ROOT / "data" / "raw" / "williams1954_cmv_counts.tsv"
OUT = ROOT / "data" / "processed" / "williams1954_cmv_coords.tsv"

LAYERS = ["PHILLIPS.SHP", "GEOFILE.SHP", "PFGOUT.SHP"]
LAT_RANGE = (33.0, 38.0)
LON_RANGE = (-92.0, -88.0)


def _norm_grid(s) -> str:
    """Match `28_macro_boundary._norm_grid` exactly."""
    return s.strip().upper().replace(" ", "") if isinstance(s, str) else ""


def _norm_name(s) -> str:
    return " ".join(str(s).strip().upper().split()) if s is not None else ""


def point_table() -> pd.DataFrame:
    """Pooled digitized points, with both join keys and coordinates."""
    frames = []
    for f in LAYERS:
        path = CMV / f
        if not path.exists():
            raise FileNotFoundError(
                f"{path.relative_to(ROOT)} is missing. The three digitized point "
                "layers (PHILLIPS, GEOFILE, PFGOUT, with their .SHX/.DBF/.PRJ "
                "sidecars) are tracked; if they are absent the working tree is "
                "incomplete rather than the script being wrong."
            )
        g = gpd.read_file(path).to_crs(4326)
        g["grid"] = g["NAME1_"].map(_norm_grid)
        g["name"] = g["NAME2_"].map(_norm_name)
        g["lat"] = g.geometry.y
        g["lon"] = g.geometry.x
        frames.append(g[["grid", "name", "lat", "lon"]])
    G = pd.concat(frames, ignore_index=True)
    G = G[G.lat.between(*LAT_RANGE) & G.lon.between(*LON_RANGE)]
    return G


def main() -> int:
    counts = pd.read_csv(COUNTS, sep="\t")
    sites = (counts[["site_name", "site_desig", "region"]]
             .drop_duplicates()
             .reset_index(drop=True))

    G = point_table()
    by_grid = G[G.grid != ""].drop_duplicates("grid").set_index("grid")
    by_name = G[G.name != ""].drop_duplicates("name").set_index("name")

    rows = []
    for r in sites.itertuples(index=False):
        key = _norm_grid(str(r.site_desig))
        if key in by_grid.index:
            hit, how = by_grid.loc[key], "grid"
        else:
            nkey = _norm_name(r.site_name)
            if nkey not in by_name.index:
                print(f"  UNMATCHED: {r.site_name!r} ({r.site_desig})")
                continue
            hit, how = by_name.loc[nkey], "name"
        rows.append({"site_name": r.site_name, "site_desig": r.site_desig,
                     "region": r.region, "lat": hit.lat, "lon": hit.lon,
                     "match": how})

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)

    n_grid = int((out["match"] == "grid").sum())
    n_name = int((out["match"] == "name").sum())
    print(f"{len(sites)} Williams sites; {len(out)} geocoded "
          f"({n_grid} by grid designation, {n_name} by site name)")
    if n_name:
        for r in out[out["match"] == "name"].itertuples(index=False):
            print(f"  name-matched: {r.site_name} ({r.site_desig})")
    if len(out) < len(sites):
        print(f"  WARNING: {len(sites) - len(out)} site(s) unmatched")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
