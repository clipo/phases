#!/usr/bin/env python3
"""Clip HydroRIVERS to the study area, so the map figures are reproducible.

Figures 1, 9 and 10 all draw hydrography through `make_map.draw_hydrorivers`,
which reads `data/HydroRIVERS_v10_na.gdb.zip` [Lehner & Grill 2013]. That source
is 70 MB of continental North America, and it is not in the repository: a clean
clone cannot regenerate three cited figures.

Tracking 70 MB to draw rivers across two degrees of the lower Mississippi valley
is the wrong trade. This writes the clip instead, which is the option
`output/findings/basin_membership_provenance.md` records as option 2 for exactly
this situation: commit the small derived intermediate, with the cut clearly
visible, rather than the bulk source.

The clip keeps every reach of flow order <= MAX_ORD intersecting the study
bounding box, in the source CRS, with the ORD_FLOW attribute the drawing code
uses for line weight. `draw_hydrorivers` prefers the clip when it is present and
falls back to the full source, so nothing breaks for anyone who has the 70 MB
file and the figures are byte-comparable either way.

A NOTE ON REPRODUCIBILITY. Re-running this does NOT reproduce the committed
file byte for byte, because GeoPackage embeds a creation timestamp in the
container. The CONTENT is identical and was checked that way on 2026-09-04:
same 9,927 reaches, identical attribute table, and geometries exactly equal
under `geom_equals_exact(tolerance=0)`. Compare content, not checksums.

Usage: .venv/bin/python scripts/clip_hydrorivers.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "HydroRIVERS_v10_na.gdb.zip"
OUT = ROOT / "data" / "Shapefiles" / "hydrorivers_lmv_cmv.gpkg"

# Study bounding box in lon/lat, padded well beyond every figure extent so a
# later change of extent does not silently clip the hydrography. The basin
# spans about 128 km; this covers the whole lower valley from the Cairo Lowland
# to below Vicksburg.
BBOX = (-92.5, 31.5, -87.5, 38.5)
MAX_ORD = 6  # the highest flow order any figure draws


def main() -> int:
    if not SRC.exists():
        print(f"missing {SRC.relative_to(ROOT)}")
        print("Obtain HydroRIVERS v1.0 (North America) from "
              "https://www.hydrosheds.org/products/hydrorivers and place it "
              "there. It is deliberately not tracked; the clip this script "
              "writes is what the repository carries.")
        return 1

    riv = gpd.read_file(SRC, bbox=BBOX)
    before = len(riv)
    riv = riv[riv["ORD_FLOW"] <= MAX_ORD]
    keep = [c for c in ("ORD_FLOW", "ORD_STRA", "MAIN_RIV", "geometry")
            if c in riv.columns]
    riv = riv[keep]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    riv.to_file(OUT, driver="GPKG")
    size_mb = OUT.stat().st_size / 1e6
    print(f"{before} reaches in bbox -> {len(riv)} at flow order <= {MAX_ORD}")
    print(f"wrote {OUT.relative_to(ROOT)} ({size_mb:.1f} MB), CRS {riv.crs}")
    print(f"columns: {list(riv.columns)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
