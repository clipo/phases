"""16_basin_membership.py — canonical drainage-based basin membership.

Defines the St. Francis basin (the Parkin-phase analysis unit) HYDROLOGICALLY:
assemblages/sites within DRAINAGE_KM of the St. Francis / Tyronza / L'Anguille
drainage (LMVHydrology). This replaces the earlier latitude cut (lat >= 34.5),
which admitted Mississippi-River-proximal sites well outside the drainage.

Measured distances to the St. Francis system, 2026-08-31 (this script's own
geometry, DRAINAGE_KM = 20):

    Walls          35.25 km   excluded
    Chuccalissa    41.25 km   excluded
    Hollywood      18.50 km   RETAINED
    Upper_Nodena   18.24 km   RETAINED
    Parkin          0.19 km   retained

An earlier version of this docstring named all four of Walls, Chuccalissa,
Hollywood and Upper Nodena as "30+ km from the St. Francis drainage" and so as
sites the hydrological rule removes. That is wrong for Hollywood and Upper
Nodena, which sit at roughly 18 km and are inside the corridor. Corrected here
against measurement rather than recollection (Verification Regime rules 1 and
2). The substantive point stands for Walls and Chuccalissa. The stricter nearest-drainage watershed rule (closer to the
St. Francis system than to the Mississippi) gives a smaller set (n=19 curated);
the 20 km corridor (n=29) is used as primary for adequate sample size, with the
watershed set available as a robustness check.

Writes the member lists to data/processed/basin_members_curated.txt and
data/processed/basin_members_broad.txt (one id per line), which the analysis
pipeline reads. Read-only on the manuscript.

Usage: .venv/bin/python analyses/16_basin_membership.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import geopandas as gpd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
from mls_emergence.dataio.pfg import load_pfg_counts  # noqa: E402
from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv  # noqa: E402
grid = importlib.import_module("12_sensitivity_grid")  # noqa: E402

PROC = ROOT / "data" / "processed"
SHP = ROOT / "data" / "Shapefiles"
DRAINAGE = {"Saint Francis River", "Tyronza River", "L'Anguille River", "Saint Francis Floodway"}
DRAINAGE_KM = 20.0
UTM = 26915


def sf_geom():
    hydro = gpd.read_file(SHP / "LMVHydrology.shp")
    return hydro[hydro["NAME"].isin(DRAINAGE)].to_crs(epsg=UTM).geometry.union_all()


def mississippi_geom():
    """The Mississippi channel, for the nearest-drainage watershed rule.

    It is not available from LMVHydrology: that layer maps large rivers as
    unnamed bank and shoreline features, so `NAME` carries no Mississippi entry
    at all. LMVMajorRivers carries it as a named polygon, which is why the
    watershed rule needs that layer and the corridor rule does not.
    """
    mr = gpd.read_file(SHP / "LMVMajorRivers.shp")
    sel = mr[mr["RIVER_NAME"] == "Mississippi River"]
    if sel.empty:
        raise ValueError("no 'Mississippi River' feature in LMVMajorRivers.shp")
    return sel.to_crs(epsg=UTM).geometry.union_all()


def main():
    PROC.mkdir(parents=True, exist_ok=True)
    sf = sf_geom()

    # curated decorated set (coords are WGS84 lat/long)
    counts, coords = grid.load_curated_full()
    coords = coords.dropna()
    cp = gpd.GeoDataFrame(
        coords.copy(),
        geometry=gpd.points_from_xy(coords["Longitude"], coords["Latitude"]),
        crs="EPSG:4326").to_crs(epsg=UTM)
    cd = cp.geometry.distance(sf) / 1000.0
    cur_members = sorted(coords.index[cd <= DRAINAGE_KM])
    (PROC / "basin_members_curated.txt").write_text("\n".join(cur_members) + "\n")

    # The stricter watershed set, used as the robustness check the main text
    # cites. An assemblage qualifies when it is inside the corridor AND nearer
    # to the St. Francis system than to the Mississippi. The docstring above has
    # asserted n=19 for this set for some time without any code computing it,
    # so the manuscript's robustness claim rested on a number the repository
    # could not produce (rule 1). It is derived and written here instead. The
    # count reproduces: 29 corridor, 19 watershed.
    md = cp.geometry.distance(mississippi_geom()) / 1000.0
    ws_members = sorted(coords.index[(cd <= DRAINAGE_KM) & (cd < md)])
    (PROC / "basin_members_watershed.txt").write_text("\n".join(ws_members) + "\n")

    # broad settlement set (coords are UTM Easting/Northing)
    broad = load_pfg_counts(ROOT / "data" / "raw" / "PFGData.xlsx")
    if not broad.index.is_unique:
        broad = broad.groupby(level=0).sum()
    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    joined, _ = join_pfg_to_lmv(broad, lmv)
    bm = joined.dropna(subset=["Easting", "Northing"]).copy()
    bp = gpd.GeoDataFrame(
        bm.copy(),
        geometry=gpd.points_from_xy(bm["Easting"], bm["Northing"]),
        crs=f"EPSG:{UTM}")
    bd = bp.geometry.distance(sf) / 1000.0
    broad_members = sorted(str(i) for i in bm.index[bd <= DRAINAGE_KM])
    (PROC / "basin_members_broad.txt").write_text("\n".join(broad_members) + "\n")

    print(f"drainage corridor <= {DRAINAGE_KM:.0f} km of the St. Francis system")
    print(f"  curated basin: {len(cur_members)} assemblages")
    print(f"  watershed subset (nearer St. Francis than Mississippi): "
          f"{len(ws_members)} assemblages")
    print(f"  broad basin:   {len(broad_members)} sites")
    print(f"  wrote {PROC/'basin_members_curated.txt'}, basin_members_watershed.txt and basin_members_broad.txt")
    print("  curated members:", ", ".join(cur_members))


if __name__ == "__main__":
    main()
