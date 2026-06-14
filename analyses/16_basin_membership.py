"""16_basin_membership.py — canonical drainage-based basin membership.

Defines the St. Francis basin (the Parkin-phase analysis unit) HYDROLOGICALLY:
assemblages/sites within DRAINAGE_KM of the St. Francis / Tyronza / L'Anguille
drainage (LMVHydrology). This replaces the earlier latitude cut (lat >= 34.5),
which admitted Mississippi-River-proximal sites (Walls, Chuccalissa, Hollywood,
Upper Nodena) that are 30+ km from the St. Francis drainage and within a few km
of the Mississippi. The stricter nearest-drainage watershed rule (closer to the
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

    # broad settlement set (coords are UTM Easting/Northing)
    broad = load_pfg_counts(ROOT / "data" / "raw" / "PFGData_sherds.csv")
    if not broad.index.is_unique:
        broad = broad.groupby(level=0).sum()
    lmv = load_lmv(ROOT / "data" / "LMVData_locations.csv")
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
    print(f"  broad basin:   {len(broad_members)} sites")
    print(f"  wrote {PROC/'basin_members_curated.txt'} and basin_members_broad.txt")
    print("  curated members:", ", ".join(cur_members))


if __name__ == "__main__":
    main()
