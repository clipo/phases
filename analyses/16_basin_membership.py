"""16_basin_membership.py — canonical basin membership, by phase.

WHAT THE ANALYSIS UNIT IS (author ruling, 2026-09-20). The paper tests the
phases of the St. Francis basin, so membership is phase membership: an
assemblage is in the curated set when Mainfort (1996, Figure 1) assigns it to
one of the five phases the Data section names — Parkin, Nodena, Kent, Walls or
Parchman. The phases as drawn spill beyond the drainage, and the unit under
test is the phase scheme, not the watershed.

This replaces a 20 km corridor around the St. Francis / Tyronza / L'Anguille
drainage, which was a geographic proxy for the same idea and disagreed with it
both ways: it admitted five assemblages Mainfort places in no phase, and
excluded nine that carry a phase label. The corridor distance is still measured
and reported below, because how far the phases reach beyond the drainage is
worth stating rather than assuming.

The older, drainage-based note follows, kept because the measurements in it are
still the record of what the corridor rule did.

16_basin_membership.py — canonical drainage-based basin membership.

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

import numpy as np
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

# The five phases the Data section names for the basin. Tipton and Jones Bayou,
# the other two in Figure 1, sit north and east of this scheme and are not part
# of the unit under test.
# Author rulings, 2026-09-21. Parchman is out on geography: its three
# assemblages lie across the Mississippi in Coahoma County (Parchman itself is
# SW1/4 NW1/4 S30 T29N R3W, the Mississippi survey grid), 30 to 46 km from the
# nearest other assemblage and 77 km from the centroid of the rest. It is not a
# phase of the St. Francis basin under any definition, and its 236 decorated
# sherds carried a third of the drift residual on their own.
ST_FRANCIS_PHASES = ("Parkin", "Nodena", "Kent", "Walls")
# A minimum decorated-sherd count, because a between-cluster term computed on
# a few dozen sherds is mostly noise and an "excess over drift" built on it is
# a ratio of two small numbers. 75 removes Cheatham (19), Connor (38),
# Norfolk (49), Notgrass (38), Pouncey (62) and Upper Nodena (47).
MIN_DECORATED = 75
# A phase left with fewer members than this cannot carry a between-group term
# and, as a singleton, would be handed its own cluster and inflate F_ST for an
# arithmetic reason. Nodena is the case: the 75-sherd rule leaves Carson Lake
# alone, and it is treated as unassigned rather than as a phase of one.
MIN_PHASE_MEMBERS = 2
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

    # Membership is Mainfort's phase label, not distance to water.
    ph = importlib.import_module("36_canonical_phase_map")
    labels, derived = ph.assign_phases_by_territory(
        [str(i) for i in coords.index], coords[["Latitude", "Longitude"]].to_numpy(float))
    labels = np.asarray(labels)
    in_phase = np.isin(labels, ST_FRANCIS_PHASES)
    n_sherds = counts.reindex(coords.index).sum(axis=1).to_numpy(float)
    too_small = n_sherds < MIN_DECORATED
    print(f"  under {MIN_DECORATED} decorated sherds, dropped: "
          f"{', '.join(f'{a} ({int(k)})' for a, k in zip(coords.index[in_phase & too_small], n_sherds[in_phase & too_small]))}")
    in_phase &= ~too_small
    for phase in ST_FRANCIS_PHASES:
        members = coords.index[in_phase & (labels == phase)]
        if 0 < len(members) < MIN_PHASE_MEMBERS:
            print(f"  {phase} is left with {len(members)} member ({', '.join(members)}); "
                  f"fewer than {MIN_PHASE_MEMBERS}, so it is not carried as a phase")
            in_phase &= ~(labels == phase)
    n_derived = int((in_phase & derived).sum())
    print(f"  phase labels: {int(in_phase.sum()) - n_derived} from Mainfort's map, "
          f"{n_derived} by territory containment "
          f"({', '.join(sorted(coords.index[in_phase & derived]))})")
    cur_members = sorted(coords.index[in_phase])
    (PROC / "basin_members_curated.txt").write_text("\n".join(cur_members) + "\n")

    # How far the phase set reaches beyond the old corridor, stated rather than
    # assumed: these are the assemblages on which the two rules disagree.
    added = sorted(coords.index[in_phase & (cd > DRAINAGE_KM)])
    dropped = sorted(coords.index[~in_phase & (cd <= DRAINAGE_KM)])
    print(f"  phase rule: {len(cur_members)} assemblages in "
          f"{', '.join(ST_FRANCIS_PHASES)}")
    print(f"    beyond the {DRAINAGE_KM:.0f} km corridor but in a phase ({len(added)}): "
          f"{', '.join(added) if added else 'none'}")
    print(f"    inside the corridor but in no phase ({len(dropped)}): "
          f"{', '.join(dropped) if dropped else 'none'}")
    if len(added):
        far = float(cd[in_phase].max())
        print(f"    farthest phase member from the drainage: {far:.1f} km")

    # The stricter watershed set, used as the robustness check the main text
    # cites. An assemblage qualifies when it is inside the corridor AND nearer
    # to the St. Francis system than to the Mississippi. The docstring above has
    # asserted n=19 for this set for some time without any code computing it,
    # so the manuscript's robustness claim rested on a number the repository
    # could not produce (rule 1). It is derived and written here instead. The
    # count reproduces: 29 corridor, 19 watershed.
    md = cp.geometry.distance(mississippi_geom()) / 1000.0
    # The robustness set must be a SUBSET of the unit under test, or it tests
    # something else. It is the phase members that are also hydrologically
    # core: nearer the St. Francis system than the Mississippi. Under the old
    # corridor rule this line ignored phase membership, which after the rule
    # changed left it holding assemblages the primary set no longer contains.
    ws_members = sorted(coords.index[in_phase & (cd < md)])
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
