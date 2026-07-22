"""36_canonical_phase_map.py — the canonical "phase = bounded area" assumption.

Reproduces the standard culture-history map in which the central Mississippi
Valley is partitioned into bounded phase territories, using the published
site-to-phase assignments of Mainfort (1996: Figure 1; the same scheme appears,
for a subset of the phases, in Azar & Steponaitis 2022:28). Each Mainfort-PFG
assemblage is given its Mainfort phase (Parkin, Nodena, Jones Bayou, Tipton,
Kent, Walls, Parchman); assemblages Mainfort does not place on that map are left
'unassigned' and drawn as open symbols. The phase territories are a Voronoi
partition of the assigned deposits, dissolved by phase and clipped to a loose
envelope around the sites, giving the bounded phase-territory picture the phase
concept assumes, the foil for the drift-generated graded membership of Figure 8
and Figure S8. Each territory is shrunk inward by a fixed margin so adjacent
phases are separated by a visible gap rather than sharing a Voronoi edge.

Writes figures/fig1_phases.png (manuscript Figure 1).

Usage: PYTHONPATH=src python3 analyses/36_canonical_phase_map.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import geopandas as gpd  # noqa: E402
from shapely.geometry import Point, MultiPoint, box  # noqa: E402
from shapely.ops import unary_union, voronoi_diagram  # noqa: E402
import make_figures as mf  # noqa: E402
import make_map as mm  # noqa: E402
m35 = importlib.import_module("35_basin_pullout")

OUT_FIG = ROOT / "figures" / "fig1_phases.png"

# Each phase territory is shrunk inward by this many meters before plotting, so
# adjacent phases show a clear gap instead of sharing a Voronoi boundary line.
PHASE_GAP_M = 2500.0

PHASES = ["Parkin", "Nodena", "Jones Bayou", "Tipton", "Kent", "Walls",
          "Parchman"]
PHASE_COL = {"Parkin": "#0072B2", "Nodena": "#E69F00", "Jones Bayou": "#D55E00",
             "Tipton": "#56B4E9", "Kent": "#009E73", "Walls": "#CC79A7",
             "Parchman": "#F0E442"}
# A faint gray per phase: the territory fills are kept very light so the
# hydrology and county lines underneath stay visible; the per-phase marker shape
# and the territory label carry the phase identity.
PHASE_GRAY = {"Parkin": 0.88, "Nodena": 0.95, "Jones Bayou": 0.91,
              "Tipton": 0.965, "Kent": 0.90, "Walls": 0.93, "Parchman": 0.87}
# A distinct closed marker per phase so each deposit's phase reads at a glance in
# grayscale (the Parkin type-site is overlaid with a star separately).
PHASE_MARKER = {"Parkin": "o", "Nodena": "s", "Jones Bayou": "^", "Tipton": "v",
                "Kent": "D", "Walls": "P", "Parchman": "X"}

# Published site-to-phase assignments transcribed from Mainfort (1996: Figure 1,
# "Late period sites and phases in the Central Mississippi Valley"). Assemblages
# Mainfort does not place on that map are left 'unassigned'.
MAINFORT = {
    "Turnbow": "Parkin", "Fortune": "Parkin", "Neeleys_Ferry": "Parkin",
    "Barton_Ranch": "Parkin", "Williamson": "Parkin", "Vernon_Paul": "Parkin",
    "Parkin": "Parkin", "Rose_Mound": "Parkin",
    "Upper_Nodena": "Nodena", "Carson_Lake": "Nodena", "Notgrass": "Nodena",
    "40LA007": "Jones Bayou", "Porter": "Jones Bayou", "Jones_Bayou": "Jones Bayou",
    "Fullen": "Jones Bayou",
    "Graves_Lake": "Tipton", "Hatchie": "Tipton", "Richardsons_Landing": "Tipton",
    "Wilder": "Tipton", "Rast": "Tipton", "Jeter": "Tipton",
    "Castile_Landing": "Kent", "Cramor_Place": "Kent", "Nickel": "Kent",
    "Davis": "Kent", "Clay_Hill": "Kent", "Kent_Place": "Kent", "Soudan": "Kent",
    "Grant": "Kent", "Starkley": "Kent",
    "Mound_Place": "Walls", "Young": "Walls", "Belle_Meade": "Walls",
    "Walls": "Walls", "Chuccalissa": "Walls", "Woodlyn": "Walls", "Beck": "Walls",
    "Irby": "Walls", "Lake_Cormorant": "Walls", "Commerce": "Walls",
    "Hollywood": "Walls",
    "West_Mounds": "Parchman", "Salomon": "Parchman", "Parchman": "Parchman",
}


def assign_phases(names, coords=None):
    """Mainfort (1996: Figure 1) phase label per assemblage; assemblages Mainfort
    does not place on that map are returned as 'unassigned'. (coords is accepted
    for a stable call signature but not used.)"""
    return np.array([MAINFORT.get(str(nm), "unassigned") for nm in names])


def main():
    counts_df, coords_df = m35.load_full()
    names = [str(x) for x in counts_df.index]
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    basin = mf._basin_members("curated")

    labels = assign_phases(names)
    lab_idx = np.array([PHASES.index(p) if p in PHASES else -1 for p in labels])
    assigned = lab_idx >= 0

    lon, lat = coords[:, 1], coords[:, 0]
    gp = gpd.GeoSeries(gpd.points_from_xy(lon, lat),
                       crs="EPSG:4326").to_crs("EPSG:26915")
    E, Nm = gp.x.to_numpy(), gp.y.to_numpy()
    margin = 10_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)
    is_basin = np.array([nm in basin for nm in names])

    fig, ax = plt.subplots(figsize=(5.8, 6.8))
    # No land tone behind the map: a white background makes the light-gray phase
    # territories stand out clearly.
    mm.basin_basemap(ax, ext, geology=False, grayscale=True,
                     show_counties=False, show_states=False, draw_rivers=False)

    # Consistent HydroRIVERS hydrography (shared with Figures 9 and 10): rivers
    # to flow order 6, width scaled by order, Mississippi main stem in blue.
    mm.draw_hydrorivers(ax, ext, max_ord=6, main_ord=3, grayscale=False, zorder=4.4)

    # State labels only, with no political boundary lines: the Mississippi River
    # carries the eastern edge as water rather than as a border, and county
    # lines are omitted. Labels are placed by hand in clear areas (the
    # assemblages fall in Arkansas, Tennessee, and Mississippi; Missouri lies
    # just north of the mapped area).
    st = mm._clip_gdf(mm._load_states(), ext).dissolve(by="STATE")
    STATE_LBL = {"AR": (0.20, 0.46), "TN": (0.92, 0.74), "MS": (0.40, 0.16)}
    for abbr, (fx, fy) in STATE_LBL.items():
        if abbr not in st.index:
            continue
        ax.text(fx, fy, abbr, transform=ax.transAxes, fontsize=11,
                fontweight="bold", color="0.3", ha="center", va="center",
                zorder=12, path_effects=[pe.withStroke(linewidth=3.0,
                                                       foreground="white")])

    # Phase territories as a non-overlapping partition: a Voronoi tessellation of
    # the assigned deposits, dissolved by phase, then clipped to a loose envelope
    # around the deposits so the territories hug the sites instead of running to
    # the map edge. Each phase gets a distinct gray (American Antiquity prints
    # without color); the territories do not overlap and the labels name them.
    asg = np.where(assigned)[0]
    pts = [Point(E[i], Nm[i]) for i in asg]
    cells = list(voronoi_diagram(MultiPoint(pts),
                                 envelope=box(ext[0], ext[2], ext[1], ext[3])).geoms)
    cell_phase = []
    for cell in cells:
        ph_here = -1
        for j, p in enumerate(pts):
            if cell.intersects(p):
                ph_here = int(lab_idx[asg[j]])
                break
        cell_phase.append(ph_here)
    envelope = unary_union([p.buffer(16_000) for p in pts])
    for k, ph in enumerate(PHASES):
        member = [cells[c] for c in range(len(cells)) if cell_phase[c] == k]
        if not member:
            continue
        poly = unary_union(member).intersection(envelope)
        # Shrink the territory inward so neighboring phases are visibly separated.
        poly = poly.buffer(-PHASE_GAP_M)
        if poly.is_empty:
            continue
        geoms = poly.geoms if poly.geom_type == "MultiPolygon" else [poly]
        for g in geoms:
            if g.is_empty or g.geom_type != "Polygon":
                continue
            x, y = g.exterior.xy
            ax.fill(x, y, facecolor=str(PHASE_GRAY[ph]), edgecolor="0.2",
                    linewidth=0.8, zorder=1.5)

    # Deposits, one closed marker shape per phase (assemblages Mainfort does not
    # place are dropped). Markers are black with a thin white edge so they read on
    # the light fills and over the hydrology.
    for k, ph in enumerate(PHASES):
        m = lab_idx == k
        if not m.any():
            continue
        ax.scatter(E[m], Nm[m], s=26, c="black", marker=PHASE_MARKER[ph],
                   edgecolor="white", linewidth=0.4, zorder=10)
    for k, ph in enumerate(PHASES):
        m = lab_idx == k
        if m.sum() == 0:
            continue
        ax.text(E[m].mean(), Nm[m].mean(), ph, fontsize=8, fontweight="bold",
                ha="center", va="center", zorder=13, color="black",
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    pk = names.index("Parkin")
    ax.scatter([E[pk]], [Nm[pk]], marker="*", s=210, c="white", edgecolor="black",
               linewidth=0.8, zorder=14)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])

    # Legend: the marker shape for each phase, plus the basin set and type-site.
    handles = [Line2D([0], [0], marker=PHASE_MARKER[ph], color="none",
                      markerfacecolor="black", markeredgecolor="white",
                      markeredgewidth=0.4, markersize=6, label=ph) for ph in PHASES]
    handles += [
        Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
               markeredgecolor="black", markersize=12, label="Parkin (type-site)"),
    ]
    ax.legend(handles=handles, fontsize=6.4, loc="lower center",
              bbox_to_anchor=(0.5, 1.005), ncol=3, framealpha=0.0,
              handletextpad=0.4, columnspacing=1.2, labelspacing=0.4,
              borderpad=0.2)

    # 25 km scale bar in the lower-left corner (the lower-right is reserved for
    # the locator inset).
    bar = 25_000.0 / (ext[1] - ext[0])
    x0, y0 = 0.06, 0.05
    ax.plot([x0, x0 + bar], [y0, y0], transform=ax.transAxes, color="black", lw=2)
    ax.text(x0 + bar / 2, y0 + 0.015, "25 km", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=7)

    # North America locator inset, placed in the blank space east of the
    # Mississippi (lower-right of the data axes). Position is taken from the
    # finalized axes box so it tracks the equal-aspect map under bbox_inches.
    fig.canvas.draw()
    pos = ax.get_position()
    inset_rect = [pos.x0 + 0.605 * pos.width, pos.y0 + 0.02 * pos.height,
                  0.38 * pos.width, 0.30 * pos.height]
    mm._add_na_inset(fig, ext, rect=inset_rect, grayscale=True)

    mf.save_all(fig, OUT_FIG)
    plt.close(fig)
    counts = {p: int((labels == p).sum()) for p in PHASES}
    basin_phases = {p: int(((labels == p) & is_basin).sum()) for p in PHASES}
    print("phase counts (all 55):", counts)
    print("basin-set (29) by phase:", {k: v for k, v in basin_phases.items() if v})
    print(f"wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
