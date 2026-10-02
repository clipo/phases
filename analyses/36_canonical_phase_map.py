"""36_canonical_phase_map.py — the canonical "phase = bounded area" assumption.

Reproduces the standard culture-history map in which the central Mississippi
Valley is partitioned into bounded phase territories, using the published
site-to-phase assignments of Mainfort (1996: Figure 1; the same scheme appears,
for a subset of the phases, in Azar & Steponaitis 2022:28). Each Mainfort-PFG
assemblage is given its Mainfort phase (Parkin, Nodena, Jones Bayou, Tipton,
Kent, Walls, Parchman); assemblages Mainfort does not place on that map are left
'unassigned' and are not drawn. Jones Bayou and Tipton have no assemblage with
coordinates in the matrix, so they appear in the legend but draw no territory. The phase territories are a Voronoi
partition of the assigned deposits, dissolved by phase and clipped to a loose
envelope around the sites, giving the bounded phase-territory picture the phase
concept assumes, the foil for the graded membership that copying across
distance generates. Each territory is shrunk inward by a fixed margin so adjacent
phases are separated by a visible gap rather than sharing a Voronoi edge.

Writes figures/fig1_phases.png. This was the manuscript's Figure 1 until
2026-10-01. Figure 1 is now the Phillips (1970) phase areas
(93_phillips1970_phase_map.py), and the territory map shown in the paper is the
one built on Mainfort's (2003) assignments for the contrast case
(94_mainfort2003_phase_map.py, Figure 10). This map is kept as the standalone
view of the territories the primary analysis assigns by, whose edges Figures 4B
and 7A overlay; both later maps draw with the helpers defined here.

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


def assign_phases_by_territory(names, coords):
    """Phase label per assemblage, filling Mainfort's gaps by territory.

    Mainfort's map covers 44 of the 55 curated assemblages. The other eleven are
    not placed outside the phases; they are simply not on his map, several of
    them because they were collected later (Holden Lake is one of Lipo's 1996-97
    collections and has no PFG site number). Dropping them would exclude
    assemblages on the ground the phases claim, which is the ground the test is
    about.

    So an unmapped assemblage takes the phase of the territory it falls in,
    where territories are built exactly as this module's map builds them: every point
    belongs to its nearest assemblage, and those areas are merged by phase. The
    reasoning is the paper's own. A phase is drawn as a bounded area, and the
    claim under test is that the area was a community. Assemblages inside it
    therefore belong in the test: if the phase is a community, they should
    resemble one another, and if it is not, they should not. Leaving them out
    because a 1996 map did not list them would decide part of the question by
    omission.

    Returns (labels, derived), where `derived` marks the labels this function
    supplied rather than Mainfort. Quote the two separately when it matters: a
    territory label is our construction, not his determination (rule 13).
    """
    labels = assign_phases(names)
    xy = np.asarray(coords, dtype=float)
    mapped = labels != "unassigned"
    derived = ~mapped
    if not mapped.any():
        return labels, derived
    for i in np.where(derived)[0]:
        d = np.hypot((xy[mapped, 0] - xy[i, 0]) * 111.32,
                     (xy[mapped, 1] - xy[i, 1]) * 111.32
                     * np.cos(np.radians(xy[i, 0])))
        labels[i] = labels[mapped][int(np.argmin(d))]
    return labels, derived


# Phillips's (1970) own phase areas, as redrawn by Lipo (2001: Figure 2.6) and
# read into coordinates by scripts/build_phillips1970_outlines.py.
PHILLIPS_OUTLINES = ROOT / "data" / "raw" / "phillips1970_phase_outlines.csv"
PHILLIPS_PHASES = ["Nodena", "Parkin", "Walls", "Kent"]
# The frame of the source figure, LMS columns N to Q by rows 10 to 14:
# (lon west, lon east, lat south, lat north).
PHILLIPS_FRAME = (-90.75, -89.75, 34.50, 35.75)


def load_phillips_outlines():
    """{phase: shapely Polygon in UTM 15N} for the four Phillips (1970) phase
    areas, each closed first vertex to last."""
    import pandas as pd
    from shapely.geometry import Polygon
    if not PHILLIPS_OUTLINES.exists():
        raise SystemExit(f"{PHILLIPS_OUTLINES} is missing; run "
                         "scripts/build_phillips1970_outlines.py --pdf <Lipo 2001>")
    df = pd.read_csv(PHILLIPS_OUTLINES, comment="#")
    missing = [p for p in PHILLIPS_PHASES if p not in set(df["phase"])]
    if missing:
        raise SystemExit(f"{PHILLIPS_OUTLINES.name} has no outline for {missing}")
    out = {}
    for ph in PHILLIPS_PHASES:
        v = df[df["phase"] == ph].sort_values("vertex")
        x, y = to_utm(v["longitude"].to_numpy(), v["latitude"].to_numpy())
        poly = Polygon(np.column_stack([x, y]))
        if not poly.is_valid:
            raise SystemExit(f"the {ph} outline is not a simple polygon")
        out[ph] = poly
    return out


def assign_phases_phillips(names, coords):
    """Phase label per assemblage from Phillips's (1970) phase areas.

    An assemblage takes the phase of the outline it lies in. One that lies in
    none, but inside the frame of the source figure, takes the phase of the
    NEAREST outline, on the same reasoning as `assign_phases_by_territory`:
    the areas are drawn as bounded ground, and leaving out an assemblage on or
    beside that ground because the line misses it would decide part of the
    question by omission. An assemblage outside the frame (the Coahoma County
    sites to the south) is 'unassigned': the figure does not cover it.

    Returns (labels, derived), `derived` marking the labels supplied by
    nearness rather than by containment, with `km`, the distance to the
    outline of the phase assigned (0 inside), available from
    `phillips_outline_distance`. Quote the two kinds separately (rule 13).
    `names` is accepted for a stable call signature and is not used: the
    assignment is by position alone.
    """
    xy = np.asarray(coords, dtype=float)
    outline = load_phillips_outlines()
    E, N = to_utm(xy[:, 1], xy[:, 0])
    w, e, s, n = PHILLIPS_FRAME
    in_frame = (xy[:, 1] > w) & (xy[:, 1] < e) & (xy[:, 0] > s) & (xy[:, 0] < n)
    labels, derived = [], []
    for i in range(len(xy)):
        pt = Point(E[i], N[i])
        inside = [ph for ph in PHILLIPS_PHASES if outline[ph].contains(pt)]
        if len(inside) > 1:
            raise ValueError(f"{names[i]} lies inside {inside}: the outlines overlap")
        if inside:
            labels.append(inside[0]); derived.append(False)
        elif in_frame[i]:
            labels.append(min(PHILLIPS_PHASES, key=lambda ph: outline[ph].exterior.distance(pt)))
            derived.append(True)
        else:
            labels.append("unassigned"); derived.append(True)
    return np.array(labels), np.array(derived)


def assign_primary_phases(names, coords):
    """The phase labels the primary analysis tests: Phillips's (1970) phases.

    Author ruling, 2026-10-01: the primary analysis follows Lipo (2001), the
    survey counts plus his collections under the phases Phillips formulated,
    and is independent of Mainfort. Every analysis takes its labels from here,
    so the scheme is chosen in one place. Until that date the pipeline used
    `assign_phases_by_territory` (Mainfort 1996: Figure 1), which differs from
    Phillips at five of the 28 analyzed assemblages; it is kept for
    `95_phillips_assignment_impact.py`, which records what the change did.

    Returns (labels, derived) exactly as `assign_phases_phillips` does.
    """
    return assign_phases_phillips(names, coords)


def phillips_outline_distance(coords, labels):
    """Kilometers from each assemblage to the outline of the Phillips phase in
    `labels` (0 when inside it, NaN when the label names no outline)."""
    xy = np.asarray(coords, dtype=float)
    outline = load_phillips_outlines()
    E, N = to_utm(xy[:, 1], xy[:, 0])
    return np.array([outline[l].distance(Point(E[i], N[i])) / 1000.0 if l in outline
                     else np.nan for i, l in enumerate(labels)])


# ---------------------------------------------------------------------------
# Map furniture shared by the phase maps (this script, 93 and 94), so the three
# are drawn on one base and cannot drift apart.
# ---------------------------------------------------------------------------
def to_utm(lon, lat):
    """Longitude/latitude arrays -> UTM zone 15N easting, northing (meters)."""
    gp = gpd.GeoSeries(gpd.points_from_xy(lon, lat),
                       crs="EPSG:4326").to_crs("EPSG:26915")
    return gp.x.to_numpy(), gp.y.to_numpy()


def draw_base(ax, ext, state_lbl):
    """White base, HydroRIVERS hydrography and hand-placed state labels.

    `state_lbl` maps a state abbreviation to its (x, y) position in axes
    fractions; a state not inside `ext` is skipped."""
    # No land tone behind the map: a white background makes the light-gray phase
    # territories stand out clearly.
    mm.basin_basemap(ax, ext, geology=False, grayscale=True,
                     show_counties=False, show_states=False, draw_rivers=False)

    # Consistent HydroRIVERS hydrography (one layer for every phase map): rivers
    # to flow order 6, width scaled by order, Mississippi main stem in blue.
    mm.draw_hydrorivers(ax, ext, max_ord=6, main_ord=3, grayscale=False, zorder=4.4)

    # State labels only, with no political boundary lines: the Mississippi River
    # carries the eastern edge as water rather than as a border, and county
    # lines are omitted. Labels are placed by hand in clear areas.
    st = mm._clip_gdf(mm._load_states(), ext).dissolve(by="STATE")
    for abbr, (fx, fy) in state_lbl.items():
        if abbr not in st.index:
            continue
        ax.text(fx, fy, abbr, transform=ax.transAxes, fontsize=11,
                fontweight="bold", color="0.3", ha="center", va="center",
                zorder=12, path_effects=[pe.withStroke(linewidth=3.0,
                                                       foreground="white")])


def voronoi_territories(E, Nm, lab_idx, phases, ext):
    """Phase territories as a non-overlapping partition: a Voronoi tessellation
    of the assigned deposits (lab_idx >= 0), dissolved by phase, clipped to a
    loose envelope around the deposits so the territories hug the sites instead
    of running to the map edge, and shrunk inward by PHASE_GAP_M so neighboring
    phases are visibly separated. Returns {phase: polygon}; a phase with no
    assigned deposit has no entry."""
    asg = np.where(lab_idx >= 0)[0]
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
    territory = {}
    for k, ph in enumerate(phases):
        member = [cells[c] for c in range(len(cells)) if cell_phase[c] == k]
        if not member:
            continue
        poly = unary_union(member).intersection(envelope)
        poly = poly.buffer(-PHASE_GAP_M)
        if poly.is_empty:
            continue
        territory[ph] = poly
    return territory


def fill_polygon(ax, poly, **kw):
    """Fill each Polygon part of `poly` (holes are not drawn)."""
    geoms = poly.geoms if poly.geom_type == "MultiPolygon" else [poly]
    for g in geoms:
        if g.is_empty or g.geom_type != "Polygon":
            continue
        x, y = g.exterior.xy
        ax.fill(x, y, **kw)


def place_phase_labels(fig, ax, ext, E, Nm, anchors, inside, spill=()):
    """Write each phase name at the clear position nearest its anchor.

    `anchors` maps phase -> (x, y), the position the label would ideally take
    (the phase's mean deposit position); `inside` maps phase -> the polygon the
    label box must stay within (None for no constraint). Drawn at the anchor a
    label sat on its own markers, so each goes to the point, nearest the anchor,
    whose text box clears every marker in (E, Nm), every label already on the
    map (rivers, states, earlier phases), and the map edge. Raises if a phase
    has no such position, unless it is named in `spill`: a territory too small
    to hold its own name then takes the nearest clear position outside it."""
    fig.canvas.draw()
    inv = ax.transData.inverted()
    _r = fig.canvas.get_renderer()
    # River and state labels already on the map are obstacles too.
    obstacles = [tx.get_window_extent(_r).transformed(inv) for tx in ax.texts]
    for ph, (cx, cy) in anchors.items():
        t = ax.text(cx, cy, ph, fontsize=8, fontweight="bold",
                    ha="center", va="center", zorder=13, color="black",
                    path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
        bb = t.get_window_extent(fig.canvas.get_renderer()).transformed(inv)
        hw, hh = bb.width / 2 + 1_000.0, bb.height / 2 + 1_000.0

        def search(poly):
            best = None
            for gx in np.arange(ext[0] + hw, ext[1] - hw, 1_000.0):
                for gy in np.arange(ext[2] + hh + 0.08 * (ext[3] - ext[2]),
                                    ext[3] - hh, 1_000.0):
                    if poly is not None and not poly.contains(box(gx - hw, gy - hh, gx + hw, gy + hh)):
                        continue
                    if np.any((np.abs(E - gx) < hw + 800.0) & (np.abs(Nm - gy) < hh + 800.0)):
                        continue
                    if any(gx - hw < o.x1 and gx + hw > o.x0 and gy - hh < o.y1 and gy + hh > o.y0
                           for o in obstacles):
                        continue
                    d = np.hypot(gx - cx, gy - cy)
                    if best is None or d < best[0]:
                        best = (d, gx, gy)
            return best

        best = search(inside.get(ph))
        if best is None and ph in spill:
            best = search(None)
        if best is None:
            raise SystemExit(f"no clear label position for phase {ph}")
        t.set_position((best[1], best[2]))
        obstacles.append(t.get_window_extent(_r).transformed(inv))


def add_scale_bar(ax, ext, x0=0.05, y0=0.955):
    """25 km scale bar at (x0, y0) in axes fractions."""
    bar = 25_000.0 / (ext[1] - ext[0])
    ax.plot([x0, x0 + bar], [y0, y0], transform=ax.transAxes, color="black", lw=2,
            zorder=15)
    ax.text(x0 + bar / 2, y0 + 0.008, "25 km", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=7, zorder=15,
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])


def move_label(ax, text, de=0.0, dn=0.0, remove=False):
    """Shift (or remove) a map label the basemap placed automatically.

    River names are put at the centroid of the river's geometry, which on some
    maps is on top of an assemblage or a boundary line. `de` and `dn` are in
    metres, east and north."""
    for t in list(ax.texts):
        if t.get_text() == text:
            if remove:
                t.remove()
            else:
                x, y = t.get_position()
                t.set_position((x + de, y + dn))


def add_locator_inset(fig, ax, ext, rect=(0.605, 0.02, 0.38, 0.30)):
    """North America locator inset. `rect` is (left, bottom, width, height) in
    fractions of the data axes; the position is taken from the finalized axes
    box so it tracks the equal-aspect map under bbox_inches."""
    fig.canvas.draw()
    pos = ax.get_position()
    inset_rect = [pos.x0 + rect[0] * pos.width, pos.y0 + rect[1] * pos.height,
                  rect[2] * pos.width, rect[3] * pos.height]
    mm._add_na_inset(fig, ext, rect=inset_rect, grayscale=True)


def main():
    counts_df, coords_df = m35.load_full()
    names = [str(x) for x in counts_df.index]
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    basin = mf._basin_members("curated")

    labels = assign_phases(names)
    lab_idx = np.array([PHASES.index(p) if p in PHASES else -1 for p in labels])

    E, Nm = to_utm(coords[:, 1], coords[:, 0])
    margin = 10_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)
    is_basin = np.array([nm in basin for nm in names])

    fig, ax = plt.subplots(figsize=(5.8, 6.8))
    # The assemblages fall in Arkansas, Tennessee, and Mississippi; Missouri
    # lies just north of the mapped area.
    draw_base(ax, ext, {"AR": (0.20, 0.46), "TN": (0.92, 0.74), "MS": (0.40, 0.16)})

    # Each phase gets a distinct gray (American Antiquity prints without color);
    # the territories do not overlap and the labels name them.
    territory = voronoi_territories(E, Nm, lab_idx, PHASES, ext)
    for ph, poly in territory.items():
        fill_polygon(ax, poly, facecolor=str(PHASE_GRAY[ph]), edgecolor="0.2",
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
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    anchors = {ph: (E[lab_idx == k].mean(), Nm[lab_idx == k].mean())
               for k, ph in enumerate(PHASES) if (lab_idx == k).any()}
    place_phase_labels(fig, ax, ext, E, Nm, anchors, territory)
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

    # 25 km scale bar in the upper-left corner, north of every territory. In the
    # lower-left it crossed the Parchman territory and the Mississippi, and the
    # lower-right is reserved for the locator inset, in the blank space east of
    # the Mississippi.
    add_scale_bar(ax, ext)
    add_locator_inset(fig, ax, ext)

    mf.save_all(fig, OUT_FIG)
    plt.close(fig)
    counts = {p: int((labels == p).sum()) for p in PHASES}
    basin_phases = {p: int(((labels == p) & is_basin).sum()) for p in PHASES}
    print(f"phase counts (all {len(names)} with coordinates):", counts)
    print(f"basin-set ({int(is_basin.sum())}) by phase:", {k: v for k, v in basin_phases.items() if v})
    print(f"wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
