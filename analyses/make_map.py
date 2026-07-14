"""make_map.py — geographic study-area map (Figure 1).

Produces figures/fig1_studyarea.png: site locations on local LMV GIS layers
showing real river geometry (St. Francis, Tyronza, Mississippi) and
valley geomorphology (LMVgeology3), scoped to the drainage-defined St. Francis
basin (within DRAINAGE_KM of the St. Francis / Tyronza / L'Anguille drainage,
matching analyses/16_basin_membership.py; Parkin-phase focus).

Map extent:
  The main map is bounded by the drainage-member site cluster, the basin
  definition used throughout the analysis. A 12 km margin is added around the
  basin site bounding box in UTM15N. This excludes the lower Yazoo / Winterville
  area and the Mississippi-proximal centers from the view.

Base-layer strategy:
  - LMVgeology3.shp  (POLYGON, EPSG:26915): valley structure background
  - LMVcounties4.shp (POLYGON, EPSG:26915): county outlines
  - LMVstates.shp    (POLYGON, EPSG:26915): state outlines + labels
  - LMVMajorRivers.shp (POLYGON, EPSG:26915): major rivers as filled polygons
  - LMVHydrology.shp (POLYLINE, geographic NAD27 ESRI:104000 -> reprojected
                       to EPSG:26915): named stream centerlines

Site coordinates:
  - Broad PFG set: drainage-basin sites only (within DRAINAGE_KM of the
    St. Francis / Tyronza / L'Anguille drainage), with UTM Zone 15N
    (Easting/Northing) from mls_emergence.dataio. Mound/ditch/fortification
    flags from LMVData-22March2006.csv. Marker size scaled by
    Max Mound Height (ft) where available (min size 6 for sites with no
    recorded mound height), making Parkin's 23 ft mound visually dominant.
  - Curated decorated set: drainage-basin sites from
    mainfort-pfg-cplXY.txt, converted to UTM15N.
  Parkin (11-N-1 / 'Parkin') marked as a large star + bold label; it is the
  within-basin primate center.

Parkin-phase delineation:
  - Sites with PFG Type == "St. Francis" from LMVData-22March2006.csv,
    restricted to the drainage basin (the same DRAINAGE_KM rule drops the
    southern outliers and the Mississippi-proximal centers). The retained
    cluster is delineated by a buffer-union: each point buffered by 9 km then
    union'd. This produces a compact outline that hugs the actual cluster.
  - St-Francis-type sites (retained set) are marked with filled triangles.

North America locator inset:
  - Small cartopy inset in the lower-left corner showing CONUS with state
    boundaries and a red rectangle indicating the NEW basin-scoped extent.

Data policy: this script does NOT print raw coordinates to stdout. It reads
gitignored data files but writes only to figures/ (also gitignored).

Usage:
    .venv/bin/python analyses/make_map.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
from scipy.spatial import cKDTree
from pyproj import Transformer
from shapely.geometry import box, MultiPoint, Point
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

from figstyle import save, OI_BLUE, OI_ORANGE, OI_GREEN, OI_VERMIL

warnings.filterwarnings("ignore")

from mls_emergence.dataio.pfg import load_pfg_counts
from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv, normalize_grid

DATA = ROOT / "data"
SHP = DATA / "Shapefiles"
UTM15N = "EPSG:26915"
PARKIN_BROAD = "11-N-1"
PARKIN_CUR = "Parkin"
# Hydrological basin definition (matches analyses/16_basin_membership.py):
# sites within DRAINAGE_KM of the St. Francis / Tyronza / L'Anguille drainage.
DRAINAGE_NAMES = {"Saint Francis River", "Tyronza River", "L'Anguille River",
                  "Saint Francis Floodway"}
DRAINAGE_KM = 20.0

# Geology grouping: DESCRIPT prefix/code -> class label
# H* = Holocene; P* = Pleistocene/older
# Three classes keep the background readable
_GEO_MEANDER = {  # Holocene meander belt, point bars, channel fills
    "Hpm1", "Hpm2", "Hpm3", "Hpm4", "Hpm5", "Hpm6",
    "Hchm", "Hcom", "Had", "Hal",
    "Hpa1", "Hpa2", "Hpa3", "Hpa4", "Hpa5", "Hpa6", "Hpa7",
    "Hpr1", "Hpr2", "Hpr3", "Hpr5", "Hpr6",
}
_GEO_BACKSWAMP = {  # Holocene backswamp / swamp / floodplain
    "HB", "Hpu", "Hps",
}
# All P* codes -> Pleistocene terrace/upland
_GEO_COLORS = {
    "meander": "#E8DFC8",    # warm tan
    "backswamp": "#D4E8D4",  # pale green
    "pleistocene": "#DCDCCC", # light gray-tan
}
_GEO_LABELS = {
    "meander": "Holocene meander belt & alluvium",
    "backswamp": "Holocene backswamp & floodplain",
    "pleistocene": "Pleistocene terrace & upland",
}


def _classify_geo(descript: str) -> str:
    if descript in _GEO_MEANDER:
        return "meander"
    if descript in _GEO_BACKSWAMP:
        return "backswamp"
    return "pleistocene"


# ---------------------------------------------------------------------------
# Load site data (UTM15N throughout)
# ---------------------------------------------------------------------------
def _load_broad_sites() -> pd.DataFrame:
    """Load broad PFG/LMV set; attach mound/ditch/fortification flags.

    Northing/Easting are already UTM Zone 15N (EPSG:26915).
    """
    broad_counts = load_pfg_counts(DATA / "raw" / "PFGData_sherds.csv")
    if not broad_counts.index.is_unique:
        broad_counts = broad_counts.groupby(level=0).sum()
    lmv = load_lmv(DATA / "LMVData_locations.csv")
    joined, _ = join_pfg_to_lmv(broad_counts, lmv)
    bm = joined.dropna(subset=["Easting", "Northing"]).copy()

    # Attach binary features from LMVData-22March2006.csv
    lmv2 = pd.read_csv(DATA / "LMVData-22March2006.csv")
    lmv2 = lmv2.dropna(subset=["Number"]).copy()
    lmv2["_k"] = lmv2["Number"].astype(str).map(normalize_grid)
    lmv2 = lmv2.drop_duplicates(subset=["_k"], keep="first").set_index("_k")
    norm_ids = pd.Index([normalize_grid(str(i)) for i in bm.index])
    ext = lmv2.reindex(norm_ids)
    ext.index = bm.index

    def to_bin(series: pd.Series) -> pd.Series:
        v = pd.to_numeric(
            series.astype(str).str.replace("?", "", regex=False), errors="coerce"
        )
        return v.fillna(0) > 0

    bm["mound"] = to_bin(ext["Mound"])
    bm["ditch"] = to_bin(ext["Ditch"])
    bm["stfr"] = to_bin(ext["St Francis"])
    bm["max_mound_ht"] = pd.to_numeric(
        ext["Max Mound Height (ft)"], errors="coerce"
    ).fillna(0.0)

    # Corrected Parkin record
    if PARKIN_BROAD in bm.index:
        bm.loc[PARKIN_BROAD, "mound"] = True
        bm.loc[PARKIN_BROAD, "ditch"] = True
        bm.loc[PARKIN_BROAD, "stfr"] = True

    # Add latitude column (UTM15N -> WGS84) for basin filtering
    t_to_ll = Transformer.from_crs(UTM15N, "EPSG:4326", always_xy=True)
    _lons, _lats = t_to_ll.transform(bm["Easting"].values, bm["Northing"].values)
    bm["lat"] = _lats

    return bm


def _load_curated_sites() -> pd.DataFrame:
    """Load curated decorated-set; convert lat/lon to UTM15N.

    The 'lat' column is retained on the returned DataFrame, but basin
    membership is applied by the caller via the drainage rule (_within_drainage).
    """
    xy = pd.read_csv(DATA / "raw" / "mainfort-pfg-cplXY.txt", sep="\t")
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    xy["lat"] = pd.to_numeric(xy["Latitude"], errors="coerce")
    xy["lon"] = pd.to_numeric(xy["Longitude"], errors="coerce")
    xy = xy.dropna(subset=["lat", "lon"])

    t = Transformer.from_crs("EPSG:4326", UTM15N, always_xy=True)
    easting, northing = t.transform(xy["lon"].values, xy["lat"].values)
    xy["Easting_utm"] = easting
    xy["Northing_utm"] = northing
    return xy


# ---------------------------------------------------------------------------
# GIS layer loaders — all return EPSG:26915
# ---------------------------------------------------------------------------
def _load_geology() -> gpd.GeoDataFrame:
    geo = gpd.read_file(SHP / "LMVgeology3.shp")
    # Already EPSG:26915
    geo["geo_class"] = geo["DESCRIPT"].fillna("").map(_classify_geo)
    return geo


def _load_hydrology() -> gpd.GeoDataFrame:
    """LMVHydrology: geographic NAD27 (ESRI:104000); set to EPSG:4267, reproject."""
    hydro = gpd.read_file(SHP / "LMVHydrology.shp")
    hydro = hydro.set_crs("EPSG:4267", allow_override=True)
    return hydro.to_crs(UTM15N)


def _load_major_rivers() -> gpd.GeoDataFrame:
    rivers = gpd.read_file(SHP / "LMVMajorRivers.shp")
    # Already EPSG:26915
    return rivers


def _load_states() -> gpd.GeoDataFrame:
    return gpd.read_file(SHP / "LMVstates.shp")  # EPSG:26915


def _load_counties() -> gpd.GeoDataFrame:
    return gpd.read_file(SHP / "LMVcounties4.shp")  # EPSG:26915


def _sf_drainage():
    """Union of the named St. Francis-system drainage features (EPSG:26915).

    Mirrors analyses/16_basin_membership.py: the basin is defined hydrologically
    as sites within DRAINAGE_KM of this geometry, which replaces the old latitude
    cut. Uses _load_hydrology so the corridor shares the reprojected CRS of the
    rivers drawn on the map.
    """
    hydro = _load_hydrology()
    return hydro[hydro["NAME"].isin(DRAINAGE_NAMES)].geometry.union_all()


def _within_drainage(df: pd.DataFrame, ekey: str, nkey: str, sf) -> np.ndarray:
    """Boolean mask: rows whose (easting, northing) lie within DRAINAGE_KM of sf."""
    pts = gpd.GeoSeries(gpd.points_from_xy(df[ekey], df[nkey]), crs=UTM15N)
    return (pts.distance(sf) / 1000.0 <= DRAINAGE_KM).to_numpy()


# ---------------------------------------------------------------------------
# Parkin-phase sites: PFG Type == "St. Francis" from LMVData-22March2006.csv
# ---------------------------------------------------------------------------
def _load_parkin_phase_sites() -> tuple[pd.DataFrame, int]:
    """Return all St-Francis-type sites; the caller applies drainage membership.

    All sites with PFG Type == 'St. Francis' are loaded (UTM Zone 15N,
    EPSG:26915), with a 'lat' column attached via pyproj (EPSG:26915 ->
    EPSG:4326). Drainage membership (within DRAINAGE_KM of the St. Francis
    system) is applied by make_map via _within_drainage, which drops the
    southern outliers and the Mississippi-proximal centers.

    Returns
    -------
    stfr_core : pd.DataFrame  (retained sites, with a 'lat' column added)
    n_all     : int           (total St-Francis-type count before the cut)
    """
    lmv2 = pd.read_csv(DATA / "LMVData-22March2006.csv")
    lmv2 = lmv2.dropna(subset=["Number"]).copy()
    stfr = lmv2[
        lmv2["PFG Type"].astype(str).str.contains(r"St\.\s*Francis", case=False, na=False)
    ].dropna(subset=["Easting", "Northing"]).copy()
    n_all = len(stfr)

    # Compute latitude for each site (UTM15N -> WGS84)
    t_to_ll = Transformer.from_crs(UTM15N, "EPSG:4326", always_xy=True)
    _lons, _lats = t_to_ll.transform(stfr["Easting"].values, stfr["Northing"].values)
    stfr = stfr.copy()
    stfr["lat"] = _lats

    # Drainage membership is applied by the caller (_within_drainage); return all.
    return stfr, n_all


# ---------------------------------------------------------------------------
# Study-area extent: basin site bbox + 12 km margin (in UTM15N meters)
#
# The extent is computed from the drainage-basin sites only, matching the
# empirical analysis scope. A 12 km margin keeps a small frame
# of geographic context (rivers, county lines) around the site cluster without
# adding unnecessary blank space southward toward the Yazoo / Winterville area.
# ---------------------------------------------------------------------------
def _study_extent(bm: pd.DataFrame, margin_m: float = 12_000) -> tuple[float, float, float, float]:
    """Compute map extent from the passed (drainage-filtered) sites + margin_m."""
    basin = bm
    e_min = basin["Easting"].min() - margin_m
    e_max = basin["Easting"].max() + margin_m
    n_min = basin["Northing"].min() - margin_m
    n_max = basin["Northing"].max() + margin_m
    return e_min, e_max, n_min, n_max


def _clip_gdf(gdf: gpd.GeoDataFrame, extent: tuple) -> gpd.GeoDataFrame:
    """Clip a GeoDataFrame to a (e_min, e_max, n_min, n_max) UTM bounding box."""
    e_min, e_max, n_min, n_max = extent
    clip_box = box(e_min, n_min, e_max, n_max)
    return gdf.clip(clip_box)


# ---------------------------------------------------------------------------
# Scale bar and north arrow (axes-fraction coordinates, UTM15N)
# ---------------------------------------------------------------------------
def _add_scale_bar(ax: plt.Axes, extent: tuple, bar_km: int = 50) -> None:
    """Draw a scale bar of bar_km kilometers.

    Placed at bottom-left of the axes in axes-fraction coordinates.
    In UTM15N, 1 m = 1 m so we convert km -> m directly.
    """
    e_min, e_max, n_min, n_max = extent
    width_m = e_max - e_min
    bar_m = bar_km * 1000
    frac = bar_m / width_m   # fraction of axes width

    x0, y0 = 0.05, 0.05
    ax.annotate(
        "",
        xy=(x0 + frac, y0),
        xytext=(x0, y0),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-", lw=2.0, color="black"),
    )
    # Tick ends
    for xf in (x0, x0 + frac):
        ax.annotate(
            "",
            xy=(xf, y0 + 0.012),
            xytext=(xf, y0 - 0.012),
            xycoords="axes fraction",
            arrowprops=dict(arrowstyle="-", lw=1.5, color="black"),
        )
    ax.text(
        x0 + frac / 2, y0 + 0.022,
        f"{bar_km} km",
        ha="center", va="bottom", fontsize=7,
        transform=ax.transAxes,
    )


def _add_north_arrow(ax: plt.Axes, x0: float = 0.93, y0: float = 0.10) -> None:
    ax.annotate(
        "N",
        xy=(x0, y0 + 0.055),
        xytext=(x0, y0),
        xycoords="axes fraction",
        ha="center", va="center",
        fontsize=9, fontweight="bold",
        arrowprops=dict(arrowstyle="-|>", lw=1.5, color="black"),
    )


# ---------------------------------------------------------------------------
# River label placement: centroid of clipped geometry
# ---------------------------------------------------------------------------
def _label_river(ax: plt.Axes, gdf: gpd.GeoDataFrame, text: str,
                 color: str, fontsize: float = 6.5, rotation: float = 0,
                 offset_e: float = 0, offset_n: float = 0) -> None:
    """Place a river label at the centroid of the union of gdf geometries."""
    if gdf.empty:
        return
    union = gdf.union_all() if hasattr(gdf, "union_all") else gdf.unary_union
    cx = union.centroid.x + offset_e
    cy = union.centroid.y + offset_n
    ax.text(
        cx, cy, text,
        fontsize=fontsize, color=color, rotation=rotation,
        ha="center", va="center",
        path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        zorder=12,
    )


# ---------------------------------------------------------------------------
# Reusable basin basemap (rivers + geology), for figures beyond the study map
# ---------------------------------------------------------------------------
def basin_basemap(ax: plt.Axes, extent: tuple, geology: bool = True,
                  grayscale: bool = False, show_counties: bool = True,
                  show_states: bool = True, draw_rivers: bool = True) -> None:
    """Draw the St. Francis basin river and geology basemap (UTM15N) onto ax,
    clipped to extent=(e_min, e_max, n_min, n_max). Sets equal aspect, the axis
    limits, and a thin frame; labels the Mississippi, St. Francis, and Tyronza.
    Reused by the emergence figure so it shares fig1's geographic base.

    geology=False skips the surficial-geology fill (which covers only the
    Arkansas/Missouri side) so the whole map reads as uniform land, leaving the
    county and state lines, which span all five states, to carry the base.
    grayscale=True renders the water and rivers in gray tones for figures that
    must print without color (American Antiquity does not publish color)."""
    geo = _clip_gdf(_load_geology(), extent) if geology else None
    hydro = _clip_gdf(_load_hydrology(), extent)
    major_rivers = _clip_gdf(_load_major_rivers(), extent)
    states = _clip_gdf(_load_states(), extent)
    counties = _clip_gdf(_load_counties(), extent)
    stfr_lines = hydro[hydro["NAME"] == "Saint Francis River"]
    tyronza_lines = hydro[hydro["NAME"] == "Tyronza River"]
    ms_poly = major_rivers[major_rivers["RIVER_NAME"] == "Mississippi River"]
    stfr_poly = major_rivers[major_rivers["RIVER_NAME"] == "St. Francis River"]

    if geology:
        for cls, color in _GEO_COLORS.items():
            subset = geo[geo["geo_class"] == cls]
            if not subset.empty:
                subset.plot(ax=ax, facecolor=color, edgecolor="none", linewidth=0,
                            zorder=0)
        if not geo.empty:
            geo.plot(ax=ax, facecolor="none", edgecolor="#C8C0B0", linewidth=0.15,
                     zorder=1)
    if show_counties and not counties.empty:
        counties.plot(ax=ax, facecolor="none", edgecolor="#AAAAAA", linewidth=0.4, zorder=2)
    if show_states and not states.empty:
        states.plot(ax=ax, facecolor="none", edgecolor="#888888", linewidth=1.0, zorder=3)

    if grayscale:
        water_blue, water_blue_dark, line_blue = "0.66", "0.25", "0.4"
        stream_col = "0.55"
    else:
        water_blue, water_blue_dark, line_blue = "#A8D0E8", "#4A8AB0", "#5BA3C9"
        stream_col = "#90C0D8"
    edge_col = "0.5" if grayscale else "#6EB5D8"
    # draw_rivers=False suppresses the local LMVHydrology lines so a caller can
    # supply a consistent hydrography layer (draw_hydrorivers) instead; the river
    # labels below are kept either way.
    if draw_rivers:
        if not ms_poly.empty:
            ms_poly.plot(ax=ax, facecolor=water_blue, edgecolor=edge_col, linewidth=0.4, zorder=4)
        if not stfr_poly.empty:
            stfr_poly.plot(ax=ax, facecolor=water_blue, edgecolor=edge_col, linewidth=0.3, zorder=4)
        if not hydro.empty:
            other = hydro[~hydro["NAME"].isin(["Saint Francis River", "Tyronza River"])
                          & (hydro["FEATURE"] == "Stream")]
            if not other.empty:
                other.plot(ax=ax, color=stream_col, linewidth=0.5, zorder=5)
        if not stfr_lines.empty:
            stfr_lines.plot(ax=ax, color=line_blue, linewidth=1.3, zorder=6)
        if not tyronza_lines.empty:
            tyronza_lines.plot(ax=ax, color=line_blue, linewidth=1.0, zorder=6)

    _label_river(ax, ms_poly, "Mississippi R.", water_blue_dark, rotation=90)
    _label_river(ax, stfr_poly, "St. Francis R.", water_blue_dark, rotation=70)
    _label_river(ax, tyronza_lines, "Tyronza R.", water_blue_dark, rotation=0)

    e_min, e_max, n_min, n_max = extent
    ax.set_xlim(e_min, e_max)
    ax.set_ylim(n_min, n_max)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.5)
        sp.set_edgecolor("0.6")


# ---------------------------------------------------------------------------
# Consistent HydroRIVERS hydrography (shared by Figures 1, 9, 10)
# ---------------------------------------------------------------------------
HYDRORIVERS = DATA / "HydroRIVERS_v10_na.gdb.zip"

# Line width (points) per HydroRIVERS flow order. Lower ORD_FLOW = larger river,
# so the Mississippi (order ~2) reads as the heaviest line and small tributaries
# taper. Used identically across the three regional maps for a uniform look.
_ORD_LW = {1: 2.4, 2: 2.0, 3: 1.3, 4: 0.9, 5: 0.6, 6: 0.42, 7: 0.30,
           8: 0.22, 9: 0.16, 10: 0.13}


def draw_hydrorivers(ax, extent, *, max_ord: int = 6, main_ord: int = 3,
                     grayscale: bool = False, lw_scale: float = 1.0,
                     zorder: float = 4.4) -> None:
    """Draw HydroRIVERS [Lehner & Grill 2013] hydrography onto ax, UTM15N.

    A single consistent hydrography layer for the regional maps: rivers down to
    flow order ``max_ord`` are drawn with line width scaled by flow order, so the
    Mississippi main stem is the heaviest line and small tributaries taper. Rivers
    of order <= ``main_ord`` (the Mississippi and the largest trunks) are drawn in
    blue; smaller tributaries in gray. ``extent`` is (e_min, e_max, n_min, n_max)
    in UTM15N; it is reprojected to a lat/lon bbox for a fast spatial read.
    """
    from shapely.geometry import box as _box
    t = Transformer.from_crs(UTM15N, "EPSG:4326", always_xy=True)
    xs = (extent[0], extent[1], extent[0], extent[1])
    ys = (extent[2], extent[2], extent[3], extent[3])
    lons, lats = t.transform(xs, ys)
    pad = 0.15  # degrees, so lines crossing the edge are not clipped early
    bbox = (min(lons) - pad, min(lats) - pad, max(lons) + pad, max(lats) + pad)

    riv = gpd.read_file(HYDRORIVERS, bbox=bbox)
    riv = riv[riv["ORD_FLOW"] <= max_ord]
    if riv.empty:
        return
    riv = riv.to_crs(UTM15N)
    riv = riv.clip(_box(extent[0], extent[2], extent[1], extent[3]))
    if riv.empty:
        return

    trib_col = "0.62"
    main_col = "0.32" if grayscale else "#5B9BC9"
    # draw smallest first so larger trunks sit on top
    for ordv in sorted(riv["ORD_FLOW"].unique(), reverse=True):
        sub = riv[riv["ORD_FLOW"] == ordv]
        lw = _ORD_LW.get(int(ordv), 0.15) * lw_scale
        is_main = ordv <= main_ord
        sub.plot(ax=ax, color=(main_col if is_main else trib_col),
                 linewidth=lw, zorder=zorder + (0.3 if is_main else 0.0),
                 capstyle="round")


# ---------------------------------------------------------------------------
# River-network (along-waterway) distance between sites
# ---------------------------------------------------------------------------
def river_distance_matrix(coords_latlon, tol: float = 250.0,
                          margin: float = 25_000.0):
    """Along-waterway shortest-path distance (km) between sites.

    coords_latlon: N x 2 array of (latitude, longitude). Builds a graph from the
    LMVHydrology centerlines clipped to the site bounding box + margin, repairs
    topology by bridging vertices within ``tol`` meters of each other, keeps the
    largest connected component, snaps each site to the nearest network node
    (adding the snap distance as an access cost at both ends), and returns
    (D_km[N,N], site_xy_utm[N,2], access_km[N], info_dict)."""
    arr = np.asarray(coords_latlon, float)
    s = gpd.GeoSeries(gpd.points_from_xy(arr[:, 1], arr[:, 0]),
                      crs="EPSG:4326").to_crs(UTM15N)
    SX, SY = s.x.to_numpy(), s.y.to_numpy()
    ext = (SX.min() - margin, SX.max() + margin,
           SY.min() - margin, SY.max() + margin)
    h = _clip_gdf(_load_hydrology(), ext)

    G = nx.Graph()

    def add_line(coords):
        for (ax, ay), (bx, by) in zip(coords[:-1], coords[1:]):
            a = (round(ax, 1), round(ay, 1))
            b = (round(bx, 1), round(by, 1))
            if a == b:
                continue
            w = float(np.hypot(a[0] - b[0], a[1] - b[1]))
            if G.has_edge(a, b):
                if w < G[a][b]["weight"]:
                    G[a][b]["weight"] = w
            else:
                G.add_edge(a, b, weight=w)

    for geom in h.geometry:
        if geom is None:
            continue
        if geom.geom_type == "LineString":
            add_line(list(geom.coords))
        elif geom.geom_type == "MultiLineString":
            for ln in geom.geoms:
                add_line(list(ln.coords))

    # repair topology: bridge vertices within tol meters across separate lines
    nodes = list(G.nodes)
    if len(nodes) > 1:
        xy = np.array(nodes, float)
        for i, j in cKDTree(xy).query_pairs(tol):
            a, b = nodes[i], nodes[j]
            if not G.has_edge(a, b):
                G.add_edge(a, b, weight=float(np.hypot(a[0] - b[0], a[1] - b[1])))

    comp = max(nx.connected_components(G), key=len)
    H = G.subgraph(comp).copy()
    hnodes = list(H.nodes)
    htree = cKDTree(np.array(hnodes, float))

    access = np.empty(len(SX))
    snap = []
    for k in range(len(SX)):
        d, idx = htree.query([SX[k], SY[k]])
        access[k] = d
        snap.append(hnodes[idx])

    N = len(SX)
    D = np.zeros((N, N))
    for i in range(N):
        L = nx.single_source_dijkstra_path_length(H, snap[i], weight="weight")
        for j in range(N):
            D[i, j] = L.get(snap[j], np.inf) + access[i] + access[j]
    np.fill_diagonal(D, 0.0)

    info = {
        "n_nodes": G.number_of_nodes(),
        "n_components": nx.number_connected_components(G),
        "largest_component": H.number_of_nodes(),
        "max_access_km": float(access.max() / 1000.0),
        "n_unreachable": int(np.isinf(D).sum()),
    }
    return D / 1000.0, np.column_stack([SX, SY]), access / 1000.0, info


# ---------------------------------------------------------------------------
# North America locator inset
# ---------------------------------------------------------------------------
def _add_na_inset(fig: plt.Figure, extent_utm: tuple, rect=None,
                  grayscale: bool = False) -> None:
    """Add a small North America locator inset.

    extent_utm: (e_min, e_max, n_min, n_max) in UTM15N meters.
    rect: figure-fraction [left, bottom, width, height] for the inset axes;
        defaults to the lower-left corner.
    grayscale: render land, water, and the study-area box in gray tones for
        figures that print without color.
    """
    e_min, e_max, n_min, n_max = extent_utm
    if rect is None:
        rect = [0.02, 0.02, 0.19, 0.26]

    if grayscale:
        land_c, ocean_c, river_c = "#DDDDDD", "#F0F0F0", "0.45"
        box_face, box_edge = "0.55", "0.15"
    else:
        land_c, ocean_c, river_c = "#E8E4DC", "#C8DCE8", "#2E6DA4"
        box_face, box_edge = "#CC000033", "red"

    # Convert UTM extent corners to geographic (lon/lat) for the study-area box
    t_to_ll = Transformer.from_crs(UTM15N, "EPSG:4326", always_xy=True)
    lon0, lat0 = t_to_ll.transform(e_min, n_min)
    lon1, lat1 = t_to_ll.transform(e_max, n_max)

    inset_ax = fig.add_axes(
        rect,
        projection=ccrs.AlbersEqualArea(central_longitude=-96, central_latitude=37.5),
    )

    # Focus on CONUS + southern Canada / northern Mexico
    inset_ax.set_extent([-130, -65, 22, 52], crs=ccrs.PlateCarree())

    # Natural Earth features
    inset_ax.add_feature(
        cfeature.NaturalEarthFeature("physical", "land", "110m",
                                     facecolor=land_c, edgecolor="none"),
        zorder=0,
    )
    inset_ax.add_feature(
        cfeature.NaturalEarthFeature("physical", "ocean", "110m",
                                     facecolor=ocean_c, edgecolor="none"),
        zorder=0,
    )
    inset_ax.add_feature(
        cfeature.NaturalEarthFeature("cultural", "admin_1_states_provinces_lines", "50m",
                                     facecolor="none", edgecolor="#AAAAAA", linewidth=0.4),
        zorder=1,
    )
    inset_ax.add_feature(
        cfeature.NaturalEarthFeature("cultural", "admin_0_countries", "110m",
                                     facecolor="none", edgecolor="#808080", linewidth=0.7),
        zorder=2,
    )
    inset_ax.coastlines(resolution="110m", linewidth=0.5, color="#505050", zorder=2)

    # Mississippi River (and major tributaries), so the study area reads as a
    # location on the river rather than an isolated point.
    inset_ax.add_feature(
        cfeature.NaturalEarthFeature("physical", "rivers_lake_centerlines", "50m",
                                     facecolor="none", edgecolor=river_c, linewidth=0.7),
        zorder=2,
    )
    inset_ax.text(-90.0, 41.5, "Mississippi R.", transform=ccrs.PlateCarree(),
                  fontsize=4.5, color=river_c, style="italic", rotation=62,
                  ha="center", va="center", zorder=4,
                  path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])

    # Rectangle for study-area extent
    study_box = box(lon0, lat0, lon1, lat1)
    xs, ys = study_box.exterior.xy
    inset_ax.fill(xs, ys, transform=ccrs.PlateCarree(),
                  facecolor=box_face, edgecolor=box_edge, linewidth=1.2, zorder=3)

    # Thin border around inset
    for spine in inset_ax.spines.values():
        spine.set_edgecolor("0.5")
        spine.set_linewidth(0.6)


# ---------------------------------------------------------------------------
# Main map
# ---------------------------------------------------------------------------
def make_map() -> None:
    # --- Load site data ---
    bm_all = _load_broad_sites()
    cur_all = _load_curated_sites()
    parkin_phase, n_parkin_all = _load_parkin_phase_sites()

    # Restrict broad, curated, and St-Francis-type sets to the drainage-defined
    # basin: sites within DRAINAGE_KM of the St. Francis / Tyronza / L'Anguille
    # drainage. This matches analyses/16_basin_membership.py, the analysis scope
    # throughout the paper, and replaces the earlier latitude cut.
    sf = _sf_drainage()
    bm = bm_all[_within_drainage(bm_all, "Easting", "Northing", sf)].copy()
    cur = cur_all[_within_drainage(cur_all, "Easting_utm", "Northing_utm", sf)].copy()
    parkin_phase = parkin_phase[_within_drainage(parkin_phase, "Easting", "Northing", sf)].copy()

    # --- Study-area extent: basin site bbox + 12 km margin (UTM15N) ---
    extent = _study_extent(bm, margin_m=12_000)
    e_min, e_max, n_min, n_max = extent

    # --- Load and clip GIS layers ---
    geo = _clip_gdf(_load_geology(), extent)
    hydro = _clip_gdf(_load_hydrology(), extent)
    major_rivers = _clip_gdf(_load_major_rivers(), extent)
    states = _clip_gdf(_load_states(), extent)
    counties = _clip_gdf(_load_counties(), extent)

    # Subset hydrology by name for labeling
    stfr_lines = hydro[hydro["NAME"] == "Saint Francis River"]
    tyronza_lines = hydro[hydro["NAME"] == "Tyronza River"]
    ms_poly = major_rivers[major_rivers["RIVER_NAME"] == "Mississippi River"]
    stfr_poly = major_rivers[major_rivers["RIVER_NAME"] == "St. Francis River"]

    # Parkin-phase delineation: ONE inclusive region around all drainage-member
    # St-Francis-type deposits. We take a (concave) hull of the deposit points and
    # pad it, so the shaded area is a single continuous polygon rather than the
    # separate blobs a per-point buffer-union produces. Sites more than DRAINAGE_KM
    # from the St. Francis system (southern outliers, Mississippi-proximal centers)
    # were already excluded by the drainage filter above.
    n_parkin_phase = len(parkin_phase)
    n_parkin_dropped = n_parkin_all - n_parkin_phase
    BUFFER_M = 7_000  # padding around the deposit hull
    parkin_pts = [
        Point(row["Easting"], row["Northing"])
        for _, row in parkin_phase.iterrows()
    ]
    _mp = MultiPoint(parkin_pts)
    try:
        from shapely import concave_hull as _concave_hull
        _hull = _concave_hull(_mp, ratio=0.35)
    except Exception:
        _hull = _mp.convex_hull
    parkin_region = _hull.buffer(BUFFER_M)

    # Verify exclusion: Winterville (~lat 33.48, id 19-L-1) must NOT be inside.
    lmv2_all = pd.read_csv(DATA / "LMVData-22March2006.csv")
    lmv2_all = lmv2_all.dropna(subset=["Number"]).copy()
    winterville_rows = lmv2_all[
        lmv2_all["Name"].astype(str).str.contains("Winterville", case=False, na=False)
    ].dropna(subset=["Easting", "Northing"])
    for _, wv in winterville_rows.iterrows():
        wv_pt = Point(wv["Easting"], wv["Northing"])
        inside = parkin_region.contains(wv_pt)
        print(f"Winterville inside Parkin-phase region: {inside}")
        assert not inside, (
            "Winterville is still inside the Parkin-phase region. "
            "Tighten the lat cutoff or buffer distance."
        )

    # Report (no raw coordinates)
    t_to_ll = Transformer.from_crs(UTM15N, "EPSG:4326", always_xy=True)
    _lon0, _lat0 = t_to_ll.transform(e_min, n_min)
    _lon1, _lat1 = t_to_ll.transform(e_max, n_max)
    print(f"MAP LAYER REPORT — drainage basin scope (<= {DRAINAGE_KM:.0f} km of drainage)")
    print(f"  Extent lat range:  {_lat0:.2f} to {_lat1:.2f}")
    print(f"  Extent lon range:  {_lon0:.2f} to {_lon1:.2f}")
    print(f"  Saint Francis River lines: {len(stfr_lines)} segments")
    print(f"  Tyronza River lines:       {len(tyronza_lines)} segments")
    print(f"  Mississippi River polygon: {len(ms_poly)} polygons")
    print(f"  St. Francis River polygon: {len(stfr_poly)} polygons")
    print(f"  Geology units in view:     {len(geo)} polygons")
    print(f"  Basin broad sites (drainage): {len(bm)}")
    print(f"  Basin curated sites (drainage): {len(cur)}")
    print(f"  St-Francis-type sites total: {n_parkin_all}, "
          f"within drainage: {n_parkin_phase}, "
          f"dropped (beyond drainage): {n_parkin_dropped}")
    print(f"  Parkin-phase extent: buffer-union, buffer = {BUFFER_M/1000:.0f} km per point")
    print("Parkin-phase extent = single inclusive hull around the drainage-member "
          "St-Francis-type deposits, approximate.")
    print("Winterville (~lat 33.5) is excluded from both the site set and the map extent.")

    # --- Set up figure ---
    fig_height = 7 * (n_max - n_min) / (e_max - e_min)
    fig, ax = plt.subplots(figsize=(7, fig_height))
    ax.set_xlim(e_min, e_max)
    ax.set_ylim(n_min, n_max)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_edgecolor("0.6")

    # -----------------------------------------------------------------------
    # Draw order
    # -----------------------------------------------------------------------

    # 1. Geology background
    for cls, color in _GEO_COLORS.items():
        subset = geo[geo["geo_class"] == cls]
        if not subset.empty:
            subset.plot(ax=ax, facecolor=color, edgecolor="none", linewidth=0, zorder=0)

    # Thin geology unit boundaries (very light)
    if not geo.empty:
        geo.plot(ax=ax, facecolor="none", edgecolor="#C8C0B0", linewidth=0.15, zorder=1)

    # 2. County outlines (thin gray)
    if not counties.empty:
        counties.plot(ax=ax, facecolor="none", edgecolor="#AAAAAA", linewidth=0.4, zorder=2)

    # 3. State outlines (heavier)
    if not states.empty:
        states.plot(ax=ax, facecolor="none", edgecolor="#888888", linewidth=1.0, zorder=3)

    # 4. Major rivers as filled blue polygons
    water_blue = "#A8D0E8"
    water_blue_dark = "#4A8AB0"
    if not ms_poly.empty:
        ms_poly.plot(ax=ax, facecolor=water_blue, edgecolor="#6EB5D8", linewidth=0.4, zorder=4)
    if not stfr_poly.empty:
        stfr_poly.plot(ax=ax, facecolor=water_blue, edgecolor="#6EB5D8", linewidth=0.3, zorder=4)

    # 5. Hydrology centerlines
    line_blue = "#5BA3C9"
    if not hydro.empty:
        # All other streams thin
        other_hydro = hydro[
            ~hydro["NAME"].isin(["Saint Francis River", "Tyronza River"])
            & (hydro["FEATURE"] == "Stream")
        ]
        if not other_hydro.empty:
            other_hydro.plot(ax=ax, color="#90C0D8", linewidth=0.5, zorder=5)

    # St. Francis: draw centerline over polygon fill
    if not stfr_lines.empty:
        stfr_lines.plot(ax=ax, color=line_blue, linewidth=1.3, zorder=6)

    # Tyronza: real geometry from LMVHydrology
    if not tyronza_lines.empty:
        tyronza_lines.plot(ax=ax, color=line_blue, linewidth=1.0, zorder=6)

    # 6. River labels
    # Mississippi label
    if not ms_poly.empty:
        _label_river(ax, ms_poly, "Mississippi R.", water_blue_dark,
                     fontsize=7.5, rotation=-20, offset_e=-5000, offset_n=10000)

    # St. Francis label
    if not stfr_lines.empty:
        _label_river(ax, stfr_lines, "St. Francis R.", water_blue_dark,
                     fontsize=7.0, rotation=15, offset_e=8000, offset_n=0)
    elif not stfr_poly.empty:
        _label_river(ax, stfr_poly, "St. Francis R.", water_blue_dark,
                     fontsize=7.0, rotation=15)

    # Tyronza label
    if not tyronza_lines.empty:
        _label_river(ax, tyronza_lines, "Tyronza R.", water_blue_dark,
                     fontsize=6.5, rotation=5, offset_e=-10000, offset_n=8000)

    # 7. State labels
    # Within the basin extent, only Arkansas has its centroid in view.
    # Missouri and Tennessee are visible at the northern / eastern map edges
    # respectively; place their abbreviations at fixed map-edge positions.
    # Missouri and Tennessee are placed by the edge-label block below (their
    # clipped centroids fall near the frame edge); the centroid loop handles only
    # states whose body is well inside the view, to avoid duplicate labels.
    state_label_positions = {
        "ARKANSAS": (-12000, -30000),
        "MISSISSIPPI": (-5000, -30000),
    }
    for _, row in states.iterrows():
        name = row["NAME"]
        if name not in state_label_positions:
            continue
        centroid = row.geometry.centroid
        dx, dy = state_label_positions[name]
        # Label states whose centroid falls in the main view
        if e_min < centroid.x < e_max and n_min < centroid.y < n_max:
            ax.text(
                centroid.x + dx, centroid.y + dy,
                name.title(),
                fontsize=9, color="0.45", alpha=0.85,
                ha="center", va="center", fontweight="bold", zorder=3,
            )

    # Edge labels: Missouri (top edge, left side) and Tennessee (right edge, upper)
    # placed at fixed axes-fraction positions since their centroids lie outside view.
    _edge_state_labels = [
        # (x_frac, y_frac, label, ha, va)
        (0.22, 0.975, "Missouri", "center", "top"),
        (0.985, 0.82, "Tennessee", "right", "center"),
    ]
    for xf, yf, lbl, ha, va in _edge_state_labels:
        ax.text(xf, yf, lbl, fontsize=9, color="0.45", alpha=0.85,
                ha=ha, va=va, fontweight="bold",
                transform=ax.transAxes, zorder=3)

    # 8. Parkin-phase delineation (drawn beneath site markers, above geology)
    # A single inclusive (concave-hull + padding) region covering all the
    # drainage-member St-Francis-type deposits, shaded as one continuous area.
    PARKIN_PHASE_FILL = "#CC79A7"   # reddish-purple (OI_PURPLE, colorblind-safe)

    def _plot_polygon(poly):
        """Fill and outline a single shapely polygon on ax."""
        xs, ys = poly.exterior.xy
        ax.fill(
            np.array(xs), np.array(ys),
            facecolor=PARKIN_PHASE_FILL, alpha=0.13, edgecolor="none",
            zorder=3,
        )
        ax.plot(
            np.array(xs), np.array(ys),
            linestyle="--", color=PARKIN_PHASE_FILL, linewidth=1.2, alpha=0.75,
            zorder=7,
        )

    # parkin_region may be a Polygon or MultiPolygon
    if hasattr(parkin_region, "geoms"):
        for _poly in parkin_region.geoms:
            _plot_polygon(_poly)
    else:
        _plot_polygon(parkin_region)

    # Label at region centroid
    hc = parkin_region.centroid
    ax.text(
        hc.x, hc.y, "Parkin phase",
        fontsize=8, color=PARKIN_PHASE_FILL, alpha=0.9,
        ha="center", va="center", fontweight="bold",
        path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        zorder=13,
    )

    # 9a. St-Francis-type sites: filled triangles (distinct phase-membership marker)
    # Draw on top of the hull fill, beneath other site symbols.
    ax.scatter(
        parkin_phase["Easting"].values, parkin_phase["Northing"].values,
        s=28, color=PARKIN_PHASE_FILL, alpha=0.9, marker="^",
        edgecolors="white", linewidths=0.4,
        zorder=9, label="St-Francis type (Parkin phase)",
    )

    # 9b. Site symbols (basin broad set), marker size scaled by Max Mound Height (ft).
    # Scale: s = S_MIN + (ht / HT_REF) * S_SCALE, so that the max height site
    # (Parkin, 23 ft) receives a large marker and sites with no recorded mound
    # height receive the minimum. The star for Parkin (step 10) is plotted on top.
    S_MIN = 6.0        # minimum scatter marker area (pt^2) — no-mound sites
    HT_REF = 23.0      # reference height = Parkin / max in dataset (ft)
    S_SCALE = 130.0    # additional area at HT_REF

    def _mound_size(ht_series: pd.Series) -> np.ndarray:
        """Map Max Mound Height (ft) -> scatter marker area (pt^2)."""
        ht = np.maximum(ht_series.values, 0.0)
        return S_MIN + (ht / HT_REF) * S_SCALE

    no_mound = bm[~bm["mound"]]
    mound_only = bm[bm["mound"] & ~bm["ditch"] & ~bm["stfr"]]
    fortified = bm[bm["ditch"] | bm["stfr"]]

    if len(no_mound) > 0:
        ax.scatter(
            no_mound["Easting"], no_mound["Northing"],
            s=_mound_size(no_mound["max_mound_ht"]),
            color="0.55", alpha=0.5, marker="o", edgecolors="none",
            zorder=7, label="No mound",
        )
    if len(mound_only) > 0:
        ax.scatter(
            mound_only["Easting"], mound_only["Northing"],
            s=_mound_size(mound_only["max_mound_ht"]),
            color=OI_BLUE, alpha=0.8, marker="o",
            edgecolors="white", linewidths=0.3,
            zorder=8, label="Mound present",
        )
    if len(fortified) > 0:
        ax.scatter(
            fortified["Easting"], fortified["Northing"],
            s=_mound_size(fortified["max_mound_ht"]),
            color=OI_ORANGE, alpha=0.9, marker="D",
            edgecolors="white", linewidths=0.4,
            zorder=9, label="Mound + ditch/fortified",
        )

    # Mound-height size legend (small inset key: 0 ft / 10 ft / 23 ft)
    _size_ht_labels = [0, 10, 23]
    _size_handles = [
        mlines.Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="0.55",
            markeredgecolor="0.55",
            markersize=np.sqrt(S_MIN + (ht / HT_REF) * S_SCALE),
            label=f"{ht} ft" if ht > 0 else "0 ft (no height recorded)",
        )
        for ht in _size_ht_labels
    ]

    # 9. Curated decorated set (open green circles; basin only)
    cur_not_parkin = cur[cur.index != PARKIN_CUR]
    if len(cur_not_parkin) > 0:
        ax.scatter(
            cur_not_parkin["Easting_utm"], cur_not_parkin["Northing_utm"],
            s=18, facecolors="none", edgecolors=OI_GREEN,
            linewidths=0.9, marker="o",
            zorder=10, label="Curated decorated set",
        )

    # 10. Parkin: large star + bold label — within-basin primate center.
    # Star size (s=400) is substantially larger than all other markers to
    # signal Parkin's rank-1 position within the basin.
    if PARKIN_CUR in cur.index:
        pk = cur.loc[PARKIN_CUR]
        ax.scatter(
            [float(pk["Easting_utm"])], [float(pk["Northing_utm"])],
            s=400, color=OI_VERMIL, marker="*",
            edgecolors="white", linewidths=1.0,
            zorder=11, label="Parkin (11-N-1)",
        )
        ax.annotate(
            "Parkin",
            (float(pk["Easting_utm"]), float(pk["Northing_utm"])),
            xytext=(10, 10), textcoords="offset points",
            fontsize=9, fontweight="bold", color=OI_VERMIL,
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
            zorder=12,
        )

    # 11. Scale bar (25 km for basin extent) + north arrow
    _add_scale_bar(ax, extent, bar_km=25)
    _add_north_arrow(ax, x0=0.93, y0=0.08)

    # 13. Legend
    geo_patches = [
        mpatches.Patch(facecolor=_GEO_COLORS["meander"], edgecolor="#C0B090",
                       linewidth=0.5, label=_GEO_LABELS["meander"]),
        mpatches.Patch(facecolor=_GEO_COLORS["backswamp"], edgecolor="#A0C0A0",
                       linewidth=0.5, label=_GEO_LABELS["backswamp"]),
        mpatches.Patch(facecolor=_GEO_COLORS["pleistocene"], edgecolor="#B0B0A0",
                       linewidth=0.5, label=_GEO_LABELS["pleistocene"]),
    ]
    site_handles = [
        mlines.Line2D([0], [0], marker="o", color="w", markerfacecolor="0.55",
                      markeredgecolor="none", markersize=5, label="No mound"),
        mlines.Line2D([0], [0], marker="o", color="w", markerfacecolor=OI_BLUE,
                      markeredgecolor="white", markersize=7, label="Mound present"),
        mlines.Line2D([0], [0], marker="D", color="w", markerfacecolor=OI_ORANGE,
                      markeredgecolor="white", markersize=7, label="Mound + ditch/fortified"),
        mlines.Line2D([0], [0], marker="^", color="w", markerfacecolor=PARKIN_PHASE_FILL,
                      markeredgecolor="white", markersize=7,
                      label=f"St-Francis type ({n_parkin_phase} sites, drainage basin)"),
        mlines.Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                      markeredgecolor=OI_GREEN, markersize=7, markeredgewidth=0.9,
                      label="Curated decorated set"),
        mlines.Line2D([0], [0], marker="*", color="w", markerfacecolor=OI_VERMIL,
                      markeredgecolor="white", markersize=14, label="Parkin (11-N-1)"),
    ]
    size_legend_handles = [
        mlines.Line2D([], [], linestyle="none", label="Marker size = max mound ht.:"),
    ] + _size_handles
    river_handles = [
        mpatches.Patch(facecolor=water_blue, edgecolor="#6EB5D8",
                       linewidth=0.5, label="Major rivers (polygon)"),
        mlines.Line2D([0], [0], color=line_blue, linewidth=1.3, label="St. Francis R. (centerline)"),
        mlines.Line2D([0], [0], color=line_blue, linewidth=1.0, label="Tyronza R. (centerline)"),
    ]
    parkin_phase_handle = [
        mpatches.Patch(facecolor=PARKIN_PHASE_FILL, alpha=0.20, edgecolor=PARKIN_PHASE_FILL,
                       linestyle="--", linewidth=1.0, label="Parkin-phase extent (approx.)"),
    ]
    all_handles = (geo_patches + parkin_phase_handle + site_handles
                   + size_legend_handles + river_handles)
    ax.legend(
        handles=all_handles, loc="lower right",
        fontsize=5.5, frameon=True, framealpha=0.9,
        edgecolor="0.7", borderpad=0.6, handlelength=1.5,
        ncol=1,
    )

    # 14. North America locator inset (lower-left corner)
    _add_na_inset(fig, extent)

    save(fig, "fig1_studyarea")
    print("fig1_studyarea.png written")


if __name__ == "__main__":
    make_map()
