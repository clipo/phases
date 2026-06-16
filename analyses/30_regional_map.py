"""30_regional_map.py — regional two-scale overview map (Figure 1A).

Sets the geographic frame for the two-scale study: the LMV St. Francis basin
(Parkin phase, Phillips-Ford-Griffin curated assemblages) and the CMV southeast-
Missouri survey of Williams (1954), the ~48 km gap between them, the named
physiographic / phase regions, Parkin, and Powers Fort (Powers phase).

NOTE on layers: the project's local GIS (LMVHydrology, major_rivers, state
polygons) stops at ~lat 36.6 N, but the CMV deposits run to 37.1 N, so the
northern SE-Missouri area has no local hydrography or state lines. We therefore
draw hydrography and state boundaries from Natural Earth (full extent),
overlaying the detailed local Mississippi / St. Francis geometry where it exists.
A North America locator inset is reused from make_map.

Powers-phase is shown only as Powers Fort (a point); no Powers-phase polygon is
held by the project, so no area is drawn.

Read-only on the manuscript. Writes figures/fig1_regional.png.

Usage: .venv/bin/python analyses/30_regional_map.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import cartopy.io.shapereader as shpreader
from pyproj import Transformer
from shapely.geometry import box

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
import make_figures as mf  # noqa: E402
mm = importlib.import_module("make_map")

UTM = "EPSG:26915"
OUT = ROOT / "figures" / "fig1_regional.png"
C_LMV = "0.20"
C_CMV = "0.55"
RIVER = "0.5"
LL2U = Transformer.from_crs("EPSG:4326", UTM, always_xy=True)
LON0, LON1, LAT0, LAT1 = -91.3, -88.8, 34.8, 37.4


def to_utm(lon, lat):
    return LL2U.transform(np.asarray(lon, float), np.asarray(lat, float))


def ne(resolution, category, name):
    """Natural Earth layer as a UTM GeoDataFrame."""
    fn = shpreader.natural_earth(resolution=resolution, category=category, name=name)
    return gpd.read_file(fn).to_crs(UTM)


def main():
    # --- assemblages first, to center and zoom the extent on the data ---
    cdf, codf = mf._load_curated()
    cc = pd.read_csv(ROOT / "data" / "processed" / "williams1954_cmv_coords.tsv", sep="\t")
    all_lon = np.concatenate([codf["Longitude"].to_numpy(float), cc["lon"].to_numpy(float)])
    all_lat = np.concatenate([codf["Latitude"].to_numpy(float), cc["lat"].to_numpy(float)])
    clon, clat = float(all_lon.mean()), float(all_lat.mean())
    HALF_LON, HALF_LAT = 1.75, 1.55          # zoomed-out half-widths (degrees)
    lon0, lon1 = clon - HALF_LON, clon + HALF_LON
    lat0, lat1 = clat - HALF_LAT, clat + HALF_LAT
    e0, n0 = LL2U.transform(lon0, lat0)
    e1, n1 = LL2U.transform(lon1, lat1)
    e_min, e_max = min(e0, e1), max(e0, e1)
    n_min, n_max = min(n0, n1), max(n0, n1)
    ext = (e_min, e_max, n_min, n_max)
    clip = box(e_min, n_min, e_max, n_max)

    lmv_e, lmv_n = to_utm(codf["Longitude"].to_numpy(float), codf["Latitude"].to_numpy(float))
    cmv_e, cmv_n = to_utm(cc["lon"].to_numpy(float), cc["lat"].to_numpy(float))
    par_e, par_n = to_utm(-90.548, 35.282)

    # --- layers ---
    ne_states = ne("10m", "cultural", "admin_1_states_provinces_lines").clip(clip)
    parts = []
    for nm in ("rivers_north_america", "rivers_lake_centerlines"):
        try:
            parts.append(ne("10m", "physical", nm).clip(clip))
        except Exception as exc:
            print(f"NE {nm} skipped:", exc)
    ne_rivers = (gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=UTM)
                 if parts else gpd.GeoDataFrame())
    # detailed local hydrology (LMV + most of the CMV up to ~36.6 N)
    loc_hydro = mm._load_hydrology().clip(clip)
    loc_rivers = mm._load_major_rivers().clip(clip)
    ms_poly = loc_rivers[loc_rivers["RIVER_NAME"] == "Mississippi River"]

    # --- figure ---
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8})
    fig, ax = plt.subplots(figsize=(6.2, 7.2))
    ax.set_facecolor("#fbfaf7")

    if not loc_hydro.empty:
        loc_hydro.plot(ax=ax, color=RIVER, linewidth=0.4, alpha=0.85, zorder=2)
    if not ne_rivers.empty:
        ne_rivers.plot(ax=ax, color=RIVER, linewidth=0.7, zorder=2)
    if not ms_poly.empty:
        ms_poly.plot(ax=ax, facecolor="#bcd6ea", edgecolor=RIVER, linewidth=0.3, zorder=2)
    if not ne_states.empty:
        ne_states.plot(ax=ax, color="#8a8a8a", linewidth=1.0, zorder=3)

    ax.scatter(cmv_e, cmv_n, s=22, c=C_CMV, edgecolor="white", linewidth=0.35,
               zorder=6, label="CMV assemblages (Williams 1954)")
    ax.scatter(lmv_e, lmv_n, s=22, c=C_LMV, edgecolor="white", linewidth=0.35,
               zorder=6, label="LMV St. Francis basin (Parkin phase)")
    ax.scatter([par_e], [par_n], marker="*", s=180, c="black", edgecolor="white",
               linewidth=0.6, zorder=7)

    def lab(lon, lat, text, color="#333333", fs=7.5, style="italic", weight="normal"):
        e, n = to_utm(lon, lat)
        ax.text(float(e), float(n), text, fontsize=fs, color=color, style=style,
                weight=weight, ha="center", va="center", zorder=9,
                path_effects=[pe.withStroke(linewidth=2.4, foreground="white")])

    lab(-90.45, 35.30, "Parkin", color="black", style="normal", weight="bold", fs=8)
    lab(-89.30, 36.92, "Cairo Lowland", color=C_CMV)
    lab(-89.62, 36.55, "Sikeston Ridge", color=C_CMV)
    lab(-90.18, 36.40, "Malden Plain", color=C_CMV)
    lab(-89.70, 36.10, "Little River", color=C_CMV)
    lab(-90.10, 35.15, "St. Francis basin", color=C_LMV)
    lab(-90.30, 34.95, "Mississippi R.", color=RIVER, fs=6.5)
    lab(-89.95, 36.85, "St. Francis R.", color=RIVER, fs=6.0)
    lab(-91.05, 35.6, "ARKANSAS", color="#777777", fs=8.5, style="normal")
    lab(-91.05, 36.9, "MISSOURI", color="#777777", fs=8.5, style="normal")
    lab(-89.05, 35.45, "TENNESSEE", color="#777777", fs=8.5, style="normal")

    ax.set_xlim(e_min, e_max); ax.set_ylim(n_min, n_max)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=6.5, frameon=True, framealpha=0.92)

    # scale bar (top-left, clear of the lower-left inset)
    bar_m = 50_000
    x0 = e_min + 0.06 * (e_max - e_min); y0 = n_max - 0.06 * (n_max - n_min)
    ax.plot([x0, x0 + bar_m], [y0, y0], color="black", lw=2, zorder=10)
    ax.text(x0 + bar_m / 2, y0 + 0.012 * (n_max - n_min), "50 km",
            ha="center", va="bottom", fontsize=7, zorder=10)

    # North America locator inset (reused from make_map)
    try:
        mm._add_na_inset(fig, ext)
    except Exception as exc:
        print("inset skipped:", exc)

    mf.save_all(fig, OUT)
    print(f"CMV n={len(cc)} lat {cc['lat'].min():.2f}-{cc['lat'].max():.2f}; LMV n={len(codf)}")
    print(f"local hydro={len(loc_hydro)} NE rivers={len(ne_rivers)} states={len(ne_states)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
