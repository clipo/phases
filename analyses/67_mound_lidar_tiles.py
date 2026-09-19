#!/usr/bin/env python3
"""Lidar relief tiles for basin mound sites, for hand delineation.

The settlement cross-check currently ranks sites by `Max Mound Height (ft)` in
a 2006 compilation whose provenance varies by site, and whose mound-area field
is zero for Parkin. Lidar gives one instrument applied the same way everywhere.
Mound *area*, though, has no instrument answer: a mound grades into the field,
so the edge is a judgment about where the constructed surface ends. This script
prepares the evidence for that judgment; a person draws the outlines in the
companion web tool, and `68_mound_areas.py` measures what they drew against the
elevation model.

Source: the USGS 3DEP elevation service, which serves the 1 m DEM built from
the state's lidar. Eastern Arkansas has two collections over the St. Francis
basin, AR_Eastern_D23 (2026) and MO_AR_2014_Lidar (2014); the service returns
the current best available, so a tile is dated by `fetched_utc` in the
manifest rather than assumed.

Each tile is a square window centered on the site's recorded coordinates. The
window is deliberately wide (default 1 km): the recorded coordinates are
survey-file centroids and can sit hundreds of metres off the feature. At
Parkin the recorded point is about 450 m south of the mound.

Two renderings are written per site, because neither alone shows an edge well:

  hillshade      relief lit from the northwest, which shows slope breaks
  local relief   elevation minus a Gaussian-smoothed copy of itself, which
                 removes the regional gradient (here the river's natural
                 levee) and leaves features at mound scale

Projection note. The service returns EPSG:3857, whose pixels are inflated by
1/cos(latitude), about 1.2249 here. A pixel is therefore 1.2249 Mercator units
and very close to 1.0 ground metre. The manifest records both so that an area
computed from pixels is convertible without re-deriving the factor.

Not in MANIFEST.md, deliberately. The tiles are written under a directory given
at run time (`--out`), so `build_manifest.py`, which resolves write paths through
module-level constants, cannot attribute them and does not list this script. They
are also untracked: `output/` is gitignored apart from `findings/`, and a set of
1 m relief renderings runs to about 1.6 MB per site. Rerun this script to rebuild
them; the 3DEP service is the source of record, not the repository.

Usage:
    python analyses/67_mound_lidar_tiles.py --sites parkin --out output/lidar
    python analyses/67_mound_lidar_tiles.py --top 12
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

SERVICE = ("https://elevation.nationalmap.gov/arcgis/rest/services/"
           "3DEPElevation/ImageServer/exportImage")
R_EARTH = 6378137.0
OUT_DEFAULT = ROOT / "output" / "lidar"
LRM_SIGMA_M = 40.0      # regional trend removed at this scale
LRM_CLIP_M = 1.5        # colour saturation, so a 6 m mound does not flatten a 1 m one


def webmerc(lat: float, lon: float) -> tuple[float, float]:
    x = R_EARTH * math.radians(lon)
    y = R_EARTH * math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    return x, y


def inv_webmerc(x: float, y: float) -> tuple[float, float]:
    lon = math.degrees(x / R_EARTH)
    lat = math.degrees(2 * math.atan(math.exp(y / R_EARTH)) - math.pi / 2)
    return lat, lon


def fetch_dem(lat: float, lon: float, span_m: float, dest: Path) -> dict:
    """A square DEM window of span_m ground metres, centered on the site."""
    x, y = webmerc(lat, lon)
    k = math.cosh(y / R_EARTH)           # Mercator inflation at this latitude
    half = span_m * k / 2
    bbox = f"{x - half},{y - half},{x + half},{y + half}"
    px = int(round(span_m))              # 1 m ground sampling
    url = (f"{SERVICE}?bbox={bbox}&bboxSR=3857&imageSR=3857&size={px},{px}"
           f"&format=tiff&pixelType=F32&interpolation=RSP_BilinearInterpolation"
           f"&f=image")
    with urllib.request.urlopen(url, timeout=180) as r:
        data = r.read()
    if len(data) < 10_000:
        raise RuntimeError(f"service returned {len(data)} bytes, not a DEM")
    dest.write_bytes(data)
    return {"bbox_3857": [x - half, y - half, x + half, y + half],
            "mercator_units_per_pixel": 2 * half / px,
            "ground_m_per_pixel": span_m / px,
            "mercator_inflation": k,
            "size_px": px,
            "center_lat": lat, "center_lon": lon,
            "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}


def render(dem: Path, out_stub: Path) -> dict:
    """Hillshade and local-relief PNGs. Returns relief statistics."""
    from PIL import Image
    from scipy.ndimage import gaussian_filter
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.cm as cm

    hs_tif = out_stub.with_suffix(".hs.tif")
    subprocess.run(["gdaldem", "hillshade", "-z", "4", "-az", "315", "-alt", "45",
                    "-compute_edges", "-q", str(dem), str(hs_tif)], check=True)
    z = np.array(Image.open(dem), dtype="float64")
    hs = np.array(Image.open(hs_tif), dtype="float64")
    lrm = z - gaussian_filter(z, sigma=LRM_SIGMA_M)

    v = np.clip(lrm, -LRM_CLIP_M, LRM_CLIP_M)
    v = (v + LRM_CLIP_M) / (2 * LRM_CLIP_M)
    rgb = cm.RdYlBu_r(v)[..., :3] * 255
    shade = np.clip(hs / 255.0, 0, 1)[..., None]
    Image.fromarray((rgb * (0.45 + 0.55 * shade)).astype("uint8")).save(
        out_stub.with_suffix(".lrm.png"), optimize=True)
    Image.fromarray(np.clip(hs, 0, 255).astype("uint8")).save(
        out_stub.with_suffix(".hs.png"), optimize=True)
    hs_tif.unlink(missing_ok=True)
    return {"elev_min_m": float(z.min()), "elev_max_m": float(z.max()),
            "lrm_min_m": float(lrm.min()), "lrm_max_m": float(lrm.max()),
            "lrm_sigma_m": LRM_SIGMA_M, "lrm_clip_m": LRM_CLIP_M}


def basin_sites() -> pd.DataFrame:
    """Basin sites with a recorded mound height, with lat/lon."""
    import importlib
    from pyproj import Transformer
    from mls_emergence.dataio.pfg import load_pfg_counts
    from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv

    b17 = importlib.import_module("17_basin_results")
    mf = importlib.import_module("make_figures")
    broad = load_pfg_counts(ROOT / "data" / "raw" / "PFGData.xlsx")
    if not broad.index.is_unique:
        broad = broad.groupby(level=0).sum()
    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    joined, _ = join_pfg_to_lmv(broad, lmv)
    bm = joined.dropna(subset=["Easting", "Northing"]).copy()
    bm = bm[[str(i) in mf._basin_members("broad") for i in bm.index]]

    lmv2 = pd.read_excel(ROOT / "data" / "LMVData-22March2006.xls",
                         sheet_name="Sheet1").dropna(subset=["Number"])
    lmv2["_k"] = lmv2["Number"].astype(str).map(b17.normalize_grid)
    lmv2 = lmv2.drop_duplicates("_k").set_index("_k")
    ext = lmv2.reindex(pd.Index([b17.normalize_grid(str(i)) for i in bm.index]))
    ext.index = bm.index
    bm["height_ft"] = pd.to_numeric(ext["Max Mound Height (ft)"], errors="coerce")
    bm["n_mounds"] = pd.to_numeric(ext["Num_Mounds"], errors="coerce")
    bm["site_name"] = ext["Site Name"] if "Site Name" in ext.columns else ""
    bm = bm[bm["height_ft"] > 0].copy()

    lat, lon = [], []
    for zone, e, n in zip(bm["Zone"], bm["Easting"], bm["Northing"]):
        epsg = 32600 + int(zone) if not pd.isna(zone) else 32615
        tr = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
        x, y = tr.transform(float(e), float(n))
        lat.append(y), lon.append(x)
    bm["lat"], bm["lon"] = lat, lon
    return bm.sort_values("height_ft", ascending=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    ap.add_argument("--span", type=float, default=1000.0, help="window, ground m")
    ap.add_argument("--top", type=int, default=0, help="N tallest recorded mounds")
    ap.add_argument("--sites", default="", help="comma-separated site numbers")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    sites = basin_sites()
    if args.sites:
        want = {s.strip().lower() for s in args.sites.split(",")}
        sites = sites[[str(i).lower() in want or str(nm).strip().lower() in want
                       for i, nm in zip(sites.index, sites["site_name"])]]
    elif args.top:
        sites = sites.head(args.top)

    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    for sid, row in sites.iterrows():
        key = str(sid).replace("/", "_")
        stub = out / key
        try:
            meta = fetch_dem(row["lat"], row["lon"], args.span,
                             stub.with_suffix(".dem.tif"))
            meta.update(render(stub.with_suffix(".dem.tif"), stub))
        except Exception as exc:                       # noqa: BLE001
            print(f"{key}: FAILED ({exc})", file=sys.stderr)
            continue
        meta.update({"site_number": str(sid),
                     "site_name": str(row.get("site_name", "")).strip(),
                     "recorded_height_ft": float(row["height_ft"]),
                     "recorded_n_mounds": (None if pd.isna(row["n_mounds"])
                                           else float(row["n_mounds"])),
                     "span_m": args.span,
                     "lrm_png": f"{key}.lrm.png", "hs_png": f"{key}.hs.png",
                     "dem_tif": f"{key}.dem.tif"})
        manifest[str(sid)] = meta
        print(f"{key}: {meta['size_px']}px, relief {meta['lrm_min_m']:+.2f} to "
              f"{meta['lrm_max_m']:+.2f} m, recorded {row['height_ft']:.0f} ft")
        time.sleep(0.5)                                # be kind to the service
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(f"\nwrote {manifest_path} ({len(manifest)} sites)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
