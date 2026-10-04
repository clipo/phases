#!/usr/bin/env python3
"""Mound height from lidar, measured the same way at every site.

The settlement ranking mixes provenances: Parkin's height is an excavator's
field measurement (Morse 1981, 1990), the other 54 come from a 2006 compilation
that records no per-site provenance, and at least one of those (13-N-3) has no
mound within 2 km of its recorded coordinate. This script replaces "the height
someone wrote down" with one procedure applied to one instrument.

THE PROCEDURE, stated once so every number here rides with it. Within a DEM
window (analyses/67_mound_lidar_tiles.py):

  1. Local relief = elevation minus a Gaussian-smoothed copy (sigma 40 m),
     which removes the regional gradient the mound sits on.
  2. The FOOTPRINT is the connected set of pixels with local relief above
     `edge_m`, containing the peak nearest the given coordinates. `edge_m` is
     the judgment this measurement cannot avoid: a mound grades into the field,
     so where it ends is a choice, not a reading.
  3. The BASELINE is a percentile of elevation in an annulus around that
     footprint, from `ring_in` to `ring_out` metres beyond its edge. A median
     resists the ditch on one side and the levee on the other.
  4. HEIGHT = summit elevation minus baseline.

Because step 2 and step 3 are choices, the script reports a sensitivity grid
rather than one number, and the caller quotes the operating point with the
value (rule 6). The spread across that grid IS the uncertainty; it is not a
confidence interval and must never be reported as one.

WHAT THIS MEASURES, and does not. It measures the mound as it stands today,
after ploughing, erosion, and in Parkin's case a state park's landscaping. At
Parkin it also runs into a definitional problem the field records share: Morse
(1990) reports the site surface itself standing 3 m above the floodplain from
occupational build-up, so a height measured against the floodplain includes
that platform, and one measured against the site surface does not. The two
answers differ by more than the correction this script was written to check,
which is why both are printed.

Usage:
    python analyses/68_mound_lidar_heights.py
    python analyses/68_mound_lidar_heights.py --site 11-N-1 --lat 35.2791 --lon -90.5582
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TILES = ROOT / "output" / "lidar"
OUT_MD = ROOT / "output" / "findings" / "mound_lidar_heights.md"

FT = 0.3048
LRM_SIGMA_M = 40.0
EDGE_GRID = (0.3, 0.5, 1.0)          # metres of local relief taken as the mound edge
RING = ((10, 40), (20, 60))          # metres beyond the footprint for the baseline
PCTL = (25, 50, 75)                  # baseline percentile in that annulus

# Mound coordinates where the site's recorded point does not sit on the mound.
# Both were located from the lidar and checked against a documentary source:
# Parkin against Morse (1990); 13-N-3 against the legal description
# SW1/4 NW1/4 S14 T3N R5E (Lee County), supplied by the author 2026-09-19.
KNOWN_MOUNDS = {
    "11-N-1": (35.27911, -90.55820),
    "13-N-3": (34.87130, -90.54478),
}


def _grids(dem_path: Path):
    from PIL import Image
    from scipy.ndimage import gaussian_filter
    z = np.array(Image.open(dem_path), dtype="float64")
    return z, z - gaussian_filter(z, sigma=LRM_SIGMA_M)


def _pixel_of(meta: dict, lat: float, lon: float, shape) -> tuple[int, int]:
    import math
    x0, y0, x1, y1 = meta["bbox_3857"]
    R = 6378137.0
    x = R * math.radians(lon)
    y = R * math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    px = (x - x0) / (x1 - x0) * shape[1]
    py = (y1 - y) / (y1 - y0) * shape[0]
    return int(round(py)), int(round(px))


def measure(z, lrm, seed_rc, mpp, edge_m, ring, pctl) -> dict | None:
    """Height of the mound containing `seed_rc`, under one set of choices."""
    from scipy.ndimage import label, binary_dilation

    lab, n = label(lrm > edge_m)
    r, c = seed_rc
    r = int(np.clip(r, 0, z.shape[0] - 1)); c = int(np.clip(c, 0, z.shape[1] - 1))
    comp = lab[r, c]
    if comp == 0:                      # the seed sits below the edge threshold:
        win = 30                       # take the nearest labelled peak instead
        r0, r1 = max(0, r - win), min(z.shape[0], r + win)
        c0, c1 = max(0, c - win), min(z.shape[1], c + win)
        near = lab[r0:r1, c0:c1]
        vals = near[near > 0]
        if vals.size == 0:
            return None
        comp = int(np.bincount(vals).argmax())
    foot = lab == comp

    grow_in = int(round(ring[0] / mpp)); grow_out = int(round(ring[1] / mpp))
    outer = binary_dilation(foot, iterations=grow_out)
    inner = binary_dilation(foot, iterations=grow_in)
    annulus = outer & ~inner
    if annulus.sum() < 50:
        return None
    baseline = float(np.percentile(z[annulus], pctl))
    summit = float(z[foot].max())
    return {"edge_m": edge_m, "ring": f"{ring[0]}-{ring[1]} m", "pctl": pctl,
            "footprint_m2": float(foot.sum()) * mpp * mpp,
            "summit_m": summit, "baseline_m": baseline,
            "height_m": summit - baseline, "height_ft": (summit - baseline) / FT}


def run(site: str, meta: dict, lat: float, lon: float) -> list[dict]:
    dem = TILES / meta["dem_tif"]
    z, lrm = _grids(dem)
    mpp = float(meta["ground_m_per_pixel"])
    seed = _pixel_of(meta, lat, lon, z.shape)
    rows = []
    for edge in EDGE_GRID:
        for ring in RING:
            for p in PCTL:
                m = measure(z, lrm, seed, mpp, edge, ring, p)
                if m:
                    m["site"] = site
                    rows.append(m)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default="")
    ap.add_argument("--lat", type=float)
    ap.add_argument("--lon", type=float)
    args = ap.parse_args()

    manifest = json.loads((TILES / "manifest.json").read_text())
    wanted = [args.site] if args.site else sorted(KNOWN_MOUNDS)
    all_rows: list[dict] = []
    lines = ["# Mound height from lidar", "",
             "One procedure, applied the same way at each site; see the header of",
             "`analyses/68_mound_lidar_heights.py` for its definition. The spread across",
             "the grid is what the choice of edge and baseline costs, not a confidence",
             "interval.", ""]

    for site in wanted:
        if site not in manifest:
            print(f"{site}: no tile in {TILES}/manifest.json", file=sys.stderr)
            continue
        meta = manifest[site]
        lat = args.lat if args.lat is not None else KNOWN_MOUNDS.get(
            site, (meta["center_lat"], meta["center_lon"]))[0]
        lon = args.lon if args.lon is not None else KNOWN_MOUNDS.get(
            site, (meta["center_lat"], meta["center_lon"]))[1]
        rows = run(site, meta, lat, lon)
        if not rows:
            print(f"{site}: no mound found at {lat}, {lon}", file=sys.stderr)
            continue
        all_rows += rows
        h = np.array([r["height_ft"] for r in rows])
        mid = [r for r in rows if r["edge_m"] == 0.5 and r["pctl"] == 50
               and r["ring"] == "10-40 m"][0]
        name = meta.get("site_name") or site
        lines += [f"## {site} ({name})", "",
                  f"- Recorded in the 2006 compilation: **{meta['recorded_height_ft']:.1f} ft**",
                  f"- Lidar, at the reference operating point (edge 0.5 m, baseline the "
                  f"median of the 10-40 m annulus): **{mid['height_ft']:.1f} ft "
                  f"({mid['height_m']:.2f} m)**, footprint {mid['footprint_m2']:.0f} m2",
                  f"- Across the {len(rows)}-cell sensitivity grid: {h.min():.1f} to "
                  f"{h.max():.1f} ft", "",
                  "| edge (m) | annulus | pctl | footprint (m2) | height (m) | height (ft) |",
                  "|---|---|---|---|---|---|"]
        for r in rows:
            lines.append(f"| {r['edge_m']:.1f} | {r['ring']} | {r['pctl']} | "
                         f"{r['footprint_m2']:.0f} | {r['height_m']:.2f} | "
                         f"{r['height_ft']:.1f} |")
        lines.append("")
        print(f"{site}: reference {mid['height_ft']:.1f} ft "
              f"({mid['height_m']:.2f} m); grid {h.min():.1f}-{h.max():.1f} ft; "
              f"recorded {meta['recorded_height_ft']:.1f} ft")

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nwrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
