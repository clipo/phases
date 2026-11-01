#!/usr/bin/env python3
"""Measure hand-drawn mound outlines against the lidar elevation model.

The outlines come from the delineation tool, where a person traced each mound's
edge on the 3DEP relief rendering (`analyses/67_mound_lidar_tiles.py`). Area is
theirs, computed from the polygon; this script adds the quantities that need
the elevation model, measured the same way at every site so the results are
comparable across the basin rather than across whoever recorded each entry in
the 2006 compilation.

Per outline:
  area        the polygon's area on the ground plane (theirs, recomputed here
              from the vertices as a check on the tool's arithmetic)
  summit      the highest elevation inside the polygon
  baseline    the median elevation in a ring from the polygon's edge outward,
              `RING_IN` to `RING_OUT` metres
  height      summit minus baseline
  volume      the elevation above that baseline, summed over the polygon

THE BASELINE IS A CHOICE, not a reading. A mound grades into the field, and the
ring's placement decides how much of the apron counts as ground. The script
therefore reports each height under three rings and prints the spread; a single
number without its ring is not a measurement. The spread is what the choice
costs; it is not a confidence interval and must not be quoted as one.

WHAT IT CANNOT SEE. The DEM is the surface as it stands: ploughed, eroded, and
at Parkin landscaped as a state park. At Parkin it also meets a definitional
problem the field records share, since the site surface itself stands about 3 m
above the floodplain from occupational build-up (Morse 1990), so a height
measured from the surrounding field includes that platform while an excavator's
figure for the mound alone does not.

Input:  output/lidar/manifest.json, the tiles it names, and the outlines saved
        by the tool (a directory of per-site JSON documents).
Output: output/findings/mound_outline_measures.md and a CSV beside it.

Usage:
    python analyses/69_mound_outline_measures.py --outlines <dir>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TILES = ROOT / "output" / "lidar"
OUT_MD = ROOT / "output" / "findings" / "mound_outline_measures.md"
OUT_CSV = ROOT / "output" / "findings" / "mound_outline_measures.csv"

FT = 0.3048
RINGS = ((10, 40), (20, 60), (5, 25))     # metres beyond the outline
REFERENCE_RING = (10, 40)                 # the one quoted as "the" height


def polygon_mask(pts, shape) -> np.ndarray:
    """Boolean mask of the polygon, in tile pixel coordinates."""
    from PIL import Image, ImageDraw
    img = Image.new("1", (shape[1], shape[0]), 0)
    ImageDraw.Draw(img).polygon([(float(x), float(y)) for x, y in pts], fill=1)
    return np.array(img, dtype=bool)


def shoelace_m2(pts, mpp: float) -> float:
    a = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % len(pts)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2 * mpp * mpp


def measure(z: np.ndarray, mask: np.ndarray, mpp: float, ring) -> dict:
    from scipy.ndimage import binary_dilation
    if mask.sum() == 0:
        return {}
    inner = binary_dilation(mask, iterations=int(round(ring[0] / mpp)))
    outer = binary_dilation(mask, iterations=int(round(ring[1] / mpp)))
    annulus = outer & ~inner
    if annulus.sum() < 30:
        return {}
    baseline = float(np.median(z[annulus]))
    summit = float(z[mask].max())
    above = np.clip(z[mask] - baseline, 0, None)
    return {"baseline_m": baseline, "summit_m": summit,
            "height_m": summit - baseline,
            "volume_m3": float(above.sum()) * mpp * mpp,
            "mean_rise_m": float(above.mean())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outlines", required=True, help="directory of per-site JSON")
    args = ap.parse_args()

    from PIL import Image
    manifest = json.loads((TILES / "manifest.json").read_text())
    docs = sorted(Path(args.outlines).glob("*.json"))
    if not docs:
        print(f"no outline documents under {args.outlines}", file=sys.stderr)
        return 2

    rows = []
    for doc in docs:
        d = json.loads(doc.read_text())
        d = d.get("data", d)
        site = str(d.get("site_number") or doc.stem)
        polys = d.get("polygons") or []
        if site not in manifest:
            print(f"{site}: no tile in the manifest, skipped", file=sys.stderr)
            continue
        meta = manifest[site]
        dem = TILES / meta["dem_tif"]
        if not dem.exists():
            print(f"{site}: {dem.name} missing, skipped", file=sys.stderr)
            continue
        z = np.array(Image.open(dem), dtype="float64")
        mpp = float(meta["ground_m_per_pixel"])
        for poly in polys:
            pts = poly.get("pts") or []
            if len(pts) < 3:
                continue
            mask = polygon_mask(pts, z.shape)
            base = {"site": site, "site_name": meta.get("site_name", ""),
                    "outline": poly.get("name", ""),
                    "recorded_height_ft": meta.get("recorded_height_ft"),
                    "area_m2": shoelace_m2(pts, mpp),
                    "area_tool_m2": poly.get("area_m2"),
                    "vertices": len(pts)}
            heights = {}
            for ring in RINGS:
                m = measure(z, mask, mpp, ring)
                if not m:
                    continue
                heights[ring] = m["height_m"]
                if ring == REFERENCE_RING:
                    base.update(m)
            if heights:
                base["height_spread_m"] = max(heights.values()) - min(heights.values())
                base["height_min_m"] = min(heights.values())
                base["height_max_m"] = max(heights.values())
            rows.append(base)

    if not rows:
        print("no measurable outlines", file=sys.stderr)
        return 2

    import pandas as pd
    df = pd.DataFrame(rows)
    df["height_ft"] = df["height_m"] / FT
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    per_site = df.groupby("site").agg(
        name=("site_name", "first"), outlines=("outline", "count"),
        recorded_ft=("recorded_height_ft", "first"),
        tallest_ft=("height_ft", "max"), total_area_m2=("area_m2", "sum"),
        total_volume_m3=("volume_m3", "sum")).sort_values("tallest_ft", ascending=False)

    L = ["# Mound outlines measured against the lidar", "",
         "Outlines drawn by hand on the 3DEP relief renderings; heights and volumes",
         "measured from the 1 m DEM. The quoted height uses a baseline of the median",
         "elevation in a 10-40 m ring outside each outline. `height_spread_m` is how far",
         "that height moves across three ring choices (10-40, 20-60, 5-25 m) and is the",
         "cost of the baseline decision, not a confidence interval.", "",
         f"{len(df)} outlines across {df['site'].nunique()} sites.", "",
         "## Per site", "",
         "| site | name | outlines | recorded ft | tallest measured ft | total area m2 | total volume m3 |",
         "|---|---|---|---|---|---|---|"]
    for site, r in per_site.iterrows():
        rec = "-" if pd.isna(r["recorded_ft"]) else f"{r['recorded_ft']:.0f}"
        L.append(f"| {site} | {r['name']} | {r['outlines']:.0f} | {rec} | "
                 f"{r['tallest_ft']:.1f} | {r['total_area_m2']:.0f} | "
                 f"{r['total_volume_m3']:.0f} |")
    L += ["", "## Every outline", "",
          "| site | outline | area m2 | height m | height ft | spread m | volume m3 |",
          "|---|---|---|---|---|---|---|"]
    for _, r in df.sort_values(["site", "outline"]).iterrows():
        L.append(f"| {r['site']} | {r['outline']} | {r['area_m2']:.0f} | "
                 f"{r.get('height_m', float('nan')):.2f} | {r['height_ft']:.1f} | "
                 f"{r.get('height_spread_m', float('nan')):.2f} | "
                 f"{r.get('volume_m3', float('nan')):.0f} |")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    print(per_site.to_string())
    print(f"\n{len(df)} outlines measured; wrote {OUT_MD.name} and {OUT_CSV.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
