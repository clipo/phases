#!/usr/bin/env python3
"""Historical USGS quadrangle crops, registered to the lidar tiles.

A scanned topographic sheet is a third witness to mound height, independent of
both the 2006 compilation and the lidar. The delta sheets carry 5-foot
contours, so a mound reads as a nest of closed contours and gives a bracket —
"between 15 and 20 feet" — rather than a number. Where the lidar shows nothing
because the ground has been levelled, an early sheet may still show what was
there, and where the compilation and the lidar disagree, the sheet breaks the
tie.

Each crop uses the SAME extent as that site's lidar tile (`bbox_3857` from
67_mound_lidar_tiles.py), so the two register pixel for pixel and an outline
drawn on the relief lands in the same place on the map.

CHOOSING A SHEET. The earliest available, then the largest scale among sheets
of that year. Two cautions, both recorded per site in the manifest so a reader
can weigh them:

  * USGS dates a sheet by its base survey, not its printing. Both "1940"
    Princedale sheets carry a 13 KV transmission line and postwar road
    casings, so they are revisions of a 1940 base. The year is a lower bound on
    the cartography's age, not its date.
  * A 1:62,500 sheet at 1 m per pixel is far past its own resolution. The crop
    is upsampled; it cannot show more than the original line work holds.

Usage:
    python analyses/70_historic_quad_tiles.py
    python analyses/70_historic_quad_tiles.py --sites 11-N-1,13-N-3
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TILES = ROOT / "output" / "lidar"
QUADS = ROOT / "output" / "quads"
CACHE = QUADS / "_sheets"
TNM = "https://tnmaccess.nationalmap.gov/api/v1/products"
DATASET = "Historical Topographic Maps"
PALETTE_COLORS = 64          # scanned line work quantises well


def sheets_for(lat: float, lon: float) -> list[dict]:
    d = 0.01
    q = urllib.parse.urlencode({
        "bbox": f"{lon - d},{lat - d},{lon + d},{lat + d}",
        "datasets": DATASET, "max": 50})
    with urllib.request.urlopen(f"{TNM}?{q}", timeout=120) as r:
        items = json.load(r).get("items", [])
    out = []
    for it in items:
        title = it.get("title", "")
        year = (it.get("publicationDate") or "")[:4]
        url = it.get("downloadURL") or ""
        if not (year.isdigit() and url.endswith(".pdf")):
            continue
        # "USGS 1:62500-scale Quadrangle for ..." — take the denominator, not
        # the year, which a looser parse happily mistakes for a scale.
        m = re.search(r"1:(\d+)[- ]scale", title)
        scale = int(m.group(1)) if m else 10 ** 9
        out.append({"title": title, "year": int(year), "scale": scale, "url": url})
    # earliest first; among equals prefer the larger scale (smaller denominator)
    out.sort(key=lambda s: (s["year"], s["scale"]))
    return out


def fetch_sheet(url: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / url.rsplit("/", 1)[-1]
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    with urllib.request.urlopen(url, timeout=300) as r:
        dest.write_bytes(r.read())
    return dest


def crop(sheet: Path, bbox_3857, size_px: int, out_png: Path) -> bool:
    """Crop to the lidar tile's extent, so the layers register."""
    from PIL import Image
    tif = out_png.with_suffix(".tif")
    x0, y0, x1, y1 = bbox_3857
    r = subprocess.run(
        ["gdalwarp", "-q", "-overwrite", "-t_srs", "EPSG:3857",
         "-te", str(x0), str(y0), str(x1), str(y1),
         "-ts", str(size_px), str(size_px), "-r", "cubic",
         str(sheet), str(tif)],
        capture_output=True, text=True)
    if r.returncode != 0 or not tif.exists():
        print(f"    gdalwarp failed: {r.stderr.strip()[:120]}", file=sys.stderr)
        return False
    im = Image.open(tif).convert("RGB")
    # A crop falling outside the sheet comes back uniform; that is not a map.
    # getcolors returns None when the image has MORE colours than the cap,
    # which is the healthy case — an earlier version read that as blank and
    # rejected every sheet it was given.
    colours = im.getcolors(maxcolors=16)
    if colours is not None and len(colours) <= 2:
        tif.unlink(missing_ok=True)
        return False
    im.convert("P", palette=Image.ADAPTIVE, colors=PALETTE_COLORS).save(
        out_png, optimize=True)
    tif.unlink(missing_ok=True)
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="")
    args = ap.parse_args()

    QUADS.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((TILES / "manifest.json").read_text())
    want = [s.strip() for s in args.sites.split(",") if s.strip()] or sorted(manifest)

    for site in want:
        meta = manifest.get(site)
        if not meta:
            print(f"{site}: no lidar tile, skipped", file=sys.stderr)
            continue
        key = site.replace("/", "_")
        out_png = QUADS / f"{key}.quad.png"
        try:
            sheets = sheets_for(meta["center_lat"], meta["center_lon"])
        except Exception as exc:                          # noqa: BLE001
            print(f"{site}: sheet query failed ({exc})", file=sys.stderr)
            continue
        placed = None
        for s in sheets[:4]:                              # try oldest first
            try:
                pdf = fetch_sheet(s["url"])
            except Exception as exc:                      # noqa: BLE001
                print(f"  {site}: download failed ({exc})", file=sys.stderr)
                continue
            if crop(pdf, meta["bbox_3857"], meta["size_px"], out_png):
                placed = s
                break
        if not placed:
            print(f"{site}: no usable sheet")
            meta.pop("quad_png", None)
            continue
        meta.update({"quad_png": out_png.name, "quad_title": placed["title"],
                     "quad_year": placed["year"], "quad_scale": placed["scale"],
                     "quad_url": placed["url"]})
        print(f"{site}: {placed['year']} 1:{placed['scale']} "
              f"{placed['title'].split('for')[-1].strip()} "
              f"({out_png.stat().st_size / 1e6:.2f} MB)")
        time.sleep(0.4)

    (TILES / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    n = sum(1 for v in manifest.values() if v.get("quad_png"))
    print(f"\n{n} of {len(manifest)} sites have a quad crop")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
