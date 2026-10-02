"""build_phillips1970_outlines.py - the Phillips (1970) phase outlines, as coordinates.

Writes `data/raw/phillips1970_phase_outlines.csv`, which
`analyses/93_phillips1970_phase_map.py` draws as manuscript Figure 1.

WHAT THE SOURCE IS. Lipo (2001: Figure 2.6, printed page 20), "Phase Areas as
Defined by Phillips (1970)", redraws Phillips's Nodena, Parkin, Walls and Kent
phase areas on the Lower Mississippi Survey quadrangle grid (columns N to Q,
rows 10 to 14). In the BAR volume's PDF that figure is vector artwork, so the
four outlines are read from the page's own paths rather than traced from a
scan. The outlines are therefore Lipo's redrawing of Phillips, not a new
reading of Phillips's plates (rule 13: the caption says so).

HOW THE OUTLINES ARE FOUND. `pdftocairo -svg` (poppler) converts the page; the
outlines are the only paths stroked at 2.1 pt. Cubic segments (the Nodena
outline) are sampled at BEZIER_STEPS points each; the others are polylines and
are kept vertex for vertex. Each outline is named by the type site it
encloses, and the script raises unless exactly four outlines are found and
each encloses exactly its own type site.

HOW THEY ARE GEOREFERENCED. The figure's frame is four 15-minute quadrangle
columns by five rows, so its corners are known: longitude -90.75 to -89.75,
latitude 34.50 to 35.75. Longitude is interpolated piecewise through the frame
and the three interior column lines the figure draws; latitude is linear
between the frame's top and bottom (the figure's interior row lines are not
single paths and are not used). The quadrangle bounds were read from the
digitized LMS point layers under data/CMV: every site designated 11-N-x falls
in latitude 35.25 to 35.50 and longitude -90.75 to -90.50, and so on.

HOW THAT IS CHECKED. The figure also plots site dots. For each analysis
assemblage inside the frame, the script finds the nearest plotted dot under
the georeference above and reports the distance; the summary goes into the
CSV header. This is a check on the georeference, not a calibration: nothing is
fitted to it.

The PDF is not in the repository (docs/references/pdfs/*.pdf is ignored), so
the path is an argument and the CSV is the tracked artifact.

Usage: .venv/bin/python scripts/build_phillips1970_outlines.py --pdf PATH [--page 32]
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from mls_emergence.dataio.coords import read_assemblage_xy  # noqa: E402

OUT = ROOT / "data" / "raw" / "phillips1970_phase_outlines.csv"
XY = ROOT / "data" / "raw" / "mainfort-pfg-cplXY.txt"

OUTLINE_STROKE_PT = 2.1
FRAME_STROKE_PT = 1.4
BEZIER_STEPS = 12
PAGE_HEIGHT_PT = 842.0
# The frame is LMS columns N to Q by rows 10 to 14, 15 minutes each.
LON_WEST, LON_EAST = -90.75, -89.75
LAT_SOUTH, LAT_NORTH = 34.50, 35.75
# Each outline is named by the type site it must enclose (project coordinates).
TYPE_SITE = {"Nodena": "Upper_Nodena", "Parkin": "Parkin", "Walls": "Walls",
             "Kent": "Kent_Place"}
# A plotted dot further than this from an assemblage is treated as "not
# plotted in the figure" rather than as a georeferencing residual.
NOT_PLOTTED_KM = 3.0


def _page_svg(pdf: Path, page: int) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "page.svg"
        subprocess.run(["pdftocairo", "-f", str(page), "-l", str(page), "-svg",
                        str(pdf), str(out)], check=True)
        return out.read_text()


def _paths(svg: str):
    """Yield (stroke_width or None, d, transform) for every <path>, with the
    path data still in its own coordinates."""
    for tag in re.findall(r"<path[^>]*/>", svg):
        d = re.search(r' d="([^"]*)"', tag)
        if d is None:
            continue
        sw = re.search(r'stroke-width="([0-9.]+)"', tag)
        tr = re.search(r'transform="matrix\(([^)]*)\)"', tag)
        mat = [float(v) for v in tr.group(1).split(",")] if tr else [1, 0, 0, 1, 0, 0]
        yield (float(sw.group(1)) if sw else None), d.group(1), mat


def _to_page(xy: np.ndarray, mat) -> np.ndarray:
    """SVG path coordinates -> PDF points with y up."""
    a, b, c, d, e, f = mat
    x = a * xy[:, 0] + c * xy[:, 1] + e
    y = b * xy[:, 0] + d * xy[:, 1] + f
    return np.column_stack([x, PAGE_HEIGHT_PT - y])


def _subpaths(d: str) -> list[np.ndarray]:
    """Flatten a path's subpaths to polylines, sampling cubic segments."""
    out = []
    for sub in re.findall(r"M[^M]*", d):
        toks = re.findall(r"[MLCZ]|-?[0-9]+\.?[0-9]*", sub)
        pts, i = [], 0
        while i < len(toks):
            t = toks[i]
            if t in "ML":
                pts.append((float(toks[i + 1]), float(toks[i + 2])))
                i += 3
            elif t == "C":
                p0 = np.array(pts[-1])
                c1, c2, p3 = (np.array([float(toks[i + 1 + 2 * k]), float(toks[i + 2 + 2 * k])])
                              for k in range(3))
                for s in np.linspace(0, 1, BEZIER_STEPS + 1)[1:]:
                    pts.append(tuple((1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * c1
                                     + 3 * (1 - s) * s ** 2 * c2 + s ** 3 * p3))
                i += 7
            elif t == "Z":
                i += 1
            else:
                raise ValueError(f"unexpected path token {t!r}")
        if len(pts) > 2:
            out.append(np.array(pts))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", type=Path, required=True,
                    help="Lipo (2001), BAR International Series 918, as PDF")
    ap.add_argument("--page", type=int, default=32,
                    help="PDF sheet carrying Figure 2.6 (printed page 20)")
    args = ap.parse_args()
    if not args.pdf.exists():
        raise SystemExit(f"{args.pdf} not found: the BAR volume is not tracked "
                         "(docs/references/pdfs/*.pdf is ignored); pass its path")

    svg = _page_svg(args.pdf, args.page)
    outlines, frame_x, frame_y, col_x, dots = [], [], [], [], []
    for sw, d, mat in _paths(svg):
        nums = np.array([float(v) for v in re.findall(r"-?[0-9]+\.?[0-9]*", d)])
        if len(nums) < 4 or len(nums) % 2:
            continue
        xy = _to_page(nums.reshape(-1, 2), mat)
        w, h = np.ptp(xy[:, 0]), np.ptp(xy[:, 1])
        if sw == OUTLINE_STROKE_PT:
            outlines += [_to_page(s, mat) for s in _subpaths(d)]
        elif sw == FRAME_STROKE_PT and len(xy) == 2:
            (frame_y if h < 0.5 else frame_x).append(xy.mean(0)[1 if h < 0.5 else 0])
        elif sw is not None and len(xy) == 2 and w < 0.5 and h > 300:
            col_x.append(xy[:, 0].mean())
        elif sw is None and 1.0 < w < 4.0 and abs(w - h) < 0.5:
            dots.append(xy.min(0) + np.array([w, h]) / 2)

    if len(frame_x) != 2 or len(frame_y) != 2 or len(col_x) != 3:
        raise SystemExit(f"frame not recovered: {len(frame_x)} vertical and "
                         f"{len(frame_y)} horizontal frame lines, {len(col_x)} "
                         "interior column lines (expected 2, 2, 3)")
    x_knots = np.array(sorted(frame_x)[:1] + sorted(col_x) + sorted(frame_x)[1:])
    lon_knots = np.linspace(LON_WEST, LON_EAST, 5)
    y0, y1 = sorted(frame_y)

    def to_lonlat(xy):
        lon = np.interp(xy[:, 0], x_knots, lon_knots)
        lat = LAT_SOUTH + (xy[:, 1] - y0) / (y1 - y0) * (LAT_NORTH - LAT_SOUTH)
        return np.column_stack([lon, lat])

    if len(outlines) != 4:
        raise SystemExit(f"expected four phase outlines, found {len(outlines)}")

    sites = read_assemblage_xy(XY, verbose=False).set_index("Assemblages")
    site_ll = {n: (float(r.Longitude), float(r.Latitude)) for n, r in sites.iterrows()}
    rows, named = [], {}
    for pts in outlines:
        ll = to_lonlat(pts)
        poly = Polygon(ll)
        if not poly.is_valid:
            poly = poly.buffer(0)
        inside = [ph for ph, s in TYPE_SITE.items() if poly.contains(Point(site_ll[s]))]
        if len(inside) != 1:
            raise SystemExit(f"an outline encloses the type sites of {inside}; "
                             "expected exactly one")
        named[inside[0]] = (pts, ll)
    if set(named) != set(TYPE_SITE):
        raise SystemExit(f"outlines found for {sorted(named)} only")
    for ph in TYPE_SITE:
        pts, ll = named[ph]
        for k, ((x, y), (lon, lat)) in enumerate(zip(pts, ll)):
            rows.append((ph, k, round(lon, 5), round(lat, 5), round(x, 2), round(y, 2)))

    # Georeference check: nearest plotted dot to each assemblage in the frame.
    dots_ll = to_lonlat(np.array(dots))
    in_frame = (dots_ll[:, 0] > LON_WEST) & (dots_ll[:, 0] < LON_EAST) \
        & (dots_ll[:, 1] > LAT_SOUTH) & (dots_ll[:, 1] < LAT_NORTH)
    dots_ll = dots_ll[in_frame]
    resid = {}
    for n, (lon, lat) in site_ll.items():
        if not (LON_WEST < lon < LON_EAST and LAT_SOUTH < lat < LAT_NORTH):
            continue
        dx = (dots_ll[:, 0] - lon) * 111.32 * np.cos(np.radians(lat))
        dy = (dots_ll[:, 1] - lat) * 111.32
        resid[n] = float(np.hypot(dx, dy).min())
    plotted = {n: r for n, r in resid.items() if r <= NOT_PLOTTED_KM}
    absent = sorted(n for n, r in resid.items() if r > NOT_PLOTTED_KM)
    r = np.array(sorted(plotted.values()))

    header = [
        "# Phase outlines of Phillips (1970) as redrawn by Lipo (2001: Figure 2.6, printed p. 20),",
        "# read from the vector paths of the BAR volume's PDF by scripts/build_phillips1970_outlines.py.",
        "# x_pt, y_pt are PDF points on that page (y up); longitude and latitude follow from the figure's",
        f"# frame, LMS columns N-Q by rows 10-14: longitude {LON_WEST} to {LON_EAST}, latitude {LAT_SOUTH} to {LAT_NORTH}.",
        "# Column lines at x_pt = " + ", ".join(f"{v:.2f}" for v in x_knots)
        + f"; frame bottom and top at y_pt = {y0:.2f}, {y1:.2f}.",
        f"# Georeference check, nothing fitted: of {len(resid)} project assemblages inside the frame, "
        f"{len(plotted)} have a plotted site dot within {NOT_PLOTTED_KM:g} km;",
        f"# for those the nearest dot is a median {np.median(r):.2f} km away (maximum {r.max():.2f} km). "
        f"No dot within {NOT_PLOTTED_KM:g} km: " + (", ".join(absent) if absent else "none") + ".",
        "# The Parkin outline is open on its west side in the source, where it runs into the frame;",
        "# readers close it with a straight segment (first vertex to last).",
    ]
    df = pd.DataFrame(rows, columns=["phase", "vertex", "longitude", "latitude", "x_pt", "y_pt"])
    with open(OUT, "w") as fh:
        fh.write("\n".join(header) + "\n")
        df.to_csv(fh, index=False)
    print("\n".join(header))
    print(df.groupby("phase", sort=False).size().to_string())
    for n in sorted(resid, key=resid.get):
        print(f"  {n:20s} {resid[n]:6.2f} km")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
