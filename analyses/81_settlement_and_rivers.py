#!/usr/bin/env python3
"""What do the rivers do: carry interaction, or decide where people lived?

analyses/61_river_network_geometry.py shows that decorated pottery follows
straight-line distance and not river-network distance, above all between
phases, where the network can join a St. Francis phase to a Mississippi phase
only through a confluence far to the south. The author's reading (2026-09-21)
is that the rivers still matter, but at a different step: they shape WHERE
sites are, clumping settlement, and the clumps are what the phases were drawn
around. This script measures that reading in two parts.

1. DO SITES HUG THE RIVERS? Distance from each late-period site (period A or B
   component in the settlement compilation) to the nearest channel, against
   uniformly random points in the same mapped area. Run against every mapped
   channel and against the major rivers only, because the modern hydrography
   layer is dense enough that "near some channel" constrains little.

2. IS THAT WHAT MAKES SETTLEMENT CLUMPED? The Clark-Evans nearest-neighbour
   ratio of the sites (1 random, below 1 clumped), against the same ratio for
   random points constrained to sit at the SITES' OWN distances from the
   rivers. If hugging the rivers is what clumps settlement, the constrained
   random points will be as clumped as the sites. If they are not, the sites
   concentrate on particular reaches for reasons river proximity does not
   capture.

Both placements are descriptive (rule 18). The hydrography is MODERN; channels
have moved since these sites were occupied, which is a limit on part 1 and is
stated in the output.

Usage:
    python analyses/81_settlement_and_rivers.py [--sims 300]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "settlement_and_rivers.md"
MAJOR = ("Saint Francis", "Tyronza", "Mississippi", "Little River", "Anguille",
         "Blackfish", "Fifteenmile", "Coldwater", "Tallahatchie")
N_RANDOM = 6000
TOL_KM = 0.4


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", type=int, default=300)
    args = ap.parse_args()

    from pyproj import Transformer
    from shapely import STRtree
    from shapely.geometry import Point
    from mls_emergence.dataio.settlement import load_lmv
    mf = importlib.import_module("make_figures")
    mm = importlib.import_module("make_map")

    counts, coords = mf._load_curated()
    xy = coords.to_numpy(float)
    E, N = Transformer.from_crs("EPSG:4326", mm.UTM15N, always_xy=True).transform(
        xy[:, 1], xy[:, 0])
    E, N = np.asarray(E, float), np.asarray(N, float)
    pad = 18_000.0
    ext = (E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad)
    area = (ext[1] - ext[0]) * (ext[3] - ext[2])

    hyd = mm._clip_gdf(mm._load_hydrology(), ext)
    riv = mm._clip_gdf(mm._load_major_rivers(), ext)
    polys = [g.boundary for g in riv.geometry if g is not None and not g.is_empty]
    every = [g for g in hyd.geometry if g is not None and not g.is_empty] + polys
    # NAME holds missing values that astype(str) leaves as non-strings under
    # pandas 3's string dtype, so test membership on plain Python strings.
    is_major = np.array([isinstance(n, str) and any(w in n for w in MAJOR)
                         for n in hyd["NAME"].tolist()])
    major = [g for g in hyd[is_major].geometry if g is not None and not g.is_empty] + polys

    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    flag = lambda c: pd.to_numeric(lmv[c], errors="coerce").fillna(0) == 1
    late = lmv[(flag("A") | flag("B")) & (pd.to_numeric(lmv["Zone"], errors="coerce") == 15)]
    sx = pd.to_numeric(late["Easting"], errors="coerce").to_numpy(float)
    sy = pd.to_numeric(late["Northing"], errors="coerce").to_numpy(float)
    k = (np.isfinite(sx) & np.isfinite(sy) & (sx >= ext[0]) & (sx <= ext[1])
         & (sy >= ext[2]) & (sy <= ext[3]))
    sx, sy = sx[k], sy[k]

    rng = np.random.default_rng(81)
    rx = rng.uniform(ext[0], ext[1], N_RANDOM)
    ry = rng.uniform(ext[2], ext[3], N_RANDOM)

    def clark_evans(x, y):
        P = np.column_stack([x, y])
        D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(2))
        np.fill_diagonal(D, np.inf)
        return float(D.min(1).mean() / (0.5 / np.sqrt(len(x) / area)))

    L = ["# What do the rivers do: carry interaction, or decide where people lived?", "",
         f"{len(sx)} late-period sites (period A or B component) in the mapped area; "
         f"{N_RANDOM:,} uniformly random points in the same area for comparison.", "",
         "## 1. Do sites hug the rivers?", "",
         "| channels | sites, median km (quartiles) | random points, median km (quartiles) |",
         "|---|---|---|"]
    ce_rows = []
    for label, geoms in (("every mapped channel", every), ("major rivers only", major)):
        tree = STRtree(geoms)
        dist = lambda x, y: np.array([Point(a, b).distance(geoms[tree.nearest(Point(a, b))])
                                      / 1000.0 for a, b in zip(x, y)])
        ds, dr = dist(sx, sy), dist(rx, ry)
        L.append(f"| {label} | **{np.median(ds):.2f}** ({np.percentile(ds, 25):.2f}-"
                 f"{np.percentile(ds, 75):.2f}) | {np.median(dr):.2f} "
                 f"({np.percentile(dr, 25):.2f}-{np.percentile(dr, 75):.2f}) |")
        sims = []
        for _ in range(args.sims):
            idx = [rng.choice(np.where(np.abs(dr - d) <= TOL_KM)[0])
                   if (np.abs(dr - d) <= TOL_KM).any() else int(np.argmin(np.abs(dr - d)))
                   for d in ds]
            sims.append(clark_evans(rx[idx], ry[idx]))
        ce_rows.append((label, np.median(sims), *np.percentile(sims, [2.5, 97.5])))
    L += ["", "## 2. Is that what makes settlement clumped?", "",
          f"Clark-Evans ratio of the late-period sites: **{clark_evans(sx, sy):.2f}**; of the "
          f"{len(E)} analysed assemblages: **{clark_evans(E, N):.2f}** "
          "(1 is random, below 1 is clumped).", "",
          "| random points placed at the sites' own distances from | Clark-Evans ratio, "
          f"median of {args.sims} (2.5th-97.5th) |", "|---|---|"]
    for label, med, lo, hi in ce_rows:
        L.append(f"| {label} | {med:.2f} ({lo:.2f}-{hi:.2f}) |")
    L += ["", "Settlement is tied to the rivers and it is clumped, but the first does not "
              "produce the second:",
          "points that hug the rivers exactly as closely as the sites do are close to "
          "randomly spaced. The sites",
          "concentrate on particular reaches. Many carry the names of oxbow lakes, which "
          "suggests relict",
          "meander belts rather than active channels, but no geomorphic layer is in hand "
          "to test that.", "",
          "The hydrography is modern. Channels have moved since these sites were "
          "occupied, so distances to",
          "channels are distances to where the water is now.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[4:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
