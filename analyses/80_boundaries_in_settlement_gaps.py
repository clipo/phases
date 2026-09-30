#!/usr/bin/env python3
"""Do the phase boundaries run through the gaps in late-period settlement?

analyses/79_phases_under_drift.py tests whether pottery made by drift, on the
real site geography, would have been sorted into the published phases. This
script tests the step that account leans on: that the phases are located where
they are because SETTLEMENT IS CLUMPED, and a boundary between two phases is a
stretch of country with few sites in it rather than a line where pottery
changes. analyses/75_groupness_surface.py already shows the second half, that
composition changes unusually slowly along the boundaries. This is the first
half, measured on evidence the ceramic analysis never touches.

THE EVIDENCE. The settlement compilation (`data/LMVData.xlsx`, 3,030 recorded
sites) flags each site's components by period. Sites with a period A or B
component are the late pre-contact settlement pattern. They are independent of
the 38 ceramic assemblages: most have no analysed collection at all. Any
compiled site within 1 km of an analysed assemblage is dropped, so the density
measured here is that of the OTHER late-period sites and an assemblage cannot
vouch for its own neighbourhood.

THE MEASUREMENT. A Gaussian kernel density of those sites (sites per 100 km^2)
is read along the internal phase boundaries, which are built exactly as Figure
1 and analyses/75 build them: Voronoi cells of the assemblages dissolved by
phase, the lines where two territories meet, sampled every 1.5 km.

THE COMPARISON. Voronoi boundaries lie between assemblages by construction,
and assemblages are themselves late-period sites, so ANY such boundary runs
through thinner settlement than the assemblages' own locations. The phase
boundaries are therefore judged against boundaries of the same construction
induced by other partitions of the same assemblages into groups of the phases'
own sizes: seeds placed at random (exact minimum-cost assignment), and
partitions drawn to be spatially compact (balanced local search), both from
analyses/74_phase_partition_test.py. If the phases are where the settlement
gaps are, their boundaries should run through LOWER density than those.

Three bandwidths are reported, because a kernel density is a function of its
bandwidth and a result that holds at one only is not a result.

Usage:
    python analyses/80_boundaries_in_settlement_gaps.py [--alt 300]
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

OUT_MD = ROOT / "output" / "findings" / "boundaries_in_settlement_gaps.md"
FIG = "fig16_settlement_gaps"
BANDWIDTHS_KM = (3.0, 5.0, 8.0)
MAIN_BW = 5.0
EXCLUDE_KM = 1.0
STEP_M = 1500.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--alt", type=int, default=300)
    args = ap.parse_args()

    import matplotlib.pyplot as plt
    from pyproj import Transformer
    from shapely.geometry import MultiPoint, Point, box
    from shapely.ops import unary_union, voronoi_diagram
    from mls_emergence.dataio.settlement import load_lmv
    fs = importlib.import_module("figstyle")
    mf = importlib.import_module("make_figures")
    mm = importlib.import_module("make_map")
    ph = importlib.import_module("36_canonical_phase_map")
    t74 = importlib.import_module("74_phase_partition_test")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    xy = coords.to_numpy(float)
    labels_ph, _ = ph.assign_phases_by_territory(names, xy)
    phases = sorted(set(labels_ph))
    pidx = np.array([phases.index(l) for l in labels_ph])
    sizes = np.bincount(pidx)
    slot = np.repeat(np.arange(len(sizes)), sizes)

    tr = Transformer.from_crs("EPSG:4326", mm.UTM15N, always_xy=True)
    E, N = tr.transform(xy[:, 1], xy[:, 0])
    E, N = np.asarray(E, float), np.asarray(N, float)
    pad = 18_000.0
    ext = (E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad)

    # Late-period settlement, independent of the analysed assemblages.
    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    flag = lambda c: pd.to_numeric(lmv[c], errors="coerce").fillna(0) == 1
    late = lmv[(flag("A") | flag("B")) & (pd.to_numeric(lmv["Zone"], errors="coerce") == 15)]
    sx = pd.to_numeric(late["Easting"], errors="coerce").to_numpy(float)
    sy = pd.to_numeric(late["Northing"], errors="coerce").to_numpy(float)
    keep = np.isfinite(sx) & np.isfinite(sy)
    keep &= (sx >= ext[0]) & (sx <= ext[1]) & (sy >= ext[2]) & (sy <= ext[3])
    sx, sy = sx[keep], sy[keep]
    n_in_extent = int(len(sx))
    near = (np.hypot(sx[:, None] - E[None, :], sy[:, None] - N[None, :]).min(1)
            <= EXCLUDE_KM * 1000.0)
    sx, sy = sx[~near], sy[~near]
    if len(sx) < 30:
        raise ValueError(f"only {len(sx)} independent late-period sites in the extent; "
                         "too few for a density")

    def density(px, py, bw_km):
        """Sites per 100 km^2 at each point, Gaussian kernel."""
        h = bw_km * 1000.0
        d2 = (px[:, None] - sx[None, :]) ** 2 + (py[:, None] - sy[None, :]) ** 2
        per_m2 = np.exp(-d2 / (2 * h * h)).sum(1) / (2 * np.pi * h * h)
        return per_m2 * 1e8

    pl = [Point(e, n) for e, n in zip(E, N)]
    cells = list(voronoi_diagram(MultiPoint(pl),
                                 envelope=box(ext[0], ext[2], ext[1], ext[3])).geoms)
    owner = []
    for c in cells:
        o = None
        for j, q in enumerate(pl):
            if c.intersects(q):
                o = j
                break
        owner.append(o)
    envelope = unary_union([q.buffer(16_000) for q in pl])

    def boundary_points(labels):
        polys = {}
        for lab in set(labels):
            mem = [cells[i] for i in range(len(cells))
                   if owner[i] is not None and labels[owner[i]] == lab]
            if mem:
                polys[lab] = unary_union(mem).intersection(envelope)
        keys, pts = list(polys), []
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                sh = polys[keys[a]].boundary.intersection(polys[keys[b]].boundary)
                for g in getattr(sh, "geoms", [sh]):
                    for ln in getattr(g, "geoms", [g]):
                        if ln.geom_type != "LineString" or ln.length == 0:
                            continue
                        for t in np.arange(0, ln.length, STEP_M):
                            q = ln.interpolate(t)
                            pts.append((q.x, q.y))
        return np.array(pts) if pts else np.empty((0, 2))

    bp = boundary_points(list(pidx))
    km = t74.km_xy(xy)
    rng = np.random.default_rng(80)
    alt_pts = {"random seeds": [], "spatially compact": []}
    for r in range(args.alt):
        seeds = km[rng.choice(len(km), size=len(sizes), replace=False)]
        lab = t74.assign_exact(km, seeds, slot)
        alt_pts["random seeds"].append(boundary_points(list(lab)))
        alt_pts["spatially compact"].append(
            boundary_points(list(t74.local_search(km, lab, sizes))))
        if (r + 1) % 50 == 0:
            print(f"  partitions {r + 1}/{args.alt}", flush=True)

    rows = []
    for bw in BANDWIDTHS_KM:
        obs = float(np.median(density(bp[:, 0], bp[:, 1], bw)))
        at_sites = float(np.median(density(E, N, bw)))
        rec = dict(bandwidth_km=bw, phase_boundaries=obs, at_assemblages=at_sites)
        for kind, plist in alt_pts.items():
            d = np.array([np.median(density(p[:, 0], p[:, 1], bw)) for p in plist if len(p)])
            rec[f"{kind}: median"] = float(np.median(d))
            rec[f"{kind}: share of partitions with DENSER boundaries"] = float(np.mean(d > obs))
            if bw == MAIN_BW:
                rec[f"_{kind}"] = d
        rows.append(rec)

    # Figure: the density surface with sites and boundaries; the matched comparison.
    main_row = next(r for r in rows if r["bandwidth_km"] == MAIN_BW)
    gx, gy = np.meshgrid(np.linspace(ext[0], ext[1], 170), np.linspace(ext[2], ext[3], 170))
    surf = density(gx.ravel(), gy.ravel(), MAIN_BW).reshape(gx.shape)
    fig = plt.figure(figsize=(7.2, 3.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.55)
    ax = fig.add_subplot(gs[0, 0])
    mm.basin_basemap(ax, ext, geology=False, grayscale=True, show_counties=False,
                     show_states=True, draw_rivers=True)
    im = ax.pcolormesh(gx, gy, surf, cmap="Greys", vmin=0,
                       vmax=float(np.percentile(surf, 99)), alpha=0.75, shading="auto",
                       zorder=1.2)
    ax.plot(sx, sy, ".", ms=1.6, color="0.25", zorder=3)
    ax.plot(E, N, "o", ms=2.6, mfc="white", mec="0.05", mew=0.6, zorder=4)
    ax.plot(bp[:, 0], bp[:, 1], ".", ms=1.4, color="0.0", zorder=5)
    ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("other late-period sites\nper 100 km$^2$", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    fs.panel_label(ax, "A")
    axB = fig.add_subplot(gs[0, 1])
    axB.hist(main_row["_random seeds"], bins=24, color="0.80", edgecolor="0.55", lw=.4,
             label="same sizes, seeds at random")
    axB.hist(main_row["_spatially compact"], bins=24, color="0.45", edgecolor="0.25", lw=.4,
             alpha=0.75, label="same sizes, spatially compact")
    axB.axvline(main_row["phase_boundaries"], color="0.05", lw=1.8, label="the phase boundaries")
    axB.set_xlabel("median settlement density along the boundaries\n(sites per 100 km$^2$)")
    axB.set_ylabel("partitions")
    axB.legend(fontsize=5.5, frameon=False)
    fs.panel_label(axB, "B")
    png = fs.save_all(fig, FIG, close=True)

    L = ["# Do the phase boundaries run through the gaps in late-period settlement?", "",
         f"{n_in_extent} compiled sites with a period A or B component fall in the mapped "
         f"area; {int(near.sum())} lie within {EXCLUDE_KM:g} km of an analysed assemblage and "
         f"are dropped, leaving {len(sx)} independent late-period sites.",
         f"Phase boundaries: {len(bp)} points at {STEP_M / 1000:g} km spacing. "
         f"{args.alt} alternative partitions of each kind, all carrying the phases' group "
         f"sizes ({', '.join(str(int(s)) for s in sizes)}).", "",
         "| kernel bandwidth | at the assemblages | along the phase boundaries | "
         "random-seed boundaries, median | share denser than the phases' | "
         "compact-partition boundaries, median | share denser than the phases' |",
         "|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['bandwidth_km']:g} km | {r['at_assemblages']:.2f} | "
                 f"**{r['phase_boundaries']:.2f}** | {r['random seeds: median']:.2f} | "
                 f"{100 * r['random seeds: share of partitions with DENSER boundaries']:.0f}% | "
                 f"{r['spatially compact: median']:.2f} | "
                 f"{100 * r['spatially compact: share of partitions with DENSER boundaries']:.0f}% |")
    L += ["", "Densities are other late-period sites per 100 km^2. \"Share denser\" is the "
              "share of alternative",
          "partitions whose boundaries run through MORE settlement than the phase "
          "boundaries do. A high share",
          "means the phase boundaries sit in unusually empty country even for boundaries of "
          "their own kind;",
          "a share near half means they are ordinary in that respect. The placement is "
          "descriptive (rule 18).", "",
          f"Figure written to {png.name} and its siblings.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:12]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
