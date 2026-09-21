#!/usr/bin/env python3
"""Is groupness a gradient or an edge, and at what scale does drift fail?

The phase scheme asserts edges: named territories with boundaries between them.
Every other measurement in this paper is taken ACROSS a partition the analysis
draws, which cannot see whether an edge exists, only how much differentiation a
given cut produces. This script asks the question the map form poses directly.

PANELS A AND B. COMPOSITIONAL TURNOVER AS A CONTINUOUS FIELD. At each point of a
grid over the basin, the ten decorated-class proportions are regressed on
position by weighted least squares, with Gaussian weights of bandwidth
BANDWIDTH_KM on the assemblages' distance to that point. The fitted gradient of
each class gives a local rate of compositional change, and their root sum of
squares is the TURNOVER RATE at that point, in proportion units per kilometre.

That quantity is what an edge would show up in. A boundary between two
interaction communities is a place where composition changes fast over a short
distance, so it appears as a ridge in this field. A gradient with no boundaries
appears as a smooth surface with no ridges. Panel A shows the field alone;
panel B draws the published phase boundaries on the same field, so the question
"do the boundaries lie on the ridges" is answered by looking, and by the number
reported beside it.

THAT NUMBER NEEDS A MATCHED COMPARISON, and the obvious one is wrong. Comparing
turnover along the boundaries against turnover over the whole mapped area gives
the phases the 19th percentile, which looks like a strong result and is an
artefact. Turnover is inversely related to local support: Spearman -0.40
between the two, median 0.0159 per km where the local weight sum is 2 to 4
against 0.0108 where it is 6 to 9, because a regression fitted to few nearby
assemblages extrapolates steeply. Voronoi boundaries lie in the INTERIOR of the
site distribution, where support is highest and the field is therefore
flattest, so any set of boundaries built this way would score low against the
whole field.

The comparison reported here instead holds the construction fixed: N_ALT
alternative partitions with the phases' own group sizes, seeds placed at
random, assigned by exact minimum-cost matching (the generator from
analyses/74_phase_partition_test.py), each dissolved into territories and
sampled for turnover exactly as the phases are. That asks whether THESE
boundaries lie on ridges, against boundaries of the same kind drawn elsewhere.

PANEL C. WHERE DRIFT STOPS ACCOUNTING FOR IT. The partition grain is converted
from a cluster count, which means nothing on the ground, to a length: the root
mean squared distance from an assemblage to its group's centroid, the group's
characteristic radius. Observed differentiation is plotted against that radius
with the calibrated spatial-drift band behind it, so the reader can see at what
spatial scale the record leaves what drift produces, and where the phases fall
on that axis.

WHAT THE MEASUREMENT CANNOT DO, stated because the panels invite more than they
support. The turnover field is an interpolation between 43 points, so it has no
information at scales finer than the spacing between them, and its ridges near
the edge of the mapped area rest on assemblages on one side only. Points whose
local weight sum falls below MIN_WEIGHT are masked rather than drawn faint.
Radius is one way to give a partition a length and is not an interaction
distance. And the whole field inherits the space-time confound: the seriation
axis derives from these same frequencies, so composition changing across space
and composition changing through time are not separated here.

Usage:
    python analyses/75_groupness_surface.py [--grid 220] [--bandwidth 15]
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

OUT_MD = ROOT / "output" / "findings" / "groupness_surface.md"
SWEEP_CSV = ROOT / "output" / "findings" / "scale_km_sweep.csv"
FIG = "fig13_groupness_surface"
BANDWIDTH_KM = 15.0
N_ALT = 200               # alternative partitions the phase boundaries are judged against
MIN_WEIGHT = 2.0          # effective assemblages behind a grid point
PHASES = ("Parkin", "Nodena", "Kent", "Walls", "Parchman")


def turnover_field(pts_km, props, gx, gy, bandwidth, min_weight):
    """Local rate of compositional change, per km, at each grid point.

    Weighted least squares of every class proportion on centred position. The
    slope terms are the gradient; their root sum of squares over classes is the
    turnover rate. Returns the field and a mask of points with enough support.
    """
    g = np.column_stack([gx.ravel(), gy.ravel()])
    out = np.full(len(g), np.nan)
    ok = np.zeros(len(g), bool)
    for i, p in enumerate(g):
        d2 = ((pts_km - p) ** 2).sum(1)
        w = np.exp(-d2 / (2.0 * bandwidth ** 2))
        if w.sum() < min_weight:
            continue
        X = np.column_stack([np.ones(len(pts_km)), pts_km - p])
        XtW = X.T * w
        try:
            beta = np.linalg.solve(XtW @ X, XtW @ props)
        except np.linalg.LinAlgError:
            continue
        out[i] = float(np.sqrt((beta[1:] ** 2).sum()))
        ok[i] = True
    return out.reshape(gx.shape), ok.reshape(gx.shape)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=220)
    ap.add_argument("--bandwidth", type=float, default=BANDWIDTH_KM)
    ap.add_argument("--alt", type=int, default=N_ALT)
    args = ap.parse_args()

    import matplotlib.pyplot as plt
    from shapely.geometry import MultiPoint, Point, box
    from shapely.ops import unary_union, voronoi_diagram
    fs = importlib.import_module("figstyle")
    t74 = importlib.import_module("74_phase_partition_test")
    mf = importlib.import_module("make_figures")
    mm = importlib.import_module("make_map")
    ph = importlib.import_module("36_canonical_phase_map")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m = counts.to_numpy(float)
    props = m / m.sum(1, keepdims=True)
    xy = coords.to_numpy(float)

    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", mm.UTM15N, always_xy=True)
    E, N = tr.transform(xy[:, 1], xy[:, 0])
    E, N = np.asarray(E, float), np.asarray(N, float)
    pts_km = np.column_stack([E, N]) / 1000.0

    labels_ph, _ = ph.assign_phases_by_territory(names, xy)

    pad = 18_000.0
    ext = (E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad)
    gx, gy = np.meshgrid(np.linspace(ext[0], ext[1], args.grid),
                         np.linspace(ext[2], ext[3], args.grid))
    field, ok = turnover_field(pts_km, props, gx / 1000.0, gy / 1000.0,
                               args.bandwidth, MIN_WEIGHT)
    masked = np.ma.masked_where(~ok, field)

    # Territories exactly as Figure 1 builds them, so the phase boundaries drawn
    # here are the same lines the phase map shows, not a second construction.
    # The same function builds the alternative partitions' boundaries, so the
    # comparison holds the construction fixed and varies only the labels.
    pl = [Point(e, n) for e, n in zip(E, N)]
    cells = list(voronoi_diagram(MultiPoint(pl),
                                 envelope=box(ext[0], ext[2], ext[1], ext[3])).geoms)
    cell_owner = []
    for c in cells:
        owner = None
        for j, q in enumerate(pl):
            if c.intersects(q):
                owner = j
                break
        cell_owner.append(owner)
    envelope = unary_union([q.buffer(16_000) for q in pl])

    def boundaries_of(labels):
        """Internal boundaries of the territories a labelling induces.

        Only where two territories MEET. The outer envelope is an artefact of
        where we stopped drawing and is not a boundary of anything.
        """
        polys = {}
        for lab in set(labels):
            member = [cells[i] for i in range(len(cells))
                      if cell_owner[i] is not None and labels[cell_owner[i]] == lab]
            if member:
                polys[lab] = unary_union(member).intersection(envelope)
        keys, internal = list(polys), []
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                shared = polys[keys[a]].boundary.intersection(polys[keys[b]].boundary)
                if not shared.is_empty:
                    internal.append(shared)
        return unary_union(internal) if internal else None

    bnd = boundaries_of(list(labels_ph))

    # Where does turnover along those boundaries sit in the field as a whole?
    def sample(geom, step=1500.0):
        lines = getattr(geom, "geoms", [geom])
        pts = []
        for ln in lines:
            if ln.geom_type not in ("LineString", "MultiLineString"):
                continue
            for sub in getattr(ln, "geoms", [ln]):
                if sub.length == 0:
                    continue
                for t in np.arange(0, sub.length, step):
                    q = sub.interpolate(t)
                    pts.append((q.x, q.y))
        return np.array(pts) if pts else np.empty((0, 2))

    def read_field(P):
        if len(P) == 0:
            return np.array([])
        ix = np.clip(np.searchsorted(gx[0], P[:, 0]) - 1, 0, gx.shape[1] - 1)
        iy = np.clip(np.searchsorted(gy[:, 0], P[:, 1]) - 1, 0, gx.shape[0] - 1)
        v = field[iy, ix]
        return v[np.isfinite(v)]

    on_b = read_field(sample(bnd)) if bnd is not None else np.array([])
    bg = field[ok]
    pct_field = 100.0 * float(np.mean(bg < np.median(on_b))) if len(on_b) else float("nan")

    # The comparison that controls for where Voronoi boundaries can fall.
    km = t74.km_xy(xy)
    sizes = np.bincount(np.array([sorted(set(labels_ph)).index(l) for l in labels_ph]))
    slot = np.repeat(np.arange(len(sizes)), sizes)
    rng = np.random.default_rng(75)
    alt_med = []
    for r in range(args.alt):
        seeds = km[rng.choice(len(km), size=len(sizes), replace=False)]
        lab = t74.assign_exact(km, seeds, slot)
        gb = boundaries_of(list(lab))
        if gb is None:
            continue
        v = read_field(sample(gb))
        if len(v):
            alt_med.append(float(np.median(v)))
        if (r + 1) % 50 == 0:
            print(f"  matched boundaries {r + 1}/{args.alt}", flush=True)
    alt_med = np.array(alt_med)
    obs_med = float(np.median(on_b)) if len(on_b) else float("nan")
    pct_alt = (100.0 * float(np.mean(alt_med < obs_med))
               if len(alt_med) and np.isfinite(obs_med) else float("nan"))

    sweep = pd.read_csv(SWEEP_CSV) if SWEEP_CSV.exists() else None

    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.95], wspace=0.06)
    vmax = float(np.nanpercentile(field[ok], 99))
    for col, (title, draw_b) in enumerate(((None, False), (None, True))):
        ax = fig.add_subplot(gs[0, col])
        mm.basin_basemap(ax, ext, geology=False, grayscale=True,
                         show_counties=False, show_states=True, draw_rivers=True)
        im = ax.pcolormesh(gx, gy, masked, cmap="Greys", vmin=0, vmax=vmax,
                           alpha=0.78, shading="auto", zorder=1.2)
        ax.plot(E, N, "o", ms=2.4, mfc="white", mec="0.1", mew=0.5, zorder=4)
        if draw_b and bnd is not None:
            for ln in getattr(bnd, "geoms", [bnd]):
                if ln.geom_type == "LineString":
                    x, y = ln.xy
                    ax.plot(x, y, color="0.05", lw=1.3, zorder=5)
        ax.set_xticks([]); ax.set_yticks([])
        fs.panel_label(ax, "AB"[col])
        if col == 0:
            mm._add_scale_bar(ax, ext, bar_km=25)
            cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
            cb.set_label("compositional turnover\n(per km)", fontsize=5.5)
            cb.ax.tick_params(labelsize=5)

    axC = fig.add_subplot(gs[0, 2])
    if sweep is not None:
        axC.fill_between(sweep.radius_km, sweep.drift_lo, sweep.drift_hi,
                         color="0.82", label="calibrated drift, 95%")
        axC.plot(sweep.radius_km, sweep.observed, color="0.1", lw=1.5,
                 marker="o", ms=3, label="observed")
        axC.set_xlabel("group radius (km)")
        axC.set_ylabel("cultural $F_{ST}$")
        axC.invert_xaxis()
        axC.legend(fontsize=5.5, frameon=False, loc="upper left")
    fs.panel_label(axC, "C")

    png = fs.save_all(fig, FIG, close=True)

    L = ["# Is groupness a gradient or an edge?", "",
         f"Basin phase set, {len(names)} assemblages, 10 decorated classes. "
         f"Turnover is the root sum of squares of the",
         f"weighted-least-squares gradient of the class proportions, Gaussian "
         f"bandwidth {args.bandwidth:g} km,",
         f"on a {args.grid} x {args.grid} grid, masked where the local weight sum "
         f"falls below {MIN_WEIGHT:g}.", "",
         "## Do the phase boundaries lie on ridges of turnover?", ""]
    if len(on_b):
        L += [f"Median turnover along the internal phase boundaries is "
              f"**{obs_med:.5f}** per km,",
              f"from {len(on_b)} points sampled at 1.5 km spacing.", "",
              "**Against boundaries of the same construction drawn elsewhere** "
              f"({len(alt_med)} partitions",
              "carrying the phases' group sizes, seeds at random, assigned by exact "
              "minimum-cost",
              f"matching, each dissolved into territories the same way): median "
              f"{np.median(alt_med):.5f} per km,",
              f"and the phase boundaries sit at the **{t74.ordinal(pct_alt)} percentile**. "
              "This is the number to read.", "",
              "**Against the whole mapped area** (median "
              f"{np.median(bg):.5f} over {int(ok.sum()):,} unmasked grid points): "
              f"the {t74.ordinal(pct_field)}",
              "percentile. **Do not read this one.** Turnover falls where local support "
              "rises",
              "(Spearman -0.40; median 0.0159 per km where the local weight sum is 2 to "
              "4, against",
              "0.0108 where it is 6 to 9), because a regression fitted to few nearby "
              "assemblages",
              "extrapolates steeply. Voronoi boundaries lie in the interior of the site "
              "distribution,",
              "where support is highest, so ANY boundaries built this way score low "
              "against the whole",
              "field. The matched comparison above removes that.", "",
              "A boundary between interaction communities is a place where composition "
              "changes fast",
              "over a short distance, so it would sit high against the matched "
              "comparison.", ""]
    else:
        L += ["No internal boundaries were recovered; nothing to report here.", ""]
    if sweep is not None:
        first = sweep[sweep.above].iloc[0] if sweep.above.any() else None
        inside = sweep[~sweep.above]
        L += ["## At what spatial scale does drift stop accounting for it?", "",
              "| group radius (km) | k | observed F_ST | drift 95% upper | above? | smallest group |",
              "|---|---|---|---|---|---|"]
        for _, r in sweep.iterrows():
            L.append(f"| {r.radius_km:.1f} | {int(r.k)} | {r.observed:.4f} | "
                     f"{r.drift_hi:.4f} | {'yes' if r.above else 'no'} | "
                     f"{int(r.min_group)} |")
        L += [""]
        if first is not None and len(inside):
            L += [f"The record sits inside the drift band at a radius of "
                  f"{inside.radius_km.max():.0f} km and leaves it by "
                  f"{first.radius_km:.0f} km.", ""]
        L += ["**The smallest-group column is the limit on this table.** Once a "
              "partition isolates a",
              "single assemblage, that group matches its own profile exactly and "
              "contributes to",
              "between-group variance for a reason that is arithmetic rather than "
              "archaeological.",
              "Rows whose smallest group is 1 or 2 are reported but should not be "
              "read as evidence",
              "of structure at that scale.", ""]
    L += [f"Figure written to {png.name} and its siblings.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"phase boundaries {obs_med:.5f}/km; matched boundaries median "
          f"{np.median(alt_med):.5f}/km -> {t74.ordinal(pct_alt)} percentile "
          f"(whole-field comparison, confounded: {t74.ordinal(pct_field)})"
          if len(on_b) else "no boundaries")
    print(f"wrote {png} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
