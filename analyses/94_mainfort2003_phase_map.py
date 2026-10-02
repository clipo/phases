"""94_mainfort2003_phase_map.py — the phases as Mainfort (2003) assigned them.

The map for the contrast case (analyses/83_mainfort_replication.py). It is the
former Figure 1 redrawn on Mainfort's (2003: Table 1) site list and his own
phase assignments instead of the 1996 assignments the primary analysis uses:
39 sites in seven phases, with Smith's (1990) Horseshoe Lake phase split out
of Walls.

The territories are built exactly as `36_canonical_phase_map.py` builds them
(every point of the map goes to its nearest site, the areas are merged by
phase, clipped to a loose envelope around the sites and shrunk so neighbors
show a gap), from all 39 sites. They are our construction from his
assignments, not lines he drew (rule 13). Because his scheme follows the
ordination rather than the map, a phase's territory can come out in more than
one piece; that is drawn as it falls.

Filled markers are the sites the contrast case tests: at least MIN_DECORATED
decorated sherds, in a phase that keeps at least two such sites (the same two
rules as analysis 83, whose constants are imported so the map cannot disagree
with the test). Open markers are his other sites.

Writes figures/fig18_mainfort2003_phases.png.

Usage: PYTHONPATH=src .venv/bin/python analyses/94_mainfort2003_phase_map.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import make_figures as mf  # noqa: E402
from mls_emergence.dataio.matrix import MAINFORT_SITES, read_mainfort_replication  # noqa: E402
m36 = importlib.import_module("36_canonical_phase_map")
m83 = importlib.import_module("83_mainfort_replication")

OUT_FIG = ROOT / "figures" / "fig18_mainfort2003_phases.png"
MIN_DECORATED = 100   # the contrast case's minimum (83's --min default)

PHASES = ["Parkin", "Nodena", "Jones Bayou", "Tipton", "Kent", "Walls",
          "Horseshoe Lake"]
# Grays and marker shapes follow 36 for the phases the two maps share;
# Horseshoe Lake takes the gray and the marker Parchman had there.
PHASE_GRAY = {**{p: m36.PHASE_GRAY[p] for p in PHASES[:-1]},
              "Horseshoe Lake": m36.PHASE_GRAY["Parchman"]}
PHASE_MARKER = {**{p: m36.PHASE_MARKER[p] for p in PHASES[:-1]},
                "Horseshoe Lake": m36.PHASE_MARKER["Parchman"]}


def main():
    sites = pd.read_csv(MAINFORT_SITES, comment="#").set_index("site")
    unknown = sorted(set(sites["phase_mainfort2003"]) - set(PHASES))
    if unknown:
        raise SystemExit(f"phases in {MAINFORT_SITES.name} this map does not draw: {unknown}")
    names = list(sites.index)
    lab_idx = np.array([PHASES.index(p) for p in sites["phase_mainfort2003"]])

    # The tested set, by analysis 83's two rules.
    _, _, ph_min = read_mainfort_replication(MIN_DECORATED)
    sizes = ph_min.value_counts()
    kept = ph_min[~ph_min.isin(sizes[sizes < m83.MIN_PHASE_MEMBERS].index)]
    tested = np.array([nm in kept.index for nm in names])

    E, Nm = m36.to_utm(sites["longitude"].to_numpy(float),
                       sites["latitude"].to_numpy(float))
    margin = 10_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)

    fig, ax = plt.subplots(figsize=(5.8, 6.8))
    m36.draw_base(ax, ext, {"AR": (0.43, 0.71), "TN": (0.88, 0.40), "MS": (0.85, 0.17)})
    # The river name is placed automatically and lands on one of his sites; move
    # it 9 km upstream, where the river is clear of markers.
    m36.move_label(ax, "Mississippi R.", dn=9_000.0)

    territory = m36.voronoi_territories(E, Nm, lab_idx, PHASES, ext)
    for ph, poly in territory.items():
        m36.fill_polygon(ax, poly, facecolor=str(PHASE_GRAY[ph]), edgecolor="0.2",
                         linewidth=0.8, zorder=1.5)

    for k, ph in enumerate(PHASES):
        m = (lab_idx == k) & tested
        ax.scatter(E[m], Nm[m], s=26, c="black", marker=PHASE_MARKER[ph],
                   edgecolor="white", linewidth=0.4, zorder=10)
        m = (lab_idx == k) & ~tested
        ax.scatter(E[m], Nm[m], s=22, facecolor="white", marker=PHASE_MARKER[ph],
                   edgecolor="black", linewidth=0.6, zorder=10)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])

    # A phase whose territory is in pieces is labeled in its largest piece.
    largest = {ph: (max(poly.geoms, key=lambda g: g.area)
                    if poly.geom_type == "MultiPolygon" else poly)
               for ph, poly in territory.items()}
    anchors = {ph: (E[lab_idx == k].mean(), Nm[lab_idx == k].mean())
               for k, ph in enumerate(PHASES)}
    # Horseshoe Lake's territory is too small to hold its own name.
    m36.place_phase_labels(fig, ax, ext, E, Nm, anchors, largest,
                           spill=("Horseshoe Lake",))

    pk = names.index("Parkin")
    ax.scatter([E[pk]], [Nm[pk]], marker="*", s=210, c="white", edgecolor="black",
               linewidth=0.8, zorder=14)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])

    handles = [Line2D([0], [0], marker=PHASE_MARKER[ph], color="none",
                      markerfacecolor="black", markeredgecolor="white",
                      markeredgewidth=0.4, markersize=6, label=ph) for ph in PHASES]
    handles += [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=0.6, markersize=5,
               label="Open: not tested"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
               markeredgecolor="black", markersize=12, label="Parkin site"),
    ]
    ax.legend(handles=handles, fontsize=6.4, loc="lower center",
              bbox_to_anchor=(0.5, 1.005), ncol=3, framealpha=0.0,
              handletextpad=0.4, columnspacing=1.2, labelspacing=0.4,
              borderpad=0.2)

    # The lower right, where 36 puts the inset, is Walls territory here. The
    # inset goes to the empty northwest corner and the scale bar to the southeast.
    m36.add_scale_bar(ax, ext, x0=0.72, y0=0.05)
    m36.add_locator_inset(fig, ax, ext, rect=(0.015, 0.74, 0.34, 0.27))

    mf.save_all(fig, OUT_FIG)
    plt.close(fig)
    pieces = {ph: (len(poly.geoms) if poly.geom_type == "MultiPolygon" else 1)
              for ph, poly in territory.items()}
    print(f"sites: {len(names)} in {len(set(lab_idx))} phases;",
          {p: int((lab_idx == k).sum()) for k, p in enumerate(PHASES)})
    print(f"tested (>= {MIN_DECORATED} decorated, phase keeps >= "
          f"{m83.MIN_PHASE_MEMBERS} such sites): {int(tested.sum())};",
          {p: int(((lab_idx == k) & tested).sum()) for k, p in enumerate(PHASES)})
    print("territory pieces:", pieces)
    print(f"wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
