"""93_phillips1970_phase_map.py — the study area and the phases as Phillips drew them.

Manuscript Figure 1. The four late-period phase areas of Phillips (1970),
Nodena, Parkin, Walls and Kent, as redrawn by Lipo (2001: Figure 2.6), over the
same hydrography as the other maps, with the decorated-ceramic assemblages of
the analysis matrix that fall inside the frame.

The outlines are the source's own lines, not a construction of ours: they are
read from the vector artwork of Lipo's figure and georeferenced from its
quadrangle frame by `scripts/build_phillips1970_outlines.py`, which writes
`data/raw/phillips1970_phase_outlines.csv` and records the check on the
georeference in that file's header. Unlike the territories of
`36_canonical_phase_map.py` (every point assigned to its nearest assemblage,
merged by phase), these areas do not tile the map: they leave ground between
them that belongs to no phase, and an assemblage can fall outside all four.

The frame is the source figure's, Lower Mississippi Survey columns N to Q by
rows 10 to 14 (longitude -90.75 to -89.75, latitude 34.50 to 35.75), with a
small margin so the outlines that run along its west edge are not drawn on the
axis line. The two assemblages of the matrix outside that frame (Parchman and
Salomon, south of it) are not drawn.

Also writes output/findings/phillips1970_outline_membership.md: for each
assemblage drawn, the phase area it lies in (the phase the analysis assigns
it since 2026-10-01), whether that is by containment or by nearness, and the
Mainfort (1996: Figure 1) label the pipeline used before.

Usage: PYTHONPATH=src .venv/bin/python analyses/93_phillips1970_phase_map.py
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
from matplotlib.patches import Patch  # noqa: E402
import make_figures as mf  # noqa: E402
m35 = importlib.import_module("35_basin_pullout")
m36 = importlib.import_module("36_canonical_phase_map")

OUT_FIG = ROOT / "figures" / "fig1_phillips1970_phases.png"
OUT_MD = ROOT / "output" / "findings" / "phillips1970_outline_membership.md"

PHASES = m36.PHILLIPS_PHASES
# The source figure's frame: LMS columns N to Q by rows 10 to 14.
LON_WEST, LON_EAST, LAT_SOUTH, LAT_NORTH = m36.PHILLIPS_FRAME
MARGIN_M = 4_000.0
# The largest distance from a known site to its plotted dot in the source
# figure is 1.21 km (the outline file's header); a site closer than this to an
# outline is not securely on either side of it.
NEAR_LINE_KM = 1.25
FILL_GRAY = "0.90"


def main():
    outline = m36.load_phillips_outlines()

    counts_df, coords_df = m35.load_full()
    names = [str(x) for x in counts_df.index]
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    analyzed = np.array([nm in mf._basin_members("curated") for nm in names])
    E, Nm = m36.to_utm(coords[:, 1], coords[:, 0])

    fx, fy = m36.to_utm(np.array([LON_WEST, LON_EAST, LON_WEST, LON_EAST]),
                        np.array([LAT_SOUTH, LAT_SOUTH, LAT_NORTH, LAT_NORTH]))
    ext = (fx.min() - MARGIN_M, fx.max() + MARGIN_M,
           fy.min() - MARGIN_M, fy.max() + MARGIN_M)
    shown = (E > ext[0]) & (E < ext[1]) & (Nm > ext[2]) & (Nm < ext[3])
    if not shown[analyzed].all():
        raise SystemExit("an analyzed assemblage falls outside the frame: "
                         + ", ".join(np.array(names)[analyzed & ~shown]))

    fig, ax = plt.subplots(figsize=(5.2, 7.2))
    m36.draw_base(ax, ext, {"AR": (0.47, 0.60), "TN": (0.88, 0.60), "MS": (0.86, 0.36)})

    # Phase areas: one light gray under the hydrography (the four do not
    # overlap, so the labels carry the identity), with the outline on top.
    for ph in PHASES:
        m36.fill_polygon(ax, outline[ph], facecolor=FILL_GRAY, edgecolor="none",
                         zorder=1.5)
        x, y = outline[ph].exterior.xy
        ax.plot(x, y, color="0.15", linewidth=1.2, zorder=6,
                solid_joinstyle="round")

    # Assemblages: the analyzed set filled, the rest of the matrix open.
    m = shown & ~analyzed
    ax.scatter(E[m], Nm[m], s=20, facecolor="white", edgecolor="black",
               linewidth=0.7, zorder=10)
    m = shown & analyzed
    ax.scatter(E[m], Nm[m], s=24, c="black", edgecolor="white", linewidth=0.4,
               zorder=10)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])

    anchors = {ph: (outline[ph].representative_point().x,
                    outline[ph].representative_point().y) for ph in PHASES}
    m36.place_phase_labels(fig, ax, ext, E[shown], Nm[shown], anchors, outline)

    pk = names.index("Parkin")
    ax.scatter([E[pk]], [Nm[pk]], marker="*", s=210, c="white", edgecolor="black",
               linewidth=0.8, zorder=14)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])

    n_an, n_other = int((shown & analyzed).sum()), int((shown & ~analyzed).sum())
    handles = [
        Patch(facecolor=FILL_GRAY, edgecolor="0.15", linewidth=1.0,
              label="Phase area (Phillips 1970)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="black",
               markeredgecolor="white", markeredgewidth=0.4, markersize=6,
               label="Analyzed assemblage"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=0.7, markersize=5,
               label="Other assemblage"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
               markeredgecolor="black", markersize=12, label="Parkin site"),
    ]
    ax.legend(handles=handles, fontsize=6.4, loc="lower center",
              bbox_to_anchor=(0.5, 1.005), ncol=2, framealpha=0.0,
              handletextpad=0.4, columnspacing=1.2, labelspacing=0.4,
              borderpad=0.2)

    m36.add_scale_bar(ax, ext, x0=0.05, y0=0.955)
    m36.add_locator_inset(fig, ax, ext, rect=(0.63, 0.015, 0.355, 0.25))

    mf.save_all(fig, OUT_FIG)
    plt.close(fig)

    # Membership: which area each assemblage lies in, and what the superseded
    # Mainfort (1996) labels said, for the record of the 2026-10-01 switch.
    labels, derived = m36.assign_primary_phases(names, coords)
    old, _ = m36.assign_phases_by_territory(names, coords)
    km = m36.phillips_outline_distance(coords, labels)
    rows = []
    for i in np.where(shown)[0]:
        rows.append((names[i], "analyzed" if analyzed[i] else "not analyzed", labels[i],
                     "nearest area" if derived[i] else "inside the area",
                     f"{km[i]:.1f}" if derived[i] else "0",
                     old[i], "" if old[i] == labels[i] else "changed"))
    tab = pd.DataFrame(rows, columns=["assemblage", "set", "Phillips phase", "assigned by",
                                      "km outside the line", "Mainfort 1996 label (superseded)",
                                      "changed"])
    an = tab[tab["set"] == "analyzed"]
    near = an[an["assigned by"] == "nearest area"]
    moved = an[an["changed"] == "changed"]
    area = {ph: outline[ph].area / 1e6 for ph in PHASES}
    lines = [
        "# The assemblages and Phillips's (1970) phase areas",
        "",
        "Written by `analyses/93_phillips1970_phase_map.py`. The areas are those of",
        "`data/raw/phillips1970_phase_outlines.csv` (Lipo 2001: Figure 2.6, after Phillips 1970).",
        "The analysis assigns each assemblage the area it lies in, or the nearest when it lies in",
        "none (`36_canonical_phase_map.assign_primary_phases`). Containment is point-in-polygon in",
        f"UTM 15N; the areas carry the source's drafting error (see the CSV header), so an assemblage",
        f"within {NEAR_LINE_KM:g} km of a line is not securely on either side of it. The last two columns",
        "record the Mainfort (1996: Figure 1) labels the pipeline used until 2026-10-01.",
        "",
        f"- Drawn in Figure 1: {n_an} analyzed assemblages and {n_other} other assemblages of the",
        f"  {len(names)}-assemblage matrix; {len(names) - n_an - n_other} fall outside the frame.",
        f"- Of the {len(an)} analyzed assemblages, {len(an) - len(near)} lie inside a phase area and {len(near)} take the",
        "  nearest: " + ("; ".join(f"{r[0]} ({r[2]}, {r[4]} km outside the line)"
                                   for r in near.itertuples(index=False)) or "none") + ".",
        "- By phase: " + ", ".join(f"{ph} {int((an['Phillips phase'] == ph).sum())}" for ph in PHASES) + ".",
        f"- {len(moved)} analyzed assemblages carry a different phase than under the 1996 labels: "
        + ("; ".join(f"{r[0]} ({r[5]} to {r[2]})" for r in moved.itertuples(index=False)) or "none") + ".",
        "- Area sizes (km2): " + ", ".join(f"{ph} {area[ph]:.0f}" for ph in PHASES) + ". No two areas overlap.",
        "",
        "| " + " | ".join(tab.columns) + " |",
        "|" + "---|" * len(tab.columns),
        *("| " + " | ".join(r) + " |" for r in tab.itertuples(index=False)),
        "",
    ]
    OUT_MD.write_text("\n".join(lines))
    print("\n".join(lines[10:18]))
    print(tab.to_string(index=False))
    print(f"wrote {OUT_FIG}\nwrote {OUT_MD}")


if __name__ == "__main__":
    main()
