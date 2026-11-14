#!/usr/bin/env python3
"""Schematic of the three models: spatial field, latent-unit and barrier.

What is real and what is illustrative. The 28 site locations are the analyzed
assemblages (centered kilometers, `dataio.basin.km_xy`). In panel C the two
edge placements are the ones the barrier model was tested with: the
two-cluster division by k-means on coordinates
(`probe.placement.kmeans2_regions`) and the straight line that best agrees
with Phillips's Parkin-versus-rest labels among the declared candidates
(`probe.placement.candidate_lines`, `best_line_for`), recomputed here by the
same calls as `analyses/probe_latent_groups/a3_contracts.py`. Everything else
is drawn to show what each model assumes and is NOT a fitted value: the smooth
surface (a fixed analytic function), the assignment of sites to three units in
panel B (one draw at random from a fixed seed), and the size of the step in
panel C.

The script asserts the counts the caption states: 30 candidate lines, and
sides of 11 and 17 sites for the k-means division.

Usage:
    .venv/bin/python analyses/103_fig3_model_schematic.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))
fs = importlib.import_module("figstyle")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from mls_emergence.dataio.basin import km_xy  # noqa: E402
from mls_emergence.probe import placement as pl  # noqa: E402

FIG = "fig3_model_schematic"
UNIT_SEED = 3            # the illustrative assignment of panel B
STEP = 0.9               # illustrative step across the edge in panel C, in the surface's arbitrary units
PAD = 12.0               # km of margin around the sites


def surface(x, y):
    """A smooth illustrative surface in arbitrary units: a tilt to the north-east with gentle curvature."""
    return 0.018 * y + 0.010 * x + 0.35 * np.sin(x / 38.0) * np.cos(y / 55.0)


def placements():
    """Site coordinates (centered km), the k-means division, the chosen line, and the number of candidate lines."""
    mf = importlib.import_module("make_figures")
    m36 = importlib.import_module("36_canonical_phase_map")
    counts, coords = mf._load_curated()
    latlon = coords.to_numpy(float)
    if len(latlon) != 28:
        raise ValueError(f"expected the 28 analyzed assemblages, found {len(latlon)}")
    labels, _ = m36.assign_primary_phases([str(i) for i in counts.index], latlon)
    labels = np.asarray(labels)
    g1 = pl.kmeans2_regions(latlon)
    cands = pl.candidate_lines(latlon)
    line = pl.best_line_for(np.where(labels == "Parkin", 1, 2), cands)
    return km_xy(latlon), g1, line, len(cands), labels


def draw_surface(ax, lim, step_line=None):
    gx, gy = np.meshgrid(np.linspace(lim[0], lim[1], 240), np.linspace(lim[2], lim[3], 240))
    z = surface(gx, gy)
    if step_line is not None:
        a = np.deg2rad(step_line["angle"])
        z = z + STEP * ((gx * np.cos(a) + gy * np.sin(a)) >= step_line["offset"])
    ax.contourf(gx, gy, z, levels=np.linspace(-1.6, 2.6, 15), cmap="Greys", vmin=-3.0, vmax=6.0)


def frame(ax, lim, letter):
    ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.6); s.set_color("0.4")
    fs.panel_label(ax, letter)
    x0, y0 = lim[0] + 5, lim[3] - 9                      # 20 km scale bar, upper left, clear of the sites
    ax.plot([x0, x0 + 20], [y0, y0], color="0", lw=1.5, solid_capstyle="butt")
    ax.text(x0 + 10, y0 + 1.8, "20 km", ha="center", va="bottom", fontsize=7)


def main() -> int:
    p, g1, line, n_cands, labels = placements()
    sides = sorted(int((g1 == k).sum()) for k in (1, 2))
    if n_cands != 30 or sides != [11, 17]:
        raise AssertionError(f"caption states 30 candidate lines and sides of 11 and 17; found {n_cands} and {sides}")
    parkin = labels == "Parkin"
    small = 1 if (g1 == 1).sum() < (g1 == 2).sum() else 2
    n_parkin_in_small = int((parkin & (g1 == small)).sum())
    if (int(parkin.sum()), n_parkin_in_small) != (12, 11):
        raise AssertionError(f"caption states eleven of twelve Parkin-area sites on one side; found {n_parkin_in_small} of {int(parkin.sum())}")
    line_sides = sorted(int((line["regions"] == k).sum()) for k in (1, 2))
    lim = (p[:, 0].min() - PAD, p[:, 0].max() + PAD, p[:, 1].min() - PAD, p[:, 1].max() + PAD)

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.5))
    fig.subplots_adjust(wspace=0.06, bottom=0.13, top=0.95)

    # (A) spatial field model: a smooth surface and nothing else
    ax = axes[0]; draw_surface(ax, lim); frame(ax, lim, "A")
    ax.scatter(p[:, 0], p[:, 1], s=20, c="0", edgecolors="1", linewidths=0.5, zorder=3)

    # (B) latent-unit model: the same surface, plus three units whose membership is not tied to location
    ax = axes[1]; draw_surface(ax, lim); frame(ax, lim, "B")
    unit = np.random.default_rng(UNIT_SEED).integers(0, 3, len(p))
    for k, (mk, fc) in enumerate((("o", "0"), ("s", "1"), ("^", "0.55"))):
        m = unit == k
        ax.scatter(p[m, 0], p[m, 1], s=30 if mk != "^" else 36, marker=mk, c=fc, edgecolors="0", linewidths=0.7, zorder=3)

    # (C) barrier model: the surface with a step at a fixed edge; both tested placements
    ax = axes[2]; draw_surface(ax, lim, step_line=line); frame(ax, lim, "C")
    a = np.deg2rad(line["angle"]); c, s = np.cos(a), np.sin(a)
    t = np.linspace(-400, 400, 2)                         # the line x cos a + y sin a = offset
    ax.plot(line["offset"] * c - t * s, line["offset"] * s + t * c, color="0", lw=1.3, ls=(0, (5, 3)), zorder=2)
    for k, fc in ((small, "0"), (3 - small, "1")):
        m = g1 == k
        ax.scatter(p[m, 0], p[m, 1], s=24, c=fc, edgecolors="0", linewidths=0.7, zorder=3)

    handles = [
        Line2D([], [], marker="o", ls="", mfc="0", mec="1", ms=5, label="A: site; shading is a smooth surface"),
        Line2D([], [], marker="s", ls="", mfc="1", mec="0", ms=5, label="B: marker shape is the unit (one assignment of many)"),
        Line2D([], [], marker="o", ls="", mfc="0", mec="0", ms=5, label="C: filled and open are the two k-means clusters"),
        Line2D([], [], color="0", lw=1.3, ls=(0, (5, 3)), label="C: straight-line edge, with a step in the surface"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False, fontsize=7, handletextpad=0.5, columnspacing=1.5)
    out = fs.save(fig, FIG)
    print(f"wrote {out}; candidate lines {n_cands}; k-means sides {sides}; Parkin-area sites in the smaller cluster "
          f"{n_parkin_in_small} of {int(parkin.sum())}; line angle {line['angle']:g}, offset {line['offset']:.1f} km, sides {line_sides}, "
          f"agreement with the Parkin split {line['agreement']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
