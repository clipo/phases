#!/usr/bin/env python3
"""Figure 6: rim attributes beside decoration on the 19 sites with rim tabulations.

Drawn from the committed tables in `output/findings/rim_r1.md` (produced by
`analyses/rims/rim_r1.py --readings`; the per-cell CSVs it read are gitignored
and absent here). Nothing is recomputed. Section A of that file supplies the
site-weighted moment F of each attribute and of the three common decorated
classes on the same 19 sites, with the neutral band at the primary operating
point, window T1 (panel A); section B the posterior rank correlation of each
site frequency with latitude under a Beta(1/2, 1/2) prior per site (panel B);
section C the rank correlation between rim and decoration pairwise differences
after linear adjustment on ranks for distance (panel C). The values the text
quotes are asserted before drawing.

Usage:
    .venv/bin/python analyses/106_fig6_rim_attributes.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))
fs = importlib.import_module("figstyle")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

FIG = "fig6_rim_attributes"
SRC = ROOT / "output" / "findings" / "rim_r1.md"
RIMS = [("interior bevel", "interior beveling"), ("exterior notching", "exterior notching"),
        ("lip flat", "lip form (flat)"), ("paste Bell (of N)", "Bell paste")]
DECOR = [("Parkin_Punctated, current counts", "Parkin Punctated"),
         ("Barton/Kent/MPI, current counts", "Barton, Kent, Mound Place Inc."),
         ("Painted, current counts", "painted types")]
POST = re.compile(r"posterior ([-0-9.]+) \[([-0-9.]+), ([-0-9.]+)\]")


def section(text: str, start: str, end: str) -> str:
    return text.split(start, 1)[1].split(end, 1)[0]


def rows_of(block: str) -> list[list[str]]:
    out = []
    for line in block.splitlines():
        if line.startswith("|") and not line.startswith("|---") and not line.startswith("| trait"):
            out.append([c.strip() for c in line.strip().strip("|").split("|")])
    return out


def main() -> int:
    text = SRC.read_text()
    secA = section(text, "### primary point, window T1", "### primary point, window 8LR")
    secB = section(text, "## B. Association with geography", "## C. ")
    secC = section(text, "## C. Rim and decoration pair differences", "## D. ")
    keys = [k for k, _ in RIMS + DECOR]

    F, band = {}, {}
    for cells in rows_of(secA):
        key = cells[0].replace(" (read separately)", "")
        if key in keys:
            F[key] = float(cells[1])
            m = re.match(r"\[([-0-9.]+), ([-0-9.]+)\] n (\d+)", cells[2])
            band[key] = (float(m.group(1)), float(m.group(2)))
    assert set(F) == set(keys), sorted(F)
    assert [round(F[k], 3) for k, _ in RIMS] == [0.105, 0.069, 0.031, 0.266], F

    lat = {}
    for cells in rows_of(secB):
        if cells[0] in keys:
            m = POST.search(cells[2])
            lat[cells[0]] = tuple(float(g) for g in m.groups())
    assert set(lat) == set(keys), sorted(lat)
    assert [round(lat[k][0], 2) for k in ("interior bevel", "exterior notching", "Parkin_Punctated, current counts",
                                          "Painted, current counts")] == [-0.58, -0.39, 0.65, -0.74], lat

    adj = {}
    for cells in rows_of(secC):
        if cells[0] in [k for k, _ in RIMS]:
            cur = re.search(r"adjusted [-0-9.]+; " + POST.pattern, cells[1])
            pfg = re.search(r"adjusted [-0-9.]+; " + POST.pattern, cells[2])
            adj[cells[0]] = (tuple(float(g) for g in cur.groups()), tuple(float(g) for g in pfg.groups()))
    assert len(adj) == 4, adj
    meds = [v[i][0] for v in adj.values() for i in (0, 1)]
    assert -0.075 <= min(meds) and max(meds) <= 0.145, meds

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.9), gridspec_kw=dict(width_ratios=[1.0, 1.0, 0.8], wspace=0.6))
    order = RIMS + DECOR
    ypos = [0, 1, 2, 3, 4.6, 5.6, 6.6]
    axA, axB, axC = axes
    for y, (k, lab) in zip(ypos, order):
        axA.plot(band[k], [y, y], color="0.75", lw=6, solid_capstyle="butt", zorder=1)
        axA.scatter([F[k]], [y], marker="o" if k in dict(RIMS) else "s", s=36, facecolor="black",
                    edgecolor="white", lw=0.6, zorder=5)
        m, lo, hi = lat[k]
        axB.errorbar([m], [y], xerr=[[m - lo], [hi - m]], fmt="o" if k in dict(RIMS) else "s", color="black",
                     ecolor="black", elinewidth=1.0, capsize=2, markersize=5)
    for ax in (axA, axB):
        ax.set_yticks(ypos)
        ax.set_ylim(7.3, -0.7)
        ax.axhline(3.8, color="0.85", lw=0.8)
    axA.set_yticklabels([lab for _, lab in order], fontsize=7)
    axB.set_yticklabels([])
    axA.set_xlabel("site-weighted moment F (19 sites)")
    axA.set_xlim(-0.01, 0.3)
    axB.axvline(0, color="0.6", lw=0.8)
    axB.set_xlabel("rank correlation with latitude\n(positive: higher to the north)")
    axB.set_xlim(-0.9, 0.9)
    for y, (k, lab) in zip(range(4), RIMS):
        (mc, loc, hic), (mp, lop, hip) = adj[k]
        axC.errorbar([mc], [y - 0.15], xerr=[[mc - loc], [hic - mc]], fmt="o", color="black", ecolor="black",
                     elinewidth=1.0, capsize=2, markersize=5)
        axC.errorbar([mp], [y + 0.15], xerr=[[mp - lop], [hip - mp]], fmt="o", color="black", ecolor="black",
                     markerfacecolor="white", elinewidth=1.0, capsize=2, markersize=5)
    axC.set_yticks(range(4))
    axC.set_yticklabels([lab for _, lab in RIMS], fontsize=7)
    axC.set_ylim(3.6, -0.6)
    axC.axvline(0, color="0.6", lw=0.8)
    axC.set_xlim(-0.3, 0.3)
    axC.set_xlabel("rank correlation, rim against\ndecoration pair differences,\nadjusted for distance")
    axC.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor="black", markeredgecolor="black",
                               ms=5, label="against current counts"),
                        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black",
                               ms=5, label="against 1940 collections")],
               fontsize=6, loc="lower left", bbox_to_anchor=(0.0, 1.02), frameon=False)
    axA.legend(handles=[Line2D([0], [0], color="0.75", lw=6, label="neutral band, primary point, window T1"),
                        Line2D([0], [0], marker="o", color="none", markerfacecolor="black", ms=5, label="rim attribute"),
                        Line2D([0], [0], marker="s", color="none", markerfacecolor="black", ms=5,
                               label="decorated class, same 19 sites")],
               fontsize=6, loc="lower left", bbox_to_anchor=(0.0, 1.02), frameon=False)
    for ax, letter in zip(axes, "ABC"):
        fs.despine(ax)
        fs.panel_label(ax, letter, y=1.20)  # above the legends anchored at 1.02
    fig.subplots_adjust(left=0.2, right=0.97, top=0.80, bottom=0.25, wspace=0.6)
    png = fs.save_all(fig, FIG, close=True)
    print(f"written {png}; F {[F[k] for k, _ in RIMS]}; lat {[lat[k][0] for k, _ in RIMS]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
