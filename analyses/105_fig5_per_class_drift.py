#!/usr/bin/env python3
"""Figure 5: the decorated classes against neutral copying with regional drift.

Drawn from the committed tables in `output/findings/per_type_drift_sensitivity.md`
(primary operating point; produced by `analyses/per_type/drift_sensitivity_figure.py`
from per-record CSVs that are gitignored and absent here). Nothing is simulated:
the figure re-plots the coverage table (neutral median, 10th to 90th percentile
band, the record) and the planted-marker tables (planted medians and the count
of planted records outside the neutral band) exactly as committed, and asserts
the values the text quotes before drawing.

(A) Site-weighted moment F of the three common classes as the amplitude F_W of
regional drift in the innovation target rises, with planted markers.
(B) Fraction of planted records outside the neutral band of the same class and
amplitude, beside the neutral cells' own fraction outside their band.
(C) Sites present for the seven rare classes at F_W = 0: neutral median and
band, the record, and the planted Walls Engraved occupancy at a multiplier of 3.

Usage:
    .venv/bin/python analyses/105_fig5_per_class_drift.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))
fs = importlib.import_module("figstyle")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

FIG = "fig5_per_class_drift"
SRC = ROOT / "output" / "findings" / "per_type_drift_sensitivity.md"
FW = ["0", "0.01", "0.03", "0.1", "0.3"]
COMMON = ["Parkin_Punctated", "Barton/Kent/MPI", "Painted"]
RARE = ["Fortune_Noded", "Ranch_Incised", "Walls_Engraved", "Wallace_Incised",
        "Rhodes_Incised", "Vernon_Paul_Applique", "Hull_Engraved"]
LABEL = {"Parkin_Punctated": "Parkin Punctated", "Barton/Kent/MPI": "Barton, Kent, Mound Place Inc.",
         "Painted": "painted types", "Fortune_Noded": "Fortune\nNoded", "Ranch_Incised": "Ranch\nIncised",
         "Walls_Engraved": "Walls\nEngraved", "Wallace_Incised": "Wallace\nIncised",
         "Rhodes_Incised": "Rhodes\nIncised", "Vernon_Paul_Applique": "Vernon Paul\nAppliqué",
         "Hull_Engraved": "Hull\nEngraved"}
CELL = re.compile(r"([-0-9.]+) \[([-0-9.]+), ([-0-9.]+)\] n (\d+): \*\*(\w+)\*\*(?:; lost (\d+))?(?:; absent (\d+))?")
MROW = re.compile(r"\| ([0-9.]+) \| ([0-9.]+) \| (\d+) of 50 \((\d+) undefined\) \| (\d+) of 50 \((\d+) undefined\) "
                  r"\| ([-0-9.]+|[a-z]+) \((eligible|ineligible)\) \| \[([-0-9.]+), ([-0-9.]+)\] n (\d+) \|")


def parse() -> tuple[dict, dict]:
    text = SRC.read_text()
    prim = text.split("## primary operating point", 1)[1].split("## secondary operating point", 1)[0]
    coverage = {}
    for line in prim.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 7 and cells[0] in COMMON + RARE:
            parsed = [CELL.match(c) for c in cells[2:]]
            assert all(parsed), line
            coverage[cells[0]] = {"record": float(cells[1]),
                                  "cells": [dict(med=float(m.group(1)), lo=float(m.group(2)), hi=float(m.group(3)),
                                                 n=int(m.group(4)), status=m.group(5)) for m in parsed]}
    assert set(coverage) == set(COMMON + RARE), sorted(coverage)
    markers = {}
    for block in prim.split("### ")[1:]:
        head = block.splitlines()[0]
        m = re.match(r"(\S+), (large|small) cluster \((\d+) sites\), statistic: (.+)", head)
        if not m:
            continue
        rows = {}
        for line in block.splitlines():
            r = MROW.match(line.strip())
            if r:
                fw, mult = r.group(1), float(r.group(2))
                med = float(r.group(7)) if r.group(8) == "eligible" and re.match(r"[-0-9.]+$", r.group(7)) else np.nan
                rows[(fw, mult)] = dict(planted_out=int(r.group(3)), neutral_out=int(r.group(5)), planted_med=med,
                                        lo=float(r.group(9)), hi=float(r.group(10)))
        assert len(rows) == 15, (head, len(rows))
        markers[(m.group(1), m.group(2))] = rows
    assert len(markers) == 8, sorted(markers)
    return coverage, markers


def main() -> int:
    cov, mk = parse()
    rec = [cov[c]["record"] for c in COMMON]
    assert [round(v, 3) for v in rec] == [0.165, 0.117, 0.199], rec
    assert [round(cov[c]["cells"][4]["med"], 3) for c in COMMON] == [0.035, 0.032, 0.029]
    assert cov["Fortune_Noded"]["cells"][0]["status"] == "outside" and cov["Hull_Engraved"]["cells"][0]["status"] == "outside"

    x = np.arange(len(FW))
    shade = {1.5: "0.72", 2.0: "0.45", 3.0: "0.0"}
    fig = plt.figure(figsize=(7.0, 8.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.25, 1.0, 1.1], hspace=0.5, wspace=0.32)
    axA = [fig.add_subplot(gs[0, j]) for j in range(3)]
    axB = [fig.add_subplot(gs[1, j]) for j in range(3)]
    axC = fig.add_subplot(gs[2, :])

    for j, c in enumerate(COMMON):
        ax = axA[j]
        cells = cov[c]["cells"]
        ax.fill_between(x, [d["lo"] for d in cells], [d["hi"] for d in cells], color="0.88", lw=0, zorder=1)
        ax.plot(x, [d["med"] for d in cells], color="0.35", lw=1.3, zorder=2)
        ax.axhline(cov[c]["record"], color="black", ls="--", lw=1.2, zorder=4)
        for clus, ls in (("large", "-"), ("small", ":")):
            rows = mk[(c, clus)]
            for mult in (1.5, 2.0, 3.0):
                ys = [rows[(fw, mult)]["planted_med"] for fw in FW]
                ax.plot(x, ys, color=shade[mult], ls=ls, lw=0.9, marker="o", ms=2.6, zorder=3)
        ax.set_yscale("log")
        ax.set_ylim(0.004, 0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(FW, fontsize=7)
        ax.set_title(LABEL[c], fontsize=8, pad=3)
        if j == 0:
            ax.set_ylabel("site-weighted moment F")
        fs.despine(ax)
        if j == 0:
            fs.panel_label(ax, "A")

    for j, c in enumerate(COMMON):
        ax = axB[j]
        for clus, ls in (("large", "-"), ("small", ":")):
            rows = mk[(c, clus)]
            for mult in (1.5, 2.0, 3.0):
                ax.plot(x, [rows[(fw, mult)]["planted_out"] / 50 for fw in FW], color=shade[mult], ls=ls,
                        lw=0.9, marker="o", ms=2.6)
        ax.plot(x, [mk[(c, "large")][(fw, 2.0)]["neutral_out"] / 50 for fw in FW], color="0.35", ls="-.", lw=1.0)
        ax.set_ylim(0, 1.04)
        ax.set_xticks(x)
        ax.set_xticklabels(FW, fontsize=7)
        ax.set_xlabel("regional drift amplitude $F_W$", fontsize=7.5)
        if j == 0:
            ax.set_ylabel("fraction of records outside\nthe neutral band")
            fs.panel_label(ax, "B")
        fs.despine(ax)

    xr = np.arange(len(RARE))
    for i, c in enumerate(RARE):
        d = cov[c]["cells"][0]
        axC.plot([i, i], [d["lo"], d["hi"]], color="0.55", lw=5, solid_capstyle="butt", zorder=1)
        axC.plot([i - 0.22, i + 0.22], [d["med"], d["med"]], color="0.2", lw=1.4, zorder=2)
        axC.scatter([i], [cov[c]["record"]], marker="D", s=42, facecolor="black", edgecolor="white", lw=0.6, zorder=5)
    iw = RARE.index("Walls_Engraved")
    for clus, dx, mkr in (("large", -0.12, "^"), ("small", 0.12, "v")):
        ymed = mk[("Walls_Engraved", clus)][("0", 3.0)]["planted_med"]
        axC.scatter([iw + dx], [ymed], marker=mkr, s=46, facecolor="white", edgecolor="black", lw=0.9, zorder=6)
    axC.set_xticks(xr)
    axC.set_xticklabels([LABEL[c] for c in RARE], fontsize=7)
    axC.set_ylabel("sites present, of 28")
    axC.set_ylim(0, 28)
    axC.set_xlim(-0.6, len(RARE) - 0.4)
    fs.despine(axC)
    fs.panel_label(axC, "C")

    handles = [Line2D([0], [0], color="0.35", lw=1.3, label="neutral median (band: 10th to 90th percentile)"),
               Line2D([0], [0], color="black", ls="--", lw=1.2, label="the record"),
               Line2D([0], [0], color="0.35", ls="-.", lw=1.0, label="neutral records outside their own band (B)"),
               Line2D([0], [0], color="0.0", ls="-", lw=0.9, marker="o", ms=2.6, label="planted, larger cluster (17 sites)"),
               Line2D([0], [0], color="0.0", ls=":", lw=0.9, marker="o", ms=2.6, label="planted, smaller cluster (11 sites)"),
               Line2D([0], [0], color="0.72", lw=2, label="multiplier 1.5"),
               Line2D([0], [0], color="0.45", lw=2, label="multiplier 2"),
               Line2D([0], [0], color="0.0", lw=2, label="multiplier 3"),
               Line2D([0], [0], marker="D", color="none", markerfacecolor="black", markeredgecolor="white", ms=6,
                      label="the record (C)"),
               Line2D([0], [0], marker="^", color="none", markerfacecolor="white", markeredgecolor="black", ms=6,
                      label="planted Walls Engraved, multiplier 3, larger / smaller cluster (C)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6.4, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(bottom=0.12, top=0.96, left=0.09, right=0.98)
    png = fs.save_all(fig, FIG, close=True)
    print(f"written {png}; record F {rec}; neutral medians at F_W 0.3 "
          f"{[cov[c]['cells'][4]['med'] for c in COMMON]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
