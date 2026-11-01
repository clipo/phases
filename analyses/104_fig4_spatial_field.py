#!/usr/bin/env python3
"""Figure 4: the spatial field (A0), drawn from committed posterior summaries.

Reads `output/findings/a0_intercept_prior.md`, the summary table of the three
A0 fits that differ only in the class-intercept prior (Normal(0, 2) inherited;
Normal(0, 10), the default; Normal centered on the pooled log-odds), and plots
the identified quantities: the standard deviation of the between-site contrast
in a class's logit share at the median nearest-neighbor spacing (4.3 km), the
median pair distance (36.6 km) and the largest pair distance (79.1 km), and,
inset, the contrast-based spatial share. The fit records (`output/gp_cache/`)
are gitignored and absent, so medians and 90 percent intervals are plotted,
not draws, and the caption says so. No model is refitted.

Usage:
    .venv/bin/python analyses/104_fig4_spatial_field.py
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

FIG = "fig4_spatial_field"
SRC = ROOT / "output" / "findings" / "a0_intercept_prior.md"
CELL = re.compile(r"([-0-9.]+) \[([-0-9.]+), ([-0-9.]+)\]")
PRIORS = [("P10", "Normal(0, 10), default", "o", "black"),
          ("P2", "Normal(0, 2)", "s", "0.45"),
          ("PC", "Normal(pooled log-odds, 2)", "^", "0.7")]


def read_table() -> tuple[list[str], dict[str, list[tuple[float, float, float]]]]:
    header, rows = None, {}
    for line in SRC.read_text().splitlines():
        if line.startswith("| quantity |"):
            header = [h.strip() for h in line.strip("|").split("|")][1:]
        elif header and line.startswith("|") and not line.startswith("|---"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            vals = [CELL.match(c) for c in cells[1:]]
            if all(vals):
                rows[cells[0]] = [tuple(float(g) for g in v.groups()) for v in vals]
    assert header and header[0].startswith("P2") and header[1].startswith("P10") and header[2].startswith("PC"), header
    return header, rows


def main() -> int:
    header, rows = read_table()
    col = {"P2": 0, "P10": 1, "PC": 2}
    dist_rows = {float(re.search(r"at ([0-9.]+) km", k).group(1)): v
                 for k, v in rows.items() if k.startswith("field contrast SD at")}
    share = next(v for k, v in rows.items() if k.startswith("contrast-based share"))
    assert sorted(dist_rows) == [4.3, 36.6, 79.1], sorted(dist_rows)
    assert [round(share[i][0], 2) for i in (col["P2"], col["P10"], col["PC"])] == [0.93, 0.93, 0.94], share

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    ds = np.array(sorted(dist_rows))
    for k, (key, label, mk, color) in enumerate(PRIORS):
        off = (k - 1) * 1.6
        med = np.array([dist_rows[d][col[key]][0] for d in ds])
        lo = np.array([dist_rows[d][col[key]][1] for d in ds])
        hi = np.array([dist_rows[d][col[key]][2] for d in ds])
        ax.errorbar(ds + off, med, yerr=[med - lo, hi - med], fmt=mk, color=color, ecolor=color,
                    elinewidth=1.0, capsize=2.5, markersize=5.5, label=label, zorder=5 - k)
    for d, txt in zip(ds, ["median nearest-\nneighbor spacing", "median pair\ndistance", "largest pair\ndistance"]):
        ax.text(d, 0.12, txt, ha="center", va="bottom", fontsize=6.5, color="0.35")
    ax.set_xlabel("distance between two sites (km)")
    ax.set_ylabel("SD of the between-site contrast\nin a class's logit share")
    ax.set_xlim(-4, 92)
    ax.set_ylim(0, 5.6)
    ax.legend(title="class-intercept prior", fontsize=7, title_fontsize=7, loc="upper left", frameon=False, bbox_to_anchor=(0.0, 0.98))
    fs.despine(ax)

    ins = ax.inset_axes([0.47, 0.68, 0.28, 0.26])
    for k, (key, label, mk, color) in enumerate(PRIORS):
        m, lo, hi = share[col[key]]
        ins.errorbar([m], [k], xerr=[[m - lo], [hi - m]], fmt=mk, color=color, ecolor=color,
                     elinewidth=1.0, capsize=2.5, markersize=5)
    ins.set_yticks(range(len(PRIORS)))
    ins.set_yticklabels([p[1].replace(", default", "") for p in PRIORS], fontsize=6)
    ins.set_xlim(0.7, 1.0)
    ins.set_xlabel("contrast-based spatial share", fontsize=6.5)
    ins.tick_params(axis="x", labelsize=6)
    ins.invert_yaxis()
    fs.despine(ins)

    fig.tight_layout()
    png = fs.save_all(fig, FIG, close=True)
    print(f"written {png}; shares " + ", ".join(f"{p[0]} {share[col[p[0]]][0]:.3f}" for p in PRIORS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
