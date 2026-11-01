#!/usr/bin/env python3
"""Figure 7: recovery of planted structure on this site geometry (A1 and A3).

Drawn from the committed recovery tables; no model is refitted.

(A) Latent-unit model A1 (`output/findings/a1_validation.md`, produced by
`analyses/probe_latent_groups/validate_a1.R`): posterior median and 90 percent
interval of the between-unit scale sigma_g against the planted value on records
simulated from the model's own generator at the 28 sites, plus the fit to a
record with no units. Fits that fail the sampler diagnostics are drawn open; they are
failed validation attempts, not estimates, and are shown so the denominator is
visible (12 of 14 recovery fits fail; the no-unit fit fails).

(B) Barrier model A3 (`output/findings/a3_validation.md`, from `validate_a3.R`;
free range, class-intercept prior Normal(0, 10)): posterior median and 90
percent interval of the barrier parameter omega against the planted value, for
the two edge placements, with the prior median marked. The six prior-by-range
combinations are tabulated in `a3_prior_comparison.md` and the supplement; this
panel shows the one named in the caption.

Usage:
    .venv/bin/python analyses/107_fig7_latent_unit_recovery.py
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

FIG = "fig7_latent_unit_recovery"
SRC_A1 = ROOT / "output" / "findings" / "a1_validation.md"
SRC_A3 = ROOT / "output" / "findings" / "a3_validation.md"
POST = re.compile(r"([-0-9.]+) \[([-0-9.]+), ([-0-9.]+)\]")
PRIOR_MEDIAN = 1.35  # half-Normal(0, 2) on sigma_g and on omega; stated in both findings files


def a1_rows() -> tuple[list[dict], dict]:
    text = SRC_A1.read_text()
    rec = text.split("## (iii) Recovery against the prior", 1)[1]
    rows = []
    for line in rec.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 10 and cells[0].startswith("sg"):
            m = POST.match(cells[4])
            rows.append(dict(record=cells[0], true=float(cells[1]), med=float(m.group(1)), lo=float(m.group(2)),
                             hi=float(m.group(3)), passed=cells[9].strip() == "pass"))
    assert len(rows) == 14, len(rows)
    null_block = text.split("Truth: sigma_g = 0, no groups", 1)[1]
    m = re.search(r"\| sigma_g \| ([-0-9.]+) \[([-0-9.]+), ([-0-9.]+)\] \|", null_block)
    gate = re.search(r"\| A1 on zero \| \d+ \| \d+ \| (pass|fail) \|", text)
    null = dict(true=0.0, med=float(m.group(1)), lo=float(m.group(2)), hi=float(m.group(3)), passed=gate.group(1) == "pass")
    return rows, null


def a3_rows() -> list[dict]:
    rows = []
    for line in SRC_A3.read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 11 and cells[0].startswith("own_"):
            tr = re.search(r"omega ([0-9.]+)", cells[3])
            m = POST.match(cells[5])
            rows.append(dict(record=cells[0], placement=cells[2], true=float(tr.group(1)), med=float(m.group(1)),
                             lo=float(m.group(2)), hi=float(m.group(3)), passed=cells[10] == "pass"))
    assert len(rows) == 24, len(rows)
    return rows


def main() -> int:
    r1, null = a1_rows()
    r3 = a3_rows()
    half = [r["med"] for r in r1 if r["true"] == 0.5]
    assert (round(min(half), 2), round(max(half), 2)) == (0.27, 0.48), half
    assert round(null["med"], 2) == 0.52 and not null["passed"], null
    assert sum(not r["passed"] for r in r1) == 12
    om0 = [r["med"] for r in r3 if r["true"] == 0.0]
    om1 = [r["med"] for r in r3 if r["true"] == 1.0]
    assert (round(min(om0), 2), round(max(om0), 2)) == (0.48, 1.21), om0
    assert (round(min(om1), 2), round(max(om1), 2)) == (0.70, 2.28), om1
    assert all(r["passed"] for r in r3)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.5))
    # A1
    for i, r in enumerate([null] + r1):
        x = r["true"] + (i % 5 - 2) * 0.03
        fc = "black" if r["passed"] else "white"
        axA.errorbar([x], [r["med"]], yerr=[[r["med"] - r["lo"]], [r["hi"] - r["med"]]], fmt="o", color="black",
                     ecolor="0.4", markerfacecolor=fc, elinewidth=0.8, capsize=1.5, markersize=4.5, zorder=5)
    lim = (-0.1, 2.2)
    axA.plot(lim, lim, color="0.6", lw=0.8, ls="-", zorder=1)
    axA.axhline(PRIOR_MEDIAN, color="0.3", lw=0.9, ls="--", zorder=1)
    axA.set_xlabel("planted between-unit scale $\\sigma_g$ (logits)")
    axA.set_ylabel("posterior median, 90% interval")
    axA.set_xlim(*lim)
    axA.set_ylim(0, 4.0)
    axA.set_xticks([0, 0.25, 0.5, 1.0, 2.0])
    axA.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor="black", markeredgecolor="black",
                               ms=5, label="passes the sampler diagnostics"),
                        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black",
                               ms=5, label="fails the sampler diagnostics"),
                        Line2D([0], [0], color="0.3", ls="--", lw=0.9, label="prior median"),
                        Line2D([0], [0], color="0.6", lw=0.8, label="recovery (planted = posterior)")],
               fontsize=6.2, loc="upper left", frameon=False)
    fs.despine(axA)
    fs.panel_label(axA, "A")
    # A3
    for r in r3:
        rep = int(r["record"].rsplit("_r", 1)[1])
        x = r["true"] + (rep - 2) * 0.03 + (-0.045 if r["placement"] == "G1" else 0.045)
        mk = "o" if r["placement"] == "G1" else "s"
        axB.errorbar([x], [r["med"]], yerr=[[r["med"] - r["lo"]], [r["hi"] - r["med"]]], fmt=mk, color="black",
                     ecolor="0.4", elinewidth=0.8, capsize=1.5, markersize=4.5, zorder=5)
    lim = (-0.1, 1.15)
    axB.plot(lim, lim, color="0.6", lw=0.8, zorder=1)
    axB.axhline(PRIOR_MEDIAN, color="0.3", lw=0.9, ls="--", zorder=1)
    axB.set_xlabel("planted barrier $\\omega$")
    axB.text(0.60, 0.98, "free range\nintercept prior Normal(0, 10)", transform=axB.transAxes, ha="center", va="top",
             fontsize=6.5, color="0.3")  # centered over omega 0.6, clear of the bars at 0.5 and 1.0
    axB.set_xlim(*lim)
    axB.set_ylim(0, 4.0)
    axB.set_xticks([0, 0.25, 0.5, 1.0])
    axB.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor="black", ms=5,
                               label="edge placement G1"),
                        Line2D([0], [0], marker="s", color="none", markerfacecolor="black", ms=5,
                               label="edge placement L"),
                        Line2D([0], [0], color="0.3", ls="--", lw=0.9, label="prior median")],
               fontsize=6.2, loc="upper left", frameon=False)
    fs.despine(axB)
    fs.panel_label(axB, "B")
    fig.tight_layout()
    png = fs.save_all(fig, FIG, close=True)
    print(f"written {png}; A1 planted 0.5 medians {min(half):.3f} to {max(half):.3f}, null {null['med']:.3f}; "
          f"A3 omega 0 {min(om0):.2f} to {max(om0):.2f}, omega 1 {min(om1):.2f} to {max(om1):.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
