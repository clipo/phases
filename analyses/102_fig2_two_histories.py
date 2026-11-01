#!/usr/bin/env python3
"""Figure 2: two histories, one record. An illustrative schematic, not data.

The St. Francis basin assemblages sit in clusters separated by stretches with
few or no collected assemblages. This figure shows why that geometry limits
what the record can say (register D24; recovery tests of A1 and A3 in
`a1_validation.md` and `a3_validation.md`): the SAME sixteen sampled values
are consistent with a gradient in a class's share across the basin (A) and
with a step at the gap between the clusters (B). The two histories differ only
where no assemblage was collected. Everything drawn here is synthetic, from a
fixed seed; the caption labels it illustrative and scopes the point to this
sample.

Usage:
    .venv/bin/python analyses/102_fig2_two_histories.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))
fs = importlib.import_module("figstyle")
import matplotlib.pyplot as plt  # noqa: E402  (after figstyle sets rcParams)

FIG = "fig2_two_histories"
SEED = 20261031
N_PER_CLUSTER = 8
GAP = (42.0, 58.0)  # km along the transect with no assemblages


def main() -> int:
    rng = np.random.default_rng(SEED)
    x_a = np.sort(rng.uniform(6.0, 40.0, N_PER_CLUSTER))
    x_b = np.sort(rng.uniform(60.0, 94.0, N_PER_CLUSTER))
    x = np.concatenate([x_a, x_b])
    slope = 0.028  # logits per km, illustrative
    gradient = lambda t: -1.4 + slope * t  # noqa: E731
    y = gradient(x) + rng.normal(0.0, 0.22, x.size)
    m_a, m_b = y[:N_PER_CLUSTER].mean(), y[N_PER_CLUSTER:].mean()

    t = np.linspace(0, 100, 400)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
    for ax, letter in zip(axes, "AB"):
        ax.axvspan(*GAP, color="0.92", zorder=0, lw=0)
        ax.text(np.mean(GAP), 2.05, "no assemblages", ha="center", va="top",
                fontsize=7, color="0.35")
        if letter == "A":
            ax.plot(t, gradient(t), color="0.2", lw=1.4, zorder=2,
                    label="gradient: copying falls with distance")
        else:
            ax.plot(t[t < 50], np.full((t < 50).sum(), m_a), color="0.2", lw=1.4, zorder=2)
            ax.plot(t[t >= 50], np.full((t >= 50).sum(), m_b), color="0.2", lw=1.4, zorder=2,
                    label="step at the gap: a bounded unit")
            ax.plot([50, 50], [m_a, m_b], color="0.2", lw=1.4, ls=":", zorder=2)
        ax.scatter(x, y, s=22, facecolor="white", edgecolor="black", lw=0.8, zorder=5,
                   label="sampled assemblages (the same 16 values in A and B)")
        ax.set_xlim(0, 100)
        ax.set_ylim(-2.1, 2.1)
        ax.set_xlabel("position along the basin (km)")
        fs.panel_label(ax, letter)
        fs.despine(ax)
    axes[0].set_ylabel("logit share of a decorated class")
    h0, l0 = axes[0].get_legend_handles_labels()
    h1, l1 = axes[1].get_legend_handles_labels()
    fig.legend(h0 + h1[:1], l0 + l1[:1], loc="lower center", ncol=3, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.08))
    fig.tight_layout()
    png = fs.save_all(fig, FIG, close=True)
    print(f"written {png} and siblings; illustrative, seed {SEED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
