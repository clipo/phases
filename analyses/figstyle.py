"""House figure style for mls-emergence manuscript figures.

Import this module in any figure script to apply the global rcParams and
access the Okabe-Ito palette. Call ``save(fig, name)`` to write
figures/{name}.png at 300 dpi with tight bounding box.

Example
-------
from analyses.figstyle import OKABE_ITO, save
import matplotlib.pyplot as plt
...
fig, ax = plt.subplots()
...
save(fig, "fig2_validation")
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Global rcParams (house style)
# ---------------------------------------------------------------------------
_SANS = ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"]

mpl.rcParams.update(
    {
        # Font
        "font.family": "sans-serif",
        "font.sans-serif": _SANS,
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        # Resolution
        "figure.dpi": 150,       # screen preview; save() overrides to 300
        "savefig.dpi": 300,
        # Spine cleanup
        "axes.spines.top": False,
        "axes.spines.right": False,
        # Figure size defaults (7 in wide = single double-column; 3.5 in = single)
        "figure.figsize": (7, 4.5),
        # Grid
        "axes.grid": False,
        # Lines
        "lines.linewidth": 1.5,
        "patch.linewidth": 0.8,
    }
)

# ---------------------------------------------------------------------------
# Grayscale palette (American Antiquity prints without color). The names keep
# the former Okabe-Ito aliases so figure scripts need no change; each maps to a
# distinct gray chosen for contrast in the common pairings (primary series dark,
# highlight/observed black). Where more than three series share a panel, the
# scripts also vary marker and line style so the figures read in grayscale.
# ---------------------------------------------------------------------------
OKABE_ITO: list[str] = [
    "0.50",   # was orange  -> medium gray (secondary series)
    "0.72",   # was sky     -> light gray (tertiary series)
    "0.38",   # was green   -> dark-medium gray
    "0.82",   # was yellow  -> very light gray
    "0.20",   # was blue    -> dark gray (primary series)
    "0.00",   # was vermillion -> black (highlight / observed)
    "0.62",   # was purple  -> medium-light gray
    "0.00",   # black
]

# Named aliases for convenience
OI_ORANGE   = OKABE_ITO[0]
OI_SKY      = OKABE_ITO[1]
OI_GREEN    = OKABE_ITO[2]
OI_YELLOW   = OKABE_ITO[3]
OI_BLUE     = OKABE_ITO[4]
OI_VERMIL   = OKABE_ITO[5]
OI_PURPLE   = OKABE_ITO[6]
OI_BLACK    = OKABE_ITO[7]


# ---------------------------------------------------------------------------
# save helper
# ---------------------------------------------------------------------------
_FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


def save(fig: "plt.Figure", name: str) -> Path:
    """Write figures/{name}.png at 300 dpi, tight bounding box.

    The figures/ directory is created if absent. Returns the Path written.
    """
    _FIGURES_DIR.mkdir(exist_ok=True)
    out = _FIGURES_DIR / f"{name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def despine(ax: "mpl.axes.Axes") -> None:
    """Remove top and right spines (applied globally via rcParams, but
    useful for axes added after the fact)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
