"""House figure style for mls-emergence manuscript figures.

Import this module in any figure script to apply the global rcParams and
access the palette. Call ``save(fig, name)`` to write figures/{name} in all
formats (PNG, PDF, SVG, TIFF) at 300 dpi with a tight bounding box.

Color policy: the palette (``OKABE_ITO`` / the ``OI_*`` aliases) is grayscale
by default, so the print (main-text) figures render without color. Set the
environment variable ``MLS_FIG_COLOR=1`` before importing this module to switch
the palette to the true Okabe-Ito colors (used by the online-only supplement
figures). The always-color ``OKABE_ITO_COLOR`` / ``OIC_*`` aliases are available
for coloring a single figure in a script that otherwise makes grayscale figures.

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

import os
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
# Palette. Main-text figures print in grayscale (American Antiquity prints
# without color), so the default OKABE_ITO is a set of distinct grays and the
# scripts also vary marker and line style where more than three series share a
# panel. The supplement is online-only and may use color: set the environment
# variable MLS_FIG_COLOR=1 to make OKABE_ITO and the OI_* aliases resolve to the
# true Okabe-Ito colors instead. The always-color palette is also exposed as
# OKABE_ITO_COLOR / OIC_* for figures that opt into color directly (used to
# color a single figure in a script that also makes grayscale print figures).
# ---------------------------------------------------------------------------
_GRAY: list[str] = [
    "0.50",   # orange      -> medium gray (secondary series)
    "0.72",   # sky         -> light gray (tertiary series)
    "0.38",   # green       -> dark-medium gray
    "0.82",   # yellow      -> very light gray
    "0.20",   # blue        -> dark gray (primary series)
    "0.00",   # vermillion  -> black (highlight / observed)
    "0.62",   # purple      -> medium-light gray
    "0.00",   # black
]

# True Okabe-Ito colorblind-safe palette.
OKABE_ITO_COLOR: list[str] = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#000000",  # black
]

_USE_COLOR = os.environ.get("MLS_FIG_COLOR", "").strip().lower() not in ("", "0", "false", "no")
OKABE_ITO: list[str] = OKABE_ITO_COLOR if _USE_COLOR else _GRAY

# Named aliases (resolve to gray or color per MLS_FIG_COLOR)
OI_ORANGE   = OKABE_ITO[0]
OI_SKY      = OKABE_ITO[1]
OI_GREEN    = OKABE_ITO[2]
OI_YELLOW   = OKABE_ITO[3]
OI_BLUE     = OKABE_ITO[4]
OI_VERMIL   = OKABE_ITO[5]
OI_PURPLE   = OKABE_ITO[6]
OI_BLACK    = OKABE_ITO[7]

# Always-color aliases (independent of MLS_FIG_COLOR)
OIC_ORANGE  = OKABE_ITO_COLOR[0]
OIC_SKY     = OKABE_ITO_COLOR[1]
OIC_GREEN   = OKABE_ITO_COLOR[2]
OIC_YELLOW  = OKABE_ITO_COLOR[3]
OIC_BLUE    = OKABE_ITO_COLOR[4]
OIC_VERMIL  = OKABE_ITO_COLOR[5]
OIC_PURPLE  = OKABE_ITO_COLOR[6]
OIC_BLACK   = OKABE_ITO_COLOR[7]


# ---------------------------------------------------------------------------
# save helper
# ---------------------------------------------------------------------------
_FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


_FORMATS = (".png", ".pdf", ".svg", ".tif")


def save_all(fig: "plt.Figure", name_or_path, dpi: int = 300,
             close: bool = False) -> Path:
    """Write a figure in every distribution format: PNG and TIFF (raster at
    ``dpi``) and PDF and SVG (vector). ``name_or_path`` may be a bare figure
    name (written under figures/) or a path; any extension is replaced, so the
    four siblings share one stem. Returns the .png Path. Set close=True to close
    the figure afterwards."""
    p = Path(name_or_path)
    if p.suffix == "":
        _FIGURES_DIR.mkdir(exist_ok=True)
        stem = _FIGURES_DIR / p.name
    else:
        p.parent.mkdir(parents=True, exist_ok=True)
        stem = p.with_suffix("")
    for ext in _FORMATS:
        fig.savefig(f"{stem}{ext}", dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return Path(f"{stem}.png")


def save(fig: "plt.Figure", name: str) -> Path:
    """Write figures/{name} in all formats (PNG, PDF, SVG, TIFF) at 300 dpi,
    tight bounding box, and close the figure. Returns the .png Path."""
    return save_all(fig, name, dpi=300, close=True)


def despine(ax: "mpl.axes.Axes") -> None:
    """Remove top and right spines (applied globally via rcParams, but
    useful for axes added after the fact)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
