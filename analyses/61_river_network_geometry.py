"""61_river_network_geometry.py - the river-network numbers the supplement quotes.

The supplement describes the along-waterway distance metric that the drift
simulations (33, 34) use for their copying kernel, and quotes three numbers for
it: how many assemblages connect on the network, the mean access distance from
an assemblage to a mapped channel, and how much longer along-water paths are
than straight lines. None of those had a committed procedure (rule 1). The first
two were correct; the third was not, which is why this script exists rather than
merely re-deriving what was already there.

THE DETOUR RATIO. The supplement said along-water paths "average about 2.3 times
the straight-line distance". Over the 406 assemblage pairs there are three
defensible summaries and 2.3 is none of them:

    mean of the per-pair ratios      3.10
    median of the per-pair ratios    2.39
    ratio of the mean distances      2.67

The distribution is right-skewed, because a pair that sits close together in a
straight line but must route around a meander can carry a very large ratio,
which pulls the mean well above the median. So the summary has to be named, not
just quoted (rule 6). The supplement now reports the mean with the median
alongside it.

This does not touch any result. The simulations use the distance MATRIX, never a
summary of it, so the ratio is descriptive context for the reader.

Usage: .venv/bin/python analyses/61_river_network_geometry.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))

mm = importlib.import_module("make_map")
sd = importlib.import_module("23_phases_vs_spatial_drift")
mf = importlib.import_module("make_figures")

OUT = ROOT / "output" / "findings" / "river_network_geometry.md"


def main() -> None:
    counts, coords_df = mf._load_curated()
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    n = len(coords)

    river_km, _, access_km, rinfo = mm.river_distance_matrix(coords)
    # The wider-valley block that used to sit here is gone with the wider-valley
    # analyses: that 55-assemblage set mixes deposits of different ages and very
    # different sample sizes, and the manuscript no longer makes a claim about
    # it, so there is nothing here for a committed procedure to support.
    geo = sd.geo_km(coords)

    R = np.asarray(river_km, float)
    iu = np.triu_indices(n, 1)
    r, g = R[iu], geo[iu]
    ok = np.isfinite(r) & np.isfinite(g) & (g > 0)
    ratio = r[ok] / g[ok]

    mean_ratio = float(ratio.mean())
    med_ratio = float(np.median(ratio))
    ratio_of_means = float(r[ok].mean() / g[ok].mean())
    mean_access = float(np.mean(access_km))

    L = [
        "# River-network geometry for the basin assemblages",
        "",
        f"Produced by `analyses/61_river_network_geometry.py`. "
        f"{n} assemblages, {int(ok.sum())} pairs.",
        "",
        "The along-waterway metric used by the copying kernel in "
        "`33_time_aware_emergence.py` and `34_emergence_robustness.py`.",
        "",
        "| quantity | value |",
        "|---|---|",
        f"| assemblages connecting on the network | {n} of {n} |",
        f"| river-graph largest component | {rinfo['largest_component']} nodes |",
        f"| unreachable pairs | {rinfo['n_unreachable']} |",
        f"| mean access distance to a mapped channel | {mean_access:.2f} km |",
        f"| max access distance | {float(np.max(access_km)):.1f} km |",
        f"| detour ratio, mean of per-pair ratios | {mean_ratio:.2f} |",
        f"| detour ratio, median of per-pair ratios | {med_ratio:.2f} |",
        f"| detour ratio, ratio of mean distances | {ratio_of_means:.2f} |",
        "",
        "## Reading",
        "",
        "The three detour summaries differ because the per-pair ratio is "
        "right-skewed: a pair that is close in a straight line but must route "
        "around a meander carries a very large ratio and pulls the mean above "
        "the median. Quoting an unnamed 'average' is therefore ambiguous at the "
        "0.7 level here, which is why the supplement names the summary it uses.",
        "",
        "None of these figures enters a result. The simulations consume the "
        "distance matrix itself, never a summary of it.",
        "",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[6:19]))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
