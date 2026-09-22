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
    # Which distance does the pottery follow? Added 2026-09-21. The supplement
    # justified the river metric as the more faithful representation of how
    # potters could interact; that is a claim about the record and this measures
    # it. Chi-square distance between decorated-class profiles, rank-correlated
    # with each metric, and each metric with the other held fixed (partial rank
    # correlation by residualising ranks). Descriptive (rule 18).
    from scipy.stats import rankdata, spearmanr
    m = counts.to_numpy(float)
    prof = m / m.sum(1, keepdims=True)
    cm = prof.mean(0)
    X = prof / np.sqrt(np.where(cm > 0, cm, 1.0))
    cer = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(2))[iu][ok]
    rr, gg = r[ok], g[ok]

    def partial(y, x, z):
        ry, rx, rz = rankdata(y), rankdata(x), rankdata(z)
        ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
        ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
        return float(np.corrcoef(ey, ex)[0, 1])

    nn = np.where(np.eye(n, dtype=bool), np.inf, geo).min(1)
    L += [
        "",
        "## Which distance does the pottery follow?",
        "",
        f"Chi-square distance between decorated-class profiles, {int(ok.sum())} pairs.",
        "",
        "| | rank correlation with ceramic distance |",
        "|---|---|",
        f"| straight-line distance | {spearmanr(cer, gg).statistic:.3f} |",
        f"| river-network distance | {spearmanr(cer, rr).statistic:.3f} |",
        f"| river, straight-line held fixed | {partial(cer, rr, gg):.3f} |",
        f"| straight-line, river held fixed | {partial(cer, gg, rr):.3f} |",
        "",
        f"The two metrics correlate at {spearmanr(rr, gg).statistic:.3f}. If the pottery "
        "followed the waterways, the river",
        "metric would carry information the straight-line one lacks. Read the third row "
        "for that. A null",
        "there does not show that waterways were unimportant: the network is MODERN "
        "hydrography, and the",
        "St. Francis, Tyronza and Mississippi have all moved since these sites were "
        "occupied, so the metric",
        "may be measuring the wrong rivers.",
    ]
    # The comparison that matters is BETWEEN phases (author, 2026-09-21). Within
    # a phase sites are close and the two metrics barely differ; whether the
    # waterways structured interaction is a question about pairs in different
    # phases, and about the phases taken as units.
    import itertools
    ph = importlib.import_module("36_canonical_phase_map")
    lab = np.array(ph.assign_phases_by_territory([str(i) for i in counts.index], coords)[0])
    same = (lab[:, None] == lab[None, :])[iu][ok]
    L += ["", "### Within phases and between them", "",
          "| pairs | n | straight-line | river | river, straight-line held fixed | "
          "straight-line, river held fixed |", "|---|---|---|---|---|---|"]
    for label, sel in (("within a phase", same), ("between phases", ~same)):
        L.append(f"| {label} | {int(sel.sum())} | {spearmanr(cer[sel], gg[sel]).statistic:.3f} | "
                 f"{spearmanr(cer[sel], rr[sel]).statistic:.3f} | "
                 f"{partial(cer[sel], rr[sel], gg[sel]):.3f} | "
                 f"{partial(cer[sel], gg[sel], rr[sel]):.3f} |")
    phases = sorted(set(lab))
    rows = []
    for a, b in itertools.combinations(phases, 2):
        ia, ib = np.where(lab == a)[0], np.where(lab == b)[0]
        pa = m[ia].sum(0) / m[ia].sum()
        pb = m[ib].sum(0) / m[ib].sum()
        w = np.sqrt(np.where(cm > 0, cm, 1.0))
        rows.append((f"{a} - {b}", float(np.sqrt((((pa - pb) / w) ** 2).sum())),
                     float(geo[np.ix_(ia, ib)].mean()), float(R[np.ix_(ia, ib)].mean())))
    rows.sort(key=lambda t: t[1])
    L += ["", "### The phases as units", "",
          "Pooled class profile per phase; distances are means over member pairs.", "",
          "| phase pair | ceramic distance | straight-line km | river km | detour |",
          "|---|---|---|---|---|"]
    for nm, cd_, g_, r_ in rows:
        L.append(f"| {nm} | {cd_:.3f} | {g_:.1f} | {r_:.1f} | {r_ / g_:.2f} |")
    cds = [t[1] for t in rows]
    L += ["", f"Across the {len(rows)} phase pairs ceramic distance rank-correlates "
              f"{spearmanr(cds, [t[2] for t in rows]).statistic:.3f} with straight-line "
              f"distance and {spearmanr(cds, [t[3] for t in rows]).statistic:.3f} with river "
              "distance. Ten pairs is indicative, no more.", "",
          "The largest detours join phases on the St. Francis to phases on the Mississippi, "
          "which the network can",
          "connect only through their confluence far to the south. Those distances describe "
          "the routing of the modern",
          "network, and they are what the drift model's copying kernel uses.",
    ]
    L += [
        "",
        "## The grain of the record",
        "",
        f"Nearest-neighbour spacing of the {n} assemblages: median {np.median(nn):.1f} km, "
        f"quartiles {np.percentile(nn, 25):.1f} to {np.percentile(nn, 75):.1f} km. No "
        "spatial statistic here",
        "can resolve structure much finer than that.",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[6:19]))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
