#!/usr/bin/env python3
"""Figure 3: what recovers the phases.

(A) The 28 analyzed assemblages on the map with Phillips's (1970) phase areas,
marked by the phase each takes, with the assemblages that Ward linkage on site
COORDINATES ALONE places with another phase ringed and named. The clustering is
recomputed here (it is deterministic and takes milliseconds) and the result is
checked against the committed findings in `output/findings/phase_recovery.md`
and `phase_recovery.csv` (rule 19): the misplaced set must be Big Eddy and
Castile Landing and the adjusted Rand index over the 26 assemblages inside a
phase area must read 0.740.

(B) Late-period settlement density along the reconstructed phase boundaries
against the boundaries of alternative divisions, read from the committed table
in `output/findings/boundaries_in_settlement_gaps.md` (produced by
`80_boundaries_in_settlement_gaps.py`). Nothing in panel B is recomputed; the
per-division densities behind the medians and shares are not committed, so the
panel shows the table's medians and shares, which is what the text quotes.

Usage:
    .venv/bin/python analyses/103_fig3_phase_recovery_map.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))
fs = importlib.import_module("figstyle")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

FIG = "fig3_phase_recovery_map"
SRC_MD = ROOT / "output" / "findings" / "phase_recovery.md"
SRC_CSV = ROOT / "output" / "findings" / "phase_recovery.csv"
SRC_GAPS = ROOT / "output" / "findings" / "boundaries_in_settlement_gaps.md"
EXPECTED_MISPLACED = {"Big_Eddy", "Castile_Landing"}
OUTSIDE_AREA = {"Cramor_Place", "Holden_Lake"}  # take the nearest area; not scored
PHASES = ("Parkin", "Kent", "Walls")
MARKERS = {"Parkin": "o", "Kent": "s", "Walls": "^"}


def committed_values() -> dict:
    """ARI values the text quotes, from the committed findings."""
    import csv
    rows = list(csv.DictReader(SRC_CSV.open()))
    pick = lambda tr, w, alg: float(next(r for r in rows if r["transform"] == tr  # noqa: E731
                                         and float(r["w"]) == w and r["algorithm"] == alg)["ari_mapped"])
    md = SRC_MD.read_text()
    m = re.search(r"\| chisq \| pottery alone \| [0-9.]+ \| \*\*([0-9.]+)\*\* \|", md)
    rnd = re.search(r"boundaries placed at random agree with the phases at a median ARI of ([0-9.]+)", md)
    assert m and rnd, "phase_recovery.md no longer carries the rows this figure reads"
    return {"ward_map": pick("raw", 0.0, "ward"), "kmeans_map": pick("raw", 0.0, "kmeans"),
            "composition": float(m.group(1)), "random": float(rnd.group(1))}


def gaps_table() -> dict[float, dict[str, float]]:
    out = {}
    for line in SRC_GAPS.read_text().splitlines():
        m = re.match(r"\| (\d+) km \| ([0-9.]+) \| \*\*([0-9.]+)\*\* \| ([0-9.]+) \| (\d+)% \| ([0-9.]+) \| (\d+)% \|", line)
        if m:
            bw, at, ph, rmed, rsh, cmed, csh = m.groups()
            out[float(bw)] = {"assemblages": float(at), "phase": float(ph), "random_med": float(rmed),
                              "random_share": int(rsh), "compact_med": float(cmed), "compact_share": int(csh)}
    assert set(out) == {3.0, 5.0, 8.0}, out
    return out


def main() -> int:
    mf = importlib.import_module("make_figures")
    m36 = importlib.import_module("36_canonical_phase_map")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    latlon = coords.to_numpy(float)
    labels, _ = m36.assign_primary_phases(names, latlon)
    labels = np.asarray(labels, dtype=object)
    E, N = m36.to_utm(latlon[:, 1], latlon[:, 0])
    E, N = np.asarray(E, float), np.asarray(N, float)
    assert len(names) == 28, len(names)

    # Ward linkage on coordinates alone, three clusters, majority phase per cluster.
    Z = linkage(np.column_stack([E, N]) / 1000.0, method="ward")
    cl = fcluster(Z, t=3, criterion="maxclust")
    inside = np.array([n not in OUTSIDE_AREA for n in names])
    majority = {}
    for c in np.unique(cl):
        members = labels[(cl == c) & inside]
        vals, cnt = np.unique(members, return_counts=True)
        majority[c] = vals[np.argmax(cnt)]
    misplaced = {n for n, c, ph, ins in zip(names, cl, labels, inside) if ins and majority[c] != ph}
    assert misplaced == EXPECTED_MISPLACED, f"Ward misplaces {sorted(misplaced)}; committed: {sorted(EXPECTED_MISPLACED)}"
    v = committed_values()
    assert abs(v["ward_map"] - 0.740) < 5e-4 and abs(v["kmeans_map"] - 0.480) < 5e-4, v
    assert abs(v["composition"] - 0.358) < 5e-4 and abs(v["random"] - 0.485) < 5e-4, v
    gaps = gaps_table()
    g5 = gaps[5.0]
    assert (g5["phase"], g5["assemblages"], g5["random_share"], g5["compact_share"]) == (0.56, 1.01, 94, 100), g5

    fig = plt.figure(figsize=(7.0, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    pad = 14_000.0
    ext = (E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad)
    m36.draw_base(ax, ext, {})
    outline = m36.load_phillips_outlines()
    for ph in PHASES + ("Nodena",):
        if ph in outline:
            m36.fill_polygon(ax, outline[ph], facecolor="0.93", edgecolor="none", zorder=1.5)
            xs, ys = outline[ph].exterior.xy
            ax.plot(xs, ys, color="0.25", lw=0.9, zorder=6)
    fills = {"Parkin": fs.OI_BLUE, "Kent": fs.OI_ORANGE, "Walls": fs.OI_GREEN}
    for ph in PHASES:
        m = labels == ph
        ax.scatter(E[m], N[m], s=30, marker=MARKERS[ph], facecolor=fills[ph], edgecolor="black",
                   lw=0.5, zorder=10, label=f"{ph} phase")
    for n, e, nn in zip(names, E, N):
        if n in misplaced:
            ax.scatter([e], [nn], s=150, facecolor="none", edgecolor=fs.OI_VERMIL, lw=1.4, zorder=12)
            # Label to the LEFT of the point: the St. Francis River label runs along
            # the right of both misplaced sites and a right-hand label sat on it.
            ax.annotate(n.replace("_", " "), (e, nn), xytext=(-9, -3), textcoords="offset points",
                        ha="right", fontsize=6.5, color=fs.OI_VERMIL, zorder=13)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    handles = [Patch(facecolor="0.93", edgecolor="0.25", lw=0.9, label="Phase area (Phillips 1970)")]
    handles += [Line2D([0], [0], marker=MARKERS[ph], color="none", markerfacecolor=fills[ph],
                       markeredgecolor="black", markeredgewidth=0.5, markersize=6, label=f"{ph} phase")
                for ph in PHASES]
    handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
                          markeredgecolor=fs.OI_VERMIL, markeredgewidth=1.4, markersize=11,
                          label="placed with another phase by\nWard linkage on coordinates alone"))
    # Below the map, not inside it: inside at lower left it covered three Kent-phase sites.
    ax.legend(handles=handles, fontsize=6.2, loc="upper left", bbox_to_anchor=(0.0, -0.01), ncol=2,
              frameon=False, borderpad=0.4, handletextpad=0.5, labelspacing=0.45, columnspacing=1.2)
    fs.panel_label(ax, "A")

    # Panel B: densities from the committed table.
    bws = [3.0, 5.0, 8.0]
    xs = np.arange(len(bws))
    series = [("assemblages", "o", "black", "black", "at the analyzed assemblages"),
              ("phase", "*", fs.OI_VERMIL, fs.OI_VERMIL, "along the phase boundaries"),
              ("random_med", "o", "white", "0.3", "random-center divisions, median"),
              ("compact_med", "s", "white", "0.3", "compact divisions, median")]
    for key, mk, fc, ec, lab in series:
        ys = [gaps[b][key] for b in bws]
        axb.plot(xs, ys, color="0.6", lw=0.6, zorder=1)
        axb.scatter(xs, ys, marker=mk, s=70 if mk == "*" else 36, facecolor=fc, edgecolor=ec, lw=0.9,
                    zorder=5, label=lab)
    axb.annotate(f"{g5['random_share']}% of random divisions\nand {g5['compact_share']}% of compact divisions\n"
                 "cross denser settlement", xy=(1, g5["phase"]), xytext=(1.08, 0.18), fontsize=6.5,
                 arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6), color="0.2")
    axb.set_xticks(xs)
    axb.set_xticklabels([f"{int(b)} km" for b in bws])
    axb.set_xlabel("kernel bandwidth")
    axb.set_ylabel("late-period sites per 100 km$^2$")
    axb.set_ylim(0, 1.15)
    axb.set_xlim(-0.4, 2.6)
    # Above the axes: inside at upper right it sat on the assemblage series near 1.0.
    axb.legend(fontsize=6.2, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, frameon=False,
               labelspacing=0.5, columnspacing=1.0, handletextpad=0.4)
    fs.despine(axb)
    fs.panel_label(axb, "B", y=1.16)

    fig.subplots_adjust(left=0.03, right=0.98, top=0.90, bottom=0.17)
    png = fs.save_all(fig, FIG, close=True)
    print(f"written {png}; Ward misplaces {sorted(misplaced)}; ARI map/Ward {v['ward_map']:.3f}, "
          f"map/k-means {v['kmeans_map']:.3f}, composition {v['composition']:.3f}, random {v['random']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
