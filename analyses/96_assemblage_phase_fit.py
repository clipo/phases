#!/usr/bin/env python3
"""Which phase does each assemblage's pottery fit best?

If the phases were bounded groups, an assemblage's decorated pottery should
look more like the rest of its own phase than like another phase. This asks
that of every assemblage directly, one at a time.

THE MEASUREMENT. For each assemblage, the chi-square distance from its
decorated-class profile to the pooled profile of each phase, with the
assemblage itself left out of its own phase's pool (so it is not compared
with itself). The distance is the one analyses 76 and 83 use: class
differences divided by the square root of the mean class proportion across
assemblages. Smaller is closer, and an assemblage "fits" the phase it is
closest to.

Sampling is carried through (rule 18): the assemblage's proportions and each
phase's pooled proportions are drawn from Dirichlet(counts + 1/2), DRAWS
times, and the share of draws in which the assemblage is closer to its own
phase than to every other is the posterior probability that its own phase
fits best. An assemblage with a probability near one half fits two phases
about equally; that is a statement about the assemblage, not a test.

The straight-line distance to the nearest member of each phase is printed
beside it, so that "fits another phase ceramically" can be read against
"lies nearer another phase on the ground".

A SENSITIVITY CASE. Two assemblages in the narrow southern arm of Phillips's
Parkin area, Big Eddy and Castile Landing, are the only ones that Ward or
average-linkage clustering on the site map places with another phase
(analysis 76; k-means misplaces four more). The last section
reruns the recovery and the phase-line comparison with those two counted as
Kent. It is a sensitivity case, not a relabeling: the primary analysis keeps
Phillips's assignments, because redrawing a scheme with the data that test it
is the circularity the paper is about.

Output: output/findings/assemblage_phase_fit.md
Usage: PYTHONPATH=src .venv/bin/python analyses/96_assemblage_phase_fit.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "assemblage_phase_fit.md"
DRAWS = 4000
SEED = 96
EVEN_BAND = (0.25, 0.75)   # a probability inside this band is called "about even"
MOVED = ("Big_Eddy", "Castile_Landing")


def chisq(p, q, w):
    return float(np.sqrt((((p - q) / w) ** 2).sum()))


def main() -> int:
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    t74 = importlib.import_module("74_phase_partition_test")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m = counts.to_numpy(float)
    xy = coords.to_numpy(float)
    labels, derived = ph.assign_primary_phases(names, xy)
    labels = np.asarray(labels)
    phases = sorted(set(labels))
    pts = t74.km_xy(xy)
    p = m / m.sum(1, keepdims=True)
    cm = p.mean(0)
    w = np.sqrt(np.where(cm > 0, cm, 1.0))
    rng = np.random.default_rng(SEED)

    rows = []
    for i, name in enumerate(names):
        pools = {q: m[[j for j in range(len(names)) if labels[j] == q and j != i]].sum(0)
                 for q in phases}
        dist = {q: chisq(p[i], pools[q] / pools[q].sum(), w) for q in phases}
        near = {q: min(float(np.linalg.norm(pts[i] - pts[j]))
                       for j in range(len(names)) if labels[j] == q and j != i) for q in phases}
        best = min(dist, key=dist.get)
        own = labels[i]
        wins = 0
        for _ in range(DRAWS):
            pi = rng.dirichlet(m[i] + 0.5)
            d = {q: chisq(pi, rng.dirichlet(pools[q] + 0.5), w) for q in phases}
            wins += min(d, key=d.get) == own
        rows.append(dict(name=name, own=own, derived=bool(derived[i]), n=int(m[i].sum()),
                         dist=dist, near=near, best=best, p_own=wins / DRAWS,
                         nearest_phase=min(near, key=near.get)))

    other = [r for r in rows if r["best"] != r["own"]]
    even = [r for r in rows if EVEN_BAND[0] <= r["p_own"] <= EVEN_BAND[1]]
    clear_own = [r for r in rows if r["p_own"] > EVEN_BAND[1]]
    clear_other = [r for r in rows if r["p_own"] < EVEN_BAND[0]]
    geo_other = [r for r in rows if r["nearest_phase"] != r["own"]]
    L = ["# Which phase does each assemblage's pottery fit best?", "",
         f"Produced by `analyses/96_assemblage_phase_fit.py`. {len(names)} assemblages, phases "
         + ", ".join(f"{q} {int((labels == q).sum())}" for q in phases)
         + " (Phillips's areas). Ceramic distance is the chi-square distance from the assemblage's "
         "decorated-class profile to each phase's pooled profile, the assemblage left out of its own "
         f"phase's pool. P(own) is the share of {DRAWS:,} draws of Dirichlet(counts + 1/2), for the "
         "assemblage and for every pool, in which its own phase is the closest. km is the "
         "straight-line distance to the nearest other member of each phase.", "",
         f"- By direct distance, {len(other)} of {len(rows)} assemblages are closer to another phase's "
         f"profile than to their own: " + (", ".join(f"{r['name']} ({r['own']}, closer to {r['best']})"
                                                     for r in other) or "none") + ".",
         f"- With sampling carried through, {len(clear_own)} fit their own phase best with probability "
         f"above {EVEN_BAND[1]}, {len(clear_other)} fit another phase best (probability of their own below "
         f"{EVEN_BAND[0]}), and {len(even)} are about even between two phases "
         f"({EVEN_BAND[0]} to {EVEN_BAND[1]}): " + (", ".join(r["name"] for r in even) or "none") + ".",
         f"- {len(geo_other)} assemblages lie nearer a member of another phase than any other member of "
         f"their own: " + (", ".join(f"{r['name']} ({r['own']}, nearest {r['nearest_phase']})"
                                     for r in geo_other) or "none") + ".",
         "",
         "| assemblage | sherds | phase | ceramic distance to " + " / ".join(phases)
         + " | closest | P(own) | km to nearest member of " + " / ".join(phases) + " |",
         "|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: r["p_own"]):
        L.append(f"| {r['name']} | {r['n']} | {r['own']}{' (nearest area)' if r['derived'] else ''} | "
                 + " / ".join(f"{r['dist'][q]:.3f}" for q in phases)
                 + f" | {r['best']} | {r['p_own']:.2f} | "
                 + " / ".join(f"{r['near'][q]:.0f}" for q in phases) + " |")
    L.append("")

    # Sensitivity: the two southern-arm assemblages counted as Kent.
    t76 = importlib.import_module("76_phase_recovery")
    t86 = importlib.import_module("86_partition_posterior")
    rev = importlib.import_module("47_revision_analysis")
    data = rev.load_sets()["basin"]
    if [str(n) for n in data["names"]] != names:
        raise RuntimeError("assemblage order differs between loaders")
    moved = labels.copy()
    for site in MOVED:
        moved[names.index(site)] = "Kent"
    G = t76.unit_scale(pts)
    C = t76.unit_scale(t76.composition_features(p, "chisq"))
    L += ["## If Big Eddy and Castile Landing are counted as Kent", "",
          "Recovery is the adjusted Rand index over all 28 assemblages (chi-square composition); the "
          "phase-line comparison is analysis 86's, with its seeds (2,000 posterior draws, 300 "
          "alternative divisions of each kind).", "",
          "| | Phillips's areas | the two counted as Kent |", "|---|---|---|"]
    cols = []
    for lab in (labels, moved):
        pl = sorted(set(lab))
        idx = np.array([pl.index(q) for q in lab])
        r = t86.compare(m, idx, pts, rev, t74, 2000, 300, 86000, dist=data["d"])
        cols.append({
            "phase sizes (Kent / Parkin / Walls)": " / ".join(str(int((lab == q).sum())) for q in pl),
            "site map alone, k-means": f"{t76.ari(idx, t76.cluster(G, 3, 'kmeans', 0, mf)):.3f}",
            "site map alone, Ward": f"{t76.ari(idx, t76.cluster(G, 3, 'ward', 0, mf)):.3f}",
            "map plus composition (weight 0.5), k-means":
                f"{t76.ari(idx, t76.cluster(np.column_stack([G, np.sqrt(0.5) * C]), 3, 'kmeans', 0, mf)):.3f}",
            "composition alone, k-means": f"{t76.ari(idx, t76.cluster(C, 3, 'kmeans', 0, mf)):.3f}",
            "F_ST, posterior median (95% CrI)":
                f"{r['phase'][0]:.4f} ({r['phase'][1]:.4f} to {r['phase'][2]:.4f})",
            "median F_ST, random-center / compact divisions": f"{r['alt'][0]:.4f} / {r['opt'][0]:.4f}",
            "P(lines separate better than random-center / compact)": f"{r['p_alt']:.2f} / {r['p_opt']:.2f}",
            "boundary excess, posterior median": f"{r['be_phase'][0]:+.1f}",
            "P(larger boundary excess than random-center / compact)": f"{r['be_p_alt']:.2f} / {r['be_p_opt']:.2f}",
        })
    for key in cols[0]:
        L.append(f"| {key} | {cols[0][key]} | {cols[1][key]} |")
    L.append("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
