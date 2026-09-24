"""60_partition_ensemble.py - is OUR partition special, or just one of many?

F15 measurement M2, per
`docs/superpowers/plans/2026-09-02-f15-partition-dependence.md`. Run because M1
met its trigger: the spread of the posterior median across k = 2 to 6 exceeds
the reported 95 percent interval width. The current values are read from M1's
own output (`output/findings/partition_sensitivity.md`) at write time rather
than typed here; the 2026-09-03 version of this docstring carried the numbers of
the superseded 34-assemblage set (0.0147 to 0.0347 against 0.0044).

THE GRAIN. K below is the grain at which the reported F_ST is defined, which is
the k the silhouette criterion selects in M1: k = 2 on the 28-assemblage matrix
(restated 2026-09-22; it was 3 on the superseded set). The k-means partition at
that grain is the one `output/bayesian_fst.md` reports.

THE QUESTION M1 LEAVES OPEN. M1 shows the GRAIN matters. It does not say whether,
at a fixed grain, the particular partition k-means found is special. If an
arbitrary partition of the same grain gives the same F_ST, then our partition
carries no information beyond how finely it chops the field, and a quantity
defined on it is measuring the chopping.

THE ENSEMBLE. Random Voronoi partitions into k cells: draw k seed points
uniformly in the bounding box of the assemblage coordinates and assign each
assemblage to its nearest seed. Voronoi is the right family because k-means
produces a Voronoi partition too, so the comparison holds the FORM of the
partition fixed and varies only which one. Cells with fewer than two assemblages
are rejected and redrawn, since F_ST is undefined on them.

NOT A TEST (rule 18). No null hypothesis, no p-value. This reports where our
partition's posterior median sits within the ensemble of medians, as a
descriptive location, and the spread of the ensemble as a sensitivity.

WHAT WOULD COUNT AS A PROBLEM, written before the numbers. If the k-means
partition's median sits unremarkably inside the ensemble, say between the 20th
and 80th percentile, then it is not a distinguished choice and the F_ST defined
on it is one draw from a family of equally defensible answers. If it sits in a
tail, the partition is doing something the alternatives do not, and the current
framing has a measured defence.

Usage: .venv/bin/python analyses/60_partition_ensemble.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from figstyle import OI_BLUE, OI_VERMIL, save  # noqa: E402

a59 = importlib.import_module("59_partition_sensitivity")

OUT_MD = ROOT / "output" / "findings" / "partition_ensemble.md"
K = 2                                # silhouette-selected grain, analysis 59 M1
N_PART = 30 if "--fast" in sys.argv else 150
PRIOR = ("beta", 1.0, 10.0)
CFG = dict(draws=600, tune=800, chains=2, target_accept=0.95)


def random_voronoi(coords, k, rng, min_per_cell=2, tries=500):
    """A random contiguous partition: nearest-seed assignment to k random seeds."""
    lo, hi = coords.min(0), coords.max(0)
    for _ in range(tries):
        seeds = rng.uniform(lo, hi, size=(k, coords.shape[1]))
        d = ((coords[:, None, :] - seeds[None, :, :]) ** 2).sum(-1)
        lab = d.argmin(1)
        if len(np.unique(lab)) == k and np.all(np.bincount(lab, minlength=k) >= min_per_cell):
            return lab
    return None


def _ordinal(x):
    n = int(round(x))
    suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf}"


def _m1_spread():
    """The across-k spread M1 measured, read from its own output (rule 1)."""
    import re
    md = (ROOT / "output" / "findings" / "partition_sensitivity.md").read_text(encoding="utf-8")
    m = re.search(r"a spread of \*\*([0-9.]+)\*\*, against a reported interval width", md)
    if not m:
        raise RuntimeError("partition_sensitivity.md does not carry M1's spread; run analysis 59 first")
    return float(m.group(1))


def main(fast=False):
    from mls_emergence.signatures.assortativity import _kmeans_labels
    coords, counts = a59.basin_inputs()
    print(f"basin: {coords.shape[0]} assemblages, k = {K}, {N_PART} random partitions")

    ours = a59.fit(a59.group_counts(counts, _kmeans_labels(coords, K, seed=7)), CFG)
    print(f"k-means partition: F_ST {ours['med']:.4f} [{ours['lo']:.4f}, {ours['hi']:.4f}]")

    rng = np.random.default_rng(11)
    meds, bad = [], 0
    for i in range(N_PART):
        lab = random_voronoi(coords, K, rng)
        if lab is None:
            bad += 1; continue
        gc = a59.group_counts(counts, lab)
        if gc.shape[0] < 2:
            bad += 1; continue
        meds.append(a59.fit(gc, CFG, seed=i)["med"])
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{N_PART} ...")
    meds = np.array(meds)
    pct = float((meds < ours["med"]).mean() * 100)
    print(f"\nensemble: n={len(meds)}, median {np.median(meds):.4f}, "
          f"range {meds.min():.4f} to {meds.max():.4f}")
    print(f"our partition sits at the {pct:.0f}th percentile of the ensemble")

    inside = 20.0 <= pct <= 80.0
    m1_spread = _m1_spread()
    L = ["# Is our partition special, or just one of many?", "",
         f"Produced by `analyses/60_partition_ensemble.py` "
         f"({'FAST' if fast else 'full'}). F15 measurement M2. Basin, "
         f"{coords.shape[0]} assemblages, k = {K}, Beta(1,10) prior, "
         f"{len(meds)} random contiguous (Voronoi) partitions"
         + (f" ({bad} rejected for empty or singleton cells)." if bad else "."), "",
         "Not a test: no null and no p-value. This locates our partition inside "
         "an ensemble of equally defensible ones and reports the spread.", "",
         "## Result", "",
         "| | cultural F_ST |", "|---|---|",
         f"| **the k-means partition we use** | **{ours['med']:.4f}** "
         f"[{ours['lo']:.4f}, {ours['hi']:.4f}] |",
         f"| ensemble median | {np.median(meds):.4f} |",
         f"| ensemble 5th to 95th percentile | {np.percentile(meds,5):.4f} to {np.percentile(meds,95):.4f} |",
         f"| ensemble full range | {meds.min():.4f} to {meds.max():.4f} |",
         f"| **our partition's position in the ensemble** | **{_ordinal(pct)} percentile** |",
         "", "## Reading", ""]
    if inside:
        L += [f"Our partition sits at the {_ordinal(pct)} percentile, which is "
              f"unremarkable. At this number of clusters the specific partition k-means found "
              f"is **one draw from a family of equally defensible partitions**, "
              f"and the F_ST defined on it is not a property the data single out.",
              "",
              f"Combined with M1, where the median moved {m1_spread:.4f} across k against "
              f"a data-driven interval width of {a59.REPORTED_WIDTH:.4f}, the picture is that the "
              f"reported F_ST is set both by **how finely the field is divided** and by "
              f"**where the cuts fall**: divisions with the same number of groups span "
              f"{np.percentile(meds,5):.4f} to {np.percentile(meds,95):.4f} (5th to 95th percentile). "
              f"A quantity this sensitive to the cutting cannot be read on its own."]
    else:
        L += [f"Our partition sits at the {_ordinal(pct)} percentile of the ensemble, "
              f"outside the 20th to 80th band set in advance. It is doing "
              f"something the alternatives do not, and the current framing has a "
              f"measured defence rather than only an assertion."]
    L += ["", "## What follows", "",
          "The partition-free alternative is already in hand and needs no such "
          "defence: the spatial GP contains no partition at all. Its predictive "
          "check on the observed F_ST is in `output/findings/gp_posterior_predictive.md` "
          "and its spatial share in `output/findings/gp_basin_fit.md`; they are "
          "not copied here. The framing decision this measurement feeds is recorded in "
          "`docs/superpowers/plans/2026-09-02-f15-partition-dependence.md`.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    ax.hist(meds, bins=25, color=OI_BLUE, alpha=0.7)
    ax.axvline(ours["med"], color=OI_VERMIL, lw=1.6,
               label=f"our k-means partition ({pct:.0f}th pct)")
    ax.set_xlabel("cultural $F_{ST}$ under a random contiguous partition")
    ax.set_ylabel("count"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); save(fig, "fig_partition_ensemble")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
