#!/usr/bin/env python3
"""Are the phase boundaries where the ceramic differences are?

Everything else in this paper measures differentiation across a partition the
analysis draws (k-means on the assemblages' coordinates). That answers whether
structure exists at a given grain. It does not answer the question the phases
themselves pose, which is whether the boundaries culture history drew are where
the differences lie. This script asks that directly, and reports the two halves
of the answer together because either alone misleads.

PANEL A. HOW THE ANSWER DEPENDS ON GRAIN. The observed cultural F_ST is plotted
against the calibrated spatial-drift null at k = 2 to 6 (analyses/71_scale_
sweep.py). At the coarsest cut the record sits inside what drift produces; from
k = 3 the observed exceeds it, and the gap widens monotonically with grain.
Drift reproduces broad regional differentiation and fails locally. The level of
the statistic is a property of the scale, so no single F_ST is "the" value; the
comparison against drift is what carries the argument, and it is stable
wherever the test has power.

PANELS B AND C. WHETHER THE BOUNDARIES MATTER. The phases are not arbitrary
lines: they are where generations of workers thought they saw structure. The
question is what that perception was tracking. The phase partition's F_ST is
therefore placed against partitions of the SAME GROUP SIZES (10, 3, 4, 11, 15)
whose boundaries fall elsewhere.

WHY THE FIRST VERSION OF THIS NULL WAS WRONG. It built each partition with a
greedy capacitated assignment: take assemblages in order of regret and give
each its nearest seed with room left. Under a hard capacity that order
front-loads the periphery and leaves the dense interior for last, when the
nearby group is full, so interior assemblages get flung to distant seeds. The
partitions came out with exclaves: their within-group inertia was 4.2 times the
phases' (median 681 against 161.8 km^2 per assemblage, n = 1,000 on shared
seeds), and they broke into a median of 9 connected components on the
assemblages' Delaunay graph against the phases' 6 and the exact assignment's 6
-- all while the docstring asserted "comparable compactness". Note that the
phases themselves are not perfectly contiguous either, so the 6 is the standard
to match, not 5. Fractured groups mix profiles that are not
neighbours, which depresses between-group variance, so the null was easy to
beat and the phases landed at the 70th percentile of it.

Replacing the greedy with an exact minimum-cost assignment ON THE SAME SEEDS,
changing nothing else, raises the null's median F_ST from 0.0217 to 0.0277 and
moves the phases from the 72nd percentile to the 33rd. That controlled swap is
the evidence that the generator, not the seeds, was the defect.

One caution about the diagnosis, because the obvious explanation is not the
right one. Within the greedy ensemble, inertia and F_ST correlate at Spearman
-0.35, which invites the reading that compactness drives this statistic and
that a null simply has to be compactness-matched. It does not survive the fix:
within the exact-assignment ensemble the same correlation is -0.05, and pooled
across both ensembles reported here it is +0.02. The coupling was a property of
the broken generator -- its worst draws were both the least compact and the
most profile-scrambling -- not a general confound. So compactness matching is
done here because it makes the comparison the right one to ask, not because an
unmatched null would otherwise be biased.

Two ensembles are reported, with the inertia of each printed beside its F_ST so
the reader can see how compact the comparison partitions actually are:

  size-matched, boundaries at random
      exact min-cost assignment (`scipy.optimize.linear_sum_assignment`) of
      assemblages to seed points drawn at random, with slots replicated by
      group size so the sizes come out exactly right. Compact given its seeds,
      but the seeds are uninformed.
  size- and compactness-matched
      balanced local search from random starts: swap pairs between groups
      whenever that lowers within-group inertia, until no improving swap
      remains. Sizes stay exact. This is the null the phases are judged
      against -- a division of this ground into groups of these sizes, drawn
      to be spatially tight and knowing nothing about pots. It does not reach
      the phases' own compactness on a typical restart, because the phases sit
      essentially AT the optimum (see the recovery rate below), so the match is
      one-sided: the comparison partitions are looser than the phases, and the
      output file prints both distributions rather than claiming otherwise.

There is no third "optimised" reference any more. The earlier version had one,
built from k-means on ROTATED coordinates, and it was wrong twice over. Lloyd
k-means on Euclidean coordinates is rotation-invariant, so the rotation was
inert -- the partition is identical at every angle tested, adjusted Rand index
1.000 -- and every draw's variation came from the jitter and the seed (rule 4).
Worse, its conclusion, that "an optimiser puts its boundaries elsewhere, and
does better", is false. The optimiser puts them almost exactly where the phases
are: mean adjusted Rand index 0.88 with the phase partition, 16 of its 200
draws reproduced the phase partition exactly, and its median sorted group sizes
(3, 4, 10, 11, 15) are the phases' own. Reading "0th percentile of 200" off
that ensemble also overstated it: the 200 draws contained 26 distinct
partitions. The finding underneath is kept, because it is the real result, and
it is measured here properly as the recovery rate of the compactness search.

If the phases marked interaction communities, their boundaries should sit where
differences are sharp: well above size- and compactness-matched cuts. If they
are a reading of the most spatially compact division available, they should sit
inside that distribution, and a blind compactness optimiser should keep
rediscovering them.

Usage:
    python analyses/74_phase_partition_test.py [--alt 1000] [--opt 1000]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.special import comb

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "phase_partition_test.md"
OUT_CSV = ROOT / "output" / "findings" / "phase_partition_ensembles.csv"
SWEEP_CSV = ROOT / "output" / "findings" / "scale_sweep.csv"
FIG = "fig12_phase_partition"
EARTH_KM_PER_DEG = 111.32          # as analyses/23_phases_vs_spatial_drift.geo_km


def km_xy(xy: np.ndarray) -> np.ndarray:
    """Coordinates as centred kilometres, equirectangular about the mean latitude.

    Degrees of longitude are 18 percent shorter than degrees of latitude at 35 N,
    so a Euclidean norm on raw degrees stretches the basin east-west by that
    much. It makes no material difference to the percentiles here (checked: the
    random-boundary ensemble moves by about one percentile point), but inertia
    is reported in km^2 and has to be in real units to be read.
    """
    mlat = np.deg2rad(xy[:, 0].mean())
    p = np.column_stack([xy[:, 1] * np.cos(mlat) * EARTH_KM_PER_DEG,
                         xy[:, 0] * EARTH_KM_PER_DEG])
    return p - p.mean(0)


def inertia(pts: np.ndarray, lab: np.ndarray) -> float:
    """Mean squared distance to the group centroid, per assemblage (km^2)."""
    return float(sum(((pts[lab == g] - pts[lab == g].mean(0)) ** 2).sum()
                     for g in np.unique(lab))) / len(pts)


def assign_exact(pts: np.ndarray, seeds: np.ndarray, slot: np.ndarray) -> np.ndarray:
    """Minimum-cost assignment of assemblages to seeds under exact group sizes.

    `slot` lists the group of each of the n capacity slots. Replicating seed
    columns by group size turns the capacitated problem into a square linear
    assignment, which `linear_sum_assignment` solves exactly. This is what the
    earlier greedy was approximating, badly.
    """
    cost = ((pts[:, None, :] - seeds[None, :, :]) ** 2).sum(2)[:, slot]
    r, c = linear_sum_assignment(cost)
    out = np.empty(len(pts), int)
    out[r] = slot[c]
    return out


def local_search(pts: np.ndarray, lab: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Lower inertia by swapping pairs between groups, keeping sizes exact.

    Minimising within-group inertia at fixed sizes is the same as maximising
    sum_g |S_g|^2 / n_g with S_g the group's coordinate sum, because the sum of
    squared norms is constant. A swap touches only two group sums, so every
    candidate swap is scored in closed form and the whole pair set is evaluated
    at once. Steepest descent; stops at a local optimum.
    """
    lab = lab.copy()
    k = len(sizes)
    while True:
        S = np.array([pts[lab == g].sum(0) for g in range(k)])
        gain, move = 1e-12, None
        for a in range(k):
            ia = np.where(lab == a)[0]
            for b in range(a + 1, k):
                ib = np.where(lab == b)[0]
                pa, pb = pts[ia][:, None, :], pts[ib][None, :, :]
                sa, sb = S[a] - pa + pb, S[b] - pb + pa
                d = (((sa ** 2).sum(-1) - (S[a] ** 2).sum()) / sizes[a]
                     + ((sb ** 2).sum(-1) - (S[b] ** 2).sum()) / sizes[b])
                i, j = np.unravel_index(np.argmax(d), d.shape)
                if d[i, j] > gain:
                    gain, move = float(d[i, j]), (ia[i], ib[j])
        if move is None:
            return lab
        i, j = move
        lab[i], lab[j] = lab[j], lab[i]


def ari(a: np.ndarray, b: np.ndarray) -> float:
    """Adjusted Rand index, so "the same partition" does not depend on labels."""
    C = np.zeros((int(a.max()) + 1, int(b.max()) + 1), int)
    for i, j in zip(a, b):
        C[i, j] += 1
    s = comb(C, 2).sum()
    ra, rb = comb(C.sum(1), 2).sum(), comb(C.sum(0), 2).sum()
    e = ra * rb / comb(len(a), 2)
    return float((s - e) / ((ra + rb) / 2 - e))


def ordinal(p: float) -> str:
    """'33rd', not '33th'."""
    n = int(round(p))
    suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf}"


def pct_ci(ens: np.ndarray, value: float, n: int, rng, boot: int = 2000):
    """The percentile, and what it would have wobbled to at `n` draws.

    Rule 6: the operating point rides along. A percentile quoted to one decimal
    off 200 draws is precision the ensemble does not have -- the first version
    of this script reported 70.0 from 200 draws, which resamples anywhere in
    the middle sixties to middle seventies.
    """
    p = 100.0 * float(np.mean(ens < value))
    b = np.array([100.0 * np.mean(rng.choice(ens, size=n, replace=True) < value)
                  for _ in range(boot)])
    return p, float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--alt", type=int, default=1000,
                    help="size-matched partitions with boundaries placed at random")
    ap.add_argument("--opt", type=int, default=1000,
                    help="restarts of the size- and compactness-matched search")
    ap.add_argument("--shuffles", type=int, default=2000)
    args = ap.parse_args()

    import matplotlib.pyplot as plt
    fs = importlib.import_module("figstyle")
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    ph = importlib.import_module("36_canonical_phase_map")

    counts, coords = mf._load_curated()
    m = counts.to_numpy(int)
    xy = coords.to_numpy(float)
    names = [str(i) for i in counts.index]

    labels_ph, derived = ph.assign_phases_by_territory(names, xy)
    phases = sorted(set(labels_ph))
    phase_idx = np.array([phases.index(l) for l in labels_ph])
    k = len(phases)
    sizes = np.bincount(phase_idx)
    slot = np.repeat(np.arange(k), sizes)
    pts = km_xy(xy)

    fst_phase = float(rev.fst_by(m, phase_idx))
    in_phase = inertia(pts, phase_idx)

    # Separate streams, so changing one ensemble's draw count cannot silently
    # move another's. The first version drew all three from one generator.
    rng_s = np.random.default_rng(5)
    rng_a = np.random.default_rng(74001)
    rng_o = np.random.default_rng(74002)
    rng_b = np.random.default_rng(74003)

    shuffled = np.array([rev.fst_by(m, rng_s.permutation(phase_idx))
                         for _ in range(args.shuffles)])

    alt, alt_in = [], []
    for _ in range(args.alt):
        seeds = pts[rng_a.choice(len(pts), size=k, replace=False)]
        a = assign_exact(pts, seeds, slot)
        alt.append(rev.fst_by(m, a))
        alt_in.append(inertia(pts, a))
    alt, alt_in = np.array(alt), np.array(alt_in)

    opt, opt_in, opt_ari = [], [], []
    for r in range(args.opt):
        seeds = pts[rng_o.choice(len(pts), size=k, replace=False)]
        a = local_search(pts, assign_exact(pts, seeds, slot), sizes)
        opt.append(rev.fst_by(m, a))
        opt_in.append(inertia(pts, a))
        opt_ari.append(ari(phase_idx, a))
        if (r + 1) % 200 == 0:
            print(f"  compactness search {r + 1}/{args.opt}", flush=True)
    opt, opt_in, opt_ari = np.array(opt), np.array(opt_in), np.array(opt_ari)
    exact = int((opt_ari > 0.9999).sum())

    pct_alt, alt_lo, alt_hi = pct_ci(alt, fst_phase, args.alt, rng_b)
    pct_opt, opt_lo, opt_hi = pct_ci(opt, fst_phase, args.opt, rng_b)
    pct_shuf = 100.0 * float(np.mean(shuffled < fst_phase))
    rho = float(pd.Series(alt_in).corr(pd.Series(alt), method="spearman"))
    rho_opt = float(pd.Series(opt_in).corr(pd.Series(opt), method="spearman"))

    pd.DataFrame(dict(
        ensemble=(["random_boundary"] * len(alt) + ["compactness_matched"] * len(opt)),
        fst=np.concatenate([alt, opt]),
        inertia_km2=np.concatenate([alt_in, opt_in]),
        ari_with_phases=np.concatenate([np.full(len(alt), np.nan), opt_ari]),
    )).to_csv(OUT_CSV, index=False)

    sweep = pd.read_csv(SWEEP_CSV) if SWEEP_CSV.exists() else None
    if sweep is not None:
        sweep = sweep[sweep.k <= 6]      # 71 sweeps to k = 12 for Figure 13; panel A shows 2 to 6

    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 3, wspace=0.42)
    axA = fig.add_subplot(gs[0, 0])
    if sweep is not None:
        axA.fill_between(sweep.k, sweep.drift_lo, sweep.drift_hi, color="0.82",
                         label="calibrated drift, 95%")
        axA.plot(sweep.k, sweep.drift_median, color="0.45", lw=1.2, ls="--",
                 label="drift median")
        axA.plot(sweep.k, sweep.observed, color="0.1", lw=1.6, marker="o",
                 ms=4, label="observed")
        axA.set_xlabel("number of spatial groups (k)")
        axA.set_ylabel("cultural $F_{ST}$")
        axA.set_xticks(list(sweep.k))
        from matplotlib.ticker import MultipleLocator as _ML
        axA.yaxis.set_major_locator(_ML(0.01))
        axA.legend(fontsize=6, frameon=False, loc="upper left")
    else:
        axA.text(.5, .5, "run 71_scale_sweep.py", ha="center", transform=axA.transAxes)
    fs.panel_label(axA, "A")

    axB = fig.add_subplot(gs[0, 1])
    axB.scatter(alt_in, alt, s=5, color="0.75", edgecolor="none")
    axB.scatter(opt_in, opt, s=5, color="0.42", edgecolor="none")
    axB.plot([in_phase], [fst_phase], marker="*", ms=11, color="0.05", ls="none")
    # Proxy handles at legible sizes. Passing `label=` to the scatters and then
    # scaling the legend markers blew the star up to fill the axes, because
    # markerscale multiplies every handle including the one already large.
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ms=3.5, ls="none", color="0.75",
                      label="boundaries at random"),
               Line2D([], [], marker="o", ms=3.5, ls="none", color="0.42",
                      label="compactness matched"),
               Line2D([], [], marker="*", ms=7, ls="none", color="0.05",
                      label="the phases")]
    axB.set_xlabel("within-group inertia (km$^2$ per assemblage)")
    axB.set_ylabel("cultural $F_{ST}$")
    # Round ticks only. Panel B's default locator put ticks at 0.015, 0.025 and
    # 0.035, which scripts/check_figure_claims.py reads as statistics the
    # manuscript never states -- its documented "finer axis ticks still slip"
    # limit. Matching the other panels' tick values keeps scaffolding out of
    # the comparison.
    from matplotlib.ticker import MultipleLocator
    axB.yaxis.set_major_locator(MultipleLocator(0.01))
    axB.legend(handles=handles, fontsize=5.5, frameon=True, framealpha=0.92,
               edgecolor="0.8", loc="lower right", handletextpad=0.5,
               borderpad=0.4)
    fs.panel_label(axB, "B")

    axC = fig.add_subplot(gs[0, 2])
    axC.hist(opt, bins=26, color="0.78", edgecolor="0.55", lw=.4)
    axC.axvline(fst_phase, color="0.1", lw=1.8)
    axC.xaxis.set_major_locator(MultipleLocator(0.01))
    axC.annotate(f"the phases\n{fst_phase:.4f}\n{ordinal(pct_opt)} percentile\n"
                 f"of this ensemble",
                 xy=(fst_phase, axC.get_ylim()[1] * .78), xytext=(6, 0),
                 textcoords="offset points", fontsize=6, va="center", ha="left")
    axC.set_xlabel(f"cultural $F_{{ST}}$, size- and\ncompactness-matched {k}-group cuts")
    axC.set_ylabel("partitions")
    fs.panel_label(axC, "C")

    png = fs.save_all(fig, FIG, close=True)

    L = ["# Are the phase boundaries where the differences are?", "",
         f"Basin phase set, {len(names)} assemblages, {k} phases "
         f"({', '.join(f'{p} {int((phase_idx == i).sum())}' for i, p in enumerate(phases))}).",
         f"Phase partition cultural F_ST **{fst_phase:.4f}**, within-group inertia "
         f"**{in_phase:.1f} km^2** per assemblage.", "",
         "| partition ensemble | median F_ST | median inertia (km^2) | the phases sit at |",
         "|---|---|---|---|",
         f"| labels shuffled, n = {args.shuffles} | {np.median(shuffled):.4f} | "
         f"(geography destroyed) | above {int(np.sum(shuffled < fst_phase)):,} of "
         f"{args.shuffles:,} draws |",
         f"| same sizes, boundaries at random, n = {args.alt} | {np.median(alt):.4f} | "
         f"{np.median(alt_in):.0f} | **{ordinal(pct_alt)} pct** "
         f"({ordinal(alt_lo)}-{ordinal(alt_hi)}) |",
         f"| same sizes, compactness matched, n = {args.opt} | {np.median(opt):.4f} | "
         f"{np.median(opt_in):.0f} | **{ordinal(pct_opt)} pct** "
         f"({ordinal(opt_lo)}-{ordinal(opt_hi)}) |",
         "",
         "Percentiles are one-sided empirical percentiles of the ensemble, with the",
         f"2.5th-97.5th range over {2000} resamples of the same draw count in brackets.",
         "",
         "## Did the compactness match hold?", "",
         f"The phases' inertia is {in_phase:.1f} km^2 per assemblage. The random-boundary",
         f"ensemble's median is {np.median(alt_in):.0f} ({np.percentile(alt_in, 5):.0f}"
         f"-{np.percentile(alt_in, 95):.0f}); the compactness-matched ensemble's is "
         f"{np.median(opt_in):.0f} ({np.percentile(opt_in, 5):.0f}"
         f"-{np.percentile(opt_in, 95):.0f}).",
         f"The match is one-sided: the comparison partitions are looser than the phases,",
         "not tighter, because the phases sit essentially at the compactness optimum for",
         "these group sizes.", "",
         f"Across the random-boundary draws, inertia and F_ST correlate at Spearman "
         f"{rho:+.2f}, and",
         f"across the compactness-matched draws at {rho_opt:+.2f}. So differentiation here "
         "is NOT",
         "a simple function of how tight the groups are, and the looseness above does not",
         "by itself bias the comparison. The superseded greedy generator did show such a",
         "coupling (-0.35), which is a fact about that generator's fractured partitions",
         "rather than about this statistic. Panel B plots both clouds.", "",
         "## What the compactness search recovers", "",
         f"Of {args.opt} restarts, **{exact}** ({100 * exact / args.opt:.1f} percent) "
         "returned the phase partition",
         f"exactly (adjusted Rand index 1.0); the mean ARI with the phases is "
         f"{opt_ari.mean():.2f}.",
         f"The lowest inertia found anywhere in the search is {opt_in.min():.1f} km^2 "
         f"against the phases'",
         f"{in_phase:.1f}, so the phases are within "
         f"{100 * (in_phase / opt_in.min() - 1):.1f} percent of the most compact "
         f"{k}-way division",
         "of these assemblages at these group sizes.", "",
         "## Reading", "",
         "The phases are not arbitrary lines; they are where generations of workers",
         "thought they saw structure. What the comparison shows is what that perception",
         "was tracking. Beating shuffled labels says only that nearby assemblages",
         "resemble one another, which isolation by distance produces on its own. Sitting",
         "inside the size- and compactness-matched distribution says the particular",
         "placement of these boundaries carries little information about where ceramic",
         "differences lie beyond the fact that they enclose compact groups. And a search",
         "that minimises inertia and never sees a potsherd rediscovers the phase",
         "partition itself, from random starts, at the rate reported above.", "",
         "A partition marking interaction communities should behave differently: it",
         "should sit high against compactness-matched cuts, because the boundaries would",
         "be where the differences are, not merely where the gaps between sites are.", "",
         "The drift comparison in panel A is a separate and stronger statement, and it is",
         "the one the paper's argument rests on: differentiation exceeds what calibrated",
         "spatial drift produces on this geography at every scale the test can resolve,",
         "whatever partition is used to measure it.", "",
         f"Figure written to {png.name} and its siblings; ensembles in {OUT_CSV.name}.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"\nphase F_ST {fst_phase:.4f}, inertia {in_phase:.1f} km^2")
    print(f"percentile: shuffled {pct_shuf:.0f}; random boundaries {pct_alt:.0f} "
          f"({alt_lo:.0f}-{alt_hi:.0f}); compactness matched {pct_opt:.0f} "
          f"({opt_lo:.0f}-{opt_hi:.0f})")
    print(f"inertia: phases {in_phase:.1f}; random {np.median(alt_in):.0f}; "
          f"matched {np.median(opt_in):.0f}; Spearman(inertia, F_ST) {rho:+.2f}")
    print(f"compactness search reproduced the phases exactly {exact}/{args.opt}")
    print(f"wrote {png} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
