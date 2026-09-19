"""59_partition_sensitivity.py - how much of the reported F_ST is the partition?

F15 measurements M1 and M3, per
`docs/superpowers/plans/2026-09-02-f15-partition-dependence.md`.

THE QUESTION. Every cultural F_ST this project reports is measured between
groups produced by k-means on coordinates, with k chosen by a silhouette score
over 2 to 6 and a fixed seed of 7. Three researcher degrees of freedom sit under
the paper's central quantity: the algorithm, the criterion that picks k, and the
seed. The paper's own thesis is that phases are arbitrary chunks of a continuous
field, so the same critique applies to its analytical partition, and it is
better measured by us than by a reviewer.

  M1  sensitivity to k: refit at k = 2..6, data, prior and seed held fixed.
  M3  sensitivity to the seed: refit at the selected k across many k-means seeds.

NOT A TEST (rule 18). "How much does the answer depend on a choice we made" is a
sensitivity question. Its Bayesian form is the spread of the posterior across
those choices, reported as a spread. No null is constructed and no p-value is
computed.

WHAT WOULD COUNT AS A PROBLEM, stated before the numbers arrive so it cannot be
adjusted afterwards. The reported F_ST is 0.0179 with a 95 percent credible
interval of about [0.0157, 0.0201], a width of 0.004. If varying k moves the
posterior median by MORE than that width, then the choice of k matters more than
the data's own uncertainty, and the partition-based quantity is carrying a
researcher choice as if it were a result.

Usage: .venv/bin/python analyses/59_partition_sensitivity.py [--fast]
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

OUT_MD = ROOT / "output" / "findings" / "partition_sensitivity.md"
K_GRID = [2, 3, 4, 5, 6]
SEEDS = list(range(1, 21))
PRIOR = ("beta", 1.0, 10.0)          # the adopted primary prior (D-31)
FULL = dict(draws=1500, tune=1500, chains=4, target_accept=0.95)
FAST = dict(draws=400, tune=600, chains=2, target_accept=0.95)
REPORTED_WIDTH = 0.0044              # width of the reported 95% CI, [0.0157, 0.0201]


def basin_inputs():
    """Basin ids, centred coordinates and counts, exactly as analysis 43 uses."""
    import make_figures as mf
    a07 = importlib.import_module("07_refined_empirical")
    inp = a07.prepare_inputs()
    members = mf._basin_members("curated")
    ids_all = [str(a) for a in inp.have_coords_ids]
    keep = [i for i, a in enumerate(ids_all) if a in members]
    coords = np.asarray(inp.cc_c)[keep]
    counts = inp.counts.loc[[ids_all[i] for i in keep]].to_numpy(float)
    return coords, counts


def group_counts(counts, labels):
    uniq = np.unique(labels)
    gc = np.array([counts[labels == c].sum(axis=0) for c in uniq])
    return gc[gc.sum(1) > 0]


def fit(gc, cfg, seed=0):
    import arviz as az
    from mls_emergence.inference import sample_fst, gini_simpson_summary
    idata = sample_fst(gc, random_seed=seed, f_prior=PRIOR, **cfg)
    g = gini_simpson_summary(idata, gc, seed=0)
    s = np.asarray(g["gst_samples"], float)
    flat = lambda d: np.concatenate([np.atleast_1d(v.values).ravel()
                                     for v in d.data_vars.values()])
    return dict(med=float(np.median(s)),
                lo=float(np.percentile(s, 2.5)), hi=float(np.percentile(s, 97.5)),
                rhat=float(np.max(flat(az.rhat(idata)))),
                ess=float(np.min(flat(az.ess(idata)))),
                ndiv=int(idata.sample_stats["diverging"].sum()),
                ntot=int(idata.sample_stats["diverging"].size))


def main(fast=False):
    cfg = FAST if fast else FULL
    seeds = SEEDS[:5] if fast else SEEDS
    from mls_emergence.signatures.assortativity import _kmeans_labels
    a07 = importlib.import_module("07_refined_empirical")
    coords, counts = basin_inputs()
    print(f"basin: {coords.shape[0]} assemblages, {counts.shape[1]} classes")

    sil = {k: a07.silhouette_mean(coords, _kmeans_labels(coords, k, seed=7))
           for k in K_GRID}
    k_sel = max(sil, key=sil.get)
    print(f"silhouette selects k = {k_sel}  "
          + ", ".join(f"k={k}:{v:.4f}" for k, v in sorted(sil.items())))

    print("\nM1: sensitivity to k")
    m1 = {}
    for k in K_GRID:
        gc = group_counts(counts, _kmeans_labels(coords, k, seed=7))
        if gc.shape[0] < 2:
            print(f"  k={k}: fewer than two non-empty groups, skipped"); continue
        m1[k] = fit(gc, cfg)
        r = m1[k]
        print(f"  k={k}: F_ST {r['med']:.4f} [{r['lo']:.4f}, {r['hi']:.4f}]  "
              f"R-hat {r['rhat']:.4f} ESS {r['ess']:.0f} div {r['ndiv']}")

    print(f"\nM3: sensitivity to the k-means seed at k = {k_sel}")
    m3 = {}
    for sd in seeds:
        gc = group_counts(counts, _kmeans_labels(coords, k_sel, seed=sd))
        if gc.shape[0] < 2:
            continue
        m3[sd] = fit(gc, cfg)
    meds3 = np.array([v["med"] for v in m3.values()])
    print(f"  {len(m3)} seeds: F_ST median ranges {meds3.min():.4f} to "
          f"{meds3.max():.4f} (spread {meds3.max()-meds3.min():.4f})")

    meds1 = np.array([m1[k]["med"] for k in sorted(m1)])
    spread1 = float(meds1.max() - meds1.min())
    spread3 = float(meds3.max() - meds3.min())

    L = ["# How much of the reported F_ST is the partition we drew?", "",
         f"Produced by `analyses/59_partition_sensitivity.py` "
         f"({'FAST' if fast else 'full'}). F15 measurements M1 and M3. "
         f"Basin, {coords.shape[0]} assemblages, Beta(1,10) prior throughout. "
         f"Not a test: this is the spread of the posterior across choices we "
         f"made, reported as a spread.", "",
         f"**Pre-stated threshold.** The reported F_ST is 0.0179 with a 95 "
         f"percent interval of width {REPORTED_WIDTH:.4f}. A choice that moves "
         f"the median by more than that matters more than the data's own "
         f"uncertainty. This was written down before the numbers were seen.", "",
         "## M1. Sensitivity to k", "",
         "| k | silhouette | F_ST median | 95% CI | R-hat | min ESS | divergences |",
         "|---|---|---|---|---|---|---|"]
    for k in sorted(m1):
        r = m1[k]
        star = " **(selected)**" if k == k_sel else ""
        L.append(f"| {k}{star} | {sil[k]:.4f} | **{r['med']:.4f}** | "
                 f"[{r['lo']:.4f}, {r['hi']:.4f}] | {r['rhat']:.4f} | "
                 f"{r['ess']:.0f} | {r['ndiv']}/{r['ntot']} |")
    L += ["",
          f"Across k = {min(m1)} to {max(m1)} the posterior median spans "
          f"{meds1.min():.4f} to {meds1.max():.4f}, a spread of "
          f"**{spread1:.4f}**, against a reported interval width of "
          f"{REPORTED_WIDTH:.4f}. "
          + (f"**That is {spread1/REPORTED_WIDTH:.1f} times the data's own "
             f"uncertainty**, so the choice of k matters more than the evidence "
             f"the data carry about F_ST."
             if spread1 > REPORTED_WIDTH else
             f"That is {spread1/REPORTED_WIDTH:.2f} times the data's own "
             f"uncertainty, so the choice of k is not the dominant term."), "",
          f"## M3. Sensitivity to the k-means seed at k = {k_sel}", "",
          f"{len(m3)} seeds. Posterior medians span {meds3.min():.4f} to "
          f"{meds3.max():.4f}, a spread of **{spread3:.4f}** "
          f"({spread3/REPORTED_WIDTH:.2f} times the reported interval width).",
          "",
          ("The seed is not a meaningful degree of freedom: the k-means helper "
           "takes the best of several initialisations, and at this k the "
           "partition is stable across seeds. One of the three degrees of "
           "freedom named in F15 can be struck."
           if spread3 < 0.2 * REPORTED_WIDTH else
           "**The seed matters.** A fixed `seed=7` is therefore carrying part "
           "of the reported number, which is not defensible and must be "
           "reported or removed."), "",
          "## What this means for M2", ""]
    L += ([f"k-dependence is the dominant term, so the partition-at-fixed-grain "
           f"ensemble (M2) is worth running: the question becomes whether the "
           f"specific partition is special among partitions of the same grain, "
           f"or just one draw."]
          if spread1 > REPORTED_WIDTH else
          [f"k-dependence is smaller than the data's own uncertainty, which "
           f"weakens the case for running M2. The partition is not obviously "
           f"carrying the result, and the honest framing may simply be to "
           f"report the dependence measured here alongside the number."])
    L += ["", "## Diagnostics (rule 16)", ""]
    for k in sorted(m1):
        r = m1[k]
        L.append(f"- k = {k}: R-hat {r['rhat']:.4f}, min ESS {r['ess']:.0f}, "
                 f"divergences {r['ndiv']}/{r['ntot']}")
    worst = max((v["rhat"] for v in m3.values()), default=float("nan"))
    L.append(f"- seed sweep: worst R-hat {worst:.4f}, "
             f"total divergences {sum(v['ndiv'] for v in m3.values())}")
    L.append("")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    ks = sorted(m1)
    ax.errorbar(ks, [m1[k]["med"] for k in ks],
                yerr=[[m1[k]["med"] - m1[k]["lo"] for k in ks],
                      [m1[k]["hi"] - m1[k]["med"] for k in ks]],
                fmt="o-", color=OI_BLUE, ms=4, capsize=3, label="posterior by k")
    ax.axhline(0.0179, color=OI_VERMIL, ls="--", lw=1, label="reported (k=3)")
    ax.set_xlabel("number of spatial clusters k"); ax.set_ylabel("cultural $F_{ST}$")
    ax.set_xticks(ks); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); save(fig, "fig_partition_sensitivity")
    print(f"\nwrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
