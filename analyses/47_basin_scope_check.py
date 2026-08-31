"""47_basin_scope_check.py - does the Bayesian F_ST fit run on the basin?

Analyses 43 and 44 call `43_bayesian_fst.basin_group_counts(inp)` on the object
returned by `07_refined_empirical.prepare_inputs()`, and group by
`inp.cluster_of` over `inp.have_coords_ids`. Neither restricts to the canonical
drainage-basin membership in `data/processed/basin_members_curated.txt`.

This script measures the consequence: it refits the Balding-Nichols model on
(A) the set analysis 43 actually uses and (B) the basin, with clusters
re-selected on the basin's own coordinates by the same silhouette rule analysis
07 uses, and reports both.

The comparison is the point. Neither fit is assumed correct in advance; the
question is how far apart they are and which one the label "St. Francis basin"
belongs to.

Usage: .venv/bin/python analyses/47_basin_scope_check.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

OUT_MD = ROOT / "output" / "findings" / "basin_scope_check.md"
SEED = 0
CFG = dict(draws=2000, tune=2000, chains=4)


def cluster_labels(coords, seed=7, kmin=2, kmax=6):
    """Analysis 07's rule: k maximizing mean silhouette over k in [kmin, kmax]."""
    # silhouette_mean is not in src/; it is copy-pasted into five analysis
    # scripts (05, 06, 07, 09, make_figures). Import analysis 07's copy, since
    # that is the one the fit under test goes through.
    from mls_emergence.signatures.assortativity import _kmeans_labels
    silhouette_mean = importlib.import_module("07_refined_empirical").silhouette_mean
    sil = {k: silhouette_mean(coords, _kmeans_labels(coords, k, seed=seed))
           for k in range(kmin, kmax + 1)}
    k_use = max(sil, key=sil.get)
    return _kmeans_labels(coords, k_use, seed=seed), k_use, sil


def group_counts(counts_arr, labels):
    uniq = np.unique(labels)
    return np.array([counts_arr[labels == c].sum(axis=0) for c in uniq])


def fit(gc, label):
    from mls_emergence.inference import sample_fst, fst_summary, gini_simpson_summary
    idata = sample_fst(gc, random_seed=SEED, **CFG)
    s = fst_summary(idata)
    g = gini_simpson_summary(idata, gc, seed=SEED)
    ndiv = int(idata.sample_stats["diverging"].sum())
    print(f"[{label}] BN median {s['fst_median']:.4f}  GS median {g['gst_median']:.4f}")
    return s, g, ndiv


def main():
    import make_figures as mf
    a07 = importlib.import_module("07_refined_empirical")
    a43 = importlib.import_module("43_bayesian_fst")

    inp = a07.prepare_inputs()
    ids_used = [str(i) for i in inp.have_coords_ids]
    basin = mf._basin_members("curated")
    ids_basin = [i for i in ids_used if i in basin]

    # (A) the AS-PUBLISHED path. This must name scope="region" explicitly: since
    # the F18 fix, basin_group_counts() defaults to the basin, so calling it
    # bare would silently make column A a duplicate of column B while the labels
    # still claimed 55 assemblages. Caught 2026-08-31 when both columns returned
    # identical numbers under mismatched headings.
    gc_a, _ = a43.basin_group_counts(inp, scope="region")
    lab_a = np.array([inp.cluster_of[a] for a in inp.have_coords_ids])
    k_a = len(np.unique(lab_a))

    # (B) the basin, clusters re-selected on the basin's own coordinates
    pos = {a: i for i, a in enumerate(ids_used)}
    rows = [pos[i] for i in ids_basin]
    coords_b = np.asarray(inp.cc_c)[rows]
    counts_b = inp.counts.loc[ids_basin].to_numpy(float)
    lab_b, k_b, sil_b = cluster_labels(coords_b)
    gc_b = group_counts(counts_b, lab_b)

    sa, ga, da = fit(gc_a, f"A: as-run, {len(ids_used)} assemblages, k={k_a}")
    sb, gb, db = fit(gc_b, f"B: basin,  {len(ids_basin)} assemblages, k={k_b}")

    extras = sorted(set(ids_used) - set(ids_basin))
    L = [
        "# Does the Bayesian F_ST fit run on the basin? Measured.", "",
        "Produced by `analyses/47_basin_scope_check.py`. Seed "
        f"{SEED}, {CFG['draws']} draws x {CFG['chains']} chains, "
        "F ~ Uniform(0,1), clusters by the analysis-07 silhouette rule "
        "(k in 2..6, kmeans seed 7).", "",
        "## The two fits", "",
        "| | A: as analysis 43 runs today | B: restricted to the basin |",
        "|---|---|---|",
        f"| assemblages | **{len(ids_used)}** | **{len(ids_basin)}** |",
        f"| spatial clusters | {k_a} | {k_b} |",
        f"| BN F_ST median | **{sa['fst_median']:.4f}** | **{sb['fst_median']:.4f}** |",
        f"| BN 95% HDI | [{sa['fst_hdi95'][0]:.4f}, {sa['fst_hdi95'][1]:.4f}] "
        f"| [{sb['fst_hdi95'][0]:.4f}, {sb['fst_hdi95'][1]:.4f}] |",
        f"| Gini-Simpson median | {ga['gst_median']:.4f} | {gb['gst_median']:.4f} |",
        f"| Gini-Simpson 95% HDI | [{ga['gst_hdi95'][0]:.4f}, {ga['gst_hdi95'][1]:.4f}] "
        f"| [{gb['gst_hdi95'][0]:.4f}, {gb['gst_hdi95'][1]:.4f}] |",
        f"| plug-in F_ST | {ga['plugin_fst']:.4f} | {gb['plugin_fst']:.4f} |",
        f"| max R-hat | {sa['rhat']:.4f} | {sb['rhat']:.4f} |",
        f"| min ESS | {sa['ess']:.0f} | {sb['ess']:.0f} |",
        f"| divergences | {da} | {db} |", "",
        "## What is in A but not in the basin", "",
        f"{len(extras)} assemblages: " + ", ".join(f"`{e}`" for e in extras) + ".", "",
        "Silhouette scores on the basin coordinates, k = 2..6: "
        + ", ".join(f"k={k}: {v:.4f}" for k, v in sorted(sil_b.items()))
        + f". k={k_b} wins clearly, which is the \"three spatial clusters\" the "
        "manuscript describes. The k=5 in column A is the silhouette optimum for "
        "the wider curated set, not for the basin.", "",
        "## Reading", "",
        "`prepare_inputs()` returns the whole curated set with coordinates. "
        "`43_bayesian_fst.basin_group_counts()` groups it by `inp.cluster_of` "
        "without restricting to `data/processed/basin_members_curated.txt`, so "
        "despite its name and despite the report headed \"Observed St. Francis "
        "basin\", column A is a regional quantity. The excluded 26 include "
        "`Walls`, `Wall`, `Chuccalissa` and `Parchman`, which are the "
        "Mississippi-proximal sites that "
        "`analyses/16_basin_membership.py` names in its own docstring as the "
        "reason the earlier latitude cut was replaced by the hydrological rule.",
        "",
        "The Gini-Simpson readout, the quantity that matches the manuscript's "
        f"estimator, moves from {ga['gst_median']:.4f} to {gb['gst_median']:.4f}, "
        f"a {100 * (1 - gb['gst_median'] / ga['gst_median']):.0f} percent "
        "reduction. The BN parameter moves less but its interval widens as it "
        "should on 29 assemblages rather than 55. Both fits are numerically "
        "healthy, so this is a scope defect and not a sampling one.", "",
        "## Scope of the consequence", "",
        "- **No manuscript number changes.** Figures 4 and 5 reach the data "
        "through `make_figures._load_curated()`, which does apply the basin "
        "membership. Only analyses 43 and 44 use the unrestricted path.",
        "- `output/bayesian_fst.md` and the `STATUS.md` entry quoting "
        "\"basin BN F = 0.067\" and \"Gini-Simpson readout 0.0327\" attribute "
        "column A's numbers to the basin. Those are column B's questions with "
        "column A's answers.",
        "- `analyses/44_bayesian_fst_validation.py` calibrates coverage, SBC and "
        "small-sample behavior at column A's geometry (5 clusters over 55) while "
        "describing it as the basin fit. The validation is not wrong, but it "
        "validates a design the paper does not use.", "",
        "## Disposition", "",
        "Restrict `basin_group_counts()` to the canonical membership and "
        "re-select k on the basin coordinates, then regenerate 43 and 44. Add a "
        "test that asserts the fitted set equals "
        "`data/processed/basin_members_curated.txt`, since a name alone did not "
        "prevent this. Scheduled as the first task of item 1.", "",
        "## Why this matters beyond the numbers", "",
        "This is the same defect class the 2026-06-10 alignment run found and "
        "fixed, where five figures were computed on the whole curated set while "
        "the text reported the basin. It recurred in code written a month later. "
        "The Verification Regime preamble in `CLAUDE.md` cites that episode as "
        "this project's own precedent for gates being necessary and never "
        "sufficient; the recurrence is the evidence that the lesson needed the "
        "test it is now getting.", "",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
