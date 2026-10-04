"""57_fst_prior_predictive.py - what does the F_ST prior imply about data?

Closes F2 in docs/CODE_REVIEW_2026-08-31.md and satisfies rule 20(a) for the
one prior in this project that has never had it: the Balding-Nichols
`F ~ Uniform(0, 1)` used by analyses 43, 44 and 50.

THE ARGUMENT ON RECORD IS ABOUT THE WRONG SCALE. The docstring of
`src/mls_emergence/inference/bayesian_fst.py` defends the flat prior on the
grounds that it "does not push the estimate toward the drift level the paper
reports". That is a statement about the parameter. Rule 20(a) asks a different
question: what DATA does this prior generate, and does that data already look
like what was observed? A prior that is flat on a parameter can be violently
informative about the observable.

THE FAMILY CONSTRAINT IS NOT A MATTER OF TASTE. Any Beta(a, b) with a > 1 has
density exactly zero at F = 0, so it asserts a priori that panmixia is
impossible. For an analysis whose comparison is against panmixia that is
disqualifying rather than merely unattractive, so the family must be Beta(1, b),
which contains Uniform(0, 1) as b = 1. Beta(2, 20) is included below only to
show what the excluded shape does.

DIRECTION (rule 20, first step). `F ~ Uniform(0, 1)` leans toward LARGE F_ST.
For the posterior that is conservative, since the paper argues F_ST is small and
a prior pulling the other way cannot manufacture the conclusion. For a Bayes
factor it is the reverse: diffuseness over values the data never support
penalises the structured model, which is sympathetic to the paper's conclusion.
Both directions are stated because a prior can lean two ways at once and rule 20
requires the direction be named before the burden can be judged.

Usage: .venv/bin/python analyses/57_fst_prior_predictive.py [--fast]
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
from figstyle import OI_BLUE, OI_VERMIL, OI_GREEN, OI_ORANGE, save  # noqa: E402

OUT_MD = ROOT / "output" / "findings" / "fst_prior_justification.md"
PRIORS = [("Uniform(0,1)", 1.0, 1.0), ("Beta(1,3)", 1.0, 3.0),
          ("Beta(1,10)", 1.0, 10.0), ("Beta(2,20)", 2.0, 20.0)]
BIG = 0.30      # "strong structure" threshold, as in ../mataa's write-up


def simulate_gst(a, b, sizes, K, n_draw, rng):
    """Induced Gini-Simpson F_ST under F ~ Beta(a, b), at a real design."""
    out = np.empty(n_draw)
    G = len(sizes)
    for i in range(n_draw):
        F = rng.beta(a, b)
        F = min(max(F, 1e-9), 1 - 1e-9)
        conc = (1.0 - F) / F
        pi = rng.dirichlet(np.ones(K))
        gc = np.empty((G, K))
        for g in range(G):
            p = rng.dirichlet(conc * pi)
            gc[g] = rng.multinomial(int(sizes[g]), p)
        w = gc.sum(1) / gc.sum()
        hs = float(np.sum(w * (1.0 - np.sum((gc / gc.sum(1, keepdims=True)) ** 2, axis=1))))
        pool = gc.sum(0) / gc.sum()
        ht = float(1.0 - np.sum(pool ** 2))
        out[i] = (ht - hs) / ht if ht > 0 else np.nan
    return out


def main(fast=False):
    n_draw = 400 if fast else 4000
    a07 = importlib.import_module("07_refined_empirical")
    a43 = importlib.import_module("43_bayesian_fst")
    inp = a07.prepare_inputs()
    gc, sizes = a43.basin_group_counts(inp)          # basin, 3 clusters
    K = gc.shape[1]
    obs = None
    from mls_emergence.signatures.variance import cultural_fst
    obs = cultural_fst(gc)
    print(f"design: {gc.shape[0]} clusters, {K} classes, sizes {sizes.astype(int).tolist()}")
    print(f"observed Gini-Simpson F_ST = {obs:.4f}")

    rows = {}
    for name, a, b in PRIORS:
        s = simulate_gst(a, b, sizes, K, n_draw, np.random.default_rng(7))
        s = s[np.isfinite(s)]
        rows[name] = s
        print(f"  {name:<13} median {np.median(s):.3f}  "
              f"95% [{np.percentile(s,2.5):.3f}, {np.percentile(s,97.5):.3f}]  "
              f"P(F_ST>{BIG}) = {float((s>BIG).mean()):.3f}")

    u = rows["Uniform(0,1)"]
    b10 = rows["Beta(1,10)"]
    L = ["# What does the F_ST prior imply about data?", "",
         f"Produced by `analyses/57_fst_prior_predictive.py` "
         f"({'FAST' if fast else 'full'}, {n_draw} draws per prior). Closes F2 "
         f"and supplies rule 20(a) for the Balding-Nichols prior.", "",
         f"Design: the St. Francis basin as analysis 43 fits it, "
         f"{gc.shape[0]} spatial clusters over {K} decorated classes at the real "
         f"per-cluster sherd totals {sizes.astype(int).tolist()}. "
         f"**Observed Gini-Simpson F_ST = {obs:.4f}.**", "",
         "## Induced prior predictive distribution of F_ST", "",
         "| prior | median | 95% predictive interval | P(F_ST > 0.30) |",
         "|---|---|---|---|"]
    for name, _, _ in PRIORS:
        s = rows[name]
        L.append(f"| {name} | {np.median(s):.3f} | [{np.percentile(s,2.5):.3f}, "
                 f"{np.percentile(s,97.5):.3f}] | **{float((s>BIG).mean()):.3f}** |")
    L += ["", "## Reading", "",
          f"Before seeing any data, `F ~ Uniform(0, 1)` expects a median "
          f"Gini-Simpson F_ST of **{np.median(u):.3f}** at this design and puts "
          f"**{100*float((u>BIG).mean()):.0f} percent** of its mass above 0.30. "
          f"The observed value is {obs:.4f}. Nobody holds the belief that "
          f"decorated-ceramic assemblages a few tens of kilometres apart in one "
          f"drainage are that strongly differentiated, so on the scale that "
          f"matters the flat prior is not uninformative; it is strongly and "
          f"wrongly informative.", "",
          f"`Beta(1, 10)` expects a median of {np.median(b10):.3f} with "
          f"{100*float((b10>BIG).mean()):.0f} percent above 0.30, which is a "
          f"defensible prior belief for this setting.", "",
          "## The family constraint", "",
          "Any `Beta(a, b)` with `a > 1` has density exactly zero at F = 0 and "
          "so asserts a priori that panmixia is impossible. Because the whole "
          "comparison in analyses 43 and 44 is against panmixia, that is "
          "disqualifying rather than merely unattractive. Beta(2, 20) is "
          "included in the table only to show what the excluded shape does; the "
          "admissible family is `Beta(1, b)`, of which Uniform(0, 1) is the "
          "b = 1 member.", "",
          "## Direction, and why the burden is asymmetric (rule 20)", "",
          "Uniform(0, 1) leans toward LARGE F_ST. For the posterior that is "
          "**conservative**: the paper argues F_ST is small, and a prior pulling "
          "the other way cannot manufacture that conclusion, which is why the "
          "reported posteriors are not in doubt. For a Bayes factor the lean "
          "reverses and becomes **sympathetic**: spreading prior mass over "
          "values the data never support penalises the structured model against "
          "panmixia, which is the direction of the paper's own reading. A prior "
          "can lean two ways at once, and rule 20 requires both be named.", "",
          "## Disposition", "",
          "The prior predictive disqualifies Uniform(0, 1) as a description of "
          "prior belief, and `../mataa` reached the same verdict independently "
          "on three of its own designs. Because the lean is conservative for the "
          "posterior, no reported posterior is overturned by this; what changes "
          "is that the prior can no longer be defended as uninformative, and any "
          "Bayes factor computed under it is reading the prior as much as the "
          "data. Adopting Beta(1, 10) as primary with Uniform reported alongside "
          "is the change this finding supports.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for (name, _, _), c in zip(PRIORS, (OI_BLUE, OI_ORANGE, OI_GREEN, OI_VERMIL)):
        ax.hist(rows[name], bins=60, range=(0, 1), density=True, histtype="step",
                color=c, label=name, lw=1.3)
    ax.axvline(obs, color="0.2", ls="--", lw=1.2)
    ax.text(obs + 0.015, ax.get_ylim()[1] * 0.9, f"observed {obs:.3f}", fontsize=6.5)
    ax.set_xlabel("prior predictive Gini-Simpson $F_{ST}$")
    ax.set_ylabel("density"); ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout(); save(fig, "fig_fst_prior_predictive")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
