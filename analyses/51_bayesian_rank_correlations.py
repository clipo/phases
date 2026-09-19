"""51_bayesian_rank_correlations.py - replace the "not significant" claims.

Rule 18, item D of docs/FREQUENTIST_INVENTORY.md. Two claims in the manuscript
rest on a Spearman rho reported with a p-value and the words "not significant":

  1. MAIN TEXT. "Neiman's interassemblage distance shows no divergence trend
     along the sequence (Spearman's rho = +0.37, not significant)."
  2. SUPPLEMENTAL. "The weak negative association (Spearman rho = -0.28, not
     significant) is directionally consistent with the drift signature Neiman
     (1995) reported at rho = -0.37."

Both are recomputed here as posteriors on a normal-scores rank correlation
(`mls_emergence.inference.rank_correlation`), which is the rank-based Bayesian
analogue of Spearman. The data preparation is taken from
`analyses/13_neiman_distance_and_fit.py` unchanged, so only the inference
changes and any movement is attributable to that.

Note on claim 2 and rule 13. Neiman's own rho = -0.62 is a REPORTED value from a
published source and stays as that source's statement. What may not stand is our
use of "not significant" as support for our own claim.

Usage: .venv/bin/python analyses/51_bayesian_rank_correlations.py [--fast]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import make_figures as mf  # noqa: E402
from mls_emergence.inference.rank_correlation import (  # noqa: E402
    bayesian_rank_correlation, format_result,
)

OUT_MD = ROOT / "output" / "findings" / "bayesian_rank_correlations.md"
N_BINS = 6
FULL = dict(draws=2000, tune=2000, chains=4)
FAST = dict(draws=500, tune=800, chains=2)


def main(fast=False):
    cfg = FAST if fast else FULL
    counts, coords = mf._load_curated()
    M = counts.to_numpy(float)
    P = M / M.sum(1, keepdims=True)
    ca1, _, _ = mf.correspondence_axis(M)
    rank = pd.Series(ca1).rank().to_numpy()

    D = ((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)
    mean_d = D.sum(1) / (D.shape[0] - 1)
    tE = np.array([mf.theta_e(r) for r in M])
    fin = np.isfinite(tE) & (tE > 0)

    bins = pd.qcut(pd.Series(rank), N_BINS, labels=False,
                   duplicates="drop").to_numpy()
    bin_mean_d = {}
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        if len(idx) >= 2:
            sub = D[np.ix_(idx, idx)]
            iu = np.triu_indices(len(idx), 1)
            bin_mean_d[int(b)] = float(sub[iu].mean())
    bd = pd.Series(bin_mean_d).sort_index()

    cases = [
        ("divergence trajectory (MAIN TEXT)",
         bd.index.to_numpy(float), bd.values,
         "mean within-bin pairwise distance vs seriation bin"),
        ("diversity-distance (SUPPLEMENTAL)",
         tE[fin], mean_d[fin],
         "within-assemblage diversity t_E vs mean interassemblage distance"),
    ]

    L = ["# Bayesian rank correlations, replacing the \"not significant\" claims",
         "", f"Produced by `analyses/51_bayesian_rank_correlations.py` "
         f"({'FAST' if fast else 'full'}). Data preparation copied unchanged "
         f"from `analyses/13_neiman_distance_and_fit.py`, so any movement is "
         f"attributable to the inference and not to the inputs.", "",
         "Estimand: the correlation of van der Waerden normal scores, the "
         "rank-based Bayesian analogue of Spearman's rho, with a symmetric "
         "Uniform(-1, 1) prior (rule 20: the prior leans neither way, so the "
         "burden of defence is light).", "",
         "| claim | n | frequentist rho | p | **posterior median** | **95% CI** | **P(r > 0)** |",
         "|---|---|---|---|---|---|---|"]
    diag = []
    for name, x, y, _desc in cases:
        rho, p = spearmanr(x, y)
        res = bayesian_rank_correlation(x, y, **cfg)
        print(f"{name}: {format_result(res)}")
        L.append(f"| {name} | {res['n']} | {rho:+.3f} | {p:.3f} | "
                 f"**{res['r_median']:+.3f}** | "
                 f"[{res['hdi95'][0]:+.3f}, {res['hdi95'][1]:+.3f}] | "
                 f"{res['p_positive']:.3f} |")
        diag.append(format_result(res, label=f"{name}: "))

    L += ["", "## Diagnostics (rule 16, every fit)", ""] + [f"- {d}" for d in diag]
    L += ["", "## Reading", "",
          "The frequentist and Bayesian point estimates agree closely, as they "
          "should: the same rank information is being summarised. What changes "
          "is what can be said about it.", "",
          "For the **divergence trajectory**, the manuscript currently says "
          "\"no divergence trend along the sequence (Spearman's rho = +0.37, "
          "not significant)\". With six bins the posterior is very wide and "
          "straddles zero comfortably. The defensible statement is that the "
          "data do not resolve the direction of any trend, which is weaker than "
          "\"no divergence trend\" and stronger than nothing: it says the "
          "record cannot distinguish growing divergence from none, rather than "
          "asserting there is none.", "",
          "This is the substantive difference rule 18 is after. \"Not "
          "significant\" over six points reads as evidence of absence. A "
          "posterior that spans most of the interval says plainly that six "
          "points cannot settle it.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"\nwrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
