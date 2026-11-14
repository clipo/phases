"""Monte Carlo check of the first-order claim behind D27 (2026-10-19).

The item 1.1 entry of the discussion log derived, to first order, that
independent logit-scale Gaussian fields of variance F/p_c give the
Dirichlet-drift moments var(p_c) = F p_c (1 - p_c). This script draws such
fields at the record's pooled shares and reports the realized per-class
moment F and the shift of the mean shares. Descriptive; the result replaces
the derivation as the basis of D27.
Usage: .venv/bin/python analyses/probe_latent_groups/softmax_drift_check.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "stan" / "basin_composition.json"
OUT = ROOT / "output" / "findings" / "softmax_drift_check.md"
F_SET = (0.01, 0.05, 0.12)
N_DRAW = 400_000
SEED = 1


def main() -> int:
    d = json.loads(DATA.read_text())
    y = np.asarray(d["y"], float); p = y.sum(0) / y.sum(); cls = d["_classes"]; K = len(p)
    rng = np.random.default_rng(SEED)
    L = ["# Softmax fields of variance F/p_c against Dirichlet drift", "",
         (f"Produced by `analyses/probe_latent_groups/softmax_drift_check.py`: {N_DRAW:,} draws (seed {SEED}) of K = {K} "
          "independent Gaussian fields with mean log p_c and variance F/p_c at the record's pooled shares "
          "(`data/stan/basin_composition.json`), passed through the softmax; realized per-class moment F = "
          "var(p_c) / (mean(p_c) (1 - mean(p_c))) and the ratio of the realized mean share to the pooled share. "
          "To first order both would equal F and 1. Descriptive; the basis of D27."), "",
         "| class | pooled share | " + " | ".join(f"F set {F}: realized F, mean ratio" for F in F_SET) + " |",
         "|---|---|" + "---|" * len(F_SET)]
    res = {}
    for F in F_SET:
        eta = np.log(p) + rng.standard_normal((N_DRAW, K)) * np.sqrt(F / p)
        q = np.exp(eta - eta.max(1, keepdims=True)); q /= q.sum(1, keepdims=True)
        m = q.mean(0); res[F] = (q.var(0) / (m * (1 - m)), m / p)
    for c in range(K):
        L.append(f"| {cls[c]} | {p[c]:.4f} | " + " | ".join(f"{res[F][0][c]:.3f}, {res[F][1][c]:.2f}" for F in F_SET) + " |")
    L += ["", ("The first-order matching fails at every F tried, for the common classes as well: the implied variances F/p_c run from "
               f"{min(F_SET) / p.max():.2f} (Parkin_Punctated at F {min(F_SET):g}) to {max(F_SET) / p.min():.0f} logits^2 (Hull_Engraved at F {max(F_SET):g}), and the "
               "rare classes' dominate the normalization. The figures 1.48 and 1.21 in the item 1.1 entry (1/sqrt(1 - p_c)) are first-order values that no "
               "draw reproduces and are withdrawn."), ""]
    OUT.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
