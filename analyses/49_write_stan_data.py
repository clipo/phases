"""49_write_stan_data.py - Stan data contracts for the spatial composition model.

Writes JSON to data/stan/ so the R/Stan layer and the Python analyses are
provably reading the same counts (program plan item 0 task 7).

Two kinds of file:

  basin_composition.json    the real St. Francis basin: 28 assemblages, 10
                            decorated classes, real counts, real pairwise
                            distances in km, with the data fingerprint
                            (hashes and geometry) under "_provenance".
  a0_recovery_share{S}_rho{R}.json
                            data simulated FROM the model for the rule-20(c)
                            recovery fits of the latent-group probe (spec 4.1):
                            spatial share 0.2 and 0.9 crossed with rho 5 and
                            60 km, one record each, truth under "_truth".
                            (Until 2026-10-16 this script wrote sim_rho{R}.json
                            for the superseded 01_recovery.R; nothing runnable
                            reads those.)

The simulator here is written in plain numpy from the model definition in
stan/spatial/composition_gp.stan and calls no Stan (Verification Regime rule 3:
expected values may not be produced by the implementation under test). The
Matern kernels are implemented independently from their mathematical form, not
translated from kernels.stanfunctions, so agreement between them is evidence
rather than tautology.

Usage: .venv/bin/python analyses/49_write_stan_data.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.dataio.basin import data_fingerprint, geo_km, load_curated

OUT = ROOT / "data" / "stan"

# Prior configuration shared by every dataset written here, so a recovery fit
# and a real fit differ only in their counts.
PRIOR = {"prior_rho_type": 2,      # uniform on log rho: no scale preferred
             "prior_rho_a": 0.0, "prior_rho_b": 0.0,
             "rho_lo": 0.5, "rho_hi": 300.0,
             "kernel": 1,              # exponential (Matern 1/2)
             "use_tau": 1,
             "prior_scale_sd": 2.0,
             "prior_alpha_sd": 10.0,    # intercept prior SD (D25, 2026-10-18: was a hard-coded 2, too narrow for rare-class log-odds of -4 to -6); mean written per payload
             "likelihood": 1}

# Rule-20(c) recovery grid for A0 (probe spec 4.1): the low share is the value
# that would embarrass the paper's headline claim. Total latent variance is
# fixed at 1 logit^2, so sigma = sqrt(share) and tau = sqrt(1 - share).
RECOVERY_SHARES = (0.2, 0.9)
RECOVERY_RHOS = (5.0, 60.0)
SEED = 20261016


def matern(D, sigma, rho, kernel):
    """Matern covariance, written from the definition. Independent of Stan."""
    D = np.asarray(D, float)
    if kernel == 1:
        C = np.exp(-D / rho)
    elif kernel == 2:
        s = np.sqrt(3.0) * D / rho
        C = (1.0 + s) * np.exp(-s)
    elif kernel == 3:
        s = np.sqrt(5.0) * D / rho
        C = (1.0 + s + s ** 2 / 3.0) * np.exp(-s)
    else:
        raise ValueError(f"kernel must be 1, 2 or 3; got {kernel}")
    return sigma ** 2 * C


def simulate_from_model(D, totals, K, rho, sigma, tau, alpha, kernel, rng):
    """Draw counts from composition_gp.stan's generative statement."""
    N = D.shape[0]
    C = K - 1
    cov = matern(D, sigma, rho, kernel) + 1e-8 * np.eye(N)
    L = np.linalg.cholesky(cov)
    f = L @ rng.standard_normal((N, C))
    u = tau * rng.standard_normal((N, C))
    eta = np.zeros((N, K))
    eta[:, 1:] = alpha[None, :] + u + f
    eta -= eta.max(axis=1, keepdims=True)
    p = np.exp(eta)
    p /= p.sum(axis=1, keepdims=True)
    return np.array([rng.multinomial(int(totals[i]), p[i]) for i in range(N)])


def _payload(y, D, **over):
    d = {"N": int(y.shape[0]), "K": int(y.shape[1]),
             "y": [[int(v) for v in row] for row in y],
             "D": [[float(v) for v in row] for row in D]}
    d.update(PRIOR)
    d["prior_alpha_mean"] = [0.0] * (int(y.shape[1]) - 1)   # default centre; D25 variants override it
    d.update(over)
    return d


def main():
    counts_df, coords_df = load_curated()
    ids = [str(i) for i in counts_df.index]
    y = counts_df.to_numpy(int)
    D = geo_km(np.asarray(coords_df, dtype=float))
    N, K = y.shape
    OUT.mkdir(parents=True, exist_ok=True)
    fp = data_fingerprint(counts_df, coords_df)

    print(f"basin: {N} assemblages x {K} classes, {y.sum()} sherds, "
          f"max distance {fp['max_pairwise_km']:.1f} km, counts sha256 {fp['sha256_counts'][:12]}")

    real = _payload(y, D)
    real["_ids"] = ids           # provenance; Stan ignores leading-underscore keys
    real["_provenance"] = fp
    real["_classes"] = list(counts_df.columns)
    (OUT / "basin_composition.json").write_text(json.dumps(real))

    rng = np.random.default_rng(SEED)
    alpha = rng.normal(0.0, 1.0, K - 1)
    totals = y.sum(axis=1)
    truth = {}
    for share in RECOVERY_SHARES:
        for rho in RECOVERY_RHOS:
            sigma, tau = float(np.sqrt(share)), float(np.sqrt(1.0 - share))
            tag = f"share{share:g}_rho{rho:g}"
            ys = simulate_from_model(D, totals, K, rho, sigma, tau, alpha, PRIOR["kernel"],
                                     np.random.default_rng(SEED + int(1000 * share) + int(rho)))
            t = {"spatial_share": share, "rho": rho, "sigma": sigma, "tau": tau,
                 "alpha": [float(a) for a in alpha], "kernel": PRIOR["kernel"],
                 "sherds": int(ys.sum()), "seed": SEED, "generator": "analyses/49_write_stan_data.simulate_from_model"}
            payload = _payload(ys, D)
            payload["_truth"] = t
            payload["_provenance"] = {"derived_from": fp, "note": "simulated counts at the real geometry and totals"}
            (OUT / f"a0_recovery_{tag}.json").write_text(json.dumps(payload))
            truth[tag] = t
            print(f"  recovery {tag}: sigma {sigma:.3f} tau {tau:.3f} -> {ys.sum()} sherds")

    (OUT / "a0_recovery_truth.json").write_text(json.dumps({"truth": truth, "seed": SEED}, indent=2))
    print(f"wrote {OUT}/basin_composition.json and {len(truth)} recovery sets")


if __name__ == "__main__":
    main()
