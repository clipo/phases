"""49_write_stan_data.py - Stan data contracts for the spatial composition model.

Writes JSON to data/stan/ so the R/Stan layer and the Python analyses are
provably reading the same counts (program plan item 0 task 7).

Two kinds of file:

  basin_composition.json    the real St. Francis basin: 29 assemblages, 10
                            decorated classes, real counts, real pairwise
                            distances in km.
  sim_rho{R}.json           data simulated FROM the model at a known rho, for
                            the recovery run that gates the whole direction.

The simulator here is written in plain numpy from the model definition in
stan/spatial/composition_gp.stan and calls no Stan (Verification Regime rule 3:
expected values may not be produced by the implementation under test). The
Matern kernels are implemented independently from their mathematical form, not
translated from kernels.stanfunctions, so agreement between them is evidence
rather than tautology.

Usage: .venv/bin/python analyses/49_write_stan_data.py
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

OUT = ROOT / "data" / "stan"

# Prior configuration shared by every dataset written here, so a recovery fit
# and a real fit differ only in their counts.
PRIOR = dict(prior_rho_type=2,      # uniform on log rho: no scale preferred
             prior_rho_a=0.0, prior_rho_b=0.0,
             rho_lo=0.5, rho_hi=300.0,
             kernel=1,              # exponential (Matern 1/2)
             use_tau=1,
             prior_scale_sd=2.0,
             likelihood=1)

TRUE_RHO = (5.0, 10.0, 20.0, 40.0, 80.0)   # 80 probes above screen 48's ceiling
TRUE_SIGMA = 1.0
TRUE_TAU = 0.5
SEED = 20260831


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
    d = dict(N=int(y.shape[0]), K=int(y.shape[1]),
             y=[[int(v) for v in row] for row in y],
             D=[[float(v) for v in row] for row in D])
    d.update(PRIOR)
    d.update(over)
    return d


def main():
    a23 = importlib.import_module("23_phases_vs_spatial_drift")
    import make_figures as mf

    counts_df, coords_df = mf._load_curated()
    ids = [str(i) for i in counts_df.index]
    y = counts_df.to_numpy(int)
    D = a23.geo_km(np.asarray(coords_df, dtype=float))
    N, K = y.shape
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"basin: {N} assemblages x {K} classes, {y.sum()} sherds, "
          f"max distance {D.max():.1f} km")

    real = _payload(y, D)
    real["_ids"] = ids           # provenance; Stan ignores leading-underscore keys
    (OUT / "basin_composition.json").write_text(json.dumps(real))

    rng = np.random.default_rng(SEED)
    alpha = rng.normal(0.0, 1.0, K - 1)
    totals = y.sum(axis=1)
    truth = {}
    for rho in TRUE_RHO:
        ys = simulate_from_model(D, totals, K, rho, TRUE_SIGMA, TRUE_TAU,
                                 alpha, PRIOR["kernel"], np.random.default_rng(int(rho) + SEED))
        (OUT / f"sim_rho{int(rho)}.json").write_text(json.dumps(_payload(ys, D)))
        truth[str(rho)] = dict(rho=rho, sigma=TRUE_SIGMA, tau=TRUE_TAU,
                               spatial_share=TRUE_SIGMA ** 2 / (TRUE_SIGMA ** 2 + TRUE_TAU ** 2),
                               kernel=PRIOR["kernel"], sherds=int(ys.sum()))
        print(f"  sim rho={rho:5.1f} km -> {ys.sum()} sherds")

    (OUT / "sim_truth.json").write_text(json.dumps(
        dict(truth=truth, alpha=[float(a) for a in alpha], seed=SEED), indent=2))
    print(f"wrote {OUT}/basin_composition.json and {len(TRUE_RHO)} simulated sets")


if __name__ == "__main__":
    main()
