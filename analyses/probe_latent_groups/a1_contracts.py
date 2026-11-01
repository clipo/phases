"""Own-generator records for the A1 validation, and the real-data A1 contract (spec 4.1c, 6.6; R2-1).

Latent-group probe, pilot step 6. Exploratory design validation (rule 1).

Data are simulated from A1's own generative statement in numpy, sharing no
code with Stan (rule 3): per class a GP field with the exponential kernel at
the pinned sigma and rho plus an assemblage term with SD tau, centered per
class across sites; an intercept equal to the A0 refit's posterior mean of
the identified sum alpha_c + mean_i f_ic (a0_pinned_nuisances.json); group
shifts delta = sigma_g times the sum-to-zero projection of iid standard
normals (the same distribution as Stan's sigma_g Q z); membership z the
Phillips partition (12/11/5); counts multinomial at the real sherd totals.

Records: "zero" (sigma_g = 0; used for the code-equivalence test with the
group scale fixed at zero and for the null calibration), and five records at
each of sigma_g = 0.25 and 0.5. Each file carries "_truth" (sigma_g, delta,
z, Delta at the true group weights, the realized between-true-group plug-in
F_ST of the counts) and the A1 data fields (G = 3, prior_group_sd = 2,
Gabriel edges, group scale free). The real-data A1 contract is the A0
contract plus the same A1 fields.

Seeds: record r of level L uses SeedSequence(2026101700 + 100 L + r).
Usage: .venv/bin/python analyses/probe_latent_groups/a1_contracts.py
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.dataio.basin import data_fingerprint, geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import readouts as ro
from mls_emergence.signatures.variance import cultural_fst

OUT = ROOT / "data" / "stan"
PIN = ROOT / "output" / "probe_latent_groups" / "a0_pinned_nuisances.json"
G = 3
PRIOR_GROUP_SD = 2.0
LEVELS = {0: [0.0], 1: [0.25] * 5, 2: [0.5] * 5, 3: [1.0] * 2, 4: [2.0] * 2}   # levels 3 and 4 added 2026-10-17: large-signal software check
# Level 5 (2026-10-17, mode-check follow-up): the same sizes 12/11/5 drawn as a SCATTERED partition (uniform over
# labelings), at sigma_g 0.5 and 1.0, to test whether A1 recovers groups the smooth field cannot mimic.
SCATTERED = {5: [0.5, 1.0]}
SEED_BASE = 2026101700


def a1_fields(edges):
    return {"G": G, "prior_group_sd": PRIOR_GROUP_SD, "fix_group_scale": 0, "group_scale_fixed": 0.0,
            "M": int(len(edges)), "edges": [[int(i) + 1, int(j) + 1] for i, j in edges]}


def simulate(D, totals, alpha_sum, sigma, tau, rho, sigma_g, z, rng):
    N = len(D); K = len(alpha_sum) + 1; C = K - 1
    cov = sigma ** 2 * np.exp(-D / rho) + 1e-8 * np.eye(N)
    L = np.linalg.cholesky(cov)
    f = L @ rng.standard_normal((N, C)) + tau * rng.standard_normal((N, C))
    f -= f.mean(axis=0, keepdims=True)                 # the constant mode is carried by the pinned intercept
    u = rng.standard_normal((G, C)); delta = sigma_g * (u - u.mean(axis=0, keepdims=True))
    eta = np.zeros((N, K)); eta[:, 1:] = alpha_sum[None, :] + f + delta[z]
    eta -= eta.max(axis=1, keepdims=True); p = np.exp(eta); p /= p.sum(axis=1, keepdims=True)
    y = np.array([rng.multinomial(int(totals[i]), p[i]) for i in range(N)])
    return y, f, delta


def main() -> int:
    counts, coords = load_curated(); fp = data_fingerprint(counts, coords)
    y_real = counts.to_numpy(int); totals = y_real.sum(1); D = geo_km(coords.to_numpy(float))
    pin = json.loads(PIN.read_text())
    m36 = importlib.import_module("36_canonical_phase_map")
    labels, _ = m36.assign_primary_phases([str(i) for i in counts.index], coords.to_numpy(float))
    names = sorted(set(labels)); z = np.array([names.index(v) for v in labels])
    sizes = np.bincount(z, minlength=G); pi_true = sizes / sizes.sum()
    edges = gen.gabriel_edges(D)
    base = json.loads((OUT / "basin_composition.json").read_text())
    a0_fields = {k: v for k, v in base.items() if not k.startswith("_") and k not in ("y",)}
    alpha_sum = np.asarray(pin["alpha_plus_field_mean"], float)
    sigma, tau, rho = float(pin["sigma"]), float(pin["tau"]), float(pin["rho"])
    OUT.mkdir(parents=True, exist_ok=True)
    truth_all = {}
    for level, sgs in LEVELS.items():
        for r, sg in enumerate(sgs):
            rng = np.random.default_rng(np.random.SeedSequence(SEED_BASE + 100 * level + r))
            y, _f, delta = simulate(D, totals, alpha_sum, sigma, tau, rho, sg, z, rng)
            tag = "zero" if sg == 0 else f"sg{sg:g}_r{r + 1}"
            truth = {"sigma_g": sg, "delta": delta.tolist(), "z": z.tolist(), "z_names": names, "pi_true": pi_true.tolist(),
                     "Delta_true": ro.delta_contrast(delta, pi_true), "sigma": sigma, "tau": tau, "rho": rho,
                     "alpha_plus_field_mean": alpha_sum.tolist(),
                     "realized_fst_between_true_groups": float(cultural_fst(np.array([y[z == g].sum(0) for g in range(G)]))),
                     "field_centered_per_class": True, "generator": "analyses/probe_latent_groups/a1_contracts.simulate",
                     "seed": SEED_BASE + 100 * level + r}
            payload = dict(a0_fields); payload.update({"y": y.tolist()}); payload.update(a1_fields(edges))
            payload["_truth"] = truth; payload["_provenance"] = {"derived_from": fp, "note": "simulated from A1 at the real geometry and totals"}
            (OUT / f"a1_{tag}.json").write_text(json.dumps(payload)); truth_all[tag] = {k: v for k, v in truth.items() if k not in ("delta",)}
            print(f"  {tag}: sigma_g {sg:g}, Delta_true {truth['Delta_true']:.3f}, realized between-group F_ST {truth['realized_fst_between_true_groups']:.4f}, sherds {int(y.sum())}")
    for level, sgs in SCATTERED.items():
        for r, sg in enumerate(sgs):
            rng = np.random.default_rng(np.random.SeedSequence(SEED_BASE + 100 * level + r))
            z_sc = rng.permutation(z)                       # same sizes, uniformly scattered over sites
            y, _f, delta = simulate(D, totals, alpha_sum, sigma, tau, rho, sg, z_sc, rng)
            tag = f"scattered_sg{sg:g}"
            truth = {"sigma_g": sg, "delta": delta.tolist(), "z": z_sc.tolist(), "z_names": names, "pi_true": pi_true.tolist(),
                     "Delta_true": ro.delta_contrast(delta, pi_true), "sigma": sigma, "tau": tau, "rho": rho,
                     "alpha_plus_field_mean": alpha_sum.tolist(), "partition": "scattered (uniform permutation of the Phillips sizes)",
                     "coherence_true_hard": ro.coherence(ro.hard_comembership(z_sc), edges),
                     "realized_fst_between_true_groups": float(cultural_fst(np.array([y[z_sc == g].sum(0) for g in range(G)]))),
                     "field_centered_per_class": True, "generator": "analyses/probe_latent_groups/a1_contracts.simulate", "seed": SEED_BASE + 100 * level + r}
            payload = dict(a0_fields); payload.update({"y": y.tolist()}); payload.update(a1_fields(edges))
            payload["_truth"] = truth; payload["_provenance"] = {"derived_from": fp, "note": "simulated from A1 at the real geometry and totals; scattered groups"}
            (OUT / f"a1_{tag}.json").write_text(json.dumps(payload)); truth_all[tag] = {k: v for k, v in truth.items() if k != "delta"}
            print(f"  {tag}: sigma_g {sg:g}, Delta_true {truth['Delta_true']:.3f}, true hard coherence {truth['coherence_true_hard']:.3f}, realized between-group F_ST {truth['realized_fst_between_true_groups']:.4f}")
    real = dict(base); real.update(a1_fields(edges)); real["_edges_note"] = "Gabriel graph on geo_km distances (D14)"
    (OUT / "basin_composition_a1.json").write_text(json.dumps(real))
    (OUT / "a1_truth.json").write_text(json.dumps(truth_all, indent=1))
    print(f"wrote {len(truth_all)} own-generator records and basin_composition_a1.json to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
