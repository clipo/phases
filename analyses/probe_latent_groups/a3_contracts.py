"""Records for the A3 identification test (design section 4; D23). Step 6b.

Exploratory design validation (rule 1). Three kinds of record, every one with
the A3 data fields (P = 1, the true or supplied placement; prior_group_sd 2;
omega free) and "_truth":

1. Own-generator from A3 (numpy, no Stan): nuisances pinned to the A0
   refit's identified quantities as for A1 (field centered per class), a
   region-level shift per class with SD omega shared by every site in a
   region, omega in {0, 0.25, 0.5, 1.0}, at two placements (G1, the k-means-2
   split; L, the candidate line that best agrees with Parkin against the
   rest, used only to choose a test placement), three records each.
1b. Step against slope: omega = 0 with the field's range at 30 and 60 km
   (pinned sigma and tau), three each; fitted with the G1 placement supplied.
2. Copying-simulator records at the primary operating point (D22) under T1
   spans at the measured time scales: S0 (no barrier), S1 (Phillips's areas,
   leak 0.1), S2 (leak 0.03), three each; the supplied placement for S1 and
   S2 is Phillips's three areas (their truth), for S0 the G1 split.

Seeds: SeedSequence(2026101720 + 1000 * kind + 100 * level + r).
Usage: .venv/bin/python analyses/probe_latent_groups/a3_contracts.py [--alpha-sd 2]
With --alpha-sd 2 the records are simulated at the Normal(0, 2) basin fit's
values and written to data/stan/alphaN0-2/ with that prior (the
prior-sensitivity arm). The simulator records (kind 2) do not depend on the
prior and are the same in both arms apart from the prior field.
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.dataio.basin import data_fingerprint, geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import placement as pl
from mls_emergence.probe import spans
from mls_emergence.signatures.variance import cultural_fst
from mls_emergence.transmission.spatial import drift_record

BASE = ROOT / "data" / "stan"
PINDIR = ROOT / "output" / "probe_latent_groups"
# The Normal(0, 10) arm's records directory and pinned values. main() uses the arm's own pair
# (intercept_arm); these module-level names are what scripts/build_manifest.py resolves.
OUT = ROOT / "data" / "stan"
PIN = ROOT / "output" / "probe_latent_groups" / "a0_pinned_nuisances.json"
POINTS = ROOT / "analyses" / "probe_latent_groups" / "operating_points.json"
TS = ROOT / "output" / "probe_latent_groups" / "time_scales_primary.json"
PRIOR_GROUP_SD = 2.0
OMEGAS = (0.0, 0.25, 0.5, 1.0)
SLOPE_RHOS = (30.0, 60.0)
LEAKS = {"S0": 1.0, "S1": 0.1, "S2": 0.03}
N_REC = 3
SEED_BASE = 2026101720


def intercept_arm(alpha_sd: float) -> dict:
    """Records directory, pinned-values file and prior SD for one intercept-prior arm.

    10 is the paper's Normal(0, 10) prior (the default paths). 2 is the
    prior-sensitivity arm: records simulated at the Normal(0, 2) basin fit's
    values, written to data/stan/alphaN0-2/ with prior_alpha_sd 2. Any other
    value raises.
    """
    if alpha_sd == 10:
        return {"sd": 10.0, "out": BASE, "pin": PINDIR / "a0_pinned_nuisances.json"}
    if alpha_sd == 2:
        return {"sd": 2.0, "out": BASE / "alphaN0-2", "pin": PINDIR / "a0_pinned_nuisances_alphaN0-2.json"}
    raise ValueError("alpha_sd takes 10 (the paper's prior) or 2 (the prior-sensitivity arm)")


def arm_fields(base: dict, pin: dict, arm: dict) -> dict:
    """The basin contract's fields with the arm's intercept prior; refuses pinned values from the other prior's fit."""
    if float(base["prior_alpha_sd"]) != 10.0:
        raise ValueError(f"basin_composition.json carries prior_alpha_sd {base['prior_alpha_sd']}, expected the paper's 10")
    if float(pin.get("prior_alpha_sd", arm["sd"])) != arm["sd"]:
        raise ValueError(f"the pinned values were written for Normal(0, {pin['prior_alpha_sd']:g}), not Normal(0, {arm['sd']:g})")
    if arm["sd"] != 10.0 and "prior_alpha_sd" not in pin:
        raise ValueError("the pinned values do not state their prior; rerun a0_pinned_nuisances.R --alpha-sd=2")
    f = {k: v for k, v in base.items() if not k.startswith("_") and k != "y"}
    f["prior_alpha_sd"] = arm["sd"]
    return f


def a3_fields(regions):
    return {"P": 1, "region": [[int(r) for r in regions]], "prior_group_sd": PRIOR_GROUP_SD, "fix_omega": 0}


def simulate_a3(D, totals, alpha_sum, sigma, tau, rho, omega, regions, rng):
    N = len(D); K = len(alpha_sum) + 1; C = K - 1
    L = np.linalg.cholesky(sigma ** 2 * np.exp(-D / rho) + 1e-8 * np.eye(N))
    f = L @ rng.standard_normal((N, C)) + tau * rng.standard_normal((N, C))
    f -= f.mean(axis=0, keepdims=True)
    reg = np.asarray(regions); labs = sorted(set(reg.tolist()))
    shift = omega * rng.standard_normal((len(labs), C))
    eta = np.zeros((N, K)); eta[:, 1:] = alpha_sum[None, :] + f + shift[[labs.index(r) for r in reg]]
    eta -= eta.max(axis=1, keepdims=True); p = np.exp(eta); p /= p.sum(axis=1, keepdims=True)
    y = np.array([rng.multinomial(int(totals[i]), p[i]) for i in range(N)])
    return y, shift


def between_fst(y, reg):
    reg = np.asarray(reg)
    return float(cultural_fst(np.array([y[reg == g].sum(0) for g in np.unique(reg)])))


def choose_line(labels, cands) -> dict:
    """Placement L: the candidate line that best agrees with Phillips's Parkin-versus-rest labels."""
    return pl.best_line_for(np.where(np.asarray(labels) == "Parkin", 1, 2), cands)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--alpha-sd", type=float, default=10.0, help="intercept prior SD of the arm: 10 (default, the paper's prior) or 2")
    arm = intercept_arm(ap.parse_args(argv).alpha_sd); OUT, PIN = arm["out"], arm["pin"]; OUT.mkdir(parents=True, exist_ok=True)
    counts, coords = load_curated(); fp = data_fingerprint(counts, coords)
    y_real = counts.to_numpy(int); totals = y_real.sum(1); xy = coords.to_numpy(float); D = geo_km(xy)
    pin = json.loads(PIN.read_text()); alpha_sum = np.asarray(pin["alpha_plus_field_mean"], float)
    sigma, tau, rho = float(pin["sigma"]), float(pin["tau"]), float(pin["rho"])
    m36 = importlib.import_module("36_canonical_phase_map")
    labels, _ = m36.assign_primary_phases([str(i) for i in counts.index], xy); labels = np.asarray(labels)
    names = sorted(set(labels)); phillips = np.array([names.index(v) + 1 for v in labels])
    g1 = pl.kmeans2_regions(xy)
    cands = pl.candidate_lines(xy)
    line = choose_line(labels, cands)
    placements = {"G1": g1, "L": line["regions"]}
    base = json.loads((BASE / "basin_composition.json").read_text())
    a0_fields = arm_fields(base, pin, arm)
    truth_all = {}

    def write(tag, y, regions, truth):
        payload = dict(a0_fields); payload["y"] = y.tolist(); payload.update(a3_fields(regions))
        payload["_truth"] = truth; payload["_provenance"] = {"derived_from": fp, "note": truth.get("kind")}
        (OUT / f"a3_{tag}.json").write_text(json.dumps(payload)); truth_all[tag] = {k: v for k, v in truth.items() if k != "shift"}
        print(f"  {tag}: {truth.get('kind')}, omega {truth.get('omega', 'n/a')}, between-region F_ST {truth['realized_fst_between_regions']:.4f}")

    # 1. own-generator at two placements
    for pname, reg in placements.items():
        for li, om in enumerate(OMEGAS):
            for r in range(N_REC):
                rng = np.random.default_rng(np.random.SeedSequence(SEED_BASE + 100 * li + 10 * list(placements).index(pname) + r))
                y, shift = simulate_a3(D, totals, alpha_sum, sigma, tau, rho, om, reg, rng)
                write(f"own_{pname}_om{om:g}_r{r + 1}", y, reg,
                      {"kind": "own-generator A3", "placement": pname, "omega": om, "shift": shift.tolist(), "sigma": sigma, "tau": tau, "rho": rho,
                       "regions": [int(v) for v in reg], "realized_fst_between_regions": between_fst(y, reg)})
    # 1b. step against slope
    for li, rr in enumerate(SLOPE_RHOS):
        for r in range(N_REC):
            rng = np.random.default_rng(np.random.SeedSequence(SEED_BASE + 1000 + 100 * li + r))
            y, _ = simulate_a3(D, totals, alpha_sum, sigma, tau, rr, 0.0, g1, rng)
            write(f"slope_rho{rr:g}_r{r + 1}", y, g1,
                  {"kind": "step-against-slope (no barrier, moderate range)", "placement": "G1", "omega": 0.0, "sigma": sigma, "tau": tau, "rho": rr,
                   "regions": [int(v) for v in g1], "realized_fst_between_regions": between_fst(y, g1)})
    # 2. copying-simulator records at the primary operating point, T1
    op = json.loads(POINTS.read_text())["primary"]; ts = json.loads(TS.read_text())
    L_d, S, burnin = int(ts["L_d"]), int(ts["probe_stored_length"]), int(ts["probe_burnin"])
    pooled = y_real.astype(float).sum(0); pooled /= pooled.sum()
    for si, (sc, leak) in enumerate(LEAKS.items()):
        w = gen.family_weights(op["kernel"], op["length"], D, labels=None if leak == 1.0 else labels, leak=leak)
        for r in range(N_REC):
            key = SEED_BASE + 2000 + 100 * si + r
            ss = np.random.SeedSequence(key); kids = ss.spawn(4)
            rec = drift_record(w, k=y_real.shape[1], n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"], steps=S, n_records=S,
                               burnin=burnin, seed=int(np.random.default_rng(kids[0]).integers(0, 2**31 - 1)), initial="uniform", target=pooled)
            lengths = spans.draw_lengths(np.random.default_rng(kids[1]), len(D), L_d)
            ends = spans.ends_t1(np.random.default_rng(kids[2]), len(D), S, L_d)
            y = spans.sample_spans(rec, ends, lengths, totals, np.random.default_rng(kids[3]))
            reg = phillips if leak < 1.0 else g1
            write(f"sim_{sc}_r{r + 1}", y, reg,
                  {"kind": f"copying simulator {sc}, leak {leak}, T1, primary operating point", "placement": "Phillips" if leak < 1.0 else "G1",
                   "leak": leak, "operating_point": op, "L_d": L_d, "S": S, "burnin": burnin, "regions": [int(v) for v in reg],
                   "realized_fst_between_regions": between_fst(y, reg), "realized_fst_between_phillips": between_fst(y, phillips)})
    (OUT / "a3_truth.json").write_text(json.dumps(truth_all, indent=1))
    (OUT / "a3_placements.json").write_text(json.dumps({"G1": [int(v) for v in g1], "L": {"angle": line["angle"], "offset": line["offset"], "regions": [int(v) for v in line["regions"]], "agreement_with_parkin_split": line["agreement"]},
                                                         "Phillips": [int(v) for v in phillips], "names": names, "candidates": [{"angle": c["angle"], "offset": c["offset"], "regions": [int(v) for v in c["regions"]]} for c in cands]}, indent=1))
    print(f"wrote {len(truth_all)} records to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
