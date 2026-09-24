"""89_partition_power.py - could the phase-line comparison see a copying boundary if there were one?

The main text argues from analysis 86 that the phase lines separate the pottery
no better than arbitrary lines of the same sizes. That argument has force only
if the comparison can detect lines that DO mark restricted copying. This script
measures its power (a rule 20c-style recovery check, 2026-09-23). It rebuilds
simulated records of analysis 84's grid, with the same seeds, at cells with and
without a copying boundary at the phase lines, and asks of each simulated
record the question 86 asks of the real one: what share of same-size
alternative divisions (around random centers, and made compact) do the phase
lines beat, on between-group F_ST and on the boundary excess? The same is asked
of the Parkin-versus-rest line, which is part of the phase boundary and so is
restricted whenever the phase lines are.

Plug-in statistics on each simulated record: the simulated counts are the truth
here, so no posterior draw of the counts is needed. The observed record's value
of each share is reported alongside (plug-in, 86's alternatives).

Output: output/findings/partition_power.md
Usage: .venv/bin/python analyses/89_partition_power.py [--reps 100] [--alt 50]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
OUT_MD = ROOT / "output" / "findings" / "partition_power.md"
CELLS = [(1.0, 0.0), (0.1, 0.0), (0.03, 0.0), (1.0, 0.2), (0.03, 0.2)]  # (copying factor, local innovation)


def alternatives(t74, pts, labels, n_alt, seed):
    sizes = np.bincount(labels); k = len(sizes); slot = np.repeat(np.arange(k), sizes)
    rng = np.random.default_rng(seed)
    alt, opt = [], []
    for _ in range(n_alt):
        a = t74.assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)
        alt.append(a); opt.append(t74.local_search(pts, a, sizes))
    return alt, opt


def shares(m, labels, alt, opt, rev, sd, d):
    f = rev.fst_by(m, labels); b = sd.boundary_excess_labeled(m, d, labels)
    fa = np.array([rev.fst_by(m, a) for a in alt]); fo = np.array([rev.fst_by(m, a) for a in opt])
    ba = np.array([sd.boundary_excess_labeled(m, d, a) for a in alt])
    bo = np.array([sd.boundary_excess_labeled(m, d, a) for a in opt])
    return dict(f_alt=np.mean(f > fa), f_opt=np.mean(f > fo), b_alt=np.mean(b > ba), b_opt=np.mean(b > bo))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--alt", type=int, default=50)
    args = ap.parse_args()
    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    sd = importlib.import_module("23_phases_vs_spatial_drift")
    t74 = importlib.import_module("74_phase_partition_test")
    a84 = importlib.import_module("84_phases_as_groups")
    import pandas as pd
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record

    data = rev.load_sets()["basin"]
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    if names != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    labels_ph, _ = ph.assign_phases_by_territory(names, coords.to_numpy(float))
    plist = sorted(set(labels_ph)); phase = np.array([plist.index(l) for l in labels_ph])
    parkin = np.array([1 if l == "Parkin" else 0 for l in labels_ph])
    pts = t74.km_xy(coords.to_numpy(float))
    alt_ph, opt_ph = alternatives(t74, pts, phase, args.alt, 89000)
    alt_pk, opt_pk = alternatives(t74, pts, parkin, args.alt, 89200)
    m_obs = data["m"]; d = data["d"]; totals = m_obs.sum(1)
    rates = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates[(rates.region == "basin") & (rates.model == "pooled")].iloc[0]
    cell = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))
    obs_ph = shares(m_obs.astype(float), phase, alt_ph, opt_ph, rev, sd, d)
    obs_pk = shares(m_obs.astype(float), parkin, alt_pk, opt_pk, rev, sd, d)

    L = ["# Could the phase-line comparison see a copying boundary if there were one?", "",
         f"Produced by `analyses/89_partition_power.py`: {args.reps} simulated records per cell, rebuilt with "
         f"analysis 84's seeds at its calibrated cell (N {cell['n_ind']}, innovation {cell['innovation']}, "
         f"mixing {cell['mixing']}, 24 km along the rivers); {args.alt} same-size alternative divisions of each kind. "
         "Each entry is the median over simulated records of the share of alternatives the tested line beats, "
         "with the share of records in which it beats at least 90 percent of them in parentheses.", "",
         "| copying factor at the phase lines | local innovation | line tested | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |",
         "|---|---|---|---|---|---|---|"]
    fmt = lambda v: f"{np.median(v):.2f} ({np.mean(np.array(v) >= 0.9) * 100:.0f}%)"
    for leak, strength in CELLS:
        li, si = a84.LEAKS.index(leak), a84.STRENGTHS.index(strength)
        w = copying_weights(d, 24.0, labels=phase, leak=leak)
        res = {"phases": [], "Parkin vs rest": []}
        for r in range(args.reps):
            seed = 84000 + 10000 * li + 1000 * si + r      # analysis 84's seeds
            rng = np.random.default_rng(seed + 5_000_000)
            tgt = a84.group_targets_by(data["pooled"], phase, strength, rng)
            rec = drift_record(w, k=m_obs.shape[1], n_ind=cell["n_ind"], seed=seed,
                               innovation=cell["innovation"], mixing=cell["mixing"],
                               burnin=1200, initial="uniform", target=tgt)
            m = sample_record(rec, data["ranks"], totals, rng).astype(float)
            res["phases"].append(shares(m, phase, alt_ph, opt_ph, rev, sd, d))
            res["Parkin vs rest"].append(shares(m, parkin, alt_pk, opt_pk, rev, sd, d))
        for line, rr in res.items():
            L.append(f"| {leak} | {strength} | {line} | " + " | ".join(fmt([x[c] for x in rr]) for c in ("f_alt", "f_opt", "b_alt", "b_opt")) + " |")
            print(L[-1], flush=True)
    L += ["", "The observed record (plug-in, these alternatives):", "",
          "| line | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |",
          "|---|---|---|---|---|",
          "| phases | " + " | ".join(f"{obs_ph[c]:.2f}" for c in ("f_alt", "f_opt", "b_alt", "b_opt")) + " |",
          "| Parkin vs rest | " + " | ".join(f"{obs_pk[c]:.2f}" for c in ("f_alt", "f_opt", "b_alt", "b_opt")) + " |", "",
          "## Reading", "",
          "Compare the rows with a copying boundary (factor below 1) against the row without. If the tested line "
          "beats most alternatives only when copying is restricted at it, the comparison has power, and the "
          "observed record's failure to beat them is evidence against a restriction there. If the rows look alike, "
          "the comparison cannot see a copying boundary and the argument from it should be dropped.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
