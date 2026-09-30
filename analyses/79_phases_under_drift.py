#!/usr/bin/env python3
"""Would a culture historian have found these phases in pottery made by drift?

Earlier work showed that the late-period phases of this region do not hold
together as units (Fox 1998; Lipo 2001; Mainfort 2003), and this paper agrees.
What none of that work explains is why careful workers saw social groups here
in the first place. This script tests one account of that directly.

THE ACCOUNT. Similarity in decorated pottery decays with distance, which
unbounded copying among neighbours produces without any group. Settlement is
clumped. Sampling a smooth gradient at clumped points gives assemblages that
are alike within a clump and different across the gap to the next one. A
procedure that groups assemblages by resemblance and by where they are, which
is what the culture-historical method did, will then return groups, and they
will look like the published phases because the clumps are where they are.

THE TEST. Run the calibrated neutral drift model on the real site geography.
It contains no groups, no boundaries and no restricted copying. Sample its
pottery at the observed sherd totals. Then treat each simulated record exactly
as analyses/76_phase_recovery.py treats the real one: cluster the assemblages
in the joint space of position and chi-square composition, at several weights
on composition, and score agreement with the published phases by adjusted Rand
index on the assemblages Mainfort himself assigned.

WHAT THE OUTCOMES WOULD MEAN.

  Drift pottery recovers the published phases about as well as the real
  pottery does. Then the appearance of phases needs nothing beyond distance
  decay and clumped settlement, and the social reading added nothing the null
  does not supply.

  The real pottery recovers them clearly better than drift pottery does. Then
  the real record carries phase-aligned structure that drift on this geography
  lacks, and the account above is incomplete.

The weight-zero column is the site map alone and is identical in both by
construction; it is printed as a check. The informative columns are the ones
where pottery enters, above all pottery alone, where geography cannot help.

Drift's pottery is scored over many realizations, so the observed value is
read against a distribution and reported with where it falls in it. That
placement is descriptive (rule 18): it says what the design manufactures under
the null, and no inference is drawn from a tail probability.

Usage:
    python analyses/79_phases_under_drift.py [--reps 300] [--seeds 10]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "phases_under_drift.md"
OUT_CSV = ROOT / "output" / "findings" / "phases_under_drift.csv"
FIG = "fig15_phases_under_drift"
MODEL = "pooled"
WEIGHTS = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, np.inf]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--seeds", type=int, default=10,
                    help="k-means initializations per clustering; the median is scored")
    args = ap.parse_args()

    import matplotlib.pyplot as plt
    fs = importlib.import_module("figstyle")
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    ph = importlib.import_module("36_canonical_phase_map")
    t74 = importlib.import_module("74_phase_partition_test")
    t76 = importlib.import_module("76_phase_recovery")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m_obs = counts.to_numpy(float)
    xy = coords.to_numpy(float)
    labels_ph, derived = ph.assign_phases_by_territory(names, xy)
    mapped = ~np.array([bool(d) for d in derived])
    phases = sorted(set(labels_ph))
    pidx = np.array([phases.index(l) for l in labels_ph])
    k = len(phases)
    G = t76.unit_scale(t74.km_xy(xy))

    data = rev.load_sets()["basin"]
    if [str(x) for x in data["names"]] != names:
        raise ValueError("simulation set and curated set list assemblages differently")
    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == MODEL)].iloc[0]
    rates = dict(innovation=float(row.innovation), mixing=float(row.mixing),
                 n_ind=int(row.n_ind))
    totals = data["m"].sum(1)

    def recover(m: np.ndarray) -> list[float]:
        """Agreement with the published phases at each weight, for one record."""
        tot = m.sum(1, keepdims=True)
        p = m / np.where(tot > 0, tot, 1.0)
        C = t76.unit_scale(t76.composition_features(p, "chisq"))
        out = []
        for w in WEIGHTS:
            F = G if w == 0 else (C if not np.isfinite(w) else
                                  np.column_stack([G, np.sqrt(w) * C]))
            a = [t76.ari(pidx[mapped], mf._kmeans_labels(F, k, seed=s)[mapped])
                 for s in range(args.seeds)]
            out.append(float(np.median(a)))
        return out

    obs = recover(m_obs)
    sims = []
    for seed in range(args.reps):
        f = rev.simulate(data, 79000 + seed, MODEL, **rates)
        m = rev.sample_record(f, data["ranks"], totals,
                              np.random.default_rng(79500 + seed)).astype(float)
        sims.append(recover(m))
        if (seed + 1) % 25 == 0:
            print(f"  {seed + 1}/{args.reps}", flush=True)
    sims = np.array(sims)

    rows = []
    for j, w in enumerate(WEIGHTS):
        d = sims[:, j]
        rows.append(dict(weight="pottery alone" if not np.isfinite(w) else w,
                         observed=obs[j], drift_median=float(np.median(d)),
                         drift_lo=float(np.percentile(d, 2.5)),
                         drift_hi=float(np.percentile(d, 97.5)),
                         share_of_drift_at_or_above_observed=float(np.mean(d >= obs[j]))))
    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    x = np.arange(len(WEIGHTS))
    ax.fill_between(x, df.drift_lo, df.drift_hi, color="0.82", label="pottery made by drift, 95%")
    ax.plot(x, df.drift_median, color="0.45", lw=1.2, ls="--", label="drift median")
    ax.plot(x, df.observed, color="0.1", lw=1.5, marker="o", ms=3.5, label="the real pottery")
    ax.set_xticks(x)
    ax.set_xticklabels(["map\nalone" if w == 0 else ("pottery\nalone" if not np.isfinite(w)
                        else f"{w:g}") for w in WEIGHTS], fontsize=6)
    ax.set_xlabel("weight on composition relative to geography")
    ax.set_ylabel("agreement with the published phases\n(adjusted Rand index)")
    ax.legend(fontsize=5.5, frameon=False, loc="upper right")
    png = fs.save_all(fig, FIG, close=True)

    L = ["# Would a culture historian have found these phases in pottery made by drift?", "",
         f"Basin phase set, {len(names)} assemblages, {k} phases; agreement scored on the "
         f"{int(mapped.sum())} assemblages Mainfort assigned.",
         f"{args.reps} realizations of the calibrated {MODEL} drift model "
         f"({rates['n_ind']} learners, innovation {rates['innovation']}, mixing "
         f"{rates['mixing']}), which contains no groups,",
         "sampled at the observed sherd totals and clustered exactly as the real record is in "
         "`76_phase_recovery.py`",
         f"(chi-square composition, k-means, median of {args.seeds} initializations).", "",
         "| weight on composition | real pottery | drift pottery, median (95 percent) | "
         "share of drift records at or above the real one |", "|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['weight']} | **{r['observed']:.3f}** | {r['drift_median']:.3f} "
                 f"({r['drift_lo']:.3f} to {r['drift_hi']:.3f}) | "
                 f"{100 * r['share_of_drift_at_or_above_observed']:.0f}% |")
    L += ["", "The weight-zero row is the site map alone and is the same in both columns by "
              "construction.", "",
          "Read the rows where pottery enters. Where the real pottery sits inside what drift "
          "pottery gives,",
          "the published phases are as recoverable from a record with no social groups in it "
          "as from the",
          "real one, and the appearance of phases needs nothing beyond distance decay and "
          "where the sites are.",
          "Where the real pottery sits above it, the record carries phase-aligned structure "
          "that drift on",
          "this geography does not produce.", "",
          f"Figure written to {png.name} and its siblings.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n" + df.to_string(index=False))
    print(f"\nwrote {png} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
