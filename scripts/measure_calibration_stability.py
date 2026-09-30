"""measure_calibration_stability.py - is the "cell most favorable to drift" identifiable?

Analysis 47 calibrates the drift model on a 297-cell grid, six seeds per cell,
keeps the cells whose three diversity summaries fall within 10 percent of the
observed values, and takes as the baseline the matched cell with the HIGHEST
six-seed median partition F_ST. That is the maximum of 50-odd noisy medians.
On 2026-09-22 a 209 m datum shift in 13 of 28 coordinates, which left every
partition and plug-in F_ST unchanged, moved the winner from (120 learners,
innovation 0.024, mixing 0.01) to (2,000, 0.001, 0.01), and with it every
downstream drift comparison. This script re-simulates every matched cell with
REPS fresh seeds and reports whether the six-seed ranking survives, which cell
wins on the larger sample, and the envelope across matched cells. It changes
nothing in the pipeline; it writes `output/findings/calibration_stability.md`.

Usage: .venv/bin/python scripts/measure_calibration_stability.py [--reps 50]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analyses"))
a47 = importlib.import_module("47_revision_analysis")

CAL = ROOT / "output" / "revision_2026_09" / "calibration.csv"
OUT = ROOT / "output" / "findings" / "calibration_stability.md"
SEED0 = 90000   # disjoint from every seed family 47 uses (60000/61000 calibration, 0.. baseline, 30000, 50000)


def main(reps):
    sets = a47.load_sets()
    data = sets["basin"]
    part = data["parkin"] if data["parkin"] is not None else data["labels"]
    totals = data["m"].sum(1)
    obs = a47.fst_by(data["m"], part)
    cal = pd.read_csv(CAL)
    cal = cal[(cal.region == "basin") & (cal.model == "pooled")]
    matched = cal[cal.matched].copy().reset_index(drop=True)
    rows = []
    for i, r in matched.iterrows():
        f = []
        for s in range(reps):
            field = a47.simulate(data, SEED0 + 1000 * i + s, "pooled", innovation=float(r.innovation),
                                 mixing=float(r.mixing), n_ind=int(r.n_ind))
            m = a47.sample_record(field, data["ranks"], totals, np.random.default_rng(SEED0 + 500000 + 1000 * i + s))
            f.append(a47.fst_by(m, part))
        f = np.asarray(f, float); f = f[np.isfinite(f)]
        rows.append(dict(n_ind=int(r.n_ind), innovation=float(r.innovation), mixing=float(r.mixing),
                         six_seed_median=float(r.fst_median), median=float(np.median(f)),
                         lo=float(np.percentile(f, 2.5)), hi=float(np.percentile(f, 97.5)),
                         covers_observed=bool(np.percentile(f, 2.5) <= obs <= np.percentile(f, 97.5)),
                         share_reaching=float(np.mean(f >= obs)), n=int(len(f))))
        print(f"{i+1}/{len(matched)} n={r.n_ind:>5.0f} inn={r.innovation:.4f} mix={r.mixing:.3f} "
              f"six-seed {r.fst_median:.4f} -> {reps}-seed {np.median(f):.4f}", flush=True)
    df = pd.DataFrame(rows)
    rho = df[["six_seed_median", "median"]].corr(method="spearman").iloc[0, 1]
    six_win = df.loc[df.six_seed_median.idxmax()]
    big_win = df.loc[df["median"].idxmax()]
    top5_six = set(df.sort_values("six_seed_median", ascending=False).head(5).index)
    top5_big = set(df.sort_values("median", ascending=False).head(5).index)
    L = ["# Is the cell most favorable to drift identifiable?", "",
         f"Produced by `scripts/measure_calibration_stability.py`, {reps} fresh seeds per matched "
         f"cell (seed family {SEED0}), pooled-profile model, basin, {len(df)} diversity-matched cells "
         f"from `output/revision_2026_09/calibration.csv`. Partition F_ST is the statistic analysis 47 "
         f"ranks cells by; observed {obs:.4f}.", "",
         "## Does the six-seed ranking survive?", "",
         f"- Spearman rank correlation between the six-seed medians and the {reps}-seed medians across "
         f"matched cells: **{rho:+.2f}**.",
         f"- Six-seed winner: ({six_win.n_ind:.0f} learners, innovation {six_win.innovation}, mixing "
         f"{six_win.mixing}); six-seed median {six_win.six_seed_median:.4f}, {reps}-seed median "
         f"**{six_win['median']:.4f}** [{six_win.lo:.4f}, {six_win.hi:.4f}].",
         f"- {reps}-seed winner: ({big_win.n_ind:.0f} learners, innovation {big_win.innovation}, mixing "
         f"{big_win.mixing}); median **{big_win['median']:.4f}** [{big_win.lo:.4f}, {big_win.hi:.4f}].",
         f"- Overlap of the top five cells under the two rankings: {len(top5_six & top5_big)} of 5.", "",
         "## The envelope across matched cells", "",
         f"- {reps}-seed medians run {df['median'].min():.4f} to {df['median'].max():.4f} "
         f"(median of medians {df['median'].median():.4f}).",
         f"- Cells whose 95 percent interval covers the observed {obs:.4f}: "
         f"**{int(df.covers_observed.sum())} of {len(df)}**.",
         f"- Share of realizations reaching the observed value, across cells: median "
         f"{df.share_reaching.median():.1%}, max {df.share_reaching.max():.1%}.", "",
         "## Every matched cell", "",
         f"| learners | innovation | mixing | six-seed median | {reps}-seed median | 95% | covers observed | share reaching |",
         "|---|---|---|---|---|---|---|---|"]
    for _, r in df.sort_values("median", ascending=False).iterrows():
        L.append(f"| {r.n_ind:.0f} | {r.innovation} | {r.mixing} | {r.six_seed_median:.4f} | {r['median']:.4f} | "
                 f"[{r.lo:.4f}, {r.hi:.4f}] | {'yes' if r.covers_observed else 'no'} | {r.share_reaching:.1%} |")
    L += ["", "## Reading", "",
          "If the rank correlation is low and the winners differ, the argmax over six-seed medians is "
          "selecting sampling noise, and the single-cell drift comparison that the manuscript reports "
          "is not a property of the calibration but of the seed. The envelope rows above are what the "
          "matched set as a whole says.", ""]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:20]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--reps", type=int, default=50)
    main(p.parse_args().reps)
