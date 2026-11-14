"""The two-cluster F_ST of every diversity-matched screen cell, with no selection on F_ST.

The operating points (D22) were chosen among cells matched on the three
diversity summaries (two local, one pooled) by the smallest four-summary loss,
whose fourth term is the squared relative error of the two-cluster F_ST
median, and the stage-2 recheck (resemblance_stage2.py) rechecked only the
matched cells whose F_ST median lay within 50 percent of the observed
value. The two-cluster F_ST therefore entered the choice of operating point
twice, so the main-text comparison of the observed value against the 30
cells at the two operating points is a comparison against cells selected
for relatively high, observed-like F_ST (neither point has the highest), not
against the screen at large. This script reports the comparison the selection did not make: for
every cell of the widened screen (resemblance_screen.py: 2,970 cells, 6
seeds each, each at its own measured time scales) that is matched on the
three diversity summaries and has an accepted time-scale measurement, where
the observed two-cluster F_ST lies relative to the cell's six seeds.

Six seeds per cell give a coarse spread (the 90th percentile of six values
is between the fifth and sixth largest), so the per-cell position is read
as above or not above that spread, and the cell counts are the result. No
new simulation; the screen table is read as written. Descriptive; no
inference rests on it. Exploratory design validation, not evidence about
the record.

Input: output/probe_latent_groups/resemblance_screen.csv (gitignored; written
by resemblance_screen.py). The observed value is recomputed from the counts at
full precision and rounded only for display: comparing against the rounded
0.0062 miscounted two cells (Codex audit 2026-11-03, F5). Output: output/findings/screen_fst_unselected.md.
Usage: .venv/bin/python analyses/probe_latent_groups/screen_fst_unselected.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses", Path(__file__).resolve().parent):
    sys.path.insert(0, str(_p))
SCREEN_CSV = ROOT / "output" / "probe_latent_groups" / "resemblance_screen.csv"
STAGE2_MD = ROOT / "output" / "findings" / "resemblance_stage2.md"
OUT = ROOT / "output" / "findings" / "screen_fst_unselected.md"
PRIMARY = ("gaussian-6km", 2000, 0.002, 0.4)
SECONDARY = ("exponential-6km", 120, 0.024, 0.4)


def observed_fst() -> float:
    from a0_predictive_summaries import spatial_clusters
    from mls_emergence.dataio.basin import load_curated
    from mls_emergence.probe import resemblance as rs
    m47 = importlib.import_module("47_revision_analysis")
    counts, coords = load_curated()
    return float(rs.summaries(counts.to_numpy(int), spatial_clusters(coords.to_numpy(float)), m47.diversity)["fst_clusters"])


def stage2_candidates() -> int:
    m = re.search(r"within 50% of the observed value in the 6-seed screen \((\d+)\)", STAGE2_MD.read_text())
    if m is None:
        raise ValueError(f"{STAGE2_MD}: candidate count not found")
    return int(m[1])


def main() -> None:
    if not SCREEN_CSV.exists():
        raise FileNotFoundError(f"{SCREEN_CSV} missing; run resemblance_screen.py first")
    obs = observed_fst()
    d = pd.read_csv(SCREEN_CSV)
    need = {"family", "n_ind", "innovation", "mixing", "matched", "scales_accepted", "fst_median", "fst_q90",
            "fst_frac_ge_observed", "fst_n_finite"}
    missing = need - set(d.columns)
    if missing:
        raise ValueError(f"{SCREEN_CSV}: missing columns {sorted(missing)}")
    m = d[d.matched.astype(bool) & d.scales_accepted.astype(bool)].copy()
    if (m.fst_n_finite < 6).any():
        raise ValueError("a matched cell has fewer than six finite F_ST values")
    above_q90 = m.fst_q90 < obs                      # observed above the cell's 90th percentile of six seeds
    above_med = m.fst_median < obs
    no_seed = m.fst_frac_ge_observed <= 0            # no seed of six at or above the observed value
    n_stage2 = stage2_candidates()

    def is_point(r, pt):
        return r.family == pt[0] and int(r.n_ind) == pt[1] and np.isclose(r.innovation, pt[2]) and np.isclose(r.mixing, pt[3])
    pts = {name: m[[is_point(r, pt) for r in m.itertuples()]] for name, pt in (("primary", PRIMARY), ("secondary", SECONDARY))}
    for name, sub in pts.items():
        if len(sub) != 1:
            raise ValueError(f"{name} operating point matched {len(sub)} screen rows")
    rank = {name: int((m.fst_median > float(sub.fst_median.iloc[0])).sum()) + 1 for name, sub in pts.items()}

    q = m.fst_median.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(4)
    fam = m.groupby("family").agg(n=("fst_median", "size"), median_of_medians=("fst_median", "median"),
                                  observed_above_q90=("fst_q90", lambda x: int((x < obs).sum())),
                                  median_at_or_above_observed=("fst_median", lambda x: int((x >= obs).sum())))
    L = ["# The two-cluster F_ST across every diversity-matched screen cell, without the F_ST selection", "",
         f"Produced by `analyses/probe_latent_groups/screen_fst_unselected.py` from `resemblance_screen.csv` (2,970 cells, six seeds each). "
         f"Observed two-cluster F_ST {obs:.4f}, compared at full precision. Descriptive; exploratory design validation, not evidence about the record.", "",
         "## What the selection did",
         "",
         f"Of {len(d)} screen cells, {int(d.matched.astype(bool).sum())} match the three diversity summaries within 10 percent and "
         f"{len(m)} of those have an accepted time-scale measurement. The stage-2 recheck took the {n_stage2} of these whose six-seed F_ST median "
         f"lay within 50 percent of the observed value, and the operating points were chosen among them by the four-summary loss, whose fourth term is "
         f"the squared relative error of that median (D22). By six-seed median the primary point ranks {rank['primary']} of {len(m)} matched cells "
         f"and the secondary {rank['secondary']} of {len(m)}.", "",
         "## Where the observed value lies in the unselected matched cells", "",
         "| quantity | value |", "|---|---|",
         f"| matched cells with accepted time scales | {len(m)} |",
         f"| six-seed F_ST median: 10th, 25th, 50th, 75th, 90th percentile over cells | {q[0.1]:.4f}, {q[0.25]:.4f}, {q[0.5]:.4f}, {q[0.75]:.4f}, {q[0.9]:.4f} |",
         f"| cells whose six-seed median is at or above the observed value | {int((~above_med).sum())} |",
         f"| cells whose six-seed 90th percentile is below the observed value (observed above that percentile) | {int(above_q90.sum())} ({100 * above_q90.mean():.0f} percent) |",
         f"| cells with no seed of six at or above the observed value | {int(no_seed.sum())} ({100 * no_seed.mean():.0f} percent) |",
         f"| mean fraction of seeds at or above the observed value, over cells | {m.fst_frac_ge_observed.mean():.3f} |",
         "", "## By copying-network family", "",
         "| family | matched cells | median of six-seed medians | cells with observed above the 90th percentile | cells with median at or above observed |",
         "|---|---|---|---|---|"]
    L += [f"| {f} | {int(r.n)} | {r.median_of_medians:.4f} | {int(r.observed_above_q90)} | {int(r.median_at_or_above_observed)} |" for f, r in fam.iterrows()]
    top = m.sort_values("fst_median", ascending=False).head(10)
    L += ["", "## The ten matched cells with the largest six-seed F_ST median", "",
          "| family | learners | innovation | mixing | F_ST median | 90th percentile | seeds at or above observed |", "|---|---|---|---|---|---|---|"]
    L += [f"| {r.family} | {int(r.n_ind)} | {r.innovation:g} | {r.mixing:g} | {r.fst_median:.4f} | {r.fst_q90:.4f} | {r.fst_frac_ge_observed:.2f} |" for r in top.itertuples()]
    L += ["", "Reading: a comparison of the observed two-cluster F_ST against the full set of diversity-matched cells, rather than against the cells "
          "selected for their closeness to it, places the observed value above the six-seed 90th percentile of most cells. The 50-seed bands at the two operating "
          "points (17 of 30 cells containing the observed value) describe cells selected for relatively high, observed-like F_ST, although neither point has the highest. "
          "Six seeds per cell resolve the position coarsely, so the counts above are stated as counts of cells, not as a probability."]
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L[4:22]))


if __name__ == "__main__":
    main()
