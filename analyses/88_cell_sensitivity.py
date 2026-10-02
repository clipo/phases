"""88_cell_sensitivity.py - does the drift shortfall depend on which matched cell is used?

The cell most favorable to drift is not stable across seed families (analysis
47's fifty-run ranking and `measure_calibration_stability.py`'s fresh-seed
ranking put different cells first, with medians a few ten-thousandths apart).
The main text reports the drift comparison at one cell, so this script scores
the same statistics at the six cells that lead analysis 47's fifty-run ranking
among those that passed its diversity confirmation: between-cluster F_ST at
k = 2 to 5 (the divisions of analysis 71), the between-phase F_ST, and the
Parkin-versus-rest F_ST, 300 realizations per cell.

THE SIX ARE READ FROM THE CALIBRATION, not listed here (rule 4). Until
2026-10-02 they were a hand-kept list from the calibration of 2026-09-23, with
a seventh cell that had matched on the six-run screen and failed the fifty-run
check. After the Cummins coordinate was corrected the calibration ranked
different cells, three of the listed six were no longer among its leaders, the
seventh no longer matched even on the screen, and the script went on scoring
the old list without complaint.

Output: output/findings/cell_sensitivity.md
Usage: .venv/bin/python analyses/88_cell_sensitivity.py [--reps 300]
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
OUT_MD = ROOT / "output" / "findings" / "cell_sensitivity.md"
CALIBRATION = ROOT / "output" / "revision_2026_09" / "calibration.csv"
N_CELLS = 6


def leading_cells(n: int = N_CELLS):
    """The `n` cells leading analysis 47's fifty-run ranking among those that
    passed its diversity confirmation (basin, pooled-profile model), the
    calibration's own selection first. Raises if the calibration is missing or
    its selection is not its own first-ranked confirmed cell."""
    import pandas as pd
    if not CALIBRATION.exists():
        raise SystemExit(f"{CALIBRATION} is missing; run analyses/47_revision_analysis.py")
    c = pd.read_csv(CALIBRATION)
    b = c[(c["region"] == "basin") & (c["model"] == "pooled") & c["stage2_matched"].astype(bool)]
    b = b.sort_values("stage2_fst_median", ascending=False)
    if len(b) < n:
        raise SystemExit(f"only {len(b)} confirmed cells in the calibration; {n} wanted")
    if not bool(b.iloc[0]["selected"]):
        raise SystemExit("the calibration's selected cell is not its first-ranked confirmed cell")
    return [(int(r.n_ind), float(r.innovation), float(r.mixing)) for r in b.head(n).itertuples()]


K_RANGE = (2, 3, 4, 5)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()
    mf = importlib.import_module("make_figures")
    rev = importlib.import_module("47_revision_analysis")
    ph = importlib.import_module("36_canonical_phase_map")
    counts, coords = mf._load_curated()
    data = rev.load_sets()["basin"]
    names = [str(i) for i in counts.index]
    if [str(n) for n in data["names"]] != names:
        raise RuntimeError("assemblage order differs between loaders")
    centred = coords.to_numpy(float) - coords.to_numpy(float).mean(0)
    labels = {f"k={k}": mf._kmeans_labels(centred, k, seed=7) for k in K_RANGE}
    labs, _ = ph.assign_primary_phases(names, coords.to_numpy(float))
    plist = sorted(set(labs))
    labels["phases"] = np.array([plist.index(l) for l in labs])
    labels["Parkin vs rest"] = np.array([1 if l == "Parkin" else 0 for l in labs])
    m_obs = counts.to_numpy(int)
    obs = {k: rev.fst_by(m_obs, v) for k, v in labels.items()}
    totals = data["m"].sum(1)
    L = ["# Does the drift shortfall depend on which matched cell is used?", "",
         f"Produced by `analyses/88_cell_sensitivity.py`, {args.reps} realizations per cell, "
         "pooled-profile model, basin, 28 assemblages. The cells are the six that lead analysis 47's "
         "fifty-run ranking among those passing its diversity confirmation, the selection first. "
         "Entries: drift median [2.5th, 97.5th percentile]; share of runs "
         "at or above the observed value.", "",
         "| cell (learners, innovation, mixing) | " + " | ".join(f"{k} (obs {obs[k]:.4f})" for k in labels) + " |",
         "|---|" + "---|" * len(labels)]
    for cell in leading_cells():
        n_ind, inn, mix = cell
        sims = {k: [] for k in labels}
        for seed in range(args.reps):
            f = rev.simulate(data, 90000 + seed, "pooled", innovation=inn, mixing=mix, n_ind=n_ind)
            m = rev.sample_record(f, data["ranks"], totals, np.random.default_rng(91000 + seed))
            for k, v in labels.items():
                sims[k].append(rev.fst_by(m, v))
        cells = []
        for k in labels:
            d = np.array(sims[k], float); d = d[np.isfinite(d)]
            lo, med, hi = np.percentile(d, [2.5, 50, 97.5])
            cells.append(f"{med:.4f} [{lo:.4f}, {hi:.4f}]; {np.mean(d >= obs[k]) * 100:.0f}%")
        L.append(f"| {n_ind}, {inn:g}, {mix:g} | " + " | ".join(cells) + " |")
        print(L[-1], flush=True)
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
