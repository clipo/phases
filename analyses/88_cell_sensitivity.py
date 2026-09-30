"""88_cell_sensitivity.py - does the drift shortfall depend on which matched cell is used?

The cell most favorable to drift is not stable across seed families (analysis
47's stage-2 winner is mixing 0.02; `measure_calibration_stability.py`'s
50-seed winner is mixing 0.01; their medians differ by about 0.0005). The main
text reports the drift comparison at one cell, so this script scores the same
statistics at the top six cells of analysis 47's stage-2 ranking, and at one
cell that matched on the six-run screen but failed the fifty-run diversity
check, reported because it gives the drift model its highest values: between-cluster F_ST at k = 2
to 5 (the divisions of analysis 71), the between-phase F_ST, and the
Parkin-versus-rest F_ST, 300 realizations per cell.

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
CELLS = [  # (learners, innovation, mixing): the top six of analysis 47's stage-2 ranking
    (2000, 0.001, 0.02),     # the selection
    (10000, 0.0005, 0.002),
    (120, 0.008, 0.1),
    (10000, 0.0002, 0.005),
    (2000, 0.002, 0.005),
    (2000, 0.0005, 0.05),
    (2000, 0.001, 0.01),     # NOT confirmed: matched on the six-run screen, failed the fifty-run
                             # diversity check; first in measure_calibration_stability.py's F_ST ranking
]
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
    labs, _ = ph.assign_phases_by_territory(names, coords.to_numpy(float))
    plist = sorted(set(labs))
    labels["phases"] = np.array([plist.index(l) for l in labs])
    labels["Parkin vs rest"] = np.array([1 if l == "Parkin" else 0 for l in labs])
    m_obs = counts.to_numpy(int)
    obs = {k: rev.fst_by(m_obs, v) for k, v in labels.items()}
    totals = data["m"].sum(1)
    L = ["# Does the drift shortfall depend on which matched cell is used?", "",
         f"Produced by `analyses/88_cell_sensitivity.py`, {args.reps} realizations per cell, "
         "pooled-profile model, basin, 28 assemblages. Each cell is a diversity-matched cell of "
         "analysis 47's grid. Entries: drift median [2.5th, 97.5th percentile]; share of runs "
         "at or above the observed value.", "",
         "| cell (learners, innovation, mixing) | " + " | ".join(f"{k} (obs {obs[k]:.4f})" for k in labels) + " |",
         "|---|" + "---|" * len(labels)]
    for cell in CELLS:
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
        L.append(f"| {n_ind}, {inn}, {mix} | " + " | ".join(cells) + " |")
        print(L[-1], flush=True)
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
