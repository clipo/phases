"""Choose the S0-stress cell within the primary family with a 50-seed recheck (D22; spec 3.2, 3.3).

Latent-group probe, pilot step 5 follow-up. Exploratory design validation
(rule 1). Candidates are the primary family's cells matched on diversity in
the 6-seed screen, ordered by F_ST median; the five lowest are rerun with 50
fresh seeds at the primary point's time scales (resemblance_stage2 machinery);
the stress cell is the rechecked cell that is still matched with the smallest
F_ST median, ties by loss then grid order.

Usage: .venv/bin/python analyses/probe_latent_groups/stress_cell.py
Output: output/probe_latent_groups/stress_cell.json, output/findings/stress_cell.md
"""
from __future__ import annotations

import json
import multiprocessing as mp
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import resemblance_stage2 as st

from mls_emergence.probe import resemblance as rs

N_LOWEST = 5


def main() -> int:
    points = json.loads((Path(__file__).resolve().parent / "operating_points.json").read_text())
    prim = points["primary"]; family = f"{prim['kernel']}-{prim['length']:g}km"
    if len(sys.argv) > 1:          # a family name overrides the file (provisional runs before a decision is logged)
        family = sys.argv[1]
    st.SEED_BASE = 2026101680
    # Reuse stage 2's setup by calling its module-level pieces.
    import importlib
    m47 = importlib.import_module("47_revision_analysis")
    from a0_predictive_summaries import spatial_clusters

    from mls_emergence.dataio.basin import geo_km, load_curated
    from mls_emergence.probe import generators as gen
    counts, coords = load_curated(); xy = coords.to_numpy(float); D = geo_km(xy)
    st.G.update(N=len(counts), K=counts.shape[1], weights={gen.family_name(k, l): gen.family_weights(k, l, D) for k, l in gen.FAMILIES},
                pooled=counts.to_numpy(float).sum(0) / counts.to_numpy(float).sum(), totals=counts.to_numpy(int).sum(1),
                labels=spatial_clusters(xy), diversity=m47.diversity)
    st.screen.G.update(st.G)
    st.G["observed"] = rs.summaries(counts.to_numpy(int), st.G["labels"], m47.diversity)
    screen = pd.read_csv(st.OUTD / "resemblance_screen.csv")
    cand = screen[(screen.family == family) & screen.matched].sort_values(["fst_median", "loss"]).head(N_LOWEST)
    args = [(c, r.family, int(r.n_ind), float(r.innovation), float(r.mixing)) for c, r in enumerate(cand.itertuples())]
    with mp.get_context("fork").Pool(len(args)) as pool:
        rows = pool.map(st.run_candidate, args)
    for r, g in zip(rows, cand.grid_index):
        r["screen_grid_index"] = int(g)
    pick = rs.select_stress_cell(rows)
    tab = pd.DataFrame(rows)
    out = {"family": family, "candidates": rows, "stress_cell": pick}
    (st.OUTD / "stress_cell.json").write_text(json.dumps(out, indent=1, default=float))
    L = ["# S0-stress cell within the primary family (D22)", "",
         (f"Produced by `analyses/probe_latent_groups/stress_cell.py`: the {N_LOWEST} lowest-F_ST matched cells of the {family} "
          f"family in the 6-seed screen, rerun with {st.N_SEEDS} fresh seeds, each at its own measured time scales. "
          "Exploratory design validation, not evidence about the record."), "",
         "| n_ind | innovation | mixing | hs (se) | rich (se) | ht (se) | matched | F_ST median [10%, 90%] | loss |", "|---|---|---|---|---|---|---|---|---|",
         *[f"| {int(r.n_ind)} | {r.innovation:g} | {r.mixing:g} | {r.mean_hs:.4f} ({r.se_hs:.4f}) | {r.mean_rich:.3f} ({r.se_rich:.3f}) | "
           f"{r.mean_ht:.4f} ({r.se_ht:.4f}) | {'yes' if r.matched else 'no'} | {r.fst_median:.4f} [{r.fst_q10:.4f}, {r.fst_q90:.4f}] | {r.loss:.4f} |"
           for r in tab.itertuples()], ""]
    L.append("No rechecked cell stayed matched; stop for discussion." if pick is None else
             f"**S0-stress = n_ind {int(pick['n_ind'])}, innovation {pick['innovation']:g}, mixing {pick['mixing']:g}** "
             f"(F_ST median {pick['fst_median']:.4f}, loss {pick['loss']:.4f}).")
    (ROOT / "output" / "findings" / "stress_cell.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[4:]))
    return 0 if pick is not None else 3


if __name__ == "__main__":
    raise SystemExit(main())
