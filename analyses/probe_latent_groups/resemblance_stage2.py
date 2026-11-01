"""Stage-2 recheck of candidate operating points with 50 fresh seeds (D21).

Latent-group probe, pilot step 5 follow-up. Exploratory design validation
(rule 1). The widened screen (resemblance_screen.py) uses 6 seeds per cell,
and the inherited operating point moved across the 10 percent richness
tolerance between two seed draws (0.034 against 0.106), so 6 seeds are too
few to decide an operating point. This script re-runs, with 50 fresh seeds
each, every cell that was matched on the three diversity summaries and whose
two-cluster F_ST median lay within 25 percent of the observed value, plus the
inherited operating point, and reports means with their across-seed standard
errors, the F_ST distribution, and the four-summary loss (three diversity
relative errors plus the F_ST relative error, squared and summed). It makes
no selection: the operating point is a logged decision.

Seeds: candidate c uses SeedSequence(2026101670 + c) with 50 run children,
three stream children and five time-scale children, as in the screen. Each
candidate measures its own time scales (screen.cell_scales) first.
Usage: .venv/bin/python analyses/probe_latent_groups/resemblance_stage2.py
Output: output/probe_latent_groups/resemblance_stage2.csv, output/findings/resemblance_stage2.md
"""
from __future__ import annotations

import importlib
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses")); sys.path.insert(0, str(Path(__file__).resolve().parent))

import resemblance_screen as screen

from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import resemblance as rs
from mls_emergence.probe import spans
from mls_emergence.transmission.spatial import drift_record

OUTD = ROOT / "output" / "probe_latent_groups"
FIND = ROOT / "output" / "findings" / "resemblance_stage2.md"
SEED_BASE = 2026101670
N_SEEDS = 50
FST_WINDOW = 0.50
INHERITED = ("exponential-24km", 2000, 0.001, 0.02)
G: dict = {}


def run_candidate(args):
    """Same cell-level time-scale measurement as the screen, then 50 seeds at those scales."""
    c, family, n_ind, innovation, mixing = args
    ss = np.random.SeedSequence(SEED_BASE + c); kids = ss.spawn(N_SEEDS + 3 + 5)
    w = G["weights"][family]
    sc = screen.cell_scales(w, n_ind, innovation, mixing, kids[N_SEEDS + 3:])
    L_d, S, burnin, scales_from = sc["L_d"], sc["S"], sc["burnin"], ("own, accepted" if sc["accepted"] else "own, NOT accepted")
    rng_len, rng_end, rng_sub = (np.random.default_rng(k) for k in kids[N_SEEDS:N_SEEDS + 3])
    rows = []
    for j in range(N_SEEDS):
        rec = drift_record(w, k=G["K"], n_ind=n_ind, innovation=innovation, mixing=mixing, steps=S, n_records=S,
                           burnin=burnin, seed=int(np.random.default_rng(kids[j]).integers(0, 2**31 - 1)), initial="uniform", target=G["pooled"])
        y = spans.sample_spans(rec, spans.ends_t1(rng_end, G["N"], S, L_d), spans.draw_lengths(rng_len, G["N"], L_d), G["totals"], rng_sub)
        rows.append(rs.summaries(y, G["labels"], G["diversity"]))
    df = pd.DataFrame(rows); obs = G["observed"]
    mean = {k: float(df[k].mean()) for k in ("hs", "rich", "ht")}
    se = {k: float(df[k].std(ddof=1) / np.sqrt(len(df))) for k in ("hs", "rich", "ht")}
    fst = df["fst_clusters"].to_numpy(float); fst = fst[np.isfinite(fst)]
    cell = rs.screen_cell(mean, obs)
    fst_rel = abs(float(np.median(fst)) - obs["fst_clusters"]) / obs["fst_clusters"]
    return {"candidate": c, "family": family, "n_ind": n_ind, "innovation": innovation, "mixing": mixing, "n_seeds": N_SEEDS,
            "L_d": L_d, "S": S, "burnin": burnin, "time_scales_from": scales_from,
            **{f"mean_{k}": v for k, v in mean.items()}, **{f"se_{k}": v for k, v in se.items()},
            "fst_median": float(np.median(fst)), "fst_q10": float(np.quantile(fst, .1)), "fst_q90": float(np.quantile(fst, .9)),
            "fst_frac_ge_observed": float(np.mean(fst >= obs["fst_clusters"])), "fst_rel": fst_rel,
            **cell, "loss4": cell["loss"] + fst_rel ** 2}


def main() -> int:
    m47 = importlib.import_module("47_revision_analysis")
    from a0_predictive_summaries import spatial_clusters
    counts, coords = load_curated(); xy = coords.to_numpy(float); D = geo_km(xy)
    G.update(N=len(counts), K=counts.shape[1], weights={gen.family_name(k, l): gen.family_weights(k, l, D) for k, l in gen.FAMILIES},
             pooled=counts.to_numpy(float).sum(0) / counts.to_numpy(float).sum(), totals=counts.to_numpy(int).sum(1),
             labels=spatial_clusters(xy), diversity=m47.diversity)
    screen.G.update(G)
    G["observed"] = rs.summaries(counts.to_numpy(int), G["labels"], m47.diversity); obs = G["observed"]
    screen_csv = pd.read_csv(OUTD / "resemblance_screen.csv")
    screen_tab = screen_csv
    cand = screen_tab[screen_tab.matched & screen_tab.scales_accepted
                      & ((screen_tab.fst_median - obs["fst_clusters"]).abs() / obs["fst_clusters"] <= FST_WINDOW)]
    cells = [(r.family, int(r.n_ind), float(r.innovation), float(r.mixing)) for r in cand.itertuples()]
    if INHERITED not in cells:
        cells.append(INHERITED)
    args = [(c, *cell) for c, cell in enumerate(cells)]
    t0 = time.time()
    with mp.get_context("fork").Pool(min(28, mp.cpu_count())) as pool:
        rows = pool.map(run_candidate, args)
    secs = time.time() - t0
    tab = pd.DataFrame(rows).sort_values("loss4"); tab.to_csv(OUTD / "resemblance_stage2.csv", index=False)
    L = ["# Stage-2 recheck of candidate operating points, 50 fresh seeds each (D21)", "",
         ("Produced by `analyses/probe_latent_groups/resemblance_stage2.py` after the widened screen; same functions and spans as the "
          "screen. Exploratory design validation, not evidence about the record. No selection is made here."), "",
         (f"Candidates: cells matched on the three diversity summaries, with an accepted time-scale measurement, and with F_ST median within {FST_WINDOW:.0%} of the observed value in "
          f"the 6-seed screen ({len(cand)}), plus the inherited operating point. {len(cells)} cells, {N_SEEDS} seeds each, {secs:.0f} s. "
          f"Observed: hs {obs['hs']:.4f}, rich {obs['rich']:.3f}, ht {obs['ht']:.4f}, two-cluster F_ST {obs['fst_clusters']:.4f}. "
          "Matched = three diversity means within 10 percent; loss4 adds the squared relative error of the F_ST median."), "",
         "| family | n_ind | innovation | mixing | L_d, S (from) | hs (se) | rich (se) | ht (se) | matched | F_ST median [10%, 90%] | seeds at or above observed | loss (3) | loss4 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in tab.itertuples():
        L.append(f"| {r.family} | {int(r.n_ind)} | {r.innovation:g} | {r.mixing:g} | {r.L_d}, {r.S} ({r.time_scales_from}) | {r.mean_hs:.4f} ({r.se_hs:.4f}) | {r.mean_rich:.3f} ({r.se_rich:.3f}) | "
                 f"{r.mean_ht:.4f} ({r.se_ht:.4f}) | {'yes' if r.matched else 'no'} | {r.fst_median:.4f} [{r.fst_q10:.4f}, {r.fst_q90:.4f}] | "
                 f"{r.fst_frac_ge_observed:.2f} | {r.loss:.4f} | {r.loss4:.4f} |")
    L += ["", "Full table: `output/probe_latent_groups/resemblance_stage2.csv`."]
    FIND.write_text("\n".join(L) + "\n"); print("\n".join(L[4:])); print(f"\nwrote {FIND}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
