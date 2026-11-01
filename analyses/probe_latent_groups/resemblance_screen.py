"""The restricted resemblance screen under the probe's pipeline (spec 3.3; R2-6).

Latent-group probe, pilot step 5. Exploratory design validation (rule 1).

Script 47's stage-1 grid (3 population sizes by 11 innovation rates by 9
mixing rates from its CONFIG; pooled-profile innovation only; 6 seeds per
cell) is re-run with straight-line distance, every step stored, and
paired-span T1 sampling at the measured L_d (time_scales.json), for every
copying-network family of mls_emergence.probe.generators (D21, 2026-10-16:
the inherited 24 km exponential kernel is re-screened rather than accepted;
exponential, Gaussian and power-law kernels at 6, 12 and 24 km, and
Gabriel-neighbor copying with no length). Each cell first measures its own
L_d and L_R (5 seeds, stored length 5,000 after a burn-in of 5,000, at most
one doubling; mls_emergence.probe.timescales.cell_time_scales) and draws its
spans and stored length from them, because the time scales depend on the
rates (found 2026-10-16: a single L_d for every cell made the screen
inconsistent with the chosen point's own spans). Each cell's
burn-in is the larger of the probe's burn-in and 5 over the cell's own
innovation rate (declared: the regional time scale of a cell is about one
over its innovation rate). The four summaries are the screen's own
(47.diversity; plug-in cultural_fst on the step-4 spatial clusters). Means
across seeds are compared with the observed values at the screen's 10
percent tolerance; the S0-stress cell is the matched cell with the smallest
F_ST median, ties by loss then grid order; the stop rule is applied to the
operating point in the same table.

Seeds: cell c (global grid order over families and rates) uses
SeedSequence(2026101650 + c) with children run (per seed j: child j),
lengths, ends, subsample, and five time-scale runs, so no stream depends on
another's consumption.

Usage: .venv/bin/python analyses/probe_latent_groups/resemblance_screen.py
Outputs: output/probe_latent_groups/resemblance_screen.csv,
  output/findings/resemblance_screen.md
"""
from __future__ import annotations

import importlib
import itertools
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses")); sys.path.insert(0, str(Path(__file__).resolve().parent))

from mls_emergence.dataio.basin import data_fingerprint, geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import resemblance as rs
from mls_emergence.probe import spans
from mls_emergence.probe import timescales as ts
from mls_emergence.transmission.spatial import drift_record

OUTD = ROOT / "output" / "probe_latent_groups"
FIND = ROOT / "output" / "findings" / "resemblance_screen.md"
OPERATING_POINT = (2000, 0.001, 0.02)
INHERITED_FAMILY = "exponential-24km"
SEED_BASE = 2026101650
N_SEEDS = 6

G: dict = {}   # shared with forked workers


def cell_scales(w, n_ind, innovation, mixing, seed_kids):
    def run(seed, S):
        return drift_record(w, k=G["K"], n_ind=n_ind, innovation=innovation, mixing=mixing, steps=S, n_records=S,
                            burnin=int(max(5000, np.ceil(5.0 / innovation))), seed=seed, initial="uniform", target=G["pooled"])
    seeds = [int(np.random.default_rng(k).integers(0, 2**31 - 1)) for k in seed_kids]
    return ts.cell_time_scales(run, seeds, innovation=innovation)


def run_cell(args):
    c, family, n_ind, innovation, mixing = args
    w = G["weights"][family]
    ss = np.random.SeedSequence(SEED_BASE + c)
    kids = ss.spawn(N_SEEDS + 3 + 5)
    sc = cell_scales(w, n_ind, innovation, mixing, kids[N_SEEDS + 3:])
    S, L_d, burnin = sc["S"], sc["L_d"], sc["burnin"]
    rng_len, rng_end, rng_sub = (np.random.default_rng(k) for k in kids[N_SEEDS:N_SEEDS + 3])
    rows = []
    for j in range(N_SEEDS):
        rec = drift_record(w, k=G["K"], n_ind=n_ind, innovation=innovation, mixing=mixing, steps=S, n_records=S,
                           burnin=burnin, seed=int(np.random.default_rng(kids[j]).integers(0, 2**31 - 1)),
                           initial="uniform", target=G["pooled"])
        lengths = spans.draw_lengths(rng_len, G["N"], L_d)
        ends = spans.ends_t1(rng_end, G["N"], S, L_d)
        y = spans.sample_spans(rec, ends, lengths, G["totals"], rng_sub)
        rows.append(rs.summaries(y, G["labels"], G["diversity"]))
    df = pd.DataFrame(rows)
    mean = {k: float(df[k].mean()) for k in ("hs", "rich", "ht")}
    fst = df["fst_clusters"].to_numpy(float); fst = fst[np.isfinite(fst)]
    out = {"grid_index": c, "family": family, "n_ind": n_ind, "innovation": innovation, "mixing": mixing, "burnin": burnin,
           "L_d": L_d, "L_R": sc["L_R"], "S": S, "scales_accepted": sc["accepted"],
           **{f"mean_{k}": v for k, v in mean.items()},
           "fst_median": float(np.median(fst)) if len(fst) else float("nan"), "fst_mean": float(fst.mean()) if len(fst) else float("nan"),
           "fst_q90": float(np.quantile(fst, 0.9)) if len(fst) else float("nan"),
           "fst_frac_ge_observed": float(np.mean(fst >= G["observed"]["fst_clusters"])) if len(fst) else float("nan"),
           "fst_n_finite": int(len(fst))}
    out.update(rs.screen_cell(mean, G["observed"]))
    return out


def main() -> int:
    m47 = importlib.import_module("47_revision_analysis")
    from a0_predictive_summaries import spatial_clusters
    cfg = m47.CONFIG
    counts, coords = load_curated()
    fp = data_fingerprint(counts, coords)
    xy = coords.to_numpy(float)
    D = geo_km(xy)
    G.update(N=len(counts), K=counts.shape[1],
             weights={gen.family_name(k, l): gen.family_weights(k, l, D) for k, l in gen.FAMILIES},
             pooled=counts.to_numpy(float).sum(0) / counts.to_numpy(float).sum(), totals=counts.to_numpy(int).sum(1),
             labels=spatial_clusters(xy), diversity=m47.diversity)
    G["observed"] = rs.summaries(counts.to_numpy(int), G["labels"], m47.diversity)
    families = [gen.family_name(k, l) for k, l in gen.FAMILIES]
    grid = [(c, f, n, i, m) for c, (f, n, i, m) in enumerate(itertools.product(families, cfg["calib_n_ind"], cfg["calib_innovations"], cfg["calib_mixings"]))]
    t0 = time.time()
    with mp.get_context("fork").Pool(min(28, mp.cpu_count())) as pool:
        rows = pool.map(run_cell, grid, chunksize=4)
    secs = time.time() - t0
    rows.sort(key=lambda r: r["grid_index"])
    tab = pd.DataFrame(rows)
    OUTD.mkdir(parents=True, exist_ok=True); tab.to_csv(OUTD / "resemblance_screen.csv", index=False)
    stress = rs.select_stress_cell([r for r in rows if r["family"] == INHERITED_FAMILY])
    op = next(r for r in rows if r["family"] == INHERITED_FAMILY and (r["n_ind"], r["innovation"], r["mixing"]) == OPERATING_POINT)
    verdict = rs.stop_rule({k: op[f"mean_{k}"] for k in ("hs", "rich", "ht")}, G["observed"])
    matched = tab[tab.matched]
    obs = G["observed"]

    def cell_str(r):
        return f"{r['family']}: n_ind {int(r['n_ind'])}, innovation {r['innovation']:g}, mixing {r['mixing']:g}"
    per_family = []
    for f in families:
        tf = tab[tab.family == f]; mf = tf[tf.matched]
        best = mf.sort_values(["fst_median", "loss"], ascending=[False, True]).iloc[0] if len(mf) else None
        per_family.append((f, len(mf), best))
    L = ["# Restricted resemblance screen under the probe's pipeline (pilot step 5)", "",
         ("Produced by `analyses/probe_latent_groups/resemblance_screen.py`; rules in `mls_emergence.probe.resemblance` "
          "(spec 3.3, R2-6). Exploratory design validation, not evidence about the record."), "",
         (f"Grid: {len(families)} copying-network families by {len(cfg['calib_n_ind'])} population sizes by {len(cfg['calib_innovations'])} innovation rates "
          f"by {len(cfg['calib_mixings'])} mixing rates = {len(grid)} cells, {N_SEEDS} seeds each, {len(grid) * N_SEEDS} runs in {secs:.0f} s on "
          f"{min(28, mp.cpu_count())} processes. Each cell at its own measured time scales (L_d from {int(tab.L_d.min())} to {int(tab.L_d.max())}, "
          f"stored length {int(tab.S.min())} to {int(tab.S.max())}; {int((~tab.scales_accepted).sum())} cells with an unaccepted measurement, flagged). "
          f"Observed counts sha256 `{fp['sha256_counts'][:12]}`; spatial clusters as in step 4 ({len(np.unique(G['labels']))} clusters)."), "",
         "| summary | observed |", "|---|---|",
         *[f"| {k} | {obs[k]:.4f} |" for k in ("hs", "rich", "ht", "fst_clusters")], "",
         "## Families: matched cells and the most spatially structured matched cell of each (D21)", "",
         "| family | matched cells | best matched cell by F_ST median | L_d, S | mean hs | mean rich | mean ht | F_ST median | F_ST 90th pct | seeds at or above observed | loss |",
         "|---|---|---|---|---|---|---|---|---|---|---|",
         *[(f"| {f} | {n} | " + (f"n_ind {int(b['n_ind'])}, innovation {b['innovation']:g}, mixing {b['mixing']:g} | {int(b['L_d'])}, {int(b['S'])} | {b['mean_hs']:.4f} | {b['mean_rich']:.3f} | "
                                 f"{b['mean_ht']:.4f} | {b['fst_median']:.4f} | {b['fst_q90']:.4f} | {b['fst_frac_ge_observed']:.2f} | {b['loss']:.4f} |"
                                 if b is not None else "none | | | | | | | | |")) for f, n, b in per_family], "",
         "## Stop rule on the inherited operating point (exponential 24 km; 2,000; 0.001; 0.02)", "",
         (f"Means across {N_SEEDS} seeds: hs {op['mean_hs']:.4f}, rich {op['mean_rich']:.4f}, ht {op['mean_ht']:.4f}; "
          f"relative errors hs {verdict['rel_hs']:.3f}, rich {verdict['rel_rich']:.3f}, ht {verdict['rel_ht']:.3f} at tolerance "
          f"{verdict['tolerance']}; two-cluster F_ST median {op['fst_median']:.4f} against observed {obs['fst_clusters']:.4f} (no tolerance). "
          f"**Verdict: {'PASS' if verdict['pass'] else 'FAIL: the probe stops for discussion'}.**"), "",
         "## Matched cells, all families, ten most structured", "",
         f"{len(matched)} of {len(grid)} cells match all three diversity summaries within 10 percent.", "",
         "| cell | mean hs | mean rich | mean ht | F_ST median | seeds at or above observed | loss |", "|---|---|---|---|---|---|---|",
         *[f"| {cell_str(r)} | {r['mean_hs']:.4f} | {r['mean_rich']:.3f} | {r['mean_ht']:.4f} | {r['fst_median']:.4f} | {r['fst_frac_ge_observed']:.2f} | {r['loss']:.4f} |"
           for _, r in matched.sort_values("fst_median", ascending=False).head(10).iterrows()], "",
         "## S0-stress cell within the inherited family (the stress cell is re-chosen once the operating family is decided)", ""]
    if stress is None:
        L.append("No cell matched: the stress cell cannot be chosen; stop for discussion.")
    else:
        L.append(f"**{cell_str(stress)}** (grid index {stress['grid_index']}): the matched cell with the smallest F_ST median "
                 f"({stress['fst_median']:.4f}; ties by loss {stress['loss']:.4f}, then grid order).")
    L += ["", "Full table: `output/probe_latent_groups/resemblance_screen.csv`."]
    FIND.parent.mkdir(parents=True, exist_ok=True); FIND.write_text("\n".join(L) + "\n")
    print("\n".join(L[6:]))
    (OUTD / "resemblance_screen_summary.json").write_text(json.dumps(
        {"observed": obs, "operating_point": op, "stop_rule": verdict, "stress_cell": stress, "n_matched": int(len(matched)),
         "seconds": secs}, indent=1, default=float))
    print(f"\nwrote {FIND}")
    return 0 if verdict["pass"] and stress is not None else 3


if __name__ == "__main__":
    raise SystemExit(main())
