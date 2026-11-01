"""Per-class resemblance of neutral operating points to the record (D24 step 1, first measurement; 2026-10-19).

The resemblance screen (resemblance_screen.py, resemblance_stage2.py)
matched four whole-record summaries. This script asks which neutral
operating points reproduce the record's per-class pattern: site-weighted
method-of-moments F per class (sampling-corrected) and the number of sites
at which each class is present. Site-level F_ST (1 - hs/ht) is added to the
screen's stored summaries without new runs; candidates are then rerun at 50
seeds through the probe's own span sampling, each at its own measured time
scales as the screen does.

Candidates, fixed before running: the two D22 operating points; every
diversity-matched cell of the 6-seed screen with accepted time scales whose
site-level F_ST lies within 25 percent of the observed value; and the D22
primary family and rates at 600, 300 and 150 learners (populations outside
the screen's grid; population is not identified, D22).

Reading rule (Beau, 2026-10-19, stated before the run): a candidate covers a
common class when the record's moment F for that class lies inside the
candidate's 10th to 90th percentile across seeds, and covers a rare class
when the record's number of sites present does. The table reports how many
of the ten classes each candidate covers; it makes no selection.

Exploratory design validation (rule 1), not evidence about the record.
Seeds: candidate c uses SeedSequence(2026101900 + c), children as in stage 2.
Usage: .venv/bin/python analyses/probe_latent_groups/per_class_resemblance.py
Outputs: output/probe_latent_groups/per_class_resemblance.csv (per seed),
  output/probe_latent_groups/screen_site_fst.csv,
  output/findings/per_class_resemblance.md
"""
from __future__ import annotations

import importlib
import json
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
FIND = ROOT / "output" / "findings" / "per_class_resemblance.md"
POINTS = ROOT / "analyses" / "probe_latent_groups" / "operating_points.json"
SEED_BASE = 2026101900
N_SEEDS = 50
SITE_FST_WINDOW = 0.25
EXTRA_POPULATIONS = (600, 300, 150)
COMMON = 3   # the first three classes carry F; the rest are read by occupancy
G: dict = {}


def moment_f(y: np.ndarray) -> np.ndarray:
    """Site-weighted method-of-moments F per class, binomial sampling variance removed; nan where a class is absent everywhere."""
    y = np.asarray(y, float); n = y.sum(1); P = y / n[:, None]; pb = P.mean(0)
    v = P.var(0, ddof=1); samp = np.mean(P * (1 - P) / (n - 1)[:, None], axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (v - samp) / (pb * (1 - pb))


def site_fst(y: np.ndarray) -> float:
    d = G["diversity"](y); return float(1 - d["hs"] / d["ht"])


def per_class_rows(y: np.ndarray) -> dict:
    f = moment_f(y); occ = (np.asarray(y) > 0).sum(0); s = rs.summaries(y, G["labels"], G["diversity"])
    return {**{f"F_{c}": f[c] for c in range(G["K"])}, **{f"occ_{c}": int(occ[c]) for c in range(G["K"])},
            "site_fst": site_fst(y), "fst_clusters": s["fst_clusters"], "rich": s["rich"], "hs": s["hs"], "ht": s["ht"]}


def run_candidate(args):
    c, label, family, n_ind, innovation, mixing = args
    ss = np.random.SeedSequence(SEED_BASE + c); kids = ss.spawn(N_SEEDS + 3 + 5)
    w = G["weights"][family]
    sc = screen.cell_scales(w, n_ind, innovation, mixing, kids[N_SEEDS + 3:])
    L_d, S, burnin = sc["L_d"], sc["S"], sc["burnin"]
    rng_len, rng_end, rng_sub = (np.random.default_rng(k) for k in kids[N_SEEDS:N_SEEDS + 3])
    rows = []
    for j in range(N_SEEDS):
        rec = drift_record(w, k=G["K"], n_ind=n_ind, innovation=innovation, mixing=mixing, steps=S, n_records=S,
                           burnin=burnin, seed=int(np.random.default_rng(kids[j]).integers(0, 2**31 - 1)), initial="uniform", target=G["pooled"])
        y = spans.sample_spans(rec, spans.ends_t1(rng_end, G["N"], S, L_d), spans.draw_lengths(rng_len, G["N"], L_d), G["totals"], rng_sub)
        rows.append({"candidate": c, "label": label, "family": family, "n_ind": n_ind, "innovation": innovation, "mixing": mixing,
                     "seed": j, "L_d": L_d, "S": S, "burnin": burnin, "scales_accepted": sc["accepted"], **per_class_rows(y)})
    return rows


def q(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if len(x) else float("nan")


def main() -> int:
    m47 = importlib.import_module("47_revision_analysis")
    from a0_predictive_summaries import spatial_clusters
    counts, coords = load_curated(); xy = coords.to_numpy(float); D = geo_km(xy); cls = list(counts.columns)
    y_real = counts.to_numpy(int)
    G.update(N=len(counts), K=counts.shape[1], weights={gen.family_name(k, l): gen.family_weights(k, l, D) for k, l in gen.FAMILIES},
             pooled=y_real.astype(float).sum(0) / y_real.sum(), totals=y_real.sum(1), labels=spatial_clusters(xy), diversity=m47.diversity)
    screen.G.update(G)
    G["observed"] = rs.summaries(y_real, G["labels"], m47.diversity)
    real = per_class_rows(y_real); obs_site = real["site_fst"]

    # 1. site-level F_ST added to the stored screen, no runs
    sc = pd.read_csv(OUTD / "resemblance_screen.csv"); sc["site_fst"] = 1 - sc.mean_hs / sc.mean_ht
    sc.to_csv(OUTD / "screen_site_fst.csv", index=False)
    m = sc[sc.matched & sc.scales_accepted]
    near = m[(m.site_fst - obs_site).abs() / obs_site <= SITE_FST_WINDOW].sort_values("site_fst", ascending=False)

    # 2. candidate set
    pts = json.loads(POINTS.read_text())
    cands = []
    for key in ("primary", "secondary"):
        p = pts[key]; cands.append((f"D22 {key}", gen.family_name(p["kernel"], p["length"]), int(p["n_ind"]), float(p["innovation"]), float(p["mixing"])))
    for r in near.itertuples():
        cands.append((f"screen, site F_ST {r.site_fst:.3f}", r.family, int(r.n_ind), float(r.innovation), float(r.mixing)))
    p = pts["primary"]
    for n in EXTRA_POPULATIONS:
        cands.append((f"D22 primary family at {n} learners", gen.family_name(p["kernel"], p["length"]), n, float(p["innovation"]), float(p["mixing"])))
    args = [(c, *cand) for c, cand in enumerate(cands)]
    t0 = time.time()
    with mp.get_context("fork").Pool(min(28, mp.cpu_count(), len(args))) as pool:
        rows = [r for rr in pool.map(run_candidate, args) for r in rr]
    secs = time.time() - t0
    per_seed = pd.DataFrame(rows); per_seed.to_csv(OUTD / "per_class_resemblance.csv", index=False)

    # 3. summary table with the reading rule
    K = G["K"]
    def fmt_f(x): return "absent" if not np.isfinite(x) else f"{x:.3f}"
    hdr = ["candidate", "family", "n_ind", "innov.", "mixing", "L_d, S", "site F_ST", "2-cluster F_ST", "rich"] + \
          [f"F {cls[c]}" for c in range(COMMON)] + [f"sites {cls[c]}" for c in range(COMMON, K)] + ["covered of 10"]
    L = ["# Per-class resemblance of neutral operating points to the record", "",
         (f"Produced by `analyses/probe_latent_groups/per_class_resemblance.py`: {len(cands)} candidates, {N_SEEDS} seeds each through the "
          f"probe's span sampling at each candidate's own time scales, {secs:.0f} s. Cells: median [10th, 90th percentile] across seeds. "
          "F is the site-weighted method-of-moments F per class with binomial sampling variance removed; sites is the number of the 28 "
          "sites at which the class is present. Reading rule (Beau, 2026-10-19): a common class is covered when the record's F lies inside "
          "the candidate's 10th to 90th percentile, a rare class when the record's site count does; the last column counts covered classes "
          "and no selection is made. Exploratory design validation, not evidence about the record."), "",
         (f"Site-level F_ST (1 - hs/ht) of the {len(m)} diversity-matched cells with accepted time scales in the 6-seed screen: "
          f"{m.site_fst.min():.4f} to {m.site_fst.max():.4f}, median {m.site_fst.median():.4f}; observed {obs_site:.4f}; {len(near)} within "
          f"{SITE_FST_WINDOW:.0%} (`output/probe_latent_groups/screen_site_fst.csv`)."), "",
         "| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    rec_row = ["**record**", "", "", "", "", "", f"{real['site_fst']:.4f}", f"{real['fst_clusters']:.4f}", f"{real['rich']:.2f}"] + \
              [fmt_f(real[f"F_{c}"]) for c in range(COMMON)] + [str(real[f"occ_{c}"]) for c in range(COMMON, K)] + [""]
    L.append("| " + " | ".join(rec_row) + " |")
    summary = []
    for c, cand in enumerate(cands):
        d = per_seed[per_seed.candidate == c]; covered = 0; cells = []
        for k in range(COMMON):
            lo, me, hi = q(d[f"F_{k}"], .1), q(d[f"F_{k}"], .5), q(d[f"F_{k}"], .9)
            hit = np.isfinite(lo) and lo <= real[f"F_{k}"] <= hi; covered += int(hit)
            cells.append(f"{fmt_f(me)} [{fmt_f(lo)}, {fmt_f(hi)}]" + (" c" if hit else ""))
        for k in range(COMMON, K):
            lo, me, hi = q(d[f"occ_{k}"], .1), q(d[f"occ_{k}"], .5), q(d[f"occ_{k}"], .9)
            hit = lo <= real[f"occ_{k}"] <= hi; covered += int(hit)
            cells.append(f"{me:.0f} [{lo:.0f}, {hi:.0f}]" + (" c" if hit else ""))
        row = [cand[0], cand[1], str(cand[2]), f"{cand[3]:g}", f"{cand[4]:g}", f"{int(d.L_d.iloc[0])}, {int(d.S.iloc[0])}" + ("" if d.scales_accepted.iloc[0] else " (not accepted)"),
               f"{q(d.site_fst, .5):.4f} [{q(d.site_fst, .1):.4f}, {q(d.site_fst, .9):.4f}]",
               f"{q(d.fst_clusters, .5):.4f} [{q(d.fst_clusters, .1):.4f}, {q(d.fst_clusters, .9):.4f}]",
               f"{q(d.rich, .5):.2f} [{q(d.rich, .1):.2f}, {q(d.rich, .9):.2f}]"] + cells + [str(covered)]
        L.append("| " + " | ".join(row) + " |"); summary.append((cand[0], covered))
    L += ["", "A trailing \"c\" marks a covered class. Per-seed values: `output/probe_latent_groups/per_class_resemblance.csv`.", ""]
    FIND.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {FIND}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
