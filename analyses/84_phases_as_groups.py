"""84_phases_as_groups.py - do the published phases behave as groups?

Two questions the earlier comparisons left open (2026-09-23).

1. THE PHASES THEMSELVES. Every earlier drift comparison scored between-group
   difference on k-means spatial clusters. The hypothesis under test is about
   the published phase assignments, so here the calibrated neutral copying
   model is scored on the phase labels: the between-phase cultural F_ST and
   the boundary excess at phase lines.

2. A FAIR GROUPS MODEL. The bounded-groups model of analysis 47 draws new
   variants from one regional pool, which no boundary interrupts, so it is
   handicapped against the one account that closes the finer-scale gap, local
   innovation (analysis 65, departure C). Here both ingredients are crossed,
   with the PHASES as the groups: copying across phase lines is multiplied by
   a factor (1 = no copying boundary) and each phase draws new variants from a
   partly reordered regional profile (share reordered 0 = regional pool). A
   bounded group in the full sense restricts copying AND innovates locally.

INFERENCE (rule 18). No test. The grid cells get a uniform prior and the
posterior over cells is estimated by rejection approximate Bayesian
computation on standardized summaries (phase F_ST, spatial F_ST at the
silhouette-selected clusters, boundary excess at phase lines, and the three
diversity summaries). The prior puts 3/4 of its mass on some copying boundary,
which leans toward the groups hypothesis, the conservative direction for a
paper that finds none (rule 20). Recovery against the prior (rule 20c):
pseudo-observations drawn from boundary cells are run through the same
posterior to show a boundary is recoverable when one is there.

Outputs: output/findings/phases_as_groups.md, output/phases_as_groups_runs.csv.
Usage: .venv/bin/python analyses/84_phases_as_groups.py [--reps 300]
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

OUT_MD = ROOT / "output" / "findings" / "phases_as_groups.md"
OUT_CSV = ROOT / "output" / "phases_as_groups_runs.csv"
LEAKS = [1.0, 0.5, 0.1, 0.03]
STRENGTHS = [0.0, 0.1, 0.2, 0.3]
SUMMARIES = ["phase_fst", "spatial_fst", "be_phase", "hs", "rich", "ht"]
ACCEPT = 0.05          # rejection-ABC acceptance fraction of all runs
TOL = 0.10             # the calibration's diversity tolerance


def group_targets_by(pooled, labels, strength, rng):
    """Per-node innovation profiles: a partial permutation of the regional
    profile per group, exactly as analysis 65 builds them (so richness and
    within-profile diversity are held fixed and only which classes are
    favored changes)."""
    pooled = np.asarray(pooled, float)
    if strength <= 0:
        return pooled
    k = pooled.size
    out = np.empty((len(labels), k))
    for g in np.unique(labels):
        prof = pooled.copy()
        idx = np.flatnonzero(rng.random(k) < strength)
        if idx.size > 1:
            prof[idx] = prof[rng.permutation(idx)]
        out[labels == g] = prof
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()

    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    sd = importlib.import_module("23_phases_vs_spatial_drift")
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record

    data = rev.load_sets()["basin"]
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    if names != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    labels_ph, _ = ph.assign_phases_by_territory(names, coords.to_numpy(float))
    phases = sorted(set(labels_ph))
    phase = np.array([phases.index(l) for l in labels_ph])
    m_obs = data["m"]; d = data["d"]; totals = m_obs.sum(1)

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == "pooled")].iloc[0]
    cell = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))

    def summarize(m):
        div = rev.diversity(m)
        return dict(phase_fst=float(rev.fst_by(m, phase)),
                    spatial_fst=float(rev.fst_by(m, data["labels"])),
                    be_phase=float(sd.boundary_excess_labeled(m, d, phase)),
                    **div)

    obs = summarize(m_obs)
    print(f"observed: phase F_ST {obs['phase_fst']:.4f}, spatial F_ST {obs['spatial_fst']:.4f}, "
          f"phase boundary excess {obs['be_phase']:+.1f}; phases {dict(zip(phases, np.bincount(phase)))}")

    rows = []
    for li, leak in enumerate(LEAKS):
        for si, strength in enumerate(STRENGTHS):
            w = copying_weights(d, 24.0, labels=phase, leak=leak)
            for r in range(args.reps):
                seed = 84000 + 10000 * li + 1000 * si + r
                rng = np.random.default_rng(seed)
                tgt = group_targets_by(data["pooled"], phase, strength, rng)
                rec = drift_record(w, k=m_obs.shape[1], n_ind=cell["n_ind"], seed=seed,
                                   innovation=cell["innovation"], mixing=cell["mixing"],
                                   burnin=1200, initial="uniform", target=tgt)
                m = sample_record(rec, data["ranks"], totals, rng)
                rows.append(dict(leak=leak, strength=strength, rep=r, **summarize(m)))
            print(f"  leak {leak}, local innovation {strength}: done", flush=True)
    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    # ---- per-cell summaries
    def matched(g):
        return all(abs(g[k].median() - obs[k]) <= TOL * abs(obs[k]) for k in ("hs", "rich", "ht"))
    cells = []
    for (leak, strength), g in df.groupby(["leak", "strength"]):
        cells.append(dict(
            leak=leak, strength=strength,
            fst_med=g.phase_fst.median(), fst_lo=g.phase_fst.quantile(.025), fst_hi=g.phase_fst.quantile(.975),
            fst_reach=float((g.phase_fst >= obs["phase_fst"]).mean()),
            be_med=g.be_phase.median(), be_lo=g.be_phase.quantile(.025), be_hi=g.be_phase.quantile(.975),
            be_reach=float((g.be_phase >= obs["be_phase"]).mean()),
            sp_med=g.spatial_fst.median(),
            matched=matched(g)))
    cdf = pd.DataFrame(cells)

    # ---- rejection ABC posterior over cells (uniform prior)
    X = df[SUMMARIES].to_numpy(float)
    sdv = np.nanstd(X, axis=0); sdv[sdv == 0] = 1.0

    def posterior(target_vec, exclude=None):
        dist = np.sqrt(np.nansum(((X - target_vec) / sdv) ** 2, axis=1))
        if exclude is not None:
            dist[exclude] = np.inf
        n_acc = max(1, int(ACCEPT * len(df)))
        acc = np.argsort(dist)[:n_acc]
        sub = df.iloc[acc]
        post = sub.groupby(["leak", "strength"]).size() / n_acc
        return post.reindex(pd.MultiIndex.from_product([LEAKS, STRENGTHS])).fillna(0.0)

    obs_vec = np.array([obs[k] for k in SUMMARIES], float)
    post = posterior(obs_vec)
    p_boundary = float(post[[l for l in LEAKS if l < 1]].sum())
    p_local = float(post.loc[(slice(None), [s for s in STRENGTHS if s > 0])].sum())
    by_leak = post.groupby(level=0).sum()
    by_strength = post.groupby(level=1).sum()

    # ---- recovery (rule 20c): pseudo-observations from boundary cells
    rng = np.random.default_rng(84999)
    rec_rows = []
    for leak in [0.1, 0.03]:
        for strength in [0.0, 0.2]:
            idx = df.index[(df.leak == leak) & (df.strength == strength)].to_numpy()
            pb = []
            for i in rng.choice(idx, size=min(40, len(idx)), replace=False):
                pp = posterior(X[i], exclude=np.array([i]))
                pb.append(float(pp[[l for l in LEAKS if l < 1]].sum()))
            rec_rows.append(dict(leak=leak, strength=strength, p_boundary_median=float(np.median(pb)),
                                 p_boundary_lo=float(np.percentile(pb, 10))))
    idx0 = df.index[(df.leak == 1.0) & (df.strength == 0.0)].to_numpy()
    pb0 = [float(posterior(X[i], exclude=np.array([i]))[[l for l in LEAKS if l < 1]].sum())
           for i in rng.choice(idx0, size=min(40, len(idx0)), replace=False)]

    L = ["# Do the published phases behave as groups?", "",
         f"Produced by `analyses/84_phases_as_groups.py`, {args.reps} runs per cell, calibrated "
         f"pooled-profile cell ({cell['n_ind']} learners, innovation {cell['innovation']}, mixing "
         f"{cell['mixing']}), interaction length 24 km along rivers. Groups are the published phases "
         f"({', '.join(f'{p} {n}' for p, n in zip(phases, np.bincount(phase)))}).", "",
         f"Observed: between-phase cultural F_ST **{obs['phase_fst']:.4f}**; boundary excess at phase "
         f"lines **{obs['be_phase']:+.1f}**; F_ST at the {len(np.unique(data['labels']))} spatial clusters "
         f"{obs['spatial_fst']:.4f}.", "",
         "## 1. Every cell", "",
         "Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the "
         "share of classes reordered in each phase's source of new variants (0 = one regional pool).", "",
         "| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |",
         "|---|---|---|---|---|---|---|"]
    for _, c in cdf.iterrows():
        L.append(f"| {c.leak:g} | {c.strength:g} | {c.fst_med:.4f} [{c.fst_lo:.4f}, {c.fst_hi:.4f}] | "
                 f"{c.fst_reach:.0%} | {c.be_med:+.1f} [{c.be_lo:+.1f}, {c.be_hi:+.1f}] | {c.be_reach:.0%} | "
                 f"{'yes' if c.matched else 'no'} |")
    base = cdf[(cdf.leak == 1.0) & (cdf.strength == 0.0)].iloc[0]
    L += ["", "## 2. The phases under neutral copying alone (leak 1, regional pool)", "",
          f"Between-phase F_ST: observed {obs['phase_fst']:.4f} against a median of {base.fst_med:.4f} "
          f"(95 percent range {base.fst_lo:.4f} to {base.fst_hi:.4f}); {base.fst_reach:.0%} of runs reach it. "
          f"Boundary excess at phase lines: observed {obs['be_phase']:+.1f} against {base.be_med:+.1f} "
          f"({base.be_lo:+.1f} to {base.be_hi:+.1f}); {base.be_reach:.0%} of runs reach it.", "",
          "## 3. Posterior over the grid (rejection ABC, uniform prior over cells)", "",
          f"Accepted the closest {ACCEPT:.0%} of {len(df)} runs on standardized "
          f"{', '.join(SUMMARIES)}.", "",
          f"- P(some copying boundary at phase lines, leak < 1) = **{p_boundary:.2f}** (prior 0.75).",
          f"- P(local innovation, share > 0) = **{p_local:.2f}** (prior 0.75).", "",
          "| leak | posterior | | local innovation | posterior |", "|---|---|---|---|---|"]
    for (lk, pv), (st, sv) in zip(by_leak.items(), by_strength.items()):
        L.append(f"| {lk:g} | {pv:.2f} | | {st:g} | {sv:.2f} |")
    L += ["", "## 4. Recovery (rule 20c)", "",
          "Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior "
          "(each excluded from its own reference table). If the design can see a boundary, "
          "P(leak < 1) should rise well above the prior of 0.75.", "",
          "| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |", "|---|---|---|---|"]
    for r in rec_rows:
        L.append(f"| {r['leak']:g} | {r['strength']:g} | {r['p_boundary_median']:.2f} | {r['p_boundary_lo']:.2f} |")
    L += [f"| 1 (no boundary) | 0 | {np.median(pb0):.2f} | {np.percentile(pb0, 10):.2f} |", "",
          "## Reading", "",
          "Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, "
          "or near the no-boundary row, says the data do not ask for a copying boundary at the phase "
          "lines once local innovation is available; a value near the boundary rows says they do.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
