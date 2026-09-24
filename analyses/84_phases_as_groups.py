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
# Parkin-specific scores (2026-09-23, author framing: does the Parkin phase,
# treated as a polity, show anything in its type frequencies beyond drift?).
# Reported per cell; not used in the ABC so the posterior is unchanged.
PARKIN = ["parkin_fst", "be_parkin"]
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


def _pct(x: float) -> str:
    """A share as a percentage; below 1 percent keep one decimal so 1 of 300 is not printed as 0%."""
    return f"{x:.1%}" if 0 < x < 0.01 else f"{x:.0%}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--placebo", type=int, default=None,
                    help="replace the phases with a same-size division around random centers drawn with this seed "
                         "(a placebo: if it gives the same posterior, the result is a property of the map)")
    args = ap.parse_args()
    if args.reps > 1000:
        raise ValueError("--reps above 1000 would collide seed families between cells")

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
    global OUT_MD, OUT_CSV
    placebo_note = ""
    if args.placebo is not None:
        t74 = importlib.import_module("74_phase_partition_test")
        pts = t74.km_xy(coords.to_numpy(float))
        sizes = np.bincount(phase); k = len(sizes); slot = np.repeat(np.arange(k), sizes)
        rngp = np.random.default_rng(args.placebo)
        a = t74.assign_exact(pts, pts[rngp.choice(len(pts), size=k, replace=False)], slot)
        placebo = a  # around random centers, not made compact: varied placebos
        placebo_note = (f"PLACEBO {args.placebo}: the groups are NOT the phases but a division with the phases' sizes "
                        f"around random centers (adjusted Rand with the phases {t74.ari(phase, placebo):.3f}); "
                        f"'phase' below means this division.")
        print(placebo_note)
        phase = placebo
        OUT_MD = OUT_MD.with_name(f"phases_as_groups_placebo{args.placebo}.md")
        OUT_CSV = OUT_CSV.with_name(f"phases_as_groups_runs_placebo{args.placebo}.csv")
    m_obs = data["m"]; d = data["d"]; totals = m_obs.sum(1)

    rates_all = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates_all[(rates_all.region == "basin") & (rates_all.model == "pooled")].iloc[0]
    cell = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))

    parkin = np.array([1 if l == "Parkin" else 0 for l in labels_ph])

    def summarize(m):
        div = rev.diversity(m)
        return dict(phase_fst=float(rev.fst_by(m, phase)),
                    spatial_fst=float(rev.fst_by(m, data["labels"])),
                    be_phase=float(sd.boundary_excess_labeled(m, d, phase)),
                    parkin_fst=float(rev.fst_by(m, parkin)),
                    be_parkin=float(sd.boundary_excess_labeled(m, d, parkin)),
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
                # separate streams for the simulator and for innovation
                # profiles and sampling (they shared one seed until 2026-09-23)
                rng = np.random.default_rng(seed + 5_000_000)
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
            pk_med=g.parkin_fst.median(), pk_lo=g.parkin_fst.quantile(.025), pk_hi=g.parkin_fst.quantile(.975),
            pk_reach=float((g.parkin_fst >= obs["parkin_fst"]).mean()),
            bp_med=g.be_parkin.median(), bp_lo=g.be_parkin.quantile(.025), bp_hi=g.be_parkin.quantile(.975),
            bp_reach=float((g.be_parkin >= obs["be_parkin"]).mean()),
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

    # How much does the leak factor actually cut cross-phase copying? Weights
    # are renormalized per site, so a factor of 0.03 does not cut the share of
    # a site's between-site copying that crosses a phase line to 3 percent.
    cross = {}
    for leak in LEAKS:
        w = copying_weights(d, 24.0, labels=phase, leak=leak)
        off = w * (1 - np.eye(len(w)))
        share = np.array([off[i, phase != phase[i]].sum() / off[i].sum() for i in range(len(w))])
        cross[leak] = (float(share.mean()), float(share.max()))

    # Sensitivity of the posterior to the approximation's own settings.
    def p_bound(target_vec, acc, cols):
        Xs = df[cols].to_numpy(float); sds = np.nanstd(Xs, axis=0); sds[sds == 0] = 1.0
        tv = np.array([obs[c] for c in cols], float)
        dist = np.sqrt(np.nansum(((Xs - tv) / sds) ** 2, axis=1))
        sub = df.iloc[np.argsort(dist)[:max(1, int(acc * len(df)))]]
        return float((sub.leak < 1).mean()), float((sub.strength > 0).mean())
    sens = []
    for acc in (0.01, 0.02, 0.05, 0.10):
        for label, cols in [("all six summaries", SUMMARIES), ("phase F_ST and boundary excess only", ["phase_fst", "be_phase"])]:
            pb, pl = p_bound(obs_vec, acc, cols)
            sens.append((acc, label, pb, pl))
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

    L = ["# Do the published phases behave as groups?" + (f" (placebo {args.placebo})" if args.placebo is not None else ""), "",
         *([placebo_note, ""] if args.placebo is not None else []),
         f"Produced by `analyses/84_phases_as_groups.py`, {args.reps} runs per cell, calibrated "
         f"pooled-profile cell ({cell['n_ind']} learners, innovation {cell['innovation']}, mixing "
         f"{cell['mixing']}), interaction length 24 km along rivers. Groups are " + ("the placebo division" if args.placebo is not None else "the published phases") + " "
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
                 f"{_pct(c.fst_reach)} | {c.be_med:+.1f} [{c.be_lo:+.1f}, {c.be_hi:+.1f}] | {_pct(c.be_reach)} | "
                 f"{'yes' if c.matched else 'no'} |")
    base = cdf[(cdf.leak == 1.0) & (cdf.strength == 0.0)].iloc[0]
    L += ["", "## The Parkin phase against the rest of the basin", "",
          f"Observed: Parkin-versus-rest cultural F_ST **{obs['parkin_fst']:.4f}**; boundary excess at the "
          f"Parkin line **{obs['be_parkin']:+.1f}**. The shares are of 300 runs per cell reaching the observed value.", "",
          "| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |",
          "|---|---|---|---|---|---|"]
    for _, c in cdf.iterrows():
        L.append(f"| {c.leak:g} | {c.strength:g} | {c.pk_med:.4f} [{c.pk_lo:.4f}, {c.pk_hi:.4f}] | {_pct(c.pk_reach)} | "
                 f"{c.bp_med:+.1f} [{c.bp_lo:+.1f}, {c.bp_hi:+.1f}] | {_pct(c.bp_reach)} |")
    L += ["", "## 2. The phases under neutral copying alone (leak 1, regional pool)", "",
          f"Between-phase F_ST: observed {obs['phase_fst']:.4f} against a median of {base.fst_med:.4f} "
          f"(95 percent range {base.fst_lo:.4f} to {base.fst_hi:.4f}); {_pct(base.fst_reach)} of runs reach it. "
          f"Boundary excess at phase lines: observed {obs['be_phase']:+.1f} against {base.be_med:+.1f} "
          f"({base.be_lo:+.1f} to {base.be_hi:+.1f}); {_pct(base.be_reach)} of runs reach it.", "",
          "## 3. Posterior over the grid (rejection ABC, uniform prior over cells)", "",
          f"Accepted the closest {ACCEPT:.0%} of {len(df)} runs on standardized "
          f"{', '.join(SUMMARIES)}.", "",
          f"- P(some copying boundary at phase lines, leak < 1) = **{p_boundary:.2f}** (prior 0.75).",
          f"- P(local innovation, share > 0) = **{p_local:.2f}** (prior 0.75).", "",
          "Sensitivity to the approximation's settings (prior 0.75 for both):", "",
          "| acceptance | summaries | P(copying boundary) | P(local innovation) |", "|---|---|---|---|"]
    for acc, label, pb, pl in sens:
        L.append(f"| {acc:.0%} | {label} | {pb:.2f} | {pl:.2f} |")
    L += ["", "What the copying factor does to cross-phase copying (share of each site's between-site "
          "copying that crosses a phase line, mean and maximum over sites):", "",
          "| copying factor | mean share | maximum share |", "|---|---|---|"]
    for leak, (mn, mx) in cross.items():
        L.append(f"| {leak:g} | {mn:.3f} | {mx:.3f} |")
    L += ["",
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
