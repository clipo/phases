"""97_design_to_separate.py - what record would separate local innovation from restricted copying?

Analysis 84 leaves two explanations of the basin's additional variation
standing: new styles arising locally with copying unrestricted, and copying
restricted at the phase lines. Its posterior probability of a restriction is
0.75 against a prior of 0.75, and on pseudo-records it returns 0.79 to 0.94
with a restriction built in and 0.66 without one. The record cannot tell them
apart. This script asks what record could.

THE QUESTION, MADE CONCRETE. Take the two explanations at the settings that
each reach the observed between-phase F_ST about as often (analysis 84):

  restricted copying   copying across the phase lines cut to 3 percent,
                       new variants from one regional pool
  local innovation     copying unrestricted, each phase drawing new variants
                       from a profile with a fifth of its classes reordered

Simulate records under each, put every one through analysis 84's own
procedure (the same sixteen-setting grid, the same six summaries, the same
5 percent rejection ABC, the record left out of its own reference table), and
read the posterior probability of a restriction it returns. A design separates
the two explanations when that probability is high for records made with a
restriction and low for records made without one.

THE DESIGNS. One change at a time from the record as it is, then together:

  sherds     every assemblage four times its observed decorated count
  classes    each of the ten decorated classes split into two or four finer
             classes (shares in geometric proportion, so each split has rarer
             members and the ten classes are recovered by adding them back up).
             Local innovation still reorders the ten classes, so the process
             is the same one read at a finer grain
  collections  the existing Phillips-Ford-Griffin collections that a lower
             sherd minimum would admit (analysis 98's inventory), at their own
             locations, phases and decorated-sherd counts: those with at least
             50, 25 and 10 decorated sherds. This is the design the authors can
             actually reach, since the record will not supply new assemblages
  sites      assemblages added at the late-period sites of the settlement
             compilation that have none (the sites analysis 80 uses), half of
             them and then all, each given the phase of the Phillips area it
             lies in or nearest, the sequence position of its nearest analyzed
             assemblage, and a sherd count drawn from the observed counts

READING THE RESULT. Two descriptive numbers per design (rule 18: no test):
the median posterior probability of a restriction under each explanation, and
the share of (restricted, local-innovation) pairs of records in which the
restricted one gets the higher probability. That share is 0.5 when the
procedure cannot tell the explanations apart and 1 when it always can.

WHAT THIS IS NOT. It is the calibrated model's answer, at the calibrated
combination, for one pair of settings. The copying model is not recalibrated
to the finer classes or the larger site set, so the designs are compared with
each other, not with the record. Finer classes here are a split of the same
frequencies; real attribute classes need not behave that way. And the added
sites are where late-period sites are recorded, not where collections exist.

CHECK ON THE MACHINERY. The first design is the record as it is, run with
analysis 84's seeds; its summaries must equal analysis 84's saved runs, and
the script stops if they do not.

Outputs: output/findings/design_to_separate.md, output/design_to_separate_runs.csv
Usage: PYTHONPATH=src .venv/bin/python analyses/97_design_to_separate.py [--reps 300] [--procs 16]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import importlib
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "design_to_separate.md"
OUT_CSV = ROOT / "output" / "design_to_separate_runs.csv"
RUNS_84 = ROOT / "output" / "phases_as_groups_runs.csv"

COLLECTIONS = ROOT / "output" / "findings" / "small_collections.csv"
COLLECTION_MINIMA = (50, 25, 10)
SHERD_FACTOR = 4
SITE_SEED = 97           # which added sites, their sherd counts
RESTRICTED = (0.03, 0.0)  # (copying factor at the phase lines, share of classes reordered)
LOCAL = (1.0, 0.2)
PHASES = ("Kent", "Parkin", "Walls")

_G: dict = {}   # per-design inputs, set before the worker pool forks


def split_profile(pooled: np.ndarray, parts: int) -> np.ndarray:
    """Split each class into `parts` finer classes with shares 1 : 1/2 : 1/4 ...

    The finer profile sums back to the original class by class, so a design
    with finer classes is the same regional frequencies read at a finer grain,
    with rarer classes among them."""
    pooled = np.asarray(pooled, float)
    if parts < 1:
        raise ValueError("parts must be at least 1")
    if parts == 1:
        return pooled.copy()
    share = 0.5 ** np.arange(parts)
    share = share / share.sum()
    return (pooled[:, None] * share[None, :]).ravel()


def _cell(job):
    """All runs of one (copying factor, local innovation) setting for the current design."""
    li, si, leak, strength, reps, seed_base = job
    g = _G
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record
    w = copying_weights(g["d"], 24.0, labels=g["phase"], leak=leak)
    rows = []
    for r in range(reps):
        seed = seed_base + 10000 * li + 1000 * si + r
        rng = np.random.default_rng(seed + 5_000_000)
        # Local innovation reorders the TEN classes, as in analysis 84, and the
        # reordered profile is then split. Reordering the finer classes instead
        # would change the process with the grain of the classification: most
        # finer classes are rare, so the same share reordered moves far less
        # frequency, and the two explanations stop being comparable.
        tgt = g["a84"].group_targets_by(g["pooled10"], g["phase"], strength, rng)
        tgt = (split_profile(tgt, g["parts"]) if tgt.ndim == 1 else
               np.array([split_profile(t, g["parts"]) for t in tgt]))
        rec = drift_record(w, k=10 * g["parts"], n_ind=g["n_ind"], seed=seed,
                           innovation=g["innovation"], mixing=g["mixing"],
                           burnin=1200, initial="uniform", target=tgt)
        m = sample_record(rec, g["ranks"], g["totals"], rng)
        rows.append(dict(sherds=1, leak=leak, strength=strength, rep=r, **g["summarize"](m)))
        rng4 = np.random.default_rng(seed + 9_000_000)
        m4 = sample_record(rec, g["ranks"], g["totals"] * SHERD_FACTOR, rng4)
        rows.append(dict(sherds=SHERD_FACTOR, leak=leak, strength=strength, rep=r, **g["summarize"](m4)))
    return rows


def p_restriction(df: pd.DataFrame, summaries, accept: float, truth) -> np.ndarray:
    """Posterior P(copying factor < 1) for every run of one setting, each run
    left out of its own reference table (analysis 84's rejection ABC)."""
    X = df[summaries].to_numpy(float)
    sdv = np.nanstd(X, axis=0)
    sdv[sdv == 0] = 1.0
    Z = X / sdv
    leak = df["leak"].to_numpy(float)
    idx = np.flatnonzero((leak == truth[0]) & (df["strength"].to_numpy(float) == truth[1]))
    n_acc = max(1, int(accept * len(df)))
    out = np.empty(len(idx))
    for j, i in enumerate(idx):
        dist = np.sqrt(np.nansum((Z - Z[i]) ** 2, axis=1))
        dist[i] = np.inf
        acc = np.argpartition(dist, n_acc)[:n_acc]
        out[j] = float((leak[acc] < 1).mean())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--procs", type=int, default=16)
    args = ap.parse_args()
    if not 20 <= args.reps <= 1000:
        raise ValueError("--reps must be between 20 and 1000 (seed families collide above 1000)")

    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    mm = importlib.import_module("make_map")
    ph = importlib.import_module("36_canonical_phase_map")
    sd = importlib.import_module("23_phases_vs_spatial_drift")
    a84 = importlib.import_module("84_phases_as_groups")
    a80 = importlib.import_module("80_boundaries_in_settlement_gaps")
    from mls_emergence.dataio.settlement import load_lmv
    from pyproj import Transformer

    data = rev.load_sets()["basin"]
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    if names != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    xy = coords.to_numpy(float)                      # latitude, longitude
    labels_ph, _ = ph.assign_primary_phases(names, xy)
    if sorted(set(labels_ph)) != sorted(PHASES):
        raise RuntimeError(f"expected phases {PHASES}, got {sorted(set(labels_ph))}")
    phase28 = np.array([PHASES.index(l) for l in labels_ph])
    totals28 = data["m"].sum(1)

    rates = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates[(rates.region == "basin") & (rates.model == "pooled")].iloc[0]
    cell = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))

    # ---- the late-period sites without an analyzed assemblage (analysis 80's selection)
    tr = Transformer.from_crs("EPSG:4326", mm.UTM15N, always_xy=True)
    E, N = (np.asarray(v, float) for v in tr.transform(xy[:, 1], xy[:, 0]))
    pad = 18_000.0
    ext = (E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad)
    lmv = load_lmv(ROOT / "data" / "LMVData.xlsx")
    flag = lambda c: pd.to_numeric(lmv[c], errors="coerce").fillna(0) == 1
    late = lmv[(flag("A") | flag("B")) & (pd.to_numeric(lmv["Zone"], errors="coerce") == 15)]
    sx = pd.to_numeric(late["Easting"], errors="coerce").to_numpy(float)
    sy = pd.to_numeric(late["Northing"], errors="coerce").to_numpy(float)
    keep = np.isfinite(sx) & np.isfinite(sy)
    keep &= (sx >= ext[0]) & (sx <= ext[1]) & (sy >= ext[2]) & (sy <= ext[3])
    sx, sy = sx[keep], sy[keep]
    near = np.hypot(sx[:, None] - E[None, :], sy[:, None] - N[None, :]).min(1) <= a80.EXCLUDE_KM * 1000.0
    sx, sy = sx[~near], sy[~near]
    n_independent = len(sx)
    back = Transformer.from_crs(mm.UTM15N, "EPSG:4326", always_xy=True)
    lon, lat = back.transform(sx, sy)
    add_xy = np.column_stack([lat, lon])
    add_lab, _ = ph.assign_primary_phases([f"site{i}" for i in range(len(add_xy))], add_xy)
    in_scheme = np.array([l in PHASES for l in add_lab])
    add_xy, add_lab = add_xy[in_scheme], np.asarray(add_lab)[in_scheme]
    sx, sy = sx[in_scheme], sy[in_scheme]
    n_added_all = len(add_xy)
    rng_sites = np.random.default_rng(SITE_SEED)
    order_added = rng_sites.permutation(n_added_all)
    add_totals = rng_sites.choice(totals28, size=n_added_all, replace=True)
    nearest = np.hypot(sx[:, None] - E[None, :], sy[:, None] - N[None, :]).argmin(1)
    add_ranks = data["ranks"][nearest]
    add_phase = np.array([PHASES.index(l) for l in add_lab])

    # ---- the existing PFG collections a lower sherd minimum would admit (analysis 98)
    if not COLLECTIONS.exists():
        raise FileNotFoundError(f"{COLLECTIONS} is missing; run analyses/98_small_collections.py first")
    coll = pd.read_csv(COLLECTIONS)
    cE, cN = (np.asarray(v, float) for v in tr.transform(coll["lon"].to_numpy(float), coll["lat"].to_numpy(float)))
    coll_ranks = data["ranks"][np.hypot(cE[:, None] - E[None, :], cN[:, None] - N[None, :]).argmin(1)]
    coll_phase = np.array([PHASES.index(l) for l in coll["phase"]])
    coll_xy = coll[["lat", "lon"]].to_numpy(float)
    coll_totals = coll["decorated"].to_numpy(int)

    def site_design(spec):
        """Inputs for the 28 analyzed assemblages plus added ones.

        `spec` is 0 (nothing added), a count of settlement-compilation sites to
        add (stand-ins, with sherd counts drawn from the observed ones), or
        ("collections", minimum): the existing PFG collections with at least
        that many decorated sherds, at their own locations and counts."""
        if spec == 0:
            return dict(d=data["d"], phase=phase28, ranks=data["ranks"], totals=totals28,
                        spatial=data["labels"])
        if isinstance(spec, tuple):
            pick = np.flatnonzero(coll_totals >= spec[1])
            x_add, p_add, r_add, t_add = coll_xy[pick], coll_phase[pick], coll_ranks[pick], coll_totals[pick]
        else:
            pick = order_added[:spec]
            x_add, p_add, r_add, t_add = add_xy[pick], add_phase[pick], add_ranks[pick], add_totals[pick]
        allxy = np.vstack([xy, x_add])
        d = mm.river_distance_matrix(allxy)[0]
        centred = allxy - allxy.mean(0)
        return dict(d=d, phase=np.concatenate([phase28, p_add]),
                    ranks=np.concatenate([data["ranks"], r_add]),
                    totals=np.concatenate([totals28, t_add]),
                    spatial=mf._kmeans_labels(centred, len(np.unique(data["labels"])), seed=7))

    half = n_added_all // 2
    designs = [
        ("the record as it is", 0, 1),
        ("classes split in two", 0, 2),
        ("classes split in four", 0, 4),
        (f"{half} sites added", half, 1),
        (f"{n_added_all} sites added", n_added_all, 1),
        (f"{n_added_all} sites added, classes split in four", n_added_all, 4),
    ]
    for mn in COLLECTION_MINIMA:
        n_c = int((coll_totals >= mn).sum())
        designs.append((f"the {n_c} excluded collections with at least {mn} decorated sherds added",
                        ("collections", mn), 1))
    n_c = int((coll_totals >= COLLECTION_MINIMA[-1]).sum())
    designs.append((f"the {n_c} excluded collections with at least {COLLECTION_MINIMA[-1]} added, "
                    f"classes split in four", ("collections", COLLECTION_MINIMA[-1]), 4))

    summaries = a84.SUMMARIES
    jobs_of = lambda base: [(li, si, leak, s, args.reps, base)
                            for li, leak in enumerate(a84.LEAKS) for si, s in enumerate(a84.STRENGTHS)]
    all_rows, results = [], []
    for di, (label, n_add, parts) in enumerate(designs):
        sdz = site_design(n_add)
        phase, d, spatial = sdz["phase"], sdz["d"], sdz["spatial"]

        def summarize(m, phase=phase, d=d, spatial=spatial):
            return dict(phase_fst=float(rev.fst_by(m, phase)),
                        spatial_fst=float(rev.fst_by(m, spatial)),
                        be_phase=float(sd.boundary_excess_labeled(m, d, phase)),
                        **rev.diversity(m))

        _G.clear()
        _G.update(d=d, phase=phase, ranks=sdz["ranks"], totals=sdz["totals"], a84=a84,
                  pooled10=np.asarray(data["pooled"], float), parts=parts, summarize=summarize, **cell)
        # Design 0 uses analysis 84's seed family so its runs can be checked
        # against 84's saved ones; the others get their own.
        base = 84000 if di == 0 else 970000 + 200000 * di
        with mp.get_context("fork").Pool(args.procs) as pool:
            rows = [r for chunk in pool.map(_cell, jobs_of(base)) for r in chunk]
        df = pd.DataFrame(rows)
        df.insert(0, "design", label)
        all_rows.append(df)

        if di == 0 and RUNS_84.exists():
            ref = pd.read_csv(RUNS_84)
            mine = df[df.sherds == 1].merge(ref, on=["leak", "strength", "rep"], suffixes=("", "_84"))
            if len(mine) == 0:
                raise RuntimeError("no runs in common with analysis 84's saved runs")
            worst = max(float(np.nanmax(np.abs(mine[c] - mine[c + "_84"]))) for c in summaries)
            if worst > 1e-9:
                raise RuntimeError(f"the record-as-it-is design does not reproduce analysis 84's runs "
                                   f"(largest difference {worst:.3g}); the designs below would not be "
                                   f"comparable with the paper's posterior")
            print(f"  check: {len(mine)} runs equal analysis 84's saved runs", flush=True)

        # Quadrupled counts are a design of their own only where the counts are
        # hypothetical; the excluded collections are scored at the counts they have.
        for sherds in ((1,) if isinstance(n_add, tuple) else (1, SHERD_FACTOR)):
            sub = df[df.sherds == sherds].reset_index(drop=True)
            pr = p_restriction(sub, summaries, a84.ACCEPT, RESTRICTED)
            pl = p_restriction(sub, summaries, a84.ACCEPT, LOCAL)
            auc = float((pr[:, None] > pl[None, :]).mean() + 0.5 * (pr[:, None] == pl[None, :]).mean())
            cellmask = lambda t: (sub.leak == t[0]) & (sub.strength == t[1])
            results.append(dict(
                design=label + ("" if sherds == 1 else f", sherds x{SHERD_FACTOR}"),
                sites=len(phase), classes=10 * parts, sherds=sherds,
                p_restricted=np.percentile(pr, [10, 50, 90]), p_local=np.percentile(pl, [10, 50, 90]),
                auc=auc,
                fst_restricted=float(sub[cellmask(RESTRICTED)].phase_fst.median()),
                fst_local=float(sub[cellmask(LOCAL)].phase_fst.median()),
                be_restricted=float(sub[cellmask(RESTRICTED)].be_phase.median()),
                be_local=float(sub[cellmask(LOCAL)].be_phase.median())))
            r = results[-1]
            print(f"  {r['design']}: P(restriction) {r['p_restricted'][1]:.2f} with one, "
                  f"{r['p_local'][1]:.2f} without; ordered correctly in {auc:.2f} of pairs", flush=True)

    runs = pd.concat(all_rows, ignore_index=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    runs.to_csv(OUT_CSV, index=False)

    prior = float(np.mean([l < 1 for l in a84.LEAKS]))
    L = ["# What record would separate local innovation from restricted copying?", "",
         f"Produced by `analyses/97_design_to_separate.py`. {args.reps} runs per setting on analysis 84's "
         f"sixteen-setting grid, at the calibrated combination (N {cell['n_ind']}, innovation "
         f"{cell['innovation']:g}, mixing {cell['mixing']:g}), for each design. Each run made with "
         f"**restricted copying** (copying across the phase lines cut to {RESTRICTED[0]:g}, one regional "
         f"pool) or with **local innovation** (copying unrestricted, {LOCAL[1]:g} of classes reordered by "
         f"phase) is put through analysis 84's rejection ABC ({a84.ACCEPT:.0%} of runs kept, six "
         f"summaries), left out of its own reference table, and the posterior probability of a "
         f"restriction at the phase lines is read off. The prior is {prior:.2f}.", "",
         f"Added sites: {n_independent} late-period sites of the settlement compilation lie more than "
         f"{a80.EXCLUDE_KM:g} km from every analyzed assemblage (analysis 80's selection); {n_added_all} of "
         f"them fall in or nearest the Kent, Parkin or Walls areas and are used. Phases of the "
         f"{28 + n_added_all}-site design: "
         + ", ".join(f"{p} {int((np.concatenate([phase28, add_phase]) == i).sum())}" for i, p in enumerate(PHASES))
         + f". Sherd counts for added sites are drawn from the observed counts (seed {SITE_SEED}). "
         "The excluded collections are those of `output/findings/small_collections.md`, each at its own "
         "location, Phillips phase and decorated-sherd count.", "",
         "| design | sites | classes | P(restriction), made with one: median (10th to 90th) | made without one: "
         "median (10th to 90th) | pairs ordered correctly | between-phase F_ST, medians (with / without) | "
         "boundary excess, medians (with / without) |",
         "|---|---|---|---|---|---|---|---|"]
    for r in results:
        L.append(f"| {r['design']} | {r['sites']} | {r['classes']} | "
                 f"{r['p_restricted'][1]:.2f} ({r['p_restricted'][0]:.2f} to {r['p_restricted'][2]:.2f}) | "
                 f"{r['p_local'][1]:.2f} ({r['p_local'][0]:.2f} to {r['p_local'][2]:.2f}) | {r['auc']:.2f} | "
                 f"{r['fst_restricted']:.4f} / {r['fst_local']:.4f} | "
                 f"{r['be_restricted']:+.1f} / {r['be_local']:+.1f} |")
    base_r = results[0]
    best = max(results, key=lambda r: r["auc"])
    L += ["", "## Reading", "",
          "\"Pairs ordered correctly\" is the share of (restricted, local-innovation) pairs of records in "
          "which the restricted one gets the higher posterior probability of a restriction: 0.5 means the "
          "procedure cannot tell the two explanations apart, 1 that it always can. It is a description of "
          "the simulated records, not a test (rule 18).", "",
          f"On the record as it is, the share is {base_r['auc']:.2f}, and the posterior is "
          f"{base_r['p_restricted'][1]:.2f} for records made with a restriction against "
          f"{base_r['p_local'][1]:.2f} for records made without one (prior {prior:.2f}). The highest "
          f"share among the designs is {best['auc']:.2f}, for \"{best['design']}\".", "",
          "The model is not recalibrated to the finer classes or the larger site set, so compare the "
          "designs with each other and not with the record. Finer classes are a geometric split of the "
          "ten classes' frequencies; real attribute classes need not behave that way. The added sites "
          "are where late-period sites are recorded, which is not where collections exist.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[8:]))
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
