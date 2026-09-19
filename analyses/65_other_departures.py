"""65_other_departures.py - the three remaining departures from the drift null.

`64_unequal_populations.py` tested the first of the four departures the
manuscript lists and found that it does not close the gap: inequality raises
differentiation but lowers the harmonic mean population, so it moves ALONG the
diversity-against-differentiation trade-off rather than off it. This script tests
the other three, on the same terms.

Every cell reports the same two things, because the standing result is that no
parameterization reproducing the assemblages' diversity reproduces their
differentiation. A cell that reaches the observed F_ST by breaking the diversity
match has not closed the gap.

A. TRANSPORT GEOGRAPHY. The copying kernel runs on river-network distance with an
   exponential decay at a fixed interaction length. Two things are assumed there,
   not measured: that movement follows the mapped channels, and that influence
   decays exponentially. The sweep swaps in straight-line distance and varies the
   interaction length, which together ask whether any reasonable geography
   produces the observed differentiation.

B. UNEQUAL ACCUMULATION SPANS. Every deposit is sampled over the same eight-slice
   trailing window. Deposits do not accumulate over equal stretches of time. A
   short window samples a narrow slice of the sequence and so looks distinct; a
   long one averages more of it away. The sweep draws per-site windows with the
   mean held at eight and rising spread.

C. INNOVATION THAT RESPECTS BOUNDARIES. This is the departure the manuscript
   flags as a limitation of its own comparison. Innovation is drawn from a
   region-wide pooled class profile, an input no boundary interrupts, so BOTH
   models understate what a boundary would do. Here each spatial group innovates
   from its own profile instead.

   The obvious way to do that is circular. Using each group's OBSERVED pooled
   profile as its innovation target feeds the answer in: the groups would differ
   because they were told to differ by exactly the amount they differ. So the
   group profiles are instead random Dirichlet perturbations of the regional
   profile, with a concentration parameter setting how far they drift from it.
   The question is then how strong an innovation boundary must be to produce the
   observed differentiation, which is answerable, rather than whether one fitted
   to the data can reproduce the data, which is not.

   This is the one departure that could move OFF the trade-off, because it raises
   differentiation without raising the drift rate.

Usage: .venv/bin/python analyses/65_other_departures.py [--fast] [--only A|B|C] [--report-only]
       [--regions basin,cmv] [--strengths 0.1,0.12,...] [--tag NAME]
       [--save-draws]
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))

r47 = importlib.import_module("47_revision_analysis")
sd = importlib.import_module("23_phases_vs_spatial_drift")
from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record  # noqa: E402

OUT_MD = ROOT / "output" / "findings" / "other_departures.md"
OUT_JSON = ROOT / "output" / "other_departures.json"

FAST = "--fast" in sys.argv
REPORT_ONLY = "--report-only" in sys.argv
SAVE_DRAWS = "--save-draws" in sys.argv
N_REAL = 60 if FAST else 250
ONLY = None
if "--only" in sys.argv:
    ONLY = sys.argv[sys.argv.index("--only") + 1].upper()
# The St. Francis basin is the analytical unit: comparable deposits and the best
# sample sizes, which is what makes it the place to test a hypothesis about
# phases. The southeast-Missouri set is a COMPARISON, included because Williams
# defined definitive phases there, not as a replication. The 55-assemblage
# valley set is deliberately excluded: it mixes deposits of different ages and
# very different sample sizes, so a differentiation measured across it confounds
# chronology with spatial process and no calibrated spatial model should be
# asked to reproduce it.
REGIONS = ["basin", "cmv"]
TOL = 0.10


def _opt(flag, cast=str):
    """Comma-separated CLI override, e.g. --regions cmv --strengths 0.10,0.11."""
    if flag not in sys.argv:
        return None
    return [cast(x) for x in sys.argv[sys.argv.index(flag) + 1].split(",")]

LENGTHS = [6.0, 12.0, 24.0, 48.0, 96.0]          # A: interaction length, km
WINDOW_SPREADS = [0, 2, 4, 6]                     # B: half-range of per-site window
# C: innovation-boundary strength. Finer below 0.4 because that is where the
# gap closes; above it the model overshoots the observed differentiation.
BOUNDARY_STRENGTHS = [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0]

# The default grid above is spaced for the BASIN, where the transition is
# gradual between 0.1 and 0.4. On the southeast-Missouri set the same grid is
# far too coarse: 0.1 is diversity-matched and 1.8x short, 0.15 already breaks
# the match, and 0.2 overshoots the observation sevenfold. One grid point
# separates "does not close it" from "overshoots", which is not enough to
# conclude that no boundary closes it there. `--strengths` refines that
# interval without disturbing the committed basin grid.
_ov = _opt("--regions")
if _ov:
    REGIONS = _ov
_ov = _opt("--strengths", float)
if _ov:
    BOUNDARY_STRENGTHS = _ov
# Keep a distinct output when a refinement is run, so a partial sweep never
# overwrites the full grid the manuscript's tables are checked against.
_TAG = "--tag" in sys.argv and sys.argv[sys.argv.index("--tag") + 1] or ""
if _TAG:
    OUT_MD = OUT_MD.with_name(f"other_departures_{_TAG}.md")
    OUT_JSON = OUT_JSON.with_name(f"other_departures_{_TAG}.json")


def matched(ach, obs):
    return all(abs(ach[k] - obs[k]) <= TOL * abs(obs[k]) for k in ("hs", "rich", "ht"))


def group_targets(pooled, labels, strength, rng):
    """One innovation profile per node, perturbed from the regional profile.

    ``strength`` 0 gives every node the regional profile, which is the model the
    paper reports. 1 gives each spatial group a strongly idiosyncratic profile.
    The perturbation is random, not read off the observed group differences,
    which is what keeps the test from assuming its own conclusion.
    """
    pooled = np.asarray(pooled, float)
    if strength <= 0:
        return pooled
    # PERMUTE, do not perturb. A first version drew each group's profile as a
    # Dirichlet around the regional one, which made the profiles sparser as well
    # as different: mean class richness fell from 5.31 to 3.48 and the diversity
    # match broke for that reason rather than because of the boundary. A
    # permutation of the regional profile holds the multiset of frequencies
    # exactly, so every group keeps the same richness and the same
    # within-profile diversity and differs only in WHICH classes it favors,
    # which is what an innovation boundary means.
    # A PARTIAL permutation, not a mixture. Blending the regional profile with a
    # permuted copy looked like an interpolation but flattened the profile: the
    # real one is highly skewed, two classes holding nine tenths of all sherds,
    # and the average of a skewed profile and a permutation of itself is much
    # more even. Diversity then came out far too HIGH (hs 0.73 against an
    # observed 0.51) for a reason that has nothing to do with boundaries.
    #
    # Permuting a random subset of the classes keeps the multiset of frequencies
    # exactly, at every strength, so richness and within-profile diversity are
    # held fixed by construction and only the assignment of frequency to class
    # changes. That is the whole content of "innovation respects the boundary".
    k = pooled.size
    out = np.empty((len(labels), k))
    for g in np.unique(labels):
        prof = pooled.copy()
        idx = np.flatnonzero(rng.random(k) < strength)
        if idx.size > 1:
            prof[idx] = prof[rng.permutation(idx)]
        out[labels == g] = prof
    return out


def run(data, cell, *, dist=None, length=24.0, window=8, strength=0.0,
        kernel="exponential", n_real=N_REAL):
    part = data["parkin"] if data["parkin"] is not None else data["labels"]
    d = data["d"] if dist is None else dist
    totals = data["m"].sum(1)
    k = data["m"].shape[1]
    fsts, divs = [], []
    for s in range(n_real):
        rng = np.random.default_rng(70000 + s)
        w = copying_weights(d, length, labels=part, leak=1.0, kernel=kernel)
        tgt = group_targets(data["pooled"], data["labels"], strength, rng)
        rec = drift_record(w, k=k, n_ind=cell["n_ind"], seed=70000 + s,
                           innovation=cell["innovation"], mixing=cell["mixing"],
                           burnin=1200, initial="uniform", target=tgt)
        m = sample_record(rec, data["ranks"], totals, rng, window=window)
        fsts.append(r47.fst_by(m, part))
        divs.append(r47.diversity(m))
    fsts = np.asarray(fsts, float)
    ach = {key: float(np.median([x[key] for x in divs])) for key in ("hs", "rich", "ht")}
    return fsts, ach


def run_region(sets, region, cell_by_region, rows):
    """Every departure, on one region. Rows carry their region so the report
    can separate a result that replicates from one that does not."""
    data = sets[region]
    cell = cell_by_region[region]
    part = data["parkin"] if data["parkin"] is not None else data["labels"]
    obs_fst = float(r47.fst_by(data["m"], part))
    obs_div = r47.diversity(data["m"])
    straight = sd.geo_km(data["coords"])
    n_sites = len(data["m"])
    print(f"{region}: observed F_ST {obs_fst:.4f}; calibrated N {cell['n_ind']}, "
          f"innovation {cell['innovation']}, mixing {cell['mixing']}\n")


    def record(dep, label, fsts, ach):
        # `--save-draws` keeps every realization, not just the quantiles. The
        # summaries answer "where does the model sit"; only the draws support a
        # posterior over the departure's strength given the observation, which
        # is what turns "the median passes through the observed value" into a
        # statement with an interval on it (rule 18).
        r = dict(departure=dep, setting=label, observed=obs_fst,
                 draws=[float(x) for x in fsts] if SAVE_DRAWS else None,
                 fst_median=float(np.median(fsts)),
                 fst_lo=float(np.percentile(fsts, 2.5)),
                 fst_hi=float(np.percentile(fsts, 97.5)),
                 reach=float(np.mean(fsts >= obs_fst)),
                 diversity_matched=bool(matched(ach, obs_div)), achieved=ach)
        rows.append(dict(r, region=region))
        print(f"  {dep} {label:<28} F_ST {r['fst_median']:.4f} "
              f"[{r['fst_lo']:.4f}, {r['fst_hi']:.4f}]  shortfall "
              f"{obs_fst / max(r['fst_median'], 1e-9):>5.1f}x  "
              f"diversity {'matched' if r['diversity_matched'] else 'BROKEN'}")

    if ONLY in (None, "A"):
        print("A. transport geography")
        for dist, dname in ((None, "river"), (straight, "straight-line")):
            for L in LENGTHS:
                f, a = run(data, cell, dist=dist, length=L)
                record("A", f"{dname}, exponential, {L:g} km", f, a)
        # Vary the FORM of the decay, not only its scale. The exponential is an
        # assumption the model made silently; a Gaussian falls away faster at
        # range and isolates distant nodes about eleven times more, which is the
        # shape most likely to raise differentiation through geography alone.
        for kern in ("gaussian", "power"):
            for L in LENGTHS:
                f, a = run(data, cell, length=L, kernel=kern)
                record("A", f"river, {kern}, {L:g} km", f, a)
        print()
    if ONLY in (None, "B"):
        print("B. unequal accumulation spans")
        for spread in WINDOW_SPREADS:
            rng = np.random.default_rng(4242)
            win = (np.full(n_sites, 8) if spread == 0
                   else np.clip(rng.integers(8 - spread, 8 + spread + 1, n_sites), 1, 20))
            f, a = run(data, cell, window=win)
            record("B", f"window 8 +/- {spread}", f, a)
        print()
    if ONLY in (None, "C"):
        print("C. innovation that respects boundaries")
        for s in BOUNDARY_STRENGTHS:
            f, a = run(data, cell, strength=s)
            record("C", f"boundary strength {s}", f, a)
        print()



def main() -> int:
    summary = json.loads((ROOT / "docs" / "manuscript" / "revisions" / "2026-09-09"
                          / "results" / "summary.json").read_text())
    cell_by_region = {c["region"]: c for c in summary["calibration"]
                      if c["model"] == "pooled"}
    if REPORT_ONLY:
        # Rewrite the report from the saved draws. The sweep costs hours; the
        # wording of its verdict should not. Same flag as 64.
        rows = json.loads(OUT_JSON.read_text())
        rows = [r for r in rows if r.get("region") in REGIONS]
    else:
        sets = r47.load_sets()
        rows = []
        for region in REGIONS:
            run_region(sets, region, cell_by_region, rows)
        OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
        OUT_JSON.write_text(json.dumps(rows, indent=2))

    L = ["# The three remaining departures from the drift null", "",
         f"Produced by `analyses/65_other_departures.py`{' (--fast)' if FAST else ''}, "
         f"{N_REAL} realizations per cell, each region at its own calibrated cell.", "",
         "The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a "
         "comparison, not a replication, because those deposits are on the whole earlier "
         "and were collected under a different sampling regime.", ""]
    for reg in sorted({r.get("region") for r in rows}):
        c = cell_by_region[reg]
        o = next(r["observed"] for r in rows if r.get("region") == reg)
        L.append(f"- **{reg}**: observed $F_{{ST}}$ {o:.4f}; calibrated cell N {c['n_ind']}, "
                 f"innovation {c['innovation']}, mixing {c['mixing']}.")
    L += ["",
          "A cell closes the gap only if it reaches its region's observed value AND stays "
          "inside the calibration's 10 percent diversity tolerance.", "",
          "| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |",
          "|---|---|---|---:|---|---:|---:|:---:|"]
    names = {"A": "transport geography", "B": "accumulation spans",
             "C": "innovation boundary"}
    regions = []
    for r in rows:
        if r.get("region") not in regions:
            regions.append(r.get("region"))
    for r in rows:
        L.append(f"| {r.get('region','')} | {names[r['departure']]} | {r['setting']} | {r['fst_median']:.4f} | "
                 f"[{r['fst_lo']:.4f}, {r['fst_hi']:.4f}] | "
                 f"{r['observed'] / max(r['fst_median'], 1e-9):.1f}x | "
                 f"{100 * r['reach']:.1f}% | "
                 f"{'matched' if r['diversity_matched'] else 'broken'} |")
    L += ["", "## Reading", ""]
    for region in regions:
        L.append(f"\n**{region}**\n")
        for dep in ("A", "B", "C"):
            rs = [r for r in rows
                  if r["departure"] == dep and r.get("region") == region]
            if not rs:
                continue
            ok = [r for r in rs if r["diversity_matched"]]
            short = lambda r: r["observed"] / max(r["fst_median"], 1e-9)
            best = min(rs, key=short)
            best_ok = min(ok, key=short) if ok else None
            L.append(
                f"- **{names[dep]}**: best shortfall {short(best):.1f}x at "
                f"{best['setting']}"
                + (f"; the best cell that keeps diversity matched is "
                   f"{best_ok['setting']} at {short(best_ok):.1f}x."
                   if best_ok else ", and no cell keeps diversity matched.")
            )
    L.append("")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
