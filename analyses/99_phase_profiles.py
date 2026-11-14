"""99_phase_profiles.py - how different are the phases, in sherds?

The paper reports the difference between phases as a cultural F_ST of 0.0123
against a model median of 0.0024. Those are correct and hard to picture. This
script states the same comparison in the units a reader of a sherd table uses:
what share of each phase's decorated sherds belongs to each class, and what
share of one phase's sherds would have to be moved to another class to turn its
profile into a neighbor's.

THE MEASURES.
  profile        each phase's decorated sherds pooled, as percentages by class.
  difference     between two phases, half the summed absolute difference of
                 their profiles, in percent: the share of sherds that would
                 have to change class to make the two profiles identical. It
                 is 0 for identical profiles and 100 for profiles with no class
                 in common. The mean over the three pairs of phases is reported.
  leading class  the range, across phases, of the share of the basin's
                 commonest class, in percentage points.

THE MODEL. The calibrated neutral copying model with no groups, exactly as
analysis 84 runs it for its no-boundary, regional-pool setting (same seeds, so
these are the same simulated records), sampled at the real sherd counts and
pooled by the same phases.

Descriptive (rule 18): observed values beside the model's median and 95 percent
range, and the share of runs that reach the observed value.

Output: output/findings/phase_profiles.md
Usage: PYTHONPATH=src .venv/bin/python analyses/99_phase_profiles.py [--reps 300]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "phase_profiles.md"
PARKIN_SEED_86 = 86200   # the seed analysis 86 gives its Parkin-against-the-rest comparison


def phase_profiles(m: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """Percent of each phase's pooled sherds in each class (phases x classes)."""
    m = np.asarray(m, float)
    pooled = np.array([m[phase == g].sum(0) for g in np.unique(phase)])
    return 100.0 * pooled / pooled.sum(1, keepdims=True)


def mean_pair_difference(prof: np.ndarray) -> float:
    """Mean over pairs of phases of half the summed absolute profile difference (percent)."""
    return float(np.mean([0.5 * np.abs(prof[a] - prof[b]).sum()
                          for a, b in combinations(range(len(prof)), 2)]))


def pair_differences(prof: np.ndarray) -> dict:
    """Half the summed absolute profile difference (percent) for every pair of phases."""
    return {(a, b): float(0.5 * np.abs(prof[a] - prof[b]).sum())
            for a, b in combinations(range(len(prof)), 2)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()
    if not 20 <= args.reps <= 1000:
        raise ValueError("--reps must be between 20 and 1000 (analysis 84's seed family)")

    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record

    data = rev.load_sets()["basin"]
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    if names != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    labels, _ = ph.assign_primary_phases(names, coords.to_numpy(float))
    phases = sorted(set(labels))
    phase = np.array([phases.index(l) for l in labels])
    classes = [str(c).replace("_", " ") for c in counts.columns]
    m_obs = data["m"]
    totals = m_obs.sum(1)

    obs = phase_profiles(m_obs, phase)
    lead = int(np.argmax(m_obs.sum(0)))
    obs_diff = mean_pair_difference(obs)
    obs_lead = float(obs[:, lead].max() - obs[:, lead].min())
    obs_fst = float(rev.fst_by(m_obs, phase))

    rates = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates[(rates.region == "basin") & (rates.model == "pooled")].iloc[0]
    w = copying_weights(data["d"], 24.0, labels=phase, leak=1.0)
    diffs, leads, fsts, bypair = [], [], [], []
    for r in range(args.reps):
        seed = 84000 + r                                  # analysis 84, no boundary, regional pool
        rng = np.random.default_rng(seed + 5_000_000)
        rec = drift_record(w, k=m_obs.shape[1], n_ind=int(row["n_ind"]), seed=seed,
                           innovation=float(row["innovation"]), mixing=float(row["mixing"]),
                           burnin=1200, initial="uniform", target=data["pooled"])
        m = sample_record(rec, data["ranks"], totals, rng)
        prof = phase_profiles(m, phase)
        ml = int(np.argmax(m.sum(0)))                     # the run's own commonest class
        diffs.append(mean_pair_difference(prof))
        bypair.append(pair_differences(prof))
        leads.append(float(prof[:, ml].max() - prof[:, ml].min()))
        fsts.append(float(rev.fst_by(m, phase)))
    diffs, leads, fsts = np.array(diffs), np.array(leads), np.array(fsts)

    q = lambda a: np.percentile(a, [2.5, 50, 97.5])
    L = ["# How different are the phases, in sherds?", "",
         f"Produced by `analyses/99_phase_profiles.py`. {len(names)} assemblages, phases "
         + ", ".join(f"{p} {int((phase == i).sum())} ({int(m_obs[phase == i].sum()):,} decorated sherds)"
                     for i, p in enumerate(phases))
         + f". Model: the calibrated copying model with no groups, {args.reps} runs with analysis 84's seeds "
         "(its no-boundary, regional-pool setting), sampled at the real sherd counts.", "",
         "## Each phase's decorated sherds, percent by class", "",
         "| class | " + " | ".join(phases) + " | all |", "|---|" + "---:|" * (len(phases) + 1)]
    allp = 100.0 * m_obs.sum(0) / m_obs.sum()
    for j in np.argsort(-allp):
        L.append(f"| {classes[j]} | " + " | ".join(f"{obs[i, j]:.1f}" for i in range(len(phases)))
                 + f" | {allp[j]:.1f} |")
    L += ["", "## The record against the model", "",
          "| measure | observed | model median (95 percent range) | runs reaching observed |",
          "|---|---:|---|---:|"]
    for label, o, a, fmt in [
            ("sherds that would have to change class to make two phases alike, mean over pairs (percent)",
             obs_diff, diffs, "{:.1f}"),
            (f"range across phases in the share of the commonest class, {classes[lead]} (percentage points)",
             obs_lead, leads, "{:.1f}"),
            ("between-phase cultural F_ST", obs_fst, fsts, "{:.4f}")]:
        lo, med, hi = q(a)
        L.append(f"| {label} | {fmt.format(o)} | {fmt.format(med)} ({fmt.format(lo)} to {fmt.format(hi)}) | "
                 f"{int((a >= o).sum())} of {len(a)} ({(a >= o).mean():.0%}) |")
    obs_pairs = pair_differences(obs)
    L += ["", "## Pair by pair", "",
          "The share of sherds that would have to change class to make the two phases' profiles alike.", "",
          "| pair | observed (percent) | model median (95 percent range) | runs reaching observed |",
          "|---|---:|---|---:|"]
    for (a, b), o in obs_pairs.items():
        arr = np.array([d[(a, b)] for d in bypair])
        lo, med, hi = q(arr)
        L.append(f"| {phases[a]} and {phases[b]} | {o:.1f} | {med:.1f} ({lo:.1f} to {hi:.1f}) | {int((arr >= o).sum())} of {len(arr)} ({(arr >= o).mean():.0%}) |")
    # ---- Is each phase's difference specific to its own line?
    # The model comparison above says the phases differ more than copying with
    # no groups makes them. It does not say the difference follows the phase
    # lines. Analysis 86's comparison does: each phase against the rest, beside
    # divisions of the same sizes around random centers and made compact.
    t86 = importlib.import_module("86_partition_posterior")
    t74 = importlib.import_module("74_phase_partition_test")
    pts = t74.km_xy(coords.to_numpy(float))
    gaps = {(x, y): float(min(np.linalg.norm(pts[i] - pts[j]) for i in np.flatnonzero(phase == x)
                              for j in np.flatnonzero(phase == y)))
            for x, y in combinations(range(len(phases)), 2)}
    L += ["", "Straight-line distance between the nearest assemblages of two phases: "
          + "; ".join(f"{phases[x]} and {phases[y]}, {g:.0f} km" for (x, y), g in gaps.items()) + "."]
    L += ["", "## Each phase against the rest, beside same-size divisions", "",
          ("Analysis 86's comparison (2,000 posterior draws of the counts, 300 alternative divisions of each "
          "kind), with one phase against the other two taken together, at the seed analysis 86 uses for the "
          "Parkin phase, so the Parkin row is that analysis's own result."), "",
          ("| division | assemblages | F_ST, posterior median (95% CrI) | median F_ST, random-center / compact "
          "divisions | P(separates better than random-center / compact) | boundary excess | "
          "P(larger than random-center / compact) |"),
          "|---|---|---|---|---|---|---|"]
    for i, pname in enumerate(phases):
        idx = (phase == i).astype(int)
        r = t86.compare(m_obs.astype(float), idx, pts, rev, t74, 2000, 300, PARKIN_SEED_86, dist=data["d"])
        if pname == "Parkin":
            # One definition (rule 6): the paper already reports Parkin against the
            # rest from analysis 86. This row must be that number, not a near copy.
            ref = (ROOT / "output" / "findings" / "partition_posterior.md").read_text()
            sec = ref[ref.index("## The Parkin phase against the rest of the basin"):]
            want = [f"= **{r['p_alt']:.2f}**", f"= **{r['p_opt']:.2f}**"]
            if not all(w in sec for w in want):
                raise RuntimeError("the Parkin row does not reproduce analysis 86's reported probabilities "
                                   f"({r['p_alt']:.2f}, {r['p_opt']:.2f}); check the seed and the labels")
        L.append(f"| {pname} against the rest | {int(idx.sum())} of {len(idx)} | {r['phase'][0]:.4f} "
                 f"({r['phase'][1]:.4f} to {r['phase'][2]:.4f}) | {r['alt'][0]:.4f} / {r['opt'][0]:.4f} | "
                 f"{r['p_alt']:.2f} / {r['p_opt']:.2f} | {r['be_phase'][0]:+.1f} | "
                 f"{r['be_p_alt']:.2f} / {r['be_p_opt']:.2f} |")
    L += ["", "## Reading", "",
          ("The three measures say the same thing in different units. In a run, the commonest class is "
          "whichever class is commonest in that run, so the second row compares like with like. The "
          "F_ST row repeats Table 1's comparison and is here so the three can be read side by side; it "
          "should equal analysis 84's no-boundary row."), "",
          ("The last table asks a different question from the first two. A phase can differ from its "
          "neighbors by more than the copying model produces and still be no more distinct than any "
          "other group of neighboring assemblages of its size; a probability near one half in that "
          "table says exactly that."), ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
