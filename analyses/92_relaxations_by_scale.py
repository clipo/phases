"""92_relaxations_by_scale.py - the relaxations of the copying model, scored at every scale that fails.

Analyses 64 and 65 relax four simplifications of the copying model and score
each variant at TWO clusters, where the calibrated baseline already covers the
observation within its predictive range (12 percent of runs reach it). The
failure that motivates the relaxations is at three and four clusters and at the
phase lines, where 0 to 3 percent of baseline runs reach the observation. The
review of 2026-09-24 (item 3) points out that a relaxation judged only at two
clusters, and only by whether its MEDIAN reaches the observation, is judged by
a different standard from the baseline, which is credited because its RANGE
covers the observation.

This script reruns the main settings of each relaxation with 64's and 65's
seeds and scores every run at two, three and four k-means clusters (analysis
71's labels) and at the phase lines, reporting for each the median, the 95
percent range, and the share of runs reaching the observed value, beside the
diversity match. The same three readouts are given for the baseline, so every
account is read by one standard: median discrepancy, predictive coverage, and
joint fit to diversity.

THE DOSE OF LOCAL INNOVATION. `65_other_departures.group_targets` reorders an
area's profile only when two or more classes are drawn (each with probability
equal to the strength) and the reordering is not the identity, so the same
strength changes far fewer profiles than its name suggests, and a row with
many small areas is not the same dose as a row with two large ones. The last
section of the output MEASURES the dose of every local-innovation row: the
share of areas whose profile differs from the regional one at all, and the
share of the profile's frequency that sits in a different class. Read a
local-innovation row beside its dose before reading it against another row.

All settings run at the calibrated combination; none is recalibrated. A setting
that breaks the diversity match here might fit after recalibration, so a
negative result is a statement about one-at-a-time perturbations of the
calibrated model, not about the relaxed model in general.

Output: output/findings/relaxations_by_scale.md
Usage: .venv/bin/python analyses/92_relaxations_by_scale.py [--reps 300]
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
OUT_MD = ROOT / "output" / "findings" / "relaxations_by_scale.md"
TOL = 0.10          # the calibration's diversity tolerance
K_SCORED = (2, 3, 4)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    args = ap.parse_args()
    sys_argv = sys.argv; sys.argv = [sys.argv[0]]       # 65 parses sys.argv at import
    r47 = importlib.import_module("47_revision_analysis")
    sd = importlib.import_module("23_phases_vs_spatial_drift")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    a64 = importlib.import_module("64_unequal_populations")
    a65 = importlib.import_module("65_other_departures")
    sys.argv = sys_argv
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record

    data = r47.load_sets()["basin"]
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    if names != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    centred = coords.to_numpy(float) - coords.to_numpy(float).mean(0)
    parts = {f"{k} clusters": mf._kmeans_labels(centred, k, seed=7) for k in K_SCORED}
    labels_ph, _ = ph.assign_primary_phases(names, coords.to_numpy(float))
    plist = sorted(set(labels_ph))
    parts["phases"] = np.array([plist.index(l) for l in labels_ph])
    a, b = np.asarray(parts["2 clusters"]), np.asarray(data["labels"])
    if len({(x, y) for x, y in zip(a, b)}) != len(set(a)) or len(set(a)) != len(set(b)):   # same division up to relabeling
        raise RuntimeError("two-cluster labels differ from the loader's spatial clusters")
    m_obs = data["m"]
    obs = {p: float(r47.fst_by(m_obs, lab)) for p, lab in parts.items()}
    obs_div = r47.diversity(m_obs)
    totals = m_obs.sum(1)
    rates = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates[(rates.region == "basin") & (rates.model == "pooled")].iloc[0]
    cell = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))
    straight = sd.geo_km(data["coords"])
    part_copy = data["parkin"] if data["parkin"] is not None else data["labels"]

    def run_65(dist=None, length=24.0, strength=0.0, kernel="exponential", areas=None):
        """Analysis 65's run(), with its seeds, scored at every partition.

        `areas` labels the areas within which new variants share one reordered
        profile; by default the two spatial clusters, as in analysis 65."""
        area_labels = data["labels"] if areas is None else np.asarray(areas)
        d = data["d"] if dist is None else dist
        out = {p: [] for p in parts}; divs = []
        pooled = np.asarray(data["pooled"], float)
        changed, moved = [], []
        for s in range(args.reps):
            rng = np.random.default_rng(70000 + s)
            w = copying_weights(d, length, labels=part_copy, leak=1.0, kernel=kernel)
            tgt = a65.group_targets(data["pooled"], area_labels, strength, rng)
            if strength > 0:
                # one row per area: every node of an area carries the same profile
                rows = np.array([tgt[np.flatnonzero(area_labels == g)[0]] for g in np.unique(area_labels)])
                diff = np.abs(rows - pooled).sum(1) / (2.0 * pooled.sum())
                changed.append(float((diff > 0).mean()))
                moved.append(float(diff.mean()))
            rec = drift_record(w, k=m_obs.shape[1], n_ind=cell["n_ind"], seed=70000 + s,
                               innovation=cell["innovation"], mixing=cell["mixing"],
                               burnin=1200, initial="uniform", target=tgt)
            m = sample_record(rec, data["ranks"], totals, rng)
            for p, lab in parts.items():
                out[p].append(r47.fst_by(m, lab))
            divs.append(r47.diversity(m))
        if strength > 0:
            dose.append((len(np.unique(area_labels)), strength, float(np.mean(changed)), float(np.mean(moved))))
        return out, divs

    def run_64(cv):
        """Analysis 64's loop, with its seeds, scored at every partition."""
        out = {p: [] for p in parts}; divs = []
        for s in range(args.reps):
            rng = np.random.default_rng(90000 + s)
            n_ind = a64.populations(len(m_obs), cell["n_ind"], cv, rng)
            rec = r47.simulate(data, 90000 + s, "pooled", innovation=cell["innovation"],
                               mixing=cell["mixing"], n_ind=n_ind)
            m = r47.sample_record(rec, data["ranks"], totals, rng)
            for p, lab in parts.items():
                out[p].append(r47.fst_by(m, lab))
            divs.append(r47.diversity(m))
        return out, divs

    dose = []   # (areas, strength, share of areas changed, share of frequency moved), in run order
    hood5 = mf._kmeans_labels(centred, 5, seed=7)
    hood8 = mf._kmeans_labels(centred, 8, seed=7)
    by_site = np.arange(len(m_obs))
    settings = [
        ("baseline", "calibrated copying model", lambda: run_65()),
        ("unequal populations", "CV 0.5", lambda: run_64(0.5)),
        ("unequal populations", "CV 1.0", lambda: run_64(1.0)),
        ("unequal populations", "CV 1.5", lambda: run_64(1.5)),
        ("transport", "straight-line, exponential, 6 km", lambda: run_65(dist=straight, length=6.0)),
        ("transport", "river, exponential, 6 km", lambda: run_65(length=6.0)),
        ("transport", "river, exponential, 12 km", lambda: run_65(length=12.0)),
        ("transport", "river, gaussian, 24 km", lambda: run_65(length=24.0, kernel="gaussian")),
        ("transport", "river, gaussian, 48 km", lambda: run_65(length=48.0, kernel="gaussian")),
        ("local innovation", "strength 0.2", lambda: run_65(strength=0.2)),
        ("local innovation", "strength 0.25", lambda: run_65(strength=0.25)),
        ("local innovation", "strength 0.3", lambda: run_65(strength=0.3)),
        # The same mechanism at finer grain (2026-10-02): new variants shared
        # within a neighborhood of sites, or arising site by site.
        ("local innovation, five neighborhoods", "strength 0.2", lambda: run_65(strength=0.2, areas=hood5)),
        ("local innovation, five neighborhoods", "strength 0.25", lambda: run_65(strength=0.25, areas=hood5)),
        ("local innovation, eight neighborhoods", "strength 0.2", lambda: run_65(strength=0.2, areas=hood8)),
        ("local innovation, eight neighborhoods", "strength 0.25", lambda: run_65(strength=0.25, areas=hood8)),
        ("local innovation, each site its own", "strength 0.1", lambda: run_65(strength=0.1, areas=by_site)),
        ("local innovation, each site its own", "strength 0.2", lambda: run_65(strength=0.2, areas=by_site)),
    ]
    L = ["# The relaxations scored at every scale", "",
         f"Produced by `analyses/92_relaxations_by_scale.py`: {args.reps} runs per setting at the calibrated "
         f"combination (N {cell['n_ind']}, innovation {cell['innovation']}, mixing {cell['mixing']}), with "
         "analysis 64's and 65's seeds; local innovation areas are the two spatial clusters, as in 65, unless the "
         "row names neighborhoods (the five- and eight-cluster k-means divisions of the site map) or single sites. Each "
         "cell gives the median between-group F_ST, its 95 percent range, and the share of runs at or above "
         "the observed value. Diversity is matched when the medians of all three within-assemblage summaries "
         "lie within 10 percent of the observed values. No setting is recalibrated.", "",
         "| account | setting | " + " | ".join(f"{p} (obs {obs[p]:.4f})" for p in parts) + " | diversity |",
         "|---|---|" + "---|" * len(parts) + "---|"]
    dose_rows = []
    for acc, label, fn in settings:
        n_dose = len(dose)
        out, divs = fn()
        if len(dose) > n_dose:
            areas_n, _, ch, mv = dose[-1]
            dose_rows.append(f"| {acc} | {label} | {areas_n} | {ch * 100:.0f}% | {mv * 100:.1f}% |")
        ach = {k: float(np.median([x[k] for x in divs])) for k in ("hs", "rich", "ht")}
        ok = all(abs(ach[k] - obs_div[k]) <= TOL * abs(obs_div[k]) for k in ("hs", "rich", "ht"))
        cells = []
        for p in parts:
            v = np.asarray(out[p], float); v = v[np.isfinite(v)]
            lo, med, hi = np.percentile(v, [2.5, 50, 97.5])
            cells.append(f"{med:.4f} [{lo:.4f}, {hi:.4f}]; {np.mean(v >= obs[p]) * 100:.0f}%")
        L.append(f"| {acc} | {label} | " + " | ".join(cells) + f" | {'matched' if ok else 'broken'} |")
        print(L[-1], flush=True)
    L += ["", "## Reading", "",
          "Read each row on three things at once: whether the median reaches the observation, whether the "
          "95 percent range covers it, and whether diversity stays matched. The baseline row is the reference; "
          "whether its range covers a finer-scale value depends on the seeds (compare Table 1), so read the "
          "share of runs reaching it rather than the range edge. A relaxation improves on the baseline where it "
          "raises that share with diversity matched. Settings that break the diversity match here were not "
          "recalibrated.", "",
          "## The dose behind each local-innovation row", "",
          "An area's profile is reordered only when two or more of its classes are drawn, each with "
          "probability equal to the strength, and the reordering is not the identity. Measured over the "
          f"same {args.reps} runs: the mean share of areas whose profile differs from the regional one at all, "
          "and the mean share of a profile's frequency that sits in a different class (half the summed "
          "absolute difference from the regional profile, averaged over areas, unchanged ones included). "
          "Two rows are the same dose only when these agree; a row that reaches the finer-scale values no "
          "more often than the baseline at a small dose has not tested the grain it names.", "",
          "| account | setting | areas | areas whose profile differs | frequency moved, mean |",
          "|---|---|---|---|---|", *dose_rows, ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
