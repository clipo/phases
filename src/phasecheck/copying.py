"""What copying across distance produces on this map, with no groups, and what the record can detect.

Four steps, in the order a reader needs them:

  calibrate            find settings of the copying model that reproduce the
                       record's own diversity, and among them the one that
                       gives the phases the largest difference (the no-groups
                       model at its most favorable)
  compare_with_model   does the record differ between phases by more than the
                       model produces?
  sorting_under_model  would records made with no groups sort into the phases?
  power                could the comparison with alternative lines detect a
                       restriction on copying at the phase boundaries here?

Every share reported is a count of simulated runs. None is a significance level.

These answer the report's questions 6 (calibrate), 7 (compare_with_model),
8 (sorting_under_model) and 9 (power). The model itself is in `model.py`.
Conventions: distances and the copying length in km; F_ST is the
Gini-Simpson cultural F_ST; the boundary excess is within-phase minus
between-phase similarity at matched distance (positive: the phase lines
separate assemblages more than distance alone predicts).

Seeds. Every simulated record has a drift seed and a sampling seed, and each
step uses its own family of seeds so no setting is judged on the runs it was
chosen on. Defaults here are the paper's (calibration 60,000 to about
100,000; `compare_with_model` and `power` 84,000; `sorting_under_model`
79,000); `report.run` passes its own families from 2,000,000 upward. Within a
family, run r uses base + r, so `MAX_RUNS` (5,000) keeps a family below the
10,000 spacing of `power`'s leak settings. When no sampling seed is given the
sampling seed is the drift seed + 5,000,000 (`simulate`).
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from . import divisions as dv
from . import measures as ms
from . import model as md
from . import questions as qs
from .data import Dataset

# The calibration grid, the paper's: 3 x 11 x 9 = 297 combinations.
N_IND = (120, 2000, 10000)                         # learners per assemblage
INNOVATIONS = (.0002, .0005, .001, .002, .004, .008, .012, .024, .048, .1, .2)
MIXINGS = (.001, .002, .005, .01, .02, .05, .1, .2, .4)
# Copying across a phase line, as a share of its strength. The position in this
# tuple sets each leak's seed family in `power`, so the order must not change.
LEAKS = (1.0, 0.5, 0.1, 0.03)
TOLERANCE = 0.10                       # a match: every diversity summary within 10 percent of the record
LENGTH_SHARE = 0.2                     # default copying length, as a share of the largest distance


def default_length(data: Dataset) -> float:
    """The copying length used when none is given: a fifth of the largest distance between assemblages.

    This is a choice, not an estimate. The paper used 24 km on a basin whose
    assemblages lie up to 268 km apart along the rivers (79 km in a straight
    line), so this default is not the paper's value for the paper's record.

    Parameters
    ----------
    data : Dataset
        Its `dist` (straight-line or supplied, km) sets the scale.

    Returns
    -------
    float
        `LENGTH_SHARE` times the largest pairwise distance, in km.
    """
    return float(LENGTH_SHARE * data.dist.max())


def diversity(m: np.ndarray) -> dict:
    """Mean diversity within assemblages, mean number of classes present, and diversity of the pooled record.

    These three summaries, none involving the phases, are what calibration
    matches (question 6).

    Parameters
    ----------
    m : array_like, shape (n, classes)
        Counts; every row must have a sherd.

    Returns
    -------
    dict
        within : mean Gini-Simpson diversity of the assemblages.
        classes : mean number of classes with at least one sherd.
        pooled : Gini-Simpson diversity of all sherds pooled.
    """
    m = np.asarray(m, float)
    p = m / m.sum(1, keepdims=True)
    pooled = m.sum(0) / m.sum()
    return {"within": float((1 - (p ** 2).sum(1)).mean()),
            "classes": float((m > 0).sum(1).mean()),
            "pooled": float(1 - (pooled ** 2).sum())}


MAX_RUNS = 5000          # runs per family; seed families are spaced 10,000 apart in `power`


def _fst(m: np.ndarray, labels: np.ndarray) -> float:
    """F_ST of a SIMULATED table; 0 when the run has lost all diversity (one class everywhere).

    The measure is undefined there, and for an observed record that is refused.
    A simulated run with one class has no difference between any groups, so it
    counts as zero and stays in the tally. Dropping such runs would keep only
    the ones that differ and overstate what the model produces.
    """
    if (np.asarray(m).sum(0) > 0).sum() < 2:
        return 0.0
    return ms.fst_by(m, labels)


def _order(data: Dataset) -> np.ndarray:
    """Sequence positions for `model.sample_record`; all ones (contemporaneous, latest) without an order."""
    return np.ones(len(data.names)) if data.order is None else data.order


def simulate(data: Dataset, rates: dict, seed: int, *, leak: float = 1.0, sample_seed=None) -> np.ndarray:
    """One simulated count table at the real sherd counts. `leak` below 1 restricts copying across phase lines.

    Parameters
    ----------
    data : Dataset
        Supplies distances, phases, sherd totals, sequence order and the
        pooled class proportions used as the innovation profile.
    rates : dict
        n_ind (learners), innovation, mixing, length (km).
    seed : int
        Drift seed.
    leak : float
        0 to 1; copying between assemblages of different phases is
        multiplied by it. 1 leaves the phases out of the model entirely.
    sample_seed : int, optional
        Seed of the sherd sampling; default `seed + 5_000_000`, far enough
        above every drift seed in use that the two never coincide.

    Returns
    -------
    numpy.ndarray of int, shape (n, classes)
    """
    w = md.copying_weights(data.dist, rates["length"], data.phase_index if leak < 1 else None, leak)
    rec = md.drift_record(w, k=data.counts.shape[1], n_ind=rates["n_ind"], mixing=rates["mixing"],
                          innovation=rates["innovation"], target=data.counts.sum(0) / data.counts.sum(),
                          seed=seed)
    rng = np.random.default_rng(seed + 5_000_000 if sample_seed is None else sample_seed)
    return md.sample_record(rec, _order(data), data.counts.sum(1), rng)


# ------------------------------------------------------------------ calibration
def _cell(args):
    """Score one combination of settings: diversity of its runs and the phases' F_ST in them.

    Takes one tuple so it can go through `ProcessPoolExecutor.map`. Returns
    the mean of each diversity summary over the runs and the array of F_ST
    between the groups of `labels`, one per run.
    """
    data, n_ind, innovation, mixing, length, seeds, sample_seeds, labels = args
    rates = {"n_ind": n_ind, "innovation": innovation, "mixing": mixing, "length": length}
    div, fst = [], []
    for s, ss in zip(seeds, sample_seeds):
        m = simulate(data, rates, s, sample_seed=ss)
        div.append(diversity(m))
        fst.append(_fst(m, labels))
    mean = {k: float(np.mean([d[k] for d in div])) for k in div[0]}
    return mean, np.asarray(fst, float)


def calibrate(data: Dataset, *, length=None, n_ind=N_IND, innovations=INNOVATIONS, mixings=MIXINGS,
              reps: int = 6, confirm_reps: int = 50, tolerance: float = TOLERANCE, jobs: int = 1,
              select_labels=None, progress=None) -> dict:
    """Tune population size, innovation rate and mixing rate to the record's own diversity.

    A combination matches when the mean of its runs is within `tolerance` of
    the record on all three diversity summaries. None of the three involves the
    phases. Matches found on `reps` runs are confirmed on `confirm_reps` fresh
    runs. Among the confirmed, the combination whose runs give the largest
    median F_ST between the groups of `select_labels` (the phases by default)
    is selected: the model with no groups at its most favorable, so that a
    shortfall is not the result of a poorly chosen setting.

    Raises if nothing matches. A model that cannot reproduce the record's
    diversity must not be used to judge its differences.

    Parameters
    ----------
    data : Dataset
    length : float, optional
        Copying length in km; default `default_length(data)`.
    n_ind, innovations, mixings : sequences
        The grid; default the paper's 297 combinations.
    reps : int
        Screening runs per combination, at least 2.
    confirm_reps : int
        Confirming runs per matched combination, 10 to 100 (each
        combination has 100 seeds of its own).
    tolerance : float
        Largest relative difference allowed on each diversity summary.
    jobs : int
        Worker processes; 1 runs in this process.
    select_labels : array_like, optional
        Groups whose F_ST picks among the confirmed combinations; default
        the phases. The paper selected by two spatial clusters.
    progress : callable, optional
        Called with a one-line message at each stage.

    Returns
    -------
    dict
        rates : the selected n_ind, innovation, mixing and length.
        table : DataFrame, one row per combination: settings, `matched`,
            `loss` (sum of squared relative differences), `sim_within`,
            `sim_classes`, `sim_pooled` (screening means), `confirmed`,
            `fst_median` (over the confirming runs), `selected`.
        observed : the record's `diversity`.
        n_cells, n_matched, n_confirmed : counts of combinations.
        reps, confirm_reps, tolerance : the settings used.
        ranked : the confirmed combinations by `fst_median`, largest first
            (ties kept in grid order).

    Raises
    ------
    ValueError
        If `reps` or `confirm_reps` is out of range, or no combination is
        confirmed (the message gives the nearest one).

    Seeds: screening runs s = 0..reps-1 use drift seed 60,000 + s and
    sampling seed 61,000 + s, the same for every combination, so screening
    compares settings on common random numbers. Confirming runs of the
    combination in row j use 70,000 + 100 j + s and 71,000 + 100 j + s,
    fresh seeds that the screening never saw.
    """
    if reps < 2 or not 10 <= confirm_reps <= 100:
        raise ValueError("calibration needs at least 2 screening runs and 10 to 100 confirming runs per "
                         "combination (each combination has 100 seeds of its own)")
    length = default_length(data) if length is None else float(length)
    labels = data.phase_index if select_labels is None else np.asarray(select_labels)
    obs = diversity(data.counts)
    cells = [(n, i, x) for n in n_ind for i in innovations for x in mixings]
    say = progress or (lambda *_: None)

    def run(jobs_args):
        """Score every combination, in worker processes when jobs > 1; results keep the input order."""
        if jobs > 1:
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                return list(ex.map(_cell, jobs_args, chunksize=4))
        return [_cell(a) for a in jobs_args]

    say(f"calibration: screening {len(cells)} combinations, {reps} runs each")
    screen = run([(data, n, i, x, length, [60000 + s for s in range(reps)],
                   [61000 + s for s in range(reps)], labels) for n, i, x in cells])
    rows = []
    for (n, i, x), (mean, fst) in zip(cells, screen):
        rel = {k: abs(mean[k] - obs[k]) / obs[k] for k in obs}      # relative miss on each summary
        rows.append(dict(n_ind=n, innovation=i, mixing=x, matched=max(rel.values()) <= tolerance,
                         loss=sum(v ** 2 for v in rel.values()), **{f"sim_{k}": v for k, v in mean.items()}))
    frame = pd.DataFrame(rows)
    frame["confirmed"], frame["fst_median"] = False, np.nan
    idx = list(frame.index[frame.matched])
    say(f"calibration: {len(idx)} matched; confirming on {confirm_reps} runs each")
    confirm = run([(data, int(frame.n_ind[j]), float(frame.innovation[j]), float(frame.mixing[j]), length,
                    [70000 + 100 * j + s for s in range(confirm_reps)],
                    [71000 + 100 * j + s for s in range(confirm_reps)], labels) for j in idx])
    for j, (mean, fst) in zip(idx, confirm):
        rel = {k: abs(mean[k] - obs[k]) / obs[k] for k in obs}
        frame.loc[j, "confirmed"] = bool(max(rel.values()) <= tolerance)
        # F_ST is defined on every run (_fst returns 0 for a run with one class), so the
        # NaN guard is a safeguard only.
        frame.loc[j, "fst_median"] = float(np.nanmedian(fst)) if np.isfinite(fst).any() else np.nan
    ok = frame[frame.confirmed & frame.fst_median.notna()]
    if not len(ok):
        best = frame.loc[frame.loss.idxmin()]
        raise ValueError(
            "no setting of the copying model reproduces this record's diversity within "
            f"{100 * tolerance:.0f} percent ({int(frame.matched.sum())} matched on screening, none confirmed). "
            f"The nearest gives within-assemblage diversity {best.sim_within:.3f}, classes present "
            f"{best.sim_classes:.1f} and pooled diversity {best.sim_pooled:.3f}, against {obs['within']:.3f}, "
            f"{obs['classes']:.1f} and {obs['pooled']:.3f} in the record. The model comparison is not "
            "reported, because a model that misses the diversity cannot be used to judge the differences.")
    # The largest median F_ST among the confirmed: the no-groups model at its most favorable.
    best = ok.loc[ok.fst_median.idxmax()]
    frame["selected"] = frame.index == best.name
    rates = {"n_ind": int(best.n_ind), "innovation": float(best.innovation),
             "mixing": float(best.mixing), "length": length}
    return {"rates": rates, "table": frame, "observed": obs, "n_cells": len(cells),
            "n_matched": int(frame.matched.sum()), "n_confirmed": len(ok), "reps": reps,
            "confirm_reps": confirm_reps, "tolerance": tolerance,
            "ranked": ok.sort_values("fst_median", ascending=False, kind="stable")[
                ["n_ind", "innovation", "mixing", "fst_median"]].reset_index(drop=True)}


# ------------------------------------------------------------------ the record against the model
def _row(label, observed, sims, fmt):
    """One line of the question 7 table: observed value, the runs' median and 95 percent range, runs reaching it."""
    a = np.asarray(sims, float)
    a = a[np.isfinite(a)]
    return {"measure": label, "observed": float(observed), "median": float(np.median(a)),
            "lo": float(np.percentile(a, 2.5)), "hi": float(np.percentile(a, 97.5)),
            "reaching": int((a >= observed).sum()), "runs": len(a), "fmt": fmt}


def compare_with_model(data: Dataset, rates: dict, *, reps: int = 300, seed: int = 84000) -> dict:
    """The record beside `reps` runs of the copying model with no groups.

    For each measure: the observed value, the model's median and 95 percent
    range, and the number of runs that reach the observed value. A measure the
    runs often reach is one the model accounts for. A measure they rarely
    reach is a difference copying across distance, as modeled, does not
    produce; it is not by itself evidence of a boundary.

    Parameters
    ----------
    data : Dataset
    rates : dict
        n_ind, innovation, mixing, length (km).
    reps : int
        Runs, 50 to `MAX_RUNS`.
    seed : int
        Run r uses drift seed `seed + r` and sampling seed
        `seed + r + 5_000_000`. The default, 84,000, is the paper's; with
        the paper's settings it reproduces 7 of 300 runs reaching the
        observed between-phase F_ST.

    Returns
    -------
    dict
        rows : list of dicts with measure, observed, median, lo, hi (2.5th
            and 97.5th percentiles of the runs), reaching (runs at or above
            the observed value), runs, fmt (format for the report). In
            order: F_ST between the phases; boundary excess at the phase
            lines; mean profile difference over pairs of phases; then, with
            three or more phases, each pair's profile difference and each
            phase's F_ST against the rest.
        runs : `reps`.
        flat : runs that held a single class (counted as F_ST 0, left out
            of the decay correlation).
        rates : the settings used.
        decay : correlation of pairwise similarity with distance, observed
            (None when undefined), and the runs' median, lo and hi.
            Negative means similarity falls with distance.

    Raises
    ------
    ValueError
        If `reps` is out of range, or fewer than 10 runs keep enough
        diversity for the decay correlation.
    """
    if not 50 <= reps <= MAX_RUNS:
        raise ValueError(f"between 50 and {MAX_RUNS} runs are needed to describe the model's range")
    idx, names, counts = data.phase_index, data.phase_names, data.counts.astype(float)
    iu = np.triu_indices(len(counts), 1)
    def decay(m):
        """Pearson correlation of pairwise Brainerd-Robinson similarity with distance (NaN if all pairs equal)."""
        s = ms.similarity_matrix(m)[iu]
        return float("nan") if np.ptp(s) == 0 else float(np.corrcoef(s, data.dist[iu])[0, 1])

    prof = lambda m: ms.pooled_profiles(m, idx)
    pairs = [(a, b) for a in range(len(names)) for b in range(a + 1, len(names))]
    sims = {k: [] for k in ["fst", "excess", "decay", "diff"]}
    rest = {p: [] for p in names}
    pair = {ab: [] for ab in pairs}
    flat = 0
    for r in range(reps):
        m = simulate(data, rates, seed + r).astype(float)
        sims["fst"].append(_fst(m, idx))
        sims["excess"].append(ms.boundary_excess(m, data.dist, idx))
        if (m.sum(0) > 0).sum() < 2:
            flat += 1                            # one class everywhere: no correlation to take
        elif np.isfinite(dc := decay(m)):
            sims["decay"].append(dc)
        pr = prof(m)
        d = {ab: ms.profile_difference(pr[ab[0]], pr[ab[1]]) for ab in pairs}
        sims["diff"].append(float(np.mean(list(d.values()))))
        for ab in pairs:
            pair[ab].append(d[ab])
        if len(names) > 2:
            for i, p in enumerate(names):
                rest[p].append(_fst(m, (idx == i).astype(int)))
    po = prof(counts)
    rows = [_row("cultural F_ST between the phases", ms.fst_by(counts, idx), sims["fst"], "{:.4f}"),
            _row("boundary excess at the phase lines", ms.boundary_excess(counts, data.dist, idx),
                 sims["excess"], "{:+.1f}"),
            _row("sherds that would have to change class to make two phases alike, mean over pairs (percent)",
                 np.mean([ms.profile_difference(po[a], po[b]) for a, b in pairs]), sims["diff"], "{:.1f}")]
    if len(pairs) > 1:
        for a, b in pairs:
            rows.append(_row(f"sherds that would have to change class, {names[a]} and {names[b]} (percent)",
                             ms.profile_difference(po[a], po[b]), pair[(a, b)], "{:.1f}"))
    if len(names) > 2:
        for i, p in enumerate(names):
            rows.append(_row(f"cultural F_ST, {p} against the rest", ms.fst_by(counts, (idx == i).astype(int)),
                             rest[p], "{:.4f}"))
    obs_decay = decay(counts)
    dec = np.asarray(sims["decay"])
    if len(dec) < 10:
        raise ValueError(f"{flat} of {reps} runs of the model lost all diversity at these settings; "
                         "they cannot stand in for this record")
    return {"rows": rows, "runs": reps, "flat": flat, "rates": rates,
            "decay": {"observed": obs_decay if np.isfinite(obs_decay) else None,
                      "median": float(np.median(dec)), "lo": float(np.percentile(dec, 2.5)),
                      "hi": float(np.percentile(dec, 97.5))}}


def sorting_under_model(data: Dataset, rates: dict, *, reps: int = 300, seed: int = 79000,
                        sample_offset: int = 500, weights=(0.0, 0.25, 0.5, 1.0, np.inf)) -> dict:
    """Would assemblages made by copying with no groups sort into the phases?

    Each simulated record is clustered by Ward linkage exactly as the real one
    is, at each weight on composition, and scored against the phases. The map
    alone is the same for every record, so its row cannot differ.

    Parameters
    ----------
    data : Dataset
    rates : dict
        n_ind, innovation, mixing, length (km).
    reps : int
        Simulated records, 50 to `MAX_RUNS` and no more than
        `sample_offset`.
    seed : int
        Record r uses drift seed `seed + r` and sampling seed
        `seed + sample_offset + r`.
    sample_offset : int
        Gap between the two seed ranges; must be at least `reps` so they
        do not overlap.
    weights : sequence of float
        Weights on composition, as in question 1 (0 the map, inf
        composition alone).

    Returns
    -------
    dict
        runs : number of records.
        rows : one dict per weight with weight, observed (the real record's
            ARI with the phases), median, lo, hi over the records, and
            reaching (records at or above the observed ARI, to 1e-12).

    Raises
    ------
    ValueError
        If `reps` is out of range or exceeds `sample_offset`.
    """
    idx, k = data.phase_index, len(data.phase_names)
    g = qs._unit_scale(data.pts)

    def recover(m):
        """Ward-linkage ARI with the phases at each weight, for one count table."""
        c = qs._unit_scale(qs._chisq_features(np.asarray(m, float)))
        out = []
        for w in weights:
            f = g if w == 0 else (c if np.isinf(w) else np.column_stack([g, np.sqrt(w) * c]))
            out.append(ms.adjusted_rand(idx, qs._cluster(f, k, "ward", 0, 0)))
        return out

    if not 50 <= reps <= MAX_RUNS:
        raise ValueError(f"between 50 and {MAX_RUNS} runs are needed")
    if sample_offset < reps:
        raise ValueError(
            f"reps ({reps}) must not exceed sample_offset ({sample_offset}), "
            "or run seeds and sample seeds overlap"
        )
    obs = recover(data.counts)
    sims = [recover(simulate(data, rates, seed + r, sample_seed=seed + sample_offset + r)) for r in range(reps)]
    arr = np.asarray(sims)
    return {"runs": len(sims),
            "rows": [{"weight": float(w), "observed": obs[j], "median": float(np.median(arr[:, j])),
                      "lo": float(np.percentile(arr[:, j], 2.5)),
                      "hi": float(np.percentile(arr[:, j], 97.5)),
                      # 1e-12: a record that ties the real one counts as reaching it
                      "reaching": int((arr[:, j] >= obs[j] - 1e-12).sum())}
                     for j, w in enumerate(weights)]}


def power(data: Dataset, rates: dict, *, records: int = 300, n_alt: int = 200, leaks=(1.0, 0.1, 0.03),
          seed: int = 84000) -> dict:
    """Could the comparison with alternative lines detect a restriction at the phase boundaries here?

    Records are simulated with copying across the phase lines unrestricted
    (leak 1) and restricted. Each is scored as the real record is: the share of
    alternative divisions the phase lines beat on cultural F_ST. If records
    with a restriction score much higher than records without, the comparison
    can detect one. The real record's own score is placed in each distribution.

    The defaults are larger than the paper's 100 records and 50 alternatives.
    At those sizes the estimated response to a restriction varied, between
    random seeds, from +0.04 to +0.31 at one setting of the paper's own record
    (measured 2026-10-06); at 300 and 200 it varied from +0.16 to +0.22.

    Parameters
    ----------
    data : Dataset
    rates : dict
        n_ind, innovation, mixing, length (km).
    records : int
        Simulated records per leak, 30 to `MAX_RUNS`.
    n_alt : int
        Alternative divisions of each kind; one set, shared by every record.
    leaks : sequence of float
        Each must be in `LEAKS`. 1 is no restriction.
    seed : int
        Base seed B. The alternatives use B + 999,999; record r at the leak
        in position i of `LEAKS` uses drift seed B + 10,000 i + r (sampling
        seed 5,000,000 above that). With the defaults (B = 84,000), the
        unrestricted records are the same runs as `compare_with_model`'s.

    Returns
    -------
    dict
        observed_around, observed_compact : share of alternatives of each
            kind the real record's phase F_ST beats (strictly).
        observed_fst : the real record's F_ST between phases.
        n_alt : as given.
        rows : one dict per leak with leak, records, around and compact (5th,
            25th, 50th, 75th, 95th percentiles of the records' scores),
            around_at_or_below and compact_at_or_below (share of records
            scoring no higher than the real record), reach_fst (share of
            records with F_ST at or above the observed).
        distinct_around, distinct_compact : distinct alternatives.
        same_as_phases_around, same_as_phases_compact : alternatives that
            are the phase division itself, which the phases cannot beat.

    Raises
    ------
    ValueError
        If `records` is out of range or a leak is not in `LEAKS`.
    """
    if not 30 <= records <= MAX_RUNS:
        raise ValueError(f"between 30 and {MAX_RUNS} simulated records per setting are needed")
    idx = data.phase_index
    around, compact = dv.alternatives(data.pts, idx, n_alt, np.random.default_rng(seed + 999_999))
    own = qs._canon(idx)

    def standing(m):
        """Share of each kind of alternative the phases beat on F_ST, and the phases' F_ST."""
        m = np.asarray(m, float)
        f = _fst(m, idx)
        return (float(np.mean([f > _fst(m, a) for a in around])),
                float(np.mean([f > _fst(m, c) for c in compact])), f)

    obs_a, obs_c, obs_f = standing(data.counts)
    rows = []
    for leak in leaks:
        if leak not in LEAKS:
            raise ValueError(f"leak must be one of {LEAKS}")
        sa, sc, sf = [], [], []
        for r in range(records):
            a, c, f = standing(simulate(data, rates, seed + 10000 * LEAKS.index(leak) + r, leak=leak))
            sa.append(a); sc.append(c); sf.append(f)
        sa, sc, sf = np.asarray(sa), np.asarray(sc), np.asarray(sf)
        rows.append({"leak": leak, "records": len(sa),
                     "around": [float(np.percentile(sa, q)) for q in (5, 25, 50, 75, 95)],
                     "compact": [float(np.percentile(sc, q)) for q in (5, 25, 50, 75, 95)],
                     # 1e-12: a score equal to the real one counts as "at or below" despite rounding
                     "around_at_or_below": float(np.mean(sa <= obs_a + 1e-12)),
                     "compact_at_or_below": float(np.mean(sc <= obs_c + 1e-12)),
                     "reach_fst": float(np.mean(sf >= obs_f))})
    return {"observed_around": obs_a, "observed_compact": obs_c, "observed_fst": obs_f,
            "n_alt": n_alt, "rows": rows,
            "distinct_around": len({qs._canon(a) for a in around}),
            "distinct_compact": len({qs._canon(c) for c in compact}),
            "same_as_phases_around": int(sum(qs._canon(a) == own for a in around)),
            "same_as_phases_compact": int(sum(qs._canon(c) == own for c in compact))}
