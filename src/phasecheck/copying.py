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

N_IND = (120, 2000, 10000)
INNOVATIONS = (.0002, .0005, .001, .002, .004, .008, .012, .024, .048, .1, .2)
MIXINGS = (.001, .002, .005, .01, .02, .05, .1, .2, .4)
LEAKS = (1.0, 0.5, 0.1, 0.03)          # copying across a phase line, as a share of its strength
TOLERANCE = 0.10
LENGTH_SHARE = 0.2                     # default copying length, as a share of the largest distance


def default_length(data: Dataset) -> float:
    """The copying length used when none is given: a fifth of the largest distance between assemblages.

    This is a choice, not an estimate. The paper used 24 km on a basin whose
    assemblages lie up to 268 km apart along the rivers (79 km in a straight
    line), so this default is not the paper's value for the paper's record.
    """
    return float(LENGTH_SHARE * data.dist.max())


def diversity(m: np.ndarray) -> dict:
    """Mean diversity within assemblages, mean number of classes present, and diversity of the pooled record."""
    m = np.asarray(m, float)
    p = m / m.sum(1, keepdims=True)
    pooled = m.sum(0) / m.sum()
    return dict(within=float((1 - (p ** 2).sum(1)).mean()), classes=float((m > 0).sum(1).mean()),
                pooled=float(1 - (pooled ** 2).sum()))


MAX_RUNS = 5000          # seed families are spaced 10,000 apart


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
    return np.ones(len(data.names)) if data.order is None else data.order


def simulate(data: Dataset, rates: dict, seed: int, *, leak: float = 1.0, sample_seed=None) -> np.ndarray:
    """One simulated count table at the real sherd counts. `leak` below 1 restricts copying across phase lines."""
    w = md.copying_weights(data.dist, rates["length"], data.phase_index if leak < 1 else None, leak)
    rec = md.drift_record(w, k=data.counts.shape[1], n_ind=rates["n_ind"], mixing=rates["mixing"],
                          innovation=rates["innovation"], target=data.counts.sum(0) / data.counts.sum(),
                          seed=seed)
    rng = np.random.default_rng(seed + 5_000_000 if sample_seed is None else sample_seed)
    return md.sample_record(rec, _order(data), data.counts.sum(1), rng)


# ------------------------------------------------------------------ calibration
def _cell(args):
    """Score one combination of settings: diversity of its runs and the phases' F_ST in them."""
    data, n_ind, innovation, mixing, length, seeds, sample_seeds, labels = args
    rates = dict(n_ind=n_ind, innovation=innovation, mixing=mixing, length=length)
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
        if jobs > 1:
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                return list(ex.map(_cell, jobs_args, chunksize=4))
        return [_cell(a) for a in jobs_args]

    say(f"calibration: screening {len(cells)} combinations, {reps} runs each")
    screen = run([(data, n, i, x, length, [60000 + s for s in range(reps)],
                   [61000 + s for s in range(reps)], labels) for n, i, x in cells])
    rows = []
    for (n, i, x), (mean, fst) in zip(cells, screen):
        rel = {k: abs(mean[k] - obs[k]) / obs[k] for k in obs}
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
    best = ok.loc[ok.fst_median.idxmax()]
    frame["selected"] = frame.index == best.name
    rates = dict(n_ind=int(best.n_ind), innovation=float(best.innovation), mixing=float(best.mixing),
                 length=length)
    return dict(rates=rates, table=frame, observed=obs, n_cells=len(cells), n_matched=int(frame.matched.sum()),
                n_confirmed=int(len(ok)), reps=reps, confirm_reps=confirm_reps, tolerance=tolerance,
                ranked=ok.sort_values("fst_median", ascending=False, kind="stable")[["n_ind", "innovation", "mixing",
                                                                     "fst_median"]].reset_index(drop=True))


# ------------------------------------------------------------------ the record against the model
def _row(label, observed, sims, fmt):
    a = np.asarray(sims, float)
    a = a[np.isfinite(a)]
    return dict(measure=label, observed=float(observed), median=float(np.median(a)),
                lo=float(np.percentile(a, 2.5)), hi=float(np.percentile(a, 97.5)),
                reaching=int((a >= observed).sum()), runs=int(len(a)), fmt=fmt)


def compare_with_model(data: Dataset, rates: dict, *, reps: int = 300, seed: int = 84000) -> dict:
    """The record beside `reps` runs of the copying model with no groups.

    For each measure: the observed value, the model's median and 95 percent
    range, and the number of runs that reach the observed value. A measure the
    runs often reach is one the model accounts for. A measure they rarely
    reach is a difference copying across distance, as modeled, does not
    produce; it is not by itself evidence of a boundary.
    """
    if not 50 <= reps <= MAX_RUNS:
        raise ValueError(f"between 50 and {MAX_RUNS} runs are needed to describe the model's range")
    idx, names, counts = data.phase_index, data.phase_names, data.counts.astype(float)
    iu = np.triu_indices(len(counts), 1)
    def decay(m):
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
    return dict(rows=rows, runs=reps, flat=flat, rates=rates,
                decay=dict(observed=obs_decay if np.isfinite(obs_decay) else None, median=float(np.median(dec)), lo=float(np.percentile(dec, 2.5)),
                           hi=float(np.percentile(dec, 97.5))))


def sorting_under_model(data: Dataset, rates: dict, *, reps: int = 300, seed: int = 79000,
                        sample_offset: int = 500, weights=(0.0, 0.25, 0.5, 1.0, np.inf)) -> dict:
    """Would assemblages made by copying with no groups sort into the phases?

    Each simulated record is clustered by Ward linkage exactly as the real one
    is, at each weight on composition, and scored against the phases. The map
    alone is the same for every record, so its row cannot differ.
    """
    idx, k = data.phase_index, len(data.phase_names)
    g = qs._unit_scale(data.pts)

    def recover(m):
        c = qs._unit_scale(qs._chisq_features(np.asarray(m, float)))
        out = []
        for w in weights:
            f = g if w == 0 else (c if np.isinf(w) else np.column_stack([g, np.sqrt(w) * c]))
            out.append(ms.adjusted_rand(idx, qs._cluster(f, k, "ward", 0, 0)))
        return out

    if not 50 <= reps <= MAX_RUNS:
        raise ValueError(f"between 50 and {MAX_RUNS} runs are needed")
    if sample_offset < reps:
        raise ValueError("more than 500 runs need a larger sample_offset, or run seeds and sample seeds overlap")
    obs = recover(data.counts)
    sims = [recover(simulate(data, rates, seed + r, sample_seed=seed + sample_offset + r)) for r in range(reps)]
    arr = np.asarray(sims)
    return dict(runs=len(sims), rows=[dict(weight=float(w), observed=obs[j], median=float(np.median(arr[:, j])),
                                           lo=float(np.percentile(arr[:, j], 2.5)),
                                           hi=float(np.percentile(arr[:, j], 97.5)),
                                           reaching=int((arr[:, j] >= obs[j] - 1e-12).sum()))
                                      for j, w in enumerate(weights)])


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
    """
    if not 30 <= records <= MAX_RUNS:
        raise ValueError(f"between 30 and {MAX_RUNS} simulated records per setting are needed")
    idx = data.phase_index
    around, compact = dv.alternatives(data.pts, idx, n_alt, np.random.default_rng(seed + 999_999))
    own = qs._canon(idx)

    def standing(m):
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
        rows.append(dict(leak=leak, records=len(sa),
                         around=[float(np.percentile(sa, q)) for q in (5, 25, 50, 75, 95)],
                         compact=[float(np.percentile(sc, q)) for q in (5, 25, 50, 75, 95)],
                         around_at_or_below=float(np.mean(sa <= obs_a + 1e-12)),
                         compact_at_or_below=float(np.mean(sc <= obs_c + 1e-12)),
                         reach_fst=float(np.mean(sf >= obs_f))))
    return dict(observed_around=obs_a, observed_compact=obs_c, observed_fst=obs_f, n_alt=n_alt, rows=rows,
                distinct_around=len({qs._canon(a) for a in around}),
                distinct_compact=len({qs._canon(c) for c in compact}),
                same_as_phases_around=int(sum(qs._canon(a) == own for a in around)),
                same_as_phases_compact=int(sum(qs._canon(c) == own for c in compact)))
