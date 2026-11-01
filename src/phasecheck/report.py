"""Run the questions in order and write a report a reader can follow.

`run` answers questions 1 to 5 (`questions.py`) and, with the model, 6 to 9
(`copying.py`), then writes, to the output directory:

  report.md             the report, one section per question, each with a
                        plain-language "Reading" derived from the numbers
  recovery.csv          question 1, one row per weight on composition
  assemblage_fit.csv    question 4, one row per assemblage
  results.json          every number behind the report
  calibration.csv       question 6, one row per combination (calibrated runs)
  model_comparison.csv  question 7, one row per measure (model runs)

The readings are chosen by thresholds on the numbers (see `_read_p`,
`_reading_7`, `_reading_9`); each threshold is stated where it is applied.
The boundary excess is printed with its sign: within-phase minus
between-phase similarity at matched distance, positive when the phase lines
separate assemblages more than distance alone predicts.

Seeds: questions 2 and 3 follow `seed`; question 4 uses `fit_seed` (96);
the model sections use fixed families far above the calibration's seeds
(see `run`), so they do not move with `seed`.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import copying as cp
from . import questions as qs
from .data import Dataset

SMALL_PHASE = 5          # below this many assemblages, a phase's own comparison has little power
LOW_COUNT = 75           # the paper's minimum for its own record; reported, not enforced


def _read_p(p: float) -> str:
    """Words for a probability that the tested lines beat the alternatives (bands at 0.05, 0.25, 0.75, 0.95)."""
    if p >= 0.95:
        return "better than nearly all of the alternatives"
    if p >= 0.75:
        return "better than most of the alternatives"
    if p > 0.25:
        return "about as well as the alternatives"
    if p > 0.05:
        return "worse than most of the alternatives"
    return "worse than nearly all of the alternatives"


def _w(w: float) -> str:
    """A weight on composition as words: 0 is the map alone, infinity composition alone."""
    return "composition alone" if np.isinf(w) else ("the map alone" if w == 0 else f"{w:g}")


def _comparison_table(r: dict) -> list:
    """Markdown lines of the question 2 table from a `questions.compare_division` result."""
    return [
        ("| measure | tested lines, posterior median (95% CrI) | random-center divisions | compact divisions | "
         "P(tested lines exceed random-center / compact) |"),
        "|---|---|---|---|---|",
        (f"| cultural F_ST | {r['fst'][0]:.4f} ({r['fst'][1]:.4f} to {r['fst'][2]:.4f}) | {r['fst_around'][0]:.4f} | "
         f"{r['fst_compact'][0]:.4f} | {r['p_fst_around']:.2f} / {r['p_fst_compact']:.2f} |"),
        (f"| boundary excess | {r['excess'][0]:+.1f} ({r['excess'][1]:+.1f} to {r['excess'][2]:+.1f}) | "
         f"{r['excess_around'][0]:+.1f} | {r['excess_compact'][0]:+.1f} | "
         f"{r['p_excess_around']:.2f} / {r['p_excess_compact']:.2f} |")]


def run(data: Dataset, out_dir, *, draws: int = 2000, n_alt: int = 300, fit_draws: int = 4000,
        seed: int = 0, fit_seed: int = 96, kmeans_seeds: int = 10, n_init: int = 500,
        source: str = "", model: bool = False, rates=None, length=None, jobs: int = 1,
        model_runs: int = 300, power_records: int = 300, progress=None) -> dict:
    """Run every question and write report.md, the tables as CSV, and results.json to `out_dir`.

    With `model=True` the copying model is calibrated to the record (or run at
    the supplied `rates`, a dictionary with n_ind, innovation, mixing and
    optionally length) and four more questions are answered. Calibration takes
    minutes to tens of minutes; `jobs` spreads it over processor cores.

    Parameters
    ----------
    data : Dataset
    out_dir : str or path-like
        Created if absent; files in it are overwritten.
    draws, n_alt, seed : int
        Posterior draws, alternative divisions of each kind, and base seed
        for questions 2 and 3 (see `questions.compare_division`).
    fit_draws, fit_seed : int
        Draws and seed for question 4.
    kmeans_seeds, n_init : int
        K-means seeds and random starts for question 1.
    source : str
        Input file name, printed in the report's first paragraph.
    model : bool
        Answer questions 6 to 9.
    rates : dict, optional
        Skip calibration and use these settings (n_ind, innovation, mixing,
        and optionally length in km).
    length : float, optional
        Copying length in km when `rates` gives none; default
        `copying.default_length`.
    jobs : int
        Worker processes for the calibration.
    model_runs : int
        Runs of the model for questions 7 and 8 (50 to 5,000).
    power_records : int
        Simulated records per setting for question 9 (30 to 5,000).
    progress : callable, optional
        Called with one-line progress messages.

    Returns
    -------
    dict
        recovery, boundary_comparison, each_phase_against_rest, fit,
        profiles: the results of questions 1 to 5; report: path of
        report.md.

    Raises
    ------
    ValueError
        From the model when calibration finds no setting that reproduces
        the record's diversity, when nearly every run loses all diversity, or
        when a run count is out of range.

    Seed families of the model sections, fixed and independent of `seed`:
    question 7 runs from 2,000,000 (and 2,100,000 and 2,200,000 for the
    middle and smallest confirmed settings), question 8 from 3,000,000 with
    sampling seeds 500,000 above, question 9 from 4,000,000. Each family
    holds at most 5,000 runs, so none overlaps another or the calibration.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    names = data.phase_names
    sizes = {p: int((data.phases == p).sum()) for p in names}
    totals = data.counts.sum(1)

    rec = qs.recovery_by_location(data, kmeans_seeds=kmeans_seeds, n_init=n_init)
    cmp_all = qs.boundary_comparison(data, draws=draws, n_alt=n_alt, seed=seed)
    rest = qs.each_phase_against_rest(data, draws=draws, n_alt=n_alt, seed=seed)
    fit = qs.own_phase_fit(data, draws=fit_draws, seed=fit_seed)
    prof = qs.phase_profiles(data)
    say = progress or (lambda *_: None)
    cal = cmp_model = sorting = pw = across = None
    if model:
        if rates is None:
            cal = cp.calibrate(data, length=length, jobs=jobs, progress=progress)
            rates = cal["rates"]
        else:
            rates = {"n_ind": int(rates["n_ind"]), "innovation": float(rates["innovation"]),
                     "mixing": float(rates["mixing"]),
                     "length": float(rates.get("length") or length or cp.default_length(data))}
        say("model: comparing the record with the copying model")
        # Seed families for this run. They sit far above the calibration's
        # (60,000 to about 100,000), so no run reuses a trajectory a setting was chosen on.
        cmp_model = cp.compare_with_model(data, rates, reps=model_runs, seed=2_000_000)
        across = [("the selected setting (largest difference between phases)" if cal is not None else
                   "the supplied setting", rates, cmp_model["rows"][0])]
        # Settings that fit equally well can disagree; with three or more confirmed, also show
        # the middle and the smallest of the ranked settings, each on a seed family of its own.
        if cal is not None and cal["n_confirmed"] >= 3:
            ranked = cal["ranked"]
            for label, j, sd in (("the middle confirmed setting", len(ranked) // 2, 2_100_000),
                                 ("the confirmed setting with the smallest difference", len(ranked) - 1, 2_200_000)):
                r = {"n_ind": int(ranked.n_ind[j]), "innovation": float(ranked.innovation[j]),
                     "mixing": float(ranked.mixing[j]), "length": rates["length"]}
                across.append((label, r, cp.compare_with_model(data, r, reps=model_runs, seed=sd)["rows"][0]))
        say("model: sorting simulated records into the phases")
        sorting = cp.sorting_under_model(data, rates, reps=model_runs, seed=3_000_000, sample_offset=500_000)
        say("model: estimating what the comparison can detect")
        pw = cp.power(data, rates, records=power_records, seed=4_000_000)

    L = ["# Do these phases behave as bounded groups?", "",
         f"Produced by phasecheck {_version()}" + (f" from `{source}`" if source else "") + ". "
         f"{len(data.names)} assemblages, {len(data.classes)} classes, {int(totals.sum()):,} sherds; phases "
         + ", ".join(f"{p} ({sizes[p]})" for p in names) + f". Distance: {data.dist_kind}. "
         f"Sherds per assemblage: {int(totals.min())} to {int(totals.max())} (median {int(np.median(totals))}). "
         f"Settings: {draws} posterior draws, {n_alt} alternative divisions of each kind, seed {seed}.", "",
         "Classes read as sherd counts (check that every one is a class, not a total, an identifier or a "
         "measurement): " + "; ".join(f"{_md(c)} ({int(t):,})" for c, t in zip(data.classes, data.counts.sum(0)))
         + ".", ""]

    # ---- what the reader must know first
    L += ["## Before reading the results", "",
          ("- These tests ask whether learning was bounded where the phase lines fall. They do not test "
           "whether a society or polity existed; political boundaries need not follow pottery."),
          ("- Sections 1 to 5 compare the phase lines with other lines on the same map. Sections 6 to 9 use "
           "a model of copying across distance with no groups, calibrated to this record, to ask how much "
           "difference distance alone produces here and what the comparisons can detect."
           if model else
           "- This run did not include the copying model (`--model`). It compares the phase lines with "
           "other lines on the same map. It does not say how much difference copying across distance alone "
           "would produce here, so a difference between phases is not by itself evidence of a boundary."),
          ("- Each probability is the phase lines' standing among alternative lines: the share of comparisons, "
           "over posterior draws of the counts, in which the phase lines separate the assemblages better. It "
           "is not a significance level. An arbitrary division of a map with any spatial pattern can land "
           "anywhere between 0 and 1, depending on how its lines happen to run. A value near one half says "
           "the phase lines are unremarkable. A value near 1 is necessary for the lines to stand apart but "
           "is not sufficient: arbitrary lines reach it too, most easily against compact divisions, which "
           "are few. The second decimal place changes from seed to seed."),
          ("- The alternative divisions and the spread of each division are built from the coordinates "
           "(straight lines). A supplied distance matrix enters only the boundary excess."),
          ("- How well these comparisons can detect a real boundary depends on the record. Section 9 "
           "measures it for this one." if model else
           "- How well these comparisons can detect a real boundary depends on the record. The paper found "
           "that for its 28 assemblages in three phases the comparison discriminates only moderately, and "
           "for a single line hardly at all. Power for this record has not been measured in this run.")]
    low = [(n, int(t)) for n, t in zip(data.names, totals) if t < LOW_COUNT]
    if low:
        L.append(f"- {len(low)} assemblage(s) have fewer than {LOW_COUNT} sherds, the minimum the paper "
                 "justified for its own record: " + ", ".join(f"{n} ({t})" for n, t in low[:12])
                 + (" ..." if len(low) > 12 else "") + ". Small collections differ by sampling alone; "
                 "consider `--min-count`.")
    if data.dropped:
        L.append("- Dropped for falling below the minimum count: "
                 + ", ".join(f"{n} ({t})" for n, t in data.dropped) + ".")
    for note in data.notes:
        L.append("- " + note)
    small = [p for p in names if sizes[p] < SMALL_PHASE]
    if small:
        L.append(f"- Phases with fewer than {SMALL_PHASE} assemblages ({', '.join(small)}) give their own "
                 "comparisons little to work with.")
    L.append("")

    # ---- 1
    L += ["## 1. Does location alone recover the phases?", "",
          (f"The assemblages are clustered into {rec['k']} groups on site coordinates and class composition "
           "together, as the weight on composition rises from the map alone to composition alone. Agreement "
           "with the phases is the adjusted Rand index (1 identical, about 0 unrelated). Divisions of the "
           f"phases' sizes around random centers agree with the phases at a median of {cmp_all['ari_around']:.2f}; "
           "values near that are what proximity alone gives any such division."), "",
          "| weight on composition | k-means | Ward linkage | average linkage | misplaced by Ward |",
          "|---|---:|---:|---:|---|"]
    for r in rec["rows"]:
        mis = r["ward_misplaced"]
        L.append(f"| {_w(r['weight'])} | {r['kmeans']:.3f} | {r['ward']:.3f} | {r['average']:.3f} | "
                 + (", ".join(mis) if len(mis) <= 8 else f"{len(mis)} assemblages") + " |")
    r0, rinf = rec["rows"][0], rec["rows"][-1]
    L += ["", (f"**Reading.** On the map alone, Ward linkage agrees with the phases at {r0['ward']:.2f} and "
          f"k-means at {r0['kmeans']:.2f}; on composition alone they give {rinf['ward']:.2f} and "
          f"{rinf['kmeans']:.2f}. If the map alone recovers the phases about as well as composition does, "
          "or better, the phases follow geography, and sorting by proximity would produce them whether or "
          "not groups existed. If phases were drawn using these same class frequencies, composition "
          "recovering them is expected by construction and is not evidence of a boundary."), ""]

    # ---- 2
    L += ["## 2. Do the phase boundaries separate the assemblages better than other lines?", "",
          ("The phases are compared with divisions of the same assemblages into groups of the same sizes: "
           "around randomly placed centers, and made compact by swapping pairs of assemblages between "
           f"groups until no swap tightens them. Of the {n_alt} drawn of each kind, "
           f"{cmp_all['distinct_around']} random-center divisions and {cmp_all['distinct_compact']} compact "
           "ones are distinct. Cultural F_ST "
           "measures how much the groups differ in class mix. The boundary excess is the similarity of "
           "pairs within a group minus that of pairs in different groups, at matched distance (positive: "
           "a step in similarity at the lines beyond what distance predicts)."), ""]
    L += _comparison_table(cmp_all)
    L += ["", f"**Reading.** On cultural F_ST the phase lines separate the assemblages "
          f"{_read_p(cmp_all['p_fst_around'])} around random centers ({cmp_all['p_fst_around']:.2f}) and "
          f"{_read_p(cmp_all['p_fst_compact'])} made compact ({cmp_all['p_fst_compact']:.2f}). "
          f"On the boundary excess the two probabilities are {cmp_all['p_excess_around']:.2f} and "
          f"{cmp_all['p_excess_compact']:.2f}. " + _excess_note(cmp_all) +
          f"The phases' own spread is {cmp_all['spread']:.1f} km² per assemblage; "
          f"{100 * cmp_all['share_around_tighter']:.0f} percent of the random-center divisions and "
          f"{100 * cmp_all['share_compact_tighter']:.0f} percent of the compact ones are tighter.", ""]

    # ---- 3
    if rest:
        L += ["## 3. Does each phase stand apart from the rest?", "",
              "Each phase against all the others taken together, beside divisions of the same two sizes.", "",
              ("| phase | assemblages | F_ST, posterior median (95% CrI) | P(F_ST exceeds random-center / compact) "
               "| boundary excess | P(excess exceeds random-center / compact) | distinct alternatives, "
               "random-center / compact |"), "|---|---:|---|---|---|---|---|"]
        for p in names:
            r = rest[p]
            L.append(f"| {_md(p)} | {sizes[p]} | {r['fst'][0]:.4f} ({r['fst'][1]:.4f} to {r['fst'][2]:.4f}) | "
                     f"{r['p_fst_around']:.2f} / {r['p_fst_compact']:.2f} | {r['excess'][0]:+.1f}"
                     f"{'*' if r['excess_bins'] == 0 else ''} | "
                     f"{r['p_excess_around']:.2f} / {r['p_excess_compact']:.2f} | "
                     f"{r['distinct_around']} / {r['distinct_compact']} |")
        L += ["", "**Reading.** A single line gives these comparisons little power. Compact two-way "
              "divisions of a map are few and alike (see the last column), so a probability against them, "
              "even 1.00, can come from arbitrary lines and carries less weight than the one against "
              "random centers."
              + (" An asterisk marks a boundary excess with no distance control, because no distance bin "
                 "held enough pairs." if any(rest[p]["excess_bins"] == 0 for p in names) else ""), ""]
    else:
        L += ["## 3. Does each phase stand apart from the rest?", "",
              "With two phases, each phase against the rest is the division already tested in section 2.", ""]

    # ---- 4
    clear = [r for r in fit if r["p_own"] > 0.75]
    other = [r for r in fit if r["likeliest"] != r["phase"]]
    unclear = [r for r in fit if r["likeliest"] == r["phase"] and r["p_own"] <= 0.75]
    L += ["## 4. Does each assemblage fit its own phase?", "",
          ("Each assemblage's class profile is compared with the pooled profile of every phase, the "
           "assemblage left out of its own phase's pool. P(own) is the share of posterior draws in which "
           "its own phase is the closest."), "",
          f"- {len(clear)} of {len(fit)} fit their own phase best with probability above 0.75.",
          f"- {len(other)} fit another phase more often than their own: "
          + ("; ".join(f"{_md(r['name'])} ({_md(r['phase'])}; {_md(r['likeliest'])} in {r['p_likeliest']:.2f} "
                       "of draws)" for r in other) or "none") + ".",
          f"- {len(unclear)} fit their own phase most often but not clearly (0.75 or less): "
          + ("; ".join(_md(r["name"]) for r in unclear) or "none") + ".", "",
          ("| assemblage | sherds | phase | P(own) | likeliest phase (share of draws) | phase of the nearest "
           "other assemblage |"), "|---|---:|---|---:|---|---|"]
    for r in sorted(fit, key=lambda r: r["p_own"]):
        L.append(f"| {_md(r['name'])} | {r['sherds']} | {_md(r['phase'])} | {r['p_own']:.2f} | "
                 f"{_md(r['likeliest'])} ({r['p_likeliest']:.2f}) | {_md(r['nearest_phase'])} |")
    L.append("")

    # ---- 5
    L += ["## 5. How different are the phases, in sherds?", "",
          ("Each phase's sherds pooled, in percent by class, and the share of sherds that would have to "
           "change class to make two phases' profiles identical (0 identical, 100 nothing in common)."), "",
          "| class | " + " | ".join(map(_md, names)) + " | all |", "|---|" + "---:|" * (len(names) + 1)]
    for j in np.argsort(-prof["overall"]):
        L.append(f"| {_md(data.classes[j])} | " + " | ".join(f"{prof['profiles'][i, j]:.1f}" for i in range(len(names)))
                 + f" | {prof['overall'][j]:.1f} |")
    L += ["", "| pair of phases | sherds that would have to change class (percent) |", "|---|---:|"]
    for (a, b), v in prof["pairs"].items():
        L.append(f"| {_md(a)} and {_md(b)} | {v:.1f} |")
    L += [f"| mean over pairs | {prof['mean']:.1f} |", "",
          ("**Reading.** These differences are descriptions. Without a model of copying across this map "
           "they cannot be called large or small: nearby places differ by distance alone, and small "
           "collections differ by sampling alone."), ""]

    if model:
        L += _model_sections(data, cal, rates, cmp_model, sorting, pw, across)
        if cal is not None:
            cal["table"].to_csv(out / "calibration.csv", index=False)
        pd.DataFrame([{k: v for k, v in r.items() if k != "fmt"} for r in cmp_model["rows"]]
                     ).to_csv(out / "model_comparison.csv", index=False)

    (out / "report.md").write_text("\n".join(L), encoding="utf-8")
    pd.DataFrame([{k: v for k, v in r.items() if not k.endswith("_misplaced")} for r in rec["rows"]]
                 ).to_csv(out / "recovery.csv", index=False)
    pd.DataFrame([dict(name=r["name"], phase=r["phase"], sherds=r["sherds"], closest=r["closest"],
                       p_own=r["p_own"], likeliest=r["likeliest"], p_likeliest=r["p_likeliest"],
                       nearest_phase=r["nearest_phase"],
                       **{f"distance_to_{q}": d for q, d in r["distance"].items()}) for r in fit]
                 ).to_csv(out / "assemblage_fit.csv", index=False)
    results = {"version": _version(), "n_assemblages": len(data.names), "phases": sizes,
               "distance": data.dist_kind,
               "classes": list(map(str, data.classes)),
               "settings": {"draws": draws, "n_alt": n_alt, "fit_draws": fit_draws, "seed": seed,
                            "fit_seed": fit_seed},
               "boundary_comparison": cmp_all, "each_phase_against_rest": rest,
               "profile_difference": {f"{a} | {b}": v for (a, b), v in prof["pairs"].items()},
               "profile_difference_mean": prof["mean"],
               "model": None if not model else {
                   "rates": rates, "calibrated": cal is not None,
                   "n_matched": None if cal is None else cal["n_matched"],
                   "n_confirmed": None if cal is None else cal["n_confirmed"],
                   "comparison": [{k: v for k, v in r.items() if k != "fmt"} for r in cmp_model["rows"]],
                   "distance_decay": cmp_model["decay"], "sorting": sorting, "power": pw,
                   "across_settings": [dict(setting=a, rates=b, **{k: v for k, v in c.items() if k != "fmt"})
                                       for a, b, c in across]}}
    (out / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return {"recovery": rec, "boundary_comparison": cmp_all, "each_phase_against_rest": rest, "fit": fit,
            "profiles": prof, "report": str(out / "report.md")}


def _share(e: int, n: int) -> str:
    """'e of n (pct%)', with one decimal below 1 percent so a small count does not read as 0%."""
    if e == 0:
        return f"none of {n}"
    pct = 100 * e / n
    return f"{e} of {n} ({pct:.1f}%)" if pct < 1 else f"{e} of {n} ({pct:.0f}%)"


def _place(row: dict) -> str:
    """Where the observed value sits in the model's runs: above, inside or below their 95 percent range."""
    if row["observed"] > row["hi"]:
        return "above"
    if row["observed"] < row["lo"]:
        return "below"
    return "inside"


def _pct(share: float) -> str:
    """A share as a percentage that does not round a small non-zero value to 0."""
    v = 100 * share
    return f"{v:.1f}" if 0 < v < 1 else f"{v:.0f}"


def _reading_7(across: list, calibrated: bool) -> str:
    """The reading of question 7, from where the observed between-phase F_ST sits at each setting tried."""
    places = [_place(r) for _, _, r in across]
    top = across[0][2]
    many = len(across) > 1
    if all(p == "above" for p in places):
        out = ("The phases differ by more than the model with no groups produces"
               + (", at its most favorable setting and at the others tried" if many else " at this setting")
               + f" (runs reaching the observed value: {_share(top['reaching'], top['runs'])}). That is a "
               "difference to explain. It is not by itself evidence of a boundary at the phase lines: "
               "variation arising locally, differences of age, or differences at other lines would produce "
               "it too, and sections 2 and 3 ask whether it is tied to these lines.")
    elif all(p == "below" for p in places):
        out = ("The phases differ by LESS than the model with no groups produces at "
               + ("every setting tried" if many else "this setting")
               + ". The model does not fit this record: it makes assemblages more different from place to "
               "place than they are. "
               + ("A copying length that is too short for this map is one possible cause; rerun with a "
                  "larger `--length-km`. " if calibrated else
                  "The supplied settings, or a copying length too short for this map, may be the cause. ")
               + "Sections 8 and 9 rest on the same settings and should not be relied on until the model "
               "fits.")
    elif not many:                                   # one setting, and the record lies inside its range
        out = ("The observed difference between the phases lies inside the range of the model with no "
               f"groups (runs reaching it: {_share(top['reaching'], top['runs'])}). Copying across "
               "distance, as modeled, can produce a difference of this size without any groups.")
    elif "above" not in places:
        out = ("At no setting tried do the phases differ by more than the model with no groups produces"
               + (f" (runs reaching the observed value at the selected setting: "
                  f"{_share(top['reaching'], top['runs'])})" if places[0] == "inside" else "")
               + ". Copying across distance, as modeled, can produce a difference of this size without any "
               "groups."
               + (" At some settings the model produces more difference than the record has, so those "
                  "settings do not fit it." if "below" in places else ""))
    else:
        out = ("The answer depends on the setting. Settings that reproduce this record's diversity equally "
               "well put the observed difference " + ", ".join(
                   f"{p} the model's range at {lab}" for (lab, _, _), p in zip(across, places))
               + ". Diversity alone does not pin the model down for this record, so it cannot say whether "
               "the phases differ by more than copying across distance produces.")
    if not calibrated:
        out += (" These settings were supplied and were not checked against the record's diversity, so this "
                "reading holds only if they fit it.")
    return out


def _model_sections(data, cal, rates, cmp_model, sorting, pw, across) -> list:
    """Sections 6 to 9: the copying model with no groups, calibrated to this record."""
    obs = cp.diversity(data.counts)
    L = ["## 6. Can copying with no groups reproduce this record's diversity?", "",
         ("The copying model has no groups. Learners copy within their own assemblage's population, "
          "sometimes from other assemblages (nearer ones more often), and sometimes adopt a new variant. "
          "Its settings are first required to reproduce the record's diversity, which does not involve the "
          f"phases: mean diversity within an assemblage ({obs['within']:.3f}), mean number of classes "
          f"present ({obs['classes']:.1f}), and diversity of the whole record pooled ({obs['pooled']:.3f})."), ""]
    if cal is not None:
        L += [(f"- {cal['n_matched']} of {cal['n_cells']} combinations of settings matched all three within "
               f"{100 * cal['tolerance']:.0f} percent on {cal['reps']} runs, and {cal['n_confirmed']} still "
               f"matched on {cal['confirm_reps']} fresh runs."),
              ("- Among the combinations that fit, the one whose runs give the phases the largest difference "
               f"was selected: {rates['n_ind']:,} learners per assemblage, innovation rate "
               f"{rates['innovation']:g}, mixing rate {rates['mixing']:g}. This choice leans toward the "
               "conclusion that copying with no groups accounts for the phases, and makes it harder to find "
               "that they differ by more. It is the most favorable setting for that one comparison only.")]
        if cal["n_confirmed"] < 5:
            L.append(f"- Only {cal['n_confirmed']} combination(s) were confirmed, so the comparison rests "
                     "on a narrow part of the model's settings.")
        if cal["n_confirmed"] > cal["n_cells"] / 4:
            L.append(f"- **{cal['n_confirmed']} of {cal['n_cells']} combinations fit, so diversity "
                     "constrains the model only weakly for this record.** Settings that fit equally well "
                     "can differ widely in what they predict for the phases; section 7 shows how widely.")
    else:
        L += [(f"- Settings were supplied, not calibrated in this run: {rates['n_ind']:,} learners per "
               f"assemblage, innovation rate {rates['innovation']:g}, mixing rate {rates['mixing']:g}. "
               "Whether they reproduce the record's diversity was not checked here.")]
    L += [(f"- Copying falls off with distance over a length of {rates['length']:.1f} km ({data.dist_kind} "
           "distance). **This length is an assumption, not an estimate.** Rerun with `--length-km` at half "
           "and double this value to see what depends on it."),
          "- Assemblages are treated as "
          + ("positioned along the supplied order, each sampled from its own part of the simulated "
             "sequence. Only the order of the values is used, not their spacing."
             if data.order is not None else
             "contemporaneous. If they span a long sequence, supply an order column (`--order-col`)."),
          ("- Also assumed, and not varied here: every assemblage has the same number of learners and the "
           "same total rate of copying from outside; copying falls off exponentially with distance; new "
           "variants arise everywhere in proportion to the record's pooled class frequencies; each deposit "
           "pools the same span of time. The paper relaxes several of these for its own record."), ""]

    L += ["## 7. Does the record differ between phases by more than copying across distance produces?", "",
          (f"The record beside {cmp_model['runs']} runs of the model at the selected setting, each sampled "
           "at the real sherd counts. The last column counts the runs that reach the observed value."), "",
          "| measure | observed | model median (95 percent range) | runs reaching observed |",
          "|---|---:|---|---:|"]
    for r in cmp_model["rows"]:
        f = r["fmt"]
        L.append(f"| {_md(r['measure'])} | {f.format(r['observed'])} | {f.format(r['median'])} "
                 f"({f.format(r['lo'])} to {f.format(r['hi'])}) | {_share(r['reaching'], r['runs'])} |")
    if cmp_model["flat"]:
        L += ["", (f"In {cmp_model['flat']} of {cmp_model['runs']} runs the simulated record held a single "
              "class. Those runs count as showing no difference between phases; they are left out of the "
              "similarity-and-distance line below, which cannot be computed for them.")]
    if len(across) > 1:
        L += ["", "The between-phase F_ST at settings that fit the record's diversity equally well:", "",
              ("| setting | learners, innovation, mixing | model median (95 percent range) | runs reaching "
               "observed | observed lies |"), "|---|---|---|---:|---|"]
        for lab, r, row in across:
            L.append(f"| {lab} | {r['n_ind']:,}, {r['innovation']:g}, {r['mixing']:g} | {row['median']:.4f} "
                     f"({row['lo']:.4f} to {row['hi']:.4f}) | {_share(row['reaching'], row['runs'])} | "
                     f"{_place(row)} the range |")
    d = cmp_model["decay"]
    if d["observed"] is None:
        L += ["", ("Similarity and distance: every pair of assemblages is equally similar in this record, so "
              "the correlation between similarity and distance cannot be computed.")]
    else:
        L += ["", _decay_line(d)]
    # 3 / runs: the rule of three, the upper 95 percent bound on a frequency never observed.
    L += ["", "**Reading.** " + _reading_7(across, cal is not None)
          + f" Where no run reaches a value, its frequency in the model is not zero; with "
          f"{cmp_model['runs']} runs it can still be as high as about {_pct(3 / cmp_model['runs'])} "
          "percent.", ""]
    L += ["## 8. Would assemblages made with no groups sort into these phases?", "",
          (f"Each of {sorting['runs']} simulated records, at the selected setting, is clustered by Ward "
           "linkage as the real record was in section 1, and scored against the phases. The map alone is "
           "the same for every record."), "",
          ("| weight on composition | real record | model median (95 percent range) | runs reaching the "
           "real record |"), "|---|---:|---|---:|"]
    for r in sorting["rows"]:
        L.append(f"| {_w(r['weight'])} | {r['observed']:.3f} | {r['median']:.3f} ({r['lo']:.3f} to "
                 f"{r['hi']:.3f}) | {_share(r['reaching'], sorting['runs'])} |")
    L += ["", ("**Reading.** Where most runs reach the real record's agreement, sorting assemblages made "
          "with no groups gives these phases about as readily as the real ones. Where few do, the real "
          "pottery lines up with the phases better than simulated pottery does; if the phases were drawn "
          "from these same class frequencies, part of that is expected by construction."), ""]

    free = pw["rows"][0]                    # leak 1, the unrestricted setting, in the default order
    L += ["## 9. Could this comparison detect a restriction at the phase boundaries?", "",
          (f"Records were simulated at the selected setting with copying across the phase lines "
           f"unrestricted and restricted, {free['records']} per setting, and each was scored as the real "
           f"record is: the share of {pw['n_alt']} alternative divisions of each kind that the phase lines "
           f"beat on cultural F_ST. The real record's own scores are {pw['observed_around']:.2f} "
           f"(random-center) and {pw['observed_compact']:.2f} (compact). These are direct scores on the "
           "counts, against a smaller set of alternatives than section 2 uses, so they differ somewhat from "
           f"the probabilities there. Of the {pw['n_alt']} alternatives of each kind, "
           f"{pw['distinct_around']} random-center and {pw['distinct_compact']} compact ones are distinct, "
           f"and {pw['same_as_phases_around']} and {pw['same_as_phases_compact']} are the phases "
           "themselves, which the phases cannot beat."), "",
          ("| copying across the phase lines | score vs random-center, 5% / 25% / 50% / 75% / 95% of records "
           "| records at or below the real score | score vs compact, 5% / 25% / 50% / 75% / 95% | records at "
           "or below the real score | records reaching the observed F_ST |"), "|---|---|---:|---|---:|---:|"]
    for r in pw["rows"]:
        lab = "unrestricted" if r["leak"] == 1 else f"cut to {100 * r['leak']:g} percent"
        L.append(f"| {lab} | " + " / ".join(f"{v:.2f}" for v in r["around"]) +
                 f" | {_pct(r['around_at_or_below'])}% | " + " / ".join(f"{v:.2f}" for v in r["compact"]) +
                 f" | {_pct(r['compact_at_or_below'])}% | {_pct(r['reach_fst'])}% |")
    L += ["", "**Reading.** " + _reading_9(pw), ""]
    return L


def _decay_line(d: dict) -> str:
    """Similarity against distance, with the direction read from the sign (negative: similarity falls)."""
    return ("Similarity and distance: the correlation between pairwise similarity and distance is "
            f"{d['observed']:+.2f} in the record, so similarity "
            + ("falls with distance" if d["observed"] < -0.1 else
               "rises with distance" if d["observed"] > 0.1 else "hardly changes with distance")
            + f". In the model it is {d['median']:+.2f} ({d['lo']:+.2f} to {d['hi']:+.2f}), and the record's "
            "value lies " + ("inside" if d["lo"] <= d["observed"] <= d["hi"] else "outside") + " that range.")


def _reading_9(pw: dict) -> str:
    """The reading of question 9, from how far a restriction moves the phase lines' score.

    The gap is the median score against random-center divisions with the
    strongest restriction minus that with none (restricted minus free).
    Up to 0.05 the comparison is said not to detect a restriction; below
    0.15 it responds slightly, below 0.3 moderately, otherwise strongly. The
    real score is then weighed by how much more common a score that low is
    under one kind of record than the other (a ratio of 1.25 or more counts,
    2 or more without "weakly"), with zero counts floored at 3/n.
    """
    free = max(pw["rows"], key=lambda r: r["leak"])        # no restriction
    tight = min(pw["rows"], key=lambda r: r["leak"])       # the strongest restriction
    n = free["records"]
    gap = round(tight["around"][2] - free["around"][2], 2)
    med = f"(median {free['around'][2]:.2f} without one, {tight['around'][2]:.2f} with the strongest)"
    if gap <= 0.05:
        out = (f"A restriction does not raise the phase lines' score on this record {med}. **The comparison "
               "with alternative lines cannot detect a restriction here**, and sections 2 and 3 say nothing "
               "about one either way.")
    else:
        how = "slightly" if gap < 0.15 else "moderately" if gap < 0.3 else "strongly"
        out = (f"A restriction raises the phase lines' score against random-center divisions {how} {med}, "
               "so the comparison responds to a restriction on this record"
               + (", though weakly." if gap < 0.15 else "."))
    a, b = free["around_at_or_below"], tight["around_at_or_below"]
    floor = 3.0 / n                                  # what "none of n" still allows (rule of three)
    out += (f" A score as low as the real one occurs in {_pct(a)} percent of records simulated with no "
            f"restriction and {_pct(b)} percent with the strongest.")
    against, for_ = max(a, floor) / max(b, floor), max(b, floor) / max(a, floor)
    bound = " (a lower bound, since none of the {} records scored that low)"
    if a == 0 and b == 0:
        out += (" Every simulated record of both kinds scored above the real one, so the real score is "
                "lower than either kind of record produces and favors neither.")
    elif against >= 1.25:
        out += (f" It is about {against:.1f} times as common without a restriction, which weighs "
                + ("" if against >= 2 else "weakly ") + "against one"
                + (bound.format("restricted") if b == 0 else "") + ".")
    elif for_ >= 1.25:
        out += (f" It is about {for_:.1f} times as common with a restriction, which weighs "
                + ("" if for_ >= 2 else "weakly ") + "for one"
                + (bound.format("unrestricted") if a == 0 else "") + ".")
    else:
        out += " The two are about equally common, so the real score does not favor either."
    out += (f" The compact divisions give {_pct(free['compact_at_or_below'])} and "
            f"{_pct(tight['compact_at_or_below'])} percent"
            + (", but so few of them are distinct that this column carries little weight"
               if pw["distinct_compact"] < 5 else "")
            + f". These figures rest on {n} records per setting and one set of alternative divisions; they "
            "change by several hundredths with the random seed, so read them as rough.")
    return out


def _md(text) -> str:
    """A name made safe for a Markdown table cell."""
    return str(text).replace("|", "/").replace("\n", " ")


def _excess_note(r: dict) -> str:
    """The sentence on how many distance bins the boundary excess rests on (the default 4 bins)."""
    if r["excess_bins"] == 0:
        return ("**No distance bin held enough pairs of each kind, so the boundary excess here is the plain "
                "gap between within-group and between-group similarity, with no control for distance.** ")
    note = f"The boundary excess rests on {r['excess_bins']} of 4 distance bins. "
    if r["excess_uncontrolled"] > 0:
        note += (f"In {100 * r['excess_uncontrolled']:.0f} percent of draws no bin qualified and the "
                 "distance control was dropped. ")
    return note


def _version() -> str:
    """The package version (imported late, since the package imports this module)."""
    from . import __version__
    return __version__
