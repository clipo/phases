"""89_partition_power.py - could the phase-line comparison see a copying boundary if there were one?

The main text argues from analysis 86 that the phase lines separate the pottery
no better than arbitrary lines of the same sizes. That argument has force only
if the comparison can detect lines that DO mark restricted copying. This script
measures its power (a rule 20c-style recovery check, 2026-09-23). It rebuilds
simulated records of analysis 84's grid, with the same seeds, at cells with and
without a copying boundary at the phase lines, and asks of each simulated
record the question 86 asks of the real one: what share of same-size
alternative divisions (around random centers, and made compact) do the phase
lines beat, on between-group F_ST and on the boundary excess? The same is asked
of the Parkin-versus-rest line, which is part of the phase boundary and so is
restricted whenever the phase lines are.

Plug-in statistics on each simulated record: the simulated counts are the truth
here, so no posterior draw of the counts is needed. The observed record's value
of each share is reported alongside (plug-in, 86's alternatives). The simulated
records are drawn at the real sherd counts, so each simulated plug-in score
already carries the sampling noise the observed plug-in score carries; that is
why the observed plug-in value, not 86's posterior probability, is the one the
simulated distributions are read against.

MEDIANS ARE NOT ENOUGH (review of 2026-09-24, item 4). A rise in the median share
from the no-boundary row to the boundary rows shows the comparison RESPONDS to
restricted copying; it does not say how often a boundary-bearing record scores
as low as the observed one. Every simulated record's four shares are therefore
written out (RECORDS_CSV), and the report gives, for each cell, the quantiles of
each share and the fraction of simulated records scoring at or below the
observed value. The whole grid is also rerun at two further diversity-matched
combinations (COMBOS), because a power statement made at one calibrated
combination is a statement about that combination. The two are fixed here and
checked against the calibration at run time: the script raises if either no
longer passes the fifty-run diversity confirmation.

Output: output/findings/partition_power.md, output/findings/partition_power_records.csv,
figures/figS12_partition_power.*
Usage: .venv/bin/python analyses/89_partition_power.py [--reps 100] [--alt 50]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
OUT_MD = ROOT / "output" / "findings" / "partition_power.md"
RECORDS_CSV = ROOT / "output" / "findings" / "partition_power_records.csv"
FIG = "figS12_partition_power"
# Further diversity-matched combinations (learners, innovation, mixing). They were
# picked from the confirmed set of the 2026-09-23 calibration; in the calibration of
# 2026-10-01 both still pass the fifty-run confirmation (fourth and eighth by
# fifty-run median), which main() verifies. The calibrated combination is read
# from calibrated_rates.csv.
COMBOS = [(10000, 0.0005, 0.002), (2000, 0.002, 0.005)]
STATS = ("f_alt", "f_opt", "b_alt", "b_opt")
CELLS = [(1.0, 0.0), (0.1, 0.0), (0.03, 0.0), (1.0, 0.2), (0.03, 0.2)]  # (copying factor, local innovation)


def alternatives(t74, pts, labels, n_alt, seed):
    sizes = np.bincount(labels); k = len(sizes); slot = np.repeat(np.arange(k), sizes)
    rng = np.random.default_rng(seed)
    alt, opt = [], []
    for _ in range(n_alt):
        a = t74.assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)
        alt.append(a); opt.append(t74.local_search(pts, a, sizes))
    return alt, opt


def shares(m, labels, alt, opt, rev, sd, d):
    f = rev.fst_by(m, labels); b = sd.boundary_excess_labeled(m, d, labels)
    fa = np.array([rev.fst_by(m, a) for a in alt]); fo = np.array([rev.fst_by(m, a) for a in opt])
    ba = np.array([sd.boundary_excess_labeled(m, d, a) for a in alt])
    bo = np.array([sd.boundary_excess_labeled(m, d, a) for a in opt])
    return dict(f_alt=np.mean(f > fa), f_opt=np.mean(f > fo), b_alt=np.mean(b > ba), b_opt=np.mean(b > bo))


LABEL = {"f_alt": "F_ST vs random-center", "f_opt": "F_ST vs compact",
         "b_alt": "boundary excess vs random-center", "b_opt": "boundary excess vs compact"}


def distribution_section(df, obs, combos):
    """Quantiles of every share, and the fraction of records at or below the observed score.

    The at-or-below fraction is what a median cannot give: how often a record
    generated WITH a given boundary scores as low as the real one. A ratio of two
    such fractions (boundary cell over no-boundary cell) is the likelihood ratio
    the observed score carries between them, under this simulator and this
    combination; it is reported as a description of the simulated distributions,
    not as a test.
    """
    L = ["## The distributions behind the medians", "",
         "For each cell, the 5th, 25th, 50th, 75th and 95th percentiles of the share of alternatives the "
         "tested line beats, over simulated records, and the fraction of simulated records whose share is "
         "at or below the observed record's. A fraction of 0 means none of the simulated records scored "
         "that low, which bounds the frequency below about 1/reps, not at zero.", ""]
    for ci, (n_ind, innov, mix) in enumerate(combos):
        sub = df[(df.n_ind == n_ind) & (df.innovation == innov) & (df.mixing == mix)]
        tag = "calibrated combination" if ci == 0 else "sensitivity combination"
        L += [f"### {n_ind} learners, innovation {innov}, mixing {mix} ({tag})", "",
              "| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |",
              "|---|---|---:|---|---|---|---:|"]
        for line in ("phases", "Parkin vs rest"):
            for st in STATS:
                o = obs[line][st]
                for leak, strength in CELLS:
                    v = sub[(sub.line == line) & (sub.leak == leak) & (sub.strength == strength)][st].to_numpy()
                    q = " / ".join(f"{x:.2f}" for x in np.percentile(v, [5, 25, 50, 75, 95]))
                    L.append(f"| {line} | {LABEL[st]} | {o:.2f} | {leak} | {strength} | {q} | {np.mean(v <= o + 1e-12):.2f} |")
        L.append("")
    return L


def draw_figure(df, obs, combo):
    """Empirical CDFs of the phase lines' shares, one curve per cell, observed score as a line."""
    import matplotlib.pyplot as plt
    fs = importlib.import_module("figstyle")
    n_ind, innov, mix = combo
    sub = df[(df.n_ind == n_ind) & (df.innovation == innov) & (df.mixing == mix)]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex=True, sharey=True)
    styles = {(1.0, 0.0): ("0.0", "-"), (0.1, 0.0): ("0.45", "-"), (0.03, 0.0): ("0.7", "-"),
              (1.0, 0.2): ("0.0", ":"), (0.03, 0.2): ("0.6", ":")}
    for ax, st, letter in zip(axes.ravel(), STATS, "ABCD"):
        for leak, strength in CELLS:
            v = np.sort(sub[(sub.line == "phases") & (sub.leak == leak) & (sub.strength == strength)][st].to_numpy())
            col, ls = styles[(leak, strength)]
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=col, ls=ls, lw=1.3,
                    label=f"copying factor {leak:g}, local innovation {strength:g}")
        ax.axvline(obs["phases"][st], color="0.0", lw=1.0, ls="--")
        ax.set_title(LABEL[st].replace("F_ST", "$F_{ST}$"), fontsize=8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        fs.panel_label(ax, letter)
    for ax in axes[1]:
        ax.set_xlabel("share of alternatives the phase lines beat")
    for ax in axes[:, 0]:
        ax.set_ylabel("cumulative share of simulated records")
    axes[1, 0].legend(fontsize=6, frameon=False, loc="upper left")   # C's upper left is empty; in A it hid the observed line
    fig.tight_layout()
    fs.save_all(fig, FIG, close=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--alt", type=int, default=50)
    ap.add_argument("--figure-only", action="store_true",
                    help="redraw the figure from RECORDS_CSV and the observed shares in OUT_MD, without simulating")
    args = ap.parse_args()
    if args.figure_only:
        import pandas as pd
        import re
        rec_df = pd.read_csv(RECORDS_CSV)
        row = re.search(r"^\| phases \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \|$", OUT_MD.read_text(), re.M)
        obs = {"phases": dict(zip(STATS, (float(x) for x in row.groups())))}
        first = rec_df.iloc[0]
        draw_figure(rec_df, obs, (int(first.n_ind), float(first.innovation), float(first.mixing)))
        print(f"redrew figures/{FIG} from {RECORDS_CSV.name}")
        return 0
    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    sd = importlib.import_module("23_phases_vs_spatial_drift")
    t74 = importlib.import_module("74_phase_partition_test")
    a84 = importlib.import_module("84_phases_as_groups")
    import pandas as pd
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record

    data = rev.load_sets()["basin"]
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    if names != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    labels_ph, _ = ph.assign_primary_phases(names, coords.to_numpy(float))
    plist = sorted(set(labels_ph)); phase = np.array([plist.index(l) for l in labels_ph])
    parkin = np.array([1 if l == "Parkin" else 0 for l in labels_ph])
    pts = t74.km_xy(coords.to_numpy(float))
    alt_ph, opt_ph = alternatives(t74, pts, phase, args.alt, 89000)
    alt_pk, opt_pk = alternatives(t74, pts, parkin, args.alt, 89200)
    m_obs = data["m"]; d = data["d"]; totals = m_obs.sum(1)
    rates = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates[(rates.region == "basin") & (rates.model == "pooled")].iloc[0]
    cell = dict(innovation=float(row["innovation"]), mixing=float(row["mixing"]), n_ind=int(row["n_ind"]))
    obs_ph = shares(m_obs.astype(float), phase, alt_ph, opt_ph, rev, sd, d)
    obs_pk = shares(m_obs.astype(float), parkin, alt_pk, opt_pk, rev, sd, d)

    L = ["# Could the phase-line comparison see a copying boundary if there were one?", "",
         f"Produced by `analyses/89_partition_power.py`: {args.reps} simulated records per cell, rebuilt with "
         f"analysis 84's seeds at its calibrated cell (N {cell['n_ind']}, innovation {cell['innovation']}, "
         f"mixing {cell['mixing']}, 24 km along the rivers); {args.alt} same-size alternative divisions of each kind. "
         "Each entry is the median over simulated records of the share of alternatives the tested line beats, "
         "with the share of records in which it beats at least 90 percent of them in parentheses.", "",
         "| copying factor at the phase lines | local innovation | line tested | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |",
         "|---|---|---|---|---|---|---|"]
    fmt = lambda v: f"{np.median(v):.2f} ({np.mean(np.array(v) >= 0.9) * 100:.0f}%)"
    combos = [(cell["n_ind"], cell["innovation"], cell["mixing"])] + [c for c in COMBOS]
    _cal = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibration.csv")
    _cal = _cal[(_cal["region"] == "basin") & (_cal["model"] == "pooled")]
    for n_, inn_, mix_ in COMBOS:
        hit = _cal[(_cal["n_ind"] == n_) & np.isclose(_cal["innovation"], inn_) & np.isclose(_cal["mixing"], mix_)]
        if hit.empty or not bool(hit["stage2_matched"].iloc[0]):
            raise SystemExit(f"combination {(n_, inn_, mix_)} no longer passes the fifty-run "
                             "diversity confirmation; choose another before reporting power at it")
    rows = []
    for ci, (n_ind, innov, mix) in enumerate(combos):
        for leak, strength in CELLS:
            li, si = a84.LEAKS.index(leak), a84.STRENGTHS.index(strength)
            w = copying_weights(d, 24.0, labels=phase, leak=leak)
            res = {"phases": [], "Parkin vs rest": []}
            for r in range(args.reps):
                seed = 84000 + 10000 * li + 1000 * si + r      # analysis 84's seeds
                rng = np.random.default_rng(seed + 5_000_000)
                tgt = a84.group_targets_by(data["pooled"], phase, strength, rng)
                rec = drift_record(w, k=m_obs.shape[1], n_ind=n_ind, seed=seed,
                                   innovation=innov, mixing=mix,
                                   burnin=1200, initial="uniform", target=tgt)
                m = sample_record(rec, data["ranks"], totals, rng).astype(float)
                for line, lab, alt, opt in (("phases", phase, alt_ph, opt_ph),
                                            ("Parkin vs rest", parkin, alt_pk, opt_pk)):
                    sh = shares(m, lab, alt, opt, rev, sd, d)
                    res[line].append(sh)
                    rows.append(dict(n_ind=n_ind, innovation=innov, mixing=mix, leak=leak,
                                     strength=strength, line=line, rep=r, **sh))
            if ci == 0:     # the headline table stays at the calibrated combination
                for line, rr in res.items():
                    L.append(f"| {leak} | {strength} | {line} | " + " | ".join(fmt([x[c] for x in rr]) for c in STATS) + " |")
                    print(L[-1], flush=True)
            else:
                print(f"  combo {n_ind}/{innov}/{mix} cell {leak}/{strength} done", flush=True)
    rec_df = pd.DataFrame(rows)
    RECORDS_CSV.parent.mkdir(parents=True, exist_ok=True)
    rec_df.to_csv(RECORDS_CSV, index=False)
    obs = {"phases": obs_ph, "Parkin vs rest": obs_pk}
    L += ["", "The observed record (plug-in, these alternatives):", "",
          "| line | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |",
          "|---|---|---|---|---|",
          "| phases | " + " | ".join(f"{obs_ph[c]:.2f}" for c in ("f_alt", "f_opt", "b_alt", "b_opt")) + " |",
          "| Parkin vs rest | " + " | ".join(f"{obs_pk[c]:.2f}" for c in ("f_alt", "f_opt", "b_alt", "b_opt")) + " |", ""]
    L += distribution_section(rec_df, obs, combos)
    L += ["## Reading", "",
          "Compare the rows with a copying boundary (factor below 1) against the row without. If the tested line "
          "beats most alternatives only when copying is restricted at it, the comparison responds, and the "
          "observed record's failure to beat them weighs against a restriction there by as much as the "
          "at-or-below fractions below differ between cells, which is moderate at best. If the rows look alike, "
          "the comparison cannot see a copying boundary and the argument from it should be dropped. "
          "Read the medians together with the distributions section: a median that moves shows the "
          "comparison responds, and only the at-or-below fractions say how strongly the observed "
          "score weighs against a boundary.", ""]
    draw_figure(rec_df, obs, combos[0])
    L += [f"Figure written to figures/{FIG}.png and its siblings; every simulated record's shares are in "
          f"`{RECORDS_CSV.relative_to(ROOT)}`.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
