"""86_partition_posterior.py - the phase-line comparisons as posterior probabilities.

Analyses 74 and 83 place a phase scheme's between-group cultural F_ST within
ensembles of alternative divisions of the same sites into groups of the same
sizes, and report a percentile. That is a randomization comparison on the
plug-in statistic, the form rule 18 asks to replace with a Bayesian one where
one exists (2026-09-23). This script states the same comparison as a
posterior probability.

Each assemblage's class proportions are drawn from their posterior,
Dirichlet(counts + 1/2), 2,000 times (Jeffreys prior: it pulls each profile
slightly toward even, which raises every division's F_ST by a similar small
amount and so leans neither way in a comparison between divisions). For each
draw, F_ST is computed under the phase scheme and under alternative divisions
drawn from the two ensembles of analysis 74: same group sizes around randomly
placed centers, and same group sizes made as spatially compact as possible.
The reported quantity is

    P(the phase lines separate the pottery better than an alternative | data)

averaged over the posterior of the counts and over the alternatives. It
carries the uncertainty in the observed counts, which the percentile did not.
Two schemes: the phases of Figure 1 on the 28-assemblage basin set, and
Mainfort's (2003) phases on his table (29 sites, at least 100 decorated
sherds).

Output: output/findings/partition_posterior.md
Usage: .venv/bin/python analyses/86_partition_posterior.py [--draws 2000] [--alt 300]
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
OUT_MD = ROOT / "output" / "findings" / "partition_posterior.md"
PER_DRAW = 20      # alternatives of each kind compared per posterior draw


def compare(m, pidx, pts, rev, t74, n_draws, n_alt, seed, dist=None):
    k = len(np.unique(pidx))
    sizes = np.bincount(pidx)
    slot = np.repeat(np.arange(k), sizes)
    rng = np.random.default_rng(seed)
    alt, opt = [], []
    for _ in range(n_alt):
        a = t74.assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)
        alt.append(a)
        opt.append(t74.local_search(pts, a, sizes))
    n_a = m.sum(1)
    rng = np.random.default_rng(seed + 1)
    f_ph, win_alt, win_opt, f_alt_med, f_opt_med = [], [], [], [], []
    for _ in range(n_draws):
        md = np.array([rng.dirichlet(r + 0.5) for r in m]) * n_a[:, None]
        fp = rev.fst_by(md, pidx)
        fa = np.array([rev.fst_by(md, alt[i]) for i in rng.choice(n_alt, PER_DRAW)])
        fo = np.array([rev.fst_by(md, opt[i]) for i in rng.choice(n_alt, PER_DRAW)])
        f_ph.append(fp); win_alt.append(np.mean(fp > fa)); win_opt.append(np.mean(fp > fo))
        f_alt_med.append(np.median(fa)); f_opt_med.append(np.median(fo))
    f_ph = np.array(f_ph)
    q = lambda a: (float(np.median(a)), float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)))
    out = dict(phase=q(f_ph), alt=q(f_alt_med), opt=q(f_opt_med),
               p_alt=float(np.mean(win_alt)), p_opt=float(np.mean(win_opt)),
               plug_in=float(rev.fst_by(m, pidx)))
    # Boundary excess at the lines (2026-09-23): is the step in similarity at
    # the phase lines, beyond distance, a property of these lines or of any
    # division of the map? Same posterior draws, fewer of them (the statistic
    # builds a similarity matrix per call).
    if dist is not None:
        sd = importlib.import_module("23_phases_vs_spatial_drift")
        rng = np.random.default_rng(seed + 2)
        b_ph, b_alt, b_opt, bw_alt, bw_opt = [], [], [], [], []
        for _ in range(min(n_draws, 500)):
            md = np.array([rng.dirichlet(r + 0.5) for r in m]) * n_a[:, None]
            bp = sd.boundary_excess_labeled(md, dist, pidx)
            ba = np.array([sd.boundary_excess_labeled(md, dist, alt[i]) for i in rng.choice(n_alt, 10)])
            bo = np.array([sd.boundary_excess_labeled(md, dist, opt[i]) for i in rng.choice(n_alt, 10)])
            b_ph.append(bp); b_alt.append(np.median(ba)); b_opt.append(np.median(bo))
            bw_alt.append(np.mean(bp > ba)); bw_opt.append(np.mean(bp > bo))
        out.update(be_phase=q(b_ph), be_alt=q(b_alt), be_opt=q(b_opt),
                   be_p_alt=float(np.mean(bw_alt)), be_p_opt=float(np.mean(bw_opt)),
                   be_plug=float(sd.boundary_excess_labeled(m, dist, pidx)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--alt", type=int, default=300)
    args = ap.parse_args()

    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    t74 = importlib.import_module("74_phase_partition_test")
    a83 = importlib.import_module("83_mainfort_replication")
    from mls_emergence.dataio.matrix import read_mainfort_replication

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    xy = coords.to_numpy(float)
    labs, _ = ph.assign_phases_by_territory(names, xy)
    plist = sorted(set(labs)); pidx = np.array([plist.index(l) for l in labs])
    data = rev.load_sets()["basin"]
    if [str(n) for n in data["names"]] != names:
        raise RuntimeError("assemblage order differs between loaders")
    res_basin = compare(counts.to_numpy(float), pidx, t74.km_xy(xy), rev, t74, args.draws, args.alt, 86000,
                        dist=data["d"])
    # The Parkin phase against the rest of the basin (2026-09-23): the paper's
    # second question is whether the phase treated as a polity stands apart.
    parkin_idx = np.array([1 if l == "Parkin" else 0 for l in labs])
    res_parkin = compare(counts.to_numpy(float), parkin_idx, t74.km_xy(xy), rev, t74, args.draws, args.alt, 86200,
                         dist=data["d"])

    mc, mco, mph = read_mainfort_replication(100)
    sizes = mph.value_counts()
    keep = ~mph.isin(sizes[sizes < a83.MIN_PHASE_MEMBERS].index)
    mc, mco, mph = mc[keep], mco[keep], mph[keep]
    ml = sorted(set(mph)); midx = np.array([ml.index(p) for p in mph])
    res_mf = compare(mc.to_numpy(float), midx, t74.km_xy(mco.to_numpy(float)), rev, t74,
                     args.draws, args.alt, 86500)

    def block(title, r, n_sites, groups):
        return [f"## {title}", "",
                f"{n_sites} sites in {groups} groups. Plug-in F_ST {r['plug_in']:.4f}.", "",
                "| quantity | posterior median | 95% credible interval |", "|---|---|---|",
                f"| F_ST under the phase scheme | {r['phase'][0]:.4f} | {r['phase'][1]:.4f} to {r['phase'][2]:.4f} |",
                f"| median F_ST, same-size divisions around random centers | {r['alt'][0]:.4f} | {r['alt'][1]:.4f} to {r['alt'][2]:.4f} |",
                f"| median F_ST, same-size divisions made compact | {r['opt'][0]:.4f} | {r['opt'][1]:.4f} to {r['opt'][2]:.4f} |",
                "",
                f"- P(phase lines separate the pottery better than a division around random centers | data) = **{r['p_alt']:.2f}**",
                f"- P(phase lines separate the pottery better than a compact division | data) = **{r['p_opt']:.2f}**", ""] + (
                ["Boundary excess at the lines (similarity lost across a line beyond what distance "
                 "predicts, river distance), 500 posterior draws, 10 alternatives of each kind per draw:", "",
                 "| quantity | posterior median | 95% credible interval |", "|---|---|---|",
                 f"| boundary excess at the phase lines (plug-in {r['be_plug']:+.1f}) | {r['be_phase'][0]:+.1f} | {r['be_phase'][1]:+.1f} to {r['be_phase'][2]:+.1f} |",
                 f"| median, same-size divisions around random centers | {r['be_alt'][0]:+.1f} | {r['be_alt'][1]:+.1f} to {r['be_alt'][2]:+.1f} |",
                 f"| median, same-size divisions made compact | {r['be_opt'][0]:+.1f} | {r['be_opt'][1]:+.1f} to {r['be_opt'][2]:+.1f} |", "",
                 f"- P(larger boundary excess at the phase lines than at a random-center division | data) = **{r['be_p_alt']:.2f}**",
                 f"- P(larger boundary excess at the phase lines than at a compact division | data) = **{r['be_p_opt']:.2f}**", ""]
                if "be_phase" in r else [])

    L = ["# The phase-line comparisons as posterior probabilities", "",
         f"Produced by `analyses/86_partition_posterior.py`: {args.draws} posterior draws of every "
         f"assemblage's class proportions (Dirichlet(counts + 1/2)), {args.alt} alternative divisions of "
         f"each kind, {PER_DRAW} of each compared per draw. Replaces the ensemble percentiles of "
         f"analyses 74 and 83 as the reported quantity (rule 18).", ""]
    L += block("The phases of Figure 1, basin set", res_basin, len(names), len(plist))
    L += block("The Parkin phase against the rest of the basin", res_parkin, len(names), 2)
    k2 = importlib.import_module("74_phase_partition_test").ari(parkin_idx, np.asarray(data["labels"]))
    L += [f"Agreement (adjusted Rand index) between the Parkin-versus-rest division and the two spatial "
          f"clusters the site layout supports: {k2:.3f}.", ""]
    L += block("Mainfort's (2003) phases, his table", res_mf, len(mph), len(ml))
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
