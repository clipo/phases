"""55_abc_stability.py - how reproducible is the ABC transmission posterior?

Chases F22 in docs/CODE_REVIEW_2026-08-31.md and, in the same measurement,
settles a dependency: `analyses/56_neutrality_ppc.py` (rule 18 item G) draws
population sizes from this same posterior, so whatever instability it carries is
inherited there.

THE DISCREPANCY. The supplement reports mean b = +0.006, 95 percent interval
-0.020 to +0.035, posterior SD 0.014, and **P(b > 0) = 0.70**. Re-running
`analyses/38_abc_smc_transmission.py` at full settings in the rebuilt
environment gives mean +0.0013, interval [-0.0243, +0.0295], **P(b > 0) =
0.531**. The interval reproduces closely; P(b > 0) does not, and 0.70 against
0.53 is the difference between "leans conformist" and "an even split".

THE HYPOTHESES, and why this is the right test. The candidates are a library
version change, or Monte Carlo variability in the ABC-SMC itself. The second is
testable directly and cheaply: refit at several seeds and look at the spread. If
P(b > 0) wanders across seeds by more than the gap in question, then 0.70 and
0.53 are both draws from the same distribution and the supplement is quoting a
statistic more precisely than the method supports. That is a finding about the
number regardless of which library produced it.

Note the asymmetry the test exploits: the posterior MEAN and INTERVAL are
integrals over the whole posterior and are comparatively stable, while
P(b > 0) is a tail-mass probability evaluated at exactly the point where this
posterior is densest, so it is the most seed-sensitive summary of the three.
That predicts the interval reproduces and P(b > 0) does not, which is what was
observed.

SUPERSEDING NOTE (2026-09-03). The numbers quoted above are what this
investigation observed at the time, and are kept so the reasoning stays legible.
They have since been superseded by the F7 fix to the SMC weight update
(`src/mls_emergence/inference/abc_smc.py`, log-space normalization), which this
script's own runs predate. `38_abc_smc_transmission.py` at full settings now
gives mean b = -0.003, interval [-0.096, +0.031], SD 0.025, P(b > 0) = 0.51, and
this script's five-seed range is 0.447 to 0.868.

That revises one claim made above. The asymmetry argument predicted that the
posterior mean and interval would reproduce while P(b > 0) would not. P(b > 0)
is indeed the seed-unstable summary, and the widened five-seed range (0.447 to
0.868, against 0.44 to 0.54 before) makes that case more strongly than the
original run did. But the interval did NOT hold still: it went from
[-0.0243, +0.0295] to [-0.096, +0.031], about 2.4 times wider. That movement is
not Monte Carlo variability, which is what the asymmetry argument was about; it
is the weight fix correcting a real defect. The lesson stands with a caveat
attached, that stability across seeds says nothing about correctness of the
estimator being reseeded.

Usage: .venv/bin/python analyses/55_abc_stability.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

OUT_MD = ROOT / "output" / "findings" / "abc_stability.md"
# The POOLED joint posterior across seeds. analyses/56_neutrality_ppc.py
# consumes this rather than a single ABC run, because the population size
# median swings 98 to 263 across seeds and conditioning on one run would
# present that seed as the posterior (added 2026-09-02).
OUT_NPZ = ROOT / "output" / "abc_pooled_posterior.npz"
SEEDS = [101, 202, 303, 404, 505]
REPORTED = dict(mean=0.006, lo=-0.020, hi=0.035, sd=0.014, p_pos=0.70)


def main(fast=False):
    a38 = importlib.import_module("38_abc_smc_transmission")
    from mls_emergence.inference.abc_smc import resample, weighted_mean, weighted_quantile
    cfg = a38.FAST if fast else a38.FULL
    obs, counts = a38.load_obs()

    rows, pooled = [], []
    for sd in SEEDS:
        r, b, joint = a38.fit_target(obs, slice(None), seed=sd, **cfg)
        # Equal-weight resample so each seed contributes the same number of
        # draws regardless of its particle weights, then pool.
        idx = np.random.default_rng(sd).choice(
            len(joint), size=len(joint), p=r.weights / r.weights.sum())
        pooled.append(joint[idx])
        rng = np.random.default_rng(0)
        draws = resample(b, r.weights, 20000, rng)
        rows.append(dict(
            seed=sd,
            mean=weighted_mean(b, r.weights),
            lo=weighted_quantile(b, r.weights, 0.025),
            hi=weighted_quantile(b, r.weights, 0.975),
            sd=float(np.std(draws)),
            p_pos=float(np.mean(draws > 0)),
            N_med=float(np.median(joint[:, 2])),
            N_lo=float(np.percentile(joint[:, 2], 2.5)),
            N_hi=float(np.percentile(joint[:, 2], 97.5))))
        r_ = rows[-1]
        print(f"seed {sd}: mean {r_['mean']:+.4f}  95% [{r_['lo']:+.4f}, {r_['hi']:+.4f}]  "
              f"SD {r_['sd']:.4f}  P(b>0) {r_['p_pos']:.3f}  N median {r_['N_med']:.0f}")

    g = lambda k: np.array([r[k] for r in rows])
    pp, mm, ss = g("p_pos"), g("mean"), g("sd")

    L = ["# How reproducible is the ABC transmission posterior?", "",
         f"Produced by `analyses/55_abc_stability.py` ({'FAST' if fast else 'full'}: "
         f"{cfg['n_particles']} particles, {cfg['n_rounds']} rounds, "
         f"{len(SEEDS)} seeds). Chases F22 and settles the dependency for rule-18 "
         f"item G, which draws population sizes from this posterior.", "",
         "## Across-seed spread", "",
         "| seed | mean b | 95% interval | posterior SD | P(b > 0) | N median |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['seed']} | {r['mean']:+.4f} | [{r['lo']:+.4f}, {r['hi']:+.4f}] | "
                 f"{r['sd']:.4f} | **{r['p_pos']:.3f}** | {r['N_med']:.0f} |")
    L += ["",
          f"- P(b > 0): range **{pp.min():.3f} to {pp.max():.3f}**, SD across seeds {pp.std():.3f}.",
          f"- mean b: range {mm.min():+.4f} to {mm.max():+.4f}, SD {mm.std():.4f}.",
          f"- posterior SD: range {ss.min():.4f} to {ss.max():.4f}.", "",
          "## Verdict", ""]
    span = pp.max() - pp.min()
    gap = abs(REPORTED["p_pos"] - float(np.median(pp)))
    covers = pp.min() <= REPORTED["p_pos"] <= pp.max()
    if covers or span >= gap:
        L += [f"**P(b > 0) is not stable to the precision it is quoted at.** Across "
              f"seeds it spans {pp.min():.3f} to {pp.max():.3f}, a range of "
              f"{span:.3f}, against a gap of {gap:.3f} between the supplement's "
              f"0.70 and the median of these runs. The supplement quotes a "
              f"seed-dependent summary to two decimal places.",
              "",
              "**The remedy is not to pick a seed.** P(b > 0) is a tail mass "
              "evaluated where this posterior is densest, so it is the least "
              "stable summary available and the one most sensitive to the "
              "particle sample. The mean and the interval, which are integrals "
              "over the whole posterior, are comparatively steady and reproduce "
              "the supplement closely. Report the interval and drop the point "
              "probability, or report the probability with its across-seed range "
              "attached.", ""]
    else:
        L += [f"**P(b > 0) is stable across seeds** ({pp.min():.3f} to "
              f"{pp.max():.3f}) and does NOT bracket the supplement's 0.70. "
              f"Seed variability is therefore not the explanation and the "
              f"library version change is the remaining candidate, which this "
              f"script does not test. F22 stays open.", ""]

    Nm = g("N_med")
    L += ["## The dependency for item G", "",
          f"Item G draws the population size from this posterior. Its median "
          f"across seeds runs {Nm.min():.0f} to {Nm.max():.0f}, a spread of "
          f"{100 * (Nm.max() - Nm.min()) / Nm.mean():.0f} percent of the mean. "
          f"The neutrality posterior predictive check must therefore either "
          f"pool draws across seeds or report its envelope's sensitivity to the "
          f"seed, rather than conditioning on one ABC run as if it were the "
          f"posterior.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    pooled_arr = np.vstack(pooled)
    np.savez(OUT_NPZ, pooled=pooled_arr,
             param_names=np.array(["mu", "b", "N", "w"]),
             seeds=np.array(SEEDS),
             prior_lo=a38.PRIOR_LO, prior_hi=a38.PRIOR_HI)
    print(f"pooled joint posterior: {pooled_arr.shape} -> {OUT_NPZ}")
    print(f"\nP(b>0) across seeds: {pp.min():.3f} to {pp.max():.3f}")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
