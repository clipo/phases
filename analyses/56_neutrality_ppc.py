"""56_neutrality_ppc.py - the neutrality envelope as a posterior predictive check.

Rule 18, item G of docs/FREQUENTIST_INVENTORY.md, and the last of the eight.

WHAT IT REPLACES. `analyses/18_kandler_shennan_neutrality.py` simulates neutral
Wright-Fisher drift forward over the sequence at a RANGE OF FIXED effective
population sizes, takes the 2.5th and 97.5th percentiles of the simulated
diversity trajectories as a "neutral envelope", and checks whether the observed
trajectory falls inside it. The supplement reads the result as "the test cannot
reject neutral transmission".

WHY THAT IS NOT ALREADY A POSTERIOR PREDICTIVE CHECK. A PPC draws parameters
from a fitted posterior. This sweeps fixed values of Ne, so the band is a
sensitivity envelope, and it is then read as an acceptance region, which is the
frequentist framing rule 18 removes. The rule explicitly retains posterior
predictive checks, so the conversion is to make it into one properly.

THE POSTERIOR IS ALREADY THERE. The ABC-SMC transmission fit estimates a joint
posterior over the innovation rate, the transmission bias, **the effective
population size** and the time-averaging window. Drawing Ne from it instead of
sweeping arbitrary values turns the envelope into a genuine PPC, using a
posterior this project already has.

POOLED ACROSS SEEDS, AND THIS IS NOT OPTIONAL. `analyses/55_abc_stability.py`
measured the ABC posterior's run-to-run variability: the median of N swings from
98 to 263 across five seeds, about 90 percent of its mean. Conditioning on one
ABC run would present that run's seed as the posterior. This script therefore
consumes the POOLED joint posterior that 55 writes, and reports the envelope's
sensitivity to the seed alongside the pooled result.

WHAT THE CHECK CAN AND CANNOT SAY. A PPC that passes says the neutral model
reproduces the observed diversity trajectory. It does not say transmission was
neutral: other processes can produce the same trajectory, and Kandler & Shennan's
own caveat applies, that with an ordinal axis and an uncertain effective size the
predictive interval is wide. The reportable statement is a Bayesian predictive
p-value with its interval, not a rejection or a failure to reject.

Usage: .venv/bin/python analyses/56_neutrality_ppc.py [--fast]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
from figstyle import OI_BLUE, OI_VERMIL, OI_GREEN, save  # noqa: E402

OUT_MD = ROOT / "output" / "findings" / "neutrality_ppc.md"
POOLED = ROOT / "output" / "abc_pooled_posterior.npz"
N_BINS = 6
GEN_PER_BIN = 20


def gini(counts):
    c = np.asarray(counts, float)
    tot = c.sum()
    if tot <= 0:
        return np.nan
    p = c / tot
    return float(1.0 - np.sum(p ** 2))


def simulate_trajectory(x0, Ne, mu, rng, n_bins=N_BINS, gen=GEN_PER_BIN):
    """Wright-Fisher drift with innovation, started at the earliest-bin
    frequencies. Same construction as analysis 18 except that innovation is
    included, because the posterior supplies a rate for it."""
    K = len(x0)
    Ne = max(int(round(Ne)), 2)
    pop = rng.choice(K, size=Ne, p=x0)
    out = [gini(np.bincount(pop, minlength=K).astype(float))]
    for _ in range(n_bins - 1):
        for _g in range(gen):
            freq = np.bincount(pop, minlength=K).astype(float)
            freq = freq / freq.sum()
            pop = rng.choice(K, size=Ne, p=freq)
            if mu > 0:                       # innovation: reassign at rate mu
                hit = rng.random(Ne) < mu
                if hit.any():
                    pop = pop.copy()
                    pop[hit] = rng.integers(0, K, size=int(hit.sum()))
        out.append(gini(np.bincount(pop, minlength=K).astype(float)))
    return np.array(out)


def main(fast=False):
    n_draw = 200 if fast else 1500
    if not POOLED.exists():
        raise FileNotFoundError(
            f"{POOLED} is missing. It is the pooled ABC posterior this check "
            "draws from; generate it with\n"
            "  .venv/bin/python analyses/55_abc_stability.py")
    d = np.load(POOLED, allow_pickle=True)
    pooled, seeds = d["pooled"], d["seeds"]
    per_seed = np.array_split(pooled, len(seeds))
    print(f"pooled posterior: {pooled.shape[0]} draws over {len(seeds)} seeds")

    counts, coords = mf._load_curated()
    M = counts.to_numpy(float)
    # ORIENTATION (fixed 2026-09-03). Must use 17_basin_results.oriented_ca,
    # which flips the correspondence axis against the radiocarbon anchors so
    # that increasing = later. The raw mf.correspondence_axis sign is arbitrary
    # and on this basin runs EXACTLY BACKWARDS: Spearman(raw, oriented) = -1.000.
    # Using it reversed the sequence, which flipped the sign of any directional
    # quantity and made forward simulation start from the latest bin.
    import importlib
    res = importlib.import_module("17_basin_results")
    ca, _ = res.oriented_ca(counts)
    bins = pd.qcut(ca.rank(), N_BINS, labels=False,
                   duplicates="drop").to_numpy()
    binmat = np.array([M[bins == b].sum(0) for b in sorted(set(bins))])
    obs_div = np.array([gini(r) for r in binmat])
    x0 = binmat[0] / binmat[0].sum()
    print("observed diversity:", np.round(obs_div, 4).tolist())

    def envelope(draws, tag, seed):
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(draws), size=min(n_draw, len(draws)), replace=True)
        sims = np.array([simulate_trajectory(x0, draws[i, 2], draws[i, 0], rng)
                         for i in idx])
        lo, hi = np.percentile(sims, [2.5, 97.5], axis=0)
        within = int(np.sum((obs_div >= lo) & (obs_div <= hi)))
        # Bayesian predictive p-value on the whole trajectory, via the summed
        # squared discrepancy from the predictive mean (a global check rather
        # than six per-bin ones, so no multiplicity arises).
        mean = sims.mean(0)
        T_obs = float(np.sum((obs_div - mean) ** 2))
        T_rep = np.sum((sims - mean) ** 2, axis=1)
        pb = float(np.mean(T_rep >= T_obs))
        print(f"  {tag}: {within}/{N_BINS} bins inside, predictive p = {pb:.3f}")
        return dict(lo=lo, hi=hi, within=within, pb=pb, sims=sims)

    print("posterior predictive envelopes:")
    pooled_env = envelope(pooled, "pooled", 0)
    seed_envs = [envelope(p, f"seed {s}", 100 + i)
                 for i, (s, p) in enumerate(zip(seeds, per_seed))]

    pbs = np.array([e["pb"] for e in seed_envs])
    wins = np.array([e["within"] for e in seed_envs])
    L = ["# The neutrality envelope as a posterior predictive check", "",
         f"Produced by `analyses/56_neutrality_ppc.py` ({'FAST' if fast else 'full'}, "
         f"{n_draw} predictive draws). Replaces the fixed-Ne sweep in "
         f"`analyses/18_kandler_shennan_neutrality.py` "
         f"(`docs/FREQUENTIST_INVENTORY.md` item G).", "",
         f"Effective population size and innovation rate are drawn from the "
         f"POOLED ABC-SMC posterior ({pooled.shape[0]} draws over "
         f"{len(seeds)} seeds), not swept over fixed values and not taken from a "
         f"single run.", "",
         "## Result", "",
         "Observed per-bin Gini-Simpson diversity: "
         + ", ".join(f"{v:.3f}" for v in obs_div) + ".", "",
         "| posterior | bins inside the 95% predictive interval | Bayesian predictive p |",
         "|---|---|---|",
         f"| **pooled** | **{pooled_env['within']}/{N_BINS}** | **{pooled_env['pb']:.3f}** |"]
    for s, e in zip(seeds, seed_envs):
        L.append(f"| seed {s} alone | {e['within']}/{N_BINS} | {e['pb']:.3f} |")
    L += ["",
          f"Across single-seed posteriors the predictive p ranges "
          f"{pbs.min():.3f} to {pbs.max():.3f} and the bins-inside count ranges "
          f"{wins.min()} to {wins.max()} of {N_BINS}. That spread is the "
          f"dependency `analyses/55_abc_stability.py` measured (N median 98 to "
          f"263 across seeds) propagating into the check, and it is the reason "
          f"the pooled row is the one to report.", "",
          "## Reading", ""]
    if 0.05 < pooled_env["pb"] < 0.95:
        L += ["The neutral model reproduces the observed diversity trajectory: "
              "the predictive p sits away from both tails and the observed "
              "trajectory lies inside the predictive interval in "
              f"{pooled_env['within']} of {N_BINS} bins.", "",
              "**What this does not say.** A passing predictive check says the "
              "neutral model is adequate for this trajectory, not that "
              "transmission was neutral; other processes produce the same "
              "trajectory. Kandler and Shennan's own caveat also applies, that "
              "with an ordinal axis and an uncertain effective size the "
              "predictive interval is wide, and drawing Ne from a posterior "
              "whose median moves by a factor of nearly three across seeds "
              "widens it further rather than narrowing it. This is concordant "
              "evidence, and weak, which is what the supplement already says."]
    else:
        L += [f"**The observed trajectory sits in the tail of the posterior "
              f"predictive distribution (p = {pooled_env['pb']:.3f}).** The "
              f"neutral model does not reproduce it, which is a substantive "
              f"departure from what the supplement reports and must be chased "
              f"before anything is written from it."]
    L += ["", "## What replaces what", "",
          "| before | after |", "|---|---|",
          "| sweep over fixed Ne | Ne drawn from the ABC posterior |",
          "| 95% band read as an acceptance region | 95% posterior predictive interval |",
          "| \"the test cannot reject neutral transmission\" | a Bayesian predictive p with its across-seed range |", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    t = np.arange(N_BINS)
    ax.fill_between(t, pooled_env["lo"], pooled_env["hi"], color=OI_BLUE,
                    alpha=0.25, label="95% posterior predictive")
    ax.plot(t, pooled_env["sims"].mean(0), color=OI_BLUE, lw=1.2,
            label="predictive mean")
    ax.plot(t, obs_div, "o-", color=OI_VERMIL, ms=4, label="observed")
    ax.set_xlabel("CA seriation bin (early to late)")
    ax.set_ylabel("Gini-Simpson diversity")
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    save(fig, "fig_neutrality_ppc")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
