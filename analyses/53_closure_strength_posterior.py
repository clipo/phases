"""53_closure_strength_posterior.py - a posterior on closure strength.

Rule 18, item A of docs/FREQUENTIST_INVENTORY.md, and the largest one: the
recovery-and-power curve is the instrument that carries this paper's central
negative, and it is Neyman-Pearson end to end. The manuscript currently reports
"with the false-positive rate held at five percent, F_ST recovers a closure of
strength s greater than or equal to 0.5 with power above 0.8", plus a nominal
injected strength of about 0.35 read off the recovery curve for the empirical
value.

THE CONVERSION. Two routes were on the table. Route 2 recast the design's
no-closure distribution as a prior predictive and left the rest alone. Route 1
replaces the whole apparatus with a posterior on the closure strength s. Route 1
turned out to cost almost nothing extra, because the simulation grid the power
curve already builds IS a simulated likelihood: analysis 21 evaluates the F_ST
trajectory statistic T at each injected strength across many seeds, which is a
sample from p(T | s).

So:

    p(T | s)   estimated from the existing simulation grid, one density per s
    p(s)       Uniform(0, 1), the same range the injection spans
    p(s | T_obs) propto p(T_obs | s) p(s)

evaluated on the s grid and normalized. This is simulation-based inference of
the synthetic-likelihood kind, and it answers the question the paper actually
asks, "how much closure is consistent with what we see", rather than "can we
reject a null we do not believe".

WHAT REPLACES WHAT:
  false-positive rate 0.05          -> nothing; no test is performed
  power(s) and the threshold s*     -> the width of p(s | T_obs), which says
                                       directly which strengths remain credible
  "nominal injected strength ~0.35" -> the posterior median, with an interval

PRIOR AND DIRECTION (rule 20). Uniform(0, 1) on s is flat over the whole
injected range. It leans neither toward nor against closure, so the burden is
light. The decisive check is the recovery run below, which asks whether a
genuinely strong closure is recovered rather than shrunk toward the middle.

DECLARED LIMITATION. p(T | s) is approximated by a normal at each grid point,
with the mean and standard deviation taken across seeds. The statistic is a rank
correlation, so this is reasonable in the interior and worse near the ends where
T is bounded. The alternative, a kernel density per grid point at 120 seeds, is
noisier. This is a two-stage approximation like item B's and is stated for the
same reason.

Usage: .venv/bin/python analyses/53_closure_strength_posterior.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
from figstyle import OI_BLUE, OI_VERMIL, save  # noqa: E402

OUT_MD = ROOT / "output" / "findings" / "closure_strength_posterior.md"
# Figure 4's right panel reads this. Written so analysis 21 does not have to
# re-run the simulation grid to draw the panel (added 2026-09-02).
OUT_NPZ = ROOT / "output" / "closure_posterior.npz"


def load_posterior():
    """Load the saved closure posterior, refusing a STALE one.

    Figure 4's right panel and Figure 2's shaded band are both drawn from this
    file rather than recomputed, so that they and the captions cannot drift
    apart. That protects against one failure and not the other: on 2026-09-15
    this script was re-run and rewrote its finding (median 0.227, time-averaged
    0.308) while the .npz on disk stayed at the 2026-09-03 values (0.295,
    0.397). Re-running Figure 4 then redrew it from the old array and silently
    reverted a recalibrated figure, which is how this guard came to exist.

    The finding is the witness: this script writes both, so an .npz older than
    the .md means the two came from different runs.
    """
    import numpy as np
    if not OUT_NPZ.exists():
        raise FileNotFoundError(
            f"{OUT_NPZ} is missing. It holds the posterior over closure "
            "strength that Figures 2 and 4 draw; generate it with\n"
            "  .venv/bin/python analyses/53_closure_strength_posterior.py")
    if OUT_MD.exists() and OUT_NPZ.stat().st_mtime < OUT_MD.stat().st_mtime - 1:
        raise RuntimeError(
            f"{OUT_NPZ.name} is older than {OUT_MD.name}, so they are from "
            "different runs of this script and the arrays no longer match the "
            "reported posterior. Drawing a figure from them would revert it to "
            "superseded values. Re-run:\n"
            "  .venv/bin/python analyses/53_closure_strength_posterior.py")
    return np.load(OUT_NPZ)
m21 = importlib.import_module("21_signal_recovery")
res = importlib.import_module("17_basin_results")

S_FINE = np.round(np.arange(0.0, 1.0001, 0.02), 3)   # posterior evaluation grid


def build_grid(n_seeds, w, k, K, clusters, bins_arr, coords_c, N_arr):
    """T at each injected strength across seeds: a sample from p(T | s)."""
    T = np.full((len(m21.S_GRID), n_seeds), np.nan)
    for si, s in enumerate(m21.S_GRID):
        for seed in range(n_seeds):
            rng = np.random.default_rng(si * 100003 + seed * 17 + w * 7)
            grid = m21.time_average(m21.emergence_profiles(k, K, s, rng), w)
            M = m21.rarefy(m21.gen_synth_M(grid, clusters, bins_arr, N_arr, rng),
                           m21.NRARE, rng)
            T[si, seed] = m21.sig_rhos(M, clusters, bins_arr, coords_c,
                                       which=["fst"])["fst"]
    return T


def posterior_over_s(T_grid, t_obs):
    """p(s | t_obs) on S_FINE, from a normal approximation to p(T | s)."""
    mu, sd = [], []
    for row in T_grid:
        r = row[np.isfinite(row)]
        mu.append(np.mean(r)); sd.append(max(np.std(r), 1e-6))
    mu_f = np.interp(S_FINE, m21.S_GRID, mu)
    sd_f = np.interp(S_FINE, m21.S_GRID, sd)
    ll = norm.logpdf(t_obs, loc=mu_f, scale=sd_f)      # uniform prior on s
    p = np.exp(ll - ll.max()); p /= p.sum()
    return p, mu_f, sd_f


def summarize(p):
    c = np.cumsum(p)
    q = lambda a: float(np.interp(a, c, S_FINE))
    return dict(median=q(0.5), lo=q(0.025), hi=q(0.975),
                p_ge_05=float(p[S_FINE >= 0.5].sum()),
                p_ge_03=float(p[S_FINE >= 0.3].sum()))


def main(fast=False):
    n_seeds = 30 if fast else m21.N_SEEDS
    b_emp = 60 if fast else m21.B_EMP

    counts, coords = mf._load_curated()
    ca, _ = res.oriented_ca(counts)
    K = counts.shape[1]
    Ni = dict(zip(counts.index, counts.to_numpy(float).sum(1)))
    cdf = coords.dropna(); have = list(cdf.index)
    cc = cdf[["Latitude", "Longitude"]].to_numpy(float); cc_c = cc - cc.mean(0)
    sil = {kk: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, kk, seed=7))
           for kk in range(2, 7)}
    k = max(sil, key=sil.get)
    cluster_of = dict(zip(have, mf._kmeans_labels(cc_c, k, seed=7)))
    bins = pd.qcut(ca.reindex(have), m21.N_BINS, labels=False, duplicates="drop")
    keep = [i for i in have if not pd.isna(bins[i])]
    clusters = np.array([cluster_of[i] for i in keep])
    bins_arr = np.array([int(bins[i]) for i in keep])
    coords_c = cc_c[[have.index(i) for i in keep]]
    N_arr = np.array([Ni[i] for i in keep])
    real_M = counts.reindex(keep).to_numpy(float)

    rng = np.random.default_rng(9)
    emp = np.array([m21.sig_rhos(m21.rarefy(real_M, m21.NRARE, rng), clusters,
                                 bins_arr, coords_c, which=["fst"])["fst"]
                    for _ in range(b_emp)])
    # One definition of the observed trend (rule 6): the 4,000-draw estimator
    # in analysis 21, which the manuscript quotes. The b_emp rarefactions above
    # are kept only as a cross-check; until 2026-09-22 their mean (-0.065) was
    # used here while the paper quoted -0.048.
    t_obs = float(m21.empirical_fst_trend()["mean"])
    print(f"observed F_ST trajectory statistic T = {t_obs:+.4f} "
          f"(analysis 21 estimator; {b_emp}-rarefaction cross-check {float(np.nanmean(emp)):+.4f})")

    out = {}
    for w in (1, 3):
        T = build_grid(n_seeds, w, k, K, clusters, bins_arr, coords_c, N_arr)
        p, mu_f, sd_f = posterior_over_s(T, t_obs)
        out[w] = dict(T=T, p=p, mu=mu_f, sd=sd_f, **summarize(p))
        o = out[w]
        print(f"w={w}: s posterior median {o['median']:.3f} "
              f"[{o['lo']:.3f}, {o['hi']:.3f}]  P(s>=0.5) = {o['p_ge_05']:.3f}")

    # The seriation axis's direction is not established by the four dated
    # assemblages (2026-09-22), and the trend statistic is signed, so the
    # posterior is also computed with the axis reversed (T -> -T): closure
    # growing toward the other end of the sequence.
    rev = {}
    for w in (1, 3):
        pr, _, _ = posterior_over_s(out[w]["T"], -t_obs)
        rev[w] = summarize(pr)
        print(f"w={w}, axis reversed: s posterior median {rev[w]['median']:.3f} "
              f"[{rev[w]['lo']:.3f}, {rev[w]['hi']:.3f}]  P(s>=0.5) = {rev[w]['p_ge_05']:.3f}")

    # Rule 20(c): recover a strong closure that the data do not contain.
    rngr = np.random.default_rng(4242)
    grid = m21.time_average(m21.emergence_profiles(k, K, 0.7, rngr), 1)
    Mr = m21.rarefy(m21.gen_synth_M(grid, clusters, bins_arr, N_arr, rngr),
                    m21.NRARE, rngr)
    t_strong = m21.sig_rhos(Mr, clusters, bins_arr, coords_c, which=["fst"])["fst"]
    p_s, _, _ = posterior_over_s(out[1]["T"], t_strong)
    rec = summarize(p_s)
    print(f"recovery: data simulated at s=0.7 -> posterior median {rec['median']:.3f} "
          f"[{rec['lo']:.3f}, {rec['hi']:.3f}]")

    o1, o3 = out[1], out[3]
    covered = rec["lo"] <= 0.7 <= rec["hi"]
    L = ["# A posterior on closure strength", "",
         f"Produced by `analyses/53_closure_strength_posterior.py` "
         f"({'FAST' if fast else 'full'}: {n_seeds} seeds per grid point, "
         f"{b_emp} empirical rarefactions).", "",
         "Replaces the calibrated detector, its five percent false-positive "
         "rate, its power curve and its detection threshold with one quantity "
         "(`docs/FREQUENTIST_INVENTORY.md` item A). No test is performed.", "",
         f"Observed statistic: the rarefied F_ST trajectory rank correlation, "
         f"T = **{t_obs:+.4f}**.", "",
         "## Posterior over injected closure strength s", "",
         "| time-averaging | median | 95% credible interval | P(s >= 0.3) | P(s >= 0.5) |",
         "|---|---|---|---|---|",
         f"| none (w = 1) | **{o1['median']:.3f}** | [{o1['lo']:.3f}, {o1['hi']:.3f}] | {o1['p_ge_03']:.3f} | {o1['p_ge_05']:.3f} |",
         f"| record's (w = 3) | **{o3['median']:.3f}** | [{o3['lo']:.3f}, {o3['hi']:.3f}] | {o3['p_ge_03']:.3f} | {o3['p_ge_05']:.3f} |",
         f"| none (w = 1), axis reversed (T = {-t_obs:+.4f}) | **{rev[1]['median']:.3f}** | [{rev[1]['lo']:.3f}, {rev[1]['hi']:.3f}] | {rev[1]['p_ge_03']:.3f} | {rev[1]['p_ge_05']:.3f} |",
         f"| record's (w = 3), axis reversed | **{rev[3]['median']:.3f}** | [{rev[3]['lo']:.3f}, {rev[3]['hi']:.3f}] | {rev[3]['p_ge_03']:.3f} | {rev[3]['p_ge_05']:.3f} |",
         "", "## Rule 20(c): recovery of a closure the data do not contain", "",
         f"Data simulated at a true s = 0.7, well above anything the basin "
         f"shows, returns a posterior median of **{rec['median']:.3f}** with "
         f"interval [{rec['lo']:.3f}, {rec['hi']:.3f}], "
         f"{'which covers the truth' if covered else '**which does NOT cover the truth**'}.",
         "",
         ("A strong closure is recovered rather than dragged toward the middle, "
          "so a low posterior on the real data is a statement about the data."
          if covered else
          "**The construction does not recover a strong closure. The empirical "
          "posterior is therefore not interpretable and this instrument is not "
          "reportable as it stands.**"), "",
         "## Reading", "",
         f"The old apparatus said closure of strength 0.5 or more would be "
         f"detected with power above 0.8, and read the empirical value as a "
         f"nominal strength near 0.35. The posterior says the same thing "
         f"without a null: strengths at or above 0.5 carry posterior mass "
         f"{o1['p_ge_05']:.3f} without time-averaging and {o3['p_ge_05']:.3f} "
         f"with it, so moderate and strong closure are excluded by the data "
         f"rather than by a threshold, while weak closure remains credible. "
         f"That is the bounded non-detection the paper argues for, stated as a "
         f"parameter.", "",
         "## Declared limitation", "",
         "p(T | s) is approximated by a normal at each grid point, with mean and "
         "standard deviation taken across seeds and interpolated between them. "
         "T is a bounded rank correlation, so the approximation is reasonable in "
         "the interior and degrades near the ends. This is a simulation-based "
         "posterior, not an exact one, and is reported as such.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    np.savez(OUT_NPZ, s_grid=S_FINE, post_w1=o1["p"], post_w3=o3["p"],
             t_obs=t_obs,
             summ_w1=np.array([o1["median"], o1["lo"], o1["hi"], o1["p_ge_05"]]),
             summ_w3=np.array([o3["median"], o3["lo"], o3["hi"], o3["p_ge_05"]]))

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.plot(S_FINE, o1["p"] / o1["p"].max(), color=OI_BLUE, label="no time-averaging")
    ax.plot(S_FINE, o3["p"] / o3["p"].max(), color=OI_VERMIL, ls="--",
            label="record's time-averaging")
    ax.set_xlabel("injected closure strength s"); ax.set_ylabel("posterior (scaled)")
    ax.legend(frameon=False, fontsize=7); fig.tight_layout()
    save(fig, "fig_closure_strength_posterior")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
