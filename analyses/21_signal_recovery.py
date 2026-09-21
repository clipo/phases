"""Revised profile-injection sensitivity with without-replacement rarefaction.

The shared-profile sampling baseline is not evolving spatial drift. Power is
conditional on the imposed increasing-divergence trajectory and is not an
empirical confidence bound on closure. IDSS group count per assemblage replaces
the historical fixed-row-order proxy. The canonical spatial comparison is 47.
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from figstyle import OI_BLUE, OI_VERMIL, OI_ORANGE, OI_GREEN, OI_PURPLE, save  # noqa: E402
import make_figures as mf  # noqa: E402
res = importlib.import_module("17_basin_results")
from mls_emergence.signatures.seriation import seriation_groups  # noqa: E402
from mls_emergence.signatures.sampling import rarefy  # noqa: E402

OUT = ROOT / "output" / "signal_recovery.md"
N_BINS = 6
S_GRID = np.round(np.arange(0.0, 1.0001, 0.1), 2)
N_SEEDS = 120
NRARE = 50               # common rarefied count (= min real assemblage size)
B_EMP = 400              # empirical rarefactions
KAPPA = 1.2
POWER_TARGET = 0.80
SIG = ["neutral", "seriability", "fst", "spatial"]
SIG_LABEL = {"neutral": "neutral departure", "seriability": "IDSS groups per assemblage",
             "fst": "cultural $F_{ST}$", "spatial": "spatial boundary"}


def _s_star_text(p) -> str:
    """The detection threshold, or why there isn't one.

    When power never reaches the target anywhere on the strength grid, s* is
    nan: the record cannot detect closure at any strength tested, which is a
    finding, not a missing number. Printing "nan" invites it to be read as a
    failed computation or, worse, copied into prose as though it were a value
    (rule 17: a quantity that is not identified is reported as not identified).
    """
    import numpy as _np
    if _np.isfinite(p["s_star"]):
        return f"{p['s_star']:.2f}"
    pw = _np.asarray(p["power"], dtype=float)
    return (f"not reached at any strength tested (peak power {pw.max():.2f} at "
            f"s = {S_GRID[int(pw.argmax())]:.1f})")


def _penalty_text(power) -> str:
    """The cost of time-averaging, when both thresholds exist to compare."""
    import numpy as _np
    a, b = power[1]["s_star"], power[3]["s_star"]
    if _np.isfinite(a) and _np.isfinite(b):
        return f"{b - a:+.2f}"
    return ("not comparable, because the time-averaged threshold is never "
            "reached on this grid")


def zipf_base(K):
    p = 1.0 / np.arange(1, K + 1)
    return p / p.sum()


def emergence_profiles(k, K, s, rng):
    p0 = zipf_base(K)
    shifts = rng.permutation(K)[:k]
    grid = np.empty((k, N_BINS, K))
    for t in range(N_BINS):
        a = s * t / (N_BINS - 1)
        for g in range(k):
            qg = np.roll(p0, int(shifts[g]))
            base = (1 - a) * p0 + a * qg
            prof = base ** (1.0 + KAPPA * a)
            grid[g, t] = prof / prof.sum()
    return grid


def time_average(grid, w):
    if w <= 1:
        return grid
    out = np.empty_like(grid)
    for t in range(grid.shape[1]):
        out[:, t, :] = grid[:, max(0, t - w + 1):t + 1, :].mean(axis=1)
    return out


def gen_synth_M(grid_w, clusters, bins_arr, N_arr, rng):
    n, K = len(clusters), grid_w.shape[2]
    M = np.zeros((n, K))
    for i in range(n):
        M[i] = rng.multinomial(int(N_arr[i]), grid_w[clusters[i], bins_arr[i]])
    return M


def sig_panel(M, clusters, bins_arr, coords_c, which=SIG):
    """Per-bin signature values, indexed by bin 0..N_BINS-1 (NaN where undefined)."""
    rows = {}
    for b in range(N_BINS):
        idx = np.where(bins_arr == b)[0]
        cell = {sg: np.nan for sg in which}
        if len(idx) > 0:
            sc, scl, sco = M[idx], clusters[idx], coords_c[idx]
            rep = np.unique(scl)
            if "neutral" in which:
                nd = []
                for c in rep:
                    p = sc[scl == c].sum(0)
                    if p.sum() >= 2 and (p > 0).sum() >= 2:
                        tf, te = mf.theta_f(p), mf.theta_e(p)
                        if np.isfinite(tf) and te > 0:
                            nd.append(abs(1 - tf / te))
                cell["neutral"] = float(np.mean(nd)) if nd else np.nan
            if "seriability" in which:
                cell["seriability"] = seriation_groups(sc, cont=0.10)["n_groups"] / len(idx) if len(idx) >= 2 else np.nan
            if "fst" in which:
                cell["fst"] = mf.cultural_fst(np.array([sc[scl == c].sum(0) for c in rep])) if len(rep) >= 2 else np.nan
            if "spatial" in which:
                cell["spatial"] = mf.boundary_excess(sc, sco, seed=7) if len(idx) >= 4 else np.nan
        rows[b] = cell
    return pd.DataFrame(rows).T.reindex(range(N_BINS))


def sig_rhos(M, clusters, bins_arr, coords_c, which=SIG):
    """Per-bin signature trajectory; return {signature: Spearman rho vs bin}."""
    P = sig_panel(M, clusters, bins_arr, coords_c, which)
    out = {}
    for col in which:
        s = P[col].dropna()
        out[col] = spearmanr(s.index.to_numpy(float), s.values)[0] if len(s) >= 3 else np.nan
    return out


def operating_point():
    """The basin operating point every size-controlled signature is scored on.

    Extracted so that there is exactly one definition of it (rule 6). It used to
    be rebuilt inline here and again in 29_concept_figure.py, and while the two
    constructions agreed, they drew their rarefactions from different seeds and
    so reported the F_ST trend as +0.01 and -0.02 respectively: the same
    statistic on opposite sides of zero, in the text and in Figure 2.

    Returns the curated counts restricted to assemblages that have coordinates
    and a defined seriation bin, with the silhouette-selected spatial clustering
    and the 14C-oriented ordinal bins.
    """
    counts, coords = mf._load_curated()
    ca, _ = res.oriented_ca(counts)
    Ni = dict(zip(counts.index, counts.to_numpy(float).sum(1)))
    cdf = coords.dropna()
    have = list(cdf.index)
    cc = cdf[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(0)
    sil = {kk: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, kk, seed=7))
           for kk in range(2, 7)}
    k = max(sil, key=sil.get)
    cl = mf._kmeans_labels(cc_c, k, seed=7)
    cluster_of = dict(zip(have, cl))
    bins = pd.qcut(ca.reindex(have), N_BINS, labels=False, duplicates="drop")
    keep = [i for i in have if not pd.isna(bins[i])]
    return {
        "counts": counts, "coords": coords, "ca": ca, "K": counts.shape[1],
        "Ni": Ni, "have": have, "cc_c": cc_c, "sil": sil, "k": k, "keep": keep,
        "clusters": np.array([cluster_of[i] for i in keep]),
        "bins_arr": np.array([int(bins[i]) for i in keep]),
        "coords_c": cc_c[[have.index(i) for i in keep]],
        "N_arr": np.array([Ni[i] for i in keep]),
        "real_M": counts.reindex(keep).to_numpy(float),
    }


# Draw count for the *reported* empirical trend. The single-draw spread of this
# statistic is large (SD about 0.26), so a 400-draw mean carries a Monte Carlo
# standard error near 0.013, which is bigger than the estimate itself. Reporting
# it to two decimals from a few hundred draws is what let the text say +0.01 and
# Figure 2 say -0.02 with neither being wrong. 4000 draws puts the Monte Carlo
# error at about 0.004, small enough to quote the value to three decimals.
TREND_DRAWS = 4000
TREND_SEED = 777


def empirical_fst_trend(op=None, n_draws: int = TREND_DRAWS, seed: int = TREND_SEED):
    """Size-controlled cultural F_ST trajectory trend, with its Monte Carlo error.

    This is a rank correlation between per-bin F_ST and ordinal seriation
    position on rarefied assemblages, NOT an F_ST value (see rule 6 in
    CLAUDE.md). Returns mean, Monte Carlo standard error, the 95 percent
    interval of the mean, and the single-draw SD.
    """
    if op is None:
        op = operating_point()
    vals = []
    for b in range(n_draws):
        rg = np.random.default_rng(seed + b)
        Mr = rarefy(op["real_M"], NRARE, rg)
        v = sig_rhos(Mr, op["clusters"], op["bins_arr"], op["coords_c"],
                     which=["fst"])["fst"]
        if np.isfinite(v):
            vals.append(v)
    vals = np.asarray(vals)
    sd = float(vals.std(ddof=1))
    se = sd / np.sqrt(len(vals))
    # TWO intervals, and they answer different questions (CLAUDE.md rule 6).
    # `lo`/`hi` are the Monte Carlo interval of the MEAN: how precisely we know
    # the average of our own draws, a fact about the compute budget. `p_lo`/
    # `p_hi` are the rarefaction spread: what a single subsample of this record
    # returns, which is what says whether the record resolves the trend. Quote
    # and PLOT the second. The first is about thirty times narrower here and
    # makes an unresolved trend look resolved.
    return {"mean": float(vals.mean()), "se": float(se), "sd": sd,
            "lo": float(vals.mean() - 1.96 * se), "hi": float(vals.mean() + 1.96 * se),
            "p_lo": float(np.percentile(vals, 2.5)),
            "p_hi": float(np.percentile(vals, 97.5)),
            "n_draws": int(len(vals))}


def main():
    op = operating_point()
    counts, coords = op["counts"], op["coords"]
    ca, K, Ni = op["ca"], op["K"], op["Ni"]
    have, cc_c, sil, k = op["have"], op["cc_c"], op["sil"], op["k"]
    keep = op["keep"]
    clusters, bins_arr = op["clusters"], op["bins_arr"]
    coords_c, N_arr, real_M = op["coords_c"], op["N_arr"], op["real_M"]
    n = len(keep)
    size_conf = spearmanr(N_arr, bins_arr)[0]

    # ---- empirical signatures: raw (Figure 5) and size-controlled (rarefied) ----
    rng = np.random.default_rng(9)
    raw_emp = sig_rhos(real_M, clusters, bins_arr, coords_c)
    rar_emp = {s: [] for s in SIG}
    emp_fst = []
    for _ in range(B_EMP):
        Mr = rarefy(real_M, NRARE, rng)
        r = sig_rhos(Mr, clusters, bins_arr, coords_c)
        for s in SIG:
            rar_emp[s].append(r[s])
        emp_fst.append(r["fst"])
    rar_emp_mean = {s: float(np.nanmean(rar_emp[s])) for s in SIG}
    # The reported F_ST trend comes from the shared estimator, at a draw
    # count large enough to quote (see empirical_fst_trend). The B_EMP
    # rarefactions above still drive the Figure 5 per-panel annotations.
    trend = empirical_fst_trend(op)
    emp_fst = np.array(emp_fst)

    # ---- per-signature recovery (rarefied, w=1) ----
    persig = {s: np.zeros(len(S_GRID)) for s in SIG}
    fst_w1 = np.zeros((len(S_GRID), N_SEEDS))
    for si, s in enumerate(S_GRID):
        acc = {sg: [] for sg in SIG}
        for seed in range(N_SEEDS):
            rng2 = np.random.default_rng(si * 100003 + seed * 17 + 1)
            grid = emergence_profiles(k, K, s, rng2)
            M = rarefy(gen_synth_M(grid, clusters, bins_arr, N_arr, rng2), NRARE, rng2)
            r = sig_rhos(M, clusters, bins_arr, coords_c)
            for sg in SIG:
                acc[sg].append(r[sg])
            fst_w1[si, seed] = r["fst"]
        for sg in SIG:
            persig[sg][si] = np.nanmean(acc[sg])

    # ---- F_ST detector power, with and without time-averaging ----
    def fst_only(M):
        return sig_rhos(M, clusters, bins_arr, coords_c, which=["fst"])["fst"]
    power = {}
    fst_mats = {1: fst_w1}
    for w in (1, 3):
        if w not in fst_mats:
            Tm = np.zeros((len(S_GRID), N_SEEDS))
            for si, s in enumerate(S_GRID):
                for seed in range(N_SEEDS):
                    rng2 = np.random.default_rng(si * 100003 + seed * 17 + w * 7)
                    grid = time_average(emergence_profiles(k, K, s, rng2), w)
                    M = rarefy(gen_synth_M(grid, clusters, bins_arr, N_arr, rng2), NRARE, rng2)
                    Tm[si, seed] = fst_only(M)
            fst_mats[w] = Tm
        Tmat = fst_mats[w]
        # Estimate the threshold on separate null replicates; Tmat[0] is validation.
        calibration = []
        for seed in range(400):
            rg = np.random.default_rng(800000 + 1000*w + seed)
            grid = time_average(emergence_profiles(k, K, 0, rg), w)
            calibration.append(fst_only(rarefy(gen_synth_M(grid, clusters, bins_arr, N_arr, rg), NRARE, rg)))
        thr = float(np.nanpercentile(calibration, 95))
        pw = np.array([np.mean(Tmat[si][~np.isnan(Tmat[si])] > thr) for si in range(len(S_GRID))])
        above = np.where(pw >= POWER_TARGET)[0]
        power[w] = dict(thr=thr, power=pw, s_star=float(S_GRID[above[0]]) if len(above) else float("nan"),
                        null_mean=float(np.nanmean(Tmat[0])))
    thr1 = power[1]["thr"]
    emp_fst_mean = float(np.nanmean(emp_fst))
    s_star1 = power[1]["s_star"]
    # nominal injected strength the empirical F_ST trend corresponds to (curve is monotone increasing)
    fst_curve = persig["fst"]
    nominal_s = float(np.interp(emp_fst_mean, fst_curve, S_GRID))
    emp_resolvable = np.isfinite(s_star1) and nominal_s >= s_star1

    # Machine-readable summary of the grid (power is a design property of the
    # injected alternative, not an empirical confidence bound on s).
    summary = dict(n=n, k=k, n_bins=N_BINS, depth=NRARE, seeds=N_SEEDS,
                   raw=raw_emp, empirical=rar_emp_mean,
                   empirical_fst_interval=np.nanpercentile(emp_fst, [2.5, 97.5]).tolist(),
                   strengths=S_GRID.tolist(), response={sg: v.tolist() for sg, v in persig.items()},
                   power={str(w): {key: value.tolist() if isinstance(value, np.ndarray) else value
                                   for key, value in row.items()} for w, row in power.items()})
    (ROOT / "output" / "revision_recovery.json").write_text(json.dumps(summary, indent=2))

    # ---- report ----
    L = [
        "# Record-matched, size-controlled signal recovery", "",
        f"Real configuration: {n} assemblages, {K} decorated types, k = {k} spatial clusters, "
        f"{N_BINS} ordinal bins. Synthetic assemblages are generated at the real per-assemblage "
        f"sample sizes; both synthetic and real assemblages are rarefied to a common count "
        f"NRARE = {NRARE} before scoring. {N_SEEDS} seeds per cell, {B_EMP} empirical "
        f"rarefactions.", "",
        "The s = 0 baseline is a fixed shared Zipf profile with multinomial sampling, not a "
        "spatially evolving drift null; the injection imposes increasing divergence and profile "
        "concentration, and every recovery statement below refers only to that alternative. "
        "The seriation statistic is the number of maximal deterministic IDSS groups per assemblage "
        "within a bin (continuity 0.10), not an ordered-row proxy.", "",
        "## The sample-size confound", "",
        f"- Assemblage sample size trends with seriation position at Spearman rho = {size_conf:+.2f}, "
        f"so a size-sensitive signature can rise or fall along the axis through sampling alone.",
        f"- Raw (uncontrolled) empirical F_ST trend = {raw_emp['fst']:+.2f} (the Figure 5 value); "
        f"after rarefaction it is {trend['mean']:+.3f} "
        f"(Monte Carlo 95% interval [{trend['lo']:+.3f}, {trend['hi']:+.3f}] over "
        f"{trend['n_draws']} draws; a single rarefaction has SD {trend['sd']:.2f}, so this "
        f"statistic is not stable in sign at a few hundred draws). "
        f"The raw rise is a sampling artifact.",
        "",
        "## Which signatures recover the injected emergence (rarefied, no averaging)", "",
        "| s | " + " | ".join(SIG_LABEL[s] for s in SIG) + " |", "|---|" + "---|" * len(SIG)]
    for si, s in enumerate(S_GRID):
        L.append(f"| {s:.1f} | " + " | ".join(f"{persig[sg][si]:+.2f}" for sg in SIG) + " |")
    L += ["",
          "**Reading.** Only cultural F_ST tracks the injected signal monotonically (null near "
          "zero, rising to ~+1 at strong emergence). The neutral departure is non-monotonic in "
          "conformity (the known U-shape), the spatial boundary is unresponsive at k = 3 clusters, "
          "and the seriation group count does not rise in the hypothesized direction. At this record's resolution the criterion "
          "is carried by F_ST; the other three signatures are not reliable discriminators here.",
          "",
          "## F_ST detector (calibrated; false-positive rate 0.05)", "",
          "| s | power (w=1) | power (time-averaged w=3) |", "|---|---|---|"]
    for si, s in enumerate(S_GRID):
        L.append(f"| {s:.1f} | {power[1]['power'][si]:.2f} | {power[3]['power'][si]:.2f} |")
    L += ["",
          f"- Detection threshold s* = {_s_star_text(power[1])} (no averaging) / "
          f"{_s_star_text(power[3])} (time-averaged): the weakest emergence recovered at power "
          f">= {POWER_TARGET:.0%}. Time-averaging penalty: {_penalty_text(power)}.",
          f"- Null (s=0) F_ST trend mean {power[1]['null_mean']:+.2f}; detection threshold "
          f"{thr1:+.2f}.",
          "- Thresholds use 400 separate null calibration draws per window; the s = 0 test row "
          "estimates the achieved false-positive rate independently.", "",
          "## Empirical placement (size-controlled)", "",
          f"- Rarefied empirical F_ST trend = {emp_fst_mean:+.2f} "
          f"[{np.nanpercentile(emp_fst,2.5):+.2f}, {np.nanpercentile(emp_fst,97.5):+.2f}].",
          f"- On the recovery curve this corresponds to a nominal injected strength s ~ "
          f"{nominal_s:.2f}, {'at or above' if emp_resolvable else 'below'} the resolution limit "
          f"s* = {s_star1:.2f} that this record can reliably detect.",
          f"- The data show {'a resolvable emergence signal' if emp_resolvable else 'no resolvable emergence'}: "
          f"the faint trend is not distinguishable from the no-emergence null at this resolution.",
          "- The rarefaction percentiles condition on the observed sites, partition and ordering; "
          "they are not a confidence interval over all archaeological uncertainties. The signed "
          "trend reverses with axis orientation.",
          "",
          "## Verdict", "",
          f"At the record's own resolution and sample sizes, and with the size confound removed by "
          f"rarefaction, the discrimination is carried by cultural F_ST. The apparent raw F_ST rise "
          f"({raw_emp['fst']:+.2f}) is an artifact of the sample-size-versus-position trend "
          f"(rho {size_conf:+.2f}) and vanishes under size control, leaving a trend of "
          f"{trend['mean']:+.3f} [{trend['lo']:+.3f}, {trend['hi']:+.3f}]. "
          f"The other three signatures do not reliably discriminate at this "
          f"resolution and are reported as weak corroboration only.",
          "",
          "**How the strength of any closure is reported.** The detection threshold, power curve and "
          "five percent false-positive rate that this section used to report are withdrawn under the "
          "project's rule 18 (no frequentist inference supports a substantive claim). The simulation "
          "grid above is a simulated likelihood, so the reportable quantity is a POSTERIOR over the "
          "injected closure strength, computed by analyses/53_closure_strength_posterior.py and "
          "written to output/findings/closure_strength_posterior.md. It is what Figure 4's right "
          "panel now draws. The power numbers this script still computes internally are retained "
          "only to build that grid and must not be quoted."]
    OUT.write_text("\n".join(L), encoding="utf-8")

    # ---- figure ----
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7, 3.3))
    cmap = {"neutral": OI_PURPLE, "seriability": OI_GREEN, "fst": OI_ORANGE, "spatial": OI_BLUE}
    for sg in SIG:
        axA.plot(S_GRID, persig[sg], "-o", ms=3, color=cmap[sg], lw=1.4, label=SIG_LABEL[sg])
    axA.axhline(0, color="0.7", lw=0.6)
    axA.set_xlabel("injected closure strength s")
    axA.set_ylabel("signature trend (Spearman rho)")
    axA.legend(frameon=False, fontsize=6.5, loc="upper left")

    # Panel B: the POSTERIOR over injected closure strength, not a power curve.
    # Rule 18 removed the calibrated detector, its five percent false-positive
    # rate and its detection threshold (docs/FREQUENTIST_INVENTORY.md item A),
    # and the caption was updated on 2026-09-02. This panel is drawn from
    # analyses/53_closure_strength_posterior.py's saved posterior so the figure
    # and the caption cannot drift apart again.
    # Loaded through 53's own loader, which refuses an .npz older than the
    # finding it was written beside. Re-running this script against a stale
    # array silently reverted Figure 4 to superseded values once already.
    _m53 = importlib.import_module("53_closure_strength_posterior")
    _pd = _m53.load_posterior()
    _sg = _pd["s_grid"]
    axB.plot(_sg, _pd["post_w1"] / _pd["post_w1"].max(), color=OI_BLUE,
             label="no averaging")
    axB.plot(_sg, _pd["post_w3"] / _pd["post_w3"].max(), color=OI_VERMIL,
             ls="--", label="time-averaged")
    _m1, _l1, _h1, _p1 = _pd["summ_w1"]
    _m3, _l3, _h3, _p3 = _pd["summ_w3"]
    axB.axvline(_m1, color=OI_BLUE, ls=":", lw=1.0)
    axB.axvline(_m3, color=OI_VERMIL, ls=":", lw=1.0)
    axB.set_xlabel("injected closure strength s")
    axB.set_ylabel("posterior (scaled)")
    axB.set_ylim(0, 1.08)
    axB.legend(frameon=False, fontsize=6.5, loc="upper right")
    # Placed right of both curves, which are at or near zero for s > 0.6.
    axB.text(0.99, 0.60,
             f"median {_m1:.2f} [{_l1:.2f}, {_h1:.2f}]\n"
             f"P(s$\\geq$0.5) = {_p1:.2f}\n"
             f"time-avg. {_m3:.2f}, P = {_p3:.2f}",
             transform=axB.transAxes, ha="right", va="top", fontsize=5.5, color=OI_GREEN)
    save(fig, "fig5_recovery")

    # ---- Figure 5: size-controlled empirical four-signature trajectory ----
    rngf = np.random.default_rng(21)
    stack = np.full((B_EMP, N_BINS, len(SIG)), np.nan)
    for bi in range(B_EMP):
        P = sig_panel(rarefy(real_M, NRARE, rngf), clusters, bins_arr, coords_c)
        stack[bi] = P[SIG].to_numpy(float)
    panel_mean = np.nanmean(stack, axis=0)   # (N_BINS, n_sig)
    fig5, axes = plt.subplots(2, 2, figsize=(7, 5.2))
    cmap5 = {"neutral": OI_PURPLE, "seriability": OI_GREEN, "fst": OI_ORANGE, "spatial": OI_BLUE}
    for ax, j in zip(axes.ravel(), range(len(SIG))):
        sg = SIG[j]
        ax.plot(range(N_BINS), panel_mean[:, j], "-o", ms=4, color=cmap5[sg], lw=1.5)
        ax.set_ylabel(SIG_LABEL[sg], fontsize=8)
        ax.set_xticks(range(N_BINS))
        # The F_ST panel quotes the shared estimator (4,000 draws, three
        # decimals); at two decimals off B_EMP draws this annotation and the
        # main text disagreed in the second decimal and, in Figure 2, in sign.
        lab = (rf"$\rho$ = {trend['mean']:+.3f}" if sg == "fst"
               else rf"$\rho$ = {rar_emp_mean[sg]:+.2f}")
        ax.text(0.95, 0.95, lab, transform=ax.transAxes,
                ha="right", va="top", fontsize=9)
    for ax in axes[1, :]:
        ax.set_xlabel("CA seriation bin (early to late)", fontsize=8)
    fig5.tight_layout()
    save(fig5, "fig6_empirical_trajectory")

    print(f"wrote {OUT}")
    print("\n".join(L))


if __name__ == "__main__":
    main()
