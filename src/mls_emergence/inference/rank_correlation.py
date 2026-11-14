"""Bayesian rank correlation: the replacement for Spearman rho plus a p-value.

Rule 18 (binding here from 2026-08-31) removes p-values as a basis for
inference. Several claims in this project rest on a Spearman rho reported with
"not significant", which is a statement about a null this project does not
believe rather than about the association it is asking after.

The construction. Both variables are mapped to NORMAL SCORES (van der Waerden),
z = Phi^-1((rank - 0.5) / n), and a bivariate normal with unit marginals is
fitted to the pairs. Its correlation parameter r is the estimand. This is the
standard Bayesian analogue of Spearman's rho: it is rank-based, so it inherits
Spearman's robustness to monotone transformation and to outliers, and it returns
a posterior with an interval instead of a point estimate with a p-value.

PRIOR AND ITS DIRECTION (rule 20). r ~ Uniform(-1, 1). This is symmetric about
zero, so it leans neither toward nor against any association, and the burden of
defence is correspondingly light. It is the one place in this project where
"uninformative" is defensible without a prior predictive argument, because the
parameter is bounded, the prior is flat on the whole of its support, and the
quantity has no scale to get wrong.

WHAT A SMALL n LOOKS LIKE HERE, and why that is the point. With six seriation
bins the posterior on r is very wide. The frequentist version reported that same
weakness as "not significant", which reads as evidence of absence. A wide
posterior says what is actually true: the data barely constrain the association.
"""
from __future__ import annotations

import numpy as np


def normal_scores(x) -> np.ndarray:
    """van der Waerden scores: Phi^-1((rank - 0.5) / n). Ties averaged."""
    from scipy.stats import norm, rankdata
    x = np.asarray(x, float)
    if x.ndim != 1:
        raise ValueError(f"expected a 1-D array, got shape {x.shape}")
    ok = np.isfinite(x)
    if ok.sum() < 3:
        raise ValueError(f"need at least 3 finite values, got {int(ok.sum())}")
    out = np.full(x.shape, np.nan)
    r = rankdata(x[ok])
    out[ok] = norm.ppf((r - 0.5) / ok.sum())
    return out


def bayesian_rank_correlation(x, y, *, draws: int = 2000, tune: int = 2000,
                              chains: int = 4, random_seed: int = 0,
                              target_accept: float = 0.95) -> dict:
    """Posterior for the normal-scores rank correlation between x and y.

    Returns median, 95% credible interval, P(r > 0), and the full rule-16
    diagnostic set. Diagnostics are attached here by construction so no caller
    can report this correlation without them.
    """
    import arviz as az
    import pymc as pm

    zx, zy = normal_scores(x), normal_scores(y)
    ok = np.isfinite(zx) & np.isfinite(zy)
    z = np.column_stack([zx[ok], zy[ok]])
    n = int(ok.sum())
    if n < 3:
        raise ValueError(f"need at least 3 complete pairs, got {n}")

    with pm.Model():
        r = pm.Uniform("r", -1.0, 1.0)
        cov = pm.math.stack([pm.math.stack([1.0, r]),
                             pm.math.stack([r, 1.0])]).reshape((2, 2))
        pm.MvNormal("obs", mu=np.zeros(2), cov=cov, observed=z)
        idata = pm.sample(draws=draws, tune=tune, chains=chains, cores=1,
                          random_seed=random_seed, target_accept=target_accept,
                          progressbar=False, compute_convergence_checks=False)

    s = np.asarray(idata.posterior["r"].values).ravel()
    div = int(idata.sample_stats["diverging"].sum())
    tot = int(idata.sample_stats["diverging"].size)
    td = idata.sample_stats.get("tree_depth")
    e = np.asarray(idata.sample_stats["energy"].values)
    ebfmi = float(np.min((np.diff(e, axis=1) ** 2).mean(axis=1) / e.var(axis=1)))
    return {
        "n": n,
        "r_median": float(np.median(s)),
        "r_mean": float(np.mean(s)),
        "hdi95": [float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))],
        "p_positive": float((s > 0).mean()),
        "rhat": float(np.asarray(az.rhat(idata, var_names=["r"])["r"].values)),
        "ess_bulk": float(np.asarray(az.ess(idata, var_names=["r"])["r"].values)),
        "ess_tail": float(np.asarray(
            az.ess(idata, var_names=["r"], method="tail")["r"].values)),
        "n_div": div, "n_total": tot, "pct_div": 100.0 * div / tot,
        "n_treedepth": int((td >= 10).sum()) if td is not None else -1,
        "ebfmi": ebfmi,
        "samples": s,
    }


def format_result(res: dict, label: str = "") -> str:
    """One-line report, diagnostics included by construction (rule 16)."""
    return (f"{label}r = {res['r_median']:+.3f} "
            f"[{res['hdi95'][0]:+.3f}, {res['hdi95'][1]:+.3f}], "
            f"P(r > 0) = {res['p_positive']:.3f}, n = {res['n']} "
            f"| R-hat {res['rhat']:.4f}, bulk ESS {res['ess_bulk']:.0f}, "
            f"tail ESS {res['ess_tail']:.0f}, divergences {res['n_div']}/"
            f"{res['n_total']} ({res['pct_div']:.2f}%), "
            f"treedepth>=10 {res['n_treedepth']}, E-BFMI {res['ebfmi']:.3f}")
