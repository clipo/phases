"""Signature 3, the within- versus between-group partition of diversity.

Provides the plug-in Gini-Simpson diversity ``gini_simpson`` and the cultural
F_ST built from it, ``cultural_fst`` = (H_T - H_S) / H_T. Used by
``inference.bayesian_fst``, ``validation.harness``, ``make_figures.py``,
``scripts/measure_datum.py`` and analyses 05, 06, 07, 09, 33, 35, 44, 47, 57 and 62.
"""
from __future__ import annotations

import numpy as np


def gini_simpson(counts: np.ndarray) -> float:
    """Gini-Simpson diversity, 1 - sum(p^2).

    This is the PLUG-IN estimator and it is biased: E[sum(p_hat^2)] =
    sum(p^2) + (1 - sum(p^2))/N, so it returns D(1 - 1/N) and understates
    diversity by about D/N. The bias grows as the group shrinks, and since
    cultural_fst averages this over groups (small) and pools it (large), F_ST
    is biased upward by an amount that varies with group size.

    That is kept deliberately, and measured rather than assumed:
    `analyses/62_gini_simpson_bias.py` and
    `output/findings/gini_simpson_bias.md`. The short version is that the
    correction moves the descriptive trend from +0.004 to +0.044, but the
    inference maps the observed trend through a recovery curve scored with the
    SAME estimator, and the curve moves further than the observation does. The
    implied closure strength therefore goes DOWN, 0.37 to 0.31, both at the
    published posterior median near 0.30. The bias cannot be the reason the
    paper reaches its conclusion, because removing it strengthens that
    conclusion slightly.

    Use the unbiased form, 1 - sum(n(n-1))/(N(N-1)), for any NEW quantity
    reported as a diversity value in its own right rather than compared against
    a like-scored simulation.

    Raises on an empty assemblage rather than returning 0.0. Zero is a
    meaningful value here (one class present), so returning it for "no sherds"
    conflates an answer with the absence of one (rule 5, F4).
    """
    c = np.asarray(counts, float)
    if np.any(c < 0):
        raise ValueError("counts must be non-negative")
    N = c.sum()
    if N == 0:
        raise ValueError("gini_simpson is undefined for an empty assemblage")
    p = c / N
    return float(1.0 - np.sum(p ** 2))

def cultural_fst(group_counts: np.ndarray) -> float:
    """F_ST = (H_T - H_S)/H_T via Gini-Simpson. rows=groups, cols=types.
    H_S = size-weighted mean within-group diversity; H_T = total-pool diversity.

    Domain (rule 5, F4). Empty groups are dropped, since they carry zero weight
    in H_S and contribute nothing to the pool. Two boundaries are then refused
    rather than returned as 0.0:

      fewer than two non-empty groups   F_ST is not defined between one group.
      H_T == 0                          the pool holds a single class, so there
                                        is no diversity to partition.

    Both used to return 0.0, which is exactly the value "no differentiation"
    takes, so a degenerate input was indistinguishable from a real null and, in
    a trajectory, would have entered as a data point supporting the paper's own
    conclusion. Measured 2026-09-02: over 400 rarefactions at the analysis's
    common count, the pooled matrix is degenerate in 0 of 400, so no reported
    number was affected. The guard is for the next caller, not this one.
    """
    g = np.asarray(group_counts, float)
    if g.ndim != 2:
        raise ValueError(f"expected a 2-D groups-by-types array, got shape {g.shape}")
    if np.any(g < 0):
        raise ValueError("counts must be non-negative")
    sizes = g.sum(axis=1)
    g = g[sizes > 0]
    sizes = sizes[sizes > 0]
    if g.shape[0] < 2:
        raise ValueError(
            f"cultural_fst needs at least two non-empty groups, got {g.shape[0]}")
    H_T = gini_simpson(g.sum(axis=0))
    if H_T == 0:
        raise ValueError(
            "cultural_fst is undefined when the pooled assemblage holds a "
            "single class: there is no total diversity to partition")
    H_within = np.array([gini_simpson(row) for row in g])
    H_S = float(np.average(H_within, weights=sizes))
    return float((H_T - H_S) / H_T)


def moment_f(y: np.ndarray) -> np.ndarray:
    """Method-of-moments F per class across sites, each site weighted equally (not by its count), with the binomial sampling variance removed (per-type table, D27).

    Shares P = y / n per site; F_c = (var_i(P_ic, ddof 1) - mean_i(P_ic (1 - P_ic) / (n_i - 1))) / (pbar_c (1 - pbar_c)),
    pbar the unweighted mean share. The subtracted term is the unbiased
    estimate of the binomial sampling variance of a site share
    (E[P(1 - P)/(n - 1)] = p(1 - p)/n), so on identical site shares the
    estimate is near zero rather than near 1/n. NaN where a class is absent
    from every site or fixed at every site. Negative estimates are kept.
    Moved here from analyses/per_type/per_class_table.py on 2026-10-19 so
    the rim readings use the same definition.
    """
    y = np.asarray(y, float)
    if y.ndim != 2 or (y < 0).any() or not np.isfinite(y).all():
        raise ValueError("moment_f takes a non-negative site-by-class count table")
    n = y.sum(1)
    if (n <= 1).any():
        raise ValueError("every site needs at least two observations")
    P = y / n[:, None]; pb = P.mean(0)
    v = P.var(0, ddof=1); samp = np.mean(P * (1 - P) / (n - 1)[:, None], axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (v - samp) / (pb * (1 - pb))
