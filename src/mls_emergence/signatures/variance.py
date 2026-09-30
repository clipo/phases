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
