"""Beaumont, Zhang & Balding (2002) local-linear regression adjustment.

Corrects the residual dependence of the accepted parameters on the discrepancy
between simulated and observed summaries. Each particle is shifted along the
fitted local-linear surface to the point where the summary equals ``s_obs``,
which removes first-order ABC bias.
"""
from __future__ import annotations

import numpy as np


def _to_unbounded(x, lo, hi):
    """Logit-transform a box-constrained parameter onto the whole real line."""
    u = np.clip((x - lo) / (hi - lo), 1e-9, 1 - 1e-9)
    return np.log(u / (1.0 - u))


def _from_unbounded(z, lo, hi):
    return lo + (hi - lo) / (1.0 + np.exp(-z))


def regression_adjust(result, s_obs, param_indices=None, bounds=None) -> np.ndarray:
    thetas = result.thetas
    S = np.asarray(result.summaries, float)
    s_obs = np.asarray(s_obs, float)
    dists = np.asarray(result.distances, float)

    if param_indices is None:
        theta_cols = thetas
    else:
        theta_cols = thetas[:, list(param_indices)]

    # Epanechnikov kernel weight on the ABC distance (bandwidth = max distance).
    h = dists.max()
    if h <= 0:
        h = 1.0
    u = dists / h
    kw = np.clip(1.0 - u ** 2, 0.0, None)

    # SUPPORT. The local-linear shift is an EXTRAPOLATION and, applied on the
    # natural scale, routinely moves particles outside the prior's box. Measured
    # on this project's ABC-SMC transmission fit, 2026-09-02, 800 particles:
    # 300 of 800 innovation rates went negative, 150 of 800 population sizes
    # went negative (to -261), and 52 of 800 time-averaging windows left [5, 30].
    # Passing `bounds` performs the adjustment on a logit-transformed scale
    # instead, so every adjusted particle is inside the box by construction.
    # Omitting `bounds` preserves the original behavior.
    if bounds is not None:
        lo, hi = np.asarray(bounds[0], float), np.asarray(bounds[1], float)
        if param_indices is not None:
            lo, hi = lo[list(param_indices)], hi[list(param_indices)]
        if lo.shape[-1] != theta_cols.shape[1]:
            raise ValueError(
                f"bounds have {lo.shape[-1]} parameters but the particles have "
                f"{theta_cols.shape[1]}")
        work = _to_unbounded(theta_cols, lo, hi)
    else:
        work = theta_cols

    X = S - s_obs                                   # (n, s_dim)
    A = np.hstack([np.ones((X.shape[0], 1)), X])    # (n, 1 + s_dim)
    ATW = A.T * kw                                  # weight the normal equations
    coef, *_ = np.linalg.lstsq(ATW @ A, ATW @ work, rcond=None)
    beta = coef[1:]                                 # (s_dim, n_params)
    adjusted = work - X @ beta                      # shift to X = 0 (s == s_obs)
    if bounds is not None:
        adjusted = _from_unbounded(adjusted, lo, hi)
    return adjusted
