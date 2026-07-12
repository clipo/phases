"""Beaumont, Zhang & Balding (2002) local-linear regression adjustment.

Corrects the residual dependence of the accepted parameters on the discrepancy
between simulated and observed summaries. Each particle is shifted along the
fitted local-linear surface to the point where the summary equals ``s_obs``,
which removes first-order ABC bias.
"""
from __future__ import annotations

import numpy as np


def regression_adjust(result, s_obs, param_indices=None) -> np.ndarray:
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

    X = S - s_obs                                   # (n, s_dim)
    A = np.hstack([np.ones((X.shape[0], 1)), X])    # (n, 1 + s_dim)
    ATW = A.T * kw                                  # weight the normal equations
    coef, *_ = np.linalg.lstsq(ATW @ A, ATW @ theta_cols, rcond=None)
    beta = coef[1:]                                 # (s_dim, n_params)
    adjusted = theta_cols - X @ beta                # shift to X = 0 (s == s_obs)
    return adjusted
