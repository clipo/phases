from __future__ import annotations
import numpy as np
from scipy.optimize import brentq


def homozygosity_f(counts: np.ndarray) -> float:
    """Unbiased homozygosity F = sum n_i(n_i-1) / (N(N-1))."""
    c = np.asarray(counts, float)
    N = c.sum()
    if N < 2:
        return float("nan")
    return float(np.sum(c * (c - 1)) / (N * (N - 1)))


def theta_f(counts: np.ndarray) -> float:
    """Neiman t_F: F = 1/(1+theta) => theta=(1-F)/F."""
    F = homozygosity_f(counts)
    if not np.isfinite(F) or F <= 0:
        return float("inf")
    return (1.0 - F) / F


def _ewens_expected_k(theta: float, n: int) -> float:
    i = np.arange(n)
    return float(np.sum(theta / (theta + i)))


def theta_e(counts: np.ndarray) -> float:
    """Ewens t_E: solve E[k|theta,n] = observed k for theta."""
    c = np.asarray(counts)
    n = int(c.sum())
    k = int((c > 0).sum())
    if k <= 1:
        return 0.0
    def f(t):
        return _ewens_expected_k(t, n) - k
    return float(brentq(f, 1e-9, 1e9))


def neutrality_verdict(counts: np.ndarray) -> dict:
    """Compare t_F and t_E. t_F<t_E => conformist; t_F>t_E => anti-conformist."""
    tf, te = theta_f(counts), theta_e(counts)
    if not np.isfinite(tf):
        return {"theta_f": tf, "theta_e": te, "direction": "undefined"}
    ratio = tf / te if te > 0 else float("inf")
    direction = "neutral"
    if ratio < 0.8:
        direction = "conformist"
    elif ratio > 1.25:
        direction = "anti-conformist"
    return {"theta_f": tf, "theta_e": te, "ratio": ratio, "direction": direction}
