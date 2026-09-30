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
    """Ewens t_E: solve E[k | theta, n] = observed k for theta.

    Domain (rule 5). E[k | theta, n] is strictly increasing in theta from 1
    (as theta -> 0) toward n (as theta -> infinity), so a finite root exists
    only for 1 < k < n. The two boundaries are refused rather than returned as
    a number:

      k == 1  every sherd one class. theta is 0 in the limit, and 0.0 is
              returned because it is the correct limiting value, not a guess.
      k == n  every sherd a distinct class. E[k] approaches n only as theta
              diverges, so there is no finite root; brentq used to raise a bare
              ValueError from inside scipy with no indication of the cause.
    """
    c = np.asarray(counts)
    if np.any(c < 0):
        raise ValueError("counts must be non-negative")
    n = int(c.sum())
    k = int((c > 0).sum())
    if n < 1:
        raise ValueError("theta_e needs at least one sherd; got an empty assemblage")
    if k <= 1:
        return 0.0
    if k >= n:
        raise ValueError(
            f"theta_e is not finite when every sherd is a distinct class "
            f"(k = {k}, n = {n}): E[k | theta, n] approaches n only as theta "
            f"diverges. Caller must decide how to treat this assemblage.")

    def f(t):
        return _ewens_expected_k(t, n) - k

    return float(brentq(f, 1e-9, 1e9))


def neutrality_verdict(counts: np.ndarray) -> dict:
    """Compare t_F and t_E, DESCRIPTIVELY. t_F < t_E leans conformist.

    NOT AN INFERENCE (F5, and rule 18). The 0.8 and 1.25 cut points below have
    no source in the literature and no sampling distribution behind them: they
    are a reading aid for the ratio, not a test, and the returned ``direction``
    string must not be quoted as evidence that transmission was conformist or
    anti-conformist. The manuscript does not use it, and the recovery
    experiment shows this signature does not discriminate at the record's
    resolution.

    Where a transmission claim is actually made, it rests on the ABC-SMC
    posterior for the bias parameter (analyses/38) and the hierarchical
    increment model (analyses/52), both of which report an interval.
    """
    tf, te = theta_f(counts), theta_e(counts)
    if not np.isfinite(tf):
        return {"theta_f": tf, "theta_e": te, "direction": "undefined"}
    ratio = tf / te if te > 0 else float("inf")
    direction = "neutral"
    if ratio < 0.8:
        direction = "conformist"
    elif ratio > 1.25:
        direction = "anti-conformist"
    return {"theta_f": tf, "theta_e": te, "ratio": ratio,
            "direction": direction,
            # Carried in the payload so the caveat travels with the value and
            # cannot be dropped by quoting the dict (rule 6).
            "direction_is_descriptive": True,
            "note": ("descriptive label from an uncalibrated ratio; not "
                     "evidence of a transmission regime. See F5.")}
