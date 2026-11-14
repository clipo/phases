"""Decorrelation lengths of the copying simulator and the stationarity check (spec 3.3; Finding 8).

L_d: the local decorrelation length. For each seed, site and class, the
sample autocorrelation of the stored frequency series at lags 0..S/4 with the
standard estimator (lagged products over the full-series variance, divisor
S); series with sample variance below the floor are excluded and counted;
the included curves are averaged with equal weight; L_d for the seed is the
smallest lag at which the average is at or below 1/e, with a minimum of one
step. L_R: the same on the regional-mean series per class, averaged over
classes. Both are the median over seeds.

Window rule (correction to spec 3.3, 2026-10-16). The biased estimator
carries a taper of about (1 - h/S), so any smooth slow series crosses 1/e by
roughly 0.22 S whatever its true time scale (measured on a linear trend, a
quarter cosine and a square root: crossings at 0.21 to 0.22 S for S from 100
to 1,600). "No crossing by S/4" therefore never triggers, and a short window
would report a window artifact as the time scale. A crossing is accepted
only when it lies at or below S/10, where the taper is at most ten percent;
otherwise the stored length is doubled, up to three times.
"""
from __future__ import annotations

import math

import numpy as np

LEVEL = 1.0 / math.e
VARIANCE_FLOOR = 1e-6


def acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Sample autocorrelation at lags 0..max_lag, divisor S (biased estimator)."""
    x = np.asarray(x, float)
    S = len(x)
    if max_lag >= S:
        raise ValueError("max_lag must be smaller than the series length")
    xc = x - x.mean()
    denom = float((xc * xc).sum())
    if denom <= 0:
        raise ValueError("series has zero variance")
    return np.array([float((xc[:S - h] * xc[h:]).sum()) / denom for h in range(max_lag + 1)])


def mean_local_curve(freq: np.ndarray, floor: float = VARIANCE_FLOOR, max_lag: int | None = None):
    """Equal-weight mean ACF over (site, class) series; returns (curve, n_included, n_excluded)."""
    freq = np.asarray(freq, float)
    S, N, K = freq.shape
    max_lag = S // 4 if max_lag is None else int(max_lag)
    curves = []; n_ex = 0
    for i in range(N):
        for k in range(K):
            x = freq[:, i, k]
            if x.var() < floor:
                n_ex += 1
                continue
            curves.append(acf(x, max_lag))
    if not curves:
        return None, 0, n_ex
    return np.mean(curves, axis=0), len(curves), n_ex


def regional_curve(freq: np.ndarray, floor: float = VARIANCE_FLOOR, max_lag: int | None = None):
    """Mean ACF over classes of the regional-mean series (unweighted mean over sites)."""
    freq = np.asarray(freq, float)
    S = freq.shape[0]
    max_lag = S // 4 if max_lag is None else int(max_lag)
    reg = freq.mean(axis=1)                      # S x K
    curves = []; n_ex = 0
    for k in range(reg.shape[1]):
        x = reg[:, k]
        if x.var() < floor:
            n_ex += 1
            continue
        curves.append(acf(x, max_lag))
    if not curves:
        return None, 0, n_ex
    return np.mean(curves, axis=0), len(curves), n_ex


def first_crossing(curve, level: float = LEVEL):
    """Smallest lag at which the curve is at or below `level`, minimum one step; None if never."""
    if curve is None:
        return None
    below = np.where(np.asarray(curve)[1:] <= level)[0]
    return int(below[0]) + 1 if len(below) else None


MAX_FRACTION = 0.10


def measure_time_scales(run, seeds, S0: int = 5000, max_doublings: int = 3, floor: float = VARIANCE_FLOOR,
                        max_fraction: float = MAX_FRACTION) -> dict:
    """Measure L_d and L_R with the doubling rule.

    `run(seed, S)` returns stored frequencies of shape (S, N, K). The stored
    length starts at S0 and doubles, at most `max_doublings` times, until
    every seed's local and regional curves cross 1/e at a lag no larger than
    `max_fraction` of S (see the module docstring). Returns the per-seed
    values, medians, the curves and the stored length used; the key "status"
    is "ok", "no_crossing" or "all_excluded".
    """
    S = int(S0); doublings = 0
    accepted = lambda L, S: L is not None and L <= max_fraction * S
    while True:
        out = {"S": S, "doublings": doublings, "seeds": list(map(int, seeds)), "local": [], "regional": [],
               "n_included_local": [], "n_excluded_local": [], "n_excluded_regional": [], "curves_local": [], "curves_regional": []}
        for sd in seeds:
            f = run(int(sd), S)
            cl, n_in, n_ex = mean_local_curve(f, floor=floor)
            cr, _, n_exr = regional_curve(f, floor=floor)
            out["curves_local"].append(None if cl is None else cl.tolist())
            out["curves_regional"].append(None if cr is None else cr.tolist())
            out["n_included_local"].append(n_in); out["n_excluded_local"].append(n_ex); out["n_excluded_regional"].append(n_exr)
            out["local"].append(first_crossing(cl)); out["regional"].append(first_crossing(cr))
        if any(c is None for c in out["curves_local"]) or any(c is None for c in out["curves_regional"]):
            out["status"] = "all_excluded"; return out
        out["max_fraction"] = max_fraction
        if all(accepted(v, S) for v in out["local"]) and all(accepted(v, S) for v in out["regional"]):
            out["status"] = "ok"
            out["L_d"] = int(np.median(out["local"])); out["L_R"] = int(np.median(out["regional"]))
            return out
        if doublings >= max_doublings:
            out["status"] = "no_crossing"; return out
        S *= 2; doublings += 1


def between_site_dispersion(freq: np.ndarray) -> np.ndarray:
    """Per stored step: variance across sites of the class frequencies, averaged over classes."""
    f = np.asarray(freq, float)
    return f.var(axis=1).mean(axis=1)


def stationarity(freq_by_seed) -> dict:
    """The half-run rule of spec 3.3.

    For each seed, the mean between-site dispersion in the first and second
    halves of the stored run. Non-stationary if the absolute difference
    between the across-seed means of the two halves exceeds the larger of the
    two halves' across-seed ranges.
    """
    h1, h2 = [], []
    for f in freq_by_seed:
        d = between_site_dispersion(f); m = len(d) // 2
        h1.append(float(d[:m].mean())); h2.append(float(d[m:].mean()))
    h1, h2 = np.array(h1), np.array(h2)
    diff = abs(h1.mean() - h2.mean())
    spread = max(float(np.ptp(h1)), float(np.ptp(h2)))
    return {"first_half": h1.tolist(), "second_half": h2.tolist(), "difference_of_means": float(diff),
            "range_across_seeds": spread, "stationary": bool(diff <= spread)}


def cell_time_scales(run, seeds, S0: int = 5000, max_doublings: int = 1, burnin_floor: int = 5000, innovation: float | None = None) -> dict:
    """Time scales and derived probe settings for one generator cell (screen use).

    Measures L_d and L_R with `measure_time_scales` and returns the stored
    length max(2 L_d + L_R, 10 L_R) and the burn-in max(burnin_floor, 5 L_R,
    5 / innovation), the last term because a cell's regional time scale is
    about one over its innovation rate. When the measurement is not accepted
    within `max_doublings` doublings the status says so and the settings fall
    back to the largest lag examined, flagged by "accepted": False.
    """
    res = measure_time_scales(run, seeds, S0=S0, max_doublings=max_doublings)
    accepted = res["status"] == "ok"
    if accepted:
        L_d, L_R = res["L_d"], res["L_R"]
    else:
        L_d = L_R = int(res["S"] // 4)
    S = int(max(2 * L_d + L_R, 10 * L_R))
    burnin = int(max(burnin_floor, 5 * L_R, 0 if innovation is None else math.ceil(5.0 / innovation)))
    return {"L_d": int(L_d), "L_R": int(L_R), "S": S, "burnin": burnin, "accepted": accepted, "status": res["status"],
            "S_measured": res["S"], "doublings": res["doublings"]}
