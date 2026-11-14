"""Rim R1 statistics (spec `docs/superpowers/specs/2026-10-19-rim-attributes-r1-design.md`, sections 4, 6, 7).

Three kinds of number, kept apart: the observed count-based moment F
(`moment_f`, one definition in the signatures package); generator bands of
the same statistic (made by the runner, not here); and posterior summaries
of latent site frequencies under the sampling model of spec section 6
(Beta(k + alpha, N - k + alpha) per site, Dirichlet(counts + alpha) per
composition, independent sites, attributes independent of one another and
of the composition as a stated working model; direct conjugate draws, no
MCMC). Geography statistics are Spearman correlations; the pair
association is the Spearman of rim and decoration pair differences after
linear adjustment on ranks for distance, with decoration differences as
total variation. Nothing here is a test statistic.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import betabinom

from mls_emergence.signatures.assortativity import partial_spearman, spearman, total_variation
from mls_emergence.signatures.variance import moment_f

__all__ = [
    "composition_draws",
    "frequency_draws",
    "geography_stats",
    "moment_f",
    "pair_association",
    "posterior_summary",
    "prior_predictive_counts",
    "prior_weight",
    "recovery_program",
]
ALPHA = 0.5


def prior_weight(N) -> np.ndarray:
    """Weight of the Beta(alpha, alpha) prior mean in the posterior mean: 1 / (N + 2 alpha), one pseudo-observation at alpha 1/2."""
    return 1.0 / (np.asarray(N, float) + 1.0)


def _check_counts(k, N):
    k = np.asarray(k, float); N = np.asarray(N, float)
    if k.shape != N.shape or k.ndim != 1 or (k < 0).any() or (N <= 0).any() or (k > N).any() \
            or not np.equal(k, np.rint(k)).all() or not np.equal(N, np.rint(N)).all():
        raise ValueError("k and N must be integer vectors with 0 <= k <= N and N > 0")
    return k, N


def frequency_draws(k, N, rng: np.random.Generator, draws: int, alpha: float = ALPHA) -> np.ndarray:
    """Posterior draws of each site's attribute frequency, shape (draws, sites)."""
    k, N = _check_counts(k, N)
    if not (isinstance(draws, (int, np.integer)) and draws > 0) or alpha <= 0:
        raise ValueError("draws must be a positive integer and alpha positive")
    return rng.beta(k + alpha, N - k + alpha, size=(int(draws), len(k)))


def composition_draws(counts, rng: np.random.Generator, draws: int, alpha: float = ALPHA) -> np.ndarray:
    """Posterior draws of each site's composition, shape (draws, sites, classes)."""
    y = np.asarray(counts, float)
    if y.ndim != 2 or (y < 0).any() or not np.equal(y, np.rint(y)).all() or (y.sum(1) <= 0).any():
        raise ValueError("counts must be a non-negative integer site-by-class table with positive row totals")
    return np.stack([rng.dirichlet(y[i] + alpha, size=int(draws)) for i in range(y.shape[0])], axis=1)


def _pairs(n: int):
    return np.triu_indices(n, 1)


def geography_stats(p, D, lat, lon) -> dict:
    """Spearman of |p_i - p_j| with distance over the pairs; Spearman of p with latitude and with longitude. NaN when undefined."""
    p = np.asarray(p, float); D = np.asarray(D, float); lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    n = len(p)
    if D.shape != (n, n) or lat.shape != (n,) or lon.shape != (n,):
        raise ValueError("p, lat and lon need one value per site and D must be n by n")
    iu = _pairs(n); dp = np.abs(p[:, None] - p[None, :])[iu]
    return {"rho_pair_distance": spearman(dp, D[iu]), "rho_latitude": spearman(p, lat), "rho_longitude": spearman(p, lon)}


def pair_association(p, comp, D) -> dict:
    """Rim pair differences |p_i - p_j| against decoration total-variation differences: Spearman, and Spearman after linear rank adjustment for distance."""
    p = np.asarray(p, float); comp = np.asarray(comp, float); D = np.asarray(D, float); n = len(p)
    if comp.shape[0] != n or D.shape != (n, n):
        raise ValueError("comp needs one composition row per site and D must be n by n")
    iu = _pairs(n); rim = np.abs(p[:, None] - p[None, :])[iu]; dec = total_variation(comp)[iu]; dist = D[iu]
    return {"spearman": spearman(rim, dec), "partial_spearman": partial_spearman(rim, dec, dist),
            "rho_rim_distance": spearman(rim, dist), "rho_decoration_distance": spearman(dec, dist)}


def posterior_summary(fn, k, N, draws: int, seeds=(1, 2), alpha: float = ALPHA, counts=None, comp_alpha: float = ALPHA) -> dict:
    """Observed value and posterior median with 90 percent interval of fn(p) (or fn(p, comp) when counts are given), with the
    Monte Carlo error of each quantile from a second seed. Non-finite values of fn are dropped and counted."""
    k, N = _check_counts(k, N)
    observed = fn(k / N) if counts is None else fn(k / N, np.asarray(counts, float) / np.asarray(counts, float).sum(1, keepdims=True))
    runs = []
    for s in seeds:
        rng = np.random.default_rng(s)
        P = frequency_draws(k, N, rng, draws, alpha)
        C = composition_draws(counts, rng, draws, comp_alpha) if counts is not None else None
        vals = np.array([fn(P[d]) if C is None else fn(P[d], C[d]) for d in range(int(draws))], float)
        fin = vals[np.isfinite(vals)]
        runs.append({"n_finite": int(len(fin)), "median": float(np.median(fin)) if len(fin) else np.nan,
                     "q05": float(np.quantile(fin, .05)) if len(fin) else np.nan, "q95": float(np.quantile(fin, .95)) if len(fin) else np.nan})
    a, b = runs[0], runs[1]
    return {"observed": float(observed), "median": a["median"], "q05": a["q05"], "q95": a["q95"], "n_finite": a["n_finite"], "draws": int(draws),
            "mc_error_median": abs(a["median"] - b["median"]), "mc_error_q05": abs(a["q05"] - b["q05"]), "mc_error_q95": abs(a["q95"] - b["q95"]),
            "alpha": alpha, "seeds": tuple(seeds)}


def prior_predictive_counts(N, alpha: float = ALPHA, qs=(0.05, 0.5, 0.95)) -> np.ndarray:
    """Beta-Binomial(N, alpha, alpha) count quantiles at each site's actual N, shape (sites, len(qs)): the observable prior predictive."""
    N = np.asarray(N, int)
    return np.column_stack([betabinom(N, alpha, alpha).ppf(q) for q in qs])


def recovery_program(D, lat, lon, N, replicates: int = 20, draws: int = 1000, seed: int = 0, alpha: float = ALPHA) -> dict:
    """Spec section 7: at three true configurations on the given sites and totals, draw counts, and check that the posterior
    90 percent interval of rho_latitude and of rho_pair_distance covers the truth and that the sign of rho_latitude is recovered.

    positive: frequency rising with latitude; zero: a frequency pattern unrelated to latitude (truth near zero);
    opposite: frequency falling with latitude. Zero counts occur at small sites by construction of the rare end.
    """
    lat = np.asarray(lat, float); lon = np.asarray(lon, float); N = np.asarray(N, int); D = np.asarray(D, float)
    rng = np.random.default_rng(seed); n = len(N)
    # True frequencies carry a latitude component and a fixed non-latitude component, so the true Spearman with latitude is
    # strong but not 1 (a truth at the boundary of the statistic's range cannot be covered by any interval).
    z = (lat - lat.min()) / np.ptp(lat); offset = np.random.default_rng(seed + 999).normal(0, 0.12, n)
    configs = {"positive": np.clip(0.08 + 0.5 * z + offset, 0.02, 0.95), "opposite": np.clip(0.58 - 0.5 * z + offset, 0.02, 0.95),
               "zero": 0.05 + 0.5 * (np.argsort(np.argsort(rng.permutation(n))) / (n - 1))}
    out = {}
    for name, p_true in configs.items():
        truth = geography_stats(p_true, D, lat, lon); cov_lat = cov_pair = sign = 0
        for r in range(replicates):
            k = rng.binomial(N, p_true)
            s_lat = posterior_summary(lambda p: geography_stats(p, D, lat, lon)["rho_latitude"], k, N, draws, seeds=(seed + 100 * r + 1, seed + 100 * r + 2), alpha=alpha)
            s_pair = posterior_summary(lambda p: geography_stats(p, D, lat, lon)["rho_pair_distance"], k, N, draws, seeds=(seed + 100 * r + 3, seed + 100 * r + 4), alpha=alpha)
            cov_lat += int(s_lat["q05"] <= truth["rho_latitude"] <= s_lat["q95"])
            cov_pair += int(s_pair["q05"] <= truth["rho_pair_distance"] <= s_pair["q95"])
            sign += int(np.sign(s_lat["median"]) == np.sign(truth["rho_latitude"]))
        out[name] = {"truth_rho_latitude": truth["rho_latitude"], "truth_rho_pair_distance": truth["rho_pair_distance"],
                     "coverage": cov_lat / replicates, "coverage_pair_distance": cov_pair / replicates, "sign_recovered": sign / replicates}
    return out
