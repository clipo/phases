"""Model-agnostic ABC-SMC (Toni et al. 2009) with weighted-population helpers.

The sampler targets an approximate Bayesian posterior for a simulator with no
tractable likelihood. It replaces plain rejection ABC by propagating a weighted
particle population through a sequence of shrinking tolerances, perturbing
survivors with a Gaussian kernel whose covariance is twice the weighted
population covariance (Filippi et al. 2013 optimal-kernel rule). Regression
adjustment of the final population lives in ``regression.py``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_triangular


@dataclass
class SMCResult:
    thetas: np.ndarray       # (n_particles, d) final weighted population
    weights: np.ndarray      # (n_particles,) normalized importance weights
    summaries: np.ndarray    # (n_particles, s_dim) accepted summaries
    distances: np.ndarray    # (n_particles,) accepted distances
    eps_schedule: list       # tolerance applied per round (non-increasing)
    accept_rates: list       # acceptance rate per round
    n_sims: list             # simulator calls per round


def _weighted_cov(thetas: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Bias-corrected weighted covariance of the particle population."""
    w = weights / weights.sum()
    mean = np.average(thetas, axis=0, weights=w)
    d = thetas - mean
    cov = (w[:, None] * d).T @ d
    denom = 1.0 - np.sum(w ** 2)
    if denom <= 0:
        denom = 1.0
    cov = cov / denom
    # jitter to keep the kernel non-singular if the population collapses
    cov = cov + 1e-12 * np.eye(cov.shape[0])
    return cov


def _mvn_pdf_rows(diff: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """N(0, cov) density for each row of ``diff`` (shape (m, d))."""
    L = np.linalg.cholesky(cov)
    y = solve_triangular(L, diff.T, lower=True)      # (d, m)
    quad = np.sum(y ** 2, axis=0)                    # (m,)
    logdet = 2.0 * np.sum(np.log(np.diag(L)))
    d = cov.shape[0]
    logpdf = -0.5 * (d * np.log(2.0 * np.pi) + logdet + quad)
    return np.exp(logpdf)


def abc_smc(prior_sampler, prior_logpdf, simulator, summary, distance, s_obs,
            n_particles: int = 500, n_rounds: int = 6, quantile: float = 0.4,
            seed: int = 0, accept_rate_floor: float = 0.02,
            max_sims_per_round: int = 500_000) -> SMCResult:
    """Sequential Monte-Carlo ABC. See module docstring."""
    rng = np.random.default_rng(seed)
    s_obs = np.asarray(s_obs, float)
    if n_rounds < 1:
        raise ValueError(f"n_rounds must be >= 1, got {n_rounds}")
    d = np.asarray(prior_sampler(rng)).size

    thetas = np.empty((n_particles, d))
    summ = np.empty((n_particles, s_obs.size))
    dist = np.empty(n_particles)

    eps_schedule: list = []
    accept_rates: list = []
    n_sims_list: list = []
    prev_thetas = None
    prev_weights = None
    kernel_cov = None
    eps = np.inf

    for t in range(n_rounds):
        i = 0
        proposals = 0
        sim_calls = 0
        while i < n_particles:
            proposals += 1
            if proposals > max_sims_per_round:
                raise RuntimeError(
                    f"round {t}: exceeded max_sims_per_round={max_sims_per_round} "
                    f"with {i}/{n_particles} particles at eps={eps:.4g}")
            if t == 0:
                theta = np.asarray(prior_sampler(rng), float)
            else:
                j = rng.choice(n_particles, p=prev_weights)
                theta = rng.multivariate_normal(prev_thetas[j], kernel_cov)
                if not np.isfinite(prior_logpdf(theta)):
                    continue
            sim_calls += 1
            s = np.asarray(summary(simulator(theta, rng)), float)
            dd = float(distance(s, s_obs))
            if dd <= eps:
                thetas[i] = theta
                summ[i] = s
                dist[i] = dd
                i += 1

        if t == 0:
            weights = np.full(n_particles, 1.0 / n_particles)
        else:
            logp = np.array([prior_logpdf(th) for th in thetas])
            weights = np.empty(n_particles)
            for k in range(n_particles):
                kern = _mvn_pdf_rows(thetas[k] - prev_thetas, kernel_cov)
                weights[k] = np.exp(logp[k]) / np.sum(prev_weights * kern)
            wsum = weights.sum()
            if wsum > 0:
                weights /= wsum
            else:
                weights = np.full(n_particles, 1.0 / n_particles)

        eps_schedule.append(float(eps) if np.isfinite(eps) else float(dist.max()))
        accept_rates.append(n_particles / sim_calls)
        n_sims_list.append(sim_calls)

        # next tolerance: monotone non-increasing quantile of accepted distances
        eps = min(eps, float(np.quantile(dist, quantile)))
        prev_thetas = thetas.copy()
        prev_weights = weights.copy()
        kernel_cov = 2.0 * _weighted_cov(prev_thetas, prev_weights)

        if accept_rates[-1] < accept_rate_floor:
            break

    return SMCResult(prev_thetas, prev_weights, summ.copy(), dist.copy(),
                     eps_schedule, accept_rates, n_sims_list)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    return float(np.sum(values * weights) / np.sum(weights))


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w) - 0.5 * w
    cw /= np.sum(w)
    return float(np.interp(q, cw, v))


def resample(values: np.ndarray, weights: np.ndarray, size: int,
             rng: np.random.Generator) -> np.ndarray:
    values = np.asarray(values)
    p = np.asarray(weights, float)
    p = p / p.sum()
    idx = rng.choice(len(values), size=size, p=p)
    return values[idx]
