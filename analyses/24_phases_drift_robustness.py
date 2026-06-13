"""24_phases_drift_robustness.py — robustness of the phases-vs-spatial-drift test.

Hardens analyses/23 in two ways the first pass flagged:

  1. PARAMETER SWEEP. The Lipo-2021 spatial-drift null depends on the interaction
     length scale, the between-node copy rate, and the innovation rate. We sweep a
     grid and report, for every cell, whether the four observed statistics fall
     inside the spatial-drift null. A finding that survives the grid is robust to
     the (unobservable) interaction parameters; one that flips is not.

  2. TIME-AVERAGED NULL. Real PFG assemblages aggregate sherds over a span of
     time; the snapshot null of analyses/23 does not. Time-averaging blurs drift
     fluctuations and homogenizes assemblages, so the snapshot null may overstate
     differentiation. We add a time-averaged null (deposition pooled over a
     trailing window) and ask whether it reconciles the observed (low) modularity
     and F_ST with spatial drift.

The simulation is vectorized across nodes so the grid is tractable.

Read-only on the manuscript. Writes output/phases_drift_robustness.md.

Usage: .venv/bin/python analyses/24_phases_drift_robustness.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import make_figures as mf  # noqa: E402
sd = importlib.import_module("23_phases_vs_spatial_drift")  # geo_km, stats_for_matrix, ...

OUT = ROOT / "output" / "phases_drift_robustness.md"


def simulate_vectorized(coords, N=100, K=10, mu=0.01, length_km=12.0,
                        m_between=0.05, steps=1500, avg_window=0,
                        total_per_node=None, seed=0):
    """Vectorized Wright-Fisher neutral drift on a distance-decayed network.

    Identical model to analyses/23.simulate_spatial_drift but stepped over all
    nodes at once. If avg_window > 0, the assemblage is the type tally pooled
    over the final `avg_window` steps (time-averaging); otherwise it is the
    final-snapshot population. Output rescaled to observed sherd totals.
    """
    rng = np.random.default_rng(seed)
    n = coords.shape[0]
    D = sd.geo_km(coords)
    W = np.exp(-D / length_km)
    np.fill_diagonal(W, 0.0)
    W = W / W.sum(axis=1, keepdims=True)
    cumW = np.cumsum(W, axis=1)                       # (n, n) per-row CDF

    pop = rng.integers(0, K, size=(n, N))
    next_type = K
    rows = np.arange(n)[:, None]

    # Time-averaging: pool DEPOSITION over a trailing window. Record the
    # population at thinned steps within the window, then tally the K most
    # frequent type labels ACROSS THE POOLED WINDOW. Tallying the pooled window
    # (rather than the original label set, which innovation steadily replaces) is
    # the same top-K tabulation the snapshot path uses, applied over time.
    start_avg = steps - avg_window
    thin = max(1, avg_window // 100) if avg_window > 0 else 1
    deposits = []

    for t in range(steps):
        from_other = rng.random((n, N)) < m_between
        src_self = pop[rows, rng.integers(0, N, size=(n, N))]
        u = rng.random((n, N))
        partners = (u[:, :, None] > cumW[:, None, :]).sum(axis=2)
        partners = np.clip(partners, 0, n - 1)
        src_other = pop[partners, rng.integers(0, N, size=(n, N))]
        drawn = np.where(from_other, src_other, src_self)
        innov = rng.random((n, N)) < mu
        n_innov = int(innov.sum())
        if n_innov:
            drawn = drawn.copy()
            drawn[innov] = np.arange(next_type, next_type + n_innov)
            next_type += n_innov
        pop = drawn
        if avg_window > 0 and t >= start_avg and (steps - 1 - t) % thin == 0:
            deposits.append(pop.copy())

    record = np.concatenate(deposits, axis=1) if avg_window > 0 else pop
    M = np.zeros((n, K))
    flat = record.ravel()
    labs, fr = np.unique(flat, return_counts=True)
    keep = labs[np.argsort(-fr)[:K]]
    ki = {lab: c for c, lab in enumerate(keep)}
    for i in range(n):
        u_, c_ = np.unique(record[i], return_counts=True)
        for lab, c in zip(u_, c_):
            if lab in ki:
                M[i, ki[lab]] += c

    if total_per_node is not None:
        out = np.zeros((n, K), dtype=int)
        for i in range(n):
            tot = int(total_per_node[i])
            p = M[i] / M[i].sum() if M[i].sum() else np.ones(K) / K
            out[i] = rng.multinomial(tot, p) if tot > 0 else 0
        return out
    return M.astype(int)


def null_cell(coords, geo_dist, totals, obs, K, *, length_km, m_between, mu,
              avg_window, n_seeds=100, base_seed=0):
    """Run n_seeds sims for one parameter cell; return inside-flags + summaries."""
    keys = ["dd_r", "Q", "be", "fst"]
    vals = {k: [] for k in keys}
    for s in range(n_seeds):
        M = simulate_vectorized(coords, N=100, K=K, mu=mu, length_km=length_km,
                                m_between=m_between, steps=1200,
                                avg_window=avg_window, total_per_node=totals,
                                seed=base_seed + s)
        if M.sum() == 0:
            continue
        st = sd.stats_for_matrix(M, geo_dist, seed=s)
        for k in keys:
            if np.isfinite(st[k]):
                vals[k].append(st[k])
    res = {}
    for k in keys:
        v = np.array(vals[k], float)
        # frac of seeds that even admitted the statistic (e.g. F_ST is undefined
        # when time-averaging collapses the synthetic basin to one community)
        defined = len(v) / max(1, n_seeds)
        if len(v) < 3:
            res[k] = dict(mean=float("nan"), lo=float("nan"), hi=float("nan"),
                          inside=None, defined=defined)
            continue
        lo, hi = np.percentile(v, [2.5, 97.5])
        res[k] = dict(mean=float(v.mean()), lo=float(lo), hi=float(hi),
                      inside=bool(lo <= obs[k] <= hi), defined=defined)
    return res


def main():
    counts_df, coords_df = mf._load_curated()
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    K = counts.shape[1]
    geo_dist = sd.geo_km(coords)
    totals = counts.sum(1)
    obs = sd.stats_for_matrix(counts, geo_dist, seed=0)

    L = ["# Robustness of the phases-vs-spatial-drift test",
         "",
         f"Observed (n = {counts.shape[0]} curated assemblages): "
         f"distance-decay r = {obs['dd_r']:+.3f}, modularity Q = {obs['Q']:.3f}, "
         f"boundary excess = {obs['be']:+.1f} BR, cultural F_ST = {obs['fst']:+.3f}.",
         "",
         "Each cell = 100 seeds. 'inside' means the observed statistic falls within",
         "the spatial-drift 95% null. A result robust to interaction parameters is",
         "one whose inside/outside verdict does not flip across the grid.", ""]

    lengths = [6.0, 12.0, 24.0]
    m_betweens = [0.02, 0.05, 0.20]
    mus = [0.005, 0.01, 0.02]

    # ---- snapshot null sweep --------------------------------------------- #
    ncells = 2 * len(lengths) * len(m_betweens) * len(mus)
    done = 0
    for avg_window, tag in [(0, "Snapshot null"), (300, "Time-averaged null (window 300)")]:
        print(f"=== {tag} ===", flush=True)
        L += [f"## {tag}", "",
              "| length_km | m_between | mu | dd_r inside | Q inside | bnd-exc inside | F_ST inside |",
              "|---|---|---|---|---|---|---|"]
        for lk in lengths:
            for mb in m_betweens:
                for mu in mus:
                    r = null_cell(coords, geo_dist, totals, obs, K,
                                  length_km=lk, m_between=mb, mu=mu,
                                  avg_window=avg_window, n_seeds=100,
                                  base_seed=10_000)
                    done += 1
                    print(f"  [{done}/{ncells}] {tag[:9]} len={lk:.0f} "
                          f"m={mb:.2f} mu={mu:.3f}", flush=True)
                    def cell(k):
                        d = r[k]
                        if d["inside"] is None:
                            return f"n/a (def {d['defined']:.0%})"
                        mark = "yes" if d["inside"] else "NO"
                        return f"{mark} ({d['lo']:+.2f},{d['hi']:+.2f})"
                    L += [f"| {lk:.0f} | {mb:.2f} | {mu:.3f} | "
                          f"{cell('dd_r')} | {cell('Q')} | {cell('be')} | {cell('fst')} |"]
        L += [""]

    # ---- summary tallies -------------------------------------------------- #
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
