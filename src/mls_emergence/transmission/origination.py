"""Generator C4: neutral copying among sites with a finite regional pool in which new types originate.

Design: docs/superpowers/specs/2026-11-14-type-origination-generator-design.md
(draft 1; exploratory, not audited). The site dynamics are those of
`spatial.drift_record`. The regional innovation target is the composition of
an explicit pool of `pool_size` draws. After burn-in the pool is resampled at
every step as a Wright-Fisher population with infinite-alleles innovation:
each of its draws is a new type with probability mu = theta / (2 pool_size),
and otherwise copies a draw of the current pool at random. theta is the
parameter of the Ewens sampling formula. With theta = 0 the pool is the
exact-loss counterpart of generator C3's Dirichlet(kappa * target) redraw
with kappa = pool_size - 1: the per-step variance of a share p is
p (1 - p) / pool_size in both, and the F accumulated over W steps is
1 - (1 - 1 / pool_size)^W.

The set of types is open, so the full record is not stored. The caller names
the accumulation spans (end step and length per site, per observation
design) before the run, and the simulator returns each design's time-averaged
site profiles over the types that have any mass in them.
"""
from __future__ import annotations

import numpy as np

INITIAL = -1          # origin step recorded for the classes present at the start


def largest_remainder(shares, total: int) -> np.ndarray:
    """Integer counts summing to `total`, proportional to `shares`: floors, then the largest remainders (ties to the lower index)."""
    s = np.asarray(shares, float)
    if s.ndim != 1 or (s < 0).any() or not np.isfinite(s).all() or s.sum() <= 0:
        raise ValueError("shares must be a finite non-negative vector with a positive sum")
    if isinstance(total, (bool, np.bool_)) or not isinstance(total, (int, np.integer)) or total < 1:
        raise ValueError("total must be a positive integer")
    x = s / s.sum() * int(total); c = np.floor(x).astype(np.int64)
    short = int(total) - int(c.sum())
    if short:
        c[np.argsort(-(x - c), kind="stable")[:short]] += 1
    return c


def pool_step(counts: np.ndarray, mu: float, rng: np.random.Generator):
    """One resampling of the pool: (copied counts over the existing types, number of new types, each entering with one draw).

    The number of new draws is Binomial(N, mu). The other N - new draws are
    Multinomial over the current shares, so a type at count zero stays at
    zero. mu = 0 draws nothing for origination.
    """
    c = np.asarray(counts, np.int64); N = int(c.sum())
    if N < 2 or (c < 0).any():
        raise ValueError("the pool needs non-negative counts summing to at least 2")
    if not np.isfinite(mu) or mu < 0 or mu > 1:
        raise ValueError("mu must lie in [0, 1]")
    n_new = int(rng.binomial(N, mu)) if mu > 0 else 0
    return rng.multinomial(N - n_new, c / N), n_new


def origination_profiles(weights, *, n_ind, mixing, innovation, steps, burnin, seed, pool_seed, target, pool_size, theta,
                         spans, max_slots: int = 4096, trace: bool = False) -> dict:
    """Run C4 and return the time-averaged site profiles of each named accumulation design.

    `spans` maps a design name to (ends, lengths), integer arrays with one
    entry per site: the site's deposit averages its type shares over the
    `length` stored steps ending at step `end` (0-based, as in
    `spans.span_profile`). `target` gives the initial shares of the starting
    classes, allocated to the pool by largest remainder; the pool is fixed
    during burn-in and resampled after the learner draw of every stored step.
    Returns {"profiles": {name: (n_sites, n_types) rows summing to 1},
    "origin": {name: origin step per column, -1 for a starting class},
    "info": {...}} and, with trace=True, the per-step (type ids, shares)
    used by the brute-force test.
    """
    w = np.asarray(weights, float); n = w.shape[0]
    if w.ndim != 2 or w.shape[1] != n or (w < 0).any() or not np.allclose(w.sum(1), 1):
        raise ValueError("weights must be a row-stochastic square matrix")
    if isinstance(pool_size, (bool, np.bool_)) or not isinstance(pool_size, (int, np.integer)) or pool_size < 2:
        raise ValueError("pool_size must be an integer of at least 2")
    if not np.isfinite(theta) or theta < 0:
        raise ValueError("theta must be finite and non-negative")
    mu = float(theta) / (2.0 * int(pool_size))
    if mu > 1:
        raise ValueError("theta exceeds twice the pool size")
    if not (0 <= mixing <= 1 and 0 <= innovation <= 1) or min(int(n_ind), steps) < 1 or burnin < 0:
        raise ValueError("invalid rates, population or time configuration")
    c0 = largest_remainder(target, int(pool_size)); k0 = len(c0)
    if max_slots < k0 + 1:
        raise ValueError("max_slots must exceed the number of starting classes")
    designs = {}
    for name, (ends, lengths) in spans.items():
        e, ln = np.asarray(ends, np.int64), np.asarray(lengths, np.int64)
        if e.shape != (n,) or ln.shape != (n,) or (ln < 1).any() or (e - ln + 1 < 0).any() or (e >= steps).any():
            raise ValueError(f"design {name!r}: spans must lie inside the stored steps, one per site")
        designs[name] = (e - ln + 1, e, ln)
    rng = np.random.default_rng(seed); prng = np.random.default_rng(pool_seed)
    counts = np.zeros(max_slots, np.int64); counts[:k0] = c0
    origin = np.full(max_slots, INITIAL, np.int64); type_id = np.full(max_slots, -1, np.int64); type_id[:k0] = np.arange(k0)
    used = np.zeros(max_slots, bool); used[:k0] = True          # a slot is used while its type is live in the pool or at a site
    next_id = k0
    p = np.zeros((n, max_slots)); p[:, :k0] = c0 / c0.sum()
    acc = {name: np.zeros((n, max_slots)) for name in designs}
    done = {name: [] for name in designs}                       # (origin step, column) of freed types that hold mass
    steps_trace = []; n_orig = 0; max_live = k0; n_freed = 0
    n_pop = int(n_ind); N = int(pool_size)
    for t in range(-burnin, steps):
        idx = np.flatnonzero(used)
        pa = p[:, idx]; tgt = counts[idx] / N
        q = (1 - mixing) * pa + mixing * (w @ pa)
        q = (1 - innovation) * q + innovation * tgt
        q = np.maximum(q, 0); q /= q.sum(1, keepdims=True)
        pa = rng.multinomial(n_pop, q) / n_pop
        p[:, idx] = pa
        if t < 0:
            continue
        for name, (start, end, _ln) in designs.items():
            on = (start <= t) & (t <= end)
            if on.any():
                acc[name][np.ix_(on, idx)] += pa[on]
        if trace:
            steps_trace.append((type_id[idx].copy(), pa.copy()))
        copied, n_new = pool_step(counts[idx], mu, prng)
        counts[idx] = copied
        dead = idx[(copied == 0) & (pa.sum(0) == 0)]            # gone from the pool and from every site: the type cannot return
        if len(dead):
            for name in designs:
                for s in dead[acc[name][:, dead].any(0)]:
                    done[name].append((int(origin[s]), int(type_id[s]), acc[name][:, s].copy()))
                acc[name][:, dead] = 0.0
            used[dead] = False; p[:, dead] = 0.0; type_id[dead] = -1; n_freed += len(dead)
        if n_new:
            free = np.flatnonzero(~used)[:n_new]
            if len(free) < n_new:
                raise RuntimeError(f"no free slot for a new type at step {t} (max_slots {max_slots})")
            counts[free] = 1; origin[free] = t; type_id[free] = np.arange(next_id, next_id + n_new)
            used[free] = True; next_id += n_new; n_orig += n_new
        max_live = max(max_live, int(used.sum()))
    out = {"profiles": {}, "origin": {}, "type_id": {},
           "info": {"pool_size": N, "mu": mu, "initial_present": int((c0 > 0).sum()), "originated": n_orig,
                    "pool_types_end": int((counts > 0).sum()), "max_live_slots": max_live}}
    for name, (_start, _end, ln) in designs.items():
        live = np.flatnonzero(used & acc[name].any(0))
        cols = [acc[name][:, s] for s in live] + [c for _o, _i, c in done[name]]
        org = [int(origin[s]) for s in live] + [o for o, _i, _c in done[name]]
        ids = [int(type_id[s]) for s in live] + [i for _o, i, _c in done[name]]
        order = np.argsort(ids, kind="stable")
        prof = np.column_stack(cols)[:, order] / ln[:, None]
        if not np.allclose(prof.sum(1), 1.0, atol=1e-9):
            raise AssertionError("a site's accumulated profile does not sum to one")
        out["profiles"][name] = prof / prof.sum(1, keepdims=True)
        out["origin"][name] = np.asarray(org)[order]; out["type_id"][name] = np.asarray(ids)[order]
    if trace:
        out["trace"] = steps_trace
    return out
