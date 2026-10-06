"""The copying model: neutral copying between assemblages, falling off with distance, with no groups.

Each assemblage is a population of learners. At every step a learner keeps the
class of someone in its own population, copies someone in another population
with probability `mixing` (nearer populations more often), or adopts a new
variant with probability `innovation`, drawn in proportion to each class's
share of the whole record. Populations are finite, so class frequencies drift.
The simulated record is then pooled over a window of steps, as a deposit
accumulates, and sampled at the real sherd counts.

A restriction on copying across phase lines is the same model with copying
between assemblages of different phases multiplied by `leak` (1 is no
restriction; 0.03 cuts it to 3 percent).
"""
from __future__ import annotations

import numpy as np

STEPS, RECORDS, BURNIN, WINDOW = 1800, 90, 1200, 8


def copying_weights(dist: np.ndarray, length: float, labels=None, leak: float = 1.0) -> np.ndarray:
    """Who copies whom: weights that fall off exponentially with distance, each row summing to 1."""
    d = np.asarray(dist, float)
    if d.ndim != 2 or d.shape[0] != d.shape[1] or len(d) < 2:
        raise ValueError("distances must be a square matrix with at least two assemblages")
    if not np.isfinite(d).all() or (d < 0).any():
        raise ValueError("distances must be finite and non-negative")
    if not length > 0:
        raise ValueError("the copying length must be positive")
    if not 0 <= leak <= 1:
        raise ValueError("leak must lie between 0 (no copying across a line) and 1 (no restriction)")
    w = np.exp(-d / length)
    if labels is not None:
        labels = np.asarray(labels)
        if labels.shape != (len(d),):
            raise ValueError("one label is required per assemblage")
        w *= np.where(labels[:, None] == labels[None, :], 1.0, leak)
    elif leak != 1.0:
        raise ValueError("a leak below 1 needs labels saying which lines restrict copying")
    np.fill_diagonal(w, 0)
    sums = w.sum(1)
    alone = np.where(sums == 0)[0]            # an isolated assemblage copies itself
    w[alone, alone] = 1
    return w / w.sum(1, keepdims=True)


def drift_record(weights: np.ndarray, *, k: int, n_ind: int, mixing: float, innovation: float,
                 target: np.ndarray, seed: int, steps: int = STEPS, n_records: int = RECORDS,
                 burnin: int = BURNIN) -> np.ndarray:
    """Class frequencies through time: an array (recorded step, assemblage, class)."""
    w = np.asarray(weights, float)
    if w.ndim != 2 or w.shape[0] != w.shape[1] or (w < 0).any() or not np.allclose(w.sum(1), 1):
        raise ValueError("weights must be a square matrix whose rows sum to 1")
    if not (0 <= mixing <= 1 and 0 <= innovation <= 1):
        raise ValueError("mixing and innovation must be probabilities")
    if min(k, int(n_ind), steps, n_records) < 1 or n_records > steps or burnin < 0:
        raise ValueError("invalid population or time settings")
    target = np.asarray(target, float)
    if target.shape != (k,) or (target < 0).any() or target.sum() <= 0:
        raise ValueError("the innovation profile must be a non-negative vector with one entry per class")
    target = target / target.sum()
    rng = np.random.default_rng(seed)
    p = np.tile(target, (len(w), 1))
    record_at = set(np.linspace(0, steps - 1, n_records).astype(int))
    records = []
    for t in range(-burnin, steps):
        q = (1 - mixing) * p + mixing * (w @ p)
        q = (1 - innovation) * q + innovation * target
        q = np.maximum(q, 0)
        q /= q.sum(1, keepdims=True)
        draw = rng.multinomial(int(n_ind), q)
        p = draw / draw.sum(1, keepdims=True)
        if t in record_at:
            records.append(p.copy())
    return np.stack(records)


def sample_record(record: np.ndarray, order: np.ndarray, totals: np.ndarray, rng: np.random.Generator,
                  window: int = WINDOW) -> np.ndarray:
    """Sherd counts drawn from the simulated frequencies, pooled over a trailing window of steps.

    `order` places each assemblage along the recorded sequence, 0 (earliest) to
    1 (latest). All ones treats the assemblages as contemporaneous.
    """
    f = np.asarray(record, float)
    order, totals = np.asarray(order, float), np.asarray(totals)
    if order.shape != (f.shape[1],) or totals.shape != order.shape:
        raise ValueError("one order value and one sherd total are required per assemblage")
    if not np.isfinite(order).all() or (order < 0).any() or (order > 1).any():
        raise ValueError("order values must lie between 0 and 1")
    if not 1 <= window <= len(f):
        raise ValueError("the window must fit within the recorded sequence")
    slots = window - 1 + np.rint(order * (len(f) - window)).astype(int)
    p = np.array([f[s - window + 1:s + 1, i].mean(0) for i, s in enumerate(slots)])
    p /= p.sum(1, keepdims=True)
    return np.array([rng.multinomial(int(n), q) for n, q in zip(totals, p)])
