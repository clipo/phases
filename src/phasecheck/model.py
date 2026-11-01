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

This is the engine behind questions 6 to 9; `copying.py` calibrates and runs
it. It is the paper's own model (`mls_emergence.transmission.spatial`),
reimplemented here so the tool stands alone; `tests/phasecheck/
test_copying_model.py::test_the_model_is_the_papers_model` checks that the
two give identical weights, frequencies and sampled counts at the same seed.
Distances and the copying length are in km. Each simulation takes an integer
seed for its drift; the sampling of sherds takes its own Generator, so the
caller decides how the two seeds relate (see `copying.simulate`).
"""
from __future__ import annotations

import numpy as np

# The paper's implementation defaults: 1,200 burn-in steps from the pooled
# profile, then 1,800 steps of which 90 evenly spaced ones are recorded (about
# every 20th step), and each deposit pools 8 consecutive recorded steps.
STEPS, RECORDS, BURNIN, WINDOW = 1800, 90, 1200, 8


def copying_weights(dist: np.ndarray, length: float, labels=None, leak: float = 1.0) -> np.ndarray:
    """Who copies whom: weights that fall off exponentially with distance, each row summing to 1.

    Row i gives the share of assemblage i's outside copying that comes from
    each other assemblage, proportional to exp(-distance / length), times
    `leak` for a pair in different groups. The diagonal is zero, since
    copying within an assemblage is the `1 - mixing` share in `drift_record`;
    an assemblage whose row would be all zero (every other weight cut to 0
    by `leak` or underflow) copies itself instead.

    Parameters
    ----------
    dist : array_like, shape (n, n)
        Distances in km, finite and non-negative, n of at least 2.
    length : float
        Distance (km) over which a weight falls by a factor e; positive.
    labels : array_like, shape (n,), optional
        Group of each assemblage; needed only when `leak` is below 1.
    leak : float
        0 to 1. Multiplier on copying across a group line; 1 is no
        restriction.

    Returns
    -------
    numpy.ndarray, shape (n, n)
        Row-stochastic weights.

    Raises
    ------
    ValueError
        On a non-square, non-finite or negative distance matrix, a length not
        above 0, a leak outside 0..1, labels of the wrong shape, or a leak
        below 1 without labels.
    """
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
    """Class frequencies through time: an array (recorded step, assemblage, class).

    Every assemblage starts at `target`. At each step its expected profile is
    mixed with its neighbors' (share `mixing`, by `weights`), then with the
    innovation profile (share `innovation`), and `n_ind` learners are drawn
    from it multinomially; that drift is the only source of randomness.

    Parameters
    ----------
    weights : array_like, shape (n, n)
        Row-stochastic copying weights, as from `copying_weights`.
    k : int
        Number of classes, at least 1.
    n_ind : int
        Learners per assemblage, at least 1. Smaller means more drift.
    mixing : float
        0 to 1. Share of each step's copying that comes from other
        assemblages.
    innovation : float
        0 to 1. Share of learners who adopt a new variant, drawn from
        `target`.
    target : array_like, shape (k,)
        Innovation profile, non-negative with a positive sum; normalized
        here. The caller passes the record's pooled class proportions.
    seed : int
        Seed of the drift.
    steps, n_records, burnin : int
        Steps after burn-in (at least 1), steps recorded among them (1 to
        `steps`), and unrecorded burn-in steps (at least 0).

    Returns
    -------
    numpy.ndarray, shape (n_records, n, k)
        Class proportions of each assemblage at each recorded step.

    Raises
    ------
    ValueError
        On weights that are not square, non-negative and row-stochastic,
        rates outside 0..1, invalid population or time settings, or a
        `target` of the wrong shape or with no positive entry.
    """
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
    record_at = set(np.linspace(0, steps - 1, n_records).astype(int))   # evenly spaced, both ends kept
    records = []
    for t in range(-burnin, steps):                 # negative t is burn-in, never recorded
        q = (1 - mixing) * p + mixing * (w @ p)
        q = (1 - innovation) * q + innovation * target
        q = np.maximum(q, 0)                         # clear rounding below zero before the draw
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

    Each assemblage's profile is the mean of `window` consecutive recorded
    steps ending at its position, so the earliest assemblage still gets a
    full window rather than a shorter one. Its counts are then a multinomial
    draw at its real sherd total.

    Parameters
    ----------
    record : array_like, shape (n_records, n, k)
        Output of `drift_record`.
    order : array_like, shape (n,)
        Positions 0 to 1.
    totals : array_like, shape (n,)
        Sherd total of each assemblage.
    rng : numpy.random.Generator
        Source of the sampling.
    window : int
        Recorded steps pooled into one deposit, 1 to n_records.

    Returns
    -------
    numpy.ndarray of int, shape (n, k)
        Simulated counts with the real row totals.

    Raises
    ------
    ValueError
        On mismatched shapes, positions outside 0..1, or a window that does
        not fit the record.
    """
    f = np.asarray(record, float)
    order, totals = np.asarray(order, float), np.asarray(totals)
    if order.shape != (f.shape[1],) or totals.shape != order.shape:
        raise ValueError("one order value and one sherd total are required per assemblage")
    if not np.isfinite(order).all() or (order < 0).any() or (order > 1).any():
        raise ValueError("order values must lie between 0 and 1")
    if not 1 <= window <= len(f):
        raise ValueError("the window must fit within the recorded sequence")
    # Last recorded step of each window: order 0 -> window-1, order 1 -> the final step.
    slots = window - 1 + np.rint(order * (len(f) - window)).astype(int)
    p = np.array([f[s - window + 1:s + 1, i].mean(0) for i, s in enumerate(slots)])
    p /= p.sum(1, keepdims=True)
    return np.array([rng.multinomial(int(n), q) for n, q in zip(totals, p)])
