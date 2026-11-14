"""Finite-class Wright–Fisher copying with an explicit observation process.

The frequency-state update is distributionally equivalent to independently
choosing source learners for every copying event. Innovations draw from the
fixed archaeological class repertoire, including possible self draws. The
innovation target is uniform by default; a supplied ``target`` profile
represents uneven aggregation of variants into the measured classes.
"""
import numpy as np


def copying_weights(distances, length=24.0, labels=None, leak=1.0, kernel="exponential"):
    """Row-stochastic copying weights from a distance matrix.

    ``kernel`` sets how influence falls with distance, which is an assumption
    the exponential default makes silently. A Gaussian falls away faster at
    range and so isolates distant nodes more; a power law has a heavy tail and
    isolates them less. Varying the interaction length moves along one family;
    changing the form asks whether the family itself is the constraint.
    """
    d = np.asarray(distances, float)
    if d.ndim != 2 or d.shape[0] != d.shape[1] or len(d) < 2:
        raise ValueError("distances must be a square matrix with at least two nodes")
    if not np.isfinite(d).all() or (d < 0).any() or length <= 0 or not 0 <= leak <= 1:
        raise ValueError("invalid distances, length, or leak")
    if kernel == "exponential":
        w = np.exp(-d / length)
    elif kernel == "gaussian":
        w = np.exp(-(d / length) ** 2)
    elif kernel == "power":
        w = (1.0 + d / length) ** -2.0
    else:
        raise ValueError(
            f"kernel must be exponential, gaussian or power, got {kernel!r}")
    if labels is not None:
        labels = np.asarray(labels)
        if labels.shape != (len(d),):
            raise ValueError("one label is required per node")
        w *= np.where(labels[:, None] == labels[None, :], 1.0, leak)
    np.fill_diagonal(w, 0)
    sums = w.sum(1)
    # A completely isolated singleton copies itself on an external-copy event.
    isolated = sums == 0
    w[np.where(isolated)[0], np.where(isolated)[0]] = 1
    return w / w.sum(1, keepdims=True)


def kappa_for_window_f(F_W, W):
    """Dirichlet concentration that accumulates regional F ``F_W`` over ``W`` target updates (spec 3.3).

    F_W = 1 - (kappa/(kappa+1))^W, inverted in the stable form
    1 / expm1(-log1p(-F_W) / W). F_W = 0 returns None (no update, no randomness
    consumed). F_W must be finite in [0, 1); W a positive integer (not a bool).
    """
    if isinstance(W, bool) or not isinstance(W, (int, np.integer)) or W < 1:
        raise ValueError("W must be a positive integer number of target updates")
    if isinstance(F_W, (bool, str, bytes)) or not isinstance(F_W, (int, float, np.integer, np.floating)):
        raise ValueError("F_W must be a real number")
    F = float(F_W)
    if not np.isfinite(F) or F < 0 or F >= 1:
        raise ValueError("F_W must be finite and lie in [0, 1)")
    if F == 0:
        return None
    kappa = float(1.0 / np.expm1(-np.log1p(-F) / int(W)))
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("F_W is too small to give a finite kappa at this W")
    return kappa


def step_target(target, kappa, rng):
    """One Dirichlet(kappa * target) redraw of the regional target (spec 3.1).

    E[target' | target] = target and Cov = (diag(target) - target target^T) /
    (kappa + 1). A component at exactly zero stays zero (the sampler's
    absorbing face); a positive component can reach an exact zero in float64,
    which the caller counts as numerical target loss.
    """
    t = np.asarray(target, float)
    if t.ndim != 1:
        raise ValueError("step_target takes one 1-D target; per-node targets are not walked")
    if t.shape[0] < 2:
        raise ValueError("a target with one class cannot drift")
    if (t < 0).any() or not np.isfinite(t).all() or t.sum() <= 0:
        raise ValueError("target must be a finite non-negative vector with a positive sum")
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    return rng.dirichlet(float(kappa) * t / t.sum())


def drift_record(weights, *, k=10, n_ind=120, mixing=.02, innovation=.012,
                 steps=1800, n_records=90, burnin=1200, seed=0,
                 initial="uniform", target=None, target_kappa=None, return_targets=False):
    """Return (record time, node, class) frequencies after an explicit burn-in.

    ``target`` is the class distribution an innovation event draws from. None
    gives the uniform K-allele model. A nonnegative vector of length ``k``
    (normalized here) gives class-specific innovation rates ``innovation *
    target``; frequencies then start at the target unless ``initial`` is
    monomorphic.

    ``target_kappa`` (regional-drift spec 3.1): when set, the 1-D target is
    redrawn from Dirichlet(kappa * target) after the learner draw of every
    post-burn-in step, saved or not, by ``step_target``; the learner draw at
    step t uses the snapshot p_t and the update after it makes p_{t+1}. None
    keeps the target fixed and consumes no extra random number, so the output
    is identical to the model without the argument. ``return_targets`` also
    returns the snapshots used at the recorded steps as a second array,
    shape (n_records, k) for a shared target and (n_records, n, k) for a
    per-node target (constant when the target is fixed).
    """
    w = np.asarray(weights, float)
    if w.ndim != 2 or w.shape[0] != w.shape[1] or (w < 0).any() or not np.allclose(w.sum(1), 1):
        raise ValueError("weights must be a row-stochastic square matrix")
    # `mixing` may be one rate for every node, or one per node. The second form
    # is what analyses/73_connectivity_mixing.py needs. Because `w` is
    # row-stochastic, a scalar mixing gives every node the same total external
    # influence however isolated it is: the kernel decides WHO you copy, never
    # HOW MUCH. A per-node vector lets total copying fall with connectivity,
    # which is the assumption the scalar form makes silently.
    mix = np.asarray(mixing, float)
    if mix.ndim not in (0, 1) or (mix.ndim == 1 and mix.shape[0] != len(w)):
        raise ValueError("mixing must be a scalar or one rate per node")
    if not (np.all(mix >= 0) and np.all(mix <= 1)) or not 0 <= innovation <= 1:
        raise ValueError("mixing and innovation must be probabilities")
    mixing = mix if mix.ndim == 0 else mix[:, None]
    # n_ind may be one population for every node, or one per node. The second
    # form is what analyses/64_unequal_populations.py needs: drift rate goes as
    # 1/N, so equal populations everywhere is a substantive assumption, not a
    # convenience, and it is the assumption this model made until now.
    n_arr = np.asarray(n_ind)
    if n_arr.ndim == 0:
        n_pop = int(n_arr)
    elif n_arr.ndim == 1 and n_arr.shape[0] == w.shape[0]:
        n_pop = np.asarray(np.rint(n_arr), dtype=np.int64)
        if (n_pop < 1).any():
            raise ValueError("every per-node population must be at least 1")
    else:
        raise ValueError(
            f"n_ind must be a scalar or one value per node ({w.shape[0]}), "
            f"got shape {n_arr.shape}")
    if min(k, int(np.min(n_pop)), steps, n_records) < 1 or n_records > steps or burnin < 0:
        raise ValueError("invalid population or time configuration")
    # target may be one profile shared by every node, or one profile per node.
    # The per-node form is what a boundary that interrupts innovation looks
    # like: a region-wide pooled profile is an input no boundary can stop, so a
    # model given one understates what a boundary would do, whatever else it
    # does. See analyses/65_other_departures.py.
    if target is None:
        target = np.full(k, 1.0 / k)
    target = np.asarray(target, float)
    if target.ndim == 1:
        if target.shape != (k,) or (target < 0).any() or target.sum() <= 0:
            raise ValueError("target must be a nonnegative vector of length k")
        target = target / target.sum()
    elif target.ndim == 2:
        if target.shape != (w.shape[0], k) or (target < 0).any() or (target.sum(1) <= 0).any():
            raise ValueError(
                f"a per-node target must have shape ({w.shape[0]}, {k}), "
                f"got {target.shape}")
        target = target / target.sum(1, keepdims=True)
    else:
        raise ValueError("target must be 1-D or 2-D")
    if target_kappa is not None:
        if target.ndim != 1:
            raise ValueError("target_kappa walks one shared target; a per-node target cannot be combined with it")
        if k < 2:
            raise ValueError("target_kappa needs at least two classes")
        if isinstance(target_kappa, (bool, np.bool_)) or not np.isfinite(target_kappa) or target_kappa <= 0:
            raise ValueError("target_kappa must be finite and positive, or None")
    rng = np.random.default_rng(seed)
    p = np.tile(target, (len(w), 1)) if target.ndim == 1 else target.copy()
    if initial == "monomorphic":
        p[:] = 0
        p[:, 0] = 1
    elif initial != "uniform":
        raise ValueError("initial must be uniform or monomorphic")
    record_at = set(np.linspace(0, steps - 1, n_records).astype(int))
    records = []; targets = []
    for t in range(-burnin, steps):
        q = (1 - mixing) * p + mixing * (w @ p)
        q = (1 - innovation) * q + innovation * target
        q = np.maximum(q, 0)
        q /= q.sum(1, keepdims=True)
        # multinomial broadcasts n over rows, so a per-node vector gives each
        # node its own sampling intensity and therefore its own drift rate.
        draw = rng.multinomial(n_pop, q)
        p = draw / draw.sum(1, keepdims=True)
        if t in record_at:
            records.append(p.copy())
            if return_targets:
                targets.append(np.array(target, float, copy=True) if target.ndim == 1 else target.copy())
        if target_kappa is not None and t >= 0:
            target = step_target(target, target_kappa, rng)
    if return_targets:
        return np.stack(records), np.stack(targets)
    return np.stack(records)


def sample_record(record, ranks, totals, rng, window=8):
    """Sample observed sherd totals with the same full trailing window per site.

    Earliest rank maps to window-1, rather than to a shorter startup window.
    Ranks of one give a contemporaneous, time-averaged comparison.
    """
    f = np.asarray(record, float)
    ranks, totals = np.asarray(ranks, float), np.asarray(totals)
    # window may be one span for every site, or one per site. Deposits do not
    # accumulate over equal stretches of time, and assuming they do is a
    # substantive choice: a short window samples a narrow slice of the sequence
    # and a long one averages more of it away.
    win = np.asarray(window)
    if win.ndim == 0:
        win_arr = np.full(f.shape[1], int(win), dtype=int)
    elif win.ndim == 1 and win.shape[0] == f.shape[1]:
        win_arr = np.rint(win).astype(int)
    else:
        raise ValueError(
            f"window must be a scalar or one value per node ({f.shape[1]}), "
            f"got shape {win.shape}")
    if (win_arr < 1).any() or (win_arr > len(f)).any():
        raise ValueError("window must fit within the recorded sequence")
    if ranks.shape != (f.shape[1],) or totals.shape != ranks.shape:
        raise ValueError("one rank and total are required per node")
    if not np.isfinite(ranks).all() or (ranks < 0).any() or (ranks > 1).any():
        raise ValueError("ranks must lie in [0,1]")
    if (totals < 0).any() or not np.equal(totals, np.floor(totals)).all():
        raise ValueError("totals must be nonnegative integers")
    slots = win_arr - 1 + np.rint(ranks * (len(f) - win_arr)).astype(int)
    p = np.array([f[s - win_arr[i] + 1:s + 1, i].mean(0) for i, s in enumerate(slots)])
    p /= p.sum(1, keepdims=True)
    return np.array([rng.multinomial(int(n), q) for n, q in zip(totals, p)])
