"""The questions, each a function of the data that returns a plain dictionary.

Uncertainty in the counts is carried by drawing every assemblage's class
proportions from Dirichlet(counts + 1/2) (the Jeffreys prior, which pulls each
profile slightly toward even and so leans neither way in a comparison between
divisions of the same assemblages). Every probability reported is a share of
those draws. Nothing here is a significance test.

Which function answers which question of the report:

  1  recovery_by_location     clustering on location and composition
  2  boundary_comparison      all phases against alternative divisions
  3  each_phase_against_rest  each phase against the rest, the same way
  4  own_phase_fit            each assemblage against each phase's profile
  5  phase_profiles           pooled profiles and their differences in sherds

`compare_division` does the work of questions 2 and 3. The boundary excess is
within-group minus between-group similarity at matched distance (see
`measures`): positive means the tested lines separate assemblages more than
distance alone predicts. F_ST is the Gini-Simpson cultural F_ST.

Seeds. For a base seed S (the command line's `--seed`, default 0),
`compare_division` uses S for the alternative divisions, S + 1 for the F_ST
draws and S + 2 for the boundary-excess draws, so changing the number of
draws of one never shifts the random numbers of another.
`each_phase_against_rest` passes S + 200 to every per-phase comparison,
keeping those streams apart from the all-phases comparison's S to S + 2.
Question 1's k-means uses seeds 0 to `kmeans_seeds - 1`, and question 4 its
own seed (96 by default); neither follows S.
"""
from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from . import divisions as dv
from . import measures as ms
from .data import Dataset

PER_DRAW_FST = 20        # alternatives of each kind compared per posterior draw
PER_DRAW_EXCESS = 10     # the same for the boundary excess, which is slower to compute
MAX_EXCESS_DRAWS = 500   # the boundary excess builds a similarity matrix per call
# Weights on composition relative to location in question 1, from the map alone
# (0) upward; composition alone (infinity) is appended by recovery_by_location.
WEIGHTS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 32.0)


def _q(a):
    """Posterior median and 2.5th and 97.5th percentiles (the 95 percent credible interval)."""
    a = np.asarray(a, float)
    return (float(np.median(a)), float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)))


def _draw(counts: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One posterior draw of the count table, each row scaled back to its sherd total.

    Each row's proportions are drawn from Dirichlet(counts + 1/2), the Jeffreys
    prior, so a class with no sherds can still receive a small share.
    """
    return np.array([rng.dirichlet(r + 0.5) for r in counts]) * counts.sum(1)[:, None]


def compare_division(counts, group_index, pts, dist, *, draws: int = 2000, n_alt: int = 300,
                     seed: int = 0) -> dict:
    """Does this division separate the assemblages better than other divisions of the same sizes?

    Returns the posterior of the division's cultural F_ST and boundary excess,
    the same for alternative divisions around random centers and compact ones,
    and the posterior probability that the tested division exceeds an
    alternative of each kind. A probability near one half means the tested
    lines do no better than arbitrary ones; near one means they do better; near
    zero means arbitrary lines separate the assemblages better.

    Parameters
    ----------
    counts : array_like, shape (n, classes)
        Sherd counts.
    group_index : array_like of int, shape (n,)
        The tested division, groups numbered 0..k-1, none empty.
    pts : numpy.ndarray, shape (n, 2)
        Flat coordinates in km, for building alternatives and spreads.
    dist : numpy.ndarray, shape (n, n)
        Distances in km for the boundary excess (may be along rivers).
    draws : int
        Posterior draws for F_ST, at least 100. The boundary excess uses
        at most `MAX_EXCESS_DRAWS` of them.
    n_alt : int
        Alternative divisions of each kind, at least 1.
    seed : int
        Base seed S; S, S + 1 and S + 2 are used (see the module docstring).

    Returns
    -------
    dict
        fst, fst_around, fst_compact : (median, 2.5%, 97.5%) over draws of
            the tested division's F_ST and of the median F_ST of the
            alternatives compared in each draw.
        p_fst_around, p_fst_compact : share of (draw, alternative)
            comparisons in which the tested division's F_ST is larger.
        fst_direct : F_ST on the observed counts, without draws.
        draws, n_alt : the settings used.
        spread : the tested division's spread (km^2 per assemblage).
        share_around_tighter, share_compact_tighter : share of alternatives
            with a smaller spread than the tested division.
        ari_around : median adjusted Rand index between the tested division
            and the random-center alternatives (what proximity alone gives).
        distinct_around, distinct_compact : number of distinct alternatives.
        excess, excess_around, excess_compact, p_excess_around,
        p_excess_compact, excess_direct : the same for the boundary excess
            (positive: within-group pairs more alike at matched distance).
        excess_bins : distance bins the observed excess rests on (0: no
            distance control).
        excess_uncontrolled : share of draws in which no bin qualified.
        excess_draws : draws used for the boundary excess.

    Raises
    ------
    ValueError
        If `draws` is below 100 or `n_alt` below 1.
    """
    if draws < 100:
        raise ValueError("at least 100 posterior draws are needed for the probabilities to mean anything")
    counts = np.asarray(counts, float)
    group_index = np.asarray(group_index)
    around, compact = dv.alternatives(pts, group_index, n_alt, np.random.default_rng(seed))
    rng = np.random.default_rng(seed + 1)           # F_ST draws: a stream of their own
    f_div, win_a, win_c, med_a, med_c = [], [], [], [], []
    for _ in range(draws):
        md = _draw(counts, rng)
        fp = ms.fst_by(md, group_index)
        # Each draw is compared with a fresh random subset of the alternatives (with
        # replacement), so all n_alt contribute over the draws at a fraction of the cost.
        fa = np.array([ms.fst_by(md, around[i]) for i in rng.choice(n_alt, PER_DRAW_FST)])
        fc = np.array([ms.fst_by(md, compact[i]) for i in rng.choice(n_alt, PER_DRAW_FST)])
        f_div.append(fp); win_a.append(np.mean(fp > fa)); win_c.append(np.mean(fp > fc))
        med_a.append(np.median(fa)); med_c.append(np.median(fc))
    out = {"fst": _q(f_div), "fst_around": _q(med_a), "fst_compact": _q(med_c),
           "p_fst_around": float(np.mean(win_a)), "p_fst_compact": float(np.mean(win_c)),
           "fst_direct": float(ms.fst_by(counts, group_index)), "draws": draws, "n_alt": n_alt,
           "spread": dv.spread(pts, group_index),
           "share_around_tighter": float(np.mean([dv.spread(pts, a) < dv.spread(pts, group_index) for a in around])),
           "share_compact_tighter": float(np.mean([dv.spread(pts, c) < dv.spread(pts, group_index) for c in compact])),
           "ari_around": float(np.median([ms.adjusted_rand(group_index, a) for a in around])),
           "distinct_around": len({_canon(a) for a in around}),
           "distinct_compact": len({_canon(c) for c in compact})}
    rng = np.random.default_rng(seed + 2)           # boundary-excess draws: a stream of their own
    b_div, b_a, b_c, bw_a, bw_c, uncontrolled = [], [], [], [], [], []
    for _ in range(min(draws, MAX_EXCESS_DRAWS)):
        md = _draw(counts, rng)
        bp, used = ms.boundary_excess_detail(md, dist, group_index)
        uncontrolled.append(used == 0)
        ba = np.array([ms.boundary_excess(md, dist, around[i]) for i in rng.choice(n_alt, PER_DRAW_EXCESS)])
        bc = np.array([ms.boundary_excess(md, dist, compact[i]) for i in rng.choice(n_alt, PER_DRAW_EXCESS)])
        b_div.append(bp); b_a.append(np.median(ba)); b_c.append(np.median(bc))
        bw_a.append(np.mean(bp > ba)); bw_c.append(np.mean(bp > bc))
    out.update(excess=_q(b_div), excess_around=_q(b_a), excess_compact=_q(b_c),
               p_excess_around=float(np.mean(bw_a)), p_excess_compact=float(np.mean(bw_c)),
               excess_direct=float(ms.boundary_excess(counts, dist, group_index)),
               excess_bins=int(ms.boundary_excess_detail(counts, dist, group_index)[1]),
               excess_uncontrolled=float(np.mean(uncontrolled)),
               excess_draws=min(draws, MAX_EXCESS_DRAWS))
    return out


def _canon(labels) -> tuple:
    """A division as a tuple that does not depend on how its groups are numbered."""
    seen = {}
    return tuple(seen.setdefault(int(v), len(seen)) for v in labels)


def boundary_comparison(data: Dataset, *, draws: int = 2000, n_alt: int = 300, seed: int = 0) -> dict:
    """Question 2: all phases together against alternative divisions of the same sizes.

    Parameters
    ----------
    data : Dataset
    draws, n_alt, seed : int
        As `compare_division`. The paper's worked example uses
        draws=2000, n_alt=300, seed=86000.

    Returns
    -------
    dict
        As `compare_division`, for the division into phases.
    """
    return compare_division(data.counts, data.phase_index, data.pts, data.dist,
                            draws=draws, n_alt=n_alt, seed=seed)


def each_phase_against_rest(data: Dataset, *, draws: int = 2000, n_alt: int = 300, seed: int = 0) -> dict:
    """Each phase against all the others taken together, beside same-size divisions.

    With only two phases each comparison is the division already tested by
    `boundary_comparison`, so an empty dictionary is returned.

    Parameters
    ----------
    data : Dataset
    draws, n_alt : int
        As `compare_division`.
    seed : int
        Base seed; every phase's comparison uses seed + 200 (see the module
        docstring), so each phase is scored against the same alternative
        two-way divisions when its size matches another's.

    Returns
    -------
    dict
        Phase name -> `compare_division` result for the division into that
        phase (1) and all the others (0); empty with two phases.
    """
    names = data.phase_names
    if len(names) == 2:
        return {}
    out = {}
    for i, p in enumerate(names):
        idx = (data.phase_index == i).astype(int)
        out[p] = compare_division(data.counts, idx, data.pts, data.dist,
                                  draws=draws, n_alt=n_alt, seed=seed + 200)
    return out


def _unit_scale(x: np.ndarray) -> np.ndarray:
    """Center, then scale so the mean squared distance from the center is 1."""
    x = np.asarray(x, float)
    x = x - x.mean(0)
    s = np.sqrt((x ** 2).sum(1).mean())
    return x / s if s > 0 else x


def _chisq_features(counts: np.ndarray) -> np.ndarray:
    """Row profiles divided by the square root of the mean profile (the chi-square metric)."""
    p = counts / counts.sum(1, keepdims=True)
    cm = p.mean(0)
    return p / np.sqrt(np.where(cm > 0, cm, 1.0))


def _kmeans(x: np.ndarray, k: int, seed: int, n_init: int) -> np.ndarray:
    """Lloyd's k-means, the lowest-spread result of `n_init` random starts."""
    rng = np.random.default_rng(seed)
    n = x.shape[0]
    best, best_inertia = None, np.inf
    for _ in range(n_init):
        centers = x[rng.choice(n, size=k, replace=False)].copy()
        labels = np.full(n, -1)
        for _ in range(50):                            # Lloyd iterations; stops early on convergence
            new = ((x[:, None, :] - centers[None, :, :]) ** 2).sum(-1).argmin(1)
            if np.array_equal(new, labels):
                break
            labels = new
            for j in range(k):
                members = x[labels == j]
                if len(members):                       # an emptied cluster keeps its old center
                    centers[j] = members.mean(0)
        inertia = float(sum(((x[labels == j] - centers[j]) ** 2).sum() for j in range(k)))
        if inertia < best_inertia:
            best_inertia, best = inertia, labels
    return best


def _cluster(x: np.ndarray, k: int, how: str, seed: int, n_init: int) -> np.ndarray:
    """Labels 0..k-1 from k-means, Ward or average linkage (seed and n_init used by k-means only)."""
    if how == "kmeans":
        return _kmeans(x, k, seed, n_init)
    if how in ("ward", "average"):
        return fcluster(linkage(pdist(x), method=how), k, criterion="maxclust") - 1
    raise ValueError(f"unknown clustering method: {how!r}")


def recovery_by_location(data: Dataset, *, weights=WEIGHTS, kmeans_seeds: int = 10, n_init: int = 500) -> dict:
    """Does clustering on location, composition, or both recover the phases?

    Assemblages are clustered into as many groups as there are phases, on site
    coordinates joined to composition (chi-square metric) multiplied by the
    square root of a weight. Weight 0 is the map alone; `inf` is composition
    alone. Agreement with the phases is the adjusted Rand index. K-means is
    reported as the median over `kmeans_seeds` seeds; Ward and average linkage
    have no random element.

    Parameters
    ----------
    data : Dataset
    weights : sequence of float
        Weights on composition (0 is the map alone); composition alone
        (infinity) is always added at the end.
    kmeans_seeds : int
        K-means seeds 0..kmeans_seeds-1; the median ARI is reported.
    n_init : int
        Random starts per k-means seed.

    Returns
    -------
    dict
        rows : one dict per weight with `weight`, the ARI for `kmeans`,
            `ward` and `average`, and for each method `<method>_misplaced`,
            the assemblages outside the cluster matched to their phase (for
            k-means, from the first seed).
        k : number of clusters, the number of phases.
    """
    idx, k = data.phase_index, len(data.phase_names)
    # Each block is scaled to unit mean squared distance, so weight 1 gives location and
    # composition equal say; sqrt(w) on the features makes w the weight on squared distance.
    g = _unit_scale(data.pts)
    c = _unit_scale(_chisq_features(data.counts.astype(float)))
    rows = []
    for w in list(weights) + [np.inf]:
        f = g if w == 0 else (c if np.isinf(w) else np.column_stack([g, np.sqrt(w) * c]))
        row = {"weight": float(w)}
        for how in ("kmeans", "ward", "average"):
            seeds = range(kmeans_seeds) if how == "kmeans" else [0]
            labs = [_cluster(f, k, how, s, n_init) for s in seeds]
            row[how] = float(np.median([ms.adjusted_rand(idx, lab) for lab in labs]))
            row[how + "_misplaced"] = _misplaced(idx, labs[0], data.names)
        rows.append(row)
    return {"rows": rows, "k": k}


def _misplaced(idx: np.ndarray, lab: np.ndarray, names: list) -> list:
    """Assemblages a clustering puts outside the cluster that best matches their phase."""
    from scipy.optimize import linear_sum_assignment
    # Match clusters to phases one to one by the largest overlap (Hungarian method on
    # the negated contingency table), then list those that fall outside their match.
    k = max(idx.max(), lab.max()) + 1
    c = np.zeros((k, k), int)
    for i, j in zip(idx, lab):
        c[i, j] += 1
    r, col = linear_sum_assignment(-c)
    match = dict(zip(r, col))
    return [n for n, i, j in zip(names, idx, lab) if match[i] != j]


def own_phase_fit(data: Dataset, *, draws: int = 4000, seed: int = 96) -> list:
    """Which phase does each assemblage's class profile fit best?

    The chi-square distance from an assemblage's profile to the pooled profile
    of each phase, with the assemblage left out of its own phase's pool.
    `p_own` is the share of posterior draws in which its own phase is the
    closest. `likeliest` is the phase closest in the most draws and
    `p_likeliest` its share; with more than two phases a low `p_own` does not
    by itself mean one other phase fits better.

    The point distances (`distance`, `closest`) use the observed profiles. In
    each draw the assemblage's profile and every phase pool are redrawn from
    Dirichlet(counts + 1/2). The chi-square metric weights each class by the
    mean profile of all assemblages.

    Parameters
    ----------
    data : Dataset
    draws : int
        Posterior draws per assemblage.
    seed : int
        Seed of the draws; fixed at 96 in the report, independent of
        `--seed`.

    Returns
    -------
    list of dict, one per assemblage in table order
        name, phase, sherds; distance (phase -> chi-square distance to its
        pool); closest (phase at the smallest distance); p_own; likeliest,
        p_likeliest; km_to_nearest (phase -> straight-line km to the nearest
        other assemblage of that phase); nearest_phase.
    """
    m = data.counts.astype(float)
    labels, phases = data.phases, data.phase_names
    p = m / m.sum(1, keepdims=True)
    cm = p.mean(0)
    w = np.sqrt(np.where(cm > 0, cm, 1.0))
    pts = data.pts
    rng = np.random.default_rng(seed)
    chisq = lambda a, b: float(np.sqrt((((a - b) / w) ** 2).sum()))
    rows = []
    for i, name in enumerate(data.names):
        others = {q: [j for j in range(len(m)) if labels[j] == q and j != i] for q in phases}
        pools = {q: m[others[q]].sum(0) for q in phases}
        dist = {q: chisq(p[i], pools[q] / pools[q].sum()) for q in phases}
        near = {q: min(float(np.linalg.norm(pts[i] - pts[j])) for j in others[q]) for q in phases}
        own = labels[i]
        wins = {q: 0 for q in phases}
        for _ in range(draws):
            pi = rng.dirichlet(m[i] + 0.5)
            d = {q: chisq(pi, rng.dirichlet(pools[q] + 0.5)) for q in phases}
            wins[min(d, key=d.get)] += 1
        likeliest = max(wins, key=wins.get)
        rows.append({"name": name, "phase": own, "sherds": int(m[i].sum()), "distance": dist,
                     "closest": min(dist, key=dist.get), "p_own": wins[own] / draws,
                     "likeliest": likeliest, "p_likeliest": wins[likeliest] / draws,
                     "km_to_nearest": near, "nearest_phase": min(near, key=near.get)})
    return rows


def phase_profiles(data: Dataset) -> dict:
    """Question 5: each phase's pooled class percentages, and the difference between every pair in sherds.

    Parameters
    ----------
    data : Dataset

    Returns
    -------
    dict
        profiles : array (phases x classes) of percentages, phases in
            `data.phase_names` order.
        overall : percentages of all sherds pooled.
        pairs : (phase a, phase b) -> share of sherds (percent, 0 to 100)
            that would have to change class to make the two alike.
        mean : mean of `pairs`.
    """
    prof = ms.pooled_profiles(data.counts, data.phase_index)
    names = data.phase_names
    pairs = {(names[a], names[b]): ms.profile_difference(prof[a], prof[b])
             for a in range(len(names)) for b in range(a + 1, len(names))}
    overall = 100.0 * data.counts.sum(0) / data.counts.sum()
    return {"profiles": prof, "overall": overall, "pairs": pairs, "mean": float(np.mean(list(pairs.values())))}
