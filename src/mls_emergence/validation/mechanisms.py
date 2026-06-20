"""Shared cultural-transmission simulator and the four mechanism generators.

The simulator produces ordinal "slices" of group-by-type counts under three
controllable knobs: between-group divergence, within-group conformity, and a
spatial rule. The four generators drive these knobs over an ordinal axis
t=0..T-1 to instantiate one genuine emergence process and three mimics.

Modeling rationale (and an honest caveat)
-----------------------------------------
A naive simulator that blends probability vectors,
``p = (1-d) * shared_pool + d * group_specific`` with a flat-Dirichlet pool,
fails on the neutrality signature. Two coupled artifacts appear:

1.  A flat or near-flat assemblage reads strongly ANTI-conformist on the
    Neiman/Ewens test, because the Ewens sampling formula expects a skewed
    abundance distribution. So the c=0 baseline is not neutral, and the
    departure metric is non-monotonic in conformity (it is U-shaped, passing
    through neutrality before turning conformist).
2.  Pooling divergent groups flattens the union assemblage, so between-group
    divergence pushes the POOLED neutrality signal the wrong way and cancels
    the within-group conformity it is supposed to reveal.

The fix used here keeps the model mechanistically honest while letting each
signature measure what it is meant to:

*   Every group's abundance SHAPE is fixed to a neutral-reading profile (a Zipf
    distribution, exponent ~2, which the Ewens test reads as theta_f/theta_e ~= 1
    at c=0). Divergence is expressed by reassigning which type IDENTITIES occupy
    which abundance ranks, not by flattening probabilities. So conformity is the
    only driver of the neutrality signature, and it drives it monotonically.
*   The neutrality signature is measured WITHIN groups and averaged, because
    conformity is a within-group transmission bias. Measuring it on the pooled
    assemblage conflates it with between-group divergence (artifact 2 above).
    This is a deliberate, documented departure from a literal "pooled counts"
    reading; see harness.py.
"""
from __future__ import annotations

import numpy as np

# Default ordinal/assemblage geometry shared by all generators.
T_DEFAULT = 8
N_GROUPS_DEFAULT = 12
N_PER_GROUP_DEFAULT = 300
N_TYPES_DEFAULT = 10

# Conformity sharpening gain: c in [0,1] maps to exponent 1 + c*K_SHARP.
K_SHARP = 4.0

# Zipf exponent for the neutral-reading abundance shape. ~2 yields theta_f/theta_e
# near 1 at c=0 for K=10, N=300 (verified empirically), i.e. a genuinely neutral
# baseline against which conformity registers.
ZIPF_EXP = 2.0


def _zipf_shape(n_types: int) -> np.ndarray:
    """Normalized Zipf rank-abundance profile over n_types classes (exponent
    ZIPF_EXP): the baseline abundance shape the generators place onto type
    orderings."""
    s = 1.0 / np.arange(1, n_types + 1) ** ZIPF_EXP
    return s / s.sum()


def _profile_from_order(order: np.ndarray, shape: np.ndarray, conformity: float) -> np.ndarray:
    """Place the fixed abundance shape onto a type-identity ordering, then sharpen.

    ``order[r]`` is the type id that occupies abundance rank r (rank 0 = most
    common). Conformity raises the profile to power (1 + c*K_SHARP).
    """
    p = np.zeros_like(shape)
    p[order] = shape
    c = float(conformity)
    if c > 0:
        p = p ** (1.0 + c * K_SHARP)
    s = p.sum()
    return p / s if s > 0 else np.full(shape.shape, 1.0 / shape.size)


def _swap_order(base: np.ndarray, n_swaps: int, rng: np.random.Generator) -> np.ndarray:
    """Adjacent-rank swaps applied to a base ordering.

    A small number of adjacent swaps produces an ordering close to ``base``;
    more swaps produce a more divergent ordering. This gives a graded, smooth
    notion of divergence usable for isolation-by-distance gradients.
    """
    o = np.asarray(base).copy()
    k = len(o)
    for _ in range(int(n_swaps)):
        i = int(rng.integers(0, k - 1))
        o[i], o[i + 1] = o[i + 1], o[i]
    return o


def _interpolate_order(
    shared_order: np.ndarray,
    target_order: np.ndarray,
    n_target_ranks: int,
    n_types: int,
) -> np.ndarray:
    """Build an ordering whose top n_target_ranks identities come from target_order
    and whose remaining ranks fall back to shared_order, skipping already-used
    types. n_target_ranks=0 returns shared_order; n_target_ranks=K returns
    target_order.
    """
    out: list[int] = []
    used: set[int] = set()
    for r in range(n_types):
        src = target_order if r < n_target_ranks else shared_order
        for ty in src:
            ty = int(ty)
            if ty not in used:
                used.add(ty)
                out.append(ty)
                break
    return np.array(out)


def _group_orders(
    n_groups: int,
    n_types: int,
    between_divergence: float,
    spatial_rule: str,
    coords: np.ndarray,
    rng: np.random.Generator,
    shared_order: np.ndarray | None = None,
) -> list[np.ndarray]:
    """Per-group type-identity orderings implementing divergence and the spatial rule.

    * "none":    each group's ordering is a fully independent permutation with
                 probability d (else it uses the shared ordering); coords ignored.
    * "bounded": groups are partitioned into a few spatial clusters; all groups in
                 a cluster share one ordering distinct from other clusters. The
                 strength of cluster distinctness scales with d. Sharp boundaries.
    * "ibd":     orderings vary SMOOTHLY along the dominant spatial axis via a
                 number of adjacent swaps proportional to position and to d. No
                 sharp edges; similarity decays continuously with distance.
    """
    d = float(between_divergence)
    if shared_order is None:
        shared_order = rng.permutation(n_types)

    if spatial_rule == "none":
        orders = []
        for _ in range(n_groups):
            if rng.random() < d:
                orders.append(rng.permutation(n_types))
            else:
                orders.append(shared_order.copy())
        return orders

    if spatial_rule == "bounded":
        n_clusters = max(2, min(4, n_groups // 3))
        labels = _spatial_clusters(coords, n_clusters, rng)
        # Fixed, mutually distinct base orderings per cluster, derived
        # deterministically from the shared order by cyclic shifts. Holding these
        # fixed across ordinal slices makes the same clusters progressively
        # differentiate (rather than reshuffling which types distinguish them),
        # which keeps the spatial-boundary signature stable and monotone.
        cluster_bases = [np.roll(shared_order, shift) for shift in range(n_clusters)]
        # Interpolate each group's ordering from the shared order toward its
        # cluster base: the first round(d*K) abundance ranks take their identity
        # from the cluster base, the remainder fall back to the shared order. d=0
        # collapses every group onto the shared order (no boundary); d=1 gives the
        # full cluster ordering (sharp boundary).
        n_cluster_ranks = int(round(d * n_types))
        orders = []
        for g in range(n_groups):
            base = cluster_bases[labels[g]]
            orders.append(_interpolate_order(shared_order, base, n_cluster_ranks, n_types))
        return orders

    if spatial_rule == "ibd":
        axis = np.asarray(coords, float)[:, 0]
        lo, hi = axis.min(), axis.max()
        frac = np.zeros_like(axis) if hi == lo else (axis - lo) / (hi - lo)
        # Adjacent swaps scale with position along the transect and with d, so
        # similarity decays smoothly with distance.
        orders = []
        for g in range(n_groups):
            n_swaps = int(round(frac[g] * d * 2.0 * n_types))
            orders.append(_swap_order(shared_order, n_swaps, rng))
        return orders

    raise ValueError(f"unknown spatial_rule: {spatial_rule!r}")


def _spatial_clusters(coords: np.ndarray, n_clusters: int, rng: np.random.Generator) -> np.ndarray:
    """Assign each coordinate to one of n_clusters via a short Lloyd's iteration."""
    coords = np.asarray(coords, float)
    n = coords.shape[0]
    n_clusters = min(n_clusters, n)
    centers = coords[rng.choice(n, size=n_clusters, replace=False)].copy()
    labels = np.zeros(n, dtype=int)
    for _ in range(25):
        dist = ((coords[:, None, :] - centers[None, :, :]) ** 2).sum(-1)
        new_labels = dist.argmin(axis=1)
        if np.array_equal(new_labels, labels):
            labels = new_labels
            break
        labels = new_labels
        for k in range(n_clusters):
            members = coords[labels == k]
            if len(members):
                centers[k] = members.mean(axis=0)
    return labels


def group_profiles(
    n_groups: int,
    n_types: int,
    between_divergence: float,
    within_conformity: float,
    spatial_rule: str,
    coords: np.ndarray,
    rng: np.random.Generator,
    orders: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Latent (noise-free) per-group sampling profiles, G x K, rows sum to 1."""
    shape = _zipf_shape(n_types)
    if orders is None:
        orders = _group_orders(
            n_groups, n_types, between_divergence, spatial_rule, coords, rng
        )
    return np.array(
        [_profile_from_order(orders[g], shape, within_conformity) for g in range(n_groups)]
    )


def sample_counts(profiles: np.ndarray, n_per_group: int, rng: np.random.Generator) -> np.ndarray:
    """Draw a multinomial assemblage of n_per_group sherds for each group profile;
    returns an integer group-by-type count matrix."""
    counts = np.zeros(profiles.shape, dtype=int)
    for g in range(profiles.shape[0]):
        counts[g] = rng.multinomial(n_per_group, profiles[g])
    return counts


def simulate_slice(
    n_groups: int,
    n_per_group: int,
    n_types: int,
    between_divergence: float,
    within_conformity: float,
    spatial_rule: str,
    coords: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate one ordinal slice of group-by-type counts (G x K integers)."""
    profiles = group_profiles(
        n_groups,
        n_types,
        between_divergence,
        within_conformity,
        spatial_rule,
        coords,
        rng,
    )
    return sample_counts(profiles, n_per_group, rng)


# --------------------------------------------------------------------------- #
# Coordinate layouts                                                          #
# --------------------------------------------------------------------------- #


def _clustered_coords(n_groups: int, rng: np.random.Generator, n_clusters: int = 3) -> np.ndarray:
    """Groups drawn around a few well-separated cluster centers."""
    offsets = np.array([[0, 0], [10, 10], [0, 10], [10, 0]])[:n_clusters]
    labels = np.repeat(np.arange(n_clusters), n_groups // n_clusters + 1)[:n_groups]
    coords = offsets[labels].astype(float) + rng.normal(0, 0.4, size=(n_groups, 2))
    return coords


def _line_coords(n_groups: int) -> np.ndarray:
    """Groups evenly spaced along a transect, for IBD gradients."""
    x = np.linspace(0, 10, n_groups)
    return np.column_stack([x, np.zeros_like(x)])


def _random_coords(n_groups: int, rng: np.random.Generator) -> np.ndarray:
    """Random 2-D coordinates (uniform in [0,10]^2), one per group, for layouts
    with no spatial structure."""
    return rng.uniform(0, 10, size=(n_groups, 2))


# --------------------------------------------------------------------------- #
# Generators                                                                   #
# --------------------------------------------------------------------------- #


def gen_group_emergence(seed: int):
    """Genuine group-level emergence.

    Between-group divergence and within-group conformity rise together over the
    ordinal axis on spatially bounded groups: clusters differentiate AND sharpen.
    All four signatures should rise jointly.
    """
    rng = np.random.default_rng(seed)
    coords = _clustered_coords(N_GROUPS_DEFAULT, rng)
    shared_order = rng.permutation(N_TYPES_DEFAULT)
    schedule = np.linspace(0.0, 0.9, T_DEFAULT)
    slices = []
    for t in range(T_DEFAULT):
        orders = _group_orders(
            N_GROUPS_DEFAULT,
            N_TYPES_DEFAULT,
            between_divergence=float(schedule[t]),
            spatial_rule="bounded",
            coords=coords,
            rng=rng,
            shared_order=shared_order,
        )
        profiles = group_profiles(
            N_GROUPS_DEFAULT,
            N_TYPES_DEFAULT,
            float(schedule[t]),
            float(schedule[t]),
            "bounded",
            coords,
            rng,
            orders=orders,
        )
        slices.append(sample_counts(profiles, N_PER_GROUP_DEFAULT, rng))
    return slices, coords


def gen_aggregated_conformity(seed: int):
    """Aggregated conformity mimic.

    One shared pool (no between-group divergence) but rising conformity. Produces
    a within-group neutrality departure with NO spatial structure and NO F_ST
    growth.
    """
    rng = np.random.default_rng(seed)
    coords = _random_coords(N_GROUPS_DEFAULT, rng)
    shared_order = rng.permutation(N_TYPES_DEFAULT)
    schedule = np.linspace(0.0, 0.9, T_DEFAULT)
    slices = []
    for t in range(T_DEFAULT):
        profiles = group_profiles(
            N_GROUPS_DEFAULT,
            N_TYPES_DEFAULT,
            0.0,
            float(schedule[t]),
            "none",
            coords,
            rng,
            orders=[shared_order.copy() for _ in range(N_GROUPS_DEFAULT)],
        )
        slices.append(sample_counts(profiles, N_PER_GROUP_DEFAULT, rng))
    return slices, coords


def gen_patchiness(seed: int):
    """Static spatial patchiness mimic.

    Constant moderate divergence on spatially bounded clusters, no conformity,
    constant over the ordinal axis. Spatial structure exists but nothing trends.
    The latent orderings are fixed across slices; only sampling noise varies.
    """
    rng = np.random.default_rng(seed)
    coords = _clustered_coords(N_GROUPS_DEFAULT, rng)
    shared_order = rng.permutation(N_TYPES_DEFAULT)
    orders = _group_orders(
        N_GROUPS_DEFAULT,
        N_TYPES_DEFAULT,
        between_divergence=0.6,
        spatial_rule="bounded",
        coords=coords,
        rng=rng,
        shared_order=shared_order,
    )
    slices = []
    for _ in range(T_DEFAULT):
        profiles = group_profiles(
            N_GROUPS_DEFAULT, N_TYPES_DEFAULT, 0.6, 0.0, "bounded", coords, rng, orders=orders
        )
        slices.append(sample_counts(profiles, N_PER_GROUP_DEFAULT, rng))
    return slices, coords


def gen_drift_space(seed: int):
    """Isolation-by-distance drift mimic.

    Smooth spatial gradient (IBD) along a transect, no conformity, constant over
    the ordinal axis. Distance structure is present but SMOOTH, not bounded, and
    assemblages stay neutral. Latent orderings fixed across slices.
    """
    rng = np.random.default_rng(seed)
    coords = _line_coords(N_GROUPS_DEFAULT)
    shared_order = rng.permutation(N_TYPES_DEFAULT)
    orders = _group_orders(
        N_GROUPS_DEFAULT,
        N_TYPES_DEFAULT,
        between_divergence=0.8,
        spatial_rule="ibd",
        coords=coords,
        rng=rng,
        shared_order=shared_order,
    )
    slices = []
    for _ in range(T_DEFAULT):
        profiles = group_profiles(
            N_GROUPS_DEFAULT, N_TYPES_DEFAULT, 0.8, 0.0, "ibd", coords, rng, orders=orders
        )
        slices.append(sample_counts(profiles, N_PER_GROUP_DEFAULT, rng))
    return slices, coords
