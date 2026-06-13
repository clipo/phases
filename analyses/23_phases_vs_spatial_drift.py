"""23_phases_vs_spatial_drift.py — are the PFG "phases" real, or spatial-drift artifacts?

Tests whether the within>between similarity structure that underwrites the Parkin
phase (and Parkin as its central node) requires bounded social interaction, or
whether it is what neutral drift on a spatially structured network produces with
no groupness at all (Lipo, DiNapoli, Madsen & Hunt 2021, PLOS ONE).

Three parts, increasing in force:

  PART A  Distance-decay strength.
          Mantel of Brainerd-Robinson similarity vs geographic distance.
          How much of the ceramic patterning is geography alone?

  PART B  Excess beyond isolation-by-distance.
          Detect ceramic communities (the data-driven analog of phases) by
          modularity maximization on the BR-similarity graph, then compute the
          distance-controlled within-minus-between similarity on those
          communities. Under pure IBD this excess collapses to ~0; a genuine
          interaction boundary leaves a positive excess. Partial Mantel of
          ceramic distance vs community membership, controlling for geography.

  PART C  Generative spatial-drift null (Lipo et al. 2021).
          Simulate neutral drift on the ACTUAL coordinate layout with
          distance-decayed between-node interaction and NO imposed boundaries.
          Run the same community detection + boundary-excess + distance-decay on
          the synthetic assemblages over many seeds. If the observed values fall
          inside this null, the phase structure needs no social boundary.

Read-only on the manuscript. Writes output/phases_vs_spatial_drift.md.

Usage: .venv/bin/python analyses/23_phases_vs_spatial_drift.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import networkx as nx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import make_figures as mf  # noqa: E402
from mls_emergence.signatures.assortativity import (  # noqa: E402
    similarity_matrix, mantel,
)

OUT = ROOT / "output" / "phases_vs_spatial_drift.md"
EARTH_KM_PER_DEG = 111.32


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
def geo_km(coords: np.ndarray) -> np.ndarray:
    """Pairwise great-circle-ish distance (km) via equirectangular approximation.

    coords columns are [Latitude, Longitude]. At ~35 N over a ~100 km basin the
    equirectangular error is well under 1 percent, far below what matters here.
    """
    lat = coords[:, 0]
    lon = coords[:, 1]
    mlat = np.deg2rad(lat.mean())
    x = np.deg2rad(lon) * np.cos(mlat) * EARTH_KM_PER_DEG / np.deg2rad(1)
    y = np.deg2rad(lat) * EARTH_KM_PER_DEG / np.deg2rad(1)
    p = np.column_stack([x, y])
    d = p[:, None, :] - p[None, :, :]
    return np.sqrt((d ** 2).sum(-1))


# --------------------------------------------------------------------------- #
# community detection on the ceramic-similarity graph
# --------------------------------------------------------------------------- #
def ceramic_communities(counts: np.ndarray, seed: int = 0):
    """Greedy-modularity communities on the BR-similarity graph.

    Edge weight = Brainerd-Robinson similarity (0..200). Returns (labels, Q,
    n_communities). This is the data-driven stand-in for "phases": groups of
    assemblages whose decorated-type frequencies are mutually more similar.
    """
    S = similarity_matrix(counts)
    n = S.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            w = S[i, j]
            if w > 0:
                G.add_edge(i, j, weight=w)
    comms = nx.community.greedy_modularity_communities(G, weight="weight")
    labels = np.full(n, -1)
    for c, members in enumerate(comms):
        for m in members:
            labels[m] = c
    Q = nx.community.modularity(G, comms, weight="weight")
    return labels, float(Q), len(comms)


def boundary_excess_labeled(counts: np.ndarray, coords_km_dist: np.ndarray,
                            labels: np.ndarray, n_bins: int = 4) -> float:
    """Distance-controlled within-minus-between BR similarity for SUPPLIED labels.

    Mirrors signatures.assortativity.boundary_excess but takes an external
    partition (here, the ceramic communities) instead of clustering coordinates.
    Within each geographic-distance bin, compares mean BR similarity of
    within-community pairs to between-community pairs at the SAME distance.
    Under pure isolation-by-distance the bin-matched gap is ~0.
    """
    S = similarity_matrix(counts)
    n = S.shape[0]
    iu = np.triu_indices(n, k=1)
    d = coords_km_dist[iu]
    s = S[iu]
    same = labels[iu[0]] == labels[iu[1]]
    if same.sum() == 0 or (~same).sum() == 0:
        return 0.0
    raw_gap = float(s[same].mean() - s[~same].mean())
    if d.max() == d.min():
        return raw_gap
    edges = np.linspace(d.min(), d.max() + 1e-9, n_bins + 1)
    gaps = []
    for b in range(n_bins):
        in_bin = (d >= edges[b]) & (d < edges[b + 1])
        w = in_bin & same
        btw = in_bin & ~same
        if w.sum() == 0 or btw.sum() == 0:
            continue
        gaps.append(float(s[w].mean() - s[btw].mean()))
    return float(np.mean(gaps)) if gaps else raw_gap


def partial_mantel(cer_dist, geo_dist, memb_dist, n_perm=4999, seed=0):
    """Partial Mantel: corr(cer_dist, memb_dist | geo_dist), permutation p.

    Residualize ceramic distance and membership distance on geographic distance,
    correlate the residuals, then permute assemblage labels to get a null. Tests
    whether community membership predicts ceramic dissimilarity BEYOND geography.
    """
    iu = np.triu_indices_from(cer_dist, k=1)
    y = cer_dist[iu]
    g = geo_dist[iu]
    m = memb_dist[iu]

    def resid(v, on):
        A = np.column_stack([np.ones_like(on), on])
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    ry = resid(y, g)
    rm = resid(m, g)
    r_obs = np.corrcoef(ry, rm)[0, 1]
    rng = np.random.default_rng(seed)
    n = cer_dist.shape[0]
    count = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        mp = memb_dist[p][:, p][iu]
        rmp = resid(mp, g)
        if abs(np.corrcoef(ry, rmp)[0, 1]) >= abs(r_obs):
            count += 1
    return float(r_obs), float((count + 1) / (n_perm + 1))


# --------------------------------------------------------------------------- #
# PART C: neutral drift on the real spatial layout (Lipo et al. 2021)
# --------------------------------------------------------------------------- #
def simulate_spatial_drift(coords, N=100, K=10, mu=0.01, length_km=12.0,
                           m_between=0.05, steps=1500, total_per_node=None,
                           seed=0):
    """Wright-Fisher neutral drift on a distance-decayed interaction network.

    Each assemblage is a node holding N individuals carrying one of K type
    labels. Each step every individual copies a type either from its own node
    (prob 1 - m_between) or from another node chosen with probability
    proportional to exp(-d/length_km) (prob m_between), with innovation rate mu.
    No group boundaries are imposed: structure is purely the spatial decay of
    interaction. After `steps`, type frequencies per node are sampled to integer
    counts matching the observed per-assemblage sherd totals (so the synthetic
    matrix has the same sampling resolution as the data).

    Returns a (n_nodes x K) integer count matrix.
    """
    rng = np.random.default_rng(seed)
    n = coords.shape[0]
    D = geo_km(coords)
    W = np.exp(-D / length_km)
    np.fill_diagonal(W, 0.0)
    W = W / W.sum(axis=1, keepdims=True)

    # state: type label of each individual at each node
    pop = rng.integers(0, K, size=(n, N))
    next_type = K  # innovations append new labels; cap richness by relabeling

    for _ in range(steps):
        new = pop.copy()
        for i in range(n):
            # source node for each of N copy events
            from_other = rng.random(N) < m_between
            # within-node copies
            src_self = pop[i, rng.integers(0, N, size=N)]
            # between-node copies
            partners = rng.choice(n, size=N, p=W[i])
            src_other = pop[partners, rng.integers(0, N, size=N)]
            drawn = np.where(from_other, src_other, src_self)
            # innovation
            innov = rng.random(N) < mu
            n_innov = int(innov.sum())
            if n_innov:
                drawn = drawn.copy()
                drawn[innov] = np.arange(next_type, next_type + n_innov)
                next_type += n_innov
            new[i] = drawn
        pop = new

    # tabulate to a fixed K-column matrix by keeping the K most frequent labels
    flat = pop.ravel()
    labels, freqs = np.unique(flat, return_counts=True)
    keep = labels[np.argsort(-freqs)[:K]]
    keep_index = {lab: c for c, lab in enumerate(keep)}
    M = np.zeros((n, K))
    for i in range(n):
        u, cnt = np.unique(pop[i], return_counts=True)
        for lab, c in zip(u, cnt):
            if lab in keep_index:
                M[i, keep_index[lab]] += c
    # rescale each node to its observed sherd total, multinomial-sampled
    if total_per_node is not None:
        out = np.zeros((n, K), dtype=int)
        for i in range(n):
            p = M[i] / M[i].sum() if M[i].sum() else np.ones(K) / K
            tot = int(total_per_node[i])
            out[i] = rng.multinomial(tot, p) if tot > 0 else 0
        return out
    return M.astype(int)


def stats_for_matrix(counts, geo_dist, seed=0):
    """Distance-decay r, modularity Q, n communities, boundary excess, F_ST."""
    S = similarity_matrix(counts)
    r_dd, _ = mantel(S, geo_dist, n_perm=1, seed=seed)  # r only; p not needed here
    labels, Q, ncom = ceramic_communities(counts, seed=seed)
    be = boundary_excess_labeled(counts, geo_dist, labels)
    fst = mf.cultural_fst(np.array([counts[labels == c].sum(0)
                                    for c in np.unique(labels)
                                    if counts[labels == c].sum() > 0])) \
        if len(np.unique(labels)) >= 2 else np.nan
    return dict(dd_r=r_dd, Q=Q, ncom=ncom, be=be, fst=float(fst))


# --------------------------------------------------------------------------- #
def main():
    counts_df, coords_df = mf._load_curated()
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    names = list(counts_df.index)
    n = counts.shape[0]
    geo_dist = geo_km(coords)
    totals = counts.sum(1)

    L = [f"# Are the PFG phases real, or spatial-drift artifacts? (n = {n} curated assemblages)",
         "",
         "Test of whether the within>between decorated-type similarity structure",
         "underwriting the Parkin phase requires bounded social interaction, or is",
         "what neutral drift on a spatially structured network produces with no",
         "groupness (Lipo et al. 2021).", ""]

    # ---- PART A: distance-decay strength ---------------------------------- #
    S = similarity_matrix(counts)
    r_dd, p_dd = mantel(S, geo_dist, n_perm=4999, seed=1)
    L += ["## Part A. Distance-decay strength",
          f"- Mantel r (BR similarity vs geographic distance): {r_dd:+.3f} (p = {p_dd:.4f}).",
          f"- Geographic distance alone explains r^2 = {r_dd**2:.2f} of the pairwise",
          "  ceramic similarity variance. Negative r = closer assemblages are more",
          "  similar, the isolation-by-distance signature.", ""]

    # ---- PART B: excess beyond IBD --------------------------------------- #
    labels, Q, ncom = ceramic_communities(counts, seed=0)
    be_obs = boundary_excess_labeled(counts, geo_dist, labels)
    raw_gap = float(S[np.triu_indices(n, 1)][labels[np.triu_indices(n,1)[0]] ==
                    labels[np.triu_indices(n,1)[1]]].mean()
                    - S[np.triu_indices(n, 1)][labels[np.triu_indices(n,1)[0]] !=
                    labels[np.triu_indices(n,1)[1]]].mean())
    memb_dist = (labels[:, None] != labels[None, :]).astype(float)
    cer_dist = 200.0 - S
    r_pm, p_pm = partial_mantel(cer_dist, geo_dist, memb_dist, n_perm=4999, seed=2)

    parkin_comm = labels[names.index("Parkin")] if "Parkin" in names else None
    comm_sizes = {int(c): int((labels == c).sum()) for c in np.unique(labels)}
    L += ["## Part B. Excess beyond isolation-by-distance",
          f"- Ceramic communities (greedy modularity on BR graph): {ncom} communities, "
          f"modularity Q = {Q:.3f}; sizes {comm_sizes}.",
          f"- Raw within-minus-between BR similarity gap: {raw_gap:+.1f} BR units.",
          f"- Distance-CONTROLLED within-minus-between gap (boundary excess): "
          f"{be_obs:+.1f} BR units.",
          "  Under pure isolation-by-distance this collapses to ~0; a genuine",
          "  interaction boundary leaves a positive excess.",
          f"- Partial Mantel, ceramic distance vs community membership controlling "
          f"for geography: r = {r_pm:+.3f} (p = {p_pm:.4f}).",
          f"- Parkin falls in community {parkin_comm}.", ""]

    # ---- PART C: generative spatial-drift null --------------------------- #
    obs = stats_for_matrix(counts, geo_dist, seed=0)
    n_seeds = 200
    null = []
    for s in range(n_seeds):
        M = simulate_spatial_drift(coords, N=100, K=counts.shape[1], mu=0.01,
                                   length_km=12.0, m_between=0.05, steps=1500,
                                   total_per_node=totals, seed=1000 + s)
        if M.sum() == 0:
            continue
        null.append(stats_for_matrix(M, geo_dist, seed=s))

    def summ(key):
        v = np.array([d[key] for d in null], float)
        v = v[np.isfinite(v)]
        lo, hi = np.percentile(v, [2.5, 97.5])
        o = obs[key]
        # fraction of null at least as extreme as observed (two-sided-ish)
        frac_ge = float((v >= o).mean())
        return v.mean(), lo, hi, o, frac_ge

    L += ["## Part C. Generative spatial-drift null (Lipo et al. 2021)",
          "Neutral drift on the real coordinate layout, distance-decayed",
          "interaction (length scale 12 km, between-node copy rate 0.05,",
          "innovation 0.01), no imposed boundaries. 95% null interval over",
          f"{len(null)} seeds vs observed:", "",
          "| statistic | observed | null mean | null 95% | obs inside null? |",
          "|---|---|---|---|---|"]
    for key, lab in [("dd_r", "distance-decay r"),
                     ("Q", "modularity Q"),
                     ("be", "boundary excess (BR)"),
                     ("fst", "cultural F_ST")]:
        mean, lo, hi, o, _ = summ(key)
        inside = "yes" if lo <= o <= hi else "NO"
        L += [f"| {lab} | {o:+.3f} | {mean:+.3f} | [{lo:+.3f}, {hi:+.3f}] | {inside} |"]
    L += ["",
          "Interpretation: where observed falls INSIDE the spatial-drift null,",
          "that aspect of the phase structure is reproduced by neutral drift on",
          "geography alone and needs no social boundary. Where it falls OUTSIDE,",
          "the data carry structure beyond spatially structured drift.", ""]

    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
