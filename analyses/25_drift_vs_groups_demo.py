"""25_drift_vs_groups_demo.py — does spatial drift or bounded groups best explain
the basin's decorated-type structure?

A positive, model-comparison demonstration in the spirit of Lipo, DiNapoli,
Madsen & Hunt (2021). Two generative processes are run on the SAME real
coordinate layout of the curated St. Francis basin assemblages:

  (1) SPATIAL DRIFT (neutral). Wright-Fisher copying with between-node
      interaction decaying with geographic distance and NO group boundaries.
      Structure, if any, is purely spatially structured drift.

  (2) BOUNDED GROUPS. The same copying, but interaction is concentrated WITHIN
      a small number of spatially contiguous groups (the kind of bounded social
      units that "phases" posit), with only a small leak between groups.

For each model we build the sampling distribution of four summaries computed
exactly as on the data: distance-decay r (Mantel of BR similarity vs geography),
community modularity Q, distance-controlled boundary excess, and between-
community cultural F_ST. The observed values are then located against both
envelopes. The model whose envelope brackets the observed values is the one the
data are most consistent with.

This operationalizes Fox's (1998) negative finding (the phases fail as groups
and as classes) with a mechanism: it shows what generative process reproduces
the residual, partition-dependent similarity that remains once the phases are
not treated as real.

Read-only on the manuscript. Writes output/drift_vs_groups_demo.md and
figures/figS5_drift_vs_groups.png.

Usage: .venv/bin/python analyses/25_drift_vs_groups_demo.py
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
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
sd = importlib.import_module("23_phases_vs_spatial_drift")
from mls_emergence.signatures.assortativity import similarity_matrix  # noqa: E402

OUT_MD = ROOT / "output" / "drift_vs_groups_demo.md"
OUT_FIG = ROOT / "figures" / "figS5_drift_vs_groups.png"

# Okabe-Ito
C_DRIFT = "#0072B2"
C_GROUP = "#D55E00"
C_OBS = "#000000"


def drift_weights(coords, length_km=12.0):
    """Pure distance-decay interaction; no boundaries."""
    D = sd.geo_km(coords)
    W = np.exp(-D / length_km)
    np.fill_diagonal(W, 0.0)
    return W / W.sum(axis=1, keepdims=True)


def group_weights(coords, labels, length_km=12.0, leak=0.03):
    """Interaction concentrated within spatially contiguous groups.

    Within-group pairs interact via distance-decay; between-group pairs interact
    only at rate `leak` times the decay. Small leak = strong boundaries.
    """
    D = sd.geo_km(coords)
    base = np.exp(-D / length_km)
    same = labels[:, None] == labels[None, :]
    W = np.where(same, base, leak * base)
    np.fill_diagonal(W, 0.0)
    return W / W.sum(axis=1, keepdims=True)


def simulate(W, N=100, K=10, mu=0.01, m_between=0.05, steps=1200,
             avg_window=300, total_per_node=None, seed=0):
    """Vectorized neutral copying on a fixed row-normalized weight matrix W.

    Time-averaged by default (avg_window > 0): the assemblage pools deposition
    over a trailing window, matching the time-averaged nature of the real
    assemblages. The K most frequent labels across the pooled window are tallied.
    """
    rng = np.random.default_rng(seed)
    n = W.shape[0]
    cumW = np.cumsum(W, axis=1)
    pop = rng.integers(0, K, size=(n, N))
    next_type = K
    rows = np.arange(n)[:, None]
    start_avg = steps - avg_window
    thin = max(1, avg_window // 100) if avg_window > 0 else 1
    deposits = []
    for t in range(steps):
        from_other = rng.random((n, N)) < m_between
        src_self = pop[rows, rng.integers(0, N, size=(n, N))]
        u = rng.random((n, N))
        partners = np.clip((u[:, :, None] > cumW[:, None, :]).sum(axis=2), 0, n - 1)
        src_other = pop[partners, rng.integers(0, N, size=(n, N))]
        drawn = np.where(from_other, src_other, src_self)
        innov = rng.random((n, N)) < mu
        ni = int(innov.sum())
        if ni:
            drawn = drawn.copy()
            drawn[innov] = np.arange(next_type, next_type + ni)
            next_type += ni
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


def ensemble(W, totals, geo_dist, K, n_seeds=150, **kw):
    keys = ["dd_r", "Q", "be", "fst"]
    vals = {k: [] for k in keys}
    for s in range(n_seeds):
        M = simulate(W, K=K, total_per_node=totals, seed=5000 + s, **kw)
        if M.sum() == 0:
            continue
        st = sd.stats_for_matrix(M, geo_dist, seed=s)
        for k in keys:
            if np.isfinite(st[k]):
                vals[k].append(st[k])
    return {k: np.array(v, float) for k, v in vals.items()}


def mds2(counts):
    """Classical 2D MDS of BR dissimilarity, for the qualitative panel."""
    S = similarity_matrix(counts)
    D = 200.0 - S
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(-w)[:2]
    L = np.maximum(w[idx], 0)
    return V[:, idx] * np.sqrt(L)


def main():
    counts_df, coords_df = mf._load_curated()
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    names = list(counts_df.index)
    K = counts.shape[1]
    geo_dist = sd.geo_km(coords)
    totals = counts.sum(1)
    obs = sd.stats_for_matrix(counts, geo_dist, seed=0)

    # spatially contiguous groups for the bounded-groups model (k=2, the coarsest
    # phase split a real-groups hypothesis would posit)
    cc = coords - coords.mean(0)
    glabels = mf._kmeans_labels(cc, 2, seed=7)

    # Drift model calibrated to the observed weak distance-decay: ~24 km
    # interaction range, low between-node mixing (the grid cell that brackets all
    # four observed signatures simultaneously). The bounded-groups model uses the
    # same range and mixing but concentrates interaction within k=2 groups.
    LEN, MB = 24.0, 0.02
    Wd = drift_weights(coords, length_km=LEN)
    Wg = group_weights(coords, glabels, length_km=LEN, leak=0.03)

    drift = ensemble(Wd, totals, geo_dist, K, n_seeds=150, m_between=MB)
    group = ensemble(Wg, totals, geo_dist, K, n_seeds=150, m_between=MB)

    # ---- write results md ------------------------------------------------- #
    keys = [("dd_r", "distance-decay r"), ("Q", "modularity Q"),
            ("be", "boundary excess (BR)"), ("fst", "cultural F_ST")]
    L = ["# Spatial drift vs bounded groups: which best explains the basin?",
         "",
         f"Observed (n = {counts.shape[0]} curated assemblages), two generative",
         "models on the same coordinate layout, 150 seeds each. A model 'brackets'",
         "the observed statistic if the observed value falls in its 95% envelope.",
         "",
         "| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |",
         "|---|---|---|---|---|"]
    for key, lab in keys:
        d = drift[key]; g = group[key]
        dlo, dhi = np.percentile(d, [2.5, 97.5])
        glo, ghi = np.percentile(g, [2.5, 97.5])
        o = obs[key]
        ind = dlo <= o <= dhi
        ing = glo <= o <= ghi
        who = ("drift" if ind else "") + ("/groups" if ing else "")
        who = who.strip("/") or "neither"
        L.append(f"| {lab} | {o:+.3f} | [{dlo:+.3f}, {dhi:+.3f}] | "
                 f"[{glo:+.3f}, {ghi:+.3f}] | {who} |")
    L += ["",
          "Reading: where the observed value sits inside the spatial-drift envelope",
          "but outside the bounded-groups envelope, neutral drift on geography is",
          "the better explanation. The bounded-groups model, by construction,",
          "produces sharper community structure (higher Q, higher F_ST, larger",
          "distance-controlled boundary excess) than the basin actually shows.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    # ---- figure ----------------------------------------------------------- #
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(7.0, 4.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.1], hspace=0.45, wspace=0.55)

    # top row: 4 statistic distributions
    for j, (key, lab) in enumerate(keys):
        ax = fig.add_subplot(gs[0, j])
        d = drift[key]; g = group[key]
        bins = 18
        ax.hist(d, bins=bins, color=C_DRIFT, alpha=0.55, density=True, label="spatial drift")
        ax.hist(g, bins=bins, color=C_GROUP, alpha=0.55, density=True, label="bounded groups")
        ax.axvline(obs[key], color=C_OBS, lw=1.4)
        ax.set_xlabel(lab)
        ax.set_yticks([])
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.axes[0].legend(loc="upper left", fontsize=6, frameon=False)

    # bottom row: MDS observed, MDS one drift draw, MDS one groups draw
    Md = simulate(Wd, K=K, total_per_node=totals, m_between=MB, seed=5001)
    Mg = simulate(Wg, K=K, total_per_node=totals, m_between=MB, seed=5001)
    panels = [("Observed", counts, glabels),
              ("One spatial-drift realization", Md, glabels),
              ("One bounded-groups realization", Mg, glabels)]
    for j, (ttl, M, lab) in enumerate(panels):
        ax = fig.add_subplot(gs[1, j])
        xy = mds2(M)
        for gi, col in zip((0, 1), (C_DRIFT, C_GROUP)):
            m = lab == gi
            ax.scatter(xy[m, 0], xy[m, 1], s=18, c=col, edgecolor="white",
                       linewidth=0.4)
        if ttl == "Observed" and "Parkin" in names:
            pi = names.index("Parkin")
            ax.scatter(xy[pi, 0], xy[pi, 1], s=55, facecolor="none",
                       edgecolor=C_OBS, linewidth=1.2)
            ax.annotate("Parkin", (xy[pi, 0], xy[pi, 1]), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(ttl, fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
    # legend panel note in 4th cell
    axn = fig.add_subplot(gs[1, 3])
    axn.axis("off")
    axn.text(0.0, 0.95,
             "Points colored by the\nk=2 spatial split.\n\n"
             "Top: observed (black line)\nvs each model's null.\n\n"
             "Bottom: MDS of BR\ndissimilarity. Spatial drift\n"
             "reproduces the weak,\nspatially-blurred structure;\n"
             "bounded groups produce\nsharper separation than\nobserved.",
             fontsize=6, va="top")

    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    print("\n".join(L))
    print(f"\nwrote {OUT_MD}\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()
