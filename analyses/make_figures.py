"""make_figures.py — house-style figure pipeline for mls-emergence.

Generates nine figures (F2-F7, S2, S3, S5) to figures/ using
analyses/figstyle.py house style (Okabe-Ito, DejaVu Sans, 300 dpi, 7 in).
Reuses data-loading and computation logic from prior analysis scripts
(analyses/04-08); does not recompute from scratch where avoidable.

Data policy: data/ is gitignored and location-sensitive. This script NEVER
prints raw coordinates to stdout. Figures use centered/relative frames with
no axis-scale tick labels where required by the data policy.

Usage:
    .venv/bin/python analyses/make_figures.py
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr

# Project root and sys.path setup so `analyses` is importable as a package.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

from figstyle import save, save_all, OI_BLUE, OI_ORANGE, OI_GREEN, OI_VERMIL, OI_SKY, OI_PURPLE, OI_BLACK, OIC_BLUE, OIC_VERMIL  # noqa: F401  (save_all re-exported as mf.save_all)

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from mls_emergence.dataio.pfg import load_pfg_counts
from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv, normalize_grid
from mls_emergence.signatures.neutral import theta_f, theta_e
from mls_emergence.signatures.variance import cultural_fst
from mls_emergence.signatures.assortativity import boundary_excess, _kmeans_labels, geo_distance
from mls_emergence.signatures.seriation import seriation_groups
from mls_emergence.validation.harness import (
    run_blind, discriminates, signatures_over_axis, SIGNATURE_COLUMNS,
)
from mls_emergence.validation.mechanisms import (
    gen_group_emergence, gen_aggregated_signaling, gen_patchiness, gen_drift_space,
)

DATA = ROOT / "data"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

PARKIN_CUR = "Parkin"
PARKIN_BROAD = "11-N-1"
DECORATED_TYPES = [
    "Parkin_Punctated", "Barton/Kent/MPI", "Painted", "Fortune_Noded",
    "Ranch_Incised", "Walls_Engraved", "Wallace_Incised", "Rhodes_Incised",
    "Vernon_Paul_Applique", "Hull_Engraved",
]
CONT_PRIMARY = 0.10

# ---------------------------------------------------------------------------
# Shared helpers (duplicated from 06/07 to keep this script self-contained)
# ---------------------------------------------------------------------------

def correspondence_axis(M: np.ndarray):
    """Correspondence analysis; return (row_ordinate_dim1, row_ordinate_dim2, inertia_frac_dim1)."""
    M = np.asarray(M, float)
    total = M.sum()
    P = M / total
    r = P.sum(axis=1)
    c = P.sum(axis=0)
    keep_c = c > 0
    P = P[:, keep_c]
    c = c[keep_c]
    Dr_inv = np.diag(1.0 / np.sqrt(r))
    Dc_inv = np.diag(1.0 / np.sqrt(c))
    S = Dr_inv @ (P - np.outer(r, c)) @ Dc_inv
    U, sig, Vt = np.linalg.svd(S, full_matrices=False)
    row_scores = Dr_inv @ U @ np.diag(sig)
    inertia = sig ** 2
    frac = float(inertia[0] / inertia.sum()) if inertia.sum() > 0 else np.nan
    return row_scores[:, 0], row_scores[:, 1], frac


def parse_cal_midpoint(s) -> float:
    if pd.isna(s):
        return np.nan
    nums = [int(x) for x in re.findall(r"\d+", str(s))]
    if len(nums) < 2:
        return np.nan
    return (nums[0] + nums[-1]) / 2.0


def norm_name(s) -> str:
    base = re.sub(r"[^a-z0-9 ]", " ", str(s).lower())
    tokens = [t for t in base.split() if t]
    drop = {"place", "mound", "mounds", "lake", "landing", "ferry",
            "bayou", "ranch", "the", "site"}
    core = [t for t in tokens if t not in drop]
    return "".join(core if core else tokens)


def silhouette_mean(coords: np.ndarray, labels: np.ndarray) -> float:
    n = coords.shape[0]
    d = geo_distance(coords)
    uniq = np.unique(labels)
    if len(uniq) < 2:
        return -1.0
    s = np.zeros(n)
    for i in range(n):
        same = labels == labels[i]
        same[i] = False
        if same.sum() == 0:
            s[i] = 0.0
            continue
        a = d[i, same].mean()
        b = np.inf
        for cc in uniq:
            if cc == labels[i]:
                continue
            mask = labels == cc
            if mask.sum():
                b = min(b, d[i, mask].mean())
        s[i] = (b - a) / max(a, b) if max(a, b) > 0 else 0.0
    return float(s.mean())


def ols_slope(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2:
        return np.nan
    xc = x[m] - x[m].mean()
    if (xc @ xc) == 0:
        return np.nan
    return float((xc @ (y[m] - y[m].mean())) / (xc @ xc))


def zscore_series(s: pd.Series) -> pd.Series:
    s = s.astype(float)
    sd = float(s.std(ddof=0))
    if sd > 0 and np.isfinite(sd):
        return (s - s.mean()) / sd
    return s * 0.0


# ---------------------------------------------------------------------------
# Load shared curated data once
# ---------------------------------------------------------------------------
def _load_curated():
    cur = pd.read_csv(
        DATA / "raw" / "mainfort-pfg-cpl.csv"
    ).dropna(subset=["Assemblages"])
    cur["Assemblages"] = cur["Assemblages"].astype(str).str.strip()
    cur = cur.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    type_cols = [c for c in DECORATED_TYPES if c in cur.columns]
    counts = cur[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    row_tot = counts.sum(axis=1)
    counts = counts[row_tot > 0]
    # Coordinates
    xy = pd.read_csv(DATA / "raw" / "mainfort-pfg-cplXY.txt", sep="\t")
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    coords_ll = xy.reindex(counts.index)[["Latitude", "Longitude"]].apply(
        pd.to_numeric, errors="coerce"
    )
    # Restrict to the St. Francis basin, defined HYDROLOGICALLY as assemblages
    # within 20 km of the St. Francis / Tyronza / L'Anguille drainage (the
    # canonical member list is generated by analyses/16_basin_membership.py).
    # This replaces the earlier latitude cut, which admitted Mississippi-River
    # sites that are not in the St. Francis drainage.
    members = _basin_members("curated")
    keep = [a for a in counts.index if a in members]
    counts = counts.loc[keep]
    coords_ll = coords_ll.loc[keep]
    return counts, coords_ll


def _basin_members(which: str) -> set:
    """Canonical drainage-basin membership (analyses/16_basin_membership.py)."""
    f = DATA / "processed" / f"basin_members_{which}.txt"
    return {ln.strip() for ln in f.read_text().splitlines() if ln.strip()}


# ---------------------------------------------------------------------------
# Idealized criterion validation: four-signature trajectories under the four
# generative mechanisms, with the signature-independence audit inset (Figure S4).
# ---------------------------------------------------------------------------
def fig4_validation() -> None:
    GENERATORS = {
        "group_emergence": gen_group_emergence,
        "aggregated_signaling": gen_aggregated_signaling,
        "patchiness": gen_patchiness,
        "drift_space": gen_drift_space,
    }
    AUDIT_SEEDS = list(range(20))
    REPORT_SEED = 42

    panels = run_blind(GENERATORS, seed=REPORT_SEED)
    verdict = discriminates(panels, deriv_threshold=0.10)

    # Signature-independence audit: mean |r| among the four signatures.
    def pooled_corr(gen, seeds):
        frames = [signatures_over_axis(*gen(s)) for s in seeds]
        return pd.concat(frames, ignore_index=True).corr()

    def mean_offdiag(corr: pd.DataFrame) -> float:
        v = corr.values
        return float(np.nanmean(np.abs(v[np.triu_indices(len(v), 1)])))

    genuine_corr = pooled_corr(gen_group_emergence, AUDIT_SEEDS[:8])
    mimic_corrs = {
        "aggregated_signaling": pooled_corr(gen_aggregated_signaling, AUDIT_SEEDS[:8]),
        "patchiness": pooled_corr(gen_patchiness, AUDIT_SEEDS[:8]),
        "drift_space": pooled_corr(gen_drift_space, AUDIT_SEEDS[:8]),
    }
    genuine_meanabs = mean_offdiag(genuine_corr)
    mimic_meanabs = {k: mean_offdiag(v) for k, v in mimic_corrs.items()}

    LABELS = {
        "group_emergence": "Group emergence\n(genuine)",
        "aggregated_signaling": "Aggregated\nsignaling",
        "patchiness": "Spatial\npatchiness",
        "drift_space": "Isolation by\ndistance",
    }
    SIG_LABELS = ["Neutral departure", "Seriability", "Cultural F_ST", "Spatial boundary"]
    # Color (online-only supplement); marker and line style also vary per signature.
    SIG_COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00"]
    SIG_MARKERS = ["o", "s", "^", "D"]
    SIG_LS = ["-", "--", ":", "-."]
    MECH_ORDER = ["group_emergence", "aggregated_signaling", "patchiness", "drift_space"]

    fig = plt.figure(figsize=(7, 6.5))
    gs_main = fig.add_gridspec(2, 2, left=0.07, right=0.65, hspace=0.50, wspace=0.45)
    gs_inset = fig.add_gridspec(1, 1, left=0.72, right=0.98, top=0.88, bottom=0.18)
    axes_main = [fig.add_subplot(gs_main[i, j]) for i in range(2) for j in range(2)]

    for ax, mech in zip(axes_main, MECH_ORDER):
        panel = panels[mech]
        x = panel.index.to_numpy(float)
        for col, label, color, mk, ls in zip(SIGNATURE_COLUMNS, SIG_LABELS,
                                              SIG_COLORS, SIG_MARKERS, SIG_LS):
            ax.plot(x, panel[col], marker=mk, markersize=3, color=color,
                    label=label, linewidth=1.4, linestyle=ls)
        conv = verdict[mech]["convergent"]
        ax.set_title(
            LABELS[mech] + ("\n* CONVERGENT" if conv else ""),
            fontsize=8, pad=3, color=OIC_VERMIL if conv else "black",
        )
        ax.set_xlabel("Ordinal step", fontsize=7)
        ax.set_ylabel("Signature value", fontsize=7)
        ax.set_xticks(range(len(x)))
        ax.tick_params(labelsize=7)
        ax.axhline(0, color="0.8", linewidth=0.5, zorder=0)

    handles, lbs = axes_main[0].get_legend_handles_labels()
    fig.legend(handles, lbs, loc="lower center", bbox_to_anchor=(0.36, -0.01),
               ncol=4, fontsize=7, frameon=False, handlelength=1.5)

    ax_ins = fig.add_subplot(gs_inset[0, 0])
    bar_vals = [genuine_meanabs] + list(mimic_meanabs.values())
    bar_colors = [OIC_VERMIL] + [OIC_BLUE] * 3
    bpos = np.arange(len(bar_vals))
    ax_ins.bar(bpos, bar_vals, color=bar_colors, width=0.6, edgecolor="none")
    ax_ins.set_xticks(bpos)
    short_labels = ["Genuine\nemergence", "Agg.\nsignaling", "Patchiness", "Drift\n(IBD)"]
    ax_ins.set_xticklabels(short_labels, fontsize=6)
    ax_ins.set_ylabel("Mean |r| among\nfour signatures", fontsize=7)
    ax_ins.tick_params(labelsize=6)
    ax_ins.axhline(0, color="0.8", linewidth=0.5)

    save(fig, "figS4_validation")
    print("figS4_validation.png written")


# ---------------------------------------------------------------------------
# F3: CA ordination of the curated decorated set
# ---------------------------------------------------------------------------
def fig3_ca_ordination() -> None:
    counts, coords_ll = _load_curated()
    M = counts.to_numpy(float)
    ca1, ca2, inertia_frac1 = correspondence_axis(M)

    # Orient CA1 with 14C: larger = later
    rc = pd.read_csv(DATA / "raw" / "14CDatesFromMainfort2001.csv")
    rc = rc[rc["Provenience"].notna() & (rc["Provenience"] != "Provenience")].copy()
    rc["cal_mid"] = rc["Calibrated Date A.D. (1 Sigma)"].map(parse_cal_midpoint)
    prov_date = rc.groupby("Provenience")["cal_mid"].agg(["mean", "count"])
    cur_norm = {norm_name(a): a for a in counts.index}
    assem_date: dict[str, float] = {}
    for prov in prov_date.index:
        pn = norm_name(prov)
        hit = cur_norm.get(pn)
        if hit is None:
            for k, v in cur_norm.items():
                if pn and (pn in k or k in pn):
                    hit = v
                    break
        if hit is not None:
            assem_date.setdefault(hit, []).append(float(prov_date.loc[prov, "mean"]))  # type: ignore
    assem_date = {a: float(np.mean(v)) for a, v in assem_date.items()}  # type: ignore
    dated = [a for a in counts.index if a in assem_date]
    ca1_d = ca1[[list(counts.index).index(a) for a in dated]]
    yr_d = np.array([assem_date[a] for a in dated])
    if len(dated) >= 3:
        rho, _ = spearmanr(ca1_d, yr_d)
        if np.isfinite(rho) and rho < 0:
            ca1 = -ca1

    ca1_s = pd.Series(ca1, index=counts.index)
    ca2_s = pd.Series(ca2, index=counts.index)

    fig, ax = plt.subplots(figsize=(7, 5))

    # Shade points by CA1 rank (continuous gradient from early to late). Grayscale
    # (truncated so the earliest points are a visible mid-gray, not white).
    from matplotlib.colors import LinearSegmentedColormap
    ca1_rank = ca1_s.rank(pct=True).values
    cmap = LinearSegmentedColormap.from_list(
        "greys_t", plt.cm.Greys(np.linspace(0.25, 1.0, 256)))
    sc = ax.scatter(
        ca1_s.values, ca2_s.values,
        c=ca1_rank, cmap=cmap,
        s=28, alpha=0.9, edgecolors="black", linewidths=0.3, zorder=3,
    )
    # Parkin highlighted
    pk_x = float(ca1_s[PARKIN_CUR])
    pk_y = float(ca2_s[PARKIN_CUR])
    ax.scatter([pk_x], [pk_y], s=170, marker="*",
               facecolor="white", edgecolors="black", linewidths=0.9,
               label="Parkin", zorder=6)
    ax.annotate("Parkin", (pk_x, pk_y), fontsize=7,
                xytext=(5, 5), textcoords="offset points", color=OI_VERMIL)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("CA1 rank (0=early, 1=late)", fontsize=7)
    cbar.ax.tick_params(labelsize=7)

    ax.set_xlabel(f"CA dimension 1 ({inertia_frac1:.2f} inertia)", fontsize=9)
    _, _, frac2 = correspondence_axis(M)
    # compute dim2 inertia separately
    M2 = counts.to_numpy(float)
    total2 = M2.sum()
    P2 = M2 / total2
    r2 = P2.sum(axis=1); c2v = P2.sum(axis=0)
    P2 = P2[:, c2v > 0]; c2v = c2v[c2v > 0]
    S2 = np.diag(1/np.sqrt(r2)) @ (P2 - np.outer(r2, c2v)) @ np.diag(1/np.sqrt(c2v))
    _, sig2, _ = np.linalg.svd(S2, full_matrices=False)
    inerti2 = sig2**2
    frac2_val = float(inerti2[1] / inerti2.sum()) if inerti2.sum() > 0 else 0.0
    ax.set_ylabel(f"CA dimension 2 ({frac2_val:.2f} inertia)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0, color="0.85", linewidth=0.5, zorder=0)
    ax.axvline(0, color="0.85", linewidth=0.5, zorder=0)

    save(fig, "fig3_ca_ordination")
    print("fig3_ca_ordination.png written")


# ---------------------------------------------------------------------------
# F6: IDSS group structure — two-panel: group-size histogram + bridge-rank lollipop
# (Replaces the unreadable force-directed hairball of fig5_idss_network.)
# ---------------------------------------------------------------------------
def fig6_idss_structure() -> None:
    counts, coords_ll = _load_curated()
    M = counts.to_numpy(float)
    idx = list(counts.index)
    sg = seriation_groups(M, cont=CONT_PRIMARY)
    membership = sg["membership"]    # dict: row_idx -> [group_ids]
    n_groups_total = sg["n_groups"]
    n_assemblages = len(idx)

    # Group-to-assemblage mapping (invert membership)
    group_to_assems: dict[int, list[int]] = {}
    for ai, grp_list in membership.items():
        for gi in grp_list:
            group_to_assems.setdefault(gi, []).append(ai)

    # Group sizes
    group_sizes = sorted(len(v) for v in group_to_assems.values())
    max_group_size = max(group_sizes) if group_sizes else 0
    n_bridges = sum(1 for i in range(n_assemblages) if len(membership.get(i, [])) > 1)

    # Per-assemblage membership count (number of maximal groups each belongs to)
    memb_count_list: list[tuple[str, int]] = [
        (idx[i], len(membership.get(i, []))) for i in range(n_assemblages)
    ]
    memb_count_list.sort(key=lambda x: -x[1])
    top_n = 15
    top_assemblages = memb_count_list[:top_n]
    top_names = [x[0] for x in top_assemblages]
    top_vals = np.array([x[1] for x in top_assemblages], float)

    # Parkin rank and membership count
    parkin_rank = next(
        (i + 1 for i, (name, _) in enumerate(memb_count_list) if name == PARKIN_CUR),
        None,
    )
    parkin_memb = next(
        (v for name, v in memb_count_list if name == PARKIN_CUR),
        0,
    )
    median_memb = float(np.median([v for _, v in memb_count_list]))

    print(f"FIG5 IDSS STRUCTURE: cont={CONT_PRIMARY}")
    print(f"  n_groups={n_groups_total}, max_group_size={max_group_size}, "
          f"n_bridges={n_bridges}/{n_assemblages}")
    print(f"  Parkin rank={parkin_rank}, membership_count={parkin_memb}, "
          f"median_membership={median_memb:.1f}")

    # -----------------------------------------------------------------------
    # Two-panel figure
    # -----------------------------------------------------------------------
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(7, 4.5),
                                             gridspec_kw={"wspace": 0.42})

    # --- LEFT panel: histogram of co-seriable group sizes ---
    from collections import Counter
    size_counts = Counter(group_sizes)
    sizes_uniq = sorted(size_counts.keys())
    bars_left = [size_counts[s] for s in sizes_uniq]

    ax_left.bar(sizes_uniq, bars_left, color=OI_BLUE, edgecolor="white",
                linewidth=0.5, width=0.7, alpha=0.85)
    ax_left.set_xlabel("Co-seriable group size\n(assemblages per maximal group)")
    ax_left.set_ylabel("Number of groups")
    ax_left.set_xticks(sizes_uniq)
    ax_left.set_xticklabels([str(s) for s in sizes_uniq])
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)

    # Annotation: summary statistics
    ax_left.text(
        0.97, 0.97,
        f"n groups = {n_groups_total}\nmax size = {max_group_size}\n"
        f"n bridges = {n_bridges}/{n_assemblages}\ncont = {CONT_PRIMARY}",
        transform=ax_left.transAxes, fontsize=7, va="top", ha="right",
        color="0.35",
    )

    # --- RIGHT panel: horizontal lollipop of top-15 assemblages by membership count ---
    y_pos = np.arange(top_n)
    parkin_in_top = [i for i, n in enumerate(top_names) if n == PARKIN_CUR]

    colors_right = [OI_VERMIL if n == PARKIN_CUR else OI_BLUE for n in top_names]

    # Stems (horizontal lines from 0 to value)
    for i, (val, col) in enumerate(zip(top_vals, colors_right)):
        ax_right.plot([0, val], [i, i], color=col, linewidth=1.2, alpha=0.7)
    # Dots at the end
    ax_right.scatter(top_vals, y_pos, color=colors_right, s=40, zorder=4,
                     edgecolors="white", linewidths=0.4)

    # Median membership line
    ax_right.axvline(median_memb, color="0.5", linestyle="--", linewidth=0.9,
                     label=f"Median ({median_memb:.0f})")

    # Clean labels: replace underscores
    clean_names = [n.replace("_", " ") for n in top_names]
    ax_right.set_yticks(y_pos)
    ax_right.set_yticklabels(clean_names, fontsize=6.5)
    ax_right.invert_yaxis()   # rank 1 at top

    ax_right.set_xlabel("Maximal groups\nmembership count")
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.legend(frameon=False, fontsize=7, loc="lower right")

    # Highlight Parkin label in vermillion
    for tick_label, name in zip(ax_right.get_yticklabels(), top_names):
        if name == PARKIN_CUR:
            tick_label.set_color(OI_VERMIL)
            tick_label.set_fontweight("bold")

    # Annotate Parkin's rank
    if parkin_in_top:
        pi = parkin_in_top[0]
        ax_right.annotate(
            f"Parkin\n(rank {parkin_rank}, n={parkin_memb})",
            xy=(top_vals[pi], y_pos[pi]),
            xytext=(top_vals[pi] + 1.5, y_pos[pi] - 0.8),
            fontsize=6, color=OI_VERMIL,
            arrowprops=dict(arrowstyle="-", color=OI_VERMIL, lw=0.7),
        )

    # Legend proxy for color meaning
    from matplotlib.lines import Line2D
    handles_right = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=OI_VERMIL,
               markersize=7, label="Parkin"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=OI_BLUE,
               markersize=7, label="Other assemblage"),
    ]
    ax_right.legend(handles=handles_right, frameon=False, fontsize=7,
                    loc="lower right")

    save(fig, "fig6_idss_structure")
    print("fig6_idss_structure.png written")


# ---------------------------------------------------------------------------
# F5 (legacy, not called): IDSS group network — kept for reference
# ---------------------------------------------------------------------------
def fig5_idss_network() -> None:
    """Legacy network hairball — not called by main(). Kept for reference."""
    counts, coords_ll = _load_curated()
    M = counts.to_numpy(float)
    idx = list(counts.index)
    sg = seriation_groups(M, cont=CONT_PRIMARY)
    membership = sg["membership"]
    n_groups_total = sg["n_groups"]
    n_assemblages = len(idx)
    memb_count = {i: len(membership.get(i, [])) for i in range(n_assemblages)}
    group_to_assems: dict[int, list[int]] = {}
    for ai, grp_list in membership.items():
        for gi in grp_list:
            group_to_assems.setdefault(gi, []).append(ai)

    G = nx.Graph()
    for i in range(n_assemblages):
        G.add_node(idx[i], memberships=memb_count[i])
    for gi, members in group_to_assems.items():
        for ii in range(len(members)):
            for jj in range(ii + 1, len(members)):
                a_i = idx[members[ii]]
                a_j = idx[members[jj]]
                if G.has_edge(a_i, a_j):
                    G[a_i][a_j]["weight"] = G[a_i][a_j].get("weight", 1) + 1
                else:
                    G.add_edge(a_i, a_j, weight=1)

    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    nodes = list(G.nodes())
    pos = nx.spring_layout(G, seed=7, k=1.2 / (len(nodes) ** 0.5 + 1))

    fig, ax = plt.subplots(figsize=(7, 5.5))
    sizes = np.array([G.nodes[n].get("memberships", 1) for n in nodes], float)
    sizes_plot = 20 + 8 * sizes
    node_colors = [OI_VERMIL if n == PARKIN_CUR else OI_BLUE for n in nodes]
    node_sizes = [320 if n == PARKIN_CUR else float(s) for n, s in zip(nodes, sizes_plot)]
    edge_weights = np.array([G[u][v].get("weight", 1) for u, v in G.edges()], float)
    max_w = edge_weights.max() if len(edge_weights) > 0 else 1.0
    edge_alphas = np.clip(edge_weights / max_w * 0.5 + 0.1, 0.05, 0.5)
    nx.draw_networkx_edges(G, pos, ax=ax,
                           edge_color=[(0.6, 0.6, 0.6, float(a)) for a in edge_alphas],
                           width=0.6)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.85, linewidths=0.0)
    pk_memb = memb_count.get(idx.index(PARKIN_CUR), 0) if PARKIN_CUR in idx else 0
    top5 = sorted(nodes, key=lambda n: G.nodes[n].get("memberships", 0), reverse=True)[:5]
    label_dict = {n: n.replace("_", " ") for n in top5
                  if n == PARKIN_CUR or G.nodes[n].get("memberships", 0) >= pk_memb * 0.6}
    label_dict[PARKIN_CUR] = "Parkin"
    nx.draw_networkx_labels(G, pos, labels=label_dict, ax=ax, font_size=6, font_color="black")
    from matplotlib.lines import Line2D
    legend_els = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=OI_VERMIL,
               markeredgecolor="none", markersize=10, label="Parkin"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=OI_BLUE,
               markeredgecolor="none", markersize=7, label="Other assemblage"),
    ]
    ax.legend(handles=legend_els, frameon=False, fontsize=8, loc="lower right")
    ax.text(0.02, 0.97,
            f"cont={CONT_PRIMARY}; {n_groups_total} co-seriable groups; "
            f"node size by group-membership count",
            transform=ax.transAxes, fontsize=6.5, va="top", color="0.4")
    ax.axis("off")
    save(fig, "fig5_idss_network")
    print("fig5_idss_network.png written")


# ---------------------------------------------------------------------------
# F5: Empirical signature trajectory (from 07) with bootstrap CIs
# ---------------------------------------------------------------------------
def fig5_empirical_trajectory() -> None:
    counts, coords_ll = _load_curated()
    coords_df = coords_ll.dropna()
    have_coords_ids = list(coords_df.index)
    cc = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(axis=0)

    # CA ordinate
    M = counts.to_numpy(float)
    ca1, _, _ = correspondence_axis(M)

    # Orient by 14C
    rc = pd.read_csv(DATA / "raw" / "14CDatesFromMainfort2001.csv")
    rc = rc[rc["Provenience"].notna() & (rc["Provenience"] != "Provenience")].copy()
    rc["cal_mid"] = rc["Calibrated Date A.D. (1 Sigma)"].map(parse_cal_midpoint)
    prov_date = rc.groupby("Provenience")["cal_mid"].agg(["mean", "count"])
    cur_norm = {norm_name(a): a for a in counts.index}
    assem_date: dict[str, float] = {}
    for prov in prov_date.index:
        pn = norm_name(prov)
        hit = cur_norm.get(pn)
        if hit is None:
            for k, v in cur_norm.items():
                if pn and (pn in k or k in pn):
                    hit = v
                    break
        if hit is not None:
            assem_date.setdefault(hit, []).append(float(prov_date.loc[prov, "mean"]))  # type: ignore
    assem_date = {a: float(np.mean(v)) for a, v in assem_date.items()}  # type: ignore
    dated = [a for a in counts.index if a in assem_date]
    idx_map = list(counts.index)
    ca1_d = np.array([ca1[idx_map.index(a)] for a in dated])
    yr_d = np.array([assem_date[a] for a in dated])
    if len(dated) >= 3:
        rho, _ = spearmanr(ca1_d, yr_d)
        if np.isfinite(rho) and rho < 0:
            ca1 = -ca1
    ca = pd.Series(ca1, index=counts.index, name="ca")

    # Spatial clusters for F_ST and boundary excess
    sil = {k: silhouette_mean(cc_c, _kmeans_labels(cc_c, k, seed=7))
           for k in range(2, 7)}
    k_use = max(sil, key=sil.get)
    cl_labels = _kmeans_labels(cc_c, k_use, seed=7)
    cluster_of = dict(zip(have_coords_ids, cl_labels))

    ca_have = ca.reindex(have_coords_ids)
    counts_have = counts.reindex(have_coords_ids)

    def panel_for_bins(n_bins: int) -> pd.DataFrame:
        bins = pd.qcut(ca_have, q=n_bins, labels=False, duplicates="drop")
        bin_ids = sorted(pd.Series(bins).dropna().unique())
        rows = {}
        for b in bin_ids:
            ids = [i for i in have_coords_ids if not pd.isna(bins[i]) and bins[i] == b]
            sc = counts_have.loc[ids].to_numpy(float)
            sco = cc_c[[have_coords_ids.index(i) for i in ids]]
            scl = np.array([cluster_of[i] for i in ids])
            # neutral departure
            nd_vals = []
            for cc_id in np.unique(scl):
                pooled = sc[scl == cc_id].sum(axis=0)
                if pooled.sum() < 2 or (pooled > 0).sum() < 2:
                    continue
                tf, te = theta_f(pooled), theta_e(pooled)
                if np.isfinite(tf) and te > 0:
                    nd_vals.append(abs(1.0 - tf / te))
            nd = float(np.mean(nd_vals)) if nd_vals else np.nan
            # F_ST
            rep = np.unique(scl)
            if len(rep) >= 2:
                gc = np.array([sc[scl == c].sum(axis=0) for c in rep])
                fst = cultural_fst(gc)
            else:
                fst = np.nan
            # boundary excess
            sb = boundary_excess(sc, sco, seed=7) if len(ids) >= 4 else np.nan
            # seriation fragmentation: mean co-seriable-group membership per assemblage
            ser = float(np.mean([memb_count[i] for i in ids])) if ids else np.nan
            rows[b] = {"neutral_departure": nd, "seriation": ser,
                       "fst": fst, "spatial_boundary": sb}
        return pd.DataFrame(rows).T.sort_index()

    # IDSS membership count per assemblage (the seriation-fragmentation signature)
    _sg = seriation_groups(M, cont=CONT_PRIMARY)
    _memb = _sg["membership"]
    _idx = list(counts.index)
    memb_count = {_idx[i]: len(_memb.get(i, [])) for i in range(len(_idx))}

    panel = panel_for_bins(6)
    SIGS = ["neutral_departure", "seriation", "fst", "spatial_boundary"]
    SIG_LABELS = {"neutral_departure": "Neutral departure (θF/θE)",
                  "seriation": "Seriation fragmentation",
                  "fst": "Cultural F_ST",
                  "spatial_boundary": "Spatial boundary excess"}
    SIG_COLORS = {"neutral_departure": OI_BLUE, "seriation": OI_PURPLE,
                  "fst": OI_ORANGE, "spatial_boundary": OI_GREEN}

    # Bootstrap CIs on the OLS slope per signature
    rng = np.random.default_rng(11)

    def boot_slope(sig: str, n_boot: int = 600) -> tuple[float, float]:
        bslopes = []
        ids_all = list(have_coords_ids)
        for _ in range(n_boot):
            samp = list(rng.choice(ids_all, size=len(ids_all), replace=True))
            ca_s = ca.reindex(samp).reset_index(drop=True)
            cnt_s = counts.reindex(samp).reset_index(drop=True)
            cl_s = np.array([cluster_of[i] for i in samp])
            cco_s = cc_c[[have_coords_ids.index(i) for i in samp]]
            try:
                bb = pd.qcut(ca_s, q=6, labels=False, duplicates="drop")
            except Exception:
                continue
            bvals = []
            for bb_id in sorted(pd.Series(bb).dropna().unique()):
                mask = (bb == bb_id).to_numpy()
                sc = cnt_s.to_numpy(float)[mask]
                cco_b = cco_s[mask]
                cl_b = cl_s[mask]
                if sig == "neutral_departure":
                    nd_vals = []
                    for cc_id in np.unique(cl_b):
                        pooled = sc[cl_b == cc_id].sum(axis=0)
                        if pooled.sum() < 2 or (pooled > 0).sum() < 2:
                            continue
                        tf, te = theta_f(pooled), theta_e(pooled)
                        if np.isfinite(tf) and te > 0:
                            nd_vals.append(abs(1.0 - tf / te))
                    val = float(np.mean(nd_vals)) if nd_vals else np.nan
                elif sig == "fst":
                    rep = np.unique(cl_b)
                    if len(rep) >= 2:
                        gc = np.array([sc[cl_b == c].sum(axis=0) for c in rep])
                        val = cultural_fst(gc)
                    else:
                        val = np.nan
                else:
                    val = boundary_excess(sc, cco_b, seed=7) if mask.sum() >= 4 else np.nan
                bvals.append((bb_id, val))
            bvdf = pd.Series({k: v for k, v in bvals}).dropna()
            if len(bvdf) >= 2:
                bslopes.append(ols_slope(bvdf.index.to_numpy(float), bvdf.values))
        bslopes_arr = np.array([b for b in bslopes if np.isfinite(b)])
        if len(bslopes_arr) < 20:
            return np.nan, np.nan
        return float(np.percentile(bslopes_arr, 2.5)), float(np.percentile(bslopes_arr, 97.5))

    # Small multiples: one panel per signature on its native scale, so the
    # dissociation (each signature moving its own way) is read directly rather
    # than from a z-standardized overlay. Per-panel Spearman rho is annotated.
    x = panel.index.to_numpy(float)
    fig, axes = plt.subplots(2, 2, figsize=(7, 5.2))
    for ax, sig in zip(axes.ravel(), SIGS):
        y = panel[sig].to_numpy(float)
        ax.plot(x, y, marker="o", color=SIG_COLORS[sig], linewidth=1.6)
        ax.set_ylabel(SIG_LABELS[sig], fontsize=8)
        ax.set_xticks(sorted(int(v) for v in x))
        s = panel[sig].dropna()
        if len(s) >= 3:
            rho, _ = spearmanr(s.index.to_numpy(float), s.values)
            ax.text(0.95, 0.95, rf"$\rho$ = {rho:+.2f}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=9)
    for ax in axes[1, :]:
        ax.set_xlabel("CA seriation bin (early to late)", fontsize=8)
    fig.tight_layout()
    save(fig, "fig5_empirical_trajectory")
    print("fig5_empirical_trajectory.png written")


# ---------------------------------------------------------------------------
# F7: Settlement mound-height ranking (basin set, like-for-like field) + bar inset
# ---------------------------------------------------------------------------
def fig7_ranksize() -> None:
    broad_counts = load_pfg_counts(DATA / "raw" / "PFGData_sherds.csv")
    if not broad_counts.index.is_unique:
        broad_counts = broad_counts.groupby(level=0).sum()
    lmv = load_lmv(DATA / "LMVData_locations.csv")
    joined, _ = join_pfg_to_lmv(broad_counts, lmv)
    bmatched = joined.dropna(subset=["Easting", "Northing"]).copy()

    # Restrict the broad settlement set to the St. Francis basin (drainage
    # definition; canonical list from analyses/16_basin_membership.py).
    broad_members = _basin_members("broad")
    bmatched = bmatched[[str(i) in broad_members for i in bmatched.index]].copy()

    lmv2 = pd.read_csv(DATA / "LMVData-22March2006.csv")
    lmv2 = lmv2.dropna(subset=["Number"]).copy()
    lmv2["_k"] = lmv2["Number"].astype(str).map(normalize_grid)
    lmv2 = lmv2.drop_duplicates(subset=["_k"], keep="first").set_index("_k")
    norm_ids = pd.Index([normalize_grid(str(i)) for i in bmatched.index])
    ext = lmv2.reindex(norm_ids)
    ext.index = bmatched.index

    def to_bin(series: pd.Series) -> pd.Series:
        v = pd.to_numeric(
            series.astype(str).str.replace("?", "", regex=False), errors="coerce"
        )
        return (v.fillna(0) > 0)

    mound = to_bin(ext["Mound"])
    ditch = to_bin(ext["Ditch"])
    stfr = to_bin(ext["St Francis"])
    # Mound HEIGHT is the like-for-like settlement-scale field: Parkin's height
    # (23 ft) is recorded, whereas its Max Mound Area is under-coded as 0, so the
    # height ranking needs no substituted value (the earlier area rank-size, which
    # required substituting Parkin's ~17-acre site area into a mound-area field,
    # was not a like-for-like comparison and is not used).
    mound_ht = pd.to_numeric(ext["Max Mound Height (ft)"], errors="coerce")
    if PARKIN_BROAD in bmatched.index:
        mound.loc[PARKIN_BROAD] = True
        ditch.loc[PARKIN_BROAD] = True
        stfr.loc[PARKIN_BROAD] = True

    ht = mound_ht.dropna()
    ht = ht[ht > 0].sort_values(ascending=False)
    ranks = np.arange(1, len(ht) + 1)
    parkin_ht_rank = (
        int((ht.index == PARKIN_BROAD).argmax() + 1)
        if PARKIN_BROAD in ht.index else None
    )
    primacy_ht = float(ht.iloc[0] / ht.iloc[1]) if len(ht) >= 2 else float("nan")

    # Main panel: within-basin mound-height ranking (rank vs height)
    fig = plt.figure(figsize=(7, 4.5))
    ax_main = fig.add_axes([0.10, 0.14, 0.62, 0.80])
    ax_bar = fig.add_axes([0.77, 0.18, 0.20, 0.72])

    ax_main.plot(ranks, ht.values, "o", ms=4, color=OI_SKY, alpha=0.75,
                 markeredgecolor="none")
    if parkin_ht_rank:
        ax_main.plot(parkin_ht_rank, ht.loc[PARKIN_BROAD], "*", ms=15,
                     color=OI_BLACK, zorder=5,
                     label=f"Parkin ({ht.loc[PARKIN_BROAD]:.0f} ft, rank {parkin_ht_rank}/{len(ht)})")
    ax_main.set_xlabel("Rank")
    ax_main.set_ylabel("Maximum mound height (ft)")
    ax_main.text(0.96, 0.78,
                 f"graded distribution\ntallest / second = {primacy_ht:.2f}",
                 transform=ax_main.transAxes, fontsize=7, va="top", ha="right",
                 color="0.35")
    ax_main.legend(frameon=False, fontsize=8, loc="upper right")
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)

    # Bar inset: mound-present / ditch / St-Francis proportions
    n_tot = len(bmatched)
    bar_vals = [float(mound.sum()) / n_tot,
                float(ditch.sum()) / n_tot,
                float(stfr.sum()) / n_tot]
    bar_labels = ["Mound\npresent", "Ditch\n(defensive)", "St Francis\n(fortified)"]
    bar_colors = [OI_BLUE, OI_ORANGE, OI_GREEN]
    xpos = np.arange(3)
    ax_bar.bar(xpos, bar_vals, color=bar_colors, width=0.6, edgecolor="none")
    ax_bar.set_xticks(xpos)
    ax_bar.set_xticklabels(bar_labels, fontsize=6.5)
    ax_bar.set_ylabel("Proportion\nof basin set", fontsize=7)
    ax_bar.set_ylim(0, 1.0)
    ax_bar.tick_params(labelsize=6.5)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.axhline(0, color="0.8", linewidth=0.5)

    save(fig, "fig7_ranksize")
    print("fig7_ranksize.png written")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:

    print("Generating F3 (CA ordination)...")
    fig3_ca_ordination()

    print("Generating S6 (idealized validation)...")
    fig4_validation()  # saves figS4_validation (the idealized-data validation; moved to Supplement)

    # Fig 4 (record-matched recovery) and Fig 5 (size-controlled empirical trajectory) are
    # generated by analyses/21_signal_recovery.py, which owns the rarefaction machinery.
    # fig5_empirical_trajectory() below produces the RAW (uncontrolled) version and is retired.

    print("Generating F6 (IDSS group structure)...")
    fig6_idss_structure()

    print("Generating F7 (rank-size + inset)...")
    fig7_ranksize()


    print("\nAll figures written to figures/.")


if __name__ == "__main__":
    main()
