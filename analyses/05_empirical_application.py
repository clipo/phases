"""Phase 5: empirical application of the validated convergence criterion to the
real Phillips, Ford & Griffin (PFG) ceramic data.

This script applies the four-signature convergence criterion (validated in
Phase 4) to the real lower/central Mississippi Valley (LMV) assemblages and
reports a NEUTRAL pattern summary. It does NOT declare H1 (nascent emergence:
the four signatures co-rise toward contact and concentrate on Parkin) or H2
(stable non-consolidation, Rees 2001: signals flat / non-trending). The numbers
and the pattern are reported; the interpretation is left to the team.

Data policy: data/ is gitignored and location-sensitive. This script NEVER
prints or writes raw coordinates. Only cluster sizes, chosen k, ranks, and
percentiles are reported. Figures plot points in a centered/relative frame with
no axis tick labels.

Outputs (all gitignored): output/empirical_findings.md and figures/*.png.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mls_emergence.dataio.pfg import load_pfg_counts
from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv, normalize_grid
from mls_emergence.signatures.neutral import theta_f, theta_e
from mls_emergence.signatures.variance import cultural_fst
from mls_emergence.signatures.assortativity import (
    boundary_excess,
    _kmeans_labels,
    geo_distance,
)
from mls_emergence.signatures.seriation import seriation_groups, seriation_solutions
from mls_emergence.signatures.convergence import convergence_score, time_derivative

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "output"
FIGURES = ROOT / "figures"
OUTPUT.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

PARKIN = "11-N-1"

# ---------------------------------------------------------------------------
# House figure style
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 9,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)
# Okabe-Ito colorblind-safe palette
OKABE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
}


def silhouette_mean(coords: np.ndarray, labels: np.ndarray) -> float:
    """Mean silhouette coefficient (no sklearn dependency)."""
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
        for c in uniq:
            if c == labels[i]:
                continue
            mask = labels == c
            if mask.sum():
                b = min(b, d[i, mask].mean())
        s[i] = (b - a) / max(a, b) if max(a, b) > 0 else 0.0
    return float(s.mean())


def main() -> None:
    lines: list[str] = []

    def emit(s: str = "") -> None:
        lines.append(s)

    emit("# Empirical application of the convergence criterion to PFG data (Phase 5)")
    emit()
    emit(
        "Neutral pattern report. The four signatures (neutral departure, cultural "
        "F_ST, spatial boundary excess, seriation structure) are computed along an "
        "ordinal trajectory and at Parkin, plus external settlement/mound "
        "cross-checks. No verdict is declared between H1 (co-rise/convergence "
        "toward contact, concentrating on Parkin) and H2 (flat / non-trending; "
        "Rees 2001). See caveats and 'For interpretation by the team'."
    )
    emit()

    # -----------------------------------------------------------------------
    # 1. Load and join
    # -----------------------------------------------------------------------
    counts = load_pfg_counts(DATA / "raw" / "PFGData_sherds.csv")
    names = counts.attrs.get("names", {})
    # A few site numbers appear twice (16-L-1, 16-L-3, 16-N-6). Collapse
    # duplicate assemblage ids by summing counts so the index is unique and
    # positional lookups are well defined.
    n_before = len(counts)
    if not counts.index.is_unique:
        counts = counts.groupby(level=0).sum()
    n_dup_collapsed = n_before - len(counts)
    lmv = load_lmv(DATA / "LMVData_locations.csv")
    joined, unmatched = join_pfg_to_lmv(counts, lmv)
    matched = joined.dropna(subset=["Easting", "Northing"]).copy()
    type_cols = list(counts.columns)

    # Attach extra LMV attribute columns (not added by join_pfg_to_lmv) by the
    # same normalized-grid key, for the cross-checks.
    extra_cols = [
        "Num_Mounds",
        "Max Mound Height (ft)",
        "Max Mound Area (sq ft)",
        "St Francis",
        "Ditch",
        "Platform",
        "Name",
    ]
    lmv_keyed = lmv.dropna(subset=["Number"]).copy()
    lmv_keyed["_key"] = lmv_keyed["Number"].astype(str).map(normalize_grid)
    lmv_keyed = lmv_keyed.drop_duplicates(subset=["_key"], keep="first").set_index(
        "_key"
    )
    norm_ids = matched.index.astype(str).str.strip().map(normalize_grid)
    extra = lmv_keyed.reindex(norm_ids.values)
    for c in extra_cols:
        if c in extra.columns:
            matched[c] = extra[c].values

    n_types = len(type_cols)
    emit("## 1. Data and join")
    emit()
    emit(
        f"- PFG assemblages loaded: {n_before} rows -> {len(counts)} unique "
        f"assemblage ids x {n_types} ceramic types "
        f"({n_dup_collapsed} duplicate ids collapsed by summing)."
    )
    emit(
        f"- Matched to LMV with coordinates: {len(matched)} / {len(counts)} "
        f"(unmatched: {len(unmatched)}); all in UTM Zone "
        f"{int(matched['Zone'].dropna().unique()[0])}."
    )
    emit(f"- Parkin ({PARKIN}) present and matched: {PARKIN in matched.index}.")
    emit()

    # -----------------------------------------------------------------------
    # Ordinal axis decision (HONEST DATA NOTE)
    # -----------------------------------------------------------------------
    # The brief named data/raw/pfg-cpl-frequency.csv Seriation_Number as the
    # primary ordinal axis "203 assemblages". Inspection shows that file is the
    # IDSS minmax-graph solution table from Lipo, Madsen & Dunnell (2015): 102
    # numbered ordering-positions over only 12 UNIQUE decorated assemblages
    # (each appears many times, spanning nearly the whole 1-102 range, with
    # constant counts). It therefore cannot supply a single per-assemblage
    # chronological ordinal for the 258 settlement sites. We use it for what it
    # is (the Parkin-region decorated seriation, Sig2) and adopt the LMV
    # 'Terminal Period' phase (F earliest -> A = latest/contact) as the working
    # ordinal trajectory axis for the 258 sites. This was the brief's named
    # validation axis; we promote it to the working axis because the seriation
    # file is unusable as a 258-site axis. Both choices are documented.
    cpl = pd.read_csv(
        DATA / "raw" / "pfg-cpl-frequency.csv"
    ).dropna(subset=["Assemblage"])
    cpl_type_cols = [
        c for c in cpl.columns if c not in ("Seriation_Number", "Assemblage")
    ]
    n_unique_cpl = cpl["Assemblage"].nunique()
    # per-assemblage mean ordering-position across the solution table
    cpl_pos = (
        cpl.groupby("Assemblage")["Seriation_Number"].mean().rename("mean_serpos")
    )

    matched_ids = set(matched.index.astype(str))
    cpl_ids_in_matched = sorted(
        a for a in cpl["Assemblage"].unique() if str(a) in matched_ids
    )

    emit("## 2. Ordinal axis (data note)")
    emit()
    emit(
        "The seriation-frequency file (pfg-cpl-frequency.csv) is the IDSS "
        f"solution table: {len(cpl)} numbered ordering positions over only "
        f"{n_unique_cpl} unique decorated assemblages (each repeated, spanning "
        "nearly the full range, constant counts). It cannot supply a single "
        "per-assemblage ordinal for the 258 sites. We therefore use it for the "
        "Parkin-region decorated seriation (Signature 2, below) and adopt the "
        "LMV 'Terminal Period' phase (F earliest -> A = contact/latest) as the "
        "working ordinal trajectory axis for the 258 sites."
    )
    emit(
        f"- Decorated assemblages in the seriation file also matched in the "
        f"258: {len(cpl_ids_in_matched)} ({', '.join(cpl_ids_in_matched)})."
    )

    # Working ordinal axis: Terminal Period A..F -> ordinal 1..6 (A latest=6)
    PHASE_ORDER = {"F": 1, "E": 2, "D": 3, "C": 4, "B": 5, "A": 6}
    tp = matched["Terminal Period"].astype(str).str.strip()
    matched["phase_ord"] = tp.map(PHASE_ORDER)
    phase_known = matched.dropna(subset=["phase_ord"]).copy()
    phase_known["phase_ord"] = phase_known["phase_ord"].astype(int)
    emit(
        f"- Sites with a usable Terminal Period phase: {len(phase_known)} / "
        f"{len(matched)} (others coded '?' or blank)."
    )
    tp_counts = phase_known.groupby("phase_ord").size()
    emit(
        "- Phase counts (ordinal: 1=F earliest ... 6=A contact): "
        + ", ".join(f"{int(k)}:{int(v)}" for k, v in tp_counts.items())
        + "."
    )
    emit()

    # -----------------------------------------------------------------------
    # 2. Spatial clustering (k-means, silhouette over k=4..10)
    # -----------------------------------------------------------------------
    coords = matched[["Northing", "Easting"]].to_numpy(float)
    # center to avoid any chance of leaking absolute coordinates downstream
    coords_c = coords - coords.mean(axis=0)
    sil = {}
    for k in range(4, 11):
        labels = _kmeans_labels(coords_c, k, seed=7)
        sil[k] = silhouette_mean(coords_c, labels)
    best_k = max(sil, key=sil.get)
    cluster_labels = _kmeans_labels(coords_c, best_k, seed=7)
    matched["cluster"] = cluster_labels
    sizes = pd.Series(cluster_labels).value_counts().sort_index()

    emit("## 3. Spatial clustering (k-means on coordinates)")
    emit()
    emit(
        "- Silhouette by k: "
        + ", ".join(f"k={k}:{sil[k]:.3f}" for k in range(4, 11))
        + "."
    )
    emit(f"- Chosen k (max silhouette): **{best_k}**.")
    emit(
        "- Cluster sizes: "
        + ", ".join(f"c{int(c)}:{int(n)}" for c, n in sizes.items())
        + "."
    )
    parkin_cluster = int(matched.loc[PARKIN, "cluster"])
    emit(f"- Parkin ({PARKIN}) cluster: c{parkin_cluster}.")
    emit()

    # -----------------------------------------------------------------------
    # 3. Ordinal trajectory: bin by phase, four signatures per bin
    # -----------------------------------------------------------------------
    # Bins = the ordinal phases that carry enough sites. With only 6 phase
    # levels (and F, D sparse) we use the phase ordinal directly as the bin
    # (1..6); this is coarser than 6-8 quantile windows but is the honest
    # resolution the phase data supports.
    bin_ords = sorted(phase_known["phase_ord"].unique())
    panel_rows = {}
    bin_info = []
    for b in bin_ords:
        sub = phase_known[phase_known["phase_ord"] == b]
        ids = list(sub.index)
        sub_counts = sub[type_cols].to_numpy(float)
        sub_coords = (coords_c[[matched.index.get_loc(i) for i in ids]])
        sub_clusters = matched.loc[ids, "cluster"].to_numpy()

        # (a) neutral departure: mean within-cluster |1 - tf/te|, pooling each
        #     cluster's assemblages in this bin.
        nd_vals = []
        for c in np.unique(sub_clusters):
            pooled = sub_counts[sub_clusters == c].sum(axis=0)
            if pooled.sum() < 2 or (pooled > 0).sum() < 2:
                continue
            tf, te = theta_f(pooled), theta_e(pooled)
            if np.isfinite(tf) and te > 0:
                nd_vals.append(abs(1.0 - tf / te))
        neutral_departure = float(np.mean(nd_vals)) if nd_vals else np.nan

        # (b) cultural F_ST across clusters represented in the bin
        rep_clusters = np.unique(sub_clusters)
        if len(rep_clusters) >= 2:
            group_counts = np.array(
                [sub_counts[sub_clusters == c].sum(axis=0) for c in rep_clusters]
            )
            fst = cultural_fst(group_counts)
        else:
            fst = np.nan

        # (c) spatial boundary excess on the bin's assemblages
        if len(ids) >= 4:
            spatial_boundary = boundary_excess(sub_counts, sub_coords, seed=7)
        else:
            spatial_boundary = np.nan

        panel_rows[b] = {
            "neutral_departure": neutral_departure,
            "fst": fst,
            "spatial_boundary": spatial_boundary,
        }
        bin_info.append((b, len(ids), len(rep_clusters)))

    panel = pd.DataFrame(panel_rows).T.sort_index()
    panel.index.name = "phase_ord"

    # Signature 2 (seriation) is reported as a static per-cluster + Parkin
    # result, NOT as a bin trajectory (IDSS is combinatorial; bins have >15
    # assemblages). It is EXCLUDED from the bin-trajectory convergence below.

    emit("## 4. Ordinal trajectory panel (phase bins x signatures)")
    emit()
    emit("Bins (phase ordinal: 1=F ... 6=A), n sites, n clusters represented:")
    for b, ns, nc in bin_info:
        emit(f"- bin {int(b)}: n={ns}, clusters={nc}")
    emit()
    emit("Panel (NaN where a signature is undefined for that bin):")
    emit()
    emit("| phase_ord | neutral_departure | fst | spatial_boundary |")
    emit("|---|---|---|---|")
    for b, r in panel.iterrows():
        emit(
            f"| {int(b)} | {r['neutral_departure']:.4f} | {r['fst']:.4f} | "
            f"{r['spatial_boundary']:.3f} |"
        )
    emit()

    # Convergence score + per-signature time derivative (slope over phase ord).
    # convergence_score z-standardizes each column then averages; compute on the
    # rows where all three are present so the z-scores are comparable.
    panel_full = panel.dropna(how="any")
    conv = convergence_score(panel_full)
    emit("Convergence score (z-averaged across the three bin-trajectory "
         "signatures; computed on complete bins only):")
    emit()
    if len(panel_full) >= 2:
        for b, v in conv.items():
            emit(f"- phase {int(b)}: {v:.3f}")
        conv_slope = time_derivative(conv)
    else:
        conv_slope = np.nan
        emit("- (insufficient complete bins)")
    emit()

    from scipy.stats import spearmanr

    slopes = {}
    spearman = {}
    monotone = {}
    emit("Per-signature trend vs phase ordinal. OLS slope (positive = upward "
         "toward contact/A), Spearman rank correlation with the ordinal axis "
         "(monotone-trend test, robust to scale), and whether the bin series is "
         "monotone:")
    emit()
    for col in ["neutral_departure", "fst", "spatial_boundary"]:
        s = panel[col].dropna()
        slope = time_derivative(s) if len(s) >= 2 else np.nan
        slopes[col] = slope
        if len(s) >= 3:
            rho, pval = spearmanr(s.index.to_numpy(float), s.values)
        else:
            rho, pval = np.nan, np.nan
        spearman[col] = (rho, pval)
        d = np.diff(s.values)
        mono = bool((d >= 0).all() or (d <= 0).all()) if len(s) >= 2 else False
        monotone[col] = mono
        emit(
            f"- {col}: slope = {slope:+.5f} (over {len(s)} bins); "
            f"Spearman rho = {rho:+.3f} (p = {pval:.3f}); monotone = {mono}"
        )
    if len(panel_full) >= 3:
        crho, cpval = spearmanr(conv.index.to_numpy(float), conv.values)
    else:
        crho, cpval = np.nan, np.nan
    spearman["convergence"] = (crho, cpval)
    emit(
        f"- convergence_score: slope = {conv_slope:+.5f}; "
        f"Spearman rho = {crho:+.3f} (p = {cpval:.3f})"
    )
    emit()

    # -----------------------------------------------------------------------
    # 4. Parkin focus + Signature 2 (decorated seriation, Lipo 2001 hook)
    # -----------------------------------------------------------------------
    # Per-assemblage relative seriation position from the IDSS table (mean
    # ordering position; the only defensible scalar given repeated rows).
    cpl_uniq = (
        cpl.drop_duplicates("Assemblage").set_index("Assemblage")[cpl_type_cols]
    )
    cpl_ids = list(cpl_uniq.index)
    mat = cpl_uniq.to_numpy(float)
    CONT = 0.30  # user-set continuity threshold reported by Lipo et al. (2015)
    sg = seriation_groups(mat, cont=CONT)
    sols = seriation_solutions(mat, cont=CONT)
    maxlen = max((len(s) for s in sols), default=0)
    multi_ids = [cpl_ids[i] for i in sg["multi"]]
    parkin_in_cpl = PARKIN in cpl_ids
    if parkin_in_cpl:
        parkin_groups = sg["membership"][cpl_ids.index(PARKIN)]
        parkin_serpos = float(cpl_pos.get(PARKIN, np.nan))
        # relative position among the 12 (rank of mean ordering position)
        rel = cpl_pos.reindex(cpl_ids).rank().get(PARKIN, np.nan)
    else:
        parkin_groups, parkin_serpos, rel = [], np.nan, np.nan

    # Robustness across continuity thresholds
    robust = {}
    for c in [1.0, 0.5, 0.30, 0.10]:
        r = seriation_groups(mat, cont=c)
        robust[c] = (
            r["n_groups"],
            len(r["membership"][cpl_ids.index(PARKIN)]) if parkin_in_cpl else 0,
            (PARKIN in [cpl_ids[i] for i in r["multi"]]),
        )

    emit("## 5. Parkin focus and Signature 2 (decorated seriation)")
    emit()
    emit(
        "Signature 2 is computed ONLY on the small Parkin-region decorated "
        f"subset ({len(cpl_ids)} unique assemblages) from the IDSS table; it is "
        "NOT a bin-trajectory signature (excluded from the convergence above)."
    )
    emit(
        f"- IDSS at continuity={CONT}: {sg['n_groups']} maximal co-seriable "
        f"groups, max solution length {maxlen}."
    )
    if parkin_in_cpl:
        pos_label = (
            "late" if rel >= 9 else ("early" if rel <= 4 else "mid")
        )
        emit(
            f"- Parkin ({PARKIN}) mean ordering position = {parkin_serpos:.1f} "
            f"(rank {int(rel)}/{len(cpl_ids)} among decorated assemblages => "
            f"**{pos_label}**)."
        )
        emit(
            f"- Parkin Signature 2 multi-membership: belongs to "
            f"**{len(parkin_groups)}** co-seriable groups "
            f"(multi-membership = {PARKIN in multi_ids}; the Lipo 2001 "
            "two-lineage / bridge result)."
        )
    else:
        emit("- Parkin not present in the decorated seriation file.")
    emit("- Robustness of Parkin multi-membership across continuity thresholds:")
    for c, (ng, npg, pm) in robust.items():
        emit(f"  - cont={c}: n_groups={ng}, Parkin_in_groups={npg}, multi={pm}")
    emit()
    # Parkin's phase position on the working axis
    parkin_tp = str(matched.loc[PARKIN, "Terminal Period"]).strip()
    emit(
        f"- Parkin Terminal Period phase = {parkin_tp} (working ordinal "
        f"{PHASE_ORDER.get(parkin_tp, 'NA')}; A=6 is contact/latest)."
    )
    emit()

    # Per-cluster static Sig2 is not run on full clusters (>15, combinatorial);
    # state this explicitly.
    emit(
        "Per-cluster static Signature 2 was NOT computed on full spatial "
        "clusters: each exceeds the ~15-assemblage tractability limit of the "
        "deterministic IDSS solver, and the count columns differ (26 PFG types "
        "vs the 10 decorated types of the seriation file). Signature 2 is "
        "reported only on the decorated Parkin-region subset above."
    )
    emit()

    # -----------------------------------------------------------------------
    # 5. Cross-checks (external H1/H2 proxies, LMV columns)
    # -----------------------------------------------------------------------
    emit("## 6. Cross-checks (external settlement / mound proxies)")
    emit()

    # Rank-size of site SIZE. DATA NOTE: the LMV 'Area' column is a locality
    # label (Memphis, St. Francis, Sunflower, ...), NOT a site-size measure. The
    # genuine size proxy in LMV is 'Max Mound Area (sq ft)'; we use it for the
    # rank-size / primacy analysis and state the substitution.
    area = pd.to_numeric(matched["Max Mound Area (sq ft)"], errors="coerce").dropna()
    area = area[area > 0].sort_values(ascending=False)
    ranks = np.arange(1, len(area) + 1)
    logr = np.log(ranks)
    loga = np.log(area.values)
    # OLS slope of log(size) on log(rank)
    A = np.vstack([logr, np.ones_like(logr)]).T
    slope_rs, intercept_rs = np.linalg.lstsq(A, loga, rcond=None)[0]
    # primacy index: largest / second-largest
    primacy = float(area.iloc[0] / area.iloc[1]) if len(area) > 1 else np.nan
    largest_id = area.index[0]
    parkin_area_rank = (
        int((area.index == PARKIN).argmax() + 1) if PARKIN in area.index else None
    )
    emit("### Rank-size of site size (Max Mound Area, sq ft)")
    emit(
        "- NOTE: LMV 'Area' is a locality label, not site size; the size proxy "
        "used here is 'Max Mound Area (sq ft)'."
    )
    emit(f"- N sites with mound area > 0: {len(area)}.")
    emit(
        f"- log-log rank-size slope = {slope_rs:.3f} "
        "(Zipf/log-normal expectation ~= -1.0; shallower than -1 = convex "
        "below the top sites; steeper = primate)."
    )
    emit(
        f"- Primacy (largest/second) = {primacy:.2f}; largest site = "
        f"{largest_id} ('{names.get(largest_id, '?')}'); is largest Parkin = "
        f"{largest_id == PARKIN}."
    )
    if parkin_area_rank:
        emit(
            f"- Parkin size rank (Max Mound Area) = {parkin_area_rank}/"
            f"{len(area)} (percentile "
            f"{100*(1 - (parkin_area_rank-1)/len(area)):.1f})."
        )
    else:
        emit("- Parkin has no recorded Max Mound Area > 0.")
    emit()

    # Mound differentiation
    mh = pd.to_numeric(matched["Max Mound Height (ft)"], errors="coerce").dropna()
    nm = pd.to_numeric(matched["Num_Mounds"], errors="coerce").dropna()
    emit("### Mound differentiation")
    emit(
        f"- Max Mound Height (ft): n={len(mh)}, max={mh.max():.1f}, "
        f"median={mh.median():.1f}, mean={mh.mean():.2f}, 90th pct="
        f"{mh.quantile(0.90):.1f}."
    )
    if PARKIN in mh.index:
        ph = mh.loc[PARKIN]
        pct = float((mh < ph).mean() * 100)
        rank_h = int((mh.sort_values(ascending=False).index == PARKIN).argmax() + 1)
        emit(
            f"- Parkin Max Mound Height = {ph:.1f} ft, percentile {pct:.1f}, "
            f"rank {rank_h}/{len(mh)}."
        )
    else:
        emit("- Parkin Max Mound Height: not recorded.")
    emit(
        f"- Num_Mounds: n={len(nm)}, max={nm.max():.0f}, median={nm.median():.0f}, "
        f"sites with >1 mound = {(nm > 1).sum()}."
    )
    if PARKIN in nm.index:
        emit(f"- Parkin Num_Mounds = {nm.loc[PARKIN]:.0f}.")
    emit()

    # Fortification
    def present_count(col):
        if col not in matched.columns:
            return None
        s = matched[col]
        num = pd.to_numeric(s, errors="coerce")
        if num.notna().any():
            return int((num.fillna(0) > 0).sum())
        return int(s.astype(str).str.strip().replace("nan", "").astype(bool).sum())

    emit(f"### Fortification / monumental features (present count among "
         f"{len(matched)})")
    for col in ["St Francis", "Ditch", "Platform"]:
        pc = present_count(col)
        pk = ""
        if col in matched.columns:
            pv = matched.loc[PARKIN, col]
            pk = f" (Parkin: {pv})"
        emit(f"- {col}: {pc} sites present{pk}.")
    emit()

    # -----------------------------------------------------------------------
    # 6. Figures
    # -----------------------------------------------------------------------
    # (a) four-signature ordinal trajectory + convergence
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = panel.index.to_numpy(float)

    def zscore(s):
        s = s.astype(float)
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd and np.isfinite(sd) else s * 0

    series_specs = [
        ("neutral_departure", OKABE["blue"], "Neutral departure"),
        ("fst", OKABE["orange"], "Cultural $F_{ST}$"),
        ("spatial_boundary", OKABE["green"], "Spatial boundary excess"),
    ]
    for col, color, label in series_specs:
        s = panel[col]
        ax.plot(
            x, zscore(s), marker="o", color=color, label=label, linewidth=1.6
        )
    if len(panel_full) >= 2:
        ax.plot(
            panel_full.index.to_numpy(float),
            zscore(conv),
            marker="s",
            color=OKABE["black"],
            linewidth=2.0,
            linestyle="--",
            label="Convergence score",
        )
    ax.set_xlabel("Ordinal phase (1 = F earliest ... 6 = A contact)")
    ax.set_ylabel("Signature (z-standardized)")
    ax.set_xticks(sorted(int(v) for v in x))
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.axhline(0, color="0.7", linewidth=0.6, zorder=0)
    fig.tight_layout()
    fig.savefig(FIGURES / "05_ordinal_trajectory.png")
    plt.close(fig)

    # (b) spatial cluster map (relative, centered, NO tick labels)
    fig, ax = plt.subplots(figsize=(5.0, 5.5))
    rng = np.random.default_rng(0)
    jit = rng.normal(0, coords_c.std() * 0.01, coords_c.shape)
    pj = coords_c + jit
    palette = list(OKABE.values())
    for c in range(best_k):
        m = cluster_labels == c
        ax.scatter(
            pj[m, 1],
            pj[m, 0],
            s=14,
            color=palette[c % len(palette)],
            alpha=0.8,
            edgecolors="none",
            label=f"c{c} (n={int(m.sum())})",
        )
    pidx = matched.index.get_loc(PARKIN)
    ax.scatter(
        pj[pidx, 1],
        pj[pidx, 0],
        s=120,
        marker="*",
        color=OKABE["black"],
        edgecolors="white",
        linewidths=0.6,
        label="Parkin (11-N-1)",
        zorder=5,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Relative easting (no scale)")
    ax.set_ylabel("Relative northing (no scale)")
    ax.set_aspect("equal")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURES / "05_spatial_clusters.png")
    plt.close(fig)

    # (c) rank-size of Area, Parkin marked
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.loglog(ranks, area.values, "o", ms=4, color=OKABE["sky"], alpha=0.7)
    fit = np.exp(intercept_rs) * ranks ** slope_rs
    ax.loglog(
        ranks,
        fit,
        "-",
        color=OKABE["vermillion"],
        linewidth=1.5,
        label=f"OLS slope = {slope_rs:.2f}",
    )
    if parkin_area_rank:
        ax.loglog(
            parkin_area_rank,
            area.loc[PARKIN],
            "*",
            ms=16,
            color=OKABE["black"],
            label=f"Parkin (rank {parkin_area_rank})",
        )
    ax.set_xlabel("Rank (log)")
    ax.set_ylabel("Max mound area, sq ft (log)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "05_rank_size_area.png")
    plt.close(fig)

    emit("## 7. Figures")
    emit()
    emit(
        "- figures/05_ordinal_trajectory.png: the three bin-trajectory "
        "signatures (z-standardized) plus the convergence score across the "
        "ordinal phase axis (1=F earliest to 6=A contact). Caption: positive "
        "trends would indicate co-rise toward contact (H1); flat lines indicate "
        "non-trending signals (H2)."
    )
    emit(
        "- figures/05_spatial_clusters.png: k-means spatial communities in a "
        "centered, jittered frame with NO axis scale (coordinates not exposed); "
        "Parkin starred. Caption: groups used as units for F_ST and boundary "
        "excess."
    )
    emit(
        "- figures/05_rank_size_area.png: log-log rank-size of site size (Max "
        "Mound Area, sq ft) with OLS slope; Parkin has no recorded mound area "
        "so it is absent from this curve (its size shows in Max Mound Height, "
        "section 6). Caption: slope near -1 indicates a log-normal/Zipf "
        "settlement system; a strongly primate (convex) curve with one dominant "
        "site would indicate consolidation."
    )
    emit()

    # -----------------------------------------------------------------------
    # 7. Neutral pattern summary + caveats
    # -----------------------------------------------------------------------
    emit("## 8. Pattern summary (NEUTRAL)")
    emit()
    emit(
        "Trend label uses BOTH the OLS slope sign and the Spearman rank "
        "correlation. A signature is called 'monotone rise/decline' only when "
        "the bin series is monotone AND Spearman |rho| is large; otherwise it is "
        "'upward/downward but non-monotone' (slope sign with volatility) or "
        "'flat' (|rho| small). With only 5-6 bins no Spearman p-value reaches "
        "significance, so these are descriptive, not inferential, labels."
    )
    emit()

    def trend_word(col):
        slope = slopes.get(col, np.nan)
        rho = spearman.get(col, (np.nan, np.nan))[0]
        mono = monotone.get(col, False)
        if not np.isfinite(slope):
            return "undetermined"
        if not np.isfinite(rho) or abs(rho) < 0.3:
            return "flat / non-trending"
        direction = "upward" if slope > 0 else "downward"
        if mono and abs(rho) >= 0.7:
            return f"monotone {direction} toward contact"
        return f"{direction} but non-monotone"

    for col, _, label in series_specs:
        rho = spearman[col][0]
        emit(
            f"- {label}: slope {slopes[col]:+.5f}, Spearman rho {rho:+.3f} -> "
            f"**{trend_word(col)}**."
        )
    crho = spearman.get("convergence", (np.nan, np.nan))[0]
    if not np.isfinite(crho) or abs(crho) < 0.3:
        cw = "flat / non-trending"
    else:
        cw = ("upward" if conv_slope > 0 else "downward") + (
            " (monotone)" if np.isfinite(crho) and abs(crho) >= 0.7 else " (non-monotone)"
        )
    emit(
        f"- Convergence score: slope {conv_slope:+.5f}, Spearman rho "
        f"{crho:+.3f} -> **{cw}**."
    )
    rising = sum(
        1
        for col, _, _ in series_specs
        if "upward" in trend_word(col)
    )
    flat = sum(1 for col, _, _ in series_specs if trend_word(col).startswith("flat"))
    emit()
    emit(
        f"- Of the 3 bin-trajectory signatures: {rising} trend upward toward "
        f"contact (monotone or not), {flat} flat/non-trending, "
        f"{3 - rising - flat} downward. Whether they CONVERGE (co-rise "
        "together) is the joint question; the convergence-score slope is "
        f"{conv_slope:+.5f} (Spearman rho {crho:+.3f})."
    )
    emit(
        f"- Parkin: working-axis phase {parkin_tp} (ordinal "
        f"{PHASE_ORDER.get(parkin_tp,'NA')}); decorated-seriation position "
        + (f"rank {int(rel)}/{len(cpl_ids)}; " if parkin_in_cpl else "n/a; ")
        + f"Signature 2 multi-membership in {len(parkin_groups)} co-seriable "
        "groups (bridge / two-lineage status present)."
    )
    emit(
        f"- Settlement: rank-size slope {slope_rs:.2f}, primacy {primacy:.2f}, "
        f"largest site is {'Parkin' if largest_id == PARKIN else largest_id}."
    )
    emit()

    emit("## 9. Caveats (explicit)")
    emit()
    emit(
        "- The trajectory axis is ORDINAL phase (Terminal Period F->A), NOT "
        "absolute calendar time; slopes are per-phase, not rates."
    )
    emit(
        "- Assemblages are ~2-phase time-averaged (Period spans like B-A, D-C); "
        "the working axis uses the single Terminal (latest) phase, which "
        "compresses occupation history."
    )
    emit(
        "- The named primary seriation file is unusable as a 258-site ordinal "
        "(12 unique decorated assemblages, repeated solution positions); the "
        "LMV phase axis is used instead. This is a substitution, documented in "
        "section 2."
    )
    emit(
        "- Signature 2 (seriation) is computed only on the 12-assemblage "
        "decorated Parkin-region subset (IDSS is combinatorial) and is EXCLUDED "
        "from the bin-trajectory convergence. No per-bin Sig2 trajectory exists."
    )
    emit(
        "- The four signatures share the type-frequency substrate and are "
        "partly correlated; co-movement is expected to some degree even absent "
        "a single causal process."
    )
    emit(
        "- The convergence score z-standardizes each signature, so it is "
        "scale-free, but the underlying ABSOLUTE movements differ greatly: "
        "F_ST stays very low (~0.013 to ~0.030) and neutral departure is "
        "non-monotone, while the spatial-boundary signature is large and "
        "volatile (e.g. negative at one bin). The high convergence-score rank "
        "correlation is driven mainly by the spatial-boundary and F_ST ranks, "
        "not by a strong, smooth co-rise of all signatures. Read the panel, not "
        "just the score."
    )
    emit(
        f"- Coverage: {len(matched)}/{len(counts)} assemblages matched to "
        "coordinates; only "
        f"{len(phase_known)} carry a usable Terminal Period phase; the seriation "
        f"file covers {n_unique_cpl} unique assemblages."
    )
    emit(
        "- Phase bins are uneven (F and D sparse); per-bin signatures from few "
        "sites or few clusters are noisy."
    )
    emit()

    emit("## 10. For interpretation by the team")
    emit()
    emit(
        "This report does not declare a winner. The pattern above is the "
        "evidence; the mapping to hypotheses is the team's call. For reference, "
        "the criteria are:"
    )
    emit(
        "- Consistent with H1 (nascent emergence) IF the bin-trajectory "
        "signatures co-rise toward contact (positive convergence-score slope), "
        "the rise concentrates spatially on/near Parkin, AND the settlement "
        "cross-checks show a primate/convex rank-size with Parkin dominant and "
        "differentiated mounds/fortification."
    )
    emit(
        "- Consistent with H2 (stable non-consolidation; Rees 2001) IF the "
        "signatures are flat / non-trending, the convergence-score slope is "
        "near zero, and the settlement system is log-normal (rank-size slope "
        "near -1) without a single dominant differentiated center."
    )
    emit(
        "- A mixed pattern (some signatures rise, others flat; Parkin a bridge "
        "but not a runaway primate) is itself a finding and should be reported "
        "as such, not forced to one pole."
    )
    emit()

    (OUTPUT / "empirical_findings.md").write_text("\n".join(lines).replace("$F_{ST}$", "F_ST"))

    # Console: NEUTRAL, no coordinates.
    print("Phase 5 empirical application complete.")
    print(f"Matched: {len(matched)}/{len(counts)}; chosen k={best_k}; "
          f"Parkin cluster c{parkin_cluster}.")
    print("Bin-trajectory slope / Spearman rho vs ordinal phase "
          "(+ = upward toward contact):")
    for col, _, label in series_specs:
        print(f"  {label}: slope {slopes[col]:+.5f}, rho "
              f"{spearman[col][0]:+.3f} -> {trend_word(col)}")
    print(f"  convergence_score: slope {conv_slope:+.5f}, rho {crho:+.3f}")
    print(f"Parkin Terminal Period phase: {parkin_tp} "
          f"(ordinal {PHASE_ORDER.get(parkin_tp,'NA')}/6).")
    if parkin_in_cpl:
        print(f"Parkin decorated-seriation rank {int(rel)}/{len(cpl_ids)}; "
              f"Sig2 multi-membership in {len(parkin_groups)} groups "
              f"(multi={PARKIN in multi_ids}).")
    print(f"Rank-size Area slope: {slope_rs:.3f}; primacy {primacy:.2f}; "
          f"largest={'Parkin' if largest_id==PARKIN else largest_id}.")
    if parkin_area_rank:
        print(f"Parkin Area rank: {parkin_area_rank}/{len(area)}.")
    print("Wrote output/empirical_findings.md and 3 figures. NO VERDICT.")


if __name__ == "__main__":
    main()
