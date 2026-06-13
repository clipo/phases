"""Phase 5 v2: two-level empirical application of the validated convergence
criterion to the real lower/central Mississippi Valley (LMV) data.

TWO LEVELS, two datasets:

  TRANSMISSION LEVEL (curated decorated set, ~55 assemblages)
    A correspondence-analysis (CA) seriation axis is built from the curated
    decorated count matrix (assemblages x ~10 decorated types). The CA first
    non-trivial axis is the seriation ordinate. The 14C calendar dates anchor
    and orient that axis (so larger ordinate = later in time). The four
    signatures (neutral departure, cultural F_ST, spatial boundary excess, and
    seriation/Signature-2 structure) are computed along the CA axis.

  SETTLEMENT LEVEL (broad ~258-site PFG/LMV set)
    Mound presence (binary), ditch presence (binary, defensive), St-Francis
    fortification count, and a mound-area rank-size analysis. Whether
    mound/ditch presence concentrates spatially near Parkin is checked.

  Parkin links the two levels (curated name 'Parkin'; broad-set id '11-N-1').

INTELLECTUAL HONESTY: this script produces the actual numbers and a NEUTRAL
pattern summary. It does NOT declare a winner between H1 (nascent emergence:
co-rise/convergence toward contact, concentrating on Parkin) and H2 (stable
non-consolidation, Rees 2001: flat / non-trending). If signals are flat, it
says so. Interpretation is the team's.

DATA POLICY: data/ is gitignored and location-sensitive. This script NEVER
prints or writes raw coordinates. Figures use a centered/relative frame with no
axis scale. Outputs (output/, figures/) are gitignored; only this script is
committed.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from mls_emergence.dataio.pfg import load_pfg_counts
from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv, normalize_grid
from mls_emergence.signatures.neutral import theta_f, theta_e
from mls_emergence.signatures.variance import cultural_fst
from mls_emergence.signatures.assortativity import (
    boundary_excess,
    _kmeans_labels,
    geo_distance,
)
from mls_emergence.signatures.seriation import seriation_groups
from mls_emergence.signatures.convergence import convergence_score, time_derivative

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "output"
FIGURES = ROOT / "figures"
OUTPUT.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

PARKIN_CUR = "Parkin"  # name in the curated set
PARKIN_BROAD = "11-N-1"  # grid id in the broad set

DECORATED_TYPES = [
    "Parkin_Punctated",
    "Barton/Kent/MPI",
    "Painted",
    "Fortune_Noded",
    "Ranch_Incised",
    "Walls_Engraved",
    "Wallace_Incised",
    "Rhodes_Incised",
    "Vernon_Paul_Applique",
    "Hull_Engraved",
]

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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def correspondence_axis(M: np.ndarray):
    """Correspondence analysis; return the first non-trivial row ordinate.

    M is an assemblages x types non-negative count matrix. Implements CA from
    numpy: P = M/M.sum(); r, c row/column masses; standardized residual matrix
    S = Dr^-1/2 (P - r c^T) Dc^-1/2; SVD; row scores = Dr^-1/2 U Sigma. The
    leading non-trivial dimension is column 0 of the SVD (the trivial unit
    singular vector is removed by centering, so the first SVD component is the
    seriation ordinate). Returns (ordinate, inertia_fraction).
    """
    M = np.asarray(M, float)
    total = M.sum()
    P = M / total
    r = P.sum(axis=1)  # row masses
    c = P.sum(axis=0)  # column masses
    # Guard zero masses: drop zero-mass columns; rows are guarded by caller.
    keep_c = c > 0
    P = P[:, keep_c]
    c = c[keep_c]
    Dr_inv = np.diag(1.0 / np.sqrt(r))
    Dc_inv = np.diag(1.0 / np.sqrt(c))
    S = Dr_inv @ (P - np.outer(r, c)) @ Dc_inv
    U, sig, Vt = np.linalg.svd(S, full_matrices=False)
    row_scores = Dr_inv @ U @ np.diag(sig)
    ordinate = row_scores[:, 0]
    inertia = sig ** 2
    frac = float(inertia[0] / inertia.sum()) if inertia.sum() > 0 else np.nan
    return ordinate, frac


def parse_cal_midpoint(s) -> float:
    """1-sigma midpoint from a calibrated string like '1280 (1310, 1350) 1429'.

    The first and last integers are the 1-sigma lower and upper bounds; the
    midpoint is their mean. Returns NaN if fewer than two integers are present.
    """
    if pd.isna(s):
        return np.nan
    nums = [int(x) for x in re.findall(r"\d+", str(s))]
    if len(nums) < 2:
        return np.nan
    return (nums[0] + nums[-1]) / 2.0


def norm_name(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


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


def main() -> None:
    lines: list[str] = []

    def emit(s: str = "") -> None:
        lines.append(s)

    emit("# Two-level empirical application of the convergence criterion (Phase 5 v2)")
    emit()
    emit(
        "Neutral pattern report on two datasets. TRANSMISSION level: a "
        "correspondence-analysis (CA) seriation axis on the curated decorated "
        "assemblages (~55), anchored and oriented by 14C calendar dates, with "
        "the four signatures computed along it. SETTLEMENT level: mound and "
        "ditch presence and a mound-area rank-size on the broad PFG/LMV set "
        "(~258). Parkin links the two. No verdict is declared between H1 "
        "(co-rise/convergence toward contact, concentrating on Parkin) and H2 "
        "(flat / non-trending; Rees 2001). See caveats and the closing section."
    )
    emit()

    # =======================================================================
    # TRANSMISSION LEVEL
    # =======================================================================
    emit("## Transmission level (curated decorated set)")
    emit()

    # -- Load curated counts ------------------------------------------------
    cur = pd.read_excel(
        DATA / "raw" / "mainfort-pfg-cpl.xlsx", sheet_name="pfg-cpl-mainfort"
    ).dropna(subset=["Assemblages"])
    cur["Assemblages"] = cur["Assemblages"].astype(str).str.strip()
    cur = cur.drop_duplicates(subset=["Assemblages"], keep="first").set_index(
        "Assemblages"
    )
    type_cols = [c for c in DECORATED_TYPES if c in cur.columns]
    counts = cur[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Drop assemblages with zero decorated total (zero row mass breaks CA).
    row_tot = counts.sum(axis=1)
    dropped_zero = list(counts.index[row_tot <= 0])
    counts = counts[row_tot > 0]

    emit("### 1. Data and join")
    emit()
    emit(
        f"- Curated decorated assemblages: {len(cur)} rows, {len(type_cols)} "
        f"decorated types ({', '.join(type_cols)})."
    )
    if dropped_zero:
        emit(
            f"- Dropped {len(dropped_zero)} assemblage(s) with zero decorated "
            f"total (undefined CA row mass): {', '.join(dropped_zero)}."
        )
    emit(f"- Assemblages used for CA: {len(counts)}.")

    # -- Coordinates --------------------------------------------------------
    xy = pd.read_csv(
        DATA / "raw" / "mainfort-pfg-cplXY.txt", sep="\t"
    )
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index(
        "Assemblages"
    )
    coords_ll = xy.reindex(counts.index)[["Latitude", "Longitude"]].apply(
        pd.to_numeric, errors="coerce"
    )
    n_with_coords = int(coords_ll.dropna().shape[0])
    emit(
        f"- Curated assemblages joined to coordinates "
        f"(mainfort-pfg-cplXY.txt): {n_with_coords} / {len(counts)}."
    )

    # -- 14C dates ----------------------------------------------------------
    rc = pd.read_excel(
        DATA / "raw" / "14CDatesFromMainfort2001.xls", sheet_name="Sheet1"
    )
    rc = rc[rc["Provenience"].notna() & (rc["Provenience"] != "Provenience")].copy()
    cal_col = "Calibrated Date A.D. (1 Sigma)"
    rc["cal_mid"] = rc[cal_col].map(parse_cal_midpoint)
    n_dates_parsed = int(rc["cal_mid"].notna().sum())
    prov_date = rc.groupby("Provenience")["cal_mid"].agg(["mean", "count"])
    emit(
        f"- 14C samples: {len(rc)} with a provenience; {n_dates_parsed} parsed "
        f"to a 1-sigma calendar midpoint; aggregated to {len(prov_date)} "
        "proveniences (mean per provenience)."
    )

    # Match proveniences to curated assemblages by normalized name (exact, then
    # substring). Report the mapping.
    cur_norm = {norm_name(a): a for a in counts.index}
    prov_to_assem: dict[str, str] = {}
    for prov in prov_date.index:
        pn = norm_name(prov)
        hit = cur_norm.get(pn)
        if hit is None:
            for k, v in cur_norm.items():
                if pn and (pn in k or k in pn):
                    hit = v
                    break
        if hit is not None:
            prov_to_assem[prov] = hit
    # Assemblage -> mean 14C date
    assem_date: dict[str, float] = {}
    for prov, assem in prov_to_assem.items():
        assem_date[assem] = float(prov_date.loc[prov, "mean"])
    n_matched_14c = len(assem_date)
    emit(
        f"- 14C proveniences matched to a curated assemblage: {n_matched_14c} "
        f"(of {len(prov_date)})."
    )
    emit("  Provenience -> assemblage (mean 1-sigma calendar AD, n samples):")
    for prov, assem in sorted(prov_to_assem.items()):
        emit(
            f"  - {prov} -> {assem}: AD {assem_date[assem]:.0f} "
            f"(n={int(prov_date.loc[prov, 'count'])})"
        )
    emit()

    # -- CA seriation axis --------------------------------------------------
    M = counts.to_numpy(float)
    ordinate, inertia_frac = correspondence_axis(M)
    ca = pd.Series(ordinate, index=counts.index, name="ca")

    # Orient with 14C: if CA is negatively correlated with calendar date, flip.
    dated_assems = [a for a in ca.index if a in assem_date]
    ca_d = ca.reindex(dated_assems).to_numpy(float)
    yr_d = np.array([assem_date[a] for a in dated_assems], float)
    if len(dated_assems) >= 3:
        rho_ca_date, p_ca_date = spearmanr(ca_d, yr_d)
    else:
        rho_ca_date, p_ca_date = np.nan, np.nan
    flipped = False
    if np.isfinite(rho_ca_date) and rho_ca_date < 0:
        ca = -ca
        ordinate = -ordinate
        flipped = True
        ca_d = -ca_d
        rho_ca_date = -rho_ca_date  # sign flips with the axis

    # Linear CA -> calendar map (coarse), reported only.
    if len(dated_assems) >= 2:
        b1, b0 = np.polyfit(ca_d, yr_d, 1)
    else:
        b1, b0 = np.nan, np.nan

    emit("### 2. CA seriation axis and 14C anchor")
    emit()
    emit(
        f"- CA first non-trivial axis carries {inertia_frac:.3f} of total "
        "inertia (the seriation ordinate; larger = later after orientation)."
    )
    emit(
        f"- Orientation: Spearman(CA ordinate, 14C mean date) on the "
        f"{len(dated_assems)} dated assemblages = {rho_ca_date:+.3f} "
        f"(p = {p_ca_date:.3f}); axis "
        + ("FLIPPED so increasing ordinate = later." if flipped
           else "kept (already increasing with time).")
    )
    if np.isfinite(b1):
        emit(
            f"- Coarse linear CA->calendar map (for reporting only): "
            f"AD = {b0:.0f} + {b1:.1f} * CA_ordinate "
            f"(n={len(dated_assems)} anchors)."
        )
    emit()

    # -- Spatial clustering of curated coords -------------------------------
    coords_df = coords_ll.dropna()
    have_coords_ids = list(coords_df.index)
    cc = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(axis=0)
    sil = {}
    for k in range(2, 7):
        lbl = _kmeans_labels(cc_c, k, seed=7)
        sil[k] = silhouette_mean(cc_c, lbl)
    best_k = max(sil, key=sil.get)
    k_use = best_k
    cl_labels = _kmeans_labels(cc_c, k_use, seed=7)
    cluster_of = dict(zip(have_coords_ids, cl_labels))
    emit("### 3. Spatial clustering of curated assemblages (k-means)")
    emit()
    emit(
        "- Silhouette by k: "
        + ", ".join(f"k={k}:{sil[k]:.3f}" for k in range(2, 7))
        + f"; chosen k = {k_use}."
    )
    sizes = pd.Series(cl_labels).value_counts().sort_index()
    emit(
        "- Cluster sizes: "
        + ", ".join(f"c{int(c)}:{int(n)}" for c, n in sizes.items())
        + "."
    )
    parkin_cluster = cluster_of.get(PARKIN_CUR)
    emit(
        "- Parkin cluster: "
        + (f"c{int(parkin_cluster)}." if parkin_cluster is not None
           else "Parkin has no coordinate.")
    )
    emit()

    # -- Four signatures along the CA axis (binned) -------------------------
    n_bins = 6
    ca_have = ca.reindex(have_coords_ids)
    # quantile bins on the CA ordinate (ordered windows along the axis)
    bins = pd.qcut(ca_have, q=n_bins, labels=False, duplicates="drop")
    bin_ids = sorted(pd.Series(bins).dropna().unique())
    counts_have = counts.reindex(have_coords_ids)
    panel_rows = {}
    bin_info = []
    sig2_by_bin = {}
    for b in bin_ids:
        ids = [i for i in have_coords_ids if not pd.isna(bins[i]) and bins[i] == b]
        sub_counts = counts_have.loc[ids].to_numpy(float)
        sub_coords = (cc_c[[have_coords_ids.index(i) for i in ids]])
        sub_clusters = np.array([cluster_of[i] for i in ids])
        ca_lo = float(ca_have.loc[ids].min())
        ca_hi = float(ca_have.loc[ids].max())

        # (a) neutral departure: mean within-cluster |1 - tf/te| (pooled)
        nd_vals = []
        for c in np.unique(sub_clusters):
            pooled = sub_counts[sub_clusters == c].sum(axis=0)
            if pooled.sum() < 2 or (pooled > 0).sum() < 2:
                continue
            tf, te = theta_f(pooled), theta_e(pooled)
            if np.isfinite(tf) and te > 0:
                nd_vals.append(abs(1.0 - tf / te))
        neutral_departure = float(np.mean(nd_vals)) if nd_vals else np.nan

        # (b) cultural F_ST across clusters represented in this bin
        rep = np.unique(sub_clusters)
        if len(rep) >= 2:
            gc = np.array([sub_counts[sub_clusters == c].sum(axis=0) for c in rep])
            fst = cultural_fst(gc)
        else:
            fst = np.nan

        # (c) spatial boundary excess
        spatial_boundary = (
            boundary_excess(sub_counts, sub_coords, seed=7) if len(ids) >= 4
            else np.nan
        )

        # (d) Signature 2 (IDSS) on the bin's assemblages, if tractable
        if 3 <= len(ids) <= 14:
            sg = seriation_groups(sub_counts, cont=0.30)
            sig2_by_bin[b] = sg["n_groups"]
        else:
            sig2_by_bin[b] = np.nan

        panel_rows[b] = {
            "neutral_departure": neutral_departure,
            "fst": fst,
            "spatial_boundary": spatial_boundary,
        }
        bin_info.append((b, len(ids), len(rep), ca_lo, ca_hi))

    panel = pd.DataFrame(panel_rows).T.sort_index()
    panel.index.name = "ca_bin"

    emit("### 4. Four signatures along the CA axis")
    emit()
    emit(
        f"Curated assemblages with coordinates binned into {len(bin_ids)} "
        "ordered windows along the CA ordinate (bin 0 = earliest)."
    )
    emit()
    emit("Bins (ca_bin, n assemblages, n clusters, CA range):")
    for b, ns, nc, lo, hi in bin_info:
        emit(f"- bin {int(b)}: n={ns}, clusters={nc}, CA [{lo:+.3f}, {hi:+.3f}]")
    emit()
    emit("| ca_bin | neutral_departure | fst | spatial_boundary | sig2_n_groups |")
    emit("|---|---|---|---|---|")
    for b, r in panel.iterrows():
        s2 = sig2_by_bin.get(b, np.nan)
        s2s = f"{int(s2)}" if (isinstance(s2, (int, float)) and np.isfinite(s2)) else "NA"
        emit(
            f"| {int(b)} | {r['neutral_departure']:.4f} | {r['fst']:.4f} | "
            f"{r['spatial_boundary']:.3f} | {s2s} |"
        )
    emit()
    emit(
        "Signature 2 (IDSS n_groups) is reported per bin where the bin is "
        "tractable (3-14 assemblages); it is NOT entered into the convergence "
        "score (it is a count of co-seriable groups, not a continuous "
        "magnitude, and bins above 14 are intractable)."
    )
    emit()

    # Convergence score + per-signature slope / Spearman vs CA-bin
    panel_full = panel.dropna(how="any")
    conv = convergence_score(panel_full) if len(panel_full) >= 2 else pd.Series(dtype=float)
    conv_slope = time_derivative(conv) if len(conv) >= 2 else np.nan
    if len(conv) >= 3:
        crho, cpval = spearmanr(conv.index.to_numpy(float), conv.values)
    else:
        crho, cpval = np.nan, np.nan

    slopes, spearman_sig, monotone = {}, {}, {}
    emit("Per-signature trend along the CA axis. OLS slope per bin (positive = "
         "rising toward later/contact), Spearman rank correlation with the "
         "ordered bin index, and monotonicity:")
    emit()
    for col in ["neutral_departure", "fst", "spatial_boundary"]:
        s = panel[col].dropna()
        slope = time_derivative(s) if len(s) >= 2 else np.nan
        slopes[col] = slope
        if len(s) >= 3:
            rho, pval = spearmanr(s.index.to_numpy(float), s.values)
        else:
            rho, pval = np.nan, np.nan
        spearman_sig[col] = (rho, pval)
        d = np.diff(s.values)
        mono = bool((d >= 0).all() or (d <= 0).all()) if len(s) >= 2 else False
        monotone[col] = mono
        emit(
            f"- {col}: slope = {slope:+.5f} (over {len(s)} bins); "
            f"Spearman rho = {rho:+.3f} (p = {pval:.3f}); monotone = {mono}"
        )
    emit(
        f"- convergence_score: slope = {conv_slope:+.5f}; "
        f"Spearman rho = {crho:+.3f} (p = {cpval:.3f}) "
        f"(over {len(panel_full)} complete bins)."
    )
    spearman_sig["convergence"] = (crho, cpval)
    emit()

    # -- Parkin link --------------------------------------------------------
    emit("### 5. Parkin link (transmission level)")
    emit()
    parkin_ca = float(ca.get(PARKIN_CUR, np.nan))
    parkin_ca_rank = int(ca.rank().get(PARKIN_CUR, np.nan)) if PARKIN_CUR in ca.index else None
    n_ca = len(ca)
    if parkin_ca_rank is not None:
        frac_rank = parkin_ca_rank / n_ca
        pos_label = "late" if frac_rank >= 0.66 else ("early" if frac_rank <= 0.33 else "mid")
        emit(
            f"- Parkin CA ordinate = {parkin_ca:+.3f}, rank {parkin_ca_rank}/"
            f"{n_ca} along the (oriented) axis => **{pos_label}**."
        )
        if np.isfinite(b1):
            emit(
                f"- Parkin mapped to calendar (coarse): AD "
                f"{b0 + b1 * parkin_ca:.0f}."
            )
    parkin_14c = assem_date.get(PARKIN_CUR)
    emit(
        "- Parkin 14C mean date: "
        + (f"AD {parkin_14c:.0f}." if parkin_14c is not None
           else "no direct 14C match in the curated set.")
    )

    # Parkin Signature-2 multi-membership: run IDSS on Parkin's spatial cluster
    # (a Parkin neighborhood), tractable if <= 14.
    if parkin_cluster is not None:
        neigh_ids = [i for i in have_coords_ids if cluster_of[i] == parkin_cluster]
        emit(
            f"- Parkin neighborhood (spatial cluster c{int(parkin_cluster)}): "
            f"{len(neigh_ids)} assemblages."
        )
        if 3 <= len(neigh_ids) <= 14:
            nm = counts_have.loc[neigh_ids].to_numpy(float)
            sg = seriation_groups(nm, cont=0.30)
            pidx = neigh_ids.index(PARKIN_CUR)
            pk_groups = sg["membership"][pidx]
            is_bridge = pidx in sg["multi"]
            emit(
                f"- Signature 2 in Parkin's neighborhood: {sg['n_groups']} "
                f"co-seriable groups; Parkin belongs to {len(pk_groups)} "
                f"(multi-membership / bridge = {is_bridge})."
            )
        else:
            emit(
                f"- Parkin neighborhood not IDSS-tractable "
                f"(n={len(neigh_ids)} outside 3-14); Signature-2 bridge status "
                "reported at the bin level (section 4) instead."
            )
    emit()

    # =======================================================================
    # SETTLEMENT LEVEL
    # =======================================================================
    emit("## Settlement level (broad PFG/LMV set)")
    emit()

    broad_counts = load_pfg_counts(DATA / "raw" / "PFGData.xlsx")
    n_broad_before = len(broad_counts)
    if not broad_counts.index.is_unique:
        broad_counts = broad_counts.groupby(level=0).sum()
    lmv = load_lmv(DATA / "LMVData.xlsx")
    joined, unmatched = join_pfg_to_lmv(broad_counts, lmv)
    bmatched = joined.dropna(subset=["Easting", "Northing"]).copy()

    # Clean binary features from LMVData-22March2006.xls, joined by Number.
    lmv2 = pd.read_excel(DATA / "LMVData-22March2006.xls", sheet_name="Sheet1")
    lmv2 = lmv2.dropna(subset=["Number"]).copy()
    lmv2["_k"] = lmv2["Number"].astype(str).map(normalize_grid)
    lmv2 = lmv2.drop_duplicates(subset=["_k"], keep="first").set_index("_k")
    norm_ids = pd.Index([normalize_grid(str(i)) for i in bmatched.index])
    ext = lmv2.reindex(norm_ids)
    n_ext = int(ext["Number"].notna().sum())

    def to_bin(series) -> pd.Series:
        v = pd.to_numeric(
            series.astype(str).str.replace("?", "", regex=False), errors="coerce"
        )
        return (v.fillna(0) > 0)

    mound = to_bin(ext["Mound"])
    ditch = to_bin(ext["Ditch"])
    stfr = to_bin(ext["St Francis"])
    platform = to_bin(ext["Platform"])
    num_mounds = pd.to_numeric(ext["Num_Mounds"], errors="coerce")
    mound_area = pd.to_numeric(ext["Max Mound Area (sq ft)"], errors="coerce")
    mound_ht = pd.to_numeric(ext["Max Mound Height (ft)"], errors="coerce")
    mound.index = bmatched.index
    ditch.index = bmatched.index
    stfr.index = bmatched.index
    platform.index = bmatched.index
    num_mounds.index = bmatched.index
    mound_area.index = bmatched.index
    mound_ht.index = bmatched.index

    emit("### 6. Settlement features")
    emit()
    emit(
        f"- Broad PFG assemblages: {n_broad_before} rows -> {len(broad_counts)} "
        f"unique ids; matched to LMV coordinates: {len(bmatched)}; of those, "
        f"{n_ext} joined to LMVData-22March2006 binary features by Number."
    )
    emit(f"- Parkin ({PARKIN_BROAD}) present in broad matched set: "
         f"{PARKIN_BROAD in bmatched.index}.")
    emit(
        f"- Mound present: {int(mound.sum())} / {len(bmatched)} "
        f"({100*mound.mean():.1f}%)."
    )
    emit(
        f"- Ditch present (defensive): {int(ditch.sum())} / {len(bmatched)} "
        f"({100*ditch.mean():.1f}%)."
    )
    # Ditch context across the FULL LMV-22 table, not just PFG-matched.
    full_ditch = to_bin(lmv2["Ditch"])
    emit(
        f"- DATA NOTE: across the full LMVData-22March table "
        f"({len(lmv2)} sites) ditches number {int(full_ditch.sum())}; the "
        f"PFG-matched subset captures {int(ditch.sum())} of them. The PFG "
        "sites are a ceramic-collection subset, not a random settlement "
        "sample, and are mound-biased (see mound %)."
    )
    emit(
        f"- St-Francis (fortified) present: {int(stfr.sum())} / {len(bmatched)}."
    )
    emit(
        f"- Platform present: {int(platform.sum())} / {len(bmatched)}."
    )
    # Parkin features
    emit(
        f"- Parkin features (LMVData-22March): "
        f"Mound={bool(mound.get(PARKIN_BROAD, False))}, "
        f"Ditch={bool(ditch.get(PARKIN_BROAD, False))}, "
        f"St-Francis={bool(stfr.get(PARKIN_BROAD, False))}, "
        f"Num_Mounds={num_mounds.get(PARKIN_BROAD, np.nan):.0f}, "
        f"Max_Mound_Height={mound_ht.get(PARKIN_BROAD, np.nan):.0f} ft, "
        f"Max_Mound_Area={mound_area.get(PARKIN_BROAD, np.nan):.0f} sq ft."
    )
    emit(
        "- DATA NOTE: the LMV coding records Parkin Ditch=0 and Max Mound "
        "Area=0, even though the Parkin site is the type-site of a "
        "ditched/palisaded fortified town. The LMV ditch/area fields do not "
        "capture Parkin's known defensive ditch; this is a data-coding "
        "limitation, flagged for the team."
    )
    emit()

    # -- Mound-area rank-size ----------------------------------------------
    area = mound_area.dropna()
    area = area[area > 0].sort_values(ascending=False)
    ranks = np.arange(1, len(area) + 1)
    if len(area) >= 2:
        logr = np.log(ranks)
        loga = np.log(area.values)
        A = np.vstack([logr, np.ones_like(logr)]).T
        slope_rs, intercept_rs = np.linalg.lstsq(A, loga, rcond=None)[0]
        primacy = float(area.iloc[0] / area.iloc[1])
    else:
        slope_rs, intercept_rs, primacy = np.nan, np.nan, np.nan
    largest_id = area.index[0] if len(area) else None
    parkin_area_rank = (
        int((area.index == PARKIN_BROAD).argmax() + 1)
        if PARKIN_BROAD in area.index else None
    )
    emit("### 7. Mound-area rank-size")
    emit()
    emit(f"- Sites with Max Mound Area > 0: {len(area)}.")
    if np.isfinite(slope_rs):
        emit(
            f"- log-log rank-size slope = {slope_rs:.3f} "
            "(Zipf/log-normal expectation ~= -1; shallower = convex below the "
            "top, steeper = primate)."
        )
        emit(
            f"- Primacy (largest/second) = {primacy:.2f}; largest site = "
            f"{largest_id}; largest is Parkin = {largest_id == PARKIN_BROAD}."
        )
    if parkin_area_rank:
        emit(
            f"- Parkin mound-area rank = {parkin_area_rank}/{len(area)}."
        )
    else:
        emit("- Parkin has no recorded Max Mound Area > 0 (excluded from this "
             "curve; see height below).")
    # Mound height context for Parkin (since area=0)
    mh = mound_ht.dropna()
    if PARKIN_BROAD in mh.index:
        ph = mh.loc[PARKIN_BROAD]
        pct = float((mh < ph).mean() * 100)
        rank_h = int((mh.sort_values(ascending=False).index == PARKIN_BROAD).argmax() + 1)
        emit(
            f"- Parkin Max Mound Height = {ph:.0f} ft, percentile {pct:.1f}, "
            f"rank {rank_h}/{len(mh)} (n with recorded height)."
        )
    emit()

    # -- Spatial concentration of mound/ditch near Parkin -------------------
    emit("### 8. Spatial concentration near Parkin")
    emit()
    bcoords = bmatched[["Northing", "Easting"]].to_numpy(float)
    bcoords_c = bcoords - bcoords.mean(axis=0)
    bsil = {}
    for k in range(4, 11):
        lbl = _kmeans_labels(bcoords_c, k, seed=7)
        bsil[k] = silhouette_mean(bcoords_c, lbl)
    bk = max(bsil, key=bsil.get)
    blabels = _kmeans_labels(bcoords_c, bk, seed=7)
    bmatched["cluster"] = blabels
    pk_cl = int(bmatched.loc[PARKIN_BROAD, "cluster"]) if PARKIN_BROAD in bmatched.index else None
    in_pk = bmatched["cluster"].to_numpy() == pk_cl
    emit(
        f"- Broad k-means: chosen k={bk} (max silhouette); Parkin cluster "
        f"c{pk_cl}, size {int(in_pk.sum())}."
    )
    if pk_cl is not None:
        m_in = mound[in_pk].mean()
        m_out = mound[~in_pk].mean()
        d_in = ditch[in_pk].mean()
        d_out = ditch[~in_pk].mean()
        emit(
            f"- Mound presence: {100*m_in:.1f}% inside Parkin's cluster vs "
            f"{100*m_out:.1f}% outside."
        )
        emit(
            f"- Ditch presence: {100*d_in:.1f}% inside Parkin's cluster vs "
            f"{100*d_out:.1f}% outside "
            f"(n ditched inside = {int(ditch[in_pk].sum())}, "
            f"outside = {int(ditch[~in_pk].sum())})."
        )
    emit()

    # =======================================================================
    # FIGURES
    # =======================================================================
    def zscore(s):
        s = s.astype(float)
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd and np.isfinite(sd) else s * 0

    # (a) CA-axis trajectory of the four signatures + convergence
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = panel.index.to_numpy(float)
    specs = [
        ("neutral_departure", OKABE["blue"], "Neutral departure"),
        ("fst", OKABE["orange"], "Cultural F_ST"),
        ("spatial_boundary", OKABE["green"], "Spatial boundary excess"),
    ]
    for col, color, label in specs:
        ax.plot(x, zscore(panel[col]), marker="o", color=color, label=label, linewidth=1.6)
    s2series = pd.Series(sig2_by_bin).reindex(panel.index).astype(float)
    if s2series.notna().sum() >= 2:
        ax.plot(panel.index.to_numpy(float), zscore(s2series), marker="^",
                color=OKABE["purple"], label="Sig2 n_groups", linewidth=1.2,
                linestyle=":")
    if len(conv) >= 2:
        ax.plot(conv.index.to_numpy(float), zscore(conv), marker="s",
                color=OKABE["black"], linewidth=2.0, linestyle="--",
                label="Convergence score")
    ax.set_xlabel("CA seriation axis (bin; 0 = earliest, increasing = later)")
    ax.set_ylabel("Signature (z-standardized)")
    ax.set_xticks(sorted(int(v) for v in x))
    ax.axhline(0, color="0.7", linewidth=0.6, zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(FIGURES / "06_ca_trajectory.png")
    plt.close(fig)

    # (b) CA ordinate vs 14C date scatter
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.scatter(ca_d, yr_d, s=40, color=OKABE["blue"], edgecolors="white", linewidths=0.5)
    for i, a in enumerate(dated_assems):
        ax.annotate(a.replace("_", " "), (ca_d[i], yr_d[i]), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    if np.isfinite(b1):
        xs = np.linspace(ca_d.min(), ca_d.max(), 50)
        ax.plot(xs, b0 + b1 * xs, "-", color=OKABE["vermillion"], linewidth=1.4,
                label=f"OLS: rho={rho_ca_date:+.2f}")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("CA seriation ordinate (oriented; larger = later)")
    ax.set_ylabel("Mean 14C calendar date (AD, 1-sigma midpoint)")
    fig.tight_layout()
    fig.savefig(FIGURES / "06_ca_vs_14c.png")
    plt.close(fig)

    # (c) mound-area rank-size, Parkin marked
    if len(area) >= 2:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        ax.loglog(ranks, area.values, "o", ms=4, color=OKABE["sky"], alpha=0.7)
        fit = np.exp(intercept_rs) * ranks ** slope_rs
        ax.loglog(ranks, fit, "-", color=OKABE["vermillion"], linewidth=1.5,
                  label=f"OLS slope = {slope_rs:.2f}")
        if parkin_area_rank:
            ax.loglog(parkin_area_rank, area.loc[PARKIN_BROAD], "*", ms=16,
                      color=OKABE["black"], label=f"Parkin (rank {parkin_area_rank})")
        ax.set_xlabel("Rank (log)")
        ax.set_ylabel("Max mound area, sq ft (log)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURES / "06_mound_rank_size.png")
        plt.close(fig)

    # (d) mound/ditch map (relative, jittered, no axis scale)
    fig, ax = plt.subplots(figsize=(5.0, 5.5))
    rng = np.random.default_rng(0)
    jit = rng.normal(0, bcoords_c.std() * 0.01, bcoords_c.shape)
    pj = bcoords_c + jit
    no_m = ~mound.to_numpy()
    ax.scatter(pj[no_m, 1], pj[no_m, 0], s=8, color="0.7", alpha=0.6,
               edgecolors="none", label="no mound")
    has_m = mound.to_numpy() & ~ditch.to_numpy()
    ax.scatter(pj[has_m, 1], pj[has_m, 0], s=16, color=OKABE["blue"], alpha=0.8,
               edgecolors="none", label="mound")
    has_d = ditch.to_numpy()
    ax.scatter(pj[has_d, 1], pj[has_d, 0], s=40, color=OKABE["vermillion"],
               marker="D", edgecolors="black", linewidths=0.4, label="ditch")
    if PARKIN_BROAD in bmatched.index:
        pidx = bmatched.index.get_loc(PARKIN_BROAD)
        ax.scatter(pj[pidx, 1], pj[pidx, 0], s=140, marker="*",
                   color=OKABE["black"], edgecolors="white", linewidths=0.6,
                   label="Parkin", zorder=6)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("Relative easting (no scale)")
    ax.set_ylabel("Relative northing (no scale)")
    ax.set_aspect("equal")
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(FIGURES / "06_mound_ditch_map.png")
    plt.close(fig)

    emit("## 9. Figures")
    emit()
    emit("- figures/06_ca_trajectory.png: the four signatures (z-standardized) "
         "plus convergence score across the CA seriation axis.")
    emit("- figures/06_ca_vs_14c.png: CA ordinate vs mean 14C calendar date for "
         "the dated curated assemblages, with the orientation fit.")
    emit("- figures/06_mound_rank_size.png: log-log mound-area rank-size with "
         "OLS slope; Parkin marked (Parkin has area=0 so it is absent unless "
         "recorded).")
    emit("- figures/06_mound_ditch_map.png: relative jittered settlement map "
         "(no axis scale) of mound/ditch presence with Parkin starred.")
    emit()

    # =======================================================================
    # NEUTRAL PATTERN SUMMARY
    # =======================================================================
    emit("## 10. Pattern summary (NEUTRAL)")
    emit()

    def trend_word(col):
        slope = slopes.get(col, np.nan)
        rho = spearman_sig.get(col, (np.nan, np.nan))[0]
        mono = monotone.get(col, False)
        if not np.isfinite(slope):
            return "undetermined"
        if not np.isfinite(rho) or abs(rho) < 0.3:
            return "flat / non-trending"
        direction = "upward" if slope > 0 else "downward"
        if mono and abs(rho) >= 0.7:
            return f"monotone {direction} toward later"
        return f"{direction} but non-monotone"

    emit(
        "Trend labels use both the OLS slope sign and the Spearman rank "
        "correlation; with only a handful of bins these are descriptive, not "
        "inferential."
    )
    emit()
    for col, _, label in specs:
        rho = spearman_sig[col][0]
        emit(
            f"- {label}: slope {slopes[col]:+.5f}, Spearman rho {rho:+.3f} -> "
            f"**{trend_word(col)}**."
        )
    crho2 = spearman_sig.get("convergence", (np.nan, np.nan))[0]
    if not np.isfinite(crho2) or abs(crho2) < 0.3:
        cw = "flat / non-trending"
    else:
        cw = ("upward" if conv_slope > 0 else "downward") + (
            " (monotone)" if abs(crho2) >= 0.7 else " (non-monotone)"
        )
    emit(
        f"- Convergence score: slope {conv_slope:+.5f}, Spearman rho "
        f"{crho2:+.3f} -> **{cw}**."
    )
    rising = sum(1 for c, _, _ in specs if "upward" in trend_word(c))
    flat = sum(1 for c, _, _ in specs if trend_word(c).startswith("flat"))
    emit()
    emit(
        f"- Of the 3 continuous signatures: {rising} trend upward toward "
        f"later, {flat} flat/non-trending, {3 - rising - flat} downward. "
        f"Convergence-score slope {conv_slope:+.5f} (rho {crho2:+.3f})."
    )
    emit(
        f"- CA<->14C: Spearman {rho_ca_date:+.3f} on {len(dated_assems)} "
        f"anchors (p {p_ca_date:.3f}); CA axis "
        + ("oriented by flipping." if flipped else "already time-ordered.")
    )
    emit(
        "- Parkin (transmission): CA position "
        + (f"rank {parkin_ca_rank}/{n_ca} ({pos_label}); "
           if parkin_ca_rank is not None else "n/a; ")
        + (f"14C AD {parkin_14c:.0f}; " if parkin_14c is not None else "no 14C; ")
        + "Signature-2 bridge status as reported in section 5."
    )
    emit(
        f"- Parkin (settlement): mound present "
        f"{bool(mound.get(PARKIN_BROAD, False))}, ditch coded "
        f"{bool(ditch.get(PARKIN_BROAD, False))} (see data note), height "
        f"{mound_ht.get(PARKIN_BROAD, np.nan):.0f} ft."
    )
    emit(
        f"- Settlement system: mound % = {100*mound.mean():.1f}, ditch % = "
        f"{100*ditch.mean():.1f}, mound-area rank-size slope = {slope_rs:.2f}, "
        f"primacy = {primacy:.2f}, largest = "
        f"{'Parkin' if largest_id == PARKIN_BROAD else largest_id}."
    )
    emit()
    emit(
        "Plain statement: report the panel, not a single index. Whether the "
        "transmission signatures co-rise toward contact, whether Parkin is an "
        "early/mid/late bridge, and whether the settlement system is primate "
        "or log-normal are the observed quantities above. No pole is asserted."
    )
    emit()

    # =======================================================================
    # CAVEATS
    # =======================================================================
    emit("## 11. Caveats (explicit)")
    emit()
    emit(
        "- The CA axis is RELATIVE. Even oriented and coarsely calibrated by "
        f"{len(dated_assems)} 14C anchors, it is a seriation ordinate, not "
        "calendar time; per-bin slopes are not rates."
    )
    emit(
        "- 14C anchoring is sparse: only "
        f"{n_matched_14c} proveniences match curated assemblages, and the "
        "name-matching used exact-then-substring matching (e.g. Kent -> "
        "Kent_Place); the CA<->date correlation rests on these few points."
    )
    emit(
        "- Date parsing approximates each sample by the mean of the first and "
        "last integers in the 1-sigma calibrated string (the 1-sigma "
        "midpoint); multi-intercept samples are summarized by this midpoint, "
        "and per-provenience dates are simple means over samples."
    )
    emit(
        "- The two levels use different datasets and ID conventions: ~55 "
        "curated decorated assemblages (name keys) vs ~258 broad PFG sites "
        "(grid-id keys). Parkin is the only guaranteed cross-level link."
    )
    emit(
        "- The broad PFG set is a ceramic-collection subset, mound-biased "
        f"({100*mound.mean():.0f}% have a mound) and capturing only "
        f"{int(ditch.sum())} of {int(full_ditch.sum())} ditched sites in the "
        "full LMV table; settlement proportions are NOT a random settlement "
        "sample."
    )
    emit(
        "- The LMV coding records Parkin Ditch=0 and Mound-Area=0 despite "
        "Parkin being a known ditched town with mounds; the LMV ditch/area "
        "fields under-record fortification. Treat the ditch counts as a floor."
    )
    emit(
        "- The four signatures share the type-frequency substrate and are "
        "partly correlated; some co-movement is expected even without a single "
        "causal process."
    )
    emit(
        "- Signature 2 (IDSS) is a count of co-seriable groups, computed only "
        "where a bin or neighborhood is tractable (3-14 assemblages); it is "
        "shown for context but excluded from the convergence score."
    )
    emit(
        f"- Bins are few ({len(bin_ids)}) and uneven; per-bin signatures from "
        "small bins or few clusters are noisy."
    )
    emit()

    # =======================================================================
    # FOR TEAM INTERPRETATION (no verdict)
    # =======================================================================
    emit("## 12. For team interpretation (no verdict)")
    emit()
    emit(
        "This report does not declare a winner. The numbers above are the "
        "evidence; mapping to hypotheses is the team's call. For reference:"
    )
    emit(
        "- Consistent with H1 (nascent emergence) IF the transmission "
        "signatures co-rise toward the later end of the CA axis (positive "
        "convergence-score slope), the rise concentrates on/near Parkin, and "
        "the settlement system shows a primate/convex mound-area rank-size "
        "with Parkin dominant and differentiated."
    )
    emit(
        "- Consistent with H2 (stable non-consolidation; Rees 2001) IF the "
        "signatures are flat / non-trending along the CA axis, the "
        "convergence-score slope is near zero, and the settlement system is "
        "log-normal (rank-size slope near -1) without a single dominant "
        "differentiated center."
    )
    emit(
        "- A mixed pattern (some signatures rise, others flat; Parkin a bridge "
        "but not a runaway primate) is itself a finding and should be reported "
        "as such, not forced to one pole."
    )
    emit()

    (OUTPUT / "empirical_findings_v2.md").write_text("\n".join(lines))

    # Console: NEUTRAL, no coordinates.
    print("Phase 5 v2 two-level empirical application complete.")
    print(f"[Transmission] curated CA assemblages={len(counts)}; "
          f"with coords={n_with_coords}; with 14C match={n_matched_14c}.")
    print(f"[Transmission] CA<->14C Spearman={rho_ca_date:+.3f} "
          f"(n={len(dated_assems)}, p={p_ca_date:.3f}); flipped={flipped}.")
    print("[Transmission] signature slope / rho along CA axis:")
    for col, _, label in specs:
        print(f"  {label}: slope {slopes[col]:+.5f}, rho "
              f"{spearman_sig[col][0]:+.3f} -> {trend_word(col)}")
    print(f"  convergence_score: slope {conv_slope:+.5f}, rho {crho2:+.3f}")
    if parkin_ca_rank is not None:
        print(f"[Parkin] CA rank {parkin_ca_rank}/{n_ca} ({pos_label}); "
              f"14C={'AD %.0f'%parkin_14c if parkin_14c else 'none'}.")
    print(f"[Settlement] mound%={100*mound.mean():.1f}, ditch%={100*ditch.mean():.1f}, "
          f"rank-size slope={slope_rs:.2f}, primacy={primacy:.2f}.")
    print(f"[Settlement] Parkin ditch coded={bool(ditch.get(PARKIN_BROAD,False))} "
          f"(LMV under-records; see findings), height="
          f"{mound_ht.get(PARKIN_BROAD,np.nan):.0f} ft.")
    print("Wrote output/empirical_findings_v2.md and figures/06_*.png. NO VERDICT.")


if __name__ == "__main__":
    main()
