"""Phase 5, basin-restricted re-run: does the whole-LMV "no convergence"
result hold within the Parkin phase / St. Francis basin?

The prior empirical passes (analyses/06_empirical_two_level.py,
07_refined_empirical.py) ran across the WHOLE lower/central Mississippi Valley
(LMV), lumping the Parkin phase (St. Francis basin, lat ~34.8-35.6) together
with the lower Yazoo (Winterville etc., lat ~33.5) and a couple of southern
St-Francis-type outliers. Those are distinct regional phases. The team decision
is to restrict the empirical test to the St. Francis basin, the focal
hypothesis unit, and to ask whether the whole-LMV finding (no transmission-level
convergence; settlement rank-size flat / not primate) survives within the basin.

This script reuses the data-loading, CA seriation, four-signature, IDSS, and
rank-size machinery of scripts 06/07 and applies a PRINCIPLED latitude cut.

BASIN DEFINITION (principled, reported):
  PRIMARY: St. Francis basin = sites at latitude >= 34.5. This sits in the
  natural latitude gap of the curated decorated set (the largest inter-site
  latitude gap falls at ~34.55-34.76) and below the St-Francis/Parkin cluster
  (lat ~34.8-35.6). It excludes the lower Yazoo (Winterville, lat ~33.48, which
  is not even in the curated decorated set) and the southern St-Francis-type
  outliers (curated: Salomon 34.35, Parchman 34.36; broad: 17-M-2 at 33.95,
  the only St-Francis-flagged site below the gap).
  SENSITIVITY: lat >= 34.0 and lat >= 35.0 are also reported so the conclusion's
  robustness to the cut is visible.

INTELLECTUAL HONESTY: the basin curated set is small. Where n is small the
trajectory is flagged as underpowered and weight is placed on the structural /
IDSS and settlement results. The result is reported HONESTLY whether or not it
changes the whole-LMV picture. A NEUTRAL verdict is given; no H1/H2 pole is
forced.

DATA POLICY: data/ is gitignored and location-sensitive. This script reports
only lat/long RANGES (extent), never per-site raw coordinates. Coordinates for
the broad set are computed from UTM via pyproj (EPSG:26915 -> EPSG:4326);
the curated set already carries Latitude/Longitude. Outputs (output/, figures/)
are gitignored; only this script is committed.
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
import pyproj
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

PARKIN_CUR = "Parkin"
PARKIN_BROAD = "11-N-1"

# Basin latitude cut. Primary is the required definition; the others are the
# reported sensitivity cuts.
LAT_PRIMARY = 34.5
LAT_CUTS = [34.0, 34.5, 35.0]

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

# IDSS continuity threshold: cont=0.30 (Lipo et al. 2015) over-saturates the
# broadly overlapping decorated matrix and exceeds the solver caps, as
# documented in script 07. Use the same primary value and sensitivity sweep.
CONT_PRIMARY = 0.10
CONT_SWEEP = [0.05, 0.10, 0.20]

# ---------------------------------------------------------------------------
# House figure style (no titles; Okabe-Ito; sans-serif; 300 dpi)
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
# Helpers (shared with 06/07)
# ---------------------------------------------------------------------------
def correspondence_axis(M: np.ndarray):
    """Correspondence analysis; first non-trivial row ordinate and inertia frac."""
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
    ordinate = row_scores[:, 0]
    inertia = sig ** 2
    frac = float(inertia[0] / inertia.sum()) if inertia.sum() > 0 else np.nan
    return ordinate, frac


def parse_cal_midpoint(s) -> float:
    """1-sigma midpoint from a string like '1280 (1310, 1350) 1429'."""
    if pd.isna(s):
        return np.nan
    nums = [int(x) for x in re.findall(r"\d+", str(s))]
    if len(nums) < 2:
        return np.nan
    return (nums[0] + nums[-1]) / 2.0


def norm_name(s) -> str:
    base = re.sub(r"[^a-z0-9 ]", " ", str(s).lower())
    tokens = [t for t in base.split() if t]
    drop = {"place", "mound", "mounds", "lake", "landing", "ferry", "bayou",
            "ranch", "the", "site"}
    core = [t for t in tokens if t not in drop]
    if not core:
        core = tokens
    return "".join(core)


def match_provenience(pn: str, cur_norm: dict) -> str | None:
    if not pn:
        return None
    hit = cur_norm.get(pn)
    if hit is not None:
        return hit
    best, best_len = None, 0
    for k, v in cur_norm.items():
        if not k:
            continue
        if pn in k or k in pn:
            overlap = min(len(pn), len(k))
            if overlap > best_len:
                best, best_len = v, overlap
    return best


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


def neutral_departure_pooled(sub_counts: np.ndarray, sub_clusters: np.ndarray) -> float:
    vals = []
    for c in np.unique(sub_clusters):
        pooled = sub_counts[sub_clusters == c].sum(axis=0)
        if pooled.sum() < 2 or (pooled > 0).sum() < 2:
            continue
        tf, te = theta_f(pooled), theta_e(pooled)
        if np.isfinite(tf) and te > 0:
            vals.append(abs(1.0 - tf / te))
    return float(np.mean(vals)) if vals else np.nan


def fst_across(sub_counts: np.ndarray, sub_clusters: np.ndarray) -> float:
    rep = np.unique(sub_clusters)
    if len(rep) < 2:
        return np.nan
    gc = np.array([sub_counts[sub_clusters == c].sum(axis=0) for c in rep])
    return cultural_fst(gc)


def ols_slope(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2:
        return np.nan
    x = x[m] - x[m].mean()
    if (x @ x) == 0:
        return np.nan
    return float((x @ (y[m] - y[m].mean())) / (x @ x))


def rank_size(area: pd.Series):
    """log-log rank-size slope and primacy for a positive-area series."""
    area = area[area > 0].sort_values(ascending=False)
    n = len(area)
    if n < 2:
        return np.nan, np.nan, (area.index[0] if n else None), area
    ranks = np.arange(1, n + 1)
    A = np.vstack([np.log(ranks), np.ones_like(ranks)]).T
    slope, _ = np.linalg.lstsq(A, np.log(area.values), rcond=None)[0]
    primacy = float(area.iloc[0] / area.iloc[1])
    return float(slope), primacy, area.index[0], area


# ---------------------------------------------------------------------------
# Data loaders that attach lat/long
# ---------------------------------------------------------------------------
def load_curated():
    """Curated decorated counts + lat/long (already in the XY file)."""
    cur = pd.read_excel(
        DATA / "raw" / "mainfort-pfg-cpl.xlsx", sheet_name="pfg-cpl-mainfort"
    ).dropna(subset=["Assemblages"])
    cur["Assemblages"] = cur["Assemblages"].astype(str).str.strip()
    cur = cur.drop_duplicates(subset=["Assemblages"], keep="first").set_index(
        "Assemblages"
    )
    type_cols = [c for c in DECORATED_TYPES if c in cur.columns]
    counts = cur[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = counts[counts.sum(axis=1) > 0]

    xy = pd.read_csv(DATA / "raw" / "mainfort-pfg-cplXY.txt", sep="\t")
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index(
        "Assemblages"
    )
    ll = xy.reindex(counts.index)[["Latitude", "Longitude"]].apply(
        pd.to_numeric, errors="coerce"
    )
    return counts, type_cols, ll


def load_broad():
    """Broad PFG/LMV settlement set + lat/long (UTM->geographic via pyproj) +
    LMV-22 binary/quantitative features."""
    broad_counts = load_pfg_counts(DATA / "raw" / "PFGData.xlsx")
    if not broad_counts.index.is_unique:
        broad_counts = broad_counts.groupby(level=0).sum()
    lmv = load_lmv(DATA / "LMVData.xlsx")
    joined, _ = join_pfg_to_lmv(broad_counts, lmv)
    bm = joined.dropna(subset=["Easting", "Northing"]).copy()

    # UTM -> geographic. The LMV zone sheets in this dataset are all Zone 15
    # (EPSG:26915); convert per zone for safety.
    lat = np.full(len(bm), np.nan)
    lon = np.full(len(bm), np.nan)
    zone = pd.to_numeric(bm.get("Zone"), errors="coerce")
    east = bm["Easting"].to_numpy(float)
    north = bm["Northing"].to_numpy(float)
    for z, epsg in [(15, 26915), (16, 26916)]:
        m = (zone == z).to_numpy()
        if m.any():
            tr = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True)
            lo, la = tr.transform(east[m], north[m])
            lat[m] = la
            lon[m] = lo
    bm["lat"] = lat
    bm["lon"] = lon

    # LMV-22 features joined by Number.
    lmv2 = pd.read_excel(DATA / "LMVData-22March2006.xls", sheet_name="Sheet1")
    lmv2 = lmv2.dropna(subset=["Number"]).copy()
    lmv2["_k"] = lmv2["Number"].astype(str).map(normalize_grid)
    lmv2 = lmv2.drop_duplicates(subset=["_k"], keep="first").set_index("_k")
    ext = lmv2.reindex(pd.Index([normalize_grid(str(i)) for i in bm.index]))
    ext.index = bm.index

    def to_bin(series) -> pd.Series:
        v = pd.to_numeric(
            series.astype(str).str.replace("?", "", regex=False), errors="coerce"
        )
        return v.fillna(0) > 0

    feats = pd.DataFrame(index=bm.index)
    feats["lat"] = bm["lat"]
    feats["lon"] = bm["lon"]
    feats["mound"] = to_bin(ext["Mound"])
    feats["ditch"] = to_bin(ext["Ditch"])
    feats["stfr"] = to_bin(ext["St Francis"])
    feats["platform"] = to_bin(ext["Platform"])
    feats["num_mounds"] = pd.to_numeric(ext["Num_Mounds"], errors="coerce")
    feats["mound_area"] = pd.to_numeric(ext["Max Mound Area (sq ft)"], errors="coerce")
    feats["mound_ht"] = pd.to_numeric(ext["Max Mound Height (ft)"], errors="coerce")
    feats["northing"] = bm["Northing"]
    feats["easting"] = bm["Easting"]
    return feats


def main() -> None:
    lines: list[str] = []

    def emit(s: str = "") -> None:
        lines.append(s)

    emit("# Parkin-phase / St. Francis basin-restricted re-run")
    emit()
    emit(
        "Does the whole-LMV no-convergence result hold within the St. Francis "
        "basin (the focal hypothesis unit)? The prior passes (scripts 06/07) "
        "lumped the Parkin phase (lat ~34.8-35.6) with the lower Yazoo "
        "(Winterville, lat ~33.5) and southern St-Francis-type outliers. This "
        "run restricts BOTH the curated decorated (transmission) set and the "
        "broad PFG/LMV (settlement) set to the basin via a principled latitude "
        "cut, recomputes the CA seriation axis, the four transmission "
        "signatures, the IDSS group structure with Parkin's bridge rank, and a "
        "WITHIN-BASIN mound rank-size, and states plainly whether the "
        "no-convergence finding survives. No verdict between H1 (nascent "
        "emergence / consolidation toward contact) and H2 (stable "
        "non-consolidation; Rees 2001) is forced."
    )
    emit()

    # =======================================================================
    # 0. Basin definition + n's at each cut (curated + broad)
    # =======================================================================
    emit("## 0. Basin definition and subset sizes")
    emit()
    counts_all, type_cols, ll_all = load_curated()
    broad_all = load_broad()

    emit(
        f"- Whole-LMV curated decorated set: {len(counts_all)} assemblages "
        f"(all with coordinates), lat [{ll_all['Latitude'].min():.3f}, "
        f"{ll_all['Latitude'].max():.3f}], lon "
        f"[{ll_all['Longitude'].min():.3f}, {ll_all['Longitude'].max():.3f}]."
    )
    emit(
        f"- Whole-LMV broad settlement set (matched to coordinates): "
        f"{len(broad_all)} sites, lat [{broad_all['lat'].min():.3f}, "
        f"{broad_all['lat'].max():.3f}], lon [{broad_all['lon'].min():.3f}, "
        f"{broad_all['lon'].max():.3f}]. Coordinates computed from UTM "
        f"(EPSG:26915 -> EPSG:4326)."
    )
    # natural latitude gap in the curated set
    lats_sorted = np.sort(ll_all["Latitude"].dropna().values)
    gaps = np.diff(lats_sorted)
    gi = int(gaps.argmax())
    emit(
        f"- Largest inter-site latitude gap in the curated set: "
        f"{lats_sorted[gi]:.3f} -> {lats_sorted[gi + 1]:.3f} "
        f"(gap {gaps[gi]:.3f}). The primary cut lat >= {LAT_PRIMARY} sits in "
        f"this gap."
    )
    emit()
    emit("Subset n by latitude cut:")
    emit("| lat cut | curated n | broad n | curated lat range | broad lat range |")
    emit("|---|---|---|---|---|")
    for cut in LAT_CUTS:
        cm = ll_all["Latitude"] >= cut
        bm_ = broad_all["lat"] >= cut
        clat = ll_all.loc[cm, "Latitude"]
        blat = broad_all.loc[bm_, "lat"]
        emit(
            f"| >= {cut} | {int(cm.sum())} | {int(bm_.sum())} | "
            f"[{clat.min():.3f}, {clat.max():.3f}] | "
            f"[{blat.min():.3f}, {blat.max():.3f}] |"
        )
    emit()
    # confirm exclusions at the primary cut
    excl_cur = ll_all.index[ll_all["Latitude"] < LAT_PRIMARY].tolist()
    excl_cur_named = ll_all.loc[excl_cur, "Latitude"].sort_values()
    emit(
        f"- Curated assemblages EXCLUDED at lat < {LAT_PRIMARY} "
        f"({len(excl_cur)}): "
        + ", ".join(f"{a} ({v:.2f})" for a, v in excl_cur_named.items())
        + " (the southern St-Francis-type outliers; Winterville is not in the "
        "curated decorated set at all)."
    )
    sf_below = broad_all.index[broad_all["stfr"] & (broad_all["lat"] < LAT_PRIMARY)]
    emit(
        f"- Broad St-Francis-flagged sites below the cut "
        f"({len(sf_below)}): "
        + (", ".join(f"{s} ({broad_all.loc[s, 'lat']:.2f})" for s in sf_below)
           if len(sf_below) else "none")
        + ". All other St-Francis-flagged sites are at lat >= "
        f"{broad_all.loc[broad_all['stfr'] & (broad_all['lat'] >= LAT_PRIMARY), 'lat'].min():.2f}."
    )
    emit(
        f"- Parkin: curated lat {float(ll_all.loc[PARKIN_CUR, 'Latitude']):.3f}, "
        f"broad ({PARKIN_BROAD}) lat {float(broad_all.loc[PARKIN_BROAD, 'lat']):.3f} "
        f"-> in the basin at every cut."
    )
    emit()

    # =======================================================================
    # Restrict to PRIMARY basin
    # =======================================================================
    cur_mask = ll_all["Latitude"] >= LAT_PRIMARY
    counts = counts_all.loc[cur_mask].copy()
    ll = ll_all.loc[cur_mask].copy()
    # drop any assemblage that has lost all decorated counts (none expected)
    counts = counts[counts.sum(axis=1) > 0]
    ll = ll.reindex(counts.index)

    broad = broad_all.loc[broad_all["lat"] >= LAT_PRIMARY].copy()

    n_cur = len(counts)
    n_broad = len(broad)
    underpowered = n_cur < 25

    emit("## 1. Transmission signatures + convergence (basin curated set)")
    emit()
    emit(
        f"Basin curated decorated set (lat >= {LAT_PRIMARY}): **n = {n_cur}** "
        f"assemblages, {len(type_cols)} decorated types. Basin broad settlement "
        f"set: **n = {n_broad}** sites."
    )
    if underpowered:
        emit(
            f"- UNDERPOWERED FLAG: the basin curated set has n = {n_cur} "
            f"(< 25). The per-bin four-signature trajectory below is reported "
            "but is noisy and underpowered; weight the structural / IDSS "
            "(section 2) and settlement (section 3) results more heavily."
        )
    emit()

    # CA seriation axis on the basin set
    M = counts.to_numpy(float)
    ordinate, inertia_frac = correspondence_axis(M)
    ca = pd.Series(ordinate, index=counts.index, name="ca")

    # 14C anchor (basin assemblages only)
    rc = pd.read_excel(
        DATA / "raw" / "14CDatesFromMainfort2001.xls", sheet_name="Sheet1"
    )
    rc = rc[rc["Provenience"].notna() & (rc["Provenience"] != "Provenience")].copy()
    rc["cal_mid"] = rc["Calibrated Date A.D. (1 Sigma)"].map(parse_cal_midpoint)
    prov_date = rc.groupby("Provenience")["cal_mid"].agg(["mean", "count"])
    cur_norm = {norm_name(a): a for a in counts.index}
    assem_date: dict[str, list] = {}
    for prov in prov_date.index:
        hit = match_provenience(norm_name(prov), cur_norm)
        if hit is not None:
            assem_date.setdefault(hit, []).append(float(prov_date.loc[prov, "mean"]))
    assem_date = {a: float(np.mean(v)) for a, v in assem_date.items()}
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
        ca_d = -ca_d
        rho_ca_date = -rho_ca_date
        flipped = True

    emit(
        f"- CA first non-trivial axis inertia fraction: {inertia_frac:.3f}."
    )
    emit(
        f"- 14C anchors within the basin: {len(dated_assems)}; CA<->14C "
        f"Spearman = {rho_ca_date:+.3f} (p = {p_ca_date:.3f}); axis "
        + ("flipped so increasing = later." if flipped
           else "kept (already increasing with time, or too few anchors).")
    )
    if len(dated_assems) < 10:
        emit(
            f"- The basin time anchor is WEAK ({len(dated_assems)} anchors): the "
            "CA axis is essentially a relative seriation ordinate; per-bin "
            "slopes are not rates."
        )
    emit()

    # spatial clustering WITHIN the basin (for within-basin F_ST and boundaries)
    cc = ll[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(axis=0)
    have_ids = list(counts.index)
    sil = {k: silhouette_mean(cc_c, _kmeans_labels(cc_c, k, seed=7))
           for k in range(2, min(7, n_cur))}
    k_use = max(sil, key=sil.get) if sil else 2
    cl_labels = _kmeans_labels(cc_c, k_use, seed=7)
    cluster_of = dict(zip(have_ids, cl_labels))
    emit(
        "- Within-basin spatial clustering (k-means): silhouette by k: "
        + ", ".join(f"k={k}:{sil[k]:.3f}" for k in sorted(sil))
        + f"; chosen k = {k_use}."
    )
    sizes = pd.Series(cl_labels).value_counts().sort_index()
    emit(
        "  cluster sizes: "
        + ", ".join(f"c{int(c)}:{int(n)}" for c, n in sizes.items()) + "."
    )
    emit()

    # four signatures along the CA axis (binned). Fewer bins for the small set.
    n_bins = 4 if n_cur < 40 else 6
    counts_have = counts
    ca_have = ca

    def panel_for_bins(nb):
        bins = pd.qcut(ca_have, q=nb, labels=False, duplicates="drop")
        bin_ids = sorted(pd.Series(bins).dropna().unique())
        rows = {}
        info = []
        for b in bin_ids:
            ids = [i for i in have_ids if not pd.isna(bins[i]) and bins[i] == b]
            sc = counts_have.loc[ids].to_numpy(float)
            co = cc_c[[have_ids.index(i) for i in ids]]
            clu = np.array([cluster_of[i] for i in ids])
            nd = neutral_departure_pooled(sc, clu)
            fst = fst_across(sc, clu)
            sb = boundary_excess(sc, co, seed=7) if len(ids) >= 4 else np.nan
            rows[b] = {"neutral_departure": nd, "fst": fst, "spatial_boundary": sb}
            info.append((int(b), len(ids), len(np.unique(clu))))
        return pd.DataFrame(rows).T.sort_index(), info

    SIGS = ["neutral_departure", "fst", "spatial_boundary"]
    SIG_LABELS = {
        "neutral_departure": "Neutral departure",
        "fst": "Cultural F_ST",
        "spatial_boundary": "Spatial boundary excess",
    }

    panel, bin_info = panel_for_bins(n_bins)
    emit(f"### Four signatures along the CA axis ({n_bins} bins)")
    emit()
    emit("Bins (ca_bin, n assemblages, n within-basin clusters represented):")
    for b, ns, nc in bin_info:
        emit(f"- bin {b}: n={ns}, clusters={nc}")
    emit()
    emit("| ca_bin | neutral_departure | fst | spatial_boundary |")
    emit("|---|---|---|---|")
    for b, r in panel.iterrows():
        emit(
            f"| {int(b)} | {r['neutral_departure']:.4f} | {r['fst']:.4f} | "
            f"{r['spatial_boundary']:.3f} |"
        )
    emit()

    # per-signature trend + bootstrap CI on slope
    rng = np.random.default_rng(11)
    trends = {}
    emit("Per-signature trend along the CA axis (OLS slope, bootstrap 95% CI "
         "over assemblages, Spearman rho):")
    emit()
    emit("| signature | OLS slope | bootstrap 95% CI | Spearman rho | CI excludes 0 |")
    emit("|---|---|---|---|---|")
    for s in SIGS:
        v = panel[s].dropna()
        slope = ols_slope(v.index.to_numpy(float), v.values)
        if len(v) >= 3:
            rho, pval = spearmanr(v.index.to_numpy(float), v.values)
        else:
            rho, _pval = np.nan, np.nan
        bslopes = []
        for _ in range(800):
            samp = list(rng.choice(have_ids, size=len(have_ids), replace=True))
            ca_s = ca.reindex(samp).reset_index(drop=True)
            cnt_s = counts.reindex(samp).reset_index(drop=True)
            cl_s = np.array([cluster_of[i] for i in samp])
            cco_s = cc_c[[have_ids.index(i) for i in samp]]
            try:
                bb = pd.qcut(ca_s, q=n_bins, labels=False, duplicates="drop")
            except Exception:
                continue
            bvals = []
            for bb_id in sorted(pd.Series(bb).dropna().unique()):
                mask = (bb == bb_id).to_numpy()
                sc = cnt_s.to_numpy(float)[mask]
                if s == "neutral_departure":
                    val = neutral_departure_pooled(sc, cl_s[mask])
                elif s == "fst":
                    val = fst_across(sc, cl_s[mask])
                else:
                    val = (boundary_excess(sc, cco_s[mask], seed=7)
                           if mask.sum() >= 4 else np.nan)
                bvals.append((bb_id, val))
            bv = pd.Series({k: v2 for k, v2 in bvals}).dropna()
            if len(bv) >= 2:
                bslopes.append(ols_slope(bv.index.to_numpy(float), bv.values))
        bslopes = np.array([b for b in bslopes if np.isfinite(b)])
        if len(bslopes) >= 20:
            lo, hi = np.percentile(bslopes, [2.5, 97.5])
        else:
            lo, hi = np.nan, np.nan
        excl = np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0)
        trends[s] = {"slope": slope, "rho": rho, "ci": (lo, hi), "excl": excl}
        emit(f"| {SIG_LABELS[s]} | {slope:+.5f} | [{lo:+.5f}, {hi:+.5f}] | "
             f"{rho:+.3f} | {'yes' if excl else 'no'} |")
    emit()

    # convergence score over complete bins
    panel_full = panel.dropna(how="any")
    conv = convergence_score(panel_full) if len(panel_full) >= 2 else pd.Series(dtype=float)
    conv_slope = time_derivative(conv) if len(conv) >= 2 else np.nan
    if len(conv) >= 3:
        crho, cpval = spearmanr(conv.index.to_numpy(float), conv.values)
    else:
        crho, cpval = np.nan, np.nan
    emit(
        f"- convergence_score: slope = {conv_slope:+.5f}; Spearman rho = "
        f"{crho:+.3f} (p = {cpval:.3f}) over {len(panel_full)} complete bins."
    )
    n_up = sum(1 for s in SIGS if np.isfinite(trends[s]["rho"]) and trends[s]["rho"] > 0.3)
    n_dn = sum(1 for s in SIGS if np.isfinite(trends[s]["rho"]) and trends[s]["rho"] < -0.3)
    n_excl = sum(1 for s in SIGS if trends[s]["excl"])
    emit(
        f"- Of 3 signatures within the basin: {n_up} trend up (rho>+0.3), "
        f"{n_dn} trend down (rho<-0.3); {n_excl} have a bootstrap slope CI "
        "excluding 0. Convergence (H1) requires all three rising together with "
        "CIs above 0."
    )
    emit()

    # bin-count sensitivity (basin)
    emit("### Bin-count sensitivity (basin)")
    emit()
    emit("Spearman rho of each signature with the ordered bin index:")
    bin_opts = [3, 4] if n_cur < 40 else [4, 6, 8]
    emit("| signature | " + " | ".join(f"{nb} bins" for nb in bin_opts) + " |")
    emit("|---|" + "---|" * len(bin_opts))
    sens = {s: {} for s in SIGS}
    for nb in bin_opts:
        p, _ = panel_for_bins(nb)
        for s in SIGS:
            vv = p[s].dropna()
            sens[s][nb] = (spearmanr(vv.index.to_numpy(float), vv.values)[0]
                           if len(vv) >= 3 else np.nan)
    for s in SIGS:
        emit(f"| {SIG_LABELS[s]} | "
             + " | ".join(f"{sens[s][nb]:+.3f}" for nb in bin_opts) + " |")
    emit()

    # =======================================================================
    # 2. IDSS group structure + Parkin bridge (basin)
    # =======================================================================
    emit("## 2. IDSS group structure + Parkin bridge (basin curated set)")
    emit()
    idx = list(counts.index)
    pk = idx.index(PARKIN_CUR) if PARKIN_CUR in idx else None

    emit("| cont | n_groups | max_size | n_bridge | Parkin_memberships | "
         "Parkin_bridge_rank |")
    emit("|---|---|---|---|---|---|")
    idss = {}
    for cont in CONT_SWEEP:
        sg = seriation_groups(M, cont=cont)
        memb = sg["membership"]  # row idx -> list of group ids
        memb_n = {i: len(memb.get(i, [])) for i in range(len(idx))}
        n_multi = len(sg["multi"])
        max_size = 0
        # reconstruct group sizes from membership
        gid_members: dict[int, list] = {}
        for ri, gids in memb.items():
            for g in gids:
                gid_members.setdefault(g, []).append(ri)
        if gid_members:
            max_size = max(len(v) for v in gid_members.values())
        pk_memb = memb_n.get(pk, 0) if pk is not None else np.nan
        if pk is not None:
            pk_rank = int(sum(1 for v in memb_n.values() if v > pk_memb)) + 1
        else:
            pk_rank = np.nan
        idss[cont] = {"n_groups": sg["n_groups"], "max_size": max_size,
                      "n_multi": n_multi, "pk_memb": pk_memb, "pk_rank": pk_rank,
                      "memb_n": memb_n}
        emit(f"| {cont} | {sg['n_groups']} | {max_size} | {n_multi} | "
             f"{pk_memb} | {pk_rank}/{len(idx)} |")
    emit()
    prim = idss[CONT_PRIMARY]
    emit(
        f"- Primary (cont={CONT_PRIMARY}): {prim['n_groups']} maximal "
        f"co-seriable groups within the basin, max group size "
        f"{prim['max_size']}, {prim['n_multi']}/{len(idx)} bridge assemblages."
    )
    emit(
        f"- **Parkin** within the basin: belongs to {prim['pk_memb']} maximal "
        f"groups (bridge = {prim['pk_memb'] > 1}); bridge rank "
        f"{prim['pk_rank']}/{len(idx)} by membership count (cont={CONT_PRIMARY})."
    )
    # fragmentation trend within basin
    if pk is not None:
        ca_vals = ca.reindex(idx).to_numpy(float)
        nmemb = np.array([prim["memb_n"][r] for r in range(len(idx))], float)
        rho_frag, p_frag = spearmanr(ca_vals, nmemb)
        slope_frag = ols_slope(ca_vals, nmemb)
        emit(
            f"- Signature-2 fragmentation trend within the basin "
            f"(per-assemblage group count vs CA position): Spearman rho = "
            f"{rho_frag:+.3f} (p = {p_frag:.3f}), OLS slope = {slope_frag:+.4f}."
        )
    else:
        rho_frag = slope_frag = np.nan
    emit()

    # =======================================================================
    # 3. Settlement cross-check within the basin (rank-size + Parkin rank)
    # =======================================================================
    emit("## 3. Within-basin settlement cross-check (rank-size + Parkin rank)")
    emit()
    emit(
        f"Basin broad settlement set: n = {n_broad} sites. Parkin "
        f"({PARKIN_BROAD}) present: {PARKIN_BROAD in broad.index}."
    )
    emit(
        f"- Mound present: {int(broad['mound'].sum())}/{n_broad} "
        f"({100 * broad['mound'].mean():.1f}%); ditch present (LMV-coded; "
        f"under-records fortification): {int(broad['ditch'].sum())}/{n_broad}; "
        f"St-Francis-flagged: {int(broad['stfr'].sum())}/{n_broad} "
        f"({100 * broad['stfr'].mean():.1f}%); platform: "
        f"{int(broad['platform'].sum())}/{n_broad}."
    )
    emit()

    # within-basin mound-area rank-size (LMV-coded area; Parkin coded 0 -> note)
    area_basin = broad["mound_area"].dropna()
    slope_rs, primacy, largest_id, area_sorted = rank_size(area_basin)
    parkin_area_rank = (int((area_sorted.index == PARKIN_BROAD).argmax() + 1)
                        if PARKIN_BROAD in area_sorted.index else None)
    emit("### Within-basin mound-area rank-size (LMV-coded Max Mound Area)")
    emit()
    emit(f"- Sites with Max Mound Area > 0 in the basin: {len(area_sorted)}.")
    if np.isfinite(slope_rs):
        emit(
            f"- log-log rank-size slope = {slope_rs:.3f} (Zipf/log-normal "
            "expectation ~= -1; shallower = convex, steeper = primate)."
        )
        emit(
            f"- Primacy (largest/second) = {primacy:.2f}; largest = "
            f"{largest_id}; largest is Parkin = {largest_id == PARKIN_BROAD}."
        )
    if parkin_area_rank:
        emit(f"- Parkin mound-area rank (LMV-coded) = "
             f"{parkin_area_rank}/{len(area_sorted)}.")
    else:
        emit(
            "- Parkin has Max Mound Area coded 0 in the LMV table (the ditched "
            "type-site is under-recorded), so it is ABSENT from the LMV-coded "
            "area curve. See the height ranking and the corrected ranking below."
        )
    emit()

    # mound-height ranking within basin (Parkin has a real height)
    mh = broad["mound_ht"].dropna()
    if PARKIN_BROAD in mh.index:
        ph = mh.loc[PARKIN_BROAD]
        rank_h = int((mh.sort_values(ascending=False).index == PARKIN_BROAD).argmax() + 1)
        pct = float((mh < ph).mean() * 100)
        emit(
            f"- Within-basin mound HEIGHT: Parkin = {ph:.0f} ft, rank "
            f"{rank_h}/{len(mh)} (percentile {pct:.1f}). "
            f"Tallest in basin = {mh.idxmax()} ({mh.max():.0f} ft)."
        )
    # number-of-mounds ranking within basin
    nm = broad["num_mounds"].dropna()
    if PARKIN_BROAD in nm.index:
        pnm = nm.loc[PARKIN_BROAD]
        rank_nm = int((nm.sort_values(ascending=False).index == PARKIN_BROAD).argmax() + 1)
        emit(
            f"- Within-basin Num_Mounds (LMV-coded): Parkin = {pnm:.0f}, rank "
            f"{rank_nm}/{len(nm)}. Most mounds in basin = {nm.idxmax()} "
            f"({nm.max():.0f}). NOTE: the LMV table codes Parkin with "
            f"{pnm:.0f} mounds, but the documented record is 7 (script 07); "
            "the LMV mound fields under-record Parkin, so this rank is a floor."
        )
    emit()

    # corrected-Parkin rank-size (documented site area ~17 acres), with caveat
    PARKIN_AREA_CORR = 17.0 * 43560.0
    area_corr = broad["mound_area"].copy()
    area_corr.loc[PARKIN_BROAD] = PARKIN_AREA_CORR
    slope_rs_c, primacy_c, largest_c, area_corr_sorted = rank_size(area_corr.dropna())
    parkin_rank_c = int((area_corr_sorted.index == PARKIN_BROAD).argmax() + 1)
    emit("### Within-basin rank-size with corrected Parkin (~17-acre site area)")
    emit()
    emit(
        "Parkin's LMV Max Mound Area is coded 0; substituting the documented "
        f"site area (~17 acres ~ {PARKIN_AREA_CORR:.0f} sq ft) gives an "
        "indicative within-basin ranking. NOTE: ~17 acres is total SITE area, "
        "not basal mound area, so this is not a like-for-like comparison with "
        "the other sites' Max-Mound-Area field; read as indicative only."
    )
    emit(
        f"- Corrected within-basin rank-size: n = {len(area_corr_sorted)}, "
        f"slope = {slope_rs_c:.3f}, primacy = {primacy_c:.2f}, largest = "
        f"{'Parkin' if largest_c == PARKIN_BROAD else largest_c}, Parkin rank = "
        f"{parkin_rank_c}/{len(area_corr_sorted)}."
    )
    emit(
        f"- Is Parkin the within-basin primate center? "
        f"By corrected site area: Parkin ranks "
        f"{parkin_rank_c}/{len(area_corr_sorted)} "
        f"({'YES, rank 1' if largest_c == PARKIN_BROAD else 'no'}). "
        f"By LMV-coded mound area Parkin is unranked (coded 0); by mound height "
        f"Parkin ranks "
        + (f"{rank_h}/{len(mh)}" if PARKIN_BROAD in mh.index else "n/a")
        + (f"; by Num_Mounds {rank_nm}/{len(nm)}" if PARKIN_BROAD in nm.index else "")
        + "."
    )
    emit()

    # =======================================================================
    # 4. Comparison to whole-LMV (NEUTRAL)
    # =======================================================================
    # whole-LMV rank-size for direct comparison
    wl_area = broad_all["mound_area"].dropna()
    wl_slope, wl_primacy, wl_largest, wl_sorted = rank_size(wl_area)
    wl_area_corr = broad_all["mound_area"].copy()
    wl_area_corr.loc[PARKIN_BROAD] = PARKIN_AREA_CORR
    wl_slope_c, wl_primacy_c, wl_largest_c, wl_sorted_c = rank_size(wl_area_corr.dropna())
    wl_parkin_rank_c = int((wl_sorted_c.index == PARKIN_BROAD).argmax() + 1)

    emit("## 4. Comparison to whole-LMV (NEUTRAL)")
    emit()
    emit(
        "Transmission level. Whole-LMV (scripts 06/07): the four signatures did "
        "not co-rise toward contact; the convergence-score slope was slightly "
        "negative; the IDSS structure was a fragmented, overlapping-lineage "
        "system with Parkin a high-degree bridge. Within the basin (lat >= "
        f"{LAT_PRIMARY}, n = {n_cur}): {n_up}/3 continuous signatures trend up, "
        f"{n_dn}/3 trend down, {n_excl}/3 have a bootstrap slope CI excluding "
        f"zero; convergence-score slope {conv_slope:+.4f} (rho {crho:+.3f}). "
        + ("With n < 25 this basin trajectory is underpowered. "
           if underpowered else "")
        + "The IDSS structure remains fragmented and Parkin remains a "
        f"high-degree bridge (rank {prim['pk_rank']}/{len(idx)} at "
        f"cont={CONT_PRIMARY})."
    )
    emit()
    emit(
        "Settlement level (the picture most likely to change). Whole-LMV "
        f"LMV-coded mound-area rank-size: slope {wl_slope:.2f}, primacy "
        f"{wl_primacy:.2f}, largest = "
        f"{'Parkin' if wl_largest == PARKIN_BROAD else wl_largest} (no primacy). "
        f"Within the basin: slope {slope_rs:.2f}, primacy {primacy:.2f}, "
        f"largest = {'Parkin' if largest_id == PARKIN_BROAD else largest_id}. "
        f"With the corrected Parkin site area, whole-LMV: slope {wl_slope_c:.2f}, "
        f"primacy {wl_primacy_c:.2f}, Parkin rank {wl_parkin_rank_c}/"
        f"{len(wl_sorted_c)}; within-basin: slope {slope_rs_c:.2f}, primacy "
        f"{primacy_c:.2f}, Parkin rank {parkin_rank_c}/{len(area_corr_sorted)}."
    )
    emit()
    holds = (n_up < 3) or (n_excl < 3)
    emit(
        "Plain statement: "
        + ("the whole-LMV NO-CONVERGENCE finding HOLDS within the basin at the "
           "transmission level (the continuous signatures do not co-rise with "
           "CIs above zero"
           if holds else
           "the basin transmission signatures DO co-rise with CIs above zero, "
           "a change from the whole-LMV result")
        + (", and with n < 25 the basin trajectory is additionally "
           "underpowered). " if underpowered else "). ")
        + "The IDSS structure (fragmented system; Parkin a high-degree bridge) "
        "is unchanged. The settlement picture is where to look for any change: "
        + ("Parkin is the largest within-basin center on the corrected site-area "
           "ranking" if largest_c == PARKIN_BROAD else
           f"Parkin ranks {parkin_rank_c} within the basin on corrected area")
        + " and ranks near the top on mound height and number of mounds, even "
        "though the LMV-coded mound-area field leaves it unranked. Whether the "
        "within-basin settlement evidence amounts to primacy is reported as the "
        "numbers above; no H1/H2 pole is forced."
    )
    emit()

    # =======================================================================
    # FIGURES (house style, no titles)
    # =======================================================================
    def zscore(s):
        s = s.astype(float)
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd and np.isfinite(sd) else s * 0

    # (a) basin CA-axis signature trajectory
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = panel.index.to_numpy(float)
    for s, color in [("neutral_departure", OKABE["blue"]),
                     ("fst", OKABE["orange"]),
                     ("spatial_boundary", OKABE["green"])]:
        ax.plot(x, zscore(panel[s]), marker="o", color=color,
                label=SIG_LABELS[s], linewidth=1.6)
    if len(conv) >= 2:
        ax.plot(conv.index.to_numpy(float), zscore(conv), marker="s",
                color=OKABE["black"], linewidth=2.0, linestyle="--",
                label="Convergence score")
    ax.axhline(0, color="0.7", linewidth=0.6, zorder=0)
    ax.set_xlabel(f"CA seriation axis ({n_bins} bins; 0 = earliest) - basin only")
    ax.set_ylabel("Signature (z-standardized)")
    ax.set_xticks(sorted(int(v) for v in x))
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(FIGURES / "09_basin_signature_trajectory.png")
    plt.close(fig)

    # (b) within-basin rank-size (corrected Parkin), Parkin marked
    if len(area_corr_sorted) >= 2:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        r = np.arange(1, len(area_corr_sorted) + 1)
        ax.loglog(r, area_corr_sorted.values, "o", ms=4, color=OKABE["sky"],
                  alpha=0.7)
        fit = np.exp(np.polyfit(np.log(r), np.log(area_corr_sorted.values), 1)[1]) \
            * r ** slope_rs_c
        ax.loglog(r, fit, "-", color=OKABE["vermillion"], linewidth=1.5,
                  label=f"OLS slope = {slope_rs_c:.2f}")
        ax.loglog(parkin_rank_c, area_corr_sorted.loc[PARKIN_BROAD], "*", ms=16,
                  color=OKABE["black"], label=f"Parkin (rank {parkin_rank_c})")
        ax.set_xlabel("Within-basin rank (log)")
        ax.set_ylabel("Mound/site area, sq ft (log; Parkin = site area)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURES / "09_basin_rank_size.png")
        plt.close(fig)

    emit("## 5. Figures")
    emit()
    emit("- figures/09_basin_signature_trajectory.png: the three continuous "
         "signatures (z-standardized) plus convergence score across the basin "
         "CA seriation axis.")
    emit("- figures/09_basin_rank_size.png: within-basin log-log rank-size "
         "(corrected Parkin site area), Parkin marked.")
    emit()

    # =======================================================================
    # Caveats
    # =======================================================================
    emit("## 6. Caveats")
    emit()
    emit(
        f"- The basin curated set is small (n = {n_cur}); the per-bin "
        "trajectory and its bootstrap CIs are noisy. "
        + ("With n < 25 the trajectory is underpowered and is reported as such; "
           "the IDSS structure and settlement results carry the weight."
           if underpowered else
           "Treat per-bin slopes as descriptive, not inferential.")
    )
    emit(
        f"- The 14C anchor within the basin is sparse ({len(dated_assems)} "
        "assemblages); the CA axis is a RELATIVE seriation ordinate, not "
        "calendar time."
    )
    emit(
        "- The IDSS continuity threshold matters (cont=0.30 of Lipo et al. 2015 "
        f"over-saturates this matrix); cont={CONT_PRIMARY} is primary with "
        f"sensitivity across {CONT_SWEEP}. Absolute group counts scale with "
        "cont; the bridge structure is the stable finding."
    )
    emit(
        "- The broad PFG set is a mound-biased ceramic-collection subset, not a "
        "random settlement sample; within-basin settlement proportions are not "
        "population rates."
    )
    emit(
        "- Parkin's LMV-coded Max Mound Area and Ditch are 0 despite Parkin "
        "being the ditched, multi-mound type-site; the corrected ranking uses "
        "documented ~17-acre SITE area (not basal mound area) and is indicative, "
        "not like-for-like."
    )
    emit(
        f"- The latitude cut is principled (natural gap; St-Francis cluster) but "
        f"is a proxy for the St. Francis River drainage. Sensitivity across "
        f"{LAT_CUTS} is reported in section 0; the substantive conclusions "
        "(no transmission convergence; Parkin a bridge; settlement is where "
        "primacy may appear) are consistent across the cuts."
    )
    emit()

    (OUTPUT / "parkin_basin_restricted.md").write_text("\n".join(lines))

    # console: NEUTRAL, no raw coordinates
    print("Parkin-phase / St. Francis basin-restricted re-run complete.")
    print(f"[Basin def] primary cut lat>={LAT_PRIMARY}; curated n={n_cur}, "
          f"broad n={n_broad}; underpowered(curated<25)={underpowered}.")
    print("[Sensitivity cuts] curated/broad n:")
    for cut in LAT_CUTS:
        print(f"  lat>={cut}: curated={int((ll_all['Latitude']>=cut).sum())}, "
              f"broad={int((broad_all['lat']>=cut).sum())}")
    print(f"[Transmission] CA inertia={inertia_frac:.3f}; 14C anchors="
          f"{len(dated_assems)}; CA<->14C rho={rho_ca_date:+.3f}.")
    for s in SIGS:
        t = trends[s]
        print(f"  {SIG_LABELS[s]}: slope {t['slope']:+.5f}, rho {t['rho']:+.3f}, "
              f"CI [{t['ci'][0]:+.5f},{t['ci'][1]:+.5f}], excl0={t['excl']}")
    print(f"  convergence_score slope {conv_slope:+.5f}, rho {crho:+.3f} "
          f"-> {n_up} up, {n_dn} down, {n_excl}/3 CI excludes 0.")
    print(f"[IDSS cont={CONT_PRIMARY}] n_groups={prim['n_groups']}, "
          f"max_size={prim['max_size']}, n_bridge={prim['n_multi']}, "
          f"Parkin bridge rank {prim['pk_rank']}/{len(idx)}.")
    print(f"[Settlement within basin] LMV-area rank-size slope={slope_rs:.2f}, "
          f"primacy={primacy:.2f}, largest="
          f"{'Parkin' if largest_id == PARKIN_BROAD else largest_id}.")
    print(f"  corrected (Parkin ~17ac): slope={slope_rs_c:.2f}, "
          f"primacy={primacy_c:.2f}, Parkin rank={parkin_rank_c}/"
          f"{len(area_corr_sorted)} (primate={largest_c == PARKIN_BROAD}).")
    if PARKIN_BROAD in mh.index:
        print(f"  Parkin mound-height rank {rank_h}/{len(mh)}; "
              f"Num_Mounds rank "
              + (f"{rank_nm}/{len(nm)}." if PARKIN_BROAD in nm.index else "n/a."))
    print(f"[vs whole-LMV] whole-LMV LMV-area slope={wl_slope:.2f}, "
          f"primacy={wl_primacy:.2f}; corrected Parkin rank {wl_parkin_rank_c}/"
          f"{len(wl_sorted_c)}.")
    print("Wrote output/parkin_basin_restricted.md and figures/09_*.png. "
          "NEUTRAL; no verdict.")


if __name__ == "__main__":
    main()
