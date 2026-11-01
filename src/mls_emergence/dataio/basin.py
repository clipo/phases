"""The St. Francis basin analysis set: loader, membership, distance, fingerprint.

One definition of each (rule 6), moved here on 2026-10-16 (pilot step 2 of the
latent-group probe; inventory rows 3 and 7 of the probe spec). Before that
date the loader lived in a figure script, `analyses/make_figures.py`
(`_load_curated`, `_basin_members`, `DECORATED_TYPES`), and the distance
function existed in three copies (`analyses/23_phases_vs_spatial_drift.geo_km`,
`analyses/28_macro_boundary.geo_km` delegating to it, and
`analyses/74_phase_partition_test.km_xy`). The bodies are unchanged except
for the domain guards noted below; `tests/dataio/test_basin.py` pins their
output against the copies the Codex audit saved from the pre-move loader
(`codex/basin_counts.csv`, `codex/basin_coords.csv`, written by
`codex/audit_checks.py` at commit 9fb4f5d), which is the witness that the move
changed nothing, and checks `geo_km` against a great-circle computation that
shares no code with it.

Re-derivation note (rule 19). The loader reads the corrected analysis matrix
(`scripts/build_analysis_matrix.py`), keeps the ten decorated classes in the
fixed order below, drops rows with no decorated sherds, joins coordinates from
`data/raw/mainfort-pfg-cplXY.txt` through `read_assemblage_xy` (datum
corrections applied there), and restricts to the committed member list
`data/processed/basin_members_curated.txt`, which `analyses/16_basin_membership.py`
writes from the phase rule (Phillips's Parkin, Nodena, Kent or Walls area; at
least 75 decorated sherds; at least two members per area). The distance is an
equirectangular approximation about the mean latitude with 111.32 km per
degree; at 35 N over an 80 km basin its error against the great-circle
distance is below one percent (measured by the test).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from mls_emergence.dataio.coords import load_corrections, read_assemblage_xy
from mls_emergence.dataio.matrix import read_analysis_matrix

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data"
RAW_PFG = DATA / "raw" / "PFGData.xlsx"
COORD_FILE = DATA / "raw" / "mainfort-pfg-cplXY.txt"

# The ten decorated classes, in the order every matrix in this project uses.
DECORATED_TYPES = [
    "Parkin_Punctated", "Barton/Kent/MPI", "Painted", "Fortune_Noded",
    "Ranch_Incised", "Walls_Engraved", "Wallace_Incised", "Rhodes_Incised",
    "Vernon_Paul_Applique", "Hull_Engraved",
]

EARTH_KM_PER_DEG = 111.32


def basin_members(which: str = "curated") -> set:
    """Canonical basin membership written by analyses/16_basin_membership.py."""
    f = DATA / "processed" / f"basin_members_{which}.txt"
    return {ln.strip() for ln in f.read_text().splitlines() if ln.strip()}


def load_curated():
    """Counts (assemblage by decorated class) and WGS84 coordinates of the analysis set.

    Returns two DataFrames indexed by assemblage name in the same order:
    `counts` with the ten classes of DECORATED_TYPES, and `coords_ll` with
    columns Latitude and Longitude.
    """
    cur = read_analysis_matrix().dropna(subset=["Assemblages"])
    cur["Assemblages"] = cur["Assemblages"].astype(str).str.strip()
    cur = cur.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    type_cols = [c for c in DECORATED_TYPES if c in cur.columns]
    counts = cur[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    row_tot = counts.sum(axis=1)
    counts = counts[row_tot > 0]
    xy = read_assemblage_xy(COORD_FILE)
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    coords_ll = xy.reindex(counts.index)[["Latitude", "Longitude"]].apply(
        pd.to_numeric, errors="coerce"
    )
    members = basin_members("curated")
    keep = [a for a in counts.index if a in members]
    counts = counts.loc[keep]
    coords_ll = coords_ll.loc[keep]
    return counts, coords_ll


def _check_latlon(coords, name: str) -> np.ndarray:
    """Domain guard (rule 5): a 2-D array of finite [Latitude, Longitude] rows."""
    c = np.asarray(coords, float)
    if c.ndim != 2 or c.shape[1] != 2:
        raise ValueError(f"{name} expects an (n, 2) array of [Latitude, Longitude]; got shape {c.shape}")
    if not np.isfinite(c).all():
        raise ValueError(f"{name}: coordinates contain NaN or infinite values")
    if np.any(np.abs(c[:, 0]) > 90) or np.any(np.abs(c[:, 1]) > 180):
        raise ValueError(f"{name}: a value is outside latitude [-90, 90] or longitude [-180, 180]")
    return c


def geo_km(coords: np.ndarray) -> np.ndarray:
    """Pairwise distance in km by the equirectangular approximation about the mean latitude.

    coords columns are [Latitude, Longitude]. Body as in
    analyses/23_phases_vs_spatial_drift.geo_km (2026-06-15); the domain guard
    is the one addition.
    """
    c = _check_latlon(coords, "geo_km")
    lat = c[:, 0]
    lon = c[:, 1]
    mlat = np.deg2rad(lat.mean())
    x = np.deg2rad(lon) * np.cos(mlat) * EARTH_KM_PER_DEG / np.deg2rad(1)
    y = np.deg2rad(lat) * EARTH_KM_PER_DEG / np.deg2rad(1)
    p = np.column_stack([x, y])
    d = p[:, None, :] - p[None, :, :]
    return np.sqrt((d ** 2).sum(-1))


def km_xy(xy: np.ndarray) -> np.ndarray:
    """Coordinates as centered kilometers, equirectangular about the mean latitude.

    Body as in analyses/74_phase_partition_test.km_xy (2026-09-21); the domain
    guard is the one addition. Euclidean distances between its rows equal
    geo_km to floating-point precision (pinned by test).
    """
    c = _check_latlon(xy, "km_xy")
    mlat = np.deg2rad(c[:, 0].mean())
    p = np.column_stack([c[:, 1] * np.cos(mlat) * EARTH_KM_PER_DEG,
                         c[:, 0] * EARTH_KM_PER_DEG])
    return p - p.mean(0)


def _sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def provisional_coordinates(names) -> list:
    """Assemblages whose coordinate correction is marked provisional.

    Read from the text of the corrections table (`load_corrections`),
    not from a hard-coded name, so the flag clears when the correction is
    finalized there. As of 2026-10-16 this is Cummins (moved 22.7 km on
    2026-10-01, pending PFG's section description).
    """
    corr = load_corrections()
    text = corr.astype(str).apply(lambda row: " ".join(row.values), axis=1)
    flagged = {str(i) for i, t in text.items() if "rovisional" in t}
    return [n for n in map(str, names) if n in flagged]


def data_fingerprint(counts: pd.DataFrame, coords: pd.DataFrame, raw_pfg: Path | None = RAW_PFG) -> dict:
    """Hashes and geometry of the analysis set, for every fit record (probe spec section 8).

    Hashes are SHA-256 over canonical bytes: counts as int64 in C order after
    checking every value is an integer; assemblage IDs and class names joined
    by newlines in their stored order; coordinates and the geo_km distance
    matrix as float64 bytes (exact bytes, so a last-bit difference between
    machines shows as a mismatch rather than being hidden by rounding; the
    geometry figures beside the hashes say whether such a mismatch matters).
    The raw workbook's file hash is included when the file is present.
    Geometry figures are recomputed here, never quoted.
    """
    if list(counts.index) != list(coords.index):
        raise ValueError("counts and coords must be indexed identically, in the same order")
    y = counts.to_numpy(float)
    if not np.all(np.isfinite(y)) or not np.all(y == np.rint(y)) or np.any(y < 0):
        raise ValueError("counts must be finite non-negative integers")
    y = np.ascontiguousarray(y.astype(np.int64))
    xy = _check_latlon(coords[["Latitude", "Longitude"]].to_numpy(float), "data_fingerprint")
    D = geo_km(xy)
    off = D[~np.eye(len(D), dtype=bool)]
    nn = np.where(np.eye(len(D), dtype=bool), np.inf, D).min(axis=1)
    out = {
        "n_assemblages": int(len(y)),
        "n_classes": int(y.shape[1]),
        "n_sherds": int(y.sum()),
        "sha256_counts": _sha256(y.tobytes()),
        "sha256_ids": _sha256("\n".join(map(str, counts.index)).encode()),
        "sha256_classes": _sha256("\n".join(map(str, counts.columns)).encode()),
        "sha256_coords": _sha256(np.ascontiguousarray(xy).tobytes()),
        "sha256_distance_km": _sha256(np.ascontiguousarray(D).tobytes()),
        "max_pairwise_km": float(D.max()),
        "median_pairwise_km": float(np.median(off)),
        "nearest_neighbor_min_km": float(nn.min()),
        "nearest_neighbor_median_km": float(np.median(nn)),
        "nearest_neighbor_max_km": float(nn.max()),
        "provisional_coordinates": provisional_coordinates(counts.index),
    }
    if raw_pfg is not None and Path(raw_pfg).is_file():
        out["sha256_raw_pfg_xlsx"] = _sha256(Path(raw_pfg).read_bytes())
    return out
