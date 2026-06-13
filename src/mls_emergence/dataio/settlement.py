from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

ZONE_SHEETS = ["Locations-Zone-15", "Locations-Zone-16"]

# Matches the canonical grid format: <digits>-<token>-<digits>
# with an optional compound suffix starting with '/'
_GRID_RE = re.compile(
    r"^(\d+)-([^-/,\s]+)-(\d+)(/.*)?$",
    re.IGNORECASE,
)


def normalize_grid(site_id: str) -> str:
    """Return the canonical base grid identifier for *site_id*.

    Transformations applied (in order):
    1. Strip leading/trailing whitespace and uppercase.
    2. Strip any compound-site suffix beginning with '/' (e.g. '/A&B',
       '/B,C,D,E').
    3. In a ``<digits>-<X>-<digits>`` grid, if the middle token *X* is the
       digit ``0``, replace it with the letter ``O`` (OCR/transcription fix).
       Numeric segments are never modified.

    If *site_id* does not match the grid pattern it is returned uppercased and
    stripped (no other changes).
    """
    s = site_id.strip().upper()
    m = _GRID_RE.match(s)
    if m is None:
        return s
    left, mid, right = m.group(1), m.group(2), m.group(3)
    # Only substitute 0->O in the middle (letter) segment.
    if mid == "0":
        mid = "O"
    return f"{left}-{mid}-{right}"


def load_lmv(path: str | Path) -> pd.DataFrame:
    """Load and concatenate both UTM zone sheets from LMVData.xlsx."""
    frames = [pd.read_excel(path, sheet_name=s) for s in ZONE_SHEETS]
    return pd.concat(frames, ignore_index=True)


def join_pfg_to_lmv(counts: pd.DataFrame, lmv: pd.DataFrame):
    """Attach coordinates, size, and period to PFG assemblages.

    The LMV 'Number' column contains the PFG grid designator (e.g. '10-P-1'),
    matching the index of *counts*. 'PFG Type' is a descriptive label, not an
    identifier.

    Both PFG site ids and LMV Number values are normalised with
    :func:`normalize_grid` before matching to handle digit-0/letter-O OCR
    substitutions and compound/sub-site suffixes (e.g. '12-N-3/A&B').

    Duplicated Number entries in LMV (sites recorded in both zone sheets or
    with multiple entries) are resolved by keeping the first occurrence so that
    reindex returns a unique mapping.

    Returns
    -------
    joined : pd.DataFrame
        *counts* with spatial/attribute columns appended.
    unmatched : list[str]
        Site numbers from *counts* that had no matching LMV record
        (Easting is NaN after the join).
    """
    key_col = "Number"
    lmv_keyed = lmv.dropna(subset=[key_col]).copy()
    lmv_keyed["_key"] = (
        lmv_keyed[key_col].astype(str).map(normalize_grid)
    )
    # Drop duplicates so reindex returns a unique 1-to-1 mapping.
    lmv_keyed = lmv_keyed.drop_duplicates(subset=["_key"], keep="first")
    lmv_keyed = lmv_keyed.set_index("_key")

    # Preserve original PFG site ids for reporting; use normalised keys for lookup.
    raw_site_ids = counts.index.astype(str).str.strip()
    norm_site_ids = raw_site_ids.map(normalize_grid)

    matched = lmv_keyed.reindex(norm_site_ids.values)

    cols = [c for c in ["Area", "Northing", "Easting", "Zone", "Type",
                        "Period", "Terminal Period"] if c in matched.columns]
    joined = counts.copy()
    for c in cols:
        joined[c] = matched[c].values

    if "Easting" in joined.columns:
        unmatched = sorted(raw_site_ids[joined["Easting"].isna()].tolist())
    else:
        unmatched = sorted(raw_site_ids.tolist())

    return joined, unmatched
