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
    """Load the combined LMV site-location table (both UTM zones).

    Accepts the combined CSV export (as in the public release) or the original
    LMVData.xlsx, whose two UTM zone sheets are concatenated.
    """
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return apply_site_coordinate_corrections(pd.read_csv(path))
    frames = [pd.read_excel(path, sheet_name=s) for s in ZONE_SHEETS]
    return apply_site_coordinate_corrections(pd.concat(frames, ignore_index=True))


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


# ---------------------------------------------------------------------------
# Published corrections to the compiled settlement table
# ---------------------------------------------------------------------------
CORRECTIONS = "mound_height_corrections.csv"
SITE_CORRECTIONS = "settlement_coordinate_corrections.csv"


def load_height_corrections(path):
    """Published mound measurements that supersede LMVData-22March2006.

    The compiled table records one `Max Mound Height (ft)` per site with no
    per-site provenance, and for some sites it disagrees with the excavation
    literature. Rather than edit the compilation, corrections live in their own
    committed file, each row naming the publication it comes from, so the
    substitution is visible in the data rather than buried in a script.

    Live rows: Parkin (11-N-1) at 21.3 ft with seven mounds, from Morse
    (1981, 1990), replacing 23.0 ft and four mounds. That single row decides
    whether Parkin ranks first in the basin, so it is data, not a constant.

    Returns a frame indexed by normalized site number.
    """
    import pandas as pd
    from pathlib import Path

    df = pd.read_csv(Path(path))
    df["_key"] = df["site_number"].astype(str).map(normalize_grid)
    return df.drop_duplicates("_key").set_index("_key")


def apply_height_corrections(heights, corrections, index_is_site_id=True):
    """Substitute corrected heights into a site-indexed height series.

    `heights` is indexed by site id as the caller holds it; matching is on the
    normalized key so that 11-N-1 and 11-N-1/A meet. Returns a new series and
    the list of site ids actually changed, which the caller should report
    rather than assume (rule 1).
    """
    out = heights.copy()
    changed = []
    for site_id in list(out.index):
        key = normalize_grid(str(site_id)) if index_is_site_id else str(site_id)
        if key in corrections.index:
            new = float(corrections.loc[key, "max_mound_height_ft"])
            if out.loc[site_id] != new:
                changed.append((str(site_id), float(out.loc[site_id]), new))
            out.loc[site_id] = new
    return out, changed


def apply_site_coordinate_corrections(lmv, path=None, verbose: bool = True):
    """Replace settlement coordinates that the PFG site table contradicts.

    `LMVData` positions four basin sites well away from where Phillips, Ford
    and Griffin put them, in one case by 50 km. Each correction here is PFG's
    own UTM, kept only where PFG's stated section contains it and the
    compilation's coordinate does not (checked against the BLM PLSS cadastral
    service, 2026-09-19); 13-N-3 also agrees with the legal description the
    author supplied and with the mound the lidar shows there.

    These coordinates set basin membership and the settlement figure's site
    set. They are NOT the assemblage coordinates that the transmission analysis
    uses — those live in `mainfort-pfg-cplXY.txt` and are corrected in
    `dataio/coords.py`.
    """
    if path is None:
        path = Path(__file__).resolve().parents[3] / "data" / "raw" / SITE_CORRECTIONS
    path = Path(path)
    if not path.exists():
        return lmv
    corr = pd.read_csv(path)
    corr["_k"] = corr["site_number"].astype(str).map(normalize_grid)
    key = next((c for c in lmv.columns if "Number" in str(c)), None)
    if key is None:
        return lmv
    out = lmv.copy()
    keys = out[key].astype(str).map(lambda v: normalize_grid(v) if isinstance(v, str) else "")
    for _, row in corr.iterrows():
        sel = keys == row["_k"]
        if not sel.any():
            continue
        before = (float(out.loc[sel, "Easting"].iloc[0]),
                  float(out.loc[sel, "Northing"].iloc[0]))
        out.loc[sel, "Easting"] = float(row["easting"])
        out.loc[sel, "Northing"] = float(row["northing"])
        out.loc[sel, "Zone"] = int(row["zone"])
        if verbose:
            d = ((float(row["easting"]) - before[0]) ** 2
                 + (float(row["northing"]) - before[1]) ** 2) ** 0.5
            print(f"settlement coordinate correction: {row['site_number']} moved "
                  f"{d / 1000:.2f} km [{str(row['source'])[:60]}...]")
    return out
