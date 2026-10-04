"""The one place the decorated-class analysis matrix is read from.

Six analyses used to read it themselves, four from a sheet of
`data/raw/mainfort-pfg-cpl.xlsx` and two from `data/raw/mainfort-pfg-cpl.csv`,
which is how results came to rest on two copies of one table with nothing
enforcing that they agreed. That raw table also turned out to be wrong: it
doubles Parkin Punctated in its Mainfort component and sums two tallies of the
same sherds at 19 assemblages (scripts/check_matrix_against_mainfort2003.py).

The matrix now comes from `data/processed/analysis_matrix.csv`, which
`scripts/build_analysis_matrix.py` derives from the raw sources under a stated
rule: one source per assemblage, never a sum. The raw files are left as they
were received.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
ANALYSIS_MATRIX = ROOT / "data" / "processed" / "analysis_matrix.csv"


def read_analysis_matrix(path: Path | None = None) -> pd.DataFrame:
    """The corrected matrix, one row per assemblage, `Assemblages` as a column.

    Raises if the file is absent rather than falling back to the raw table: a
    silent fallback is exactly how the defective matrix would come back.
    """
    p = Path(path) if path is not None else ANALYSIS_MATRIX
    if not p.is_file():
        raise FileNotFoundError(
            f"{p} is missing. Build it with `python scripts/build_analysis_matrix.py`; "
            "the raw workbook is not a substitute (it double-counts).")
    return pd.read_csv(p)


MAINFORT_MATRIX = ROOT / "data" / "processed" / "mainfort2003_matrix.csv"
MAINFORT_SITES = ROOT / "data" / "processed" / "mainfort2003_sites.csv"

# Mainfort's (2003) ten decorated types onto this project's ten analysis
# classes. He tallies Barton Incised and Kent Incised separately where the
# analysis pools them with Mound Place Incised, and Old Town Red and Nodena Red
# and White separately where the analysis pools them with Avenue Polychrome;
# he does not tally Wallace Incised or Hull Engraved at all, so those two
# columns are zero in his matrix and are carried as zeros rather than dropped,
# to keep the class vocabulary identical across the two analyses.
MAINFORT_TO_CLASS = {
    "Parkin_Punctated": ["Parkin_Punctated"],
    "Barton/Kent/MPI": ["Barton_Incised", "Kent_Incised"],
    "Painted": ["Old_Town_Red", "Nodena_Red_and_White"],
    "Fortune_Noded": ["Fortune_Noded"],
    "Ranch_Incised": ["Ranch_Incised"],
    "Walls_Engraved": ["Walls_Engraved"],
    "Wallace_Incised": [],
    "Rhodes_Incised": ["Rhodes_Incised"],
    "Vernon_Paul_Applique": ["Vernon_Paul_Applique"],
    "Hull_Engraved": [],
}


def read_mainfort_replication(min_decorated: int = 100):
    """Mainfort's (2003) Table 1 as an analysis matrix, for the replication.

    Returns (counts, coords, phases): counts indexed by site with the ten
    analysis classes as columns, coordinates as a [Latitude, Longitude] frame
    on the same index, and his own phase assignment per site. Sites with fewer
    than `min_decorated` decorated sherds are dropped, and the function says
    which. Plainware never enters (author ruling, 2026-09-21: the analysis is
    limited to stylistic types).
    """
    m = pd.read_csv(MAINFORT_MATRIX, comment="#")
    sites = pd.read_csv(MAINFORT_SITES, comment="#").set_index("site")
    counts = pd.DataFrame({cls: m[cols].sum(axis=1) if cols else 0.0
                           for cls, cols in MAINFORT_TO_CLASS.items()})
    counts.index = m["site"].astype(str)
    total = counts.sum(axis=1)
    dropped = counts.index[total < min_decorated].tolist()
    keep = counts.index[total >= min_decorated]
    counts = counts.loc[keep].astype(float)
    coords = sites.loc[keep, ["latitude", "longitude"]].rename(
        columns={"latitude": "Latitude", "longitude": "Longitude"}).astype(float)
    if coords.isna().any().any():
        raise ValueError(f"missing coordinates for {coords.index[coords.isna().any(axis=1)].tolist()}")
    phases = sites.loc[keep, "phase_mainfort2003"].astype(str)
    if dropped:
        print(f"Mainfort replication: {len(dropped)} sites under {min_decorated} decorated "
              f"sherds dropped: {', '.join(dropped)}")
    return counts, coords, phases
