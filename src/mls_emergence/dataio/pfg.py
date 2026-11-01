"""Readers for the Phillips, Ford and Griffin (PFG) ceramic tables.

Provides ``load_pfg_counts`` (assemblage x type sherd counts from the
``SherdData`` sheet of ``data/raw/PFGData.xlsx``, or a CSV export of it) and
``load_pfg_attributes`` (the type attribute table on ``Sheet2``). Used by the
empirical analyses (05, 06, 07, 09, 16, 17, 27, 67), ``make_figures.py`` and
``make_map.py``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

NON_TYPE_COLS = ["Site Name", "Site Number", "Sherd Total"]

def load_pfg_counts(path: str | Path) -> pd.DataFrame:
    """Assemblage x type count matrix, indexed by PFG Site Number."""
    path = Path(path)
    raw = pd.read_csv(path) if path.suffix.lower() == ".csv" else pd.read_excel(path, sheet_name="SherdData")
    # Strip whitespace from column names (e.g., 'Old Town Red ' has trailing space)
    raw.columns = raw.columns.str.strip()
    raw = raw.dropna(subset=["Site Number"])
    type_cols = [c for c in raw.columns if c not in NON_TYPE_COLS]
    df = raw.set_index("Site Number")[type_cols].fillna(0).astype(int)
    df = df[df.sum(axis=1) > 0]
    df.attrs["names"] = raw.set_index("Site Number")["Site Name"].to_dict()
    return df

def load_pfg_attributes(path: str | Path) -> pd.DataFrame:
    """Type attribute table (temper, surface treatment, decoration)."""
    path = Path(path)
    return pd.read_csv(path) if path.suffix.lower() == ".csv" else pd.read_excel(path, sheet_name="Sheet2")
