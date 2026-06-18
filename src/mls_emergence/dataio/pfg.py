"""Load the Phillips-Ford-Griffin (1951) decorated-ceramic tables.

Reads the assemblage-by-type sherd-count matrix and the type-attribute table
(compiled by Lipo 2001) into pandas frames keyed by Lower Mississippi Survey
site number. See data/README.md for provenance.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

NON_TYPE_COLS = ["Site Name", "Site Number", "Sherd Total"]

def load_pfg_counts(path: str | Path) -> pd.DataFrame:
    """Assemblage x type count matrix, indexed by PFG Site Number."""
    raw = pd.read_csv(path)
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
    return pd.read_csv(path)
