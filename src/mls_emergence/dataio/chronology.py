"""Load period and chronology attributes from the LMV location table.

Provides the assemblage period assignments used alongside the seriation-derived
ordinal chronology. See data/README.md for provenance.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

def load_period_assignments(path: str | Path) -> pd.DataFrame:
    """Load the combined LMV location table and keep identity + period columns."""
    df = pd.read_csv(path)
    # Strip whitespace from column names to guard against source variants
    df.columns = df.columns.str.strip()
    keep = ["Number", "StateNum", "Name", "Area", "Type", "PFG Type",
            "Period", "Initial Period", "Terminal Period"]
    keep = [c for c in keep if c in df.columns]
    return df[keep]
