from __future__ import annotations
from pathlib import Path
import pandas as pd

PERIOD_SHEETS = ["Locations-Zone-15", "Locations-Zone-16"]

def load_period_assignments(path: str | Path) -> pd.DataFrame:
    """Concatenate the zone sheets and keep identity + period columns."""
    frames = [pd.read_excel(path, sheet_name=s) for s in PERIOD_SHEETS]
    df = pd.concat(frames, ignore_index=True)
    # Strip whitespace from column names to guard against source variants
    df.columns = df.columns.str.strip()
    keep = ["Number", "StateNum", "Name", "Area", "Type", "PFG Type",
            "Period", "Initial Period", "Terminal Period"]
    keep = [c for c in keep if c in df.columns]
    return df[keep]
