"""Convergence scoring across the four transmission signatures.

Standardizes each signature across the ordinal sequence and averages them, so
convergence measures the co-movement of all four rather than any single one, and
returns the time-ordered slope used to test whether the signatures rise
together.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

def convergence_score(panel: pd.DataFrame) -> pd.Series:
    """Mean of column-standardized signatures: high where all four are jointly high.
    Convergence = co-movement; z-score each signature then average across the four so
    a single extreme signature cannot dominate."""
    z = (panel - panel.mean()) / panel.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1)

def time_derivative(series: pd.Series) -> float:
    """OLS slope of the time-ordered series (the H1 rising/accelerating test)."""
    x = np.asarray(series.index, float); y = np.asarray(series.values, float)
    x = x - x.mean()
    return float((x @ (y - y.mean())) / (x @ x))
