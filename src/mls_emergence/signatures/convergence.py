"""Combining the four signatures into one convergence readout.

Provides ``convergence_score`` (the mean of column-standardized signatures in a
time-ordered panel) and ``time_derivative`` (the least-squares slope of a series
on its index). Used by ``transmission.model``, ``validation.harness`` and
analyses 05, 06 and 09.
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
