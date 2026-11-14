"""The two copies of the aggregated ceramic matrix must be the same matrix.

`data/raw/mainfort-pfg-cpl.csv` and the `pfg-cpl-mainfort` sheet of
`data/raw/mainfort-pfg-cpl.xlsx` both hold the decorated-class counts
aggregated from Phillips-Ford-Griffin, Mainfort (1996) and Lipo (2001). Two
entry points read different copies: `make_figures._load_curated` reads the csv
and feeds every basin figure and F_ST, while `07_refined_empirical` reads the
workbook and feeds the Bayesian fit. An edit to one and not the other would put
different results on different data with nothing to say so.

Rule 3: the discriminating probe is the cell-by-cell comparison, not the
shape. Two files of the same shape and the same total can still differ.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "data" / "raw" / "mainfort-pfg-cpl.csv"
XLSX = ROOT / "data" / "raw" / "mainfort-pfg-cpl.xlsx"


@pytest.mark.data
def test_csv_and_workbook_hold_the_same_counts():
    c = pd.read_csv(CSV)
    x = pd.read_excel(XLSX, sheet_name="pfg-cpl-mainfort")
    x = x.dropna(subset=[x.columns[0]])
    ck, xk = c.columns[0], x.columns[0]
    names_c = [str(v).strip() for v in c[ck]]
    names_x = [str(v).strip() for v in x[xk]]
    assert names_c == names_x, "assemblage rows differ, or differ in order"
    assert list(c.columns[1:]) == list(x.columns[1:]), "class columns differ"

    a = c.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
    b = x.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
    differ = np.argwhere(a != b)
    assert len(differ) == 0, (
        "cells differ between the csv and the workbook: "
        + "; ".join(f"{names_c[i]} / {c.columns[1 + j]}: csv {a[i, j]:g}, xlsx {b[i, j]:g}"
                    for i, j in differ[:8]))
    # Pinned so a wholesale replacement of BOTH files is also noticed.
    assert a.shape == (55, 10)
    assert int(a.sum()) == 38101
