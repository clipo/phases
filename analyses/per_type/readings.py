"""Reading rules of the regional-drift spec (section 5): accounting, bands, statuses, separation, pairing, run loading.

Shared by analyses/per_type/per_class_table.py (per-cell accounting) and
analyses/per_type/drift_sensitivity_figure.py (coverage sets and separation
tables). Descriptive summaries only: a band is the 10th to 90th percentile of
the finite values (NumPy's linear quantile rule, endpoints inclusive) and is
eligible when at least MIN_FINITE of the cell's records give a finite value.
Nothing here is a test statistic or a success level.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

MIN_FINITE = 45
STATUSES = ("absent", "fixed", "undefined")


def accounting(cell: pd.DataFrame, c: int) -> dict:
    """Counts for class c over a cell's records: total, finite F, absent, fixed, numerically lost from the target, undefined pool."""
    f = pd.to_numeric(cell[f"F_{c}"], errors="coerce").to_numpy(float)
    st = cell[f"status_{c}"].astype(str).to_numpy() if f"status_{c}" in cell else np.array([""] * len(cell))
    return {"n_total": int(len(cell)), "n_finite": int(np.isfinite(f).sum()), "absent": int((st == "absent").sum()),
            "fixed": int((st == "fixed").sum()),
            "lost": int(cell[f"lost_{c}"].astype(bool).sum()) if f"lost_{c}" in cell else 0,
            "undefined_pool": int((cell["reason"].astype(str) == "undefined_pool").sum()) if "reason" in cell else 0}


def band(values, min_finite: int = MIN_FINITE) -> dict:
    """10th to 90th percentile band and median of the finite values; eligible when n_finite >= min_finite."""
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float); v = v[np.isfinite(v)]
    n = int(len(v))
    if n == 0:
        return {"lo": np.nan, "hi": np.nan, "median": np.nan, "n_finite": 0, "eligible": False}
    return {"lo": float(np.quantile(v, .1)), "hi": float(np.quantile(v, .9)), "median": float(np.median(v)),
            "n_finite": n, "eligible": bool(n >= min_finite)}


def status(x, b: dict) -> str:
    """inside, outside, or unavailable (ineligible band or undefined value)."""
    if not b["eligible"] or x is None or not np.isfinite(x):
        return "unavailable"
    return "inside" if b["lo"] <= x <= b["hi"] else "outside"


def status3(x, b: dict) -> str:
    """below, inside, above, or unavailable (rim R1 spec 4.1: only 'above' is read as more variation than the baseline)."""
    if not b["eligible"] or x is None or not np.isfinite(x):
        return "unavailable"
    return "below" if x < b["lo"] else "above" if x > b["hi"] else "inside"


def separation(planted, b: dict, min_finite: int = MIN_FINITE) -> dict:
    """Counts of planted records outside and inside an eligible neutral band, undefined records beside them."""
    v = pd.to_numeric(pd.Series(planted), errors="coerce").to_numpy(float); fin = np.isfinite(v)
    n_undefined = int((~fin).sum())
    if not b["eligible"]:
        return {"available": False, "n_out": None, "n_in": None, "n_undefined": n_undefined, "text": "unavailable",
                "median": np.nan, "median_eligible": False}
    vf = v[fin]; out = (vf < b["lo"]) | (vf > b["hi"])
    n_out, n_in = int(out.sum()), int((~out).sum())
    return {"available": True, "n_out": n_out, "n_in": n_in, "n_undefined": n_undefined,
            "text": f"{n_out} of {n_out + n_in} ({n_undefined} undefined)",
            "median": float(np.median(vf)) if len(vf) else np.nan, "median_eligible": bool(len(vf) >= min_finite)}


def paired_changes(planted, neutral) -> dict:
    """Per-record paired change, counted only where both statistics are finite."""
    a = pd.to_numeric(pd.Series(planted), errors="coerce").to_numpy(float)
    b = pd.to_numeric(pd.Series(neutral), errors="coerce").to_numpy(float)
    if a.shape != b.shape:
        raise ValueError("paired series must have the same length and order")
    ok = np.isfinite(a) & np.isfinite(b)
    return {"increased": int((a[ok] > b[ok]).sum()), "decreased": int((a[ok] < b[ok]).sum()),
            "tied": int((a[ok] == b[ok]).sum()), "unavailable": int((~ok).sum())}


KEY = ("point", "generator", "F_W", "record")
MARKER_KEY = ("marker_class", "cluster", "m")


def load_runs(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate run CSVs, refusing duplicate record keys (point, generator, F_W, record, and the marker columns where present)."""
    df = pd.concat(list(frames), ignore_index=True)
    key = [k for k in KEY] + [k for k in MARKER_KEY if k in df.columns]
    dup = df.duplicated(subset=key, keep=False)
    if dup.any():
        rows = df.loc[dup, key].drop_duplicates().head(5).to_dict("records")
        raise ValueError(f"duplicate record keys across runs: {rows}")
    return df


def _evaluate(statuses: list[str]) -> str:
    """hit when every class is inside; miss when any is outside; otherwise unavailable."""
    if any(s == "outside" for s in statuses):
        return "miss"
    if any(s == "unavailable" for s in statuses):
        return "unavailable"
    return "hit"


def first_m_label(neutral_band: dict, planted: list) -> str:
    """Spec 5: the first tested m whose eligible planted median lies outside the neutral band; 'none on this grid' when every eligible median is inside; 'unavailable' when the neutral band or every planted median is ineligible; 'not run' when nothing was planted.

    `planted` is a list of (m, median_eligible, median) in increasing m.
    """
    if not planted:
        return "not run"
    if not neutral_band["eligible"] or not any(e for _m, e, _v in planted):
        return "unavailable"
    for m, eligible, med in planted:
        if eligible and status(med, neutral_band) == "outside":
            return f"{m:g}"
    return "none on this grid"


def coverage_sets(statuses: dict, grid, names: dict | None = None) -> dict:
    """C, R, H and their availability sets over the grid (spec 5). statuses[a] = {"common": [...], "rare": [...]}.

    C holds at a when all common classes are inside; it is unavailable when any
    is unavailable and none is outside (an outside class is a miss whatever the
    others). R likewise. H = C and R, unavailable if either is.
    """
    out = {k: [] for k in ("C", "R", "H", "E_C", "E_R", "E_H")}
    label = (lambda group, i: names[group][i]) if names else (lambda group, i: f"{group} class {i}")
    out["unavailable_cells"] = [(a, label(group, i)) for a in grid for group in ("common", "rare")
                                for i, s in enumerate(statuses[a][group]) if s == "unavailable"]
    for a in grid:
        ec, er = _evaluate(statuses[a]["common"]), _evaluate(statuses[a]["rare"])
        if ec != "unavailable":
            out["E_C"].append(a)
        if ec == "hit":
            out["C"].append(a)
        if er != "unavailable":
            out["E_R"].append(a)
        if er == "hit":
            out["R"].append(a)
        if ec != "unavailable" and er != "unavailable":
            out["E_H"].append(a)
            if ec == "hit" and er == "hit":
                out["H"].append(a)
    return out


def _set(xs) -> str:
    return "{" + ", ".join(f"{x:g}" for x in xs) + "}"


def coverage_wording(sets: dict) -> str:
    """The pre-stated wording of spec section 5 for the four coverage cases."""
    if not sets["E_H"]:
        cells = "; ".join(f"F_W {a:g}, {name}" for a, name in sets.get("unavailable_cells", []))
        return (f"no amplitude on the grid is evaluable for both sets (ineligible cells: {cells or 'none'}; common classes "
                f"evaluable at {_set(sets['E_C'])}, rare classes at {_set(sets['E_R'])})")
    if sets["H"]:
        return f"inside the bands at F_W in {_set(sets['H'])}; first tested value {min(sets['H']):g}"
    if sets["C"] and sets["R"]:
        return (f"the common classes are inside the bands at {_set(sets['C'])} and the rare classes at {_set(sets['R'])}, "
                "at no common tested amplitude")
    return (f"not inside the bands at any evaluable tested amplitude (common classes evaluable at {_set(sets['E_C'])}, "
            f"inside at {_set(sets['C'])}; rare classes evaluable at {_set(sets['E_R'])}, inside at {_set(sets['R'])})")
