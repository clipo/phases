"""The resemblance stop rule and the restricted screen (spec 3.3; R2-6).

Summaries are the operating-point screen's own: `diversity()` of
analyses/47_revision_analysis.py (mean within-assemblage Gini-Simpson "hs",
mean class richness "rich", pooled Gini-Simpson "ht") and the plug-in
`cultural_fst` on counts aggregated by the saved spatial clusters
("fst_clusters"). The screen compares MEANS across seeds with the observed
values at a relative tolerance (0.10, the screen's); the loss is the sum of
squared relative errors; the S0-stress cell is the matched cell with the
smallest F_ST median, ties by the smallest loss, then by grid order.
"""
from __future__ import annotations

import numpy as np

from mls_emergence.signatures.variance import cultural_fst

SUMMARIES = ("hs", "rich", "ht")      # the three screened on
TOLERANCE = 0.10


def summaries(counts: np.ndarray, labels: np.ndarray, diversity) -> dict:
    m = np.asarray(counts, float)
    d = dict(diversity(m))
    labels = np.asarray(labels)
    d["fst_clusters"] = float(cultural_fst(np.array([m[labels == g].sum(0) for g in np.unique(labels)])))
    return d


def relative_errors(mean: dict, observed: dict) -> dict:
    return {k: abs(mean[k] - observed[k]) / observed[k] for k in SUMMARIES}


def screen_cell(mean: dict, observed: dict, tolerance: float = TOLERANCE) -> dict:
    rel = relative_errors(mean, observed)
    return {"matched": bool(max(rel.values()) <= tolerance), "loss": float(sum(v ** 2 for v in rel.values())),
            **{f"rel_{k}": float(v) for k, v in rel.items()}}


def select_stress_cell(rows: list[dict]) -> dict | None:
    """Among matched rows: smallest fst_median, ties by smallest loss, then by grid order (list order)."""
    matched = [(i, r) for i, r in enumerate(rows) if r.get("matched")]
    if not matched:
        return None
    i, r = min(matched, key=lambda t: (t[1]["fst_median"], t[1]["loss"], t[0]))
    return dict(r, grid_index=i)


def stop_rule(mean: dict, observed: dict, tolerance: float = TOLERANCE) -> dict:
    """The operating point passes if every diversity summary's mean is within the tolerance."""
    c = screen_cell(mean, observed, tolerance)
    return {"pass": c["matched"], "tolerance": tolerance, **{k: v for k, v in c.items() if k.startswith("rel_")}}


def profile_f(rows: np.ndarray) -> np.ndarray:
    """Sampling-free per-class F of probability rows (regional-drift spec 3.3): var_i(r_ic, ddof 1) / (mean (1 - mean)).

    No binomial subtraction (the rows are probabilities, not counts, so the
    count-based `moment_f` is not applicable); NaN where the mean is 0 or 1.
    """
    r = np.asarray(rows, float)
    if r.ndim != 2 or r.shape[0] < 2:
        raise ValueError("profile_f needs at least two probability rows")
    if not np.isfinite(r).all() or (r < 0).any() or (r > 1).any() or not np.allclose(r.sum(1), 1.0, atol=1e-9, rtol=0):
        raise ValueError("profile_f takes probability rows: finite, in [0, 1], summing to 1")
    m = r.mean(0); v = r.var(0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = v / (m * (1 - m))
    out[(m <= 0) | (m >= 1)] = np.nan
    return out
