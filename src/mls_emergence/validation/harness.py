"""Blind validation harness: compute the four signatures over an ordinal axis,
run all mechanisms blind, and apply the convergence-discrimination rule.

The harness never branches on the mechanism name. It applies the identical
signature pipeline to every generator's output, so any discrimination it
achieves is a property of the criterion, not of the harness.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd

from mls_emergence.signatures import convergence, neutral, seriation, variance
from mls_emergence.signatures.assortativity import boundary_excess

SIGNATURE_COLUMNS = ["neutral_departure", "seriability", "fst", "spatial_boundary"]


def _within_conformity_departure(slice_counts: np.ndarray) -> float:
    """Mean within-group signed conformity departure, 1 - theta_f/theta_e.

    Positive => conformist (fewer effective types than Ewens neutrality predicts).
    Conformity is a within-group transmission bias; measuring it on the pooled
    assemblage conflates it with between-group divergence (pooling divergent
    groups flattens the union and cancels the signal). See mechanisms.py.
    """
    vals = []
    for row in slice_counts:
        tf = neutral.theta_f(row)
        te = neutral.theta_e(row)
        if np.isfinite(tf) and te > 0:
            vals.append(1.0 - tf / te)
    return float(np.mean(vals)) if vals else 0.0


def signatures_over_axis(slices, coords) -> pd.DataFrame:
    """Compute the four cultural-transmission signatures for each ordinal slice.

    Columns:
      neutral_departure  within-group conformity departure (Neiman/Ewens).
      seriability        a co-seriability score that RISES as the assemblage
                         approaches frequency-seriation order. Defined as the
                         negative of seriation.unimodality_violation, which falls
                         monotonically as groups differentiate into clean,
                         orderable frequency profiles. NOTE: the literal
                         seriation.n_seriation_groups proxy is unusable here:
                         with G=12 > its brute-force max_rows it saturates at G
                         for every slice (slope 0), and it is anti-oriented (more
                         groups means LESS structure, and genuine emergence
                         REDUCES the number of seriation groups). Using negative
                         violation fixes both the saturation and the orientation
                         so the convergence rule ("all four trend up") is
                         coherent. Documented deviation from the spec column name.
      fst                cultural F_ST among groups.
      spatial_boundary   boundary_excess: within-cluster minus between-cluster
                         mean similarity at matched distance. Distinguishes a
                         sharp bounded edge (large positive excess) from smooth
                         isolation-by-distance (excess near zero), which a plain
                         Mantel r cannot do (both give negative r).
    """
    rows = []
    for s in slices:
        s = np.asarray(s)
        rows.append(
            {
                "neutral_departure": _within_conformity_departure(s),
                "seriability": float(-seriation.unimodality_violation(s)),
                "fst": variance.cultural_fst(s),
                "spatial_boundary": boundary_excess(s, coords),
            }
        )
    df = pd.DataFrame(rows, columns=SIGNATURE_COLUMNS)
    df.index = range(len(df))
    return df


def run_blind(generators: dict[str, Callable], seed: int) -> dict[str, pd.DataFrame]:
    """Run each named generator and compute its signature panel, blind to name."""
    panels: dict[str, pd.DataFrame] = {}
    for name, gen in generators.items():
        slices, coords = gen(seed)
        panels[name] = signatures_over_axis(slices, coords)
    return panels


def discriminates(
    panels: dict[str, pd.DataFrame], deriv_threshold: float = 0.10
) -> dict[str, dict]:
    """Apply the convergence criterion: a mechanism is 'convergent' iff ALL FOUR
    signatures show a positive ordinal trend above ``deriv_threshold``.

    The four signatures live on very different scales (F_ST ~0-1 vs the BR
    boundary excess ~0-200), so a single raw-slope threshold is not comparable
    across them. Each signature is therefore standardized by its pooled spread
    ACROSS all mechanisms before the slope is taken, so ``deriv_threshold`` is in
    standard-deviation-of-the-signature per ordinal step. A signature that is
    flat for every mechanism (zero pooled spread) is treated as non-trending.

    Returns per-mechanism {trends: {col: {slope, slope_std, terminal}},
    convergent: bool}.
    """
    pooled = {
        col: pd.concat([p[col] for p in panels.values()], ignore_index=True)
        for col in SIGNATURE_COLUMNS
    }
    scale = {col: float(pooled[col].std(ddof=0)) for col in SIGNATURE_COLUMNS}

    result: dict[str, dict] = {}
    for name, panel in panels.items():
        trends: dict[str, dict] = {}
        for col in SIGNATURE_COLUMNS:
            raw_slope = convergence.time_derivative(panel[col])
            std_slope = raw_slope / scale[col] if scale[col] > 0 else 0.0
            trends[col] = {
                "slope": float(raw_slope),
                "slope_std": float(std_slope),
                "terminal": float(panel[col].iloc[-1]),
            }
        convergent = all(
            trends[col]["slope_std"] > deriv_threshold for col in SIGNATURE_COLUMNS
        )
        result[name] = {"trends": trends, "convergent": convergent}
    return result
