"""Reading rules R1 to R5 of the drift-extension spec (2026-11-14, draft 3, section 6).

One pure function, `read`, from one operating point's per-cell status table to
the sets E_L, L, E_S and S, the R1 lists, the R5 table by coupling and exactly
one of three statement forms. Statuses are those of `readings.status3`:
inside, above, below, unavailable. Every statement is worded as marginal
coverage: each named summary of the record lies inside its own band. Nothing
here says that a simulated record matches the record in several summaries at
once, and no count is a test statistic. Cells that share a window share
simulated histories, so counts of cells are not counts of independent
confirmations: `n_histories` gives the number of distinct histories behind a
set of cells.
"""
from __future__ import annotations

import pandas as pd

SPANS = ("Ld", "Ld/4", "1")
WINDOWS = ("8LR", "2LR", "32LR")
COUPLINGS = ("c0", "lat50N", "lat100N", "lat50S", "lat100S", "stepL", "stepS")
N_AMPLITUDES = 8
N_CELLS = len(SPANS) * len(WINDOWS) * len(COUPLINGS) * N_AMPLITUDES     # 504 per operating point
RECORDS_PER_CELL = 50
D0 = ("Ld", "8LR", "c0")                                                # the baseline design
F3 = ("F_0", "F_1", "F_2")                                              # the three common-class F values
SP = ("D_dist", "A_lat_0", "A_lat_1", "A_lat_2")                        # the four spatial readouts
RARE = tuple(f"occ_{k}" for k in range(3, 10))                          # occupancy of the seven rare classes
REQUIRED = F3 + SP + RARE
STATUSES = ("inside", "above", "below", "unavailable")
CELL = ["span", "window", "coupling", "F_W"]


def _order(cells) -> list[tuple]:
    return sorted(cells, key=lambda c: (SPANS.index(c[0]), WINDOWS.index(c[1]), COUPLINGS.index(c[2]), c[3]))


def _wide(table: pd.DataFrame) -> pd.DataFrame:
    """The table as cells by required statistics, after checking the grid, the statistics and the status words."""
    missing = [c for c in CELL + ["stat", "status"] if c not in table.columns]
    if missing:
        raise ValueError(f"status table lacks columns {missing}")
    t = table[table["stat"].isin(REQUIRED)]
    if t.duplicated(subset=CELL + ["stat"]).any():
        raise ValueError("status table has more than one row for a cell and statistic")
    bad = sorted(set(t["status"]) - set(STATUSES))
    if bad:
        raise ValueError(f"unknown status words {bad}; the rules read {STATUSES}")
    for col, allowed in (("span", SPANS), ("window", WINDOWS), ("coupling", COUPLINGS)):
        extra = sorted(set(t[col]) - set(allowed))
        if extra:
            raise ValueError(f"unknown {col} values {extra}")
    w = t.pivot(index=CELL, columns="stat", values="status")
    if len(w) != N_CELLS or t["F_W"].nunique() != N_AMPLITUDES:
        raise ValueError(f"the rules are applied to one operating point's {N_CELLS} cells; the table has {len(w)}")
    lacking = [s for s in REQUIRED if s not in w.columns or w[s].isna().any()]
    if lacking:
        raise ValueError(f"cells lack a status for {lacking}")
    return w


def read(table: pd.DataFrame) -> dict:
    """Sets, R1 lists, R2 and R3 counts, the R5 table and the statement form for one operating point."""
    w = _wide(table)
    f3, sp = w[list(F3)], w[list(SP)]
    ev_l = ~(f3 == "unavailable").any(axis=1)
    in_l = ev_l & (f3 == "inside").all(axis=1)
    ev_s = in_l & ~(sp == "unavailable").any(axis=1)
    in_s = ev_s & (sp == "inside").all(axis=1)
    any_above = ev_l & ~in_l & (f3 == "above").any(axis=1)
    below_only = ev_l & ~in_l & ~(f3 == "above").any(axis=1)        # none above, so at least one below
    cells = lambda mask: _order(list(w.index[mask]))
    out = {"n_cells": N_CELLS, "E_L": cells(ev_l), "L": cells(in_l), "E_S": cells(ev_s), "S": cells(in_s),
           "n_any_above": int(any_above.sum()), "n_below_only": int(below_only.sum()),
           "L_not_evaluable_for_space": cells(in_l & ~ev_s)}
    amps = sorted(set(w.index.get_level_values("F_W")))
    out["R1"] = {s: {st: [a for a in amps if w.loc[(*D0, a), s] == st] for st in STATUSES} for s in F3}
    coupling = w.index.get_level_values("coupling")
    out["R5"] = pd.DataFrame([{"coupling": c, "cells": int((coupling == c).sum()), "E_L": int((ev_l & (coupling == c)).sum()),
                               "L": int((in_l & (coupling == c)).sum()), "E_S": int((ev_s & (coupling == c)).sum()),
                               "S": int((in_s & (coupling == c)).sum())} for c in COUPLINGS])
    n_el, n_l, n_es, n_s = (len(out[k]) for k in ("E_L", "L", "E_S", "S"))
    if n_l == 0:
        out["form"] = "a"
        out["statement"] = ("At no evaluable tested cell do all three common-class F values of the record lie inside their bands. "
                            f"The cells evaluable for the three values number {n_el} of {N_CELLS}. Of these, {out['n_any_above']} have at "
                            f"least one value above its band, and {out['n_below_only']} have none above and at least one below.")
    elif n_s == 0:
        out["form"] = "b"
        tail = ("None of those cells is evaluable for space, so the spatial readouts are compared at none of them." if n_es == 0 else
                f"At none of the {n_es} of those cells that are evaluable for space do its spatial readouts also lie inside theirs.")
        out["statement"] = (f"The record's three common-class F values lie inside their bands at {n_l} cells of the {n_el} evaluable "
                            f"(of {N_CELLS} tested), listed below. " + tail)
    else:
        out["form"] = "c"
        out["statement"] = ("The record's three common-class F values and its four spatial readouts each lie inside their own band "
                            f"at {n_s} cells, listed below, of the {n_es} cells evaluable for both (the F values are inside at {n_l} "
                            f"of {n_el} evaluable cells, of {N_CELLS} tested).")
    return out


def n_histories(cells) -> int:
    """Distinct simulated histories behind a set of cells: one per (window, amplitude, record)."""
    return RECORDS_PER_CELL * len({(c[1], c[3]) for c in cells})


def diagnostics(table: pd.DataFrame, cell: tuple) -> dict:
    """R4 for one cell: rare classes inside and unavailable, losses from the target per common class, finite counts. Excludes nothing."""
    t = table[(table["span"] == cell[0]) & (table["window"] == cell[1]) & (table["coupling"] == cell[2]) & (table["F_W"] == cell[3])]
    t = t.set_index("stat")
    if any(s not in t.index for s in REQUIRED):
        raise ValueError(f"cell {cell} is not in the table with every required statistic")
    return {"rare_inside": int((t.loc[list(RARE), "status"] == "inside").sum()),
            "rare_unavailable": int((t.loc[list(RARE), "status"] == "unavailable").sum()),
            "n_lost": {s: int(t.loc[s, "n_lost"]) for s in F3},
            "n_finite": {s: int(t.loc[s, "n_finite"]) for s in F3 + SP}}


def _amps(xs) -> str:
    return ", ".join(f"{a:g}" for a in xs) if xs else "none"


def r1_text(r1: dict, names: dict) -> list[str]:
    """R1 as one line per common class: the amplitudes at each status at the baseline design, an empty list written as "none"."""
    return [f"{names[s]}: inside at F_W {_amps(r1[s]['inside'])}; above at {_amps(r1[s]['above'])}; below at {_amps(r1[s]['below'])}; "
            f"unavailable at {_amps(r1[s]['unavailable'])}." for s in F3]
