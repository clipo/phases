"""Reading rules R1 to R5 of the drift-extension spec (2026-11-14, draft 3, section 6) on synthetic status tables.

Each table is built here by hand with known set sizes, so the expected counts
are literals and are not computed through `read`. Rule 3: the rival each case
rejects is named beside it.
"""
import itertools
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analyses" / "per_type"))
import drift_extension_readings as R

AMPS = (0.0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9)


def table(level="below", space="inside", rare="inside", lost=0, override=None):
    """One point's 504 cells by 14 required statistics. `override(cell)` may return {stat: status} for a cell."""
    rows = []
    for span, window, coupling, fw in itertools.product(R.SPANS, R.WINDOWS, R.COUPLINGS, AMPS):
        cell = {"span": span, "window": window, "coupling": coupling, "F_W": fw}
        st = {s: level for s in R.F3}; st.update({s: space for s in R.SP}); st.update({s: rare for s in R.RARE})
        if override:
            st.update(override(cell) or {})
        for s, status in st.items():
            rows.append({**cell, "stat": s, "status": status, "n_finite": 50, "n_lost": 0 if s == "D_dist" else lost})
    return pd.DataFrame(rows)


def test_the_grid_is_504_cells_and_a_short_table_raises():
    t = table()
    assert len(t) == 504 * 14
    out = R.read(t)
    assert out["n_cells"] == 504
    short = t[~((t.span == "1") & (t.window == "32LR") & (t.coupling == "stepS") & (t.F_W == 0.9))]
    with pytest.raises(ValueError):
        R.read(short)                         # 503 cells: a silent smaller denominator is the rival
    with pytest.raises(ValueError):
        R.read(t[t.stat != "D_dist"])         # a missing required statistic
    with pytest.raises(ValueError):
        R.read(pd.concat([t, t.iloc[:1]]))    # a duplicated cell-statistic row
    bad = t.copy(); bad.loc[0, "status"] = "outside"
    with pytest.raises(ValueError):
        R.read(bad)                           # a status outside the four of status3


def test_every_evaluable_cell_below_gives_form_a_with_the_below_count():
    out = R.read(table(level="below"))
    assert out["form"] == "a" and len(out["E_L"]) == 504 and out["L"] == [] and out["S"] == []
    assert out["n_any_above"] == 0 and out["n_below_only"] == 504
    assert "At no evaluable tested cell" in out["statement"]
    # Rival rejected: counting "below" as more variation than the generator produces.


def test_above_and_below_are_counted_separately_and_above_takes_precedence():
    # One class above and one below in the same cell is counted under "at least one above".
    def ov(c):
        if c["coupling"] == "c0":
            return {"F_0": "above", "F_1": "below", "F_2": "inside"}
    out = R.read(table(level="below", override=ov))
    assert out["n_any_above"] == 72 and out["n_below_only"] == 432 and out["form"] == "a"


def test_every_cell_unavailable_gives_form_a_with_an_empty_evaluable_set():
    out = R.read(table(level="unavailable"))
    assert out["form"] == "a" and out["E_L"] == [] and out["L"] == []
    assert out["n_any_above"] == 0 and out["n_below_only"] == 0
    assert "0 of 504" in out["statement"]
    # Rival rejected: treating an unavailable cell as a miss, which would report 504 evaluable cells.


def test_one_unavailable_status_removes_the_cell_from_the_evaluable_set_even_with_one_above():
    def ov(c):
        if c["coupling"] == "stepL" and c["F_W"] == 0.9:
            return {"F_0": "above", "F_1": "unavailable"}
    out = R.read(table(level="inside", override=ov))
    assert len(out["E_L"]) == 504 - 9 and len(out["L"]) == 504 - 9


def test_r1_lists_non_adjacent_amplitudes_without_filling_the_gap():
    # At d0 (span Ld, window 8LR, coupling c0) class 0 is inside at 0.03 and 0.7 only.
    def ov(c):
        if (c["span"], c["window"], c["coupling"]) == R.D0:
            return {"F_0": "inside" if c["F_W"] in (0.03, 0.7) else "above", "F_1": "above", "F_2": "unavailable" if c["F_W"] == 0.9 else "below"}
    out = R.read(table(level="below", override=ov))
    r1 = out["R1"]
    assert r1["F_0"]["inside"] == [0.03, 0.7]
    assert r1["F_0"]["above"] == [0.0, 0.01, 0.1, 0.3, 0.5, 0.9]
    assert r1["F_0"]["below"] == [] and r1["F_1"]["above"] == list(AMPS)
    assert r1["F_2"]["unavailable"] == [0.9] and len(r1["F_2"]["below"]) == 7
    # Rival rejected: reporting a threshold "from 0.03 upward", which would include 0.1 to 0.5.
    lines = R.r1_text(r1, {"F_0": "A", "F_1": "B", "F_2": "C"})
    assert lines[0] == "A: inside at F_W 0.03, 0.7; above at 0, 0.01, 0.1, 0.3, 0.5, 0.9; below at none; unavailable at none."


def test_level_reached_with_no_cell_evaluable_for_space_gives_form_b_with_m_zero():
    def ov(c):
        if c["F_W"] == 0.9 and c["span"] == "1":
            return {"F_0": "inside", "F_1": "inside", "F_2": "inside", "D_dist": "unavailable"}
    out = R.read(table(level="above", space="inside", override=ov))
    assert out["form"] == "b" and len(out["L"]) == 21 and out["E_S"] == [] and out["S"] == []
    assert len(out["L_not_evaluable_for_space"]) == 21
    assert "21 cells" in out["statement"] and "None of those cells is evaluable for space" in out["statement"]
    # Rival rejected: reading an unavailable spatial readout as "not inside" (m would be 21).


def test_level_reached_with_space_not_inside_gives_form_b_with_m_counted():
    def ov(c):
        if c["F_W"] == 0.9 and c["span"] == "1":
            return {"F_0": "inside", "F_1": "inside", "F_2": "inside", "D_dist": "above"}
    out = R.read(table(level="above", space="inside", override=ov))
    assert out["form"] == "b" and len(out["L"]) == 21 and len(out["E_S"]) == 21 and out["S"] == []
    assert "At none of the 21 of those cells" in out["statement"]


def test_s_non_empty_only_under_c0_gives_form_c_and_the_r5_table_shows_it():
    def ov(c):
        if c["F_W"] == 0.9 and c["span"] == "1":
            return {"F_0": "inside", "F_1": "inside", "F_2": "inside", "A_lat_1": "inside" if c["coupling"] == "c0" else "below"}
    out = R.read(table(level="above", space="inside", override=ov))
    assert out["form"] == "c" and len(out["L"]) == 21 and len(out["E_S"]) == 21 and len(out["S"]) == 3
    r5 = out["R5"].set_index("coupling")
    assert list(r5.index) == list(R.COUPLINGS)
    assert (r5["cells"] == 72).all() and (r5["E_L"] == 72).all() and (r5["L"] == 3).all() and (r5["E_S"] == 3).all()
    assert r5.loc["c0", "S"] == 3 and r5.drop("c0")["S"].sum() == 0
    assert "at 3 cells" in out["statement"]
    # The statement must not collapse the table into a claim about which timing is required.
    assert "required" not in out["statement"] and "suffices" not in out["statement"]
    # Histories: the 3 cells of S are the three windows at one amplitude, 3 x 50 records.
    assert R.n_histories(out["S"]) == 150 and R.n_histories(out["L"]) == 150


def test_r4_reports_rare_class_counts_and_losses_and_excludes_nothing():
    def ov(c):
        st = {"occ_3": "unavailable", "occ_4": "above"}
        if c["F_W"] == 0.5:
            st.update({"F_0": "inside", "F_1": "inside", "F_2": "inside"})
        return st
    t = table(level="above", lost=50, override=ov)
    out = R.read(t)
    assert len(out["L"]) == 63                       # every class lost in all 50 records, and no cell excluded
    d = R.diagnostics(t, out["L"][0])
    assert d["rare_inside"] == 5 and d["rare_unavailable"] == 1
    assert d["n_lost"] == {"F_0": 50, "F_1": 50, "F_2": 50} and d["n_finite"]["D_dist"] == 50
    never = R.diagnostics(table(level="inside", lost=0), ("Ld", "8LR", "c0", 0.0))
    assert never["n_lost"] == {"F_0": 0, "F_1": 0, "F_2": 0}      # no class ever lost
