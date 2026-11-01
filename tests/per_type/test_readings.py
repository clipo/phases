"""Reading rules of the regional-drift spec, section 5 (tests 7.11, 7.12, 7.13, 7.15): accounting, eligibility, statuses, separation counts."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
pct = importlib.import_module("per_class_table")
rd = importlib.import_module("readings")


def test_summary_accounting_counts_absent_fixed_lost_and_undefined():
    labels = np.array([0, 1, 1])
    normal = np.array([[5, 3, 2], [4, 4, 2], [6, 2, 2]])
    absent1 = np.array([[5, 0, 5], [4, 0, 6], [6, 0, 4]])
    fixed1 = np.array([[0, 10, 0], [0, 10, 0], [0, 10, 0]])          # class 1 fixed at every site: the pooled sample holds one class
    rows = [pct.per_class_rows(normal, labels), pct.per_class_rows(absent1, labels), pct.per_class_rows(fixed1, labels)]
    assert rows[1]["status_1"] == "absent" and rows[1]["occ_1"] == 0
    assert rows[2]["status_1"] == "fixed" and rows[2]["occ_1"] == 3 and rows[2]["status_0"] == "absent"
    assert rows[2]["reason"] == "undefined_pool" and not np.isfinite(rows[2]["fst_clusters"]) and not np.isfinite(rows[2]["site_fst"])
    assert rows[0]["reason"] == "" and np.isfinite(rows[0]["fst_clusters"])
    assert pct.fmt_stat(float("nan"), "fixed") == "fixed" and pct.fmt_stat(float("nan"), "absent") == "absent" and pct.fmt_stat(0.1234, "finite") == "0.123"
    cell = pd.DataFrame([{"F_1": v, "status_1": s, "lost_1": l} for v, s, l in
                         [(0.1, "finite", False)] * 30 + [(float("nan"), "absent", True)] * 12 + [(float("nan"), "fixed", False)] * 8])
    acc = rd.accounting(cell, 1)
    assert acc == {"n_total": 50, "n_finite": 30, "absent": 12, "fixed": 8, "lost": 12, "undefined_pool": 0}
    band = rd.band(cell["F_1"], min_finite=45)
    assert band["n_finite"] == 30 and band["eligible"] is False
    assert rd.status(0.1, band) == "unavailable"
    good = rd.band(pd.Series(np.linspace(0.0, 1.0, 50)), min_finite=45)
    assert good["eligible"] and rd.status(0.5, good) == "inside" and rd.status(2.0, good) == "outside" and rd.status(-1.0, good) == "outside"


def test_separation_counts_exclude_undefined():
    neutral = pd.Series(list(np.linspace(0.0, 0.1, 45)) + [float("nan")] * 5)
    planted = pd.Series([0.5] * 40 + [float("nan")] * 10)
    b = rd.band(neutral, min_finite=45)
    sep = rd.separation(planted, b)
    assert (sep["n_out"], sep["n_in"], sep["n_undefined"]) == (40, 0, 10)
    assert sep["text"] == "40 of 40 (10 undefined)"
    assert sep["median_eligible"] is False                                       # 40 finite planted values, under 45
    paired = rd.paired_changes(planted, neutral)
    assert paired == {"increased": 40, "decreased": 0, "tied": 0, "unavailable": 10}
    null = rd.separation(neutral, b)
    assert null["n_out"] == 10 and null["n_in"] == 35 and null["n_undefined"] == 5   # 10 of 45 under the linear quantile rule


def test_profile_f_is_sampling_free():
    from mls_emergence.probe.resemblance import profile_f
    rows = np.array([[.5, .3, .2], [.1, .6, .3], [.25, .25, .5], [.7, .2, .1]])
    f = profile_f(rows)
    for c in range(3):
        m = rows[:, c].mean(); assert np.isclose(f[c], rows[:, c].var(ddof=1) / (m * (1 - m)))
    z = profile_f(np.array([[0., 1.], [0., 1.]]))
    assert np.isnan(z).all()                                                    # mean 0 or 1: undefined
    with pytest.raises(ValueError):
        pct.moment_f(rows)                                                        # the count-based estimator refuses probability rows (one observation per site)


def test_outputs_are_per_point_and_duplicates_refused():
    assert pct.output_stem("C3", "primary") == "per_type_table_C3_primary"
    assert pct.output_stem("M", "secondary", marker_class="Barton/Kent/MPI", cluster="small") == "per_type_table_M_secondary_Barton-Kent-MPI_small"
    a = pd.DataFrame({"point": ["primary"] * 2, "F_W": [0.1, 0.1], "record": [0, 1], "generator": ["C3"] * 2})
    b = pd.DataFrame({"point": ["primary"] * 2, "F_W": [0.1, 0.1], "record": [1, 2], "generator": ["C3"] * 2})
    with pytest.raises(ValueError):
        rd.load_runs([a, b])
    c = pd.DataFrame({"point": ["secondary"] * 2, "F_W": [0.1, 0.1], "record": [0, 1], "generator": ["C3"] * 2})
    assert len(rd.load_runs([a, c])) == 4


def test_coverage_sets_and_prestated_wording():
    grid = (0.0, 0.01, 0.03, 0.1, 0.3)
    # statuses per amplitude: common classes (3) and rare classes (7), each a list of status strings
    def mk(common, rare):
        return {a: {"common": common[i], "rare": rare[i]} for i, a in enumerate(grid)}
    ins, out, un = ["inside"] * 3, ["outside"] * 3, ["unavailable"] * 3
    rin, rout = ["inside"] * 7, ["outside"] * 7
    sets = rd.coverage_sets(mk([out, out, ins, ins, out], [rin, rin, rin, rout, rout]), grid)
    assert sets["C"] == [0.03, 0.1] and sets["R"] == [0.0, 0.01, 0.03] and sets["H"] == [0.03]
    assert sets["E_C"] == list(grid) and sets["E_H"] == list(grid)
    assert rd.coverage_wording(sets).startswith("inside the bands at F_W in {0.03}; first tested value 0.03")
    disjoint = rd.coverage_sets(mk([out, out, out, ins, ins], [rin, rin, rout, rout, rout]), grid)
    assert disjoint["H"] == [] and "at no common tested amplitude" in rd.coverage_wording(disjoint)
    empty = rd.coverage_sets(mk([out] * 5, [rin] * 5), grid)
    assert "not inside the bands at any evaluable tested amplitude" in rd.coverage_wording(empty)
    unavail = rd.coverage_sets(mk([un, un, un, un, un], [rin] * 5), grid)
    assert unavail["E_C"] == [] and unavail["E_H"] == [] and "no amplitude on the grid is evaluable for both sets" in rd.coverage_wording(unavail)
    assert len(unavail["unavailable_cells"]) == 15 and "F_W 0.3, common class 2" in rd.coverage_wording(unavail)   # the ineligible cells are named (positions when no names are given)
    mixed = rd.coverage_sets(mk([["inside", "unavailable", "inside"]] + [ins] * 4, [rin] * 5), grid)
    assert mixed["E_C"] == [0.01, 0.03, 0.1, 0.3] and mixed["C"] == [0.01, 0.03, 0.1, 0.3]   # an unavailable class makes that amplitude unevaluable, not a miss
    partly = rd.coverage_sets(mk([["outside", "unavailable", "inside"]] + [ins] * 4, [rin] * 5), grid)
    assert partly["E_C"] == list(grid) and 0.0 not in partly["C"]                              # an outside class is a miss whatever the others


def test_first_m_label_distinguishes_unavailable_from_none():
    elig = {"eligible": True, "lo": 0.0, "hi": 0.1, "median": 0.05, "n_finite": 50}
    assert rd.first_m_label(elig, [(1.5, True, 0.05), (2.0, True, 0.2)]) == "2"          # first eligible median outside
    assert rd.first_m_label(elig, [(1.5, True, 0.05), (2.0, True, 0.08)]) == "none on this grid"
    assert rd.first_m_label(elig, [(1.5, False, 0.5), (2.0, False, 0.9)]) == "unavailable"   # every planted median ineligible
    assert rd.first_m_label({**elig, "eligible": False}, [(1.5, True, 0.5)]) == "unavailable"
    assert rd.first_m_label(elig, []) == "not run"


def test_coverage_wording_names_classes_when_given():
    grid = (0.0, 0.01)
    st = {a: {"common": ["unavailable"] * 3, "rare": ["inside"] * 7} for a in grid}
    sets = rd.coverage_sets(st, grid, names={"common": ["A", "B", "C"], "rare": list("defghij")})
    assert "F_W 0, A" in rd.coverage_wording(sets) and "common class 0" not in rd.coverage_wording(sets)


def test_profile_f_rejects_rows_that_are_not_probabilities():
    from mls_emergence.probe.resemblance import profile_f
    for bad in (np.array([[1.2, -0.2], [0.2, 0.8]]), np.array([[0.5, 0.4], [0.3, 0.7]]), np.array([[np.nan, 1.0], [0.5, 0.5]])):
        with pytest.raises(ValueError):
            profile_f(bad)


def test_status3_distinguishes_below_from_above():
    b = rd.band(pd.Series(np.linspace(0.0, 1.0, 50)))
    assert rd.status3(-0.1, b) == "below" and rd.status3(0.5, b) == "inside" and rd.status3(2.0, b) == "above"
    assert rd.status3(0.5, {**b, "eligible": False}) == "unavailable" and rd.status3(float("nan"), b) == "unavailable"
    assert rd.status3(b["lo"], b) == "inside" and rd.status3(b["hi"], b) == "inside"     # endpoints inclusive
