"""Screen rules on a synthetic table (spec 3.3; R2-6)."""
from mls_emergence.probe import resemblance as rs


def test_screen_cell_and_stop_rule_use_the_tolerance_on_relative_error():
    obs = {"hs": 0.5, "rich": 6.0, "ht": 0.6}
    assert rs.screen_cell({"hs": 0.54, "rich": 6.5, "ht": 0.58}, obs)["matched"]          # all within 10 percent
    c = rs.screen_cell({"hs": 0.56, "rich": 6.5, "ht": 0.58}, obs)                        # hs off by 12 percent
    assert not c["matched"] and abs(c["rel_hs"] - 0.12) < 1e-12
    assert rs.stop_rule({"hs": 0.54, "rich": 6.5, "ht": 0.58}, obs)["pass"]
    assert not rs.stop_rule({"hs": 0.56, "rich": 6.5, "ht": 0.58}, obs)["pass"]


def test_stress_cell_selection_smallest_fst_then_loss_then_grid_order():
    rows = [
        {"matched": True, "fst_median": 0.010, "loss": 0.01, "name": "a"},
        {"matched": True, "fst_median": 0.004, "loss": 0.02, "name": "b"},
        {"matched": False, "fst_median": 0.001, "loss": 0.00, "name": "c"},   # not matched: never chosen
        {"matched": True, "fst_median": 0.004, "loss": 0.01, "name": "d"},    # ties b on fst, smaller loss
        {"matched": True, "fst_median": 0.004, "loss": 0.01, "name": "e"},    # ties d exactly: grid order keeps d
    ]
    pick = rs.select_stress_cell(rows)
    assert pick["name"] == "d" and pick["grid_index"] == 3
    assert rs.select_stress_cell([{"matched": False, "fst_median": 0, "loss": 0}]) is None
