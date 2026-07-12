import importlib


def _mod():
    import sys
    from pathlib import Path
    ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))
    return importlib.import_module("41_hierarchical_convergence_validation")


def test_pconv_calibration_directions():
    val = _mod()
    cases = [
        {"label": "all_up", "slope_panel": 1.2, "slope_ser": 1.2},
        {"label": "all_down", "slope_panel": -1.2, "slope_ser": -1.2},
    ]
    res = val.pconv_calibration(cases, T=6, draws=300, tune=300, chains=2, seed=0)
    by = {r["label"]: r["p_convergence"] for r in res}
    assert by["all_up"] > 0.8
    assert by["all_down"] < 0.2
