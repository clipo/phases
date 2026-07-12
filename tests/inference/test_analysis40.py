import importlib
import numpy as np
import pytest


@pytest.fixture(scope="module")
def a40():
    import sys
    from pathlib import Path
    ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "analyses"))
    return importlib.import_module("40_hierarchical_convergence")


def test_standardization_helpers(a40):
    # z-standardization of a trajectory returns mean~0, sd~1
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    z, mean, sd = a40._zstd(x)
    assert abs(z.mean()) < 1e-9 and abs(z.std() - 1.0) < 1e-9
    assert abs(mean - 3.5) < 1e-9 and sd > 0


@pytest.mark.data
def test_build_panel_and_ses_shapes_and_finiteness(a40):
    a07 = importlib.import_module("07_refined_empirical")
    inp = a07.prepare_inputs()
    y, se, t, labels = a40.build_panel_and_ses(inp, n_bins=6, n_boot=40, seed=0)
    assert y.shape == se.shape == (3, 6)
    assert t.shape == (6,)
    assert labels == ["neutral_departure", "fst", "spatial_boundary"]
    assert np.all(np.isfinite(y)) and np.all(se >= 0)
    b_ser, se_ser = a40.seriation_slope_and_se(inp, n_boot=40, seed=1)
    assert np.isfinite(b_ser) and se_ser >= 0
