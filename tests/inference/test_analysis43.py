import importlib
import numpy as np


def _mod():
    import sys
    from pathlib import Path
    ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))
    return importlib.import_module("43_bayesian_fst")


def test_group_counts_builder_shapes():
    a43 = _mod()

    class _Inp:
        # minimal synthetic stand-in for prepare_inputs() output
        import pandas as pd
        have_coords_ids = ["a", "b", "c", "d"]
        cluster_of = {"a": 0, "b": 0, "c": 1, "d": 1}
        counts = pd.DataFrame(
            [[10, 0, 2], [8, 1, 1], [0, 9, 3], [1, 8, 2]],
            index=["a", "b", "c", "d"], columns=["t1", "t2", "t3"],
        )

    # scope="region" is the path that groups by inp.cluster_of directly, which
    # is what this synthetic fixture exercises. The basin path re-selects
    # clusters on real coordinates and validates ids against the canonical
    # membership, so it is covered against real data in test_basin_scope.py.
    gc, sizes = a43.basin_group_counts(_Inp(), scope="region")
    assert gc.shape == (2, 3)
    assert np.allclose(gc[0], [18, 1, 3])   # cluster 0 = a+b
    assert np.allclose(sizes, gc.sum(axis=1))


def test_basin_scope_rejects_ids_outside_the_membership():
    """The guard that F18 needed. Synthetic ids are in no basin, so the basin
    path must refuse rather than silently grouping whatever it was handed."""
    import pytest
    a43 = _mod()

    class _Inp:
        import pandas as pd
        have_coords_ids = ["a", "b", "c", "d"]
        cluster_of = {"a": 0, "b": 0, "c": 1, "d": 1}
        cc_c = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        counts = pd.DataFrame(
            [[10, 0, 2], [8, 1, 1], [0, 9, 3], [1, 8, 2]],
            index=["a", "b", "c", "d"], columns=["t1", "t2", "t3"],
        )

    with pytest.raises(RuntimeError, match="canonical basin membership"):
        a43.basin_group_counts(_Inp(), scope="basin")
