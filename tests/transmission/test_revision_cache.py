"""A changed archaeological input must invalidate the spatial simulation cache."""
import copy
import importlib
from pathlib import Path
import sys

import numpy as np


def test_cache_fingerprint_binds_counts_order_and_parameters(monkeypatch):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analyses"))
    revision = importlib.import_module("47_revision_analysis")
    data = {"example": dict(m=np.array([[10, 20], [20, 10]]),
                            coords=np.array([[0., 1.], [1., 0.]]),
                            d=np.array([[0., 2.], [2., 0.]]), ranks=np.array([0., 1.]),
                            labels=np.array([0, 1]), parkin=None, names=["a", "b"],
                            pooled=np.array([0.5, 0.5]))}
    original = revision.fingerprint(data)
    assert original == revision.fingerprint(copy.deepcopy(data))
    for key in ["m", "ranks", "d", "labels", "pooled"]:
        changed = copy.deepcopy(data)
        changed["example"][key].flat[0] += 1
        assert original != revision.fingerprint(changed)
    monkeypatch.setitem(revision.CONFIG, "calib_reps", 7)
    assert original != revision.fingerprint(data)
