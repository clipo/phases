import importlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))


@pytest.mark.data
def test_prepare_inputs_shapes():
    a07 = importlib.import_module("07_refined_empirical")
    inp = a07.prepare_inputs()
    assert inp.M.shape[0] == len(inp.idx)
    assert inp.cc_c.shape[0] == len(inp.have_coords_ids)
    assert len(inp.ca) == len(inp.idx)
    assert set(inp.cluster_of.keys()) <= set(inp.have_coords_ids)
    assert inp.nmemb.shape[0] == len(inp.idx)
    assert inp.ca_vals.shape[0] == len(inp.idx)
    # panel builder returns the 3 panel signatures
    panel = a07.panel_for_bins(inp, 6)
    assert list(panel.columns) == ["neutral_departure", "fst", "spatial_boundary"]
